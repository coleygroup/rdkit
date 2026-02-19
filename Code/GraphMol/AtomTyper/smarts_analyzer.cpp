#include "smarts_analyzer.hpp"  //header file for smarts_analyzer.cpp
#include "atom_typer.hpp"
#include "atom_typer_query_reorder.hpp"
#include <GraphMol/GraphMol.h>  // Core RDKit molecule functionality
#include <GraphMol/MolOps.h>    // For molecule operations from RDKit molecule
#include <GraphMol/SmilesParse/SmilesParse.h>  // For parsing SMILES strings
#include <GraphMol/SmilesParse/SmartsWrite.h>  // For writing SMARTS strings
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/Descriptors/MolDescriptors.h>  // For molecular descriptors like the AND/OR/NOT queries
#include <GraphMol/QueryAtom.h>  // For QueryAtom and Query functionality
#include <GraphMol/QueryBond.h>
#include <GraphMol/QueryOps.h>  // For ATOM_EQUALS_QUERY and other query types
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <Query/QueryObjects.h>   // For the base Query class
#include <Query/EqualityQuery.h>  // For Equality Query functionality
#include <sstream>  // For string stream operations, work with data types better
#include <stdexcept>  // For std exceptions
#include <algorithm>  // For std algorithms like std::max, std::min, std::sort, std::unique
#include <cctype>
#include <functional>
#include <limits>
#include <map>
#include <regex>
#include <set>
#include <unordered_map>
#include <unordered_set>

#if __has_include(<tracy/Tracy.hpp>)
#include <tracy/Tracy.hpp>
#else
#define ZoneScoped
#define ZoneScopedN(x)
#endif

// stop variable name clashing with other libraries
namespace atom_typer {

/**
 * Private implementation class PIMPL, from the header file
 * This class contains the actual implementation details for SmartsAnalyzer
 */
class SmartsAnalyzer::Impl {
 public:
  using AtomQuery = Queries::Query<int, const RDKit::Atom *, true>;
  using BondQuery = Queries::Query<int, const RDKit::Bond *, true>;

  // Returns number of logical alternatives encoded by this query (for single
  // operator) So, OR(x, y, z) returns 3, AND(x, y) returns 1, NOT(x) returns 1
  int count_query_dof_logical(
      const Queries::Query<int, const RDKit::Atom *, true> *query) {
    // no query to consider, only one logical option
    if (!query) {
      return 1;
    }

    // retrieve name of query, such as "AtomAnd", "AtomOr", "AtomNot", or
    // specific atom property
    const std::string &desc = query->getDescription();

    // OR node, go through number of branches from this node
    // RDKit uses "AtomOr" for atom-level OR queries
    if (desc == "AtomOr" || desc == "Or") {
      int sum = 0;

      // each option under OR node adds combinatorially, sums up via
      // iteration/recursion
      for (auto it = query->beginChildren(); it != query->endChildren(); ++it) {
        // calls count_query_dof recursively on each child, accumulating total
        sum += count_query_dof_logical(it->get());
      }

      return sum;
    }

    // AND node: all constraints must be satisfied (no branching at this level)
    // However, children may still have OR nodes that create branches
    // RDKit uses "AtomAnd" for atom-level AND queries
    if (desc == "AtomAnd" || desc == "And") {
      // start at 1, since if no children, AND is trivially satisfied
      // multiplies if has OR nodes in children
      int product = 1;

      // Recursively count DOF in children and multiply
      // Example: C&(H1,H2) means Carbon AND (1H OR 2H) = 2 alternatives
      for (auto it = query->beginChildren(); it != query->endChildren(); ++it) {
        product *= count_query_dof_logical(it->get());
      }

      return product;
    }

    // NOT node: negation doesn't branch (simplified - just counts as 1)
    // RDKit uses "AtomNot" for atom-level NOT queries
    if (desc == "AtomNot" || desc == "Not") {
      return 1;
    }

    // If it's not a logical operator, it's a leaf node, one option
    return 1;
  }

  int count_query_dof(
      const Queries::Query<int, const RDKit::Atom *, true> *query,
      const std::vector<PatternItem> &pattern_types, unsigned int atom_idx) {
    const int logical_dof = count_query_dof_logical(query);
    for (const auto &item : pattern_types) {
      if (item.kind != PatternItemKind::Atom) {
        continue;
      }
      if (static_cast<unsigned int>(item.atom.atom_idx) != atom_idx) {
        continue;
      }

      if (item.atom.atomic_number <= 0) {
        return logical_dof;
      }

      const int valence_dof = std::max(1, item.atom.remaining_valence + 1);
      return std::min(logical_dof, valence_dof);
    }

    return logical_dof;
  }

  // class of atom constraints
  // TODO: expand constraints in a case by case basis for different atom
  // properties
  struct AtomConstraints {
    // number of valence electrons (can calculate bond connections from this)
    int min_valence = 0;
    int max_valence = 8;  // octet rule cap (temp)

    int min_h = 0;  // no hydrogens
    int max_h = 4;  // max hydrogens (temp)

    // formal charge (case by case for organic molecules)
    int min_charge = -8;
    int max_charge = 8;

    // flag for contradiction, if true, no valid atom can satisfy constraints
    bool contradiction = false;

    // Track excluded values for NOT queries
    std::vector<int> excluded_atomic_nums;
    std::vector<int> excluded_charges;
    std::vector<int> excluded_h_counts;

    // Check if constraints are satisfied
    bool is_satisfied() const { return !contradiction; }

    // Deduplicate exclusion lists
    // makes sure there are no duplicates in the exclusion vectors and sorts
    // them
    void deduplicate_exclusions() {
      std::sort(excluded_atomic_nums.begin(), excluded_atomic_nums.end());
      excluded_atomic_nums.erase(
          std::unique(excluded_atomic_nums.begin(), excluded_atomic_nums.end()),
          excluded_atomic_nums.end());

      std::sort(excluded_charges.begin(), excluded_charges.end());
      excluded_charges.erase(
          std::unique(excluded_charges.begin(), excluded_charges.end()),
          excluded_charges.end());

      std::sort(excluded_h_counts.begin(), excluded_h_counts.end());
      excluded_h_counts.erase(
          std::unique(excluded_h_counts.begin(), excluded_h_counts.end()),
          excluded_h_counts.end());
    }
  };

  // Check rules and contradictions for AND
  // returns True if there is a contradiction, and vice versa
  // merges list of child constraints into one and sets limits accordingly
  AtomConstraints merge_and(const std::vector<AtomConstraints> &children) {
    AtomConstraints result;

    // loop through the list of child constraints
    for (const auto &c : children) {
      // Make sure constraints are tightened
      result.min_valence = std::max(result.min_valence, c.min_valence);
      result.max_valence = std::min(result.max_valence, c.max_valence);

      result.min_h = std::max(result.min_h, c.min_h);
      result.max_h = std::min(result.max_h, c.max_h);

      result.min_charge = std::max(result.min_charge, c.min_charge);
      result.max_charge = std::min(result.max_charge, c.max_charge);

      // Merge exclusion lists (union of all exclusions)
      result.excluded_atomic_nums.insert(result.excluded_atomic_nums.end(),
                                         c.excluded_atomic_nums.begin(),
                                         c.excluded_atomic_nums.end());
      result.excluded_charges.insert(result.excluded_charges.end(),
                                     c.excluded_charges.begin(),
                                     c.excluded_charges.end());
      result.excluded_h_counts.insert(result.excluded_h_counts.end(),
                                      c.excluded_h_counts.begin(),
                                      c.excluded_h_counts.end());
    }

    // Contradicting constraints
    if (result.min_valence > result.max_valence ||
        result.min_h > result.max_h || result.min_charge > result.max_charge) {
      result.contradiction = true;
      return result;
    }

    // Note: Element-specific valence limits are already set by
    // apply_element_valence() when the atomic number is known from leaf
    // queries. No additional enforcement needed here.

    // Deduplicate exclusion lists
    result.deduplicate_exclusions();

    return result;
  }

  // Check rules for OR
  AtomConstraints merge_or(const std::vector<AtomConstraints> &children) {
    // OR takes the union (most permissive) of valid branches
    AtomConstraints result;
    result.contradiction = true;  // Start assuming failure

    // Loop through children to find valid branches
    for (const auto &c : children) {
      // Only consider non-contradictory branches
      if (!c.contradiction) {
        // First valid branch sets initial constraints
        if (result.contradiction) {
          result = c;  // First valid branch
          result.contradiction = false;
        } else {
          // Expand constraints to union
          result.min_valence = std::min(result.min_valence, c.min_valence);
          result.max_valence = std::max(result.max_valence, c.max_valence);
          result.min_h = std::min(result.min_h, c.min_h);
          result.max_h = std::max(result.max_h, c.max_h);
          result.min_charge = std::min(result.min_charge, c.min_charge);
          result.max_charge = std::max(result.max_charge, c.max_charge);
        }
      }
    }

    return result;
  }

  // Check rules for NOT - tracks excluded values
  AtomConstraints apply_not(const AtomConstraints &negated,
                            const std::string &negated_desc) {
    AtomConstraints result;

    // Try to cast the negated query to extract the specific value being
    // excluded For example, "NOT AtomicNum=6" means exclude carbon

    // If specific values were constrained in the negated query, we exclude them
    if (negated_desc == "AtomicNum") {
      // If negated specified a specific element, exclude it
      // This would need the actual value, which we'd get from analyzing the
      // child
      result.excluded_atomic_nums = negated.excluded_atomic_nums;
    }

    else if (negated_desc == "AtomFormalCharge") {
      // Track excluded charges
      if (negated.min_charge == negated.max_charge) {
        result.excluded_charges.push_back(negated.min_charge);
      }
    }

    else if (negated_desc == "AtomHCount" ||
             negated_desc == "AtomTotalHCount") {
      // Track excluded H counts
      if (negated.min_h == negated.max_h) {
        result.excluded_h_counts.push_back(negated.min_h);
      }
    }

    // Copy any existing exclusions from the negated constraints
    // merges existing exclusions from negated constraints
    result.excluded_atomic_nums.insert(result.excluded_atomic_nums.end(),
                                       negated.excluded_atomic_nums.begin(),
                                       negated.excluded_atomic_nums.end());
    result.excluded_charges.insert(result.excluded_charges.end(),
                                   negated.excluded_charges.begin(),
                                   negated.excluded_charges.end());
    result.excluded_h_counts.insert(result.excluded_h_counts.end(),
                                    negated.excluded_h_counts.begin(),
                                    negated.excluded_h_counts.end());

    return result;
  }

  // analyze valence query
  // extract constraints from a query node
  // recursively traverse query tree to extract constraints
  AtomConstraints analyze_query(
      const Queries::Query<int, const RDKit::Atom *, true> *q) {
    // No query means any atom, no constraints
    if (!q) {
      return {};
    }

    // get description of query, such as "AtomAnd", "AtomOr", "AtomNot", or
    // specific atom property
    const auto &desc = q->getDescription();

    // recursively handle logical operators
    // RDKit uses "AtomAnd" for atom-level AND queries
    if (desc == "AtomAnd" || desc == "And") {
      std::vector<AtomConstraints> children;

      // gather constraints from each child
      for (auto it = q->beginChildren(); it != q->endChildren(); ++it) {
        children.push_back(analyze_query(it->get()));
      }

      // for each member of children, merge constraints
      return merge_and(children);
    }

    // OR operator
    // RDKit uses "AtomOr" for atom-level OR queries
    if (desc == "AtomOr" || desc == "Or") {
      std::vector<AtomConstraints> children;

      for (auto it = q->beginChildren(); it != q->endChildren(); ++it) {
        children.push_back(analyze_query(it->get()));
      }

      // merge constraints
      return merge_or(children);
    }

    // NOT operator
    // RDKit uses "AtomNot" for atom-level NOT queries
    if (desc == "AtomNot" || desc == "Not") {
      // analyze negated child
      auto it = q->beginChildren();
      if (it != q->endChildren()) {
        AtomConstraints negated = analyze_query(it->get());
        return apply_not(negated, it->get()->getDescription());
      }
      return AtomConstraints{};
    }

    // Leaf queries - extract values from specific query types
    AtomConstraints c;

    // Try to cast to EqualityQuery to extract value
    const auto *eq = dynamic_cast<
        const Queries::EqualityQuery<int, const RDKit::Atom *, true> *>(q);

    // If it's an equality query, extract the value and set constraints
    // accordingly
    if (eq) {
      int val = eq->getVal();

      if (desc == "AtomicNum") {
        // Set valence limits based on element
        apply_element_valence(c, val);

        // Store the atomic number so NOT queries can exclude it
        c.excluded_atomic_nums.clear();

        // Note: for positive match, we don't add to exclusions
        // Exclusions are populated by NOT operations
      }

      else if (desc == "AtomHCount" || desc == "AtomTotalHCount" ||
               desc == "AtomImplicitHCount") {
        apply_hydrogen(c, val);
      }

      else if (desc == "AtomFormalCharge") {
        apply_charge(c, val);
      } else if (desc == "AtomExplicitValence" || desc == "AtomTotalValence") {
        apply_valence(c, val);
      } else if (desc == "AtomExplicitDegree" || desc == "AtomTotalDegree") {
        // Degree is closely related to valence
        c.min_valence = val;
        c.max_valence = val;
      }
    }

    return c;
  }

  // Helper to set valence based on element
  void apply_element_valence(AtomConstraints &c, int atomic_num) {
    // Set realistic valence limits based on common elements
    switch (atomic_num) {
      case 1:  // H
        c.min_valence = 0;
        c.max_valence = 1;
        break;
      case 6:  // C
        c.min_valence = 0;
        c.max_valence = 4;
        break;
      case 7:  // N
        c.min_valence = 0;
        c.max_valence = 4;  // 3 + lone pair or positively charged
        break;
      case 8:  // O
        c.min_valence = 0;
        c.max_valence = 2;
        break;
      case 9:  // F
        c.min_valence = 0;
        c.max_valence = 1;
        break;
      case 15:  // P
        c.min_valence = 0;
        c.max_valence = 5;
        break;
      case 16:  // S
        c.min_valence = 0;
        c.max_valence = 6;
        break;
      case 17:  // Cl
        c.min_valence = 0;
        c.max_valence = 7;  // Can be hypervalent (e.g., ClO4-)
        break;
      case 35:  // Br
        c.min_valence = 0;
        c.max_valence = 7;  // Can be hypervalent
        break;
      case 53:  // I
        c.min_valence = 0;
        c.max_valence = 7;  // Can be hypervalent
        break;
      default:
        // For other elements, use periodic table groups as guide
        c.min_valence = 0;
        c.max_valence = 8;
        break;
    }
  }

  // Apply valence constraints
  void apply_valence(AtomConstraints &c, int val) {
    c.min_valence = val;
    c.max_valence = val;
  }

  void apply_hydrogen(AtomConstraints &c, int val) {
    c.min_h = val;
    c.max_h = val;
  }

  void apply_charge(AtomConstraints &c, int val) {
    c.min_charge = val;
    c.max_charge = val;
  }

  // Extract SMARTS string representation for a single query node
  std::string query_to_smarts(
      const Queries::Query<int, const RDKit::Atom *, true> *q) {
    if (!q) {
      return "*";
    }

    const std::string &desc = q->getDescription();

    // Handle AND nodes by combining all child properties
    if (desc == "AtomAnd" || desc == "And") {
      std::string element = "";
      std::string properties = "";
      bool is_aromatic = false;

      // Collect all properties from AND children
      for (auto it = q->beginChildren(); it != q->endChildren(); ++it) {
        const auto *child = it->get();
        const std::string &child_desc = child->getDescription();
        const auto *eq = dynamic_cast<
            const Queries::EqualityQuery<int, const RDKit::Atom *, true> *>(
            child);

        if (eq) {
          int val = eq->getVal();

          if (child_desc == "AtomicNum") {
            // Element symbols
            if (val == 6) {
              element = is_aromatic ? "c" : "C";
            } else if (val == 7) {
              element = is_aromatic ? "n" : "N";
            } else if (val == 8) {
              element = is_aromatic ? "o" : "O";
            } else if (val == 16) {
              element = is_aromatic ? "s" : "S";
            } else if (val == 15) {
              element = is_aromatic ? "p" : "P";
            } else {
              element = "#" + std::to_string(val);
            }
          } else if (child_desc == "AtomType") {
            int atomic_num = val & 0xFF;
            is_aromatic = (val > 255);
            if (atomic_num == 6) {
              element = is_aromatic ? "c" : "C";
            } else if (atomic_num == 7) {
              element = is_aromatic ? "n" : "N";
            } else if (atomic_num == 8) {
              element = is_aromatic ? "o" : "O";
            } else if (atomic_num == 16) {
              element = is_aromatic ? "s" : "S";
            } else if (atomic_num == 15) {
              element = is_aromatic ? "p" : "P";
            } else {
              element = "#" + std::to_string(atomic_num);
            }
          } else if (child_desc == "AtomHCount" ||
                     child_desc == "AtomTotalHCount") {
            properties += "H" + std::to_string(val);
          } else if (child_desc == "AtomFormalCharge") {
            if (val > 0) {
              properties += "+" + std::to_string(val);
            } else if (val < 0) {
              properties += std::to_string(val);
            } else {
              properties += "+0";
            }
          }
        }
      }

      if (element.empty()) {
        element = "*";
      }
      return "[" + element + properties + "]";
    }

    // Try to extract value from equality queries
    const auto *eq = dynamic_cast<
        const Queries::EqualityQuery<int, const RDKit::Atom *, true> *>(q);

    if (eq) {
      int val = eq->getVal();

      if (desc == "AtomicNum") {
        // Return element symbol in brackets
        return "[#" + std::to_string(val) + "]";
      } else if (desc == "AtomType") {
        // AtomType encodes both atomic number and aromaticity
        int atomic_num = val & 0xFF;
        bool aromatic = (val > 255);

        // Common elements
        if (atomic_num == 6) {
          return aromatic ? "[c]" : "[C]";
        }
        if (atomic_num == 7) {
          return aromatic ? "[n]" : "[N]";
        }
        if (atomic_num == 8) {
          return aromatic ? "[o]" : "[O]";
        }
        if (atomic_num == 16) {
          return aromatic ? "[s]" : "[S]";
        }
        if (atomic_num == 15) {
          return aromatic ? "[p]" : "[P]";
        }

        return "[#" + std::to_string(atomic_num) + "]";
      } else if (desc == "AtomHCount" || desc == "AtomTotalHCount") {
        return "[H" + std::to_string(val) + "]";
      } else if (desc == "AtomFormalCharge") {
        if (val > 0) {
          return "[+" + std::to_string(val) + "]";
        }
        if (val < 0) {
          return "[" + std::to_string(val) + "]";
        }
        return "[+0]";
      }
    }

    // For complex queries (NOT, etc.), return generic wildcard
    return "[*]";
  }

  // Enumerate SMARTS strings from a query node
  std::vector<std::string> enumerate_query_smarts(
      const Queries::Query<int, const RDKit::Atom *, true> *q) {
    std::vector<std::string> results;

    if (!q) {
      results.push_back("[*]");
      return results;
    }

    const std::string &desc = q->getDescription();

    // OR node: recursively enumerate all branches
    if (desc == "AtomOr" || desc == "Or") {
      for (auto it = q->beginChildren(); it != q->endChildren(); ++it) {
        auto child_variants = enumerate_query_smarts(it->get());
        results.insert(results.end(), child_variants.begin(),
                       child_variants.end());
      }
      return results;
    }

    // For non-OR nodes, convert to SMARTS string
    results.push_back(query_to_smarts(q));
    return results;
  }

  bool is_or_query(const std::string &desc) const {
    return desc == "AtomOr" || desc == "Or";
  }

  bool is_and_query(const std::string &desc) const {
    return desc == "AtomAnd" || desc == "And";
  }

  bool is_not_query(const std::string &desc) const {
    return desc == "AtomNot" || desc == "Not";
  }

  std::unique_ptr<AtomQuery> normalize_query_push_not(const AtomQuery *q,
                                                      bool negate = false) {
    (void)negate;
    if (!q) {
      return nullptr;
    }
    // Preserve NOT semantics exactly as written in the original query.
    return std::unique_ptr<AtomQuery>(q->copy());
  }

  std::vector<std::vector<const AtomQuery *>> dnf_terms(const AtomQuery *q) {
    if (!q) {
      return {{}};
    }

    const std::string desc = q->getDescription();
    if (is_or_query(desc)) {
      std::vector<std::vector<const AtomQuery *>> out;
      for (auto it = q->beginChildren(); it != q->endChildren(); ++it) {
        auto childTerms = dnf_terms(it->get());
        out.insert(out.end(), childTerms.begin(), childTerms.end());
      }
      return out;
    }

    if (is_and_query(desc)) {
      std::vector<std::vector<const AtomQuery *>> acc(1);
      for (auto it = q->beginChildren(); it != q->endChildren(); ++it) {
        auto childTerms = dnf_terms(it->get());
        std::vector<std::vector<const AtomQuery *>> next;
        for (const auto &lhs : acc) {
          for (const auto &rhs : childTerms) {
            auto merged = lhs;
            merged.insert(merged.end(), rhs.begin(), rhs.end());
            next.push_back(std::move(merged));
          }
        }
        acc = std::move(next);
      }
      return acc;
    }

    return {{q}};
  }

  std::unique_ptr<AtomQuery> build_conjunction_query(
      const std::vector<const AtomQuery *> &term) {
    if (term.empty()) {
      return std::unique_ptr<AtomQuery>(RDKit::makeAtomNullQuery());
    }
    if (term.size() == 1) {
      return std::unique_ptr<AtomQuery>(term.front()->copy());
    }

    auto out = std::unique_ptr<AtomQuery>(new RDKit::ATOM_AND_QUERY);
    out->setDescription("AtomAnd");
    out->setNegation(false);
    for (const auto *primitive : term) {
      out->addChild(AtomQuery::CHILD_TYPE(primitive->copy()));
    }
    return out;
  }

  std::vector<std::unique_ptr<AtomQuery>> expand_recursive_primitive_variants(
      const AtomQuery *primitive, int max_recursive_variants = 128,
      bool carry_atom_maps = false) {
    std::vector<std::unique_ptr<AtomQuery>> out;
    if (!primitive) {
      return out;
    }

    if (primitive->getDescription() != "RecursiveStructure") {
      out.emplace_back(primitive->copy());
      return out;
    }

    // Negated recursive queries must remain opaque.  Expanding !$(A,B,C) into
    // individual !$(A), !$(B), !$(C) is semantically equivalent as an AND
    // conjunction, but when combined with other variant dimensions via the
    // Cartesian product in enumerate_query_variants it produces an exponential
    // OR explosion that the downstream consolidation steps cannot recombine.
    // Keeping the whole !$(...) unexpanded preserves both semantics and size.
    if (primitive->getNegation()) {
      out.emplace_back(primitive->copy());
      return out;
    }

    const auto *rsq =
        dynamic_cast<const RDKit::RecursiveStructureQuery *>(primitive);
    if (!rsq || !rsq->getQueryMol()) {
      out.emplace_back(primitive->copy());
      return out;
    }

    const std::string recursive_smarts =
        RDKit::MolToSmarts(*rsq->getQueryMol());
    SmartsAnalyzer nested;
    std::vector<std::string> expanded;
    try {
      expanded = nested.enumerate_variants(
          recursive_smarts, max_recursive_variants, false, carry_atom_maps);
    } catch (const std::exception &e) {
      // If expansion fails, return the original query
      std::cout << "Warning: Failed to expand recursive SMARTS '"
                << recursive_smarts << "': " << e.what() << std::endl;
      throw std::runtime_error("Failed to expand recursive SMARTS: " +
                               std::string(e.what()));
    }
    std::set<std::string> seen;
    for (const auto &sv : expanded) {
      if (!seen.insert(sv).second) {
        continue;
      }
      std::unique_ptr<RDKit::ROMol> submol(RDKit::SmartsToMol(sv));
      if (!submol) {
        continue;
      }
      auto *q = new RDKit::RecursiveStructureQuery(
          new RDKit::ROMol(*submol, true), rsq->getSerialNumber());
      q->setNegation(primitive->getNegation());

      // Preserve the recursive anchor primitive in the same branch. Without
      // this, some recursive anchors (notably organic-subset atoms like O)
      // can serialize as unconstrained anchor paths in downstream rendering
      // (e.g. $(O-[S]) -> $(-[S])), which broadens semantics.
      AtomQuery *anchor_copy = nullptr;
      if (const auto *anchor_atom = submol->getAtomWithIdx(0)) {
        if (const auto *anchor_qatom =
                dynamic_cast<const RDKit::QueryAtom *>(anchor_atom);
            anchor_qatom && anchor_qatom->getQuery()) {
          anchor_copy = anchor_qatom->getQuery()->copy();
        } else {
          const int atomic_num = anchor_atom->getAtomicNum();
          if (atomic_num > 0) {
            const int aromatic_flag = anchor_atom->getIsAromatic() ? 1 : 0;
            anchor_copy =
                RDKit::makeAtomTypeQuery(atomic_num, aromatic_flag);
          }
        }
      }

      if (anchor_copy) {
        auto *and_q = new RDKit::ATOM_AND_QUERY;
        and_q->setDescription("AtomAnd");
        and_q->addChild(AtomQuery::CHILD_TYPE(q));
        and_q->addChild(AtomQuery::CHILD_TYPE(anchor_copy));
        out.emplace_back(and_q);
      } else {
        out.emplace_back(q);
      }
    }

    if (out.empty()) {
      out.emplace_back(primitive->copy());
    }
    return out;
  }

  std::vector<std::unique_ptr<AtomQuery>> enumerate_query_variants(
      const AtomQuery *q, bool carry_atom_maps = false) {
    std::vector<std::unique_ptr<AtomQuery>> out;
    if (!q) {
      return out;
    }

    auto normalized = normalize_query_push_not(q, false);
    if (!normalized) {
      return out;
    }

    const auto terms = dnf_terms(normalized.get());
    out.reserve(terms.size());
    for (const auto &term : terms) {
      if (term.empty()) {
        out.push_back(build_conjunction_query(term));
        continue;
      }

      std::vector<std::vector<std::unique_ptr<AtomQuery>>> option_storage;
      std::vector<std::vector<const AtomQuery *>> option_ptrs;
      option_storage.reserve(term.size());
      option_ptrs.reserve(term.size());

      for (const auto *primitive : term) {
        auto options = expand_recursive_primitive_variants(primitive, 128,
                                                           carry_atom_maps);
        std::vector<const AtomQuery *> ptrs;
        ptrs.reserve(options.size());
        for (const auto &opt : options) {
          ptrs.push_back(opt.get());
        }
        option_storage.push_back(std::move(options));
        option_ptrs.push_back(std::move(ptrs));
      }

      std::vector<const AtomQuery *> chosen(option_ptrs.size(), nullptr);
      std::function<void(size_t)> expand = [&](size_t idx) {
        if (idx == option_ptrs.size()) {
          out.push_back(build_conjunction_query(chosen));
          return;
        }
        for (const auto *cand : option_ptrs[idx]) {
          chosen[idx] = cand;
          expand(idx + 1);
        }
      };
      expand(0);
    }

    return out;
  }

  RDKit::QueryAtom *hoist_recursive_anchor_primitives_for_and_query(
      const RDKit::QueryAtom *qa, bool verbose = false) {
    if (!qa || !qa->getQuery()) {
      return nullptr;
    }

    const auto *original_query = qa->getQuery();
    if (!is_and_query(original_query->getDescription())) {
      auto *result = new RDKit::QueryAtom(*qa);
      return result;
    }

    // Build a structural signature from a query tree for deduplication.
    // Walks the query nodes rather than serializing to SMARTS text.
    std::function<std::string(const AtomQuery *)> query_signature;
    query_signature = [&query_signature](const AtomQuery *q) -> std::string {
      if (!q) {
        return "";
      }
      std::string sig = q->getDescription();
      const auto *eq = dynamic_cast<const RDKit::ATOM_EQUALS_QUERY *>(q);
      if (eq) {
        sig += "=" + std::to_string(eq->getVal());
      }
      if (q->getNegation()) {
        sig += "!";
      }
      bool first_child = true;
      for (auto it = q->beginChildren(); it != q->endChildren(); ++it) {
        sig += first_child ? "(" : ",";
        sig += query_signature(it->get());
        first_child = false;
      }
      if (!first_child) {
        sig += ")";
      }
      return sig;
    };

    // Collect signatures of ALL existing non-recursive primitives (including
    // those inside nested AND/OR sub-trees) for dedup so we don't hoist a
    // primitive that already appears somewhere in the query.
    std::set<std::string> existing_sigs;
    // Also collect atomic numbers present in the outer query (from
    // AtomAtomicNum or AtomType) so that we can suppress hoisting an
    // AtomType whose atomic number is already constrained.
    std::set<int> existing_atomic_nums;
    std::function<void(const AtomQuery *)> collect_existing_sigs;
    collect_existing_sigs = [&](const AtomQuery *node) {
      if (!node) {
        return;
      }
      if (node->getDescription() == "RecursiveStructure") {
        return;  // don't descend into recursive queries
      }
      // For leaf nodes (no children), record their signature.
      if (node->beginChildren() == node->endChildren()) {
        existing_sigs.insert(query_signature(node));
        // Track atomic numbers from AtomAtomicNum and AtomType leaves.
        const auto *eq = dynamic_cast<const RDKit::ATOM_EQUALS_QUERY *>(node);
        if (eq) {
          const std::string &d = node->getDescription();
          if (d == "AtomAtomicNum") {
            existing_atomic_nums.insert(eq->getVal());
          } else if (d == "AtomType") {
            existing_atomic_nums.insert(eq->getVal() & 0xFF);
          }
        }
        return;
      }
      // For composite nodes, recurse into children.
      for (auto it = node->beginChildren(); it != node->endChildren(); ++it) {
        collect_existing_sigs(it->get());
      }
    };
    collect_existing_sigs(original_query);

    // Scan RecursiveStructure children; extract each anchor atom's query.
    std::vector<std::unique_ptr<AtomQuery>> hoisted;
    for (auto it = original_query->beginChildren();
         it != original_query->endChildren(); ++it) {
      const auto *child = it->get();
      if (!child || child->getDescription() != "RecursiveStructure") {
        continue;
      }
      // Do not hoist from negated recursive expressions: !$(C=X) does NOT
      // imply the atom is C; it only means the atom doesn't match C=X.
      if (child->getNegation()) {
        continue;
      }

      const auto *rec =
          static_cast<const RDKit::RecursiveStructureQuery *>(child);
      const RDKit::ROMol *qmol = rec->getQueryMol();
      if (!qmol || qmol->getNumAtoms() == 0) {
        continue;
      }

      const auto *anchor_atom = qmol->getAtomWithIdx(0);
      if (!anchor_atom) {
        continue;
      }

      const AtomQuery *anchor_q = nullptr;
      std::unique_ptr<AtomQuery> synthesized_anchor_q;
      if (const auto *anchor_qatom =
              dynamic_cast<const RDKit::QueryAtom *>(anchor_atom);
          anchor_qatom && anchor_qatom->getQuery()) {
        anchor_q = anchor_qatom->getQuery();
      } else {
        // Some recursive anchors are parsed as typed atoms without an explicit
        // query node (e.g. organic-subset symbols like O).  Synthesize an
        // AtomType query so anchor constraints can still be hoisted.
        const int atomic_num = anchor_atom->getAtomicNum();
        if (atomic_num > 0) {
          const int aromatic_flag = anchor_atom->getIsAromatic() ? 1 : 0;
          synthesized_anchor_q.reset(
              RDKit::makeAtomTypeQuery(atomic_num, aromatic_flag));
          anchor_q = synthesized_anchor_q.get();
        }
      }

      if (!anchor_q) {
        continue;
      }
      const std::string &desc = anchor_q->getDescription();
      // Skip trivial / wildcard anchors (no specific constraint).
      if (desc.empty()) {
        continue;
      }

      const std::string sig = query_signature(anchor_q);
      if (!existing_sigs.count(sig)) {
        // Also suppress hoisting if the anchor's atomic number is already
        // constrained in the outer query (e.g. outer has #6 / AtomAtomicNum=6,
        // anchor is C / AtomType encoding atomic_num=6).
        bool redundant = false;
        const auto *anchor_eq =
            dynamic_cast<const RDKit::ATOM_EQUALS_QUERY *>(anchor_q);
        if (anchor_eq && !anchor_q->getNegation()) {
          const std::string &adesc = anchor_q->getDescription();
          int anchor_anum = -1;
          if (adesc == "AtomAtomicNum") {
            anchor_anum = anchor_eq->getVal();
          } else if (adesc == "AtomType") {
            anchor_anum = anchor_eq->getVal() & 0xFF;
          }
          if (anchor_anum > 0 && existing_atomic_nums.count(anchor_anum)) {
            redundant = true;
          }
        }
        if (!redundant) {
          hoisted.emplace_back(anchor_q->copy());
          existing_sigs.insert(sig);
        }
      }
    }

    if (hoisted.empty()) {
      auto *result = new RDKit::QueryAtom(*qa);
      return result;
    }

    // Rebuild the AND query: original children + hoisted anchor primitives.
    auto *new_and = new RDKit::ATOM_AND_QUERY;
    new_and->setDescription("AtomAnd");
    for (auto it = original_query->beginChildren();
         it != original_query->endChildren(); ++it) {
      new_and->addChild(AtomQuery::CHILD_TYPE((*it)->copy()));
    }
    for (auto &prim : hoisted) {
      new_and->addChild(AtomQuery::CHILD_TYPE(prim.release()));
    }

    if (verbose) {
      RDKit::QueryAtom new_token;
      new_token.setQuery(new_and->copy());
      std::cout << "[enumerate_variants] hoisted recursive-anchor primitive "
                << "atom idx=" << qa->getIdx() << " token='"
                << RDKit::SmartsWrite::GetAtomSmarts(qa) << "' -> '"
                << RDKit::SmartsWrite::GetAtomSmarts(&new_token) << "'"
                << std::endl;
    }

    auto *result = new RDKit::QueryAtom();
    result->setQuery(new_and);
    return result;
  }

  RDKit::QueryAtom *normalize_atomicnum_arom_aliphatic_and_query(
      const RDKit::QueryAtom *qa, bool verbose = false) {
    if (!qa || !qa->getQuery()) {
      return nullptr;
    }

    const auto *query = qa->getQuery();
    if (!is_and_query(query->getDescription())) {
      return new RDKit::QueryAtom(*qa);
    }

    // Scan a query node (and recurse into nested AND children) to find
    // AtomAtomicNum and AtomIsAromatic/AtomIsAliphatic leaf primitives.
    // Returns true if the node itself was fully consumed (i.e. it only
    // contained the primitives we want to collapse).
    struct ScanResult {
      int atomic_num = 0;
      bool has_atomic_num = false;
      bool has_aromatic = false;
      bool has_aliphatic = false;
    };

    // Recursive scanner: walks AND sub-trees looking for the two target
    // leaf types while collecting everything else as "remaining".
    // `remaining_out` receives copies of children that should be kept.
    // Returns scan results for the primitives found.
    std::function<ScanResult(const AtomQuery *,
                             std::vector<AtomQuery::CHILD_TYPE> &)>
        scan_and_collect;
    scan_and_collect =
        [&](const AtomQuery *node,
            std::vector<AtomQuery::CHILD_TYPE> &remaining_out) -> ScanResult {
      ScanResult sr;
      for (auto it = node->beginChildren(); it != node->endChildren(); ++it) {
        const auto *child = it->get();
        if (!child) {
          continue;
        }
        const std::string &desc = child->getDescription();

        if (desc == "AtomAtomicNum" && !child->getNegation()) {
          const auto *eq = static_cast<const RDKit::ATOM_EQUALS_QUERY *>(child);
          const int val = eq->getVal();
          if (val > 0 && !sr.has_atomic_num) {
            sr.atomic_num = val;
            sr.has_atomic_num = true;
            continue;  // consume this child
          }
        } else if (desc == "AtomIsAromatic" && !child->getNegation()) {
          if (!sr.has_aromatic && !sr.has_aliphatic) {
            sr.has_aromatic = true;
            continue;  // consume this child
          }
        } else if (desc == "AtomIsAliphatic" && !child->getNegation()) {
          if (!sr.has_aromatic && !sr.has_aliphatic) {
            sr.has_aliphatic = true;
            continue;  // consume this child
          }
        } else if (is_and_query(desc)) {
          // Recurse into nested AND queries.
          std::vector<AtomQuery::CHILD_TYPE> nested_remaining;
          ScanResult nested = scan_and_collect(child, nested_remaining);
          // Merge findings.
          if (nested.has_atomic_num && !sr.has_atomic_num) {
            sr.atomic_num = nested.atomic_num;
            sr.has_atomic_num = true;
          }
          if (nested.has_aromatic && !sr.has_aromatic && !sr.has_aliphatic) {
            sr.has_aromatic = true;
          }
          if (nested.has_aliphatic && !sr.has_aromatic && !sr.has_aliphatic) {
            sr.has_aliphatic = true;
          }
          // If the nested AND still has remaining children, keep them
          // wrapped in an AND node.
          if (!nested_remaining.empty()) {
            if (nested_remaining.size() == 1) {
              remaining_out.push_back(std::move(nested_remaining[0]));
            } else {
              auto *inner_and = new RDKit::ATOM_AND_QUERY;
              inner_and->setDescription(desc);
              for (auto &rc : nested_remaining) {
                inner_and->addChild(std::move(rc));
              }
              remaining_out.push_back(AtomQuery::CHILD_TYPE(inner_and));
            }
          }
          continue;
        }

        // Not a target primitive — keep it.
        remaining_out.push_back(AtomQuery::CHILD_TYPE(child->copy()));
      }
      return sr;
    };

    std::vector<AtomQuery::CHILD_TYPE> remaining_children;
    ScanResult sr = scan_and_collect(query, remaining_children);

    // Nothing to normalize if we don't have both an atomic-number primitive
    // and an aromatic/aliphatic primitive.
    if (!sr.has_atomic_num || (!sr.has_aromatic && !sr.has_aliphatic)) {
      return new RDKit::QueryAtom(*qa);
    }

    // Determine the aromaticity flag for makeAtomTypeQuery.
    // aromatic → 1, aliphatic → 0.
    const int aromatic_flag = sr.has_aromatic ? 1 : 0;

    // Build the replacement AtomType query that encodes both atomic number
    // and aromaticity in a single leaf.
    std::unique_ptr<RDKit::ATOM_EQUALS_QUERY> type_query(
        RDKit::makeAtomTypeQuery(sr.atomic_num, aromatic_flag));
    if (!type_query) {
      return new RDKit::QueryAtom(*qa);
    }

    if (verbose) {
      std::cout << "[normalize_atomicnum_arom_aliphatic_and_query] normalized atomic number + aromatic/"
                << "aliphatic conjunction via query tree: "
                << "AtomAtomicNum=" << sr.atomic_num
                << (sr.has_aromatic ? " + aromatic" : " + aliphatic")
                << " -> AtomType query" << std::endl;
    }

    // If there are no remaining children, the whole AND reduces to just
    // the AtomType query.
    if (remaining_children.empty()) {
      auto *result = new RDKit::QueryAtom();
      result->setQuery(type_query.release());
      return result;
    }

    // Otherwise, rebuild the AND query with the AtomType query replacing
    // the two removed children.
    auto *new_and = new RDKit::ATOM_AND_QUERY;
    new_and->setDescription("AtomAnd");
    new_and->addChild(AtomQuery::CHILD_TYPE(type_query.release()));
    for (auto &child : remaining_children) {
      new_and->addChild(std::move(child));
    }

    auto *result = new RDKit::QueryAtom();
    result->setQuery(new_and);
    return result;
  }

  bool is_or_bond_query(const std::string &desc) const {
    return desc == "BondOr" || desc == "Or";
  }

  bool is_and_bond_query(const std::string &desc) const {
    return desc == "BondAnd" || desc == "And";
  }

  bool is_not_bond_query(const std::string &desc) const {
    return desc == "BondNot" || desc == "Not";
  }

  std::unique_ptr<BondQuery> normalize_bond_query_push_not(
      const BondQuery *q, bool negate = false) {
    (void)negate;
    if (!q) {
      return nullptr;
    }
    // Preserve NOT semantics exactly as written in the original query.
    return std::unique_ptr<BondQuery>(q->copy());
  }

  std::vector<std::vector<const BondQuery *>> dnf_bond_terms(
      const BondQuery *q) {
    if (!q) {
      return {{}};
    }

    const std::string desc = q->getDescription();
    if (is_or_bond_query(desc)) {
      std::vector<std::vector<const BondQuery *>> out;
      for (auto it = q->beginChildren(); it != q->endChildren(); ++it) {
        auto childTerms = dnf_bond_terms(it->get());
        out.insert(out.end(), childTerms.begin(), childTerms.end());
      }
      return out;
    }

    if (is_and_bond_query(desc)) {
      std::vector<std::vector<const BondQuery *>> acc(1);
      for (auto it = q->beginChildren(); it != q->endChildren(); ++it) {
        auto childTerms = dnf_bond_terms(it->get());
        std::vector<std::vector<const BondQuery *>> next;
        for (const auto &lhs : acc) {
          for (const auto &rhs : childTerms) {
            auto merged = lhs;
            merged.insert(merged.end(), rhs.begin(), rhs.end());
            next.push_back(std::move(merged));
          }
        }
        acc = std::move(next);
      }
      return acc;
    }

    return {{q}};
  }

  std::unique_ptr<BondQuery> build_bond_conjunction_query(
      const std::vector<const BondQuery *> &term) {
    if (term.empty()) {
      return std::unique_ptr<BondQuery>(RDKit::makeBondNullQuery());
    }
    if (term.size() == 1) {
      return std::unique_ptr<BondQuery>(term.front()->copy());
    }

    auto out = std::unique_ptr<BondQuery>(new RDKit::BOND_AND_QUERY);
    out->setDescription("BondAnd");
    out->setNegation(false);
    for (const auto *primitive : term) {
      out->addChild(BondQuery::CHILD_TYPE(primitive->copy()));
    }
    return out;
  }

  std::vector<std::unique_ptr<BondQuery>> enumerate_bond_query_variants(
      const BondQuery *q) {
    std::vector<std::unique_ptr<BondQuery>> out;
    if (!q) {
      return out;
    }

    auto normalized = normalize_bond_query_push_not(q, false);
    if (!normalized) {
      return out;
    }

    const auto terms = dnf_bond_terms(normalized.get());
    out.reserve(terms.size());
    for (const auto &term : terms) {
      out.push_back(build_bond_conjunction_query(term));
    }

    return out;
  }

  // Helper to check if a query value is in the exclusion list
  bool is_excluded_value(
      const Queries::Query<int, const RDKit::Atom *, true> *q,
      const AtomConstraints &constraints) {
    const auto *eq = dynamic_cast<
        const Queries::EqualityQuery<int, const RDKit::Atom *, true> *>(q);
    if (!eq) {
      // Not an equality query - could be complex (AND/OR/NOT)
      // For complex queries, check recursively through children
      const std::string &desc = q->getDescription();
      if (desc == "AtomAnd" || desc == "And" || desc == "AtomOr" ||
          desc == "Or") {
        for (auto it = q->beginChildren(); it != q->endChildren(); ++it) {
          if (is_excluded_value(it->get(), constraints)) {
            return true;
          }
        }
      }
      return false;
    }

    const std::string &desc = q->getDescription();
    int val = eq->getVal();

    // Check if this specific value is excluded
    if (desc == "AtomicNum") {
      return std::find(constraints.excluded_atomic_nums.begin(),
                       constraints.excluded_atomic_nums.end(),
                       val) != constraints.excluded_atomic_nums.end();
    } else if (desc == "AtomFormalCharge") {
      return std::find(constraints.excluded_charges.begin(),
                       constraints.excluded_charges.end(),
                       val) != constraints.excluded_charges.end();
    } else if (desc == "AtomHCount" || desc == "AtomTotalHCount" ||
               desc == "AtomImplicitHCount") {
      return std::find(constraints.excluded_h_counts.begin(),
                       constraints.excluded_h_counts.end(),
                       val) != constraints.excluded_h_counts.end();
    }

    return false;
  }

  void enumerate_cartesian(
      const std::vector<
          std::vector<const Queries::Query<int, const RDKit::Atom *, true> *>>
          &atom_variants,
      size_t idx, RDKit::ROMol *mol, std::vector<std::string> &out, int max) {
    const size_t max_count =
        max > 0 ? static_cast<size_t>(max)
                : std::numeric_limits<size_t>::max();
    // Base case: all atoms processed
    if (out.size() >= max_count) {
      return;
    }

    // If we've assigned queries for all atoms, generate SMARTS
    if (idx == atom_variants.size()) {
      out.push_back(RDKit::MolToSmarts(*mol));
      return;
    }

    // Process current atom
    RDKit::Atom *atom = mol->getAtomWithIdx(idx);
    RDKit::QueryAtom *qatom = dynamic_cast<RDKit::QueryAtom *>(atom);
    if (!qatom) {
      // If not a query atom, skip to next
      enumerate_cartesian(atom_variants, idx + 1, mol, out, max);
      return;
    }

    // Save original query to restore later
    auto *original_query = qatom->getQuery();

    // If no variants for this atom, just keep the original and continue
    if (atom_variants[idx].empty()) {
      enumerate_cartesian(atom_variants, idx + 1, mol, out, max);
      return;
    }

    // Track if we found any valid variants
    bool found_valid = false;
    bool modified_query = false;

    // Iterate through all variants for this atom
    for (auto *q : atom_variants[idx]) {
      // Analyze constraints for this query variant
      AtomConstraints constraints = analyze_query(q);

      // Skip if constraints are contradictory (chemically impossible)
      if (!constraints.is_satisfied()) {
        continue;
      }

      // Check if this specific query value is in the exclusion list
      if (is_excluded_value(q, constraints)) {
        continue;
      }

      found_valid = true;

      // Only modify the query if it's different from the original
      if (q != original_query) {
        // Clone the query to avoid ownership issues with sub-nodes
        auto *cloned_query = q->copy();
        qatom->setQuery(cloned_query);
        modified_query = true;
      }
      enumerate_cartesian(atom_variants, idx + 1, mol, out, max);
      if (out.size() >= max_count) {
        break;
      }
    }

    // If no valid variants were found, all were filtered out
    // Continue with original query to at least try
    if (!found_valid) {
      qatom->setQuery(original_query);
      enumerate_cartesian(atom_variants, idx + 1, mol, out, max);
      modified_query = false;  // We set it back to original
    }

    // Restore original query only if we modified it
    if (modified_query) {
      qatom->setQuery(original_query);
    }
  }
};

// SmartsAnalyzer implementation
SmartsAnalyzer::SmartsAnalyzer() : pimpl(std::make_unique<Impl>()) {}

SmartsAnalyzer::~SmartsAnalyzer() = default;

/**
 * Find the vector of all strings of SMARTS variants w/ max
 */
std::vector<std::string> SmartsAnalyzer::enumerate_variants(
    const std::string &smarts, int max, bool verbose, bool carry_atom_maps,
    bool enumerate_bond_order) {
  ZoneScopedN("SmartsAnalyzer::enumerate_variants");
  std::vector<std::string> results;
  if (verbose) {
    std::cout << "[enumerate_variants] input='" << smarts << "' max=" << max
              << std::endl;
  }

  if (smarts.empty() || max <= 0) {
    return results;
  }
  // std::cout << "here ! " << smarts << std::endl;
  std::unique_ptr<RDKit::ROMol> mol(RDKit::SmartsToMol(smarts));

  if (!mol) {
    std::cout << "[enumerate_variants] failed to parse SMARTS: '" << smarts
              << "'" << std::endl;
    return results;
  }

  std::vector<unsigned int> input_atom_maps;
  input_atom_maps.reserve(mol->getNumAtoms());
  for (const auto *atom : mol->atoms()) {
    input_atom_maps.push_back(atom ? atom->getAtomMapNum() : 0U);
  }
  if (verbose) {
    std::cout << "[enumerate_variants] parsed mol atoms=" << mol->getNumAtoms()
              << " bonds=" << mol->getNumBonds() << std::endl;
  }

  std::vector<unsigned int> variant_atom_indices;
  std::vector<std::vector<std::unique_ptr<SmartsAnalyzer::Impl::AtomQuery>>>
      variant_storage;
  std::vector<std::vector<const SmartsAnalyzer::Impl::AtomQuery *>>
      variant_ptrs;
  std::vector<unsigned int> variant_bond_indices;
  std::vector<std::vector<std::unique_ptr<SmartsAnalyzer::Impl::BondQuery>>>
      bond_variant_storage;
  std::vector<std::vector<const SmartsAnalyzer::Impl::BondQuery *>>
      bond_variant_ptrs;

  for (const auto *atom : mol->atoms()) {
    const auto *qatom = dynamic_cast<const RDKit::QueryAtom *>(atom);
    if (!qatom || !qatom->getQuery()) {
      continue;
    }

    std::unique_ptr<RDKit::QueryAtom> query_for_variant_enumeration(
        pimpl->hoist_recursive_anchor_primitives_for_and_query(qatom, verbose));
    if (verbose) {
      std::cout << "[hoist_recursive_anchor_primitives_for_and_query_pre] "
                << RDKit::SmartsWrite::GetAtomSmarts(qatom) << "\n";
      const std::string ve_atom = query_for_variant_enumeration
                                      ? RDKit::SmartsWrite::GetAtomSmarts(
                                            query_for_variant_enumeration.get())
                                      : "None";
      std::cout << "[hoist_recursive_anchor_primitives_for_and_query] "
                << ve_atom << "\n";
    }
    const auto *norm_input = query_for_variant_enumeration
                                 ? query_for_variant_enumeration.get()
                                 : qatom;
    std::unique_ptr<RDKit::QueryAtom> query_for_symbol_normalization(
        pimpl->normalize_atomicnum_arom_aliphatic_and_query(norm_input,
                                                            verbose));
    if (verbose) {
      const std::string sn_atom =
          query_for_symbol_normalization
              ? RDKit::SmartsWrite::GetAtomSmarts(
                    query_for_symbol_normalization.get())
              : "None";
      std::cout << "[normalize_atomicnum_arom_aliphatic_and_query] " << sn_atom
                << "\n";
    }
    const SmartsAnalyzer::Impl::AtomQuery *variant_query =
        query_for_symbol_normalization &&
                query_for_symbol_normalization->getQuery()
            ? query_for_symbol_normalization->getQuery()
            : (query_for_variant_enumeration &&
                       query_for_variant_enumeration->getQuery()
                   ? query_for_variant_enumeration->getQuery()
                   : qatom->getQuery());
    auto variants =
        pimpl->enumerate_query_variants(variant_query, carry_atom_maps);
    if (verbose) {
      std::cout << "[enumerate_variants] atom idx=" << atom->getIdx()
                << " query_desc='" << qatom->getQuery()->getDescription()
                << "' variants=" << variants.size() << std::endl;
    }
    if (variants.empty()) {
      if (verbose) {
        std::cout << "Warning: No variants generated for atom index "
                  << atom->getIdx() << " with query description '"
                  << qatom->getQuery()->getDescription()
                  << "'. Using original query." << std::endl;
      }
      variants.emplace_back(qatom->getQuery()->copy());
    }

    std::vector<const SmartsAnalyzer::Impl::AtomQuery *> ptrs;
    ptrs.reserve(variants.size());
    for (const auto &v : variants) {
      ptrs.push_back(v.get());
    }

    variant_atom_indices.push_back(atom->getIdx());
    variant_storage.push_back(std::move(variants));
    variant_ptrs.push_back(std::move(ptrs));
  }

  if (enumerate_bond_order) {
    for (const auto *bond : mol->bonds()) {
      if (!bond->hasQuery() || !bond->getQuery()) {
        continue;
      }

      auto variants = pimpl->enumerate_bond_query_variants(bond->getQuery());
      if (verbose) {
        std::cout << "[enumerate_variants] bond idx=" << bond->getIdx()
                  << " query_desc='" << bond->getQuery()->getDescription()
                  << "' variants=" << variants.size() << std::endl;
      }
      if (variants.empty()) {
        variants.emplace_back(bond->getQuery()->copy());
      }

      std::vector<const SmartsAnalyzer::Impl::BondQuery *> ptrs;
      ptrs.reserve(variants.size());
      for (const auto &v : variants) {
        ptrs.push_back(v.get());
      }

      variant_bond_indices.push_back(bond->getIdx());
      bond_variant_storage.push_back(std::move(variants));
      bond_variant_ptrs.push_back(std::move(ptrs));
    }
  }

  if (variant_atom_indices.empty() && variant_bond_indices.empty()) {
    results.push_back(smarts);
    return results;
  }
  std::unordered_set<std::string> seen;
  seen.reserve(static_cast<size_t>(std::max(16, max * 2)));
  std::vector<size_t> atom_choice(variant_ptrs.size(), 0);
  std::vector<size_t> bond_choice(bond_variant_ptrs.size(), 0);

  const auto ensure_bracketed = [](std::string token) {
    if (token.empty() || token.front() != '[') {
      token = "[" + token + "]";
    }
    return token;
  };

  const auto matching_paren = [](const std::string &txt, size_t open_pos) {
    int depth = 0;
    for (size_t i = open_pos; i < txt.size(); ++i) {
      if (txt[i] == '(') {
        ++depth;
      } else if (txt[i] == ')') {
        --depth;
        if (depth == 0) {
          return i;
        }
      }
    }
    return std::string::npos;
  };

  const auto get_inner = [&](const std::string &atom_token) {
    const std::string token = ensure_bracketed(atom_token);
    if (token.size() >= 2 && token.front() == '[' && token.back() == ']') {
      return token.substr(1, token.size() - 2);
    }
    return token;
  };

  // Ring-closure digits (e.g. "1", "%10") are outside atom brackets.
  // If present in recursive SMARTS, do not flatten "$(...)" because the current
  // subtree serializer is acyclic and can misserialize closures.
  const auto has_ring_closure_marker = [](const std::string &txt) {
    int bracket_depth = 0;
    for (size_t i = 0; i < txt.size(); ++i) {
      const char ch = txt[i];
      if (ch == '[') {
        ++bracket_depth;
        continue;
      }
      if (ch == ']') {
        if (bracket_depth > 0) {
          --bracket_depth;
        }
        continue;
      }
      if (bracket_depth != 0) {
        continue;
      }
      if (ch == '%' || std::isdigit(static_cast<unsigned char>(ch))) {
        return true;
      }
    }
    return false;
  };

  // Preserve recursive SMARTS that contain wildcard/dummy atoms (e.g. [*:3]).
  // These are often introduced during typing as attachment placeholders and
  // should remain scoped inside the recursive constraint instead of being
  // flattened into top-level branches.
  const auto has_dummy_atom_token = [](const std::string &txt) {
    return txt.find("[*") != std::string::npos;
  };
  std::function<std::string(const RDKit::ROMol &, unsigned int, int,
                            std::set<unsigned int> &)>
      serialize_recursive_subtree;
  serialize_recursive_subtree =
      [&](const RDKit::ROMol &m, unsigned int root, int parent,
          std::set<unsigned int> &path) -> std::string {
    if (root >= m.getNumAtoms()) {
      return std::string();
    }
    if (path.count(root)) {
      return ensure_bracketed(
          RDKit::SmartsWrite::GetAtomSmarts(m.getAtomWithIdx(root)));
    }

    path.insert(root);
    std::string out = ensure_bracketed(
        RDKit::SmartsWrite::GetAtomSmarts(m.getAtomWithIdx(root)));
    for (const auto *nbr : m.atomNeighbors(m.getAtomWithIdx(root))) {
      const unsigned int nidx = nbr->getIdx();
      if (static_cast<int>(nidx) == parent) {
        continue;
      }
      auto *bond = m.getBondBetweenAtoms(root, nidx);
      if (!bond) {
        continue;
      }
      const std::string child =
          serialize_recursive_subtree(m, nidx, static_cast<int>(root), path);
      if (!child.empty()) {
        out +=
            "(" + RDKit::SmartsWrite::GetBondSmarts(bond, root) + child + ")";
      }
    }
    path.erase(root);
    return out;
  };

  const auto flatten_atom_token = [&](const std::string &atom_token) {
    std::string inner = get_inner(atom_token);
    std::vector<std::string> appended_tails;
    std::string recursive_anchor_inner;

    size_t search_pos = 0;
    while (true) {
      const auto qpos = inner.find("$(", search_pos);
      if (qpos == std::string::npos) {
        break;
      }

      // Preserve negated recursive expressions exactly as-is (e.g. !$(...)).
      // These are normalized later via normalize_negated_recursive_queries().
      if (qpos > 0 && inner[qpos - 1] == '!') {
        search_pos = qpos + 2;
        continue;
      }

      const size_t open_pos = qpos + 1;
      const size_t close_pos = matching_paren(inner, open_pos);
      if (close_pos == std::string::npos || close_pos <= open_pos + 1) {
        break;
      }

      const std::string recursive_smarts =
          inner.substr(open_pos + 1, close_pos - open_pos - 1);

      // Preserve recursive expressions as "$(...)" when flattening is unsafe
      // or semantically undesirable:
      // 1) cyclic recursive SMARTS (ring-closure digits outside brackets),
      // 2) typed dummy/wildcard atoms (e.g. [*:3]) that should stay local to
      //    the recursive scope.
      if (has_ring_closure_marker(recursive_smarts) ||
          has_dummy_atom_token(recursive_smarts)) {
        if (verbose && has_dummy_atom_token(recursive_smarts)) {
          std::cout << "[enumerate_variants] preserving recursive dummy scope: "
                    << "$(" << recursive_smarts << ")" << std::endl;
        }
        search_pos = close_pos + 1;
        continue;
      }

      std::unique_ptr<RDKit::ROMol> recursive_mol(
          RDKit::SmartsToMol(recursive_smarts));
      if (!recursive_mol || recursive_mol->getNumAtoms() == 0) {
        break;
      }

      // Remove '&$(' or ';$(' as one unit when present.
      size_t erase_from = qpos;
      if (erase_from > 0 &&
          (inner[erase_from - 1] == '&' || inner[erase_from - 1] == ';')) {
        --erase_from;
      }
      inner.erase(erase_from, close_pos - erase_from + 1);

      auto *anchor = recursive_mol->getAtomWithIdx(0);
      if (anchor && recursive_anchor_inner.empty()) {
        recursive_anchor_inner = get_inner(
            ensure_bracketed(RDKit::SmartsWrite::GetAtomSmarts(anchor)));
        const auto colon = recursive_anchor_inner.rfind(':');
        if (colon != std::string::npos &&
            colon + 1 < recursive_anchor_inner.size()) {
          bool all_digits = true;
          for (size_t i = colon + 1; i < recursive_anchor_inner.size(); ++i) {
            if (!std::isdigit(
                    static_cast<unsigned char>(recursive_anchor_inner[i]))) {
              all_digits = false;
              break;
            }
          }
          if (all_digits) {
            recursive_anchor_inner = recursive_anchor_inner.substr(0, colon);
          }
        }
      }
      for (const auto *nbr : recursive_mol->atomNeighbors(anchor)) {
        const unsigned int nidx = nbr->getIdx();
        auto *bond = recursive_mol->getBondBetweenAtoms(0, nidx);
        if (!bond) {
          continue;
        }
        std::set<unsigned int> path;
        path.insert(0);
        const std::string subtree =
            serialize_recursive_subtree(*recursive_mol, nidx, 0, path);
        if (!subtree.empty()) {
          appended_tails.push_back(RDKit::SmartsWrite::GetBondSmarts(bond, 0) +
                                   subtree);
        }
      }

      search_pos = erase_from;
    }

    while (inner.find("&&") != std::string::npos) {
      inner.replace(inner.find("&&"), 2, "&");
    }
    while (!inner.empty() && (inner.front() == '&' || inner.front() == ';' ||
                              inner.front() == ',')) {
      inner.erase(inner.begin());
    }
    while (!inner.empty() && (inner.back() == '&' || inner.back() == ';')) {
      inner.pop_back();
    }

    if (!recursive_anchor_inner.empty() && !inner.empty() &&
        inner.front() == ':') {
      inner = recursive_anchor_inner + inner;
    } else if (!recursive_anchor_inner.empty() && inner.empty()) {
      inner = recursive_anchor_inner;
    }

    // If recursive flattening removed the only atom primitive, keep token
    // valid. Example: "$([N...]):1" -> ":1" must become "*:1".
    if (inner.empty()) {
      inner = "*";
    } else if (!inner.empty() && inner.front() == ':') {
      inner.insert(inner.begin(), '*');
    }

    std::string out = "[" + inner + "]";
    for (const auto &tail : appended_tails) {
      out += "(" + tail + ")";
    }
    return out;
  };

  const auto flatten_recursive_tokens = [&](const std::string &candidate) {
    std::string out;
    out.reserve(candidate.size() + 16);
    for (size_t i = 0; i < candidate.size();) {
      if (candidate[i] != '[') {
        out.push_back(candidate[i++]);
        continue;
      }

      int depth = 0;
      size_t j = i;
      for (; j < candidate.size(); ++j) {
        if (candidate[j] == '[') {
          ++depth;
        } else if (candidate[j] == ']') {
          --depth;
          if (depth == 0) {
            break;
          }
        }
      }
      if (j >= candidate.size()) {
        out += candidate.substr(i);
        break;
      }

      out += flatten_atom_token(candidate.substr(i, j - i + 1));
      i = j + 1;
    }
    return out;
  };

  std::function<void(size_t)> enumerate_bonds;
  std::function<void(size_t)> enumerate_atoms;

  enumerate_bonds = [&](size_t bond_pos) {
    if (results.size() >= static_cast<size_t>(max)) {
      return;
    }

    if (bond_pos == bond_variant_ptrs.size()) {
      RDKit::RWMol variant(*mol);

      for (size_t i = 0; i < variant_ptrs.size(); ++i) {
        auto *qa = dynamic_cast<RDKit::QueryAtom *>(
            variant.getAtomWithIdx(variant_atom_indices[i]));
        if (!qa) {
          continue;
        }
        qa->setQuery(variant_ptrs[i][atom_choice[i]]->copy());
      }

      for (size_t i = 0; i < bond_variant_ptrs.size(); ++i) {
        auto *qb = variant.getBondWithIdx(variant_bond_indices[i]);
        if (!qb) {
          continue;
        }
        qb->setQuery(bond_variant_ptrs[i][bond_choice[i]]->copy());
      }

      try {
        std::string candidate =
            flatten_recursive_tokens(RDKit::MolToSmarts(variant));
        if (candidate.empty()) {
          return;
        }
        std::unique_ptr<RDKit::ROMol> check(RDKit::SmartsToMol(candidate));
        if (carry_atom_maps && check) {
          const size_t n = std::min(static_cast<size_t>(check->getNumAtoms()),
                                    input_atom_maps.size());
          for (size_t i = 0; i < n; ++i) {
            auto *dst_atom =
                check->getAtomWithIdx(static_cast<unsigned int>(i));
            if (input_atom_maps[i] != 0U && dst_atom &&
                dst_atom->getAtomMapNum() == 0U) {
              dst_atom->setAtomMapNum(input_atom_maps[i]);
            }
          }
          candidate = RDKit::MolToSmarts(*check);
        }
        if (check && seen.insert(candidate).second) {
          if (verbose) {
            std::cout << "[enumerate_variants] candidate='" << candidate << "'"
                      << std::endl;
          }
          results.push_back(candidate);
        }
      } catch (...) {
        // Ignore invalid candidate
      }
      return;
    }

    for (size_t i = 0; i < bond_variant_ptrs[bond_pos].size(); ++i) {
      bond_choice[bond_pos] = i;
      enumerate_bonds(bond_pos + 1);
      if (results.size() >= static_cast<size_t>(max)) {
        break;
      }
    }
  };

  enumerate_atoms = [&](size_t atom_pos) {
    if (results.size() >= static_cast<size_t>(max)) {
      return;
    }

    if (atom_pos == variant_ptrs.size()) {
      enumerate_bonds(0);
      return;
    }

    for (size_t i = 0; i < variant_ptrs[atom_pos].size(); ++i) {
      atom_choice[atom_pos] = i;
      enumerate_atoms(atom_pos + 1);
      if (results.size() >= static_cast<size_t>(max)) {
        break;
      }
    }
  };

  enumerate_atoms(0);
  return results;
}

std::string SmartsAnalyzer::add_atom_maps(const std::string &smarts,
                                          unsigned int start_map) {
  if (smarts.empty()) {
    throw std::runtime_error("Empty SMARTS");
  }

  std::unique_ptr<RDKit::ROMol> mol(RDKit::SmartsToMol(smarts));
  if (!mol) {
    throw std::runtime_error("Invalid SMARTS");
  }

  using AtomQuery = Queries::Query<int, const RDKit::Atom *, true>;

  unsigned int next_map = start_map;
  std::set<const RDKit::ROMol *> visited;

  std::function<void(RDKit::ROMol &, bool)> assign_maps_in_mol;
  std::function<void(const AtomQuery *)> visit_query;

  visit_query = [&](const AtomQuery *q) {
    if (!q) {
      return;
    }

    const auto *rsq = dynamic_cast<const RDKit::RecursiveStructureQuery *>(q);
    if (rsq && rsq->getQueryMol()) {
      auto *submol = const_cast<RDKit::ROMol *>(rsq->getQueryMol());
      // Recursive SMARTS: keep the first token (anchor atom, idx 0) unmapped.
      assign_maps_in_mol(*submol, true);
    }

    for (auto it = q->beginChildren(); it != q->endChildren(); ++it) {
      visit_query(it->get());
    }
  };

  assign_maps_in_mol = [&](RDKit::ROMol &target,
                           bool skip_first_atom_map_in_this_mol) {
    if (visited.count(&target)) {
      return;
    }
    visited.insert(&target);

    for (auto *atom : target.atoms()) {
      if (skip_first_atom_map_in_this_mol && atom->getIdx() == 0) {
        // Do not map the first token of a recursive expression.
        atom->setAtomMapNum(0);
        // Keep numbering deterministic/stable for downstream processing.
        // ++next_map;
      } else {
        atom->setAtomMapNum(next_map++);
      }
    }

    for (auto *atom : target.atoms()) {
      const auto *qatom = dynamic_cast<const RDKit::QueryAtom *>(atom);
      if (!qatom || !qatom->getQuery()) {
        continue;
      }
      visit_query(qatom->getQuery());
    }
  };

  assign_maps_in_mol(*mol, false);
  return RDKit::MolToSmarts(*mol);
}


//TODO: remove
int SmartsAnalyzer::calculate_dof(const std::string &smarts) {
  // nothing in SMARTS, return 0
  if (smarts.empty()) {
    return 0;
  }

  else {
    // Parse SMARTS into RDKit molecule
    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmartsToMol(smarts));
    if (!mol) {
      throw std::runtime_error("Failed to parse SMARTS: " + smarts);
    }

    AtomTyper typer;
    const auto pattern_types = typer.type_pattern_from_smarts(smarts);

    // DOF calculation logic goes here
    int dof = 1;                  // combinatorial logic, starts at 1
    const int MAX_DOF = 1000000;  // Prevent overflow

    // iterate over all the atoms in the RDKit molecule
    for (const auto atom : mol->atoms()) {
      const auto *qatom = dynamic_cast<const RDKit::QueryAtom *>(atom);
      if (!qatom) {
        continue;
      }

      int atom_dof = pimpl->count_query_dof(qatom->getQuery(), pattern_types,
                                            atom->getIdx());

      // Check for overflow before multiplying
      if (dof > 0 && atom_dof > MAX_DOF / dof) {
        // Would overflow - return max value
        return MAX_DOF;
      }

      // Multiply into total DOF
      dof *= atom_dof;
    }

    return dof;
  }
}

namespace {

std::string strip_brackets_and_map(std::string token) {
  if (token.size() >= 2 && token.front() == '[' && token.back() == ']') {
    token = token.substr(1, token.size() - 2);
  }
  const auto colon = token.rfind(':');
  if (colon != std::string::npos && colon + 1 < token.size()) {
    bool all_digits = true;
    for (size_t i = colon + 1; i < token.size(); ++i) {
      if (!std::isdigit(static_cast<unsigned char>(token[i]))) {
        all_digits = false;
        break;
      }
    }
    if (all_digits) {
      token = token.substr(0, colon);
    }
  }
  return token;
}

bool is_explicit_hydrogen_atom(const RDKit::Atom *atom) {
  if (!atom) {
    return false;
  }
  if (atom->getAtomicNum() == 1) {
    return true;
  }

  const std::string token =
      strip_brackets_and_map(RDKit::SmartsWrite::GetAtomSmarts(atom));
  if (token.size() >= 2 && token[0] == '#' &&
      std::isdigit(static_cast<unsigned char>(token[1]))) {
    size_t pos = 1;
    while (pos < token.size() &&
           std::isdigit(static_cast<unsigned char>(token[pos]))) {
      ++pos;
    }
    const int anum = std::stoi(token.substr(1, pos - 1));
    if (anum == 1) {
      return true;
    }
  }
  if (token == "H") {
    return true;
  }
  if (token.size() > 1 && token[0] == 'H') {
    return std::all_of(token.begin() + 1, token.end(), [](char c) {
      return std::isdigit(static_cast<unsigned char>(c));
    });
  }
  return false;
}

bool is_h_primitive_delim(char ch) {
  return ch == ';' || ch == '&' || ch == ',';
}

bool is_degree_primitive_delim(char ch) {
  return ch == ';' || ch == '&' || ch == ',';
}

void analyze_h_primitives(const std::string &inner, bool &has_positive_h,
                          bool &has_positive_h0, bool &has_negated_h,
                          bool &has_negated_h0) {
  has_positive_h = false;
  has_positive_h0 = false;
  has_negated_h = false;
  has_negated_h0 = false;

  int paren_depth = 0;
  for (size_t i = 0; i < inner.size(); ++i) {
    const char ch = inner[i];
    if (ch == '(') {
      ++paren_depth;
      continue;
    }
    if (ch == ')') {
      --paren_depth;
      continue;
    }
    if (paren_depth != 0) {
      continue;
    }

    // Negated H primitive: !H, !H0, !H1 ...
    if (ch == '!' && i + 1 < inner.size() && inner[i + 1] == 'H') {
      size_t j = i + 2;
      while (j < inner.size() &&
             std::isdigit(static_cast<unsigned char>(inner[j]))) {
        ++j;
      }
      if (j == inner.size() || is_h_primitive_delim(inner[j])) {
        has_negated_h = true;
        const std::string tok = inner.substr(i + 1, j - (i + 1));
        if (tok == "H0") {
          has_negated_h0 = true;
        }
      }
      i = (j > i ? j - 1 : i);
      continue;
    }

    // Positive H primitive: H, H0, H1 ...
    if (ch == 'H') {
      const bool at_start_or_delim =
          (i == 0) || is_h_primitive_delim(inner[i - 1]);
      if (!at_start_or_delim || (i > 0 && inner[i - 1] == '!')) {
        continue;
      }
      size_t j = i + 1;
      while (j < inner.size() &&
             std::isdigit(static_cast<unsigned char>(inner[j]))) {
        ++j;
      }
      if (j == inner.size() || is_h_primitive_delim(inner[j])) {
        has_positive_h = true;
        const std::string tok = inner.substr(i, j - i);
        if (tok == "H0") {
          has_positive_h0 = true;
        }
      }
      i = (j > i ? j - 1 : i);
    }
  }
}

std::string replace_positive_h0_with_h_count(const std::string &inner,
                                             unsigned int h_count) {
  std::string out;
  out.reserve(inner.size() + 8);

  int paren_depth = 0;
  for (size_t i = 0; i < inner.size(); ++i) {
    const char ch = inner[i];
    if (ch == '(') {
      ++paren_depth;
      out.push_back(ch);
      continue;
    }
    if (ch == ')') {
      --paren_depth;
      out.push_back(ch);
      continue;
    }

    if (paren_depth == 0 && ch == 'H') {
      const bool at_start_or_delim =
          (i == 0) || is_h_primitive_delim(inner[i - 1]);
      if (at_start_or_delim && !(i > 0 && inner[i - 1] == '!')) {
        size_t j = i + 1;
        while (j < inner.size() &&
               std::isdigit(static_cast<unsigned char>(inner[j]))) {
          ++j;
        }
        if (j == inner.size() || is_h_primitive_delim(inner[j])) {
          const std::string tok = inner.substr(i, j - i);
          if (tok == "H0") {
            out += "H" + std::to_string(h_count);
            i = j - 1;
            continue;
          }
        }
      }
    }

    out.push_back(ch);
  }

  return out;
}

std::string decrement_positive_degree_primitives(const std::string &inner,
                                                 unsigned int decrement_by) {
  if (decrement_by == 0) {
    return inner;
  }

  std::string out;
  out.reserve(inner.size() + 8);

  int paren_depth = 0;
  for (size_t i = 0; i < inner.size(); ++i) {
    const char ch = inner[i];
    if (ch == '(') {
      ++paren_depth;
      out.push_back(ch);
      continue;
    }
    if (ch == ')') {
      --paren_depth;
      out.push_back(ch);
      continue;
    }

    if (paren_depth == 0 && ch == 'D') {
      const bool at_start_or_delim =
          (i == 0) || is_degree_primitive_delim(inner[i - 1]);
      if (at_start_or_delim && !(i > 0 && inner[i - 1] == '!')) {
        size_t j = i + 1;
        while (j < inner.size() &&
               std::isdigit(static_cast<unsigned char>(inner[j]))) {
          ++j;
        }
        if (j == inner.size() || is_degree_primitive_delim(inner[j])) {
          unsigned int d_value = 1;
          if (j > i + 1) {
            d_value = static_cast<unsigned int>(
                std::stoul(inner.substr(i + 1, j - (i + 1))));
          }
          const unsigned int updated_d =
              (d_value > decrement_by) ? (d_value - decrement_by) : 0;
          out += "D" + std::to_string(updated_d);
          i = j - 1;
          continue;
        }
      }
    }

    out.push_back(ch);
  }

  return out;
}

std::string set_or_append_lowercase_h_count(const std::string &inner,
                                            unsigned int h_count) {
  std::string out;
  out.reserve(inner.size() + 8);

  bool replaced = false;
  int paren_depth = 0;
  for (size_t i = 0; i < inner.size(); ++i) {
    const char ch = inner[i];
    if (ch == '(') {
      ++paren_depth;
      out.push_back(ch);
      continue;
    }
    if (ch == ')') {
      --paren_depth;
      out.push_back(ch);
      continue;
    }

    if (!replaced && paren_depth == 0 && ch == 'h') {
      const bool at_start_or_delim =
          (i == 0) || is_h_primitive_delim(inner[i - 1]);
      if (at_start_or_delim && !(i > 0 && inner[i - 1] == '!')) {
        size_t j = i + 1;
        while (j < inner.size() &&
               std::isdigit(static_cast<unsigned char>(inner[j]))) {
          ++j;
        }
        if (j == inner.size() || is_h_primitive_delim(inner[j])) {
          out += "h" + std::to_string(h_count);
          i = j - 1;
          replaced = true;
          continue;
        }
      }
    }

    out.push_back(ch);
  }

  if (!replaced) {
    if (!out.empty() && out.back() != ';') {
      out += ';';
    }
    out += "h" + std::to_string(h_count);
  }

  return out;
}

size_t find_matching_paren(const std::string &txt, size_t open_pos) {
  int depth = 0;
  for (size_t i = open_pos; i < txt.size(); ++i) {
    if (txt[i] == '(') {
      ++depth;
    } else if (txt[i] == ')') {
      --depth;
      if (depth == 0) {
        return i;
      }
    }
  }
  return std::string::npos;
}

std::string remove_recursive_expression_atom_maps(const std::string &smarts) {
  const std::regex atom_re("\\[[^\\]]+\\]");

  const auto strip_atom_maps_in_text = [&](const std::string &text) {
    std::string out;
    out.reserve(text.size());

    size_t cursor = 0;
    for (std::sregex_iterator it(text.begin(), text.end(), atom_re), end;
         it != end; ++it) {
      const auto &m = *it;
      const size_t start = static_cast<size_t>(m.position());
      const size_t len = static_cast<size_t>(m.length());
      out += text.substr(cursor, start - cursor);
      out += "[" + strip_brackets_and_map(m.str()) + "]";
      cursor = start + len;
    }
    out += text.substr(cursor);

    return out;
  };

  std::string out;
  out.reserve(smarts.size());

  size_t pos = 0;
  while (pos < smarts.size()) {
    const size_t qpos = smarts.find("$(", pos);
    if (qpos == std::string::npos) {
      out += smarts.substr(pos);
      break;
    }

    out += smarts.substr(pos, qpos - pos);

    const size_t open_pos = qpos + 1;  // '(' in '$('
    const size_t close_pos = find_matching_paren(smarts, open_pos);
    if (close_pos == std::string::npos || close_pos <= open_pos + 1) {
      out += smarts.substr(qpos);
      break;
    }

    const std::string inner =
        smarts.substr(open_pos + 1, close_pos - open_pos - 1);
    out += std::string("$(") + strip_atom_maps_in_text(inner) + ")";

    pos = close_pos + 1;
  }

  return out;
}

// ---------------------------------------------------------------------------
// collapse_equivalent_recursive_or_arms
//
// Walks each atom's query tree looking for OR nodes whose children include
// RecursiveStructure queries.  Two recursive arms are considered equivalent
// when their inner SMARTS are identical after stripping atom maps via
// remove_recursive_expression_atom_maps.  Duplicate arms are removed.
//
// Additionally, any !$(...) term that is present (map-insensitively) in
// EVERY arm of an OR is factored out to an AND prefix:
//   OR(A & !$R, B & !$R)  ->  !$R & OR(A, B)
// ---------------------------------------------------------------------------
std::string collapse_equivalent_recursive_or_arms(const std::string &smarts) {
  using AtomQuery = Queries::Query<int, const RDKit::Atom *, true>;

  std::unique_ptr<RDKit::ROMol> mol(RDKit::SmartsToMol(smarts));
  if (!mol) return smarts;

  const auto is_or = [](const std::string &d) {
    return d == "AtomOr" || d == "Or";
  };
  const auto is_and = [](const std::string &d) {
    return d == "AtomAnd" || d == "And";
  };

  // Map-insensitive signature for a RecursiveStructure node (prefix "!" when
  // negated so that negated and non-negated forms never collide).
  const auto recursive_sig = [](const AtomQuery *q) -> std::string {
    if (!q || q->getDescription() != "RecursiveStructure") return "";
    const auto *rsq =
        dynamic_cast<const RDKit::RecursiveStructureQuery *>(q);
    if (!rsq || !rsq->getQueryMol()) return "";
    const std::string inner = RDKit::MolToSmarts(*rsq->getQueryMol());
    const std::string map_free =
        remove_recursive_expression_atom_maps(inner);
    return std::string(q->getNegation() ? "!(" : "(") + map_free + ")";
  };

  // Flatten nested OR into a flat arm list.
  std::function<void(const AtomQuery *, std::vector<const AtomQuery *> &)>
      flatten_or;
  flatten_or = [&](const AtomQuery *q,
                   std::vector<const AtomQuery *> &out) {
    if (!q) return;
    if (is_or(q->getDescription())) {
      for (auto it = q->beginChildren(); it != q->endChildren(); ++it)
        flatten_or(it->get(), out);
    } else {
      out.push_back(q);
    }
  };

  // Flatten nested AND into a flat leaf list.
  std::function<void(const AtomQuery *, std::vector<const AtomQuery *> &)>
      flatten_and;
  flatten_and = [&](const AtomQuery *q,
                    std::vector<const AtomQuery *> &out) {
    if (!q) return;
    if (is_and(q->getDescription())) {
      for (auto it = q->beginChildren(); it != q->endChildren(); ++it)
        flatten_and(it->get(), out);
    } else {
      out.push_back(q);
    }
  };

  // Collect negated-recursive sigs from a single OR arm.
  const auto neg_rec_sigs_in_arm =
      [&](const AtomQuery *arm)
      -> std::map<std::string, const AtomQuery *> {
    std::map<std::string, const AtomQuery *> result;
    if (!arm) return result;
    const std::string d = arm->getDescription();
    if (d == "RecursiveStructure" && arm->getNegation()) {
      const std::string sig = recursive_sig(arm);
      if (!sig.empty()) result[sig] = arm;
    } else if (is_and(d)) {
      for (auto it = arm->beginChildren(); it != arm->endChildren(); ++it) {
        const auto *child = it->get();
        if (child && child->getDescription() == "RecursiveStructure" &&
            child->getNegation()) {
          const std::string sig = recursive_sig(child);
          if (!sig.empty()) result[sig] = child;
        }
      }
    }
    return result;
  };

  // Map-stripped SMARTS signature for any query node (not just pure recursive
  // nodes).  Used to detect duplicate OR arms regardless of their structure.
  const auto arm_canonical_sig = [&](const AtomQuery *arm) -> std::string {
    if (!arm) return "";
    // Fast path: pure RecursiveStructure arm already handled by recursive_sig.
    const std::string rs = recursive_sig(arm);
    if (!rs.empty()) return rs;
    // General path: round-trip through a temporary QueryAtom so that RDKit's
    // SMARTS serialiser can render the full query tree, then strip atom maps.
    RDKit::QueryAtom qa(0);
    qa.setQuery(arm->copy());
    std::string s;
    try {
      s = RDKit::SmartsWrite::GetAtomSmarts(&qa);
      s = remove_recursive_expression_atom_maps(s);
    } catch (...) {
      s = "";
    }
    return s;
  };

  // Process one OR node.  Returns a new node (caller owns) or nullptr when
  // nothing changed.
  const auto process_or = [&](const AtomQuery *or_q) -> AtomQuery * {
    std::vector<const AtomQuery *> arms;
    flatten_or(or_q, arms);
    if (arms.size() < 2) return nullptr;

    // Step 1 – deduplicate ALL OR arms by map-insensitive SMARTS signature,
    // not just arms that are purely RecursiveStructure nodes.  This catches
    // duplicates of the form AND(primitive, !$(X)) as well.
    std::vector<const AtomQuery *> deduped;
    std::set<std::string> seen;
    bool any_dup = false;
    for (const auto *arm : arms) {
      const std::string sig = arm_canonical_sig(arm);
      if (!sig.empty()) {
        if (!seen.insert(sig).second) {
          any_dup = true;
          continue;
        }
      }
      deduped.push_back(arm);
    }

    // Sort deduped arms by canonical signature for consistent ordering.
    std::sort(deduped.begin(), deduped.end(),
              [&](const AtomQuery *a, const AtomQuery *b) {
                return arm_canonical_sig(a) < arm_canonical_sig(b);
              });

    // Step 2 – find negated-recursive sigs common to ALL remaining arms.
    std::map<std::string, const AtomQuery *> common;
    bool first = true;
    for (const auto *arm : deduped) {
      auto m = neg_rec_sigs_in_arm(arm);
      if (first) {
        common = m;
        first = false;
      } else {
        for (auto it = common.begin(); it != common.end();) {
          it = m.count(it->first) ? std::next(it) : common.erase(it);
        }
      }
    }

    // Check if sorting changed the arm order.
    bool order_changed = false;
    if (!any_dup && deduped.size() == arms.size()) {
      for (size_t i = 0; i < deduped.size(); ++i) {
        if (deduped[i] != arms[i]) { order_changed = true; break; }
      }
    }

    if (!any_dup && !order_changed && common.empty()) return nullptr;

    // Step 3 – build new OR, removing common neg-recursive from each arm.
    auto *new_or = new RDKit::ATOM_OR_QUERY;
    new_or->setDescription("AtomOr");
    for (const auto *arm : deduped) {
      if (common.empty()) {
        new_or->addChild(AtomQuery::CHILD_TYPE(arm->copy()));
        continue;
      }
      const std::string d = arm->getDescription();
      if (is_and(d)) {
        auto *new_and = new RDKit::ATOM_AND_QUERY;
        new_and->setDescription(d);
        for (auto it = arm->beginChildren(); it != arm->endChildren(); ++it) {
          const auto *child = it->get();
          if (child && child->getDescription() == "RecursiveStructure" &&
              child->getNegation() &&
              common.count(recursive_sig(child))) {
            continue;  // stripped — will be re-added above the OR
          }
          new_and->addChild(AtomQuery::CHILD_TYPE(child->copy()));
        }
        const auto nc =
            std::distance(new_and->beginChildren(), new_and->endChildren());
        if (nc == 0) {
          delete new_and;
          new_or->addChild(AtomQuery::CHILD_TYPE(arm->copy()));
        } else if (nc == 1) {
          new_or->addChild(
              AtomQuery::CHILD_TYPE(new_and->beginChildren()->get()->copy()));
          delete new_and;
        } else {
          new_or->addChild(AtomQuery::CHILD_TYPE(new_and));
        }
      } else {
        new_or->addChild(AtomQuery::CHILD_TYPE(arm->copy()));
      }
    }

    if (common.empty()) return new_or;

    // Wrap: AND(common_neg_rec..., new_or)
    auto *new_and = new RDKit::ATOM_AND_QUERY;
    new_and->setDescription("AtomAnd");
    for (const auto &[sig, node] : common) {
      new_and->addChild(AtomQuery::CHILD_TYPE(node->copy()));
    }
    new_and->addChild(AtomQuery::CHILD_TYPE(new_or));
    return new_and;
  };

  bool changed = false;
  for (auto *atom : mol->atoms()) {
    auto *qatom = dynamic_cast<RDKit::QueryAtom *>(atom);
    if (!qatom || !qatom->getQuery()) continue;

    AtomQuery *q = qatom->getQuery();
    const std::string desc = q->getDescription();

    if (is_or(desc)) {
      AtomQuery *r = process_or(q);
      if (r) {
        qatom->setQuery(r);
        changed = true;
      }
    } else if (is_and(desc)) {
      std::vector<std::pair<size_t, AtomQuery *>> replacements;
      size_t idx = 0;
      for (auto it = q->beginChildren(); it != q->endChildren(); ++it, ++idx) {
        const auto *child = it->get();
        if (child && is_or(child->getDescription())) {
          AtomQuery *r = process_or(child);
          if (r) replacements.push_back({idx, r});
        }
      }
      if (!replacements.empty()) {
        auto *new_and = new RDKit::ATOM_AND_QUERY;
        new_and->setDescription(desc);
        size_t ci = 0, ri = 0;
        for (auto it = q->beginChildren(); it != q->endChildren();
             ++it, ++ci) {
          if (ri < replacements.size() && replacements[ri].first == ci) {
            new_and->addChild(
                AtomQuery::CHILD_TYPE(replacements[ri].second));
            ++ri;
          } else {
            new_and->addChild(AtomQuery::CHILD_TYPE(it->get()->copy()));
          }
        }
        qatom->setQuery(new_and);
        changed = true;
      }
    }
  }

  // --- Second pass: deduplicate RecursiveStructure children within AND
  //     nodes.  Typing + carry-through can produce duplicate $() or !$()
  //     expressions that inflate the output (and weaken negated exclusions).
  //     We flatten nested binary AND trees (produced by expandQuery) so that
  //     duplicates hidden in inner AND nodes are also caught.
  for (auto *atom : mol->atoms()) {
    auto *qatom = dynamic_cast<RDKit::QueryAtom *>(atom);
    if (!qatom || !qatom->getQuery()) continue;

    AtomQuery *q = qatom->getQuery();
    if (!is_and(q->getDescription())) continue;

    // Flatten the (possibly nested binary) AND tree into a flat leaf list.
    std::vector<const AtomQuery *> leaves;
    flatten_and(q, leaves);

    std::set<std::string> seen_rec_sigs;
    std::set<size_t> dup_indices;
    for (size_t i = 0; i < leaves.size(); ++i) {
      const auto *child = leaves[i];
      if (!child || child->getDescription() != "RecursiveStructure") continue;
      const std::string sig = recursive_sig(child);
      if (sig.empty()) continue;
      if (!seen_rec_sigs.insert(sig).second) {
        dup_indices.insert(i);  // duplicate
      }
    }

    if (!dup_indices.empty()) {
      auto *new_and = new RDKit::ATOM_AND_QUERY;
      new_and->setDescription("AtomAnd");
      for (size_t i = 0; i < leaves.size(); ++i) {
        if (!dup_indices.count(i)) {
          new_and->addChild(AtomQuery::CHILD_TYPE(leaves[i]->copy()));
        }
      }
      qatom->setQuery(new_and);
      changed = true;
    }
  }

  if (!changed) return smarts;
  return RDKit::MolToSmarts(*mol);
}

// ---------------------------------------------------------------------------
// factor_or_common_primitives
//
// For each atom in a SMARTS molecule, look for OR nodes in the query tree
// whose children (AND-terms) share a common primitive type (e.g.
// AtomHCount) across every alternative.  When removing that primitive
// leaves identical residuals for every value, the OR is factored.
//
// This version:
//  1) Flattens nested OR nodes into a single list of alternatives.
//  2) Separates alternatives into factorable (all-leaves) vs unfactorable
//     (containing RecursiveStructure / NOT / sub-OR).
//  3) Groups factorable terms by residual and factors per-group — groups
//     with only one term are left as-is rather than blocking the whole OR.
//
//   OR(+&D1&H0, +&D1&H1, +&D2&H0, +&D2&H1)
//     →  OR( AND(OR(H0,H1), +, D1),  AND(OR(H0,H1), +, D2) )
// ---------------------------------------------------------------------------
std::string factor_or_common_primitives(
  const std::string &smarts, bool verbose,
  const std::vector<std::string> &user_factoring_priority = {}) {
  using AtomQuery = Queries::Query<int, const RDKit::Atom *, true>;
  static thread_local int recursive_pass_depth = 0;

  if (verbose) {
    std::cout << "\n===== factor_or_common_primitives =====\n"
              << "  input SMARTS: " << smarts << "\n";
  }

  const auto is_and = [](const std::string &d) {
    return d == "AtomAnd" || d == "And";
  };
  const auto is_or = [](const std::string &d) {
    return d == "AtomOr" || d == "Or";
  };
  const auto canonical_factoring_category = [](const std::string &text) {
    std::string key;
    key.reserve(text.size());
    for (char c : text) {
      key.push_back(static_cast<char>(
          std::tolower(static_cast<unsigned char>(c))));
    }

    if (key == "atomtype" || key == "atomatomicnum" || key == "atomicnum" ||
        key == "symbol") {
      return std::string("symbol");
    }
    if (key == "atomaromatic" || key == "atomisaromatic" ||
        key == "atomaliphatic" || key == "atomisaliphatic" ||
        key == "aromaticity") {
      return std::string("aromaticity");
    }
    if (key == "atomformalcharge" || key == "charge") {
      return std::string("charge");
    }
    if (key == "atomtotaldegree" || key == "x") {
      return std::string("x");
    }
    if (key == "atomdegree" || key == "atomexplicitdegree" || key == "d") {
      return std::string("d");
    }
    if (key == "atomhcount" || key == "atomtotalhcount" ||
        key == "atomimplicithcount" || key == "h") {
      return std::string("h");
    }
    if (key == "recursive" || key == "recursion") {
      return std::string("recursion");
    }
    return std::string();
  };

  std::vector<std::string> category_order = {
      "symbol", "aromaticity", "charge", "x", "d", "h", "recursion"};
  if (!user_factoring_priority.empty()) {
    std::vector<std::string> custom;
    custom.reserve(category_order.size());
    for (const auto &raw : user_factoring_priority) {
      const auto cat = canonical_factoring_category(raw);
      if (cat.empty()) {
        continue;
      }
      if (std::find(custom.begin(), custom.end(), cat) == custom.end()) {
        custom.push_back(cat);
      }
    }
    if (!custom.empty()) {
      for (const auto &def_cat : category_order) {
        if (std::find(custom.begin(), custom.end(), def_cat) == custom.end()) {
          custom.push_back(def_cat);
        }
      }
      category_order = std::move(custom);
    }
  }

  std::map<std::string, int> category_rank;
  int rank = static_cast<int>(category_order.size());
  for (const auto &cat : category_order) {
    category_rank[cat] = rank--;
  }

  const auto factoring_primitive_priority = [&](const std::string &d) {
    const auto cat = canonical_factoring_category(d);
    auto it = category_rank.find(cat);
    if (it != category_rank.end()) {
      return it->second;
    }
    return 0;
  };

  // Helper: render a query sub-tree as a SMARTS atom expression.
  const auto query_to_smarts_str = [](const AtomQuery *q) -> std::string {
    if (!q) return "(null)";
    RDKit::QueryAtom tmp;
    tmp.setQuery(q->copy());
    return RDKit::SmartsWrite::GetAtomSmarts(&tmp);
  };

  // Leaf-primitive signature for structural comparison.
  struct PrimSig {
    std::string desc;
    int val;
    bool neg;
    bool operator<(const PrimSig &o) const {
      if (desc != o.desc) return desc < o.desc;
      if (val != o.val) return val < o.val;
      return static_cast<int>(neg) < static_cast<int>(o.neg);
    }
    bool operator==(const PrimSig &o) const {
      return desc == o.desc && val == o.val && neg == o.neg;
    }
    bool operator!=(const PrimSig &o) const { return !(*this == o); }
  };

  // Leaf data: signature + pointer to original query node (for cloning).
  struct LeafInfo {
    PrimSig sig;
    const AtomQuery *node;
  };

  // Flatten an AND term (or single leaf) into its leaf primitives.
  // Returns false when the term contains non-factorizable sub-queries
  // (RecursiveStructure, NOT, nested OR).
  std::function<bool(const AtomQuery *, std::vector<LeafInfo> &)>
      collect_leaves;
  collect_leaves = [&](const AtomQuery *q, std::vector<LeafInfo> &out) -> bool {
    if (!q) return false;
    const std::string &d = q->getDescription();
    if (is_and(d)) {
      for (auto it = q->beginChildren(); it != q->endChildren(); ++it) {
        if (!collect_leaves(it->get(), out)) return false;
      }
      return true;
    }
    if (is_or(d) || d == "RecursiveStructure" || d == "AtomNot" || d == "Not") {
      return false;
    }
    const auto *eq = dynamic_cast<const RDKit::ATOM_EQUALS_QUERY *>(q);
    if (!eq) return false;
    out.push_back({{d, eq->getVal(), q->getNegation()}, q});
    return true;
  };

  // Recursively flatten nested ORs into a single list of leaf alternatives.
  std::function<void(const AtomQuery *, std::vector<const AtomQuery *> &)>
      flatten_or;
  flatten_or = [&](const AtomQuery *q, std::vector<const AtomQuery *> &out) {
    if (!q) return;
    if (is_or(q->getDescription())) {
      for (auto it = q->beginChildren(); it != q->endChildren(); ++it) {
        flatten_or(it->get(), out);
      }
    } else {
      out.push_back(q);
    }
  };

  // ------------------------------------------------------------------
  // try_factor — attempt per-group factoring of an OR node.
  // Returns a new query (caller owns) or nullptr if nothing was factored.
  // ------------------------------------------------------------------
  auto try_factor = [&](const AtomQuery *or_node) -> AtomQuery * {
    // ---- Step 1: Flatten nested ORs ----
    std::vector<const AtomQuery *> flat_terms;
    flatten_or(or_node, flat_terms);

    if (verbose) {
      std::cout << "  try_factor: " << flat_terms.size()
                << " flat alternatives after OR-flattening\n";
      for (size_t i = 0; i < flat_terms.size(); ++i) {
        std::cout << "    [" << i << "] (" << flat_terms[i]->getDescription()
                  << ") " << query_to_smarts_str(flat_terms[i]) << "\n";
      }
    }

    if (flat_terms.size() < 2) {
      if (verbose) std::cout << "    < 2 flat terms, nothing to factor\n";
      return nullptr;
    }

    // ---- Step 2: Separate factorable vs unfactorable terms ----
    struct TermData {
      std::vector<LeafInfo> leaves;  // sorted by sig
      const AtomQuery *original;
    };
    std::vector<TermData> factorable;
    std::vector<const AtomQuery *> unfactorable;

    for (size_t i = 0; i < flat_terms.size(); ++i) {
      TermData td;
      td.original = flat_terms[i];
      if (collect_leaves(flat_terms[i], td.leaves) && !td.leaves.empty()) {
        std::sort(
            td.leaves.begin(), td.leaves.end(),
            [](const LeafInfo &a, const LeafInfo &b) { return a.sig < b.sig; });
        factorable.push_back(std::move(td));
        if (verbose) {
          std::cout << "      term " << i << " -> factorable, "
                    << factorable.back().leaves.size() << " leaves: {";
          for (const auto &l : factorable.back().leaves)
            std::cout << " " << l.sig.desc << "=" << l.sig.val;
          std::cout << " }\n";
        }
      } else {
        unfactorable.push_back(flat_terms[i]);
        if (verbose) {
          std::cout << "      term " << i
                    << " -> unfactorable (recursive/NOT/sub-OR), kept as-is\n";
        }
      }
    }

    if (factorable.size() < 2) {
      if (verbose) std::cout << "    < 2 factorable terms, nothing to factor\n";
      return nullptr;
    }

    // ---- Step 3: Find candidate primitive types ----
    // A candidate must: (a) appear in at least 2 factorable terms,
    // (b) have at most 1 occurrence per term, (c) have 2+ distinct values.
    std::map<std::string, std::set<int>> prim_vals;  // desc -> values seen
    std::map<std::string, size_t>
        prim_term_count;  // desc -> # terms containing it
    for (const auto &td : factorable) {
      std::map<std::string, int> counts;
      for (const auto &l : td.leaves) {
        counts[l.sig.desc]++;
        prim_vals[l.sig.desc].insert(l.sig.val);
      }
      for (const auto &[d, c] : counts) {
        if (c == 1) prim_term_count[d]++;
      }
    }

    std::vector<std::string> candidates;
    for (const auto &[desc, vals] : prim_vals) {
      if (vals.size() < 2 || prim_term_count[desc] < 2) continue;
      // Verify exactly-one-per-term for every term that has it:
      bool ok = true;
      for (const auto &td : factorable) {
        int cnt = 0;
        for (const auto &l : td.leaves)
          if (l.sig.desc == desc) ++cnt;
        if (cnt > 1) {
          ok = false;
          break;
        }
      }
      if (ok) candidates.push_back(desc);
    }

    if (verbose) {
      std::cout << "    candidate primitives to factor: {";
      for (const auto &c : candidates)
        std::cout << " " << c << "(" << prim_vals[c].size() << " vals)";
      std::cout << " }\n";
    }
    if (candidates.empty()) {
      if (verbose) std::cout << "    no candidates, nothing to factor\n";
      return nullptr;
    }

    // ---- Step 4: Score each candidate by term reduction ----
    // For each candidate, group factorable terms by residual (all prims
    // except the candidate). Groups with 2+ terms can be factored.
    // Reduction = sum over groups of (group_size - 1) for groups with 2+.

    struct GroupEntry {
      int cand_val;
      size_t term_idx;  // index into `factorable`
      const AtomQuery *cand_node;
    };
    struct ResidualGroup {
      std::vector<PrimSig> residual_sigs;
      std::vector<GroupEntry> entries;
    };

    std::string best_cand;
    int best_reduction = 0;
    int best_priority = std::numeric_limits<int>::min();

    for (const auto &cand : candidates) {
      std::map<std::vector<PrimSig>, ResidualGroup> groups;
      size_t terms_without = 0;

      for (size_t i = 0; i < factorable.size(); ++i) {
        const auto &td = factorable[i];
        bool found = false;
        int val = 0;
        const AtomQuery *cand_node = nullptr;
        std::vector<PrimSig> residual;
        for (const auto &l : td.leaves) {
          if (l.sig.desc == cand && !found) {
            val = l.sig.val;
            cand_node = l.node;
            found = true;
          } else {
            residual.push_back(l.sig);
          }
        }
        if (!found) {
          ++terms_without;
          continue;
        }
        std::sort(residual.begin(), residual.end());
        auto &g = groups[residual];
        g.residual_sigs = residual;
        g.entries.push_back({val, i, cand_node});
      }

      int reduction = 0;
      for (const auto &[res, grp] : groups) {
        if (grp.entries.size() >= 2) {
          // Count distinct values in this group.
          std::set<int> vals;
          for (const auto &e : grp.entries) vals.insert(e.cand_val);
          if (vals.size() >= 2) {
            // N terms collapse to 1; reduction = N - 1
            reduction += static_cast<int>(grp.entries.size()) - 1;
          }
        }
      }

      if (verbose) {
        std::cout << "    candidate '" << cand << "': " << groups.size()
                  << " residual groups, " << terms_without
                  << " terms without, reduction=" << reduction << "\n";
        for (const auto &[res, grp] : groups) {
          std::cout << "      group {";
          for (const auto &s : res) std::cout << " " << s.desc << "=" << s.val;
          std::cout << " } -> " << grp.entries.size() << " terms, vals={";
          for (const auto &e : grp.entries) std::cout << " " << e.cand_val;
          std::cout << " }\n";
        }
      }

      const int cand_priority = factoring_primitive_priority(cand);
      if (reduction > best_reduction ||
          (reduction == best_reduction && cand_priority > best_priority) ||
          (reduction == best_reduction && cand_priority == best_priority &&
           (best_cand.empty() || cand < best_cand))) {
        best_reduction = reduction;
        best_cand = cand;
        best_priority = cand_priority;
      }
    }

    if (best_reduction == 0) {
      if (verbose)
        std::cout << "    no candidate yields a reduction, nothing to factor\n";
      return nullptr;
    }

    if (verbose) {
      std::cout << "    >>> best candidate: '" << best_cand
                << "' with reduction=" << best_reduction << "\n";
    }

    // ---- Step 5: Build the factored OR using the best candidate ----
    // Recompute groups for the best candidate.
    std::map<std::vector<PrimSig>, ResidualGroup> groups;
    std::vector<size_t> terms_without;  // factorable indices without cand

    for (size_t i = 0; i < factorable.size(); ++i) {
      const auto &td = factorable[i];
      bool found = false;
      int val = 0;
      const AtomQuery *cand_node = nullptr;
      std::vector<PrimSig> residual;
      for (const auto &l : td.leaves) {
        if (l.sig.desc == best_cand && !found) {
          val = l.sig.val;
          cand_node = l.node;
          found = true;
        } else {
          residual.push_back(l.sig);
        }
      }
      if (!found) {
        terms_without.push_back(i);
        continue;
      }
      std::sort(residual.begin(), residual.end());
      auto &g = groups[residual];
      g.residual_sigs = residual;
      g.entries.push_back({val, i, cand_node});
    }

    // Build the new OR.
    auto *new_or = new RDKit::ATOM_OR_QUERY;
    new_or->setDescription("AtomOr");

    for (auto &[res_sigs, grp] : groups) {
      // Count distinct values.
      std::set<int> vals;
      for (const auto &e : grp.entries) vals.insert(e.cand_val);

      if (grp.entries.size() >= 2 && vals.size() >= 2) {
        // Factor this group: AND( OR(cand_vals…), residual… )
        // Deduplicate candidate values (keep first node per value).
        std::map<int, const AtomQuery *> val_to_node;
        for (const auto &e : grp.entries) {
          if (val_to_node.find(e.cand_val) == val_to_node.end())
            val_to_node[e.cand_val] = e.cand_node;
        }

        // Get residual leaf nodes from the first entry.
        const auto &first_term = factorable[grp.entries[0].term_idx];
        std::vector<const AtomQuery *> residual_nodes;
        for (const auto &l : first_term.leaves) {
          if (l.sig.desc != best_cand) {
            residual_nodes.push_back(l.node);
          }
        }

        // Keep explicit OR values as-is. Inferring a complete "universe" from
        // observed values is unsound and can broaden queries
        // (e.g. [F,Cl,Br,I] -> [*]).
        auto *or_vals = new RDKit::ATOM_OR_QUERY;
        or_vals->setDescription("AtomOr");
        for (const auto &[v, nd] : val_to_node) {
          or_vals->addChild(AtomQuery::CHILD_TYPE(nd->copy()));
        }

        AtomQuery *and_node = nullptr;
        if (residual_nodes.empty()) {
          // No residuals — just the OR of candidate values.
          and_node = or_vals;
        } else {
          auto *a = new RDKit::ATOM_AND_QUERY;
          a->setDescription("AtomAnd");
          a->addChild(AtomQuery::CHILD_TYPE(or_vals));
          for (const auto *n : residual_nodes) {
            a->addChild(AtomQuery::CHILD_TYPE(n->copy()));
          }
          and_node = a;
        }

        if (verbose) {
          std::cout << "    factored group {";
          for (const auto &s : res_sigs)
            std::cout << " " << s.desc << "=" << s.val;
          std::cout << " }: " << grp.entries.size() << " terms -> "
                    << query_to_smarts_str(and_node) << "\n";
        }

        // If the simplified node is a plain AND (no nested OR), it can
        // be added directly without a recursive wrapper.
        const bool needs_recursive_wrap = true;  // has OR>AND>OR

        if (needs_recursive_wrap) {
          // Wrap in $([…]) to avoid the OR>AND>OR serialization issue.
          RDKit::QueryAtom tmp_qa;
          tmp_qa.setQuery(and_node);  // takes ownership
          const std::string inner_smarts =
              RDKit::SmartsWrite::GetAtomSmarts(&tmp_qa);
          std::unique_ptr<RDKit::ROMol> inner_mol(
              RDKit::SmartsToMol(inner_smarts));
          if (inner_mol) {
            auto *rsq = new RDKit::RecursiveStructureQuery(
                new RDKit::ROMol(*inner_mol, true));
            rsq->setDescription("RecursiveStructure");
            new_or->addChild(AtomQuery::CHILD_TYPE(rsq));
          } else {
            // Fallback: keep original terms.
            for (const auto &e : grp.entries) {
              new_or->addChild(AtomQuery::CHILD_TYPE(
                  factorable[e.term_idx].original->copy()));
            }
          }
        } else {
          new_or->addChild(AtomQuery::CHILD_TYPE(and_node));
        }
      } else {
        // Single-term group or single value — keep original terms.
        for (const auto &e : grp.entries) {
          new_or->addChild(
              AtomQuery::CHILD_TYPE(factorable[e.term_idx].original->copy()));
        }
      }
    }

    // Add factorable terms that don't contain the candidate.
    for (size_t i : terms_without) {
      new_or->addChild(AtomQuery::CHILD_TYPE(factorable[i].original->copy()));
    }

    // Add unfactorable terms (recursive, NOT, etc.).
    for (const auto *t : unfactorable) {
      new_or->addChild(AtomQuery::CHILD_TYPE(t->copy()));
    }

    // If the new OR has only one child, return the child directly.
    if (std::distance(new_or->beginChildren(), new_or->endChildren()) == 1) {
      auto *single = new_or->beginChildren()->get()->copy();
      delete new_or;
      if (verbose) {
        std::cout << "    result (single child, unwrapped): "
                  << query_to_smarts_str(single) << "\n";
      }
      return single;
    }

    if (verbose) {
      std::cout << "    result OR: " << query_to_smarts_str(new_or) << "\n";
    }
    return new_or;
  };

  // ---- Parse, transform, serialize ----
  std::unique_ptr<RDKit::ROMol> mol(RDKit::SmartsToMol(smarts));
  if (!mol) return smarts;

  bool changed = false;
  for (auto *atom : mol->atoms()) {
    auto *qatom = dynamic_cast<RDKit::QueryAtom *>(atom);
    if (!qatom || !qatom->getQuery()) continue;

    auto *q = qatom->getQuery();
    const std::string &desc = q->getDescription();

    if (verbose) {
      std::cout << "  atom " << atom->getIdx() << ": "
                << RDKit::SmartsWrite::GetAtomSmarts(
                       static_cast<const RDKit::Atom *>(qatom))
                << "  (top-level=" << desc << ")\n";
    }

    if (is_or(desc)) {
      AtomQuery *factored = try_factor(q);
      if (factored) {
        if (verbose)
          std::cout << "    -> factored top-level OR: "
                    << query_to_smarts_str(factored) << "\n";
        qatom->setQuery(factored);
        changed = true;
      } else {
        if (verbose) std::cout << "    -> top-level OR: no factoring\n";
      }
    } else if (is_and(desc)) {
      // Check whether any OR child can be factored.
      std::vector<std::pair<size_t, AtomQuery *>> factored_children;
      size_t idx = 0;
      for (auto it = q->beginChildren(); it != q->endChildren(); ++it, ++idx) {
        const auto *child = it->get();
        if (child && is_or(child->getDescription())) {
          if (verbose)
            std::cout << "    AND child " << idx << " is OR ("
                      << query_to_smarts_str(child) << "), trying factor\n";
          AtomQuery *f = try_factor(child);
          if (f) {
            if (verbose)
              std::cout << "    -> factored AND-child OR at idx " << idx << ": "
                        << query_to_smarts_str(f) << "\n";
            factored_children.push_back({idx, f});
          } else {
            if (verbose)
              std::cout << "    -> AND-child OR at idx " << idx
                        << ": no factoring\n";
          }
        }
      }

      if (!factored_children.empty()) {
        auto *new_and = new RDKit::ATOM_AND_QUERY;
        new_and->setDescription(desc);
        size_t child_idx = 0;
        size_t fc = 0;
        for (auto it = q->beginChildren(); it != q->endChildren();
             ++it, ++child_idx) {
          if (fc < factored_children.size() &&
              factored_children[fc].first == child_idx) {
            auto *factored = factored_children[fc].second;
            // Flatten an AND result into the parent AND.
            if (is_and(factored->getDescription())) {
              for (auto fit = factored->beginChildren();
                   fit != factored->endChildren(); ++fit) {
                new_and->addChild(AtomQuery::CHILD_TYPE(fit->get()->copy()));
              }
              delete factored;
            } else {
              new_and->addChild(AtomQuery::CHILD_TYPE(factored));
            }
            ++fc;
          } else {
            const auto *child = it->get();
            if (child) {
              new_and->addChild(AtomQuery::CHILD_TYPE(child->copy()));
            }
          }
        }
        qatom->setQuery(new_and);
        changed = true;
      }
    }
  }

  // ---- Pass 2: Factor within recursive expressions & unwrap ----
  // Helper: check if a query tree contains any OR node.
  std::function<bool(const AtomQuery *)> has_or_node;
  has_or_node = [&](const AtomQuery *q) -> bool {
    if (!q) return false;
    if (is_or(q->getDescription())) return true;
    for (auto it = q->beginChildren(); it != q->endChildren(); ++it) {
      if (has_or_node(it->get())) return true;
    }
    return false;
  };

  // Recursively walk a query tree, processing RecursiveStructureQuery nodes.
  // Returns a new query (caller owns) if any change was made, or nullptr.
  std::function<AtomQuery *(const AtomQuery *)> process_recursive_children;
  process_recursive_children = [&](const AtomQuery *q) -> AtomQuery * {
    if (!q) return nullptr;
    const std::string &d = q->getDescription();

    if (d == "RecursiveStructure") {
      const auto *rsq = dynamic_cast<const RDKit::RecursiveStructureQuery *>(q);
      if (!rsq || !rsq->getQueryMol()) return nullptr;

      // Serialize inner mol, factor it, check for changes.
      const std::string inner_smarts = RDKit::MolToSmarts(*rsq->getQueryMol());
      const std::string factored_inner =
          factor_or_common_primitives(inner_smarts, verbose,
                        user_factoring_priority);

      // Parse the factored inner SMARTS.
      std::unique_ptr<RDKit::ROMol> inner_mol(
          RDKit::SmartsToMol(factored_inner));
      if (!inner_mol) return nullptr;

      // Check if we can unwrap: 1-atom mol, no OR in the anchor query.
      if (inner_mol->getNumAtoms() == 1) {
        const auto *anchor = dynamic_cast<const RDKit::QueryAtom *>(
            inner_mol->getAtomWithIdx(0));
        if (anchor && anchor->getQuery() && !has_or_node(anchor->getQuery())) {
          if (verbose) {
            std::cout << "    unwrapping recursive $(" << factored_inner
                      << ") -> inline primitives\n";
          }
          // Return a copy of the anchor's query (unwrapped).
          auto *unwrapped = anchor->getQuery()->copy();
          // Preserve negation from the original RecursiveStructureQuery.
          if (q->getNegation()) {
            unwrapped->setNegation(!unwrapped->getNegation());
          }
          return unwrapped;
        }
      }

      // Inner mol changed but can't unwrap — rebuild RecursiveStructureQuery.
      if (factored_inner != inner_smarts) {
        auto *new_rsq = new RDKit::RecursiveStructureQuery(
            new RDKit::ROMol(*inner_mol, true), rsq->getSerialNumber());
        new_rsq->setNegation(q->getNegation());
        new_rsq->setDescription("RecursiveStructure");
        if (verbose) {
          std::cout << "    factored recursive inner: $(" << inner_smarts
                    << ") -> $(" << factored_inner << ")\n";
        }
        return new_rsq;
      }
      return nullptr;
    }

    // For AND/OR/NOT nodes, recurse into children.
    if (is_and(d) || is_or(d) || d == "AtomNot" || d == "Not") {
      std::vector<std::pair<size_t, AtomQuery *>> replacements;
      size_t idx = 0;
      for (auto it = q->beginChildren(); it != q->endChildren(); ++it, ++idx) {
        AtomQuery *result = process_recursive_children(it->get());
        if (result) {
          replacements.push_back({idx, result});
        }
      }
      if (replacements.empty()) return nullptr;

      // Rebuild this node with replaced children.
      AtomQuery *new_node = nullptr;
      if (is_and(d)) {
        auto *n = new RDKit::ATOM_AND_QUERY;
        n->setDescription(d);
        new_node = n;
      } else if (is_or(d)) {
        auto *n = new RDKit::ATOM_OR_QUERY;
        n->setDescription(d);
        new_node = n;
      } else {
        // NOT or other composite node — use AND wrapper to hold single child.
        auto *n = new RDKit::ATOM_AND_QUERY;
        n->setDescription(d);
        n->setNegation(q->getNegation());
        new_node = n;
      }

      size_t child_idx = 0;
      size_t ri = 0;
      for (auto it = q->beginChildren(); it != q->endChildren();
           ++it, ++child_idx) {
        if (ri < replacements.size() && replacements[ri].first == child_idx) {
          auto *repl = replacements[ri].second;
          // If replacing a RecursiveStructure with unwrapped AND query,
          // flatten it into the parent AND.
          if (is_and(d) && is_and(repl->getDescription()) &&
              !repl->getNegation()) {
            for (auto fit = repl->beginChildren(); fit != repl->endChildren();
                 ++fit) {
              new_node->addChild(AtomQuery::CHILD_TYPE(fit->get()->copy()));
            }
            delete repl;
          } else {
            new_node->addChild(AtomQuery::CHILD_TYPE(repl));
          }
          ++ri;
        } else {
          new_node->addChild(AtomQuery::CHILD_TYPE(it->get()->copy()));
        }
      }
      return new_node;
    }

    return nullptr;
  };

  if (recursive_pass_depth == 0) {
    ++recursive_pass_depth;
    for (auto *atom : mol->atoms()) {
      auto *qatom = dynamic_cast<RDKit::QueryAtom *>(atom);
      if (!qatom || !qatom->getQuery()) continue;

      AtomQuery *result = process_recursive_children(qatom->getQuery());
      if (result) {
        qatom->setQuery(result);
        changed = true;
        if (verbose) {
          std::cout << "  atom " << atom->getIdx()
                    << " after recursive processing: "
                    << RDKit::SmartsWrite::GetAtomSmarts(
                           static_cast<const RDKit::Atom *>(qatom))
                    << "\n";
        }
      }
    }
    --recursive_pass_depth;
  }

  if (!changed) {
    if (verbose)
      std::cout << "  result: no change -> " << smarts << "\n"
                << "===== end factor_or_common_primitives =====\n\n";
    return smarts;
  }
  auto result = RDKit::MolToSmarts(*mol);
  if (verbose)
    std::cout << "  result: " << result << "\n"
              << "===== end factor_or_common_primitives =====\n\n";
  return result;
}

std::string canonicalize_smarts_basic(const std::string &smarts) {
  std::unique_ptr<RDKit::ROMol> mol(RDKit::SmartsToMol(smarts));
  if (!mol) {
    return smarts;
  }
  return RDKit::MolToSmarts(*mol);
}

std::string clear_atom_maps_in_smarts(const std::string &smarts) {
  std::unique_ptr<RDKit::ROMol> mol(RDKit::SmartsToMol(smarts));
  if (!mol) {
    return smarts;
  }
  for (auto *atom : mol->atoms()) {
    if (atom) {
      atom->setAtomMapNum(0);
    }
  }
  return RDKit::MolToSmarts(*mol);
}

std::string canonicalize_negated_recursive_inner(const std::string &inner) {
  static thread_local int recursive_canon_depth = 0;

  if (recursive_canon_depth > 0) {
    return canonicalize_smarts_basic(inner);
  }
  atom_typer::SmartsAnalyzer::StandardSmartsWorkflowOptions workflow_options;
  workflow_options.include_x_in_reserialization = false;
  workflow_options.enumerate_bond_order = false;

  ++recursive_canon_depth;
  try {
    atom_typer::SmartsAnalyzer sa;

    const auto out =
        sa.standard_smarts({inner}, false, false, false, workflow_options);
    --recursive_canon_depth;
    if (!out.empty() && !out.front().empty()) {
      return clear_atom_maps_in_smarts(out.front());
    }
    return clear_atom_maps_in_smarts(canonicalize_smarts_basic(inner));
  } catch (...) {
    --recursive_canon_depth;
    return clear_atom_maps_in_smarts(canonicalize_smarts_basic(inner));
  }
}

std::string normalize_negated_recursive_queries(const std::string &smarts) {
  using AtomQuery = Queries::Query<int, const RDKit::Atom *, true>;

  // Fast path: no negated recursive query syntax present.
  if (smarts.find("!$(") == std::string::npos) {
    return smarts;
  }

  const auto normalize_negated_recursive_in_atom_token =
      [&](const std::string &atom_token) {
        std::string token = atom_token;
        if (token.size() < 2 || token.front() != '[' || token.back() != ']') {
          return token;
        }

        std::string inner = token.substr(1, token.size() - 2);
        size_t pos = 0;
        while (true) {
          const size_t qpos = inner.find("!$(", pos);
          if (qpos == std::string::npos) {
            break;
          }

          const size_t open_pos = qpos + 2;  // '(' in '!$('
          const size_t close_pos = find_matching_paren(inner, open_pos);
          if (close_pos == std::string::npos || close_pos <= open_pos + 1) {
            break;
          }

          const std::string recursive_smarts =
              inner.substr(open_pos + 1, close_pos - open_pos - 1);
          const std::string canonical_recursive =
              canonicalize_negated_recursive_inner(recursive_smarts);

          inner.replace(open_pos + 1, close_pos - open_pos - 1,
                        canonical_recursive);
          pos = open_pos + 1 + canonical_recursive.size();
        }

        return "[" + inner + "]";
      };

  const std::function<bool(const AtomQuery *, bool)> has_negated_recursive =
      [&](const AtomQuery *q, bool inherited_negation) {
        if (!q) {
          return false;
        }

        const bool effective_negation = inherited_negation ^ q->getNegation();
        const std::string desc = q->getDescription();
        if (desc == "RecursiveStructure" && effective_negation) {
          return true;
        }

        bool child_negation = effective_negation;
        if (desc == "AtomNot" || desc == "Not") {
          child_negation = !effective_negation;
        }

        for (auto it = q->beginChildren(); it != q->endChildren(); ++it) {
          if (has_negated_recursive(it->get(), child_negation)) {
            return true;
          }
        }
        return false;
      };

  std::unique_ptr<RDKit::ROMol> parsed(RDKit::SmartsToMol(smarts));
  if (!parsed) {
    return smarts;
  }

  RDKit::RWMol mol(*parsed);
  for (auto *atom : mol.atoms()) {
    auto *qatom = dynamic_cast<RDKit::QueryAtom *>(atom);
    if (!qatom || !qatom->getQuery()) {
      continue;
    }
    if (!has_negated_recursive(qatom->getQuery(), false)) {
      continue;
    }

    const std::string atom_token = RDKit::SmartsWrite::GetAtomSmarts(qatom);
    const std::string rewritten_token =
        normalize_negated_recursive_in_atom_token(atom_token);
    if (rewritten_token == atom_token) {
      continue;
    }

    std::unique_ptr<RDKit::ROMol> one_atom(RDKit::SmartsToMol(rewritten_token));
    if (!one_atom || one_atom->getNumAtoms() != 1) {
      continue;
    }

    auto *src = dynamic_cast<RDKit::QueryAtom *>(one_atom->getAtomWithIdx(0));
    if (!src || !src->getQuery()) {
      continue;
    }
    qatom->setQuery(src->getQuery()->copy());
  }

  return RDKit::MolToSmarts(mol);
}

std::string move_explicit_hydrogens_to_not_h0(const std::string &smarts) {
  const auto likely_has_explicit_h_atom = [](const std::string &s) {
    for (size_t i = 0; i < s.size(); ++i) {
      if (s[i] != '[' || i + 1 >= s.size()) {
        continue;
      }

      const size_t j = i + 1;
      if (s[j] == 'H') {
        if (j + 1 >= s.size()) {
          return true;
        }
        const char next = s[j + 1];
        // Exclude element symbols like Hg.
        if (!std::islower(static_cast<unsigned char>(next))) {
          return true;
        }
      }

      if (s[j] == '#' && j + 1 < s.size() && s[j + 1] == '1') {
        const size_t k = j + 2;
        if (k >= s.size() || !std::isdigit(static_cast<unsigned char>(s[k]))) {
          return true;
        }
      }
    }
    return false;
  };

  // Fast path: no explicit hydrogen atom tokens to fold.
  if (!likely_has_explicit_h_atom(smarts)) {
    return smarts;
  }

  std::unique_ptr<RDKit::ROMol> parsed(RDKit::SmartsToMol(smarts));
  if (!parsed) {
    return smarts;
  }

  RDKit::RWMol mol(*parsed);
  std::set<unsigned int> hydrogen_atom_indices;
  std::map<unsigned int, unsigned int> explicit_hydrogen_count;

  for (const auto *atom : mol.atoms()) {
    if (!is_explicit_hydrogen_atom(atom)) {
      continue;
    }
    if (atom->getDegree() != 1) {
      continue;
    }

    const auto *nbr = *(mol.atomNeighbors(atom).begin());
    if (!nbr) {
      continue;
    }

    hydrogen_atom_indices.insert(atom->getIdx());
    ++explicit_hydrogen_count[nbr->getIdx()];
  }

  for (const auto &entry : explicit_hydrogen_count) {
    const auto atom_idx = entry.first;
    const auto num_explicit_h = entry.second;
    if (num_explicit_h == 0) {
      continue;
    }

    auto *dst = dynamic_cast<RDKit::QueryAtom *>(mol.getAtomWithIdx(atom_idx));
    if (!dst) {
      continue;
    }

    std::string inner =
        strip_brackets_and_map(RDKit::SmartsWrite::GetAtomSmarts(dst));
    inner = set_or_append_lowercase_h_count(inner, num_explicit_h);

    if (inner.empty()) {
      inner = "*";
    }
    const std::string rewritten = "[" + inner + "]";
    std::unique_ptr<RDKit::ROMol> probe(RDKit::SmartsToMol(rewritten));
    if (!probe || probe->getNumAtoms() != 1) {
      continue;
    }
    auto *src = dynamic_cast<RDKit::QueryAtom *>(probe->getAtomWithIdx(0));
    if (!src || !src->getQuery()) {
      continue;
    }
    dst->setQuery(src->getQuery()->copy());
  }

  std::vector<unsigned int> remove_indices(hydrogen_atom_indices.begin(),
                                           hydrogen_atom_indices.end());
  std::sort(remove_indices.begin(), remove_indices.end(), std::greater<>());
  for (const auto idx : remove_indices) {
    if (idx < mol.getNumAtoms()) {
      mol.removeAtom(idx);
    }
  }

  return RDKit::MolToSmarts(mol);
}

std::vector<std::string> perceive_tautomers(
    const std::vector<std::string> &variants, bool verbose = false) {
  struct TautomerRule {
    std::string before;
    std::string after;
  };

  // Rule list can be extended as needed.
  // Current focus: nitro-like valence representations.
  static const std::vector<TautomerRule> rules = {
      {"[#7+:1](=[#8:2])=[#8:3]", "[#7+:1](=[#8:2])-[#8-:3]"},
      {"[#7:1](=[#8:2])=[#8:3]", "[#7+:1](=[#8:2])-[#8-:3]"},
  };

  std::vector<std::string> out;
  std::set<std::string> seen;
  out.reserve(variants.size());

  for (const auto &v : variants) {
    std::string standardized = v;

    std::unique_ptr<RDKit::ROMol> current(RDKit::SmartsToMol(standardized));
    if (current) {
      for (const auto &rule : rules) {
        std::unique_ptr<RDKit::ROMol> before_q(RDKit::SmartsToMol(rule.before));
        std::unique_ptr<RDKit::ROMol> after_q(RDKit::SmartsToMol(rule.after));
        if (!before_q || !after_q) {
          continue;
        }

        std::vector<RDKit::MatchVectType> matches;
        if (!RDKit::SubstructMatch(*current, *before_q, matches) ||
            matches.empty()) {
          continue;
        }

        const auto replaced = RDKit::replaceSubstructs(
            *current, *before_q, *after_q, false, 0, false);
        if (!replaced.empty() && replaced.front()) {
          current = std::make_unique<RDKit::ROMol>(*replaced.front());
          standardized = RDKit::MolToSmarts(*current);
          if (verbose) {
            std::cout << "[perceive_tautomers] applied rule '" << rule.before
                      << "' -> '" << rule.after << "' on '" << v << "' => '"
                      << standardized << "'" << std::endl;
          }
        }
      }
    }

    if (seen.insert(standardized).second) {
      out.push_back(std::move(standardized));
    }
  }

  return out;
}

// ---------------------------------------------------------------------------
// normalize_smarts_encoding
//
// Final encoding normalization pass, controlled by StandardSmartsWorkflowOptions
// flags applied after factoring and collapsing.
//
//  remove_aa_wildcard    – Remove A,a OR sub-trees (aromaticity tautology).
//  symbol_form           – Expanded: C->[#6&A], c->[#6&a];
//                          Condensed: [#6&A]->C, [#6&a]->c.
//  fold_singleton_or     – Collapse OR(X) -> X.
//  explicit_charge_values– + -> +1, - -> -1, ++ -> +2, -- -> -2.
//  rewrite_or_primitives_to_negated + or_primitive_rewrite_atomic_nums:
//                         per-atom rewrites like H1,H2,H3,H4 -> !H0.
// ---------------------------------------------------------------------------
std::string normalize_smarts_encoding(
    const std::string &smarts,
    const SmartsAnalyzer::StandardSmartsWorkflowOptions &opts) {
  using AtomQuery = Queries::Query<int, const RDKit::Atom *, true>;

  const bool do_tree = opts.remove_aa_wildcard || opts.fold_singleton_or ||
                       opts.rewrite_or_primitives_to_negated;
  const bool do_expand =
      opts.symbol_form == SmartsAnalyzer::SymbolForm::Expanded;
  const bool do_condense =
      opts.symbol_form == SmartsAnalyzer::SymbolForm::Condensed;
  const bool do_charge = opts.explicit_charge_values;

  if (!do_tree && !do_expand && !do_condense && !do_charge) return smarts;

  // Organic-subset lookup tables (2-char entries must precede 1-char).
  static const std::vector<std::pair<std::string, int>> kAliphatic = {
      {"Cl", 17}, {"Br", 35}, {"B", 5},  {"C", 6},  {"N", 7},
      {"O", 8},   {"F", 9},   {"P", 15}, {"S", 16}, {"I", 53}};
  static const std::vector<std::pair<std::string, int>> kAromatic = {
      {"b", 5}, {"c", 6}, {"n", 7}, {"o", 8}, {"p", 15}, {"s", 16}};
  // Z -> (uppercase symbol, lowercase symbol or "").
  static const std::map<int, std::pair<std::string, std::string>> kZ2Sym = {
      {5, {"B", "b"}},   {6, {"C", "c"}},  {7, {"N", "n"}},
      {8, {"O", "o"}},   {9, {"F", ""}},   {15, {"P", "p"}},
      {16, {"S", "s"}},  {17, {"Cl", ""}}, {35, {"Br", ""}},
      {53, {"I", ""}}};

  const auto is_and = [](const std::string &d) {
    return d == "AtomAnd" || d == "And";
  };
  const auto is_or = [](const std::string &d) {
    return d == "AtomOr" || d == "Or";
  };

  const auto is_target_atomic_num = [&](int z) {
    if (!opts.rewrite_or_primitives_to_negated) return false;
    return std::find(opts.or_primitive_rewrite_atomic_nums.begin(),
                     opts.or_primitive_rewrite_atomic_nums.end(),
                     z) != opts.or_primitive_rewrite_atomic_nums.end();
  };

  const auto eq_val_if_desc = [](const AtomQuery *q,
                                 const std::set<std::string> &descs,
                                 int *out_val = nullptr) -> bool {
    const auto *eq = dynamic_cast<const RDKit::ATOM_EQUALS_QUERY *>(q);
    if (!eq) return false;
    if (!descs.count(q->getDescription())) return false;
    if (q->getNegation()) return false;
    if (out_val) *out_val = eq->getVal();
    return true;
  };

  std::function<int(const AtomQuery *)> infer_atomic_num;
  infer_atomic_num = [&](const AtomQuery *q) -> int {
    if (!q) return -1;
    int z = -1;
    if (eq_val_if_desc(q, {"AtomicNum", "AtomAtomicNum"}, &z)) return z;
    if (is_and(q->getDescription())) {
      int found = -1;
      for (auto it = q->beginChildren(); it != q->endChildren(); ++it) {
        int child_z = infer_atomic_num(it->get());
        if (child_z < 0) continue;
        if (found < 0) {
          found = child_z;
        } else if (found != child_z) {
          return -1;
        }
      }
      return found;
    }
    return -1;
  };

  std::function<bool(const AtomQuery *, std::set<int> &)> collect_h_values_from_honly_query;
  collect_h_values_from_honly_query = [&](const AtomQuery *q,
                                          std::set<int> &vals) -> bool {
    if (!q) return false;
    const std::string d = q->getDescription();
    if (is_or(d)) {
      bool had_child = false;
      for (auto it = q->beginChildren(); it != q->endChildren(); ++it) {
        had_child = true;
        if (!collect_h_values_from_honly_query(it->get(), vals)) return false;
      }
      return had_child;
    }
    int h = -1;
    if (eq_val_if_desc(q, {"AtomHCount", "AtomTotalHCount", "AtomImplicitHCount"}, &h)) {
      vals.insert(h);
      return true;
    }
    return false;
  };

  const auto is_allowed_extra_for_h4_carbon = [&](const AtomQuery *q) {
    int v = -1;
    if (eq_val_if_desc(q, {"AtomExplicitDegree", "AtomTotalDegree", "AtomDegree"}, &v)) {
      return v == 0;
    }
    if (eq_val_if_desc(q, {"AtomFormalCharge"}, &v)) {
      return v == 0;
    }
    return false;
  };

  std::function<bool(const AtomQuery *, int, std::set<int> &)> collect_h_values_from_arm;
  collect_h_values_from_arm = [&](const AtomQuery *arm,
                                  int atom_z,
                                  std::set<int> &vals) -> bool {
    if (!arm) return false;

    int h = -1;
    if (eq_val_if_desc(arm,
                       {"AtomHCount", "AtomTotalHCount", "AtomImplicitHCount"},
                       &h)) {
      vals.insert(h);
      return true;
    }

    if (is_and(arm->getDescription())) {
      int hcount_seen = -1;
      std::vector<const AtomQuery *> extras;
      for (auto it = arm->beginChildren(); it != arm->endChildren(); ++it) {
        const auto *child = it->get();
        int hv = -1;
        if (eq_val_if_desc(child,
                           {"AtomHCount", "AtomTotalHCount", "AtomImplicitHCount"},
                           &hv)) {
          if (hcount_seen >= 0 && hcount_seen != hv) return false;
          hcount_seen = hv;
          continue;
        }
        extras.push_back(child);
      }
      if (hcount_seen < 0) return false;
      for (const auto *extra : extras) {
        if (!(atom_z == 6 && hcount_seen == 4 &&
              is_allowed_extra_for_h4_carbon(extra))) {
          return false;
        }
      }
      vals.insert(hcount_seen);
      return true;
    }

    if (arm->getDescription() == "RecursiveStructure" && !arm->getNegation()) {
      const auto *rsq =
          dynamic_cast<const RDKit::RecursiveStructureQuery *>(arm);
      if (!rsq || !rsq->getQueryMol()) return false;
      const auto *qm = rsq->getQueryMol();
      if (qm->getNumAtoms() != 1) return false;
      const auto *a0 = qm->getAtomWithIdx(0);
      const auto *qa0 = dynamic_cast<const RDKit::QueryAtom *>(a0);
      if (!qa0 || !qa0->getQuery()) return false;
      std::set<int> inner_vals;
      if (!collect_h_values_from_honly_query(qa0->getQuery(), inner_vals)) return false;
      vals.insert(inner_vals.begin(), inner_vals.end());
      return true;
    }

    return false;
  };

  // Returns the aromaticity value of a sole-aromaticity leaf:
  //   1=aromatic, 0=aliphatic, -1=not a sole aromaticity leaf.
  const auto sole_aromaticity_value = [](const AtomQuery *q) -> int {
    const auto *eq = dynamic_cast<const RDKit::ATOM_EQUALS_QUERY *>(q);
    if (!eq) return -1;
    const std::string &d = eq->getDescription();
    // "a" in SMARTS → AtomIsAromatic (val=1)
    if (d == "AtomAromatic" || d == "AtomIsAromatic") {
      const bool is_arom = (eq->getVal() != 0) ^ eq->getNegation();
      return is_arom ? 1 : 0;
    }
    // "A" in SMARTS → AtomIsAliphatic (val=1)
    if (d == "AtomAliphatic" || d == "AtomIsAliphatic") {
      const bool is_ali = (eq->getVal() != 0) ^ eq->getNegation();
      return is_ali ? 0 : 1;
    }
    return -1;
  };

  // Flatten nested OR into a flat list of arm pointers (no copies).
  std::function<void(const AtomQuery *, std::vector<const AtomQuery *> &)>
      flatten_or;
  flatten_or = [&](const AtomQuery *q,
                   std::vector<const AtomQuery *> &out) {
    if (!q) return;
    if (is_or(q->getDescription())) {
      for (auto it = q->beginChildren(); it != q->endChildren(); ++it)
        flatten_or(it->get(), out);
    } else {
      out.push_back(q);
    }
  };

  // Sentinel pointer meaning "this node is always-true; drop from AND parent."
  AtomQuery *const kTautology =
      reinterpret_cast<AtomQuery *>(static_cast<uintptr_t>(1));

  // Recursive query-tree rewriter.
  // Returns: nullptr = no change, kTautology = always-true (remove from AND),
  //          otherwise a new owned node to substitute.
  std::function<AtomQuery *(const AtomQuery *)> rewrite;
  rewrite = [&](const AtomQuery *q) -> AtomQuery * {
    if (!q) return nullptr;
    const std::string &d = q->getDescription();

    if (is_or(d)) {
      std::vector<const AtomQuery *> flat;
      flatten_or(q, flat);

      // remove_aa_wildcard: if both a sole "aromatic" arm and a sole
      //  "aliphatic" arm exist, the whole OR is a tautology.
      if (opts.remove_aa_wildcard) {
        bool has_arom = false, has_ali = false;
        for (const auto *arm : flat) {
          const int av = sole_aromaticity_value(arm);
          if (av == 1) has_arom = true;
          else if (av == 0) has_ali  = true;
        }
        if (has_arom && has_ali) return kTautology;
      }

      // fold_singleton_or: OR with one arm -> arm directly.
      if (opts.fold_singleton_or && flat.size() == 1)
        return flat[0]->copy();

      // Recursively rewrite children.
      bool changed = false;
      std::vector<AtomQuery *> new_arms;
      new_arms.reserve(flat.size());
      for (const auto *arm : flat) {
        AtomQuery *rw = rewrite(arm);
        if (rw == kTautology) {
          // A tautology arm makes the whole OR a tautology.
          for (auto *tmp : new_arms) delete tmp;
          return kTautology;
        }
        if (rw) {
          new_arms.push_back(rw);
          changed = true;
        } else {
          new_arms.push_back(arm->copy());
        }
      }

      if (!changed) {
        for (auto *tmp : new_arms) delete tmp;
        return nullptr;
      }
      if (new_arms.size() == 1) return new_arms[0];

      auto *new_or = new RDKit::ATOM_OR_QUERY;
      new_or->setDescription(d);
      new_or->setNegation(q->getNegation());
      for (auto *arm : new_arms)
        new_or->addChild(AtomQuery::CHILD_TYPE(arm));
      return new_or;

    } else if (is_and(d)) {
      bool changed = false;
      std::vector<AtomQuery *> new_children;
      for (auto it = q->beginChildren(); it != q->endChildren(); ++it) {
        AtomQuery *rw = rewrite(it->get());
        if (rw == kTautology) {
          // Drop always-true child from AND.
          changed = true;
          continue;
        }
        if (rw) {
          new_children.push_back(rw);
          changed = true;
        } else {
          new_children.push_back(it->get()->copy());
        }
      }

      // Normalize explicit degree for neutral hydrogen-bearing halogens.
      // RDKit's explicit-degree semantics count explicit H isotopes as
      // neighbors, so H1-bearing halogen terms should use D1 rather than D0.
      const int atom_z = infer_atomic_num(q);
      if (atom_z == 17 || atom_z == 35 || atom_z == 53) {
        bool has_h1 = false;
        ssize_t degree_idx = -1;
        int degree_val = -1;
        for (size_t i = 0; i < new_children.size(); ++i) {
          int val = -1;
          if (eq_val_if_desc(new_children[i],
                             {"AtomHCount", "AtomTotalHCount", "AtomImplicitHCount"},
                             &val) &&
              val == 1) {
            has_h1 = true;
          }
          if (eq_val_if_desc(new_children[i],
                             {"AtomExplicitDegree", "AtomDegree", "AtomTotalDegree"},
                             &val)) {
            degree_idx = static_cast<ssize_t>(i);
            degree_val = val;
          }
        }

        if (has_h1 && degree_idx >= 0 && degree_val == 0) {
          auto *deg_eq =
              dynamic_cast<RDKit::ATOM_EQUALS_QUERY *>(new_children[degree_idx]);
          if (deg_eq) {
            deg_eq->setVal(1);
            changed = true;
          }
        }
      }

      if (opts.rewrite_or_primitives_to_negated) {
        if (is_target_atomic_num(atom_z)) {
          for (size_t i = 0; i < new_children.size(); ++i) {
            AtomQuery *child = new_children[i];
            if (!child || !is_or(child->getDescription())) continue;
            std::vector<const AtomQuery *> flat_arms;
            flatten_or(child, flat_arms);
            if (flat_arms.empty()) continue;

            std::set<int> h_vals;
            bool all_convertible = true;
            for (const auto *arm : flat_arms) {
              if (!collect_h_values_from_arm(arm, atom_z, h_vals)) {
                all_convertible = false;
                break;
              }
            }

            if (all_convertible && h_vals == std::set<int>({1, 2, 3, 4})) {
              auto *not_h0 = RDKit::makeAtomHCountQuery(0);
              not_h0->setNegation(true);
              delete new_children[i];
              new_children[i] = not_h0;
              changed = true;
            }
          }
        }
      }

      if (!changed) {
        for (auto *c : new_children) delete c;
        return nullptr;
      }
      if (new_children.empty()) return kTautology;
      if (new_children.size() == 1) return new_children[0];

      auto *new_and = new RDKit::ATOM_AND_QUERY;
      new_and->setDescription(d);
      new_and->setNegation(q->getNegation());
      for (auto *c : new_children)
        new_and->addChild(AtomQuery::CHILD_TYPE(c));
      return new_and;

    } else {
      return nullptr;  // leaf — no structural change needed
    }
  };

  std::string result = smarts;

  // ── Query-tree pass ──────────────────────────────────────────────────────
  if (do_tree) {
    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmartsToMol(result));
    if (mol) {
      bool changed = false;
      for (auto *atom : mol->atoms()) {
        auto *qa = dynamic_cast<RDKit::QueryAtom *>(atom);
        if (!qa || !qa->getQuery()) continue;
        AtomQuery *rw = rewrite(qa->getQuery());
        if (rw == kTautology) {
          // Replace with a null (match-anything) query.
          qa->setQuery(RDKit::makeAtomNullQuery());
          changed = true;
        } else if (rw) {
          qa->setQuery(rw);
          changed = true;
        }
      }
      if (changed) result = RDKit::MolToSmarts(*mol);
    }
  }

  // ── Optional primitive rewrite pass (token-level) ───────────────────────
  // Handles selected per-atom OR-primitive encodings that are equivalent to
  // a negated primitive, e.g. carbon H1,H2,H3,H4 -> !H0.
  if (opts.rewrite_or_primitives_to_negated && is_target_atomic_num(6)) {
    std::string rewritten;
    rewritten.reserve(result.size());
    size_t i = 0;
    while (i < result.size()) {
      if (result[i] != '[') {
        rewritten += result[i++];
        continue;
      }
      size_t j = i;
      int depth = 0;
      for (; j < result.size(); ++j) {
        if (result[j] == '[') ++depth;
        else if (result[j] == ']') {
          --depth;
          if (!depth) break;
        }
      }
      if (j >= result.size()) {
        rewritten += result.substr(i);
        break;
      }

      std::string inner = result.substr(i + 1, j - i - 1);
      const bool atom6 =
          inner.rfind("C", 0) == 0 ||
          inner.rfind("#6", 0) == 0;
      if (atom6) {
        auto replace_all = [](std::string &s,
                              const std::string &from,
                              const std::string &to) {
          size_t p = 0;
          while ((p = s.find(from, p)) != std::string::npos) {
            s.replace(p, from.size(), to);
            p += to.size();
          }
        };
        replace_all(inner, "$([H1,H2,H3]),H4&D0&+0", "!H0");
        replace_all(inner, "$([H1,H2,H3]),H4&+0", "!H0");
        replace_all(inner, "H1,H2,H3,H4", "!H0");
      }

      rewritten += "[" + inner + "]";
      i = j + 1;
    }
    result = std::move(rewritten);

  }


  // Helper: is c a SMARTS atom-expression delimiter?
  const auto is_delim = [](char c) {
    return c == '&' || c == ';' || c == ',' || c == '!';
  };

  // ── Expand element symbols: C->[#6&A], c->[#6&a] ────────────────────────
  if (do_expand) {
    std::map<std::string, int> ali_sym_to_z, aro_sym_to_z;
    for (const auto &[sym, z] : kAliphatic) ali_sym_to_z[sym] = z;
    for (const auto &[sym, z] : kAromatic)  aro_sym_to_z[sym] = z;

    std::string expanded;
    expanded.reserve(result.size() * 2);
    size_t i = 0;
    while (i < result.size()) {
      if (result[i] == '[') {
        // Bracket atom: find closing ']'.
        size_t j = i;
        int depth = 0;
        for (; j < result.size(); ++j) {
          if (result[j] == '[') ++depth;
          else if (result[j] == ']') { --depth; if (!depth) break; }
        }
        const std::string tok = result.substr(i, j - i + 1);
        std::string inner = tok.substr(1, tok.size() - 2);

        // Walk inner, expanding element symbols not already in #N form.
        std::string new_inner;
        new_inner.reserve(inner.size());
        size_t k = 0;
        while (k < inner.size()) {
          // Pass '#N' sequences through unchanged.
          if (inner[k] == '#') {
            new_inner += inner[k++];
            while (k < inner.size() &&
                   std::isdigit((unsigned char)inner[k]))
              new_inner += inner[k++];
            continue;
          }
          // At a position following a delimiter (or the start), check for an
          // organic-subset element symbol.
          const bool at_el_pos =
              k == 0 || is_delim(inner[k - 1]);
          if (at_el_pos) {
            // Try 2-char aliphatic (Cl, Br) first.
            bool matched = false;
            if (k + 1 < inner.size() &&
                std::islower((unsigned char)inner[k + 1])) {
              const std::string two = inner.substr(k, 2);
              auto it = ali_sym_to_z.find(two);
              if (it != ali_sym_to_z.end()) {
                new_inner += "#" + std::to_string(it->second) + "&A";
                k += 2;
                matched = true;
              }
            }
            if (!matched && std::isupper((unsigned char)inner[k])) {
              const std::string one(1, inner[k]);
              auto it = ali_sym_to_z.find(one);
              if (it != ali_sym_to_z.end()) {
                new_inner += "#" + std::to_string(it->second) + "&A";
                ++k;
                matched = true;
              }
            }
            // Lowercase aromatic symbols (b,c,n,o,p,s) — skip 'a' (aromaticity
            // primitive) and 'r' (ring primitive).
            if (!matched && std::islower((unsigned char)inner[k]) &&
                inner[k] != 'a' && inner[k] != 'r') {
              const std::string one(1, inner[k]);
              auto it = aro_sym_to_z.find(one);
              if (it != aro_sym_to_z.end()) {
                new_inner += "#" + std::to_string(it->second) + "&a";
                ++k;
                matched = true;
              }
            }
            if (!matched) new_inner += inner[k++];
          } else {
            new_inner += inner[k++];
          }
        }
        expanded += "[" + new_inner + "]";
        i = j + 1;
      } else {
        // Non-bracket context: expand standalone organic-subset atoms.
        const char c = result[i];
        if (c == '(' || c == ')' || c == '.' || c == '-' || c == '=' ||
            c == '~' || c == '@' || c == '/' || c == '\\' || c == '%' ||
            c == ':' || std::isdigit((unsigned char)c)) {
          expanded += result[i++];
          continue;
        }
        // '#' outside brackets is ring-bond order — pass through.
        if (c == '#' && (i == 0 || result[i - 1] != '[')) {
          expanded += result[i++];
          continue;
        }
        // 2-char aliphatic (Cl, Br).
        if (i + 1 < result.size() &&
            std::islower((unsigned char)result[i + 1])) {
          if (result.substr(i, 2) == "Cl" || result.substr(i, 2) == "Br") {
            const int z = (result[i + 1] == 'l') ? 17 : 35;
            expanded += "[#" + std::to_string(z) + "&A]";
            i += 2;
            continue;
          }
        }
        // 1-char aliphatic uppercase.
        if (std::isupper((unsigned char)c)) {
          const std::string one(1, c);
          auto it = ali_sym_to_z.find(one);
          if (it != ali_sym_to_z.end()) {
            expanded += "[#" + std::to_string(it->second) + "&A]";
            ++i;
            continue;
          }
        }
        // 1-char aromatic lowercase.
        if (std::islower((unsigned char)c) && c != 'a' && c != 'r') {
          const std::string one(1, c);
          auto it = aro_sym_to_z.find(one);
          if (it != aro_sym_to_z.end()) {
            expanded += "[#" + std::to_string(it->second) + "&a]";
            ++i;
            continue;
          }
        }
        expanded += result[i++];
      }
    }
    result = std::move(expanded);
  }

  // ── Condense symbols: [#6&A]->C, [#6&a]->c ──────────────────────────────
  // Walk bracket tokens and replace (#Z, A/a) primitive pairs with the
  // organic-subset element symbol wherever the condensation is unambiguous.
  if (do_condense) {
    std::string condensed;
    condensed.reserve(result.size());
    size_t i = 0;
    while (i < result.size()) {
      if (result[i] != '[') {
        condensed += result[i++];
        continue;
      }
      size_t j = i;
      int depth = 0;
      for (; j < result.size(); ++j) {
        if (result[j] == '[') ++depth;
        else if (result[j] == ']') { --depth; if (!depth) break; }
      }
      const std::string tok = result.substr(i, j - i + 1);
      const std::string inner = tok.substr(1, tok.size() - 2);

      // Scan inner for a #Z token and a standalone A/a primitive.
      size_t z_val = 0;
      bool z_found = false;
      int arom_val = -1;  // 0=aliphatic, 1=aromatic
      size_t z_start = std::string::npos, z_end = 0;
      size_t a_start = std::string::npos;

      for (size_t p = 0; p < inner.size(); ) {
        if (inner[p] == '#') {
          const size_t ns = p + 1;
          size_t ne = ns;
          while (ne < inner.size() && std::isdigit((unsigned char)inner[ne]))
            ++ne;
          if (ne > ns && !z_found) {
            z_val = std::stoul(inner.substr(ns, ne - ns));
            z_found = true;
            z_start = p;
            z_end   = ne;
          }
          p = ne;
          continue;
        }
        // Standalone A or a aromaticity primitive.
        if ((inner[p] == 'A' || inner[p] == 'a') &&
            arom_val < 0 &&
            (p == 0 || is_delim(inner[p - 1])) &&
            (p + 1 >= inner.size() || is_delim(inner[p + 1]) ||
             inner[p + 1] == ':')) {
          // Ensure it's not the start of a multi-char token.
          if (p + 1 >= inner.size() ||
              !std::isalpha((unsigned char)inner[p + 1])) {
            arom_val = (inner[p] == 'a') ? 1 : 0;
            a_start = p;
          }
        }
        ++p;
      }

      if (z_found && arom_val >= 0 &&
          kZ2Sym.count(static_cast<int>(z_val))) {
        const auto &syms = kZ2Sym.at(static_cast<int>(z_val));
        const std::string &sym =
            (arom_val == 1) ? syms.second : syms.first;
        if (!sym.empty()) {
          // Build new inner by deleting the #Z and A/a tokens (plus any
          // immediately preceding delimiter), then prepend the symbol.
          std::string ni = inner;
          // Collect the two spans to erase; erase back-to-front.
          std::vector<std::pair<size_t, size_t>> spans;
          const auto add_span = [&](size_t start, size_t end_excl) {
            size_t s = start, e = end_excl;
            if (s > 0 && is_delim(ni[s - 1])) --s;  // eat preceding delim
            spans.push_back({s, e - s});
          };
          add_span(z_start, z_end);
          add_span(a_start, a_start + 1);
          std::sort(spans.begin(), spans.end(),
                    [](const auto &x, const auto &y) {
                      return x.first > y.first;
                    });
          for (const auto &[s, l] : spans) {
            if (s < ni.size())
              ni.erase(s, std::min(l, ni.size() - s));
          }
          // Strip any leading delimiter left over.
          while (!ni.empty() && is_delim(ni[0])) ni.erase(0, 1);
          // Re-attach symbol.
          ni = sym + (ni.empty() ? "" : "&" + ni);
          condensed += "[" + ni + "]";
          i = j + 1;
          continue;
        }
      }
      condensed += tok;
      i = j + 1;
    }
    result = std::move(condensed);
  }

  // ── Explicit charge values: + -> +1, - -> -1, ++ -> +2, -- -> -2 ────────
  if (do_charge) {
    std::string charged;
    charged.reserve(result.size());
    size_t i = 0;
    while (i < result.size()) {
      if (result[i] != '[') {
        charged += result[i++];
        continue;
      }
      size_t j = i;
      int depth = 0;
      for (; j < result.size(); ++j) {
        if (result[j] == '[') ++depth;
        else if (result[j] == ']') { --depth; if (!depth) break; }
      }
      const std::string tok = result.substr(i, j - i + 1);
      const std::string inner = tok.substr(1, tok.size() - 2);
      std::string new_inner;
      new_inner.reserve(inner.size() + 4);
      size_t k = 0;
      while (k < inner.size()) {
        if (inner[k] == '+' || inner[k] == '-') {
          const char sign = inner[k];
          size_t run = 0;
          while (k + run < inner.size() && inner[k + run] == sign) ++run;
          const bool has_digit =
              (k + run < inner.size() &&
               std::isdigit((unsigned char)inner[k + run]));
          if (has_digit) {
            // Already explicit: copy the sign(s) and digit(s) as-is.
            for (size_t r = 0; r < run; ++r) new_inner += sign;
            k += run;
          } else {
            // Abbreviated: +/++ → +1/+2, -/-- → -1/-2.
            new_inner += sign;
            if (run > 1) new_inner += std::to_string(static_cast<int>(run));
            else         new_inner += '1';
            k += run;
          }
        } else {
          new_inner += inner[k++];
        }
      }
      charged += "[" + new_inner + "]";
      i = j + 1;
    }
    result = std::move(charged);
  }

  return result;
}

}  // namespace

std::vector<std::vector<std::string>> SmartsAnalyzer::generate_all_combos(
    std::vector<std::string> smarts_list, bool verbose, bool ignoreValence,
    bool catchErrors, const StandardSmartsWorkflowOptions &workflow_options,
    const StandardSmartsLogOptions &log_options) {
  ZoneScopedN("SmartsAnalyzer::generate_all_combos");
  atom_typer::SmartsAnalyzer sa;
  atom_typer::AtomTyper at;

  const auto log_enabled = [&](unsigned int flag) {
    if (verbose) {
      return true;
    }
    return log_options.enabled && (log_options.flags & flag) != 0u;
  };

  at.set_debug_level(log_enabled(LogRecanon) ? DebugLevel::Trace
                                             : DebugLevel::Off);
  query_reorder::set_comparison_trace_enabled(
      log_enabled(LogRecanonComparisons));

  std::vector<std::vector<std::string>> results;
  std::unordered_map<std::string, std::string> h_merge_cache;
  std::unordered_map<std::string, std::string> normalized_neg_recursive_cache;
  std::unordered_map<std::string, bool> valence_cache;
  for (const auto &test1 : smarts_list) {
    std::vector<std::string> these_results;
    // std::cout << "Processing SMARTS: " << test1 << std::endl;
    const std::string mapped_smarts = sa.add_atom_maps(test1);
    int max_amap = 0;
    auto mol = RDKit::SmartsToMol(mapped_smarts, 0, true);
    for (const auto *atom : mol->atoms()) {
      int amap = atom->getAtomMapNum();
      if (amap > max_amap) {
        max_amap = amap;
      }
    }
    if (log_enabled(LogMapping)) {
      std::cout << "Mapped SMARTS: " << mapped_smarts << std::endl;
    }

    const bool enumerate_verbose = log_enabled(LogVariants);
    const auto variants =
        sa.enumerate_variants(mapped_smarts, 1000, enumerate_verbose, true,
                              workflow_options.enumerate_bond_order);

    if (log_enabled(LogSummary)) {
      std::cout << "Generated " << variants.size() << " variants in r1.\n";
    }
    // const auto def = at.get_default_query_embedding();
    int intial_max_amap = max_amap;
    for (const auto &variant : variants) {
      max_amap = intial_max_amap;  // reset max_amap for next iteration
      if (log_enabled(LogVariants)) {
        std::cout << ">>> Processing variant: " << variant << std::endl;
      }
      std::string h_merged_smarts;
      {
        auto hm_it = h_merge_cache.find(variant);
        if (hm_it != h_merge_cache.end()) {
          h_merged_smarts = hm_it->second;
        } else {
          h_merged_smarts = move_explicit_hydrogens_to_not_h0(variant);
          h_merge_cache.emplace(variant, h_merged_smarts);
        }
      }
      // if (verbose) {
      //   std::cout << "Pre H-merge: " << variant << std::endl;
      //   std::cout << "Post H-merge: " << h_merged_smarts << std::endl;
      // }
      std::string extracted_smarts;
      try {
        extracted_smarts = at.type_atoms_from_smarts(
            h_merged_smarts, false, max_amap, log_enabled(LogAtomTyping),
            workflow_options.include_x_in_reserialization,
          workflow_options.enumerate_bond_order,
          workflow_options.extracted_primitives_mask);
      } catch (const std::exception &e) {
        if (log_enabled(LogErrors)) {
          std::cout << "Error processing variant '" << variant
                    << "': " << e.what() << std::endl;
        }
        if (catchErrors) {
          throw std::runtime_error("Error processing variant: " + variant +
                                   " - " + std::string(e.what()));
        } else {
          extracted_smarts = h_merged_smarts;  // Fallback to H-merged SMARTS if
                                               // processing fails.
        }
      }
      if (log_enabled(LogAtomTyping)) {
        std::cout << ">>> Extracted SMARTS: " << extracted_smarts << std::endl;
      }
      RDKit::STR_VECT variants2_taut;
      try {
        variants2_taut =
            sa.enumerate_variants(extracted_smarts, 20000, enumerate_verbose,
                                  false, workflow_options.enumerate_bond_order);
        // const auto variants2_taut = perceive_tautomers(variants2, verbose);

      } catch (const std::exception &e) {
        if (log_enabled(LogErrors)) {
          std::cout << "Error enumerating variants for '" << extracted_smarts
                    << "': " << e.what() << std::endl;
        }
        if (catchErrors) {
          throw std::runtime_error(
              "Error enumerating variants for: " + extracted_smarts + " - " +
              std::string(e.what()));
        } else {
          these_results.push_back(
              extracted_smarts);  // Fallback to extracted SMARTS if tautomer
                                  // enumeration fails.
          continue;
        }
      }

      std::vector<std::string> these_results_split;
      std::unordered_set<std::string> local_seen;
      local_seen.reserve(variants2_taut.size() * 2 + 1);
      for (const auto &v2_raw : variants2_taut) {
        std::string v2;
        {
          auto nnr_it = normalized_neg_recursive_cache.find(v2_raw);
          if (nnr_it != normalized_neg_recursive_cache.end()) {
            v2 = nnr_it->second;
          } else {
            // v2 = normalize_negated_recursive_queries(v2_raw);
            v2 = v2_raw;
            normalized_neg_recursive_cache.emplace(v2_raw, v2);
          }
        }

        if (!local_seen.insert(v2).second) {
          continue;
        }

        bool valence_ok = false;
        {
          auto vc_it = valence_cache.find(v2);
          if (vc_it != valence_cache.end()) {
            valence_ok = vc_it->second;
          } else {
            valence_ok =
                at.is_valid_valence_smarts(v2, log_enabled(LogValidation));
            valence_cache.emplace(v2, valence_ok);
          }
        }

        if (valence_ok) {
          if (log_enabled(LogValidation)) {
            std::cout << "\t" << v2 << std::endl;
          }
          // these_results.push_back(v2);
          these_results_split.push_back(v2);
        } else {
          if (log_enabled(LogValidation)) {
            std::cout << "  Invalid valence SMARTS: " << v2 << std::endl;
          }
          if (ignoreValence) {
            these_results_split.push_back(v2);
          }

        }
      }

      std::vector<std::string> first_consolidation;
      try {
        first_consolidation =
            at.consolidate_smarts_by_atom_maps_recanon(these_results_split);
      } catch (const std::exception &e) {
        if (log_enabled(LogErrors)) {
          std::cout << "Warning: recanon consolidation failed: " << e.what()
                    << " (falling back to split results)" << std::endl;
        }
        first_consolidation = these_results_split;
      }

      std::string final_smarts =
          at.consolidate_smarts_with_recursive_paths(first_consolidation);
      these_results.push_back(final_smarts);
      if (log_enabled(LogFinal)) {
        for (const auto &s : first_consolidation) {
          std::cout << "\t[first_consolidation] " << s << std::endl;
        }
        std::cout << " split console: " << final_smarts << std::endl;
      }
    }
    if (log_enabled(LogSummary)) {
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<\n\n " << std::endl;
      std::cout << "\t" << std::endl;
    }
    results.push_back(std::move(these_results));
  }
  return results;
}

std::vector<std::vector<std::string>> SmartsAnalyzer::generate_all_combos(
    std::vector<std::string> smarts_list, bool verbose,
    bool include_x_in_reserialization, bool ignoreValence, bool catchErrors,
    const StandardSmartsLogOptions &log_options) {
  StandardSmartsWorkflowOptions workflow_options;
  workflow_options.include_x_in_reserialization = include_x_in_reserialization;
  return generate_all_combos(smarts_list, verbose, ignoreValence, catchErrors,
                             workflow_options, log_options);
}

std::vector<std::string> SmartsAnalyzer::standard_smarts(
    const std::vector<std::string> &smarts_list, bool verbose,
    bool ignoreValence, bool catchErrors,
    const StandardSmartsWorkflowOptions &workflow_options,
    const StandardSmartsLogOptions &log_options) {
  ZoneScopedN("SmartsAnalyzer::standard_smarts");
  atom_typer::SmartsAnalyzer sa;
  atom_typer::AtomTyper at;

  const auto log_enabled = [&](unsigned int flag) {
    if (verbose) {
      return true;
    }
    return log_options.enabled && (log_options.flags & flag) != 0u;
  };
  std::vector<std::vector<std::string>> results;
  try {
    results =
        sa.generate_all_combos(smarts_list, verbose, ignoreValence, catchErrors,
                               workflow_options, log_options);
  } catch (const std::exception &e) {
    if (log_enabled(LogErrors)) {
      std::cout << "Error in generate_all_combos: " << e.what() << std::endl;
    }
    if (catchErrors) {
      throw std::runtime_error("Error in generate_all_combos: " +
                               std::string(e.what()));
    } else {
      return smarts_list;  // Fallback to original SMARTS list if processing
                           // fails.
    }
  }
  int idx = 0;
  std::vector<std::string> out;
  out.reserve(results.size());
  for (const auto &res : results) {
    std::string this_smarts = smarts_list[idx++];
    // std::vector<std::string> consolidated_pre_h =
    //     at.consolidate_smarts_by_atom_maps(res);
    // std::vector<std::string> h_merged;
    // h_merged.reserve(res.size());
    // for (const auto &s : res) {
    //   std::string h_merged_smarts = move_explicit_hydrogens_to_not_h0(s);
    //   h_merged.push_back(h_merged_smarts);
    //   if (verbose) {
    //     std::cout << "Pre H-merge: " << s << std::endl;
    //     std::cout << "Post H-merge: " << h_merged_smarts << std::endl;
    //   }
    // }
    if (log_enabled(LogSummary) || log_enabled(LogFinal)) {
      std::cout << ">>> input: " << this_smarts << std::endl;
      for (const auto &s : res) {
        std::cout << "\t" << s << std::endl;
      }
    }
    std::vector<std::string> non_empty;
    non_empty.reserve(res.size());
    for (const auto &s : res) {
      if (!s.empty()) {
        non_empty.push_back(s);
      }
    }

    // Keep ALL non-empty variants.  Different variants may originate from
    // different OR arms in the original SMARTS (e.g. [CH2,$([CH1][#6]),
    // $([CX4]([#6])[#6])]).  Filtering by max mapped-atom count would
    // discard shorter arms that match legitimately different structures.
    // The downstream consolidation steps (recanon + recursive_paths)
    // correctly merge variants with different topologies.
    std::vector<std::string> preferred;
    preferred.reserve(non_empty.size());
    for (const auto &s : non_empty) {
      preferred.push_back(remove_recursive_expression_atom_maps(s));
    }
    if (preferred.empty()) {
      preferred = res;
    }

    for (const auto &s : preferred) {
      if (log_enabled(LogSummary) || log_enabled(LogFinal)) {
        std::cout << "\t [preferred] " << s << std::endl;
      }
    }

    std::vector<std::string> final_smarts_1 =
        at.consolidate_smarts_by_atom_maps_recanon(preferred);
    // for (const auto &s : final_smarts_1){
    //   std::cout << "\t [test] " << s << std::endl;
    // }

    std::string final_smarts;
    if (!final_smarts_1.empty() && final_smarts_1.size() == 1) {
      final_smarts = final_smarts_1[0];
    } else {
      // fallback only when recanon cannot produce a single result
      const auto &fallback_inputs =
          final_smarts_1.empty() ? preferred : final_smarts_1;
      final_smarts = at.consolidate_smarts_with_recursive_paths(fallback_inputs);
    }
    // std::string final_smarts =
        // at.consolidate_smarts_with_recursive_paths(preferred);
    // if (max_mapped_atoms > 0 &&
    //     mapped_atom_count(final_smarts) < max_mapped_atoms) {
    //   auto fallback = at.consolidate_smarts_by_atom_maps(preferred);
    //   if (!fallback.empty()) {
    //     std::sort(fallback.begin(), fallback.end(),
    //               [&](const std::string &a, const std::string &b) {
    //                 const size_t am = mapped_atom_count(a);
    //                 const size_t bm = mapped_atom_count(b);
    //                 if (am != bm) {
    //                   return am > bm;
    //                 }
    //                 if (a.size() != b.size()) {
    //                   return a.size() < b.size();
    //                 }
    //                 return a < b;
    //               });
    //     final_smarts = fallback.front();
    //   }
    // }
    final_smarts = remove_recursive_expression_atom_maps(final_smarts);
    if (log_enabled(LogFinal)) {
      std::cout << "\t[final_pre_factor] " << final_smarts << std::endl;
    }
    try {
      final_smarts = factor_or_common_primitives(
          final_smarts, false, workflow_options.factoring_priority);
    } catch (const std::exception &e) {
      if (log_enabled(LogErrors)) {
        std::cout << "Warning: factor_or_common_primitives failed for '"
                  << final_smarts << "': " << e.what()
                  << " (keeping unfactored)" << std::endl;
      }
    } catch (...) {
      if (log_enabled(LogErrors)) {
        std::cout << "Warning: factor_or_common_primitives failed for '"
                  << final_smarts << "' (keeping unfactored)" << std::endl;
      }
    }
    try {
      final_smarts = collapse_equivalent_recursive_or_arms(final_smarts);
    } catch (const std::exception &e) {
      if (log_enabled(LogErrors)) {
        std::cout << "Warning: collapse_equivalent_recursive_or_arms failed for '"
                  << final_smarts << "': " << e.what()
                  << " (keeping as-is)" << std::endl;
      }
    } catch (...) {
      if (log_enabled(LogErrors)) {
        std::cout << "Warning: collapse_equivalent_recursive_or_arms failed for '"
                  << final_smarts << "' (keeping as-is)" << std::endl;
      }
    }

    try {
      final_smarts =
          normalize_smarts_encoding(final_smarts, workflow_options);
    } catch (const std::exception &e) {
      if (log_enabled(LogErrors)) {
        std::cout << "Warning: normalize_smarts_encoding failed for '"
                  << final_smarts << "': " << e.what()
                  << " (keeping as-is)" << std::endl;
      }
    } catch (...) {
      if (log_enabled(LogErrors)) {
        std::cout << "Warning: normalize_smarts_encoding failed for '"
                  << final_smarts << "' (keeping as-is)" << std::endl;
      }
    }

    if (log_enabled(LogFinal)) {
      // std::cout << "Input SMARTS: " << this_smarts << std::endl;
      //   std::cout << "Generated " << res.size() << " variants:\n";
      std::cout << "\t[final] " << final_smarts << std::endl;
      std::cout << std::endl;
    }
    out.push_back(final_smarts);
  }
  return out;
}

std::vector<std::string> SmartsAnalyzer::standard_smarts(
    const std::vector<std::string> &smarts_list, bool verbose,
    bool include_x_in_reserialization, bool ignoreValence, bool catchErrors,
    const StandardSmartsLogOptions &log_options) {
  StandardSmartsWorkflowOptions workflow_options;
  workflow_options.include_x_in_reserialization = include_x_in_reserialization;
  return standard_smarts(smarts_list, verbose, ignoreValence, catchErrors,
                         workflow_options, log_options);
}
}  // namespace atom_typer