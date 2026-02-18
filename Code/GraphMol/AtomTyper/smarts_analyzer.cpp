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
#include <GraphMol/Substruct/SubstructMatch.h>
#include <Query/QueryObjects.h>   // For the base Query class
#include <Query/EqualityQuery.h>  // For Equality Query functionality
#include <sstream>  // For string stream operations, work with data types better
#include <stdexcept>  // For std exceptions
#include <algorithm>  // For std algorithms like std::max, std::min, std::sort, std::unique
#include <cctype>
#include <functional>
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
      out.emplace_back(q);
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
    // Base case: all atoms processed
    if (out.size() >= max) {
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
      if (out.size() >= max) {
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

    auto variants =
        pimpl->enumerate_query_variants(qatom->getQuery(), carry_atom_maps);
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
            if (input_atom_maps[i] != 0U) {
              check->getAtomWithIdx(static_cast<unsigned int>(i))
                  ->setAtomMapNum(input_atom_maps[i]);
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

  ++recursive_canon_depth;
  try {
    atom_typer::SmartsAnalyzer sa;
    const auto out = sa.standard_smarts({inner}, false);
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

}  // namespace

std::vector<std::vector<std::string>> SmartsAnalyzer::generate_all_combos(
  std::vector<std::string> smarts_list, bool verbose, bool ignoreValence,
  bool catchErrors,
  const StandardSmartsWorkflowOptions &workflow_options,
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
    auto mol = RDKit::SmartsToMol(mapped_smarts);
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
    const auto variants = sa.enumerate_variants(
      mapped_smarts, 1000, enumerate_verbose, true,
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
            workflow_options.enumerate_bond_order);
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
      const auto variants2_taut = sa.enumerate_variants(
          extracted_smarts, 20000, enumerate_verbose, false,
          workflow_options.enumerate_bond_order);
      // const auto variants2_taut = perceive_tautomers(variants2, verbose);

      bool seen_invalid = false;
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
            v2 = normalize_negated_recursive_queries(v2_raw);
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

          seen_invalid = true;
        }
      }

      std::vector<std::string> removed_recursive_maps_consolidation;
      for (const auto &s : these_results_split) {
        std::string removed_map = remove_recursive_expression_atom_maps(s);
        if (removed_recursive_maps_consolidation.empty() ||
            std::find(removed_recursive_maps_consolidation.begin(),
                      removed_recursive_maps_consolidation.end(),
                      removed_map) ==
                removed_recursive_maps_consolidation.end()) {
          removed_recursive_maps_consolidation.push_back(removed_map);
        }
      }

      std::vector<std::string> first_consolidation =
          at.consolidate_smarts_by_atom_maps_recanon(removed_recursive_maps_consolidation);

      std::string final_smarts = at.consolidate_smarts_with_recursive_paths(
          first_consolidation);
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
  workflow_options.include_x_in_reserialization =
    include_x_in_reserialization;
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

  std::vector<std::vector<std::string>> results =
      sa.generate_all_combos(smarts_list, verbose, ignoreValence, catchErrors,
                 workflow_options, log_options);

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
    const auto mapped_atom_count = [](const std::string &s) -> size_t {
      if (s.empty()) {
        return 0;
      }
      std::unique_ptr<RDKit::ROMol> mol(RDKit::SmartsToMol(s));
      if (!mol) {
        return 0;
      }
      size_t count = 0;
      for (const auto *atom : mol->atoms()) {
        if (atom && atom->getAtomMapNum() > 0) {
          ++count;
        }
      }
      return count;
    };

    std::vector<std::string> non_empty;
    non_empty.reserve(res.size());
    for (const auto &s : res) {
      if (!s.empty()) {
        non_empty.push_back(s);
      }
    }

    size_t max_mapped_atoms = 0;
    for (const auto &s : non_empty) {
      max_mapped_atoms = std::max(max_mapped_atoms, mapped_atom_count(s));
    }

    std::vector<std::string> preferred;
    if (max_mapped_atoms > 0) {
      preferred.reserve(non_empty.size());
      for (const auto &s : non_empty) {
        if (mapped_atom_count(s) == max_mapped_atoms) {
          preferred.push_back(remove_recursive_expression_atom_maps(s));
        }
      }
    }
    if (preferred.empty()) {
      preferred = non_empty;
    }
    if (preferred.empty()) {
      preferred = res;
    }

    for (const auto &s : preferred) {
      if (log_enabled(LogSummary) || log_enabled(LogFinal)) {
        std::cout << "\t [preferred] " << s << std::endl;
      }
    }

    std::string final_smarts =
        at.consolidate_smarts_with_recursive_paths(preferred);
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
  workflow_options.include_x_in_reserialization =
    include_x_in_reserialization;
  return standard_smarts(smarts_list, verbose, ignoreValence, catchErrors,
             workflow_options, log_options);
}
}  // namespace atom_typer