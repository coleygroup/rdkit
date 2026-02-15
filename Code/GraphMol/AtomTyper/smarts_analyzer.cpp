#include "smarts_analyzer.hpp"  //header file for smarts_analyzer.cpp
#include "atom_typer.hpp"
#include <GraphMol/GraphMol.h>  // Core RDKit molecule functionality
#include <GraphMol/MolOps.h>    // For molecule operations from RDKit molecule
#include <GraphMol/SmilesParse/SmilesParse.h>  // For parsing SMILES strings
#include <GraphMol/SmilesParse/SmartsWrite.h>  // For writing SMARTS strings
#include <GraphMol/Descriptors/MolDescriptors.h>  // For molecular descriptors like the AND/OR/NOT queries
#include <GraphMol/QueryAtom.h>  // For QueryAtom and Query functionality
#include <GraphMol/QueryBond.h>
#include <GraphMol/QueryOps.h>    // For ATOM_EQUALS_QUERY and other query types
#include <Query/QueryObjects.h>   // For the base Query class
#include <Query/EqualityQuery.h>  // For Equality Query functionality
#include <sstream>  // For string stream operations, work with data types better
#include <stdexcept>  // For std exceptions
#include <algorithm>  // For std algorithms like std::max, std::min, std::sort, std::unique
#include <cctype>
#include <functional>
#include <map>
#include <set>

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
      std::cerr << "Warning: Failed to expand recursive SMARTS '"
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
    const std::string &smarts, int max, bool verbose, bool carry_atom_maps) {
  std::vector<std::string> results;
  if (verbose) {
    std::cerr << "[enumerate_variants] input='" << smarts << "' max=" << max
              << std::endl;
  }

  if (smarts.empty() || max <= 0) {
    return results;
  }

  std::unique_ptr<RDKit::ROMol> mol(RDKit::SmartsToMol(smarts));
  if (!mol) {
    throw std::runtime_error("Invalid SMARTS");
  }

  std::vector<unsigned int> input_atom_maps;
  input_atom_maps.reserve(mol->getNumAtoms());
  for (const auto *atom : mol->atoms()) {
    input_atom_maps.push_back(atom ? atom->getAtomMapNum() : 0U);
  }
  if (verbose) {
    std::cerr << "[enumerate_variants] parsed mol atoms=" << mol->getNumAtoms()
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
      std::cerr << "[enumerate_variants] atom idx=" << atom->getIdx()
                << " query_desc='" << qatom->getQuery()->getDescription()
                << "' variants=" << variants.size() << std::endl;
    }
    if (variants.empty()) {
      if (verbose) {
        std::cerr << "Warning: No variants generated for atom index "
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

  for (const auto *bond : mol->bonds()) {
    if (!bond->hasQuery() || !bond->getQuery()) {
      continue;
    }

    auto variants = pimpl->enumerate_bond_query_variants(bond->getQuery());
    if (verbose) {
      std::cerr << "[enumerate_variants] bond idx=" << bond->getIdx()
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

  if (variant_atom_indices.empty() && variant_bond_indices.empty()) {
    results.push_back(smarts);
    return results;
  }
  std::set<std::string> seen;
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

    while (true) {
      const auto qpos = inner.find("$(");
      if (qpos == std::string::npos) {
        break;
      }
      const size_t open_pos = qpos + 1;
      const size_t close_pos = matching_paren(inner, open_pos);
      if (close_pos == std::string::npos || close_pos <= open_pos + 1) {
        break;
      }

      const std::string recursive_smarts =
          inner.substr(open_pos + 1, close_pos - open_pos - 1);
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
    }

    while (inner.find("&&") != std::string::npos) {
      inner.replace(inner.find("&&"), 2, "&");
    }
    while (!inner.empty() && (inner.back() == '&' || inner.back() == ';')) {
      inner.pop_back();
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
          check.reset(RDKit::SmartsToMol(candidate));
        }
        if (check && seen.insert(candidate).second) {
          if (verbose) {
            std::cerr << "[enumerate_variants] candidate='" << candidate << "'"
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

  std::function<void(RDKit::ROMol &)> assign_maps_in_mol;
  std::function<void(const AtomQuery *)> visit_query;

  visit_query = [&](const AtomQuery *q) {
    if (!q) {
      return;
    }

    const auto *rsq = dynamic_cast<const RDKit::RecursiveStructureQuery *>(q);
    if (rsq && rsq->getQueryMol()) {
      auto *submol = const_cast<RDKit::ROMol *>(rsq->getQueryMol());
      assign_maps_in_mol(*submol);
    }

    for (auto it = q->beginChildren(); it != q->endChildren(); ++it) {
      visit_query(it->get());
    }
  };

  assign_maps_in_mol = [&](RDKit::ROMol &target) {
    if (visited.count(&target)) {
      return;
    }
    visited.insert(&target);

    for (auto *atom : target.atoms()) {
      atom->setAtomMapNum(next_map++);
    }

    for (auto *atom : target.atoms()) {
      const auto *qatom = dynamic_cast<const RDKit::QueryAtom *>(atom);
      if (!qatom || !qatom->getQuery()) {
        continue;
      }
      visit_query(qatom->getQuery());
    }
  };

  assign_maps_in_mol(*mol);
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

std::vector<std::vector<std::string>> SmartsAnalyzer::generate_all_combos(
    std::vector<std::string> smarts_list, bool verbose) {
  atom_typer::SmartsAnalyzer sa;
  atom_typer::AtomTyper at;

  std::vector<std::vector<std::string>> results;
  for (const auto &test1 : smarts_list) {
    std::vector<std::string> these_results;
    const std::string mapped_smarts = sa.add_atom_maps(test1);
    int max_amap = 0;
    auto mol = RDKit::SmartsToMol(mapped_smarts);
    for (const auto *atom : mol->atoms()) {
      int amap = atom->getAtomMapNum();
      if (amap > max_amap) {
        max_amap = amap;
      }
    }
    if (verbose) {
      std::cout << "Mapped SMARTS: " << mapped_smarts << std::endl;
    }

    const auto variants =
        sa.enumerate_variants(mapped_smarts, 1000, false, true);
    if (verbose) {
      std::cout << "Generated " << variants.size() << " variants:\n";
    }
    // const auto def = at.get_default_query_embedding();
    int intial_max_amap = max_amap;
    for (const auto &variant : variants) {
      max_amap = intial_max_amap;  // reset max_amap for next iteration
      const auto variants2 = sa.enumerate_variants(
          at.enumerate_dof_smarts(variant, false, max_amap), 1000, false,
          false);
      bool seen_invalid = false;
      for (const auto &v2 : variants2) {
        if (at.is_valid_valence_smarts(v2, false)) {
          if (verbose) {
            std::cout << "\t" << v2 << std::endl;
          }
          these_results.push_back(v2);
        } else {
          if (verbose) {
            std::cout << "  Invalid valence SMARTS: " << v2 << std::endl;
          }

          seen_invalid = true;
        }
      }
    }
    if (verbose) {
      std::cout << "\t" << std::endl;
    }
    results.push_back(std::move(these_results));
  }
  return results;
}

std::vector<std::string> SmartsAnalyzer::standard_smarts(
    const std::vector<std::string> &smarts_list, bool verbose) {
  atom_typer::SmartsAnalyzer sa;
  atom_typer::AtomTyper at;
  std::vector<std::vector<std::string>> results =
      sa.generate_all_combos(smarts_list, verbose);
  // std::cout << "Generated " << results.size() << " sets of variants:\n";
  int idx = 0;
  std::vector<std::string> out;
  out.reserve(results.size());
  for (const auto &res : results) {
    // for (const auto &s : res) {
    //   std::cout << s << std::endl;
    // }
    std::string this_smarts = smarts_list[idx++];
    std::vector<std::string> consolidated =
        at.consolidate_smarts_by_atom_maps(res);
    std::string final_smarts =
        at.consolidate_smarts_with_recursive_paths(consolidated);

    if (verbose) {
      std::cout << "Input SMARTS: " << this_smarts << std::endl;
      std::cout << "Generated " << res.size() << " variants:\n";
      std::cout << "\tConsolidated into " << consolidated.size()
                << " SMARTS:\n";
      for (const auto &s : consolidated) {
        std::cout << "\t\t" << s << std::endl;
      }
      std::cout << "\tFinal consolidated SMARTS: " << final_smarts << std::endl;
      std::cout << std::endl;
    }
    out.push_back(final_smarts);
  }
  return out;
}
}  // namespace atom_typer