#include "atom_typer.hpp"
#include "atom_typer_query_reorder.hpp"
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/RWMol.h>
#include <Query/EqualityQuery.h>
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <optional>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>

#if __has_include(<tracy/Tracy.hpp>)
#include <tracy/Tracy.hpp>
#else
#define ZoneScoped
#define ZoneScopedN(x)
#endif

namespace atom_typer {

/**
 * Private implementation class (PIMPL pattern)
 */
class AtomTyper::Impl {
 public:
  bool use_canonical = false;
  DebugLevel default_debug_level = DebugLevel::Off;
  bool recanon_comparison_trace = false;

  struct EnumerationSettings {
    int h_min = 0;
    int h_max = 4;
    int charge_min = -1;
    int charge_max = 1;
    DebugLevel debug_level = DebugLevel::Off;
    bool include_x_in_reserialization = false;
    bool enumerate_bond_order = true;
    unsigned int forced_primitives_mask = 0u;
  };

  static int debug_level_to_int(DebugLevel level) {
    return static_cast<int>(level);
  }

  struct QueryConstraints {
    std::optional<int> atomic_number;
    std::optional<int> formal_charge;
    std::optional<int> h_count;
    std::optional<int> lower_h_count;
    std::optional<int> explicit_degree_query;  // SMARTS D primitive
    std::optional<int> total_degree;
    std::optional<int> valence;
    std::optional<int> explicit_valence;
    std::optional<int> hybridization;
    std::optional<bool> aromatic;
    std::optional<bool> aliphatic;
    std::optional<bool> in_ring;
    std::optional<int> ring_count;
    std::optional<int> min_ring_size;
    int min_bonds = 0;  // explicit topology neighbors in SMARTS graph
    std::vector<int> excluded_formal_charge;
    std::vector<int> excluded_h_count;
    std::vector<int> excluded_explicit_degree_query;
    std::vector<int> excluded_total_degree;
    bool excluded_aromatic = false;
    bool excluded_aliphatic = false;
  };

  void collect_query_constraints(const RDKit::Atom::QUERYATOM_QUERY *query,
                                 QueryConstraints &constraints,
                                 bool inherited_negation = false) {
    if (!query) {
      return;
    }

    const auto &descr = query->getDescription();
    const bool self_negated = query->getNegation();
    const bool effective_negated = inherited_negation ^ self_negated;

    if (descr == "AtomNot" || descr == "Not") {
      for (auto it = query->beginChildren(); it != query->endChildren(); ++it) {
        collect_query_constraints(it->get(), constraints, !effective_negated);
      }
      return;
    }

    if (!effective_negated) {
      auto eq = dynamic_cast<const RDKit::ATOM_EQUALS_QUERY *>(query);
      if ((descr == "AtomAtomicNum" || descr == "AtomNum") && eq) {
        constraints.atomic_number = eq->getVal();
      } else if (descr == "AtomFormalCharge" && eq) {
        constraints.formal_charge = eq->getVal();
      } else if (descr == "AtomHCount" && eq) {
        if (constraints.lower_h_count.has_value()) {
          constraints.h_count = eq->getVal();
        } else {
          constraints.h_count = eq->getVal();
        }
      } else if (descr == "AtomImplicitHCount" && eq) {
        // if (constraints.h_count.has_value()) {
        //   constraints.h_count =  eq->getVal();
        // }
        constraints.lower_h_count = eq->getVal();
      } else if (descr == "AtomExplicitDegree" && eq) {
        constraints.explicit_degree_query = eq->getVal();
      } else if (descr == "AtomTotalDegree" && eq) {
        constraints.total_degree = eq->getVal();
      } else if ((descr == "AtomTotalValence" ||
                  descr == "AtomExplicitValence") &&
                 eq) {
        constraints.valence = eq->getVal();
      } else if (descr == "AtomHybridization" && eq) {
        constraints.hybridization = eq->getVal();
      } else if (descr == "AtomAromatic") {
        constraints.aromatic = true;
      } else if (descr == "AtomAliphatic") {
        constraints.aliphatic = true;
      } else if (descr == "AtomInRing") {
        constraints.in_ring = true;
      } else if ((descr == "AtomInNRings" || descr == "AtomRingCount") && eq) {
        if (eq->getVal() < 0) {
          // R without a number: val=-1 means "in any ring"
          constraints.in_ring = true;
        } else {
          constraints.ring_count = eq->getVal();
          constraints.in_ring = eq->getVal() > 0;
        }
      } else if (descr == "AtomMinRingSize" && eq) {
        constraints.min_ring_size = eq->getVal();
        constraints.in_ring = true;
      }
    } else {
      auto eq = dynamic_cast<const RDKit::ATOM_EQUALS_QUERY *>(query);
      if (descr == "AtomFormalCharge" && eq) {
        constraints.excluded_formal_charge.push_back(eq->getVal());
      } else if ((descr == "AtomHCount" || descr == "AtomImplicitHCount") &&
                 eq) {
        constraints.excluded_h_count.push_back(eq->getVal());
      } else if (descr == "AtomExplicitDegree" && eq) {
        constraints.excluded_explicit_degree_query.push_back(eq->getVal());
      } else if (descr == "AtomTotalDegree" && eq) {
        constraints.excluded_total_degree.push_back(eq->getVal());
      } else if (descr == "AtomAromatic") {
        constraints.excluded_aromatic = true;
      } else if (descr == "AtomAliphatic") {
        constraints.excluded_aliphatic = true;
      } else if (descr == "AtomInRing") {
        // Negated AtomInRing means "not in ring" → in_ring = false
        constraints.in_ring = false;
      } else if ((descr == "AtomInNRings" || descr == "AtomRingCount") && eq) {
        if (eq->getVal() < 0) {
          // Negated R (no number): val=-1 means "not in any ring"
          constraints.in_ring = false;
        } else if (eq->getVal() == 0) {
          // Negated ring-count: !R0 means "in a ring"
          constraints.in_ring = true;
        }
      }
    }

    for (auto it = query->beginChildren(); it != query->endChildren(); ++it) {
      collect_query_constraints(it->get(), constraints, effective_negated);
    }
  }

  std::string infer_hybridization_from_bonds(const AtomType &at) {
    if (at.num_triple_bonds > 0) {
      return "SP";
    }
    if (at.num_double_bonds > 0 || at.num_aromatic_bonds > 0) {
      return "SP2";
    }
    if (at.num_single_bonds > 0) {
      return "SP3";
    }
    return "UNSPECIFIED";
  }

  std::string hybridization_to_string(int hyb) {
    switch (static_cast<RDKit::Atom::HybridizationType>(hyb)) {
      case RDKit::Atom::SP:
        return "SP";
      case RDKit::Atom::SP2:
        return "SP2";
      case RDKit::Atom::SP3:
        return "SP3";
      case RDKit::Atom::SP3D:
        return "SP3D";
      case RDKit::Atom::SP3D2:
        return "SP3D2";
      default:
        return "UNSPECIFIED";
    }
  }

  std::string bond_type_to_string(RDKit::Bond::BondType bt) {
    switch (bt) {
      case RDKit::Bond::BondType::SINGLE:
        return "SINGLE";
      case RDKit::Bond::BondType::DOUBLE:
        return "DOUBLE";
      case RDKit::Bond::BondType::TRIPLE:
        return "TRIPLE";
      case RDKit::Bond::BondType::AROMATIC:
        return "AROMATIC";
      default:
        return "OTHER";
    }
  }

  BondType make_bond_type(const RDKit::Bond *bond) {
    BondType bt;
    bt.bond_idx = static_cast<int>(bond->getIdx());
    bt.begin_atom_idx = static_cast<int>(bond->getBeginAtomIdx());
    bt.end_atom_idx = static_cast<int>(bond->getEndAtomIdx());
    bt.bond_type_enum = static_cast<int>(bond->getBondType());
    bt.bond_type_name = bond_type_to_string(bond->getBondType());
    bt.smarts_pattern =
        RDKit::SmartsWrite::GetBondSmarts(bond, bt.begin_atom_idx);
    bt.is_aromatic = bond->getIsAromatic();
    bt.is_in_ring =
        bond->getOwningMol().getRingInfo()->numBondRings(bond->getIdx()) > 0;
    return bt;
  }

  double bond_order_from_type(RDKit::Bond::BondType bt) {
    switch (bt) {
      case RDKit::Bond::BondType::SINGLE:
        return 1.0;
      case RDKit::Bond::BondType::DOUBLE:
        return 2.0;
      case RDKit::Bond::BondType::TRIPLE:
        return 3.0;
      case RDKit::Bond::BondType::AROMATIC:
        return 1.5;
      default:
        return 1.0;
    }
  }

  double used_connectivity_valence(RDKit::Atom *atom, RDKit::ROMol *mol) {
    double used = static_cast<double>(atom->getNumExplicitHs());
    for (const auto &bond : mol->atomBonds(atom)) {
      used += bond_order_from_type(bond->getBondType());
    }
    return used;
  }

  int default_valence_from_api(int atomic_number, int formal_charge) {
    const auto *ptable = RDKit::PeriodicTable::getTable();
    if (atomic_number <= 0 ||
        atomic_number > static_cast<int>(ptable->getMaxAtomicNumber())) {
      return -1;
    }
    const auto &ovalens = ptable->getValenceList(atomic_number);
    // for (int ovalen : ovalens) {
    //   std::cout << ovalen << " ";
    // }
    // std::cout << std::endl;
    unsigned int effective_atomic_num =
        static_cast<unsigned int>(atomic_number);
    if (ovalens.size() > 1 || ovalens[0] != -1) {
      int effective = atomic_number - formal_charge;
      effective = std::clamp(effective, 0,
                             static_cast<int>(ptable->getMaxAtomicNumber()));
      effective_atomic_num = static_cast<unsigned int>(effective);
    }

    return ptable->getDefaultValence(effective_atomic_num);
  }

  int calculate_valence_from_connectivity(RDKit::Atom *atom,
                                          RDKit::ROMol *mol) {
    double accum = used_connectivity_valence(atom, mol);
    accum += 0.1;
    return std::max(0, static_cast<int>(std::round(accum)));
  }

  unsigned int remaining_valence_from_connectivity(int atomic_number,
                                                   int formal_charge,
                                                   RDKit::Atom *atom,
                                                   RDKit::ROMol *mol) {
    if (atomic_number <= 0) {
      return 0U;
    }
    const int default_valence =
        default_valence_from_api(atomic_number, formal_charge);
    if (default_valence < 0) {
      return 0U;
    }
    const double used = used_connectivity_valence(atom, mol);
    return static_cast<unsigned int>(
        std::max(0, static_cast<int>(std::floor(default_valence - used))));
  }

  /**
   * Extract atom type information from an RDKit atom
   */
  AtomType extract_atom_type(RDKit::Atom *atom, RDKit::ROMol *mol, int idx) {
    AtomType at;
    at.atom_idx = idx;
    at.atomic_number = atom->getAtomicNum();
    at.formal_charge = atom->getFormalCharge();
    at.num_hydrogens = atom->getTotalNumHs();
    at.min_bonds = atom->getDegree();
    const int default_valence =
        default_valence_from_api(at.atomic_number, at.formal_charge);
    at.max_valence = default_valence >= 0
                         ? default_valence
                         : calculate_valence_from_connectivity(atom, mol);
    at.is_aromatic = atom->getIsAromatic();
    at.is_aliphatic = !at.is_aromatic;
    at.is_in_ring = mol->getRingInfo()->numAtomRings(idx) > 0;

    // Initialize ring-related fields
    at.ring_size = 0;
    at.num_ring_bonds = 0;
    at.num_aliphatic_rings = 0;
    at.num_aromatic_rings = 0;
    at.ring_connectivity = 0;
    at.atom_type_enumeration =
        RDKit::makeAtomType(at.atomic_number, at.is_aromatic);
    at.num_single_bonds = 0;
    at.num_double_bonds = 0;
    at.num_triple_bonds = 0;
    at.num_aromatic_bonds = 0;
    at.remaining_valence = static_cast<int>(remaining_valence_from_connectivity(
        at.atomic_number, at.formal_charge, atom, mol));
    at.bond_details.clear();
    at.source_atom_smarts = RDKit::SmartsWrite::GetAtomSmarts(atom);
    // ring_connectivity counts neighbors that are also in rings
    for (const auto &nbr : mol->atomNeighbors(atom)) {
      if (mol->getRingInfo()->numAtomRings(nbr->getIdx()) > 0) {
        at.ring_connectivity++;
      }
    }

    // find number of aliphatic/aromatic rings
    // start both counts at 0
    if (!atom->getIsAromatic()) {
      // the atom is aliphatic, so all rings aliphatic
      at.num_aliphatic_rings = mol->getRingInfo()->numAtomRings(idx);
      at.num_aromatic_rings = 0;
    } else {
      at.num_aliphatic_rings = 0;
      at.num_aromatic_rings = 0;

      for (const auto &ring : mol->getRingInfo()->atomRings()) {
        if (std::find(ring.begin(), ring.end(), idx) != ring.end()) {
          // check if ring is aromatic
          bool is_aromatic_ring = true;
          for (int atom_id : ring) {
            RDKit::Atom *ring_atom = mol->getAtomWithIdx(atom_id);
            if (!ring_atom->getIsAromatic()) {
              is_aromatic_ring = false;
              break;
            }
          }

          // add if aromatic or aliphatic
          if (is_aromatic_ring) {
            at.num_aromatic_rings += 1;
          } else {
            at.num_aliphatic_rings += 1;
          }
        }
      }
    }

    // Get smallest ring size
    // Note: ring_size is stored as int to match RDKit's AtomType convention
    // and to maintain API consistency. Max ring size in practice is << INT_MAX.
    at.ring_size = 0;

    // declare smallest ring
    std::vector<int> smallest_ring;
    if (at.is_in_ring) {
      const auto &ring_info = mol->getRingInfo();
      for (const auto &ring : ring_info->atomRings()) {
        if (std::find(ring.begin(), ring.end(), idx) != ring.end()) {
          // find the smallest ring/size
          int current_ring_size = static_cast<int>(ring.size());
          if (at.ring_size == 0 || current_ring_size < at.ring_size) {
            at.ring_size = current_ring_size;
            smallest_ring = ring;
          }
        }
      }
    }
    // use smallest_ring to get the member atoms' indices
    for (int ind : smallest_ring) {
      at.ring_membership_list.push_back(ind);
    }

    // find number of ring bonds (bonds that are in ANY ring, not just smallest)
    at.num_ring_bonds = 0;
    for (const auto &bond : mol->atomBonds(atom)) {
      if (mol->getRingInfo()->numBondRings(bond->getIdx()) > 0) {
        at.num_ring_bonds++;
      }
    }

    // Get hybridization
    switch (atom->getHybridization()) {
      case RDKit::Atom::SP:
        at.hybridization = "SP";
        break;
      case RDKit::Atom::SP2:
        at.hybridization = "SP2";
        break;
      case RDKit::Atom::SP3:
        at.hybridization = "SP3";
        break;
      case RDKit::Atom::SP3D:
        at.hybridization = "SP3D";
        break;
      case RDKit::Atom::SP3D2:
        at.hybridization = "SP3D2";
        break;
      default:
        at.hybridization = "UNSPECIFIED";
    }

    // Get Chirality using R/S notation
    // Use RDKit's CIP assignment to get R/S labels
    std::string cip_code;
    if (atom->hasProp("_CIPCode")) {
      atom->getProp("_CIPCode", cip_code);
      at.chirality = cip_code;  // Will be "R" or "S"
    } else {
      // No CIP code assigned
      switch (atom->getChiralTag()) {
        case RDKit::Atom::CHI_TETRAHEDRAL_CW:
        case RDKit::Atom::CHI_TETRAHEDRAL_CCW:
          at.chirality = "UNSPECIFIED";  // Has chirality but no R/S assignment
          break;
        default:
          at.chirality = "UNSPECIFIED";
      }
    }

    // Get neighbors & bond types
    std::map<int, int> bond_type_counts;
    for (const auto &nbr : mol->atomNeighbors(atom)) {
      at.neighbors.push_back(nbr->getIdx());

      RDKit::Bond *bond = mol->getBondBetweenAtoms(idx, nbr->getIdx());
      if (bond) {
        int bond_type_int = static_cast<int>(bond->getBondType());
        bond_type_counts[bond_type_int]++;
        at.bond_details.push_back(make_bond_type(bond));

        switch (bond->getBondType()) {
          case RDKit::Bond::BondType::SINGLE:
            at.num_single_bonds++;
            break;
          case RDKit::Bond::BondType::DOUBLE:
            at.num_double_bonds++;
            break;
          case RDKit::Bond::BondType::TRIPLE:
            at.num_triple_bonds++;
            break;
          case RDKit::Bond::BondType::AROMATIC:
            at.num_aromatic_bonds++;
            break;
          default: {
            // For SMARTS query bonds, getBondType() may return OTHER or
            // UNSPECIFIED even for explicit = / # bonds.  Fall back to
            // the bond's SMARTS string to classify them.
            const std::string bs = RDKit::SmartsWrite::GetBondSmarts(
                bond, static_cast<int>(bond->getBeginAtomIdx()));
            if (bs == "=" || bs == "=") {
              at.num_double_bonds++;
            } else if (bs == "#") {
              at.num_triple_bonds++;
            } else if (bs == ":" || bs == "~") {
              at.num_aromatic_bonds++;
            } else {
              at.num_single_bonds++;
            }
            break;
          }
        }
      }
    }
    at.bond_types = bond_type_counts;

    // Generate SMARTS pattern for this atom type
    at.smarts_pattern = generate_smarts_pattern(at);

    return at;
  }

  /**
   * Extract atom type information directly from SMARTS query constraints.
   * This avoids inferred atom properties (e.g. [C] getting H=4).
   */
  AtomType extract_smarts_atom_type(RDKit::Atom *atom, RDKit::ROMol *mol,
                                    int idx) {
    auto at = extract_atom_type(atom, mol, idx);
    at.source_atom_smarts = RDKit::SmartsWrite::GetAtomSmarts(atom);
    if (!atom->hasQuery()) {
      return at;
    }

    QueryConstraints constraints;
    collect_query_constraints(atom->getQuery(), constraints);

    const auto infer_symbol_aromaticity_flags = [](const std::string &smarts) {
      std::string token = smarts;
      if (token.size() >= 2 && token.front() == '[' && token.back() == ']') {
        token = token.substr(1, token.size() - 2);
      }

      const auto colon = token.rfind(':');
      if (colon != std::string::npos && colon + 1 < token.size()) {
        bool digits_only = true;
        for (size_t i = colon + 1; i < token.size(); ++i) {
          if (!std::isdigit(static_cast<unsigned char>(token[i]))) {
            digits_only = false;
            break;
          }
        }
        if (digits_only) {
          token = token.substr(0, colon);
        }
      }

      const auto first_delim = token.find_first_of(";,&");
      if (first_delim != std::string::npos) {
        token = token.substr(0, first_delim);
      }

      bool was_negated = false;
      while (!token.empty() && token.front() == '!') {
        token.erase(token.begin());
        was_negated = true;
      }

      // A negated atom constraint (e.g. !#1 = "not hydrogen") could match
      // any element, so the atom can be either aromatic or aliphatic.
      if (was_negated) {
        return std::pair<std::optional<bool>, std::optional<bool>>{true, true};
      }

      // lowercase aromatic symbols
      if (token == "c" || token == "n" || token == "o" || token == "p" ||
          token == "s" || token == "as" || token == "se" || token == "a") {
        return std::pair<std::optional<bool>, std::optional<bool>>{true, false};
      }

      // uppercase element symbols => aliphatic
      if (token == "C" || token == "N" || token == "O" || token == "P" ||
          token == "S" || token == "As" || token == "Se" || token == "A" ||
          token == "F" || token == "Cl" || token == "Br" || token == "I") {
        return std::pair<std::optional<bool>, std::optional<bool>>{false, true};
      }

      // atomic number representation for C/N/O/P/S => potentially both
      if (token == "#6" || token == "#7" || token == "#8" || token == "#15" ||
          token == "#16") {
        return std::pair<std::optional<bool>, std::optional<bool>>{true, true};
      }

      // default fallback
      return std::pair<std::optional<bool>, std::optional<bool>>{false, true};
    };

    const auto inferred_arom_flags =
        infer_symbol_aromaticity_flags(at.source_atom_smarts);

    // std::cout << "Inferred aromaticity flags: " << at.source_atom_smarts
    //           << " -> " << inferred_arom_flags.first.value_or("None") << ", "
    //           << inferred_arom_flags.second.value_or("None") << std::endl;
    at.explicit_atomic_num = constraints.atomic_number;
    at.explicit_charge = constraints.formal_charge;
    // if (constraints.h_count.has_value() &&
    // constraints.lower_h_count.has_value()) {
    //   if (*constraints.h_count == *constraints.lower_h_count) {
    //     at.explicit_H = NULL;
    //   } else {
    //     at.explicit_H = *constraints.h_count - *constraints.lower_h_count;
    //   }
    // }
    at.explicit_H = constraints.h_count;
    at.explicit_lower_h = constraints.lower_h_count;
    at.explicit_D = constraints.explicit_degree_query;
    at.explicit_X = constraints.total_degree;
    at.explicit_valence = constraints.valence;
    at.explicit_hybridization = constraints.hybridization;
    at.explicit_aromatic = constraints.aromatic;
    at.explicit_aliphatic = constraints.aliphatic;
    at.inferred_aromatic = inferred_arom_flags.first;
    at.inferred_aliphatic = inferred_arom_flags.second;
    at.explicit_in_ring = constraints.in_ring;
    at.explicit_ring_count = constraints.ring_count;
    at.explicit_min_ring_size = constraints.min_ring_size;
    at.excluded_charges = constraints.excluded_formal_charge;
    at.excluded_h_counts = constraints.excluded_h_count;
    at.excluded_D_values = constraints.excluded_explicit_degree_query;
    at.excluded_X_values = constraints.excluded_total_degree;
    at.excluded_aromatic = constraints.excluded_aromatic;
    at.excluded_aliphatic = constraints.excluded_aliphatic;

    std::sort(at.excluded_charges.begin(), at.excluded_charges.end());
    at.excluded_charges.erase(
        std::unique(at.excluded_charges.begin(), at.excluded_charges.end()),
        at.excluded_charges.end());
    std::sort(at.excluded_h_counts.begin(), at.excluded_h_counts.end());
    at.excluded_h_counts.erase(
        std::unique(at.excluded_h_counts.begin(), at.excluded_h_counts.end()),
        at.excluded_h_counts.end());
    std::sort(at.excluded_D_values.begin(), at.excluded_D_values.end());
    at.excluded_D_values.erase(
        std::unique(at.excluded_D_values.begin(), at.excluded_D_values.end()),
        at.excluded_D_values.end());
    std::sort(at.excluded_X_values.begin(), at.excluded_X_values.end());
    at.excluded_X_values.erase(
        std::unique(at.excluded_X_values.begin(), at.excluded_X_values.end()),
        at.excluded_X_values.end());

    if (constraints.atomic_number.has_value()) {
      at.atomic_number = *constraints.atomic_number;
    }
    if (!at.explicit_charge.has_value()) {
      at.max_valence = default_valence_from_api(at.atomic_number, 0);
    } else {
      const int default_valence =
          default_valence_from_api(at.atomic_number, at.formal_charge);
      at.max_valence = default_valence;
    }
    at.hybridization = infer_hybridization_from_bonds(at);
    if (at.explicit_aromatic.has_value()) {
      at.is_aromatic = *at.explicit_aromatic;
    } else {
      at.is_aromatic = at.inferred_aromatic.value_or(at.is_aromatic);
    }
    if (at.explicit_aliphatic.has_value()) {
      at.is_aliphatic = *at.explicit_aliphatic;
    } else {
      at.is_aliphatic = at.inferred_aliphatic.value_or(at.is_aliphatic);
    }
    if (!at.is_aromatic && !at.is_aliphatic) {
      at.is_aliphatic = true;
    }
    // }
    // if (has_aliphatic_constraint) {
    // }

    if (constraints.in_ring.has_value()) {
      at.is_in_ring = *constraints.in_ring;
      if (!at.is_in_ring) {
        at.ring_size = 0;
        at.ring_membership_list.clear();
        at.num_ring_bonds = 0;
      }
    }
    if (constraints.min_ring_size.has_value()) {
      at.ring_size = *constraints.min_ring_size;
    }

    at.atom_type_enumeration =
        RDKit::makeAtomType(at.atomic_number, at.is_aromatic);

    at.remaining_valence = at.max_valence -
                           calculate_valence_from_connectivity(atom, mol) -
                           at.explicit_H.value_or(0);

    at.smarts_pattern = generate_smarts_pattern(at);
    return at;
  }

  /**
   * Generate a SMARTS pattern based on atom type
   */
  std::string generate_smarts_pattern(const AtomType &at) {
    std::stringstream ss;
    ss << "[";

    // Add atomic number
    if (at.atomic_number == 6) {
      ss << "C";
    } else if (at.atomic_number == 7) {
      ss << "N";
    } else if (at.atomic_number == 8) {
      ss << "O";
    } else if (at.atomic_number == 16) {
      ss << "S";
    } else if (at.atomic_number == 1) {
      ss << "H";
    } else {
      ss << "#" << at.atomic_number;
    }
    ss << ";";
    // Add degree
    // Add hydrogen count
    if (at.num_hydrogens > 0) {
      ss << "H" << at.num_hydrogens;
      ss << ";";
    }
    // Add charge
    if (at.formal_charge != 0) {
      if (at.formal_charge > 0) {
        ss << "+" << at.formal_charge;
      } else {
        ss << at.formal_charge;
      }
      ss << ";";
    }
    // Add aromaticity
    if (at.is_aromatic && at.is_aliphatic) {
      ss << "A,a";
      ss << ";";

    } else if (at.is_aromatic) {
      ss << "a";
      ss << ";";

    } else if (at.is_aliphatic) {
      ss << "A";
      ss << ";";
    }
    // Add ring membership
    if (at.is_in_ring) {
      ss << "R";
      if (at.ring_size > 0) {
        ss << at.ring_size;
      }
    }
    if (ss.str().back() == ';') {
      ss.seekp(-1, std::ios_base::end);  // Remove trailing semicolon
    }
    ss << "]";
    return ss.str();
  }

  std::string atom_symbol_for_smarts(const AtomType &at) {
    if (at.is_aromatic && !at.is_aliphatic) {
      if (at.atomic_number == 6) {
        return "c";
      }
      if (at.atomic_number == 7) {
        return "n";
      }
      if (at.atomic_number == 8) {
        return "o";
      }
      if (at.atomic_number == 16) {
        return "s";
      }
      if (at.atomic_number == 15) {
        return "p";
      }
    }

    if (at.atomic_number == 6) {
      return "C";
    }
    if (at.atomic_number == 7) {
      return "N";
    }
    if (at.atomic_number == 8) {
      return "O";
    }
    if (at.atomic_number == 16) {
      return "S";
    }
    if (at.atomic_number == 15) {
      return "P";
    }
    if (at.atomic_number == 1) {
      return "H";
    }
    return "#" + std::to_string(at.atomic_number);
  }

  std::vector<std::string> atom_symbols_for_smarts(const AtomType &at,
                                                   int aromaticity_flag) {
    if (aromaticity_flag == 1) {
      if (at.atomic_number == 6) {
        return {"c"};
      } else if (at.atomic_number == 7) {
        return {"n"};
      } else if (at.atomic_number == 8) {
        return {"o"};
      } else if (at.atomic_number == 16) {
        return {"s"};
      } else if (at.atomic_number == 15) {
        return {"p"};
      } else {
        return {atom_symbol_for_smarts(at), "a"};
      }
    }
    return {atom_symbol_for_smarts(at)};
  }

  std::string charge_token(int charge) {
    if (charge >= 0) {
      return "+" + std::to_string(charge);
    }
    return std::to_string(charge);
  }

  bool is_valence_consistent(const AtomType &at, int h_count, int charge,
                             int single_bonds, int double_bonds,
                             int triple_bonds, int aromatic_bonds, bool arom,
                             DebugLevel debug_level) {
    RDKit::RWMol temp_mol;

    // --- create the central atom ---
    auto *temp_atom = new RDKit::Atom(at.atomic_number);
    temp_atom->setFormalCharge(charge);
    temp_atom->setNoImplicit(true);
    temp_atom->setNumExplicitHs(
        static_cast<unsigned int>(std::max(0, h_count)));

    const bool want_aromatic = arom;
    temp_atom->setIsAromatic(want_aromatic);

    const unsigned int aidx = temp_mol.addAtom(temp_atom, false, true);

    // --- if aromatic, embed into a minimal aromatic ring (benzene-like) ---
    // This guarantees sanitization/aromatic valence bookkeeping is possible.
    if (want_aromatic) {
      int modeled_aromatic_bonds = aromatic_bonds;
      // Query fragments may expose only one aromatic arm in the explicit
      // SMARTS graph (e.g. terminal [c] in an aromatic path), while the
      // matched atom still has two aromatic ring neighbors in the target.
      // Model that underconstrained case as a minimal aromatic environment.
      if (modeled_aromatic_bonds == 1) {
        modeled_aromatic_bonds = 2;
      }

      if (at.atomic_number == 6){
        if (charge == 0 && h_count == 1 && modeled_aromatic_bonds == 2 && double_bonds == 0 && triple_bonds == 0) {
          // Furan-like: 2-electron lone-pair donor in a 5-member ring.
          return true;
        }
        if (charge == 0 && h_count == 0 && modeled_aromatic_bonds == 2 && double_bonds == 1 && triple_bonds == 0) {
          // *c(=O)*
          return true;
        }
        else if (charge == 0 && h_count == 1 && single_bonds == 2 && double_bonds == 0 && triple_bonds == 0) {
          // Furan-like: 2-electron lone-pair donor in a 5-member ring.
          return true;
        } else if (charge == 0 && h_count == 1 && modeled_aromatic_bonds == 2 && single_bonds == 0 && double_bonds == 0 && triple_bonds == 0) {
          // Pyrrole-like: 1-electron lone-pair donor in a 5-member ring.
          return true;
        }
      }

      // In an aromatic ring, the atom will have 2 aromatic bonds in-ring.
      // If your counting expects something else, you'll need a fused-ring
      // builder. For now, we treat ring bonds as satisfying aromaticity
      // requirements.
      if (modeled_aromatic_bonds != 0 && modeled_aromatic_bonds != 2 &&
          modeled_aromatic_bonds != 3) {
        // Conservative: refuse cases we can't model without fused systems.
        // (You can relax this if you implement fused-ring construction.)
        return false;
      }

      // For aromatic heteroatoms, a 5-member aromatic context is often a
      // better minimal probe than benzene-like 6-member embedding.
      // - [nH]-like atoms (h_count > 0): pyrrole / imidazole nitrogen
      // - Lone-pair donors with D2,H0,+0: furan oxygen, thiophene sulfur,
      //   selenophene selenium — these contribute 2 electrons to the π-system
      //   via a lone pair and naturally occur in 5-member aromatic rings.
      int aromatic_ring_size = 6;
      const bool is_lp_donor_heteroatom =
          (at.atomic_number == 8 || at.atomic_number == 16 ||
           at.atomic_number == 34);
      if ((at.atomic_number == 7 || is_lp_donor_heteroatom) &&
          h_count > 0) {
        aromatic_ring_size = 5;
      } else if (is_lp_donor_heteroatom && h_count == 0 && charge == 0 &&
                 modeled_aromatic_bonds == 2) {
        // Furan-like: 2-electron lone-pair donor in a 5-member ring.
        aromatic_ring_size = 5;
      }

      if (at.atomic_number == 8 || at.atomic_number == 7){
        if (charge == 0 && h_count == 0 && modeled_aromatic_bonds == 2 && double_bonds == 0 && triple_bonds == 0) {
          // Furan-like: 2-electron lone-pair donor in a 5-member ring.
          return true;
        } else if (charge == 0 && h_count == 1 && modeled_aromatic_bonds == 2 && single_bonds == 0 && double_bonds == 0 && triple_bonds == 0) {
          // Pyrrole-like: 1-electron lone-pair donor in a 5-member ring.
          return true;
        }
      }

      std::vector<unsigned int> ring(static_cast<size_t>(aromatic_ring_size));
      ring[0] = aidx;

      for (int i = 1; i < aromatic_ring_size; ++i) {
        auto *c = new RDKit::Atom(6);  // aromatic carbon
        c->setIsAromatic(true);
        ring[i] = temp_mol.addAtom(c, false, true);
      }

      // Connect ring with aromatic bonds.
      for (int i = 0; i < aromatic_ring_size; ++i) {
        int j = (i + 1) % aromatic_ring_size;
        temp_mol.addBond(ring[i], ring[j], RDKit::Bond::BondType::AROMATIC);
        auto *b = temp_mol.getBondBetweenAtoms(ring[i], ring[j]);
        b->setIsAromatic(true);
      }

      // Approximate fused aromatic environments where the probed atom has
      // three aromatic neighbors (e.g., bridgehead-like aromatic atoms in
      // fused systems). Add one extra aromatic neighbor to atom 0.
      if (modeled_aromatic_bonds == 3) {
        auto *c = new RDKit::Atom(6);
        c->setIsAromatic(true);
        const unsigned int fidx = temp_mol.addAtom(c, false, true);
        temp_mol.addBond(aidx, fidx, RDKit::Bond::BondType::AROMATIC);
        auto *b = temp_mol.getBondBetweenAtoms(aidx, fidx);
        if (b) {
          b->setIsAromatic(true);
        }
      }
    }

    // --- add exocyclic neighbors/bonds on the central atom ---
    auto add_nbrs = [&](int n, RDKit::Bond::BondType bt) {
      for (int i = 0; i < n; ++i) {
        auto *nbr = new RDKit::Atom(
            0);  // atomic number 0 = wildcard/query-like; safer than "*"
        unsigned int nidx = temp_mol.addAtom(nbr, false, true);
        temp_mol.addBond(aidx, nidx, bt);
      }
    };

    add_nbrs(single_bonds, RDKit::Bond::BondType::SINGLE);
    add_nbrs(double_bonds, RDKit::Bond::BondType::DOUBLE);
    add_nbrs(triple_bonds, RDKit::Bond::BondType::TRIPLE);
    add_nbrs(at.explicit_lower_h.value_or(0), RDKit::Bond::BondType::SINGLE);

    // IMPORTANT: do NOT add aromatic bonds to standalone neighbors here.
    // Aromatic bond types generally require both atoms to be aromatic *and*
    // part of an aromatic system; otherwise sanitization/property cache will
    // complain. We treat aromaticity as handled by the ring embedding above.
    // if (debug_level_to_int(debug_level) >=
    //     debug_level_to_int(DebugLevel::Verbose)) {
    //   std::cout << "Constructed probe molecule: "
    //             << RDKit::MolToSmiles(temp_mol)
    //             << " with properties: Z=" << at.atomic_number
    //             << " H=" << h_count << " q=" << charge << " s=" <<
    //             single_bonds
    //             << " d=" << double_bonds << " t=" << triple_bonds
    //             << " a=" << aromatic_bonds << " arom=" << (arom ? 1 : 0)
    //             << std::endl;
    // }
    // --- force valence computation / cache update ---
    unsigned int failedOp = 0;
    try {
      // Use only properties sanitization for these synthetic probe molecules.
      // Full sanitization on intentionally odd intermediates (especially with
      // wildcard neighbors / aromatic probes) can trigger parser/sanitizer
      // instability on some platforms.
      //   RDKit::MolOps::sanitizeMol(temp_mol, failedOp,
      //                              RDKit::MolOps::SANITIZE_PROPERTIES);
      temp_mol.updatePropertyCache(false);
    } catch (...) {
      if (debug_level_to_int(debug_level) >=
          debug_level_to_int(DebugLevel::Verbose)) {
        std::cout << "Sanitization failed for candidate: Z=" << at.atomic_number
                  << " " << failedOp << " H=" << h_count << " q=" << charge
                  << " s=" << single_bonds << " d=" << double_bonds
                  << " t=" << triple_bonds << " a=" << aromatic_bonds
                  << " arom=" << (arom ? 1 : 0) << std::endl;
      }
      return false;
    }

    if (temp_mol.getAtomWithIdx(aidx)->getTotalNumHs() != h_count) {
      // if (debug_level_to_int(debug_level) >=
      //     debug_level_to_int(DebugLevel::Verbose)) {
      //   std::cout << "[is_valence_consistent] H count mismatch for "
      //             << RDKit::MolToSmiles(temp_mol) << ": expected " << h_count
      //             << ", got " <<
      //             temp_mol.getAtomWithIdx(aidx)->getTotalNumHs()
      //             << std::endl;
      // }
      return false;
    }

    if (debug_level_to_int(debug_level) >=
            debug_level_to_int(DebugLevel::Verbose) &&
        !temp_mol.getAtomWithIdx(aidx)->hasValenceViolation()) {
      int deg = single_bonds + double_bonds + triple_bonds + aromatic_bonds;
      std::cout << "[is_valence_consistent] " << RDKit::MolToSmiles(temp_mol)
                << " " << h_count << " " << charge << " " << single_bonds << " "
                << double_bonds << " " << triple_bonds << " " << aromatic_bonds
                << " " << temp_mol.getAtomWithIdx(0)->getNumRadicalElectrons()
                << " " << deg << std::endl;
    }

    return !temp_mol.getAtomWithIdx(aidx)->hasValenceViolation();
  }
  std::string enumerate_atom_alternatives(const PatternItem &at,
                                          const EnumerationSettings &settings,
                                          bool map_new_atoms, int &max_amap) {
    if (at.atom.atomic_number <= 0) {
      if (!at.atom.source_atom_smarts.empty()) {
        return at.atom.source_atom_smarts;
      }
      return "[*]";
    }


    struct Candidate {
      int degree;
      int h_count;
      int charge;
      int single_bonds;
      int double_bonds;
      int triple_bonds;
      int aromatic_bonds;
      int arom;
    };

    if (settings.debug_level == DebugLevel::Verbose) {
      std::cout << "[enumerate_atom_alternatives] atom index "
                << at.atom.atom_idx
                << " with SMARTS pattern: " << at.atom.source_atom_smarts
                << " extracted: " << "explicit_charge: "
                << (at.atom.explicit_charge.has_value()
                        ? std::to_string(*at.atom.explicit_charge)
                        : "None")
                << "\n explicit_H: "
                << (at.atom.explicit_H.has_value()
                        ? std::to_string(*at.atom.explicit_H)
                        : "None")
                << "\n explicit_valence: "
                << (at.atom.explicit_valence.has_value()
                        ? std::to_string(*at.atom.explicit_valence)
                        : "None")
                << "\n explicit_lower_h: "
                << (at.atom.explicit_lower_h.has_value()
                        ? std::to_string(*at.atom.explicit_lower_h)
                        : "None")
                << "\n explicit_D: "
                << (at.atom.explicit_D.has_value()
                        ? std::to_string(*at.atom.explicit_D)
                        : "None")
                << "\n explicit_X: "
                << (at.atom.explicit_X.has_value()
                        ? std::to_string(*at.atom.explicit_X)
                        : "None")
                << "\n explicit_aromatic: "
                << (at.atom.explicit_aromatic.has_value()
                        ? std::to_string(*at.atom.explicit_aromatic)
                        : "None")
                << "\n explicit_aliphatic: "
                << (at.atom.explicit_aliphatic.has_value()
                        ? std::to_string(*at.atom.explicit_aliphatic)
                        : "None")
                << "\n is_aromatic: "
                << (at.atom.is_aromatic ? "true" : "false")
                << "\n is_aliphatic: "
                << (at.atom.is_aliphatic ? "true" : "false")
                << "\n explicit_in_ring: "
                << (at.atom.explicit_in_ring.has_value()
                        ? std::to_string(*at.atom.explicit_in_ring)
                        : "None")
                << " explicit_ring_count: "
                << (at.atom.explicit_ring_count.has_value()
                        ? std::to_string(*at.atom.explicit_ring_count)
                        : "None")
                << " explicit_min_ring_size: "
                << (at.atom.explicit_min_ring_size.has_value()
                        ? std::to_string(*at.atom.explicit_min_ring_size)
                        : "None")
                << "\n excluded_charges: "
                << (at.atom.excluded_charges.empty() ? "None" : "")
                << "\n excluded_h_counts: "
                << (at.atom.excluded_h_counts.empty() ? "None" : "")
                << "\n excluded_D_values: "
                << (at.atom.excluded_D_values.empty() ? "None" : "")
                << "\n excluded_X_values: "
                << (at.atom.excluded_X_values.empty() ? "None" : "")
                << "\n excluded_aromatic: "
                << (at.atom.excluded_aromatic ? "true" : "false")
                << "\n max_valence: " << at.atom.max_valence << std::endl
                << "bond_types:" << std::endl;
      for (const auto &bt : at.atom.bond_types) {
        std::cout << "  " << bt.first << ": " << bt.second << std::endl;
      }
    }

    const bool force_h_terms =
      (settings.forced_primitives_mask & (1u << 0)) != 0u;
    const bool force_degree_terms =
      (settings.forced_primitives_mask & (1u << 1)) != 0u;
    const bool force_x_terms =
      (settings.forced_primitives_mask & (1u << 2)) != 0u;
    const bool force_charge_terms =
      (settings.forced_primitives_mask & (1u << 3)) != 0u;
    const bool force_valence_terms =
      (settings.forced_primitives_mask & (1u << 4)) != 0u;
    const bool force_ring_count_terms =
      (settings.forced_primitives_mask & (1u << 5)) != 0u;
    const bool force_ring_size_terms =
      (settings.forced_primitives_mask & (1u << 6)) != 0u;
    const bool force_ring_bond_count_terms =
      (settings.forced_primitives_mask & (1u << 7)) != 0u;

    bool use_Xs = settings.include_x_in_reserialization || force_x_terms;
    if (at.atom.explicit_X.has_value() || !at.atom.excluded_X_values.empty()) {
      use_Xs = true;
    }

    std::vector<int> hs;
    if (at.atom.explicit_H.has_value()) {
      hs.push_back(*at.atom.explicit_H);
    } else {
      const int h_min = std::max(0, settings.h_min);
      const int h_max = std::max(h_min, settings.h_max);
      for (int h = h_min; h <= h_max; ++h) {
        hs.push_back(h);
      }
    }
    hs.erase(std::remove_if(hs.begin(), hs.end(),
                            [&](int h) {
                              return std::find(
                                         at.atom.excluded_h_counts.begin(),
                                         at.atom.excluded_h_counts.end(),
                                         h) != at.atom.excluded_h_counts.end();
                            }),
             hs.end());

    // Per-element charge range table.
    // Maps atomic number → {min_charge, max_charge}.
    // Elements not in the table fall back to settings.charge_min/max.
    static const std::unordered_map<int, std::pair<int, int>>
        element_charge_range = {
            {1,  {0, 1}},   // H:  0, +1
            {3,  {0, 1}},   // Li: 0, +1
            {5,  {-1, 0}},  // B:  -1, 0
            {6,  {-1, 1}},  // C:  -1, 0, +1
            {7,  {0, 1}},  // N:  -1, 0, +1
            {8,  {-1, 0}},  // O:  -1, 0
            {9,  {-1, 0}},  // F:  -1, 0
            {11, {0, 1}},   // Na: 0, +1
            {12, {0, 2}},   // Mg: 0, +2
            {14, {-1, 1}},  // Si: -1, 0, +1
            {15, {-1, 1}},  // P:  -1, 0, +1
            {16, {-1, 1}},  // S:  -1, 0, +1
            {17, {-1, 0}},  // Cl: -1, 0
            {19, {0, 1}},   // K:  0, +1
            {20, {0, 2}},   // Ca: 0, +2
            {26, {0, 3}},   // Fe: 0..+3
            {29, {0, 2}},   // Cu: 0..+2
            {30, {0, 2}},   // Zn: 0..+2
            {34, {-1, 1}},  // Se: -1, 0, +1
            {35, {-1, 0}},  // Br: -1, 0
            {53, {-1, 1}},  // I:  -1, 0, +1
        };

    std::vector<int> charges;
    if (at.atom.explicit_charge.has_value()) {
      charges.push_back(*at.atom.explicit_charge);
    } else {
      int charge_min = std::min(settings.charge_min, settings.charge_max);
      int charge_max = std::max(settings.charge_min, settings.charge_max);

      // Override from per-element table if present.
      auto ecr_it = element_charge_range.find(at.atom.atomic_number);
      if (ecr_it != element_charge_range.end()) {
        charge_min = std::max(charge_min, ecr_it->second.first);
        charge_max = std::min(charge_max, ecr_it->second.second);
      }

      for (int c = charge_min; c <= charge_max; ++c) {
        charges.push_back(c);
      }
    }
    charges.erase(
        std::remove_if(charges.begin(), charges.end(),
                       [&](int c) {
                         return std::find(at.atom.excluded_charges.begin(),
                                          at.atom.excluded_charges.end(),
                                          c) != at.atom.excluded_charges.end();
                       }),
        charges.end());

    int lower_h_count = 0;
    if (at.atom.explicit_lower_h.has_value()) {
      lower_h_count = *at.atom.explicit_lower_h;
    }

    std::vector<int> degree_values;
    if (at.atom.explicit_D.has_value()) {
      degree_values.push_back(std::max(*at.atom.explicit_D, lower_h_count));
    } else {
      for (int d = at.atom.min_bonds + lower_h_count; d <= 5; ++d) {
        degree_values.push_back(d);
      }
    }
    degree_values.erase(
        std::remove_if(degree_values.begin(), degree_values.end(),
                       [&](int d) {
                         return std::find(at.atom.excluded_D_values.begin(),
                                          at.atom.excluded_D_values.end(), d) !=
                                    at.atom.excluded_D_values.end() ||
                                std::find(at.atom.excluded_X_values.begin(),
                                          at.atom.excluded_X_values.end(),
                                          d) != at.atom.excluded_X_values.end();
                       }),
        degree_values.end());

    std::vector<int> single_bond_values;
    int min_num_single_bonds = 0;
    if (at.atom.is_aliphatic && !at.atom.is_aromatic) {
      min_num_single_bonds = at.atom.num_single_bonds;
    }
    for (int d = min_num_single_bonds; d <= 4 - lower_h_count; ++d) {
      single_bond_values.push_back(d);
    }

    std::vector<int> double_bond_values;
    for (int d = at.atom.num_double_bonds; d <= 2; ++d) {
      double_bond_values.push_back(d);
    }

    std::vector<int> triple_bond_values;
    for (int d = at.atom.num_triple_bonds; d <= 1; ++d) {
      triple_bond_values.push_back(d);
    }
    std::vector<int> aromatic_bond_values;
    std::vector<int> arom;
    // arom.push_back(false);
    if (at.atom.is_aromatic) {
      // if (!at.atom.excluded_aromatic) {
      // }
      arom.push_back(1);
      // Aromatic atoms in valid aromatic systems have at least two aromatic
      // neighbors. Without this floor, cases like "cC" can incorrectly emit
      // a D1 aromatic branch alongside D3.
      const int aromatic_floor = std::max(at.atom.num_aromatic_bonds, 2);
      for (int d = aromatic_floor; d <= 3; ++d) {
        aromatic_bond_values.push_back(d);
      }
    } else {
      aromatic_bond_values.push_back(0);
    }

    if (at.atom.is_aliphatic) {
      arom.push_back(0);
    }

    // std::cout << arom.size() << std::endl;
    // for (int a : arom) {
    //   std::cout << "Arom: " << a << std::endl;
    // }

    std::vector<Candidate> candidates;
    std::set<std::tuple<int, int, int, int, int, int, int, bool>>
        seen_candidates;
    if (settings.debug_level == DebugLevel::Verbose) {
      std::cout << "[enumerate_atom_alternatives] Candidate generation ranges:" << std::endl;
      std::cout << "  H: " << hs.front() << " to " << hs.back() << std::endl;
      std::cout << "  Charge: " << charges.front() << " to " << charges.back()
                << std::endl;
      std::cout << "  Degree: " << degree_values.front() << " to "
                << degree_values.back() << std::endl;
      std::cout << "  Single bonds: " << single_bond_values.front() << " to "
                << single_bond_values.back() << std::endl;
      std::cout << "  Double bonds: " << double_bond_values.front() << " to "
                << double_bond_values.back() << std::endl;
      std::cout << "  Triple bonds: " << triple_bond_values.front() << " to "
                << triple_bond_values.back() << std::endl;
      std::cout << "  Aromatic bonds: " << aromatic_bond_values.front()
                << " to " << aromatic_bond_values.back() << std::endl;
      // std::cout << " aromatic flag: " << arom.front() << " to " <<
      // arom.back()
      //           << std::endl;
    }
    candidates.reserve(64);
    for (int h_count : hs) {
      for (int charge : charges) {
        for (int single_bonds : single_bond_values) {
          for (int double_bonds : double_bond_values) {
            for (int triple_bonds : triple_bond_values) {
              for (int aromatic_bonds : aromatic_bond_values) {
                for (int aromflag : arom) {
                  if (settings.debug_level == DebugLevel::Verbose) {
                  std::cout
                      << "Testing candidate: "
                      << " Z=" << at.atom.atomic_number << "H=" << h_count
                      << " Charge=" << charge << " Single=" << single_bonds
                      << " Double=" << double_bonds
                      << " Triple=" << triple_bonds
                      << " Aromatic=" << aromatic_bonds << " Arom=" <<
                      aromflag
                      << std::endl;
                  }
                  if (!aromflag) {
                    aromatic_bonds = 0;
                  }
                  const int degree = single_bonds + double_bonds +
                                     triple_bonds + aromatic_bonds +
                                     lower_h_count;
                  if (aromflag){
                    if (degree+h_count < 2 || degree+h_count > 3) {
                      if (settings.debug_level == DebugLevel::Verbose) {
                        std::cout << "  Skipping due to aromatic sulfur degree constraint: " << degree << std::endl;
                      }
                      continue;
                    }
                  }

                  int lower_h = 0;
                  if (at.atom.explicit_lower_h.has_value()) {
                    lower_h = *at.atom.explicit_lower_h;
                  }
                  if (at.atom.explicit_X.has_value() &&
                      (degree + h_count) != *at.atom.explicit_X) {
                    if (settings.debug_level == DebugLevel::Verbose) {
                      std::cout << "  Skipping due to explicit X constraint:" << *at.atom.explicit_X << std::endl;
                    }
                    continue;
                  }

                  if (at.atom.explicit_valence.has_value()) {
                    int valence_from_bonds = single_bonds + 2 * double_bonds +
                                             3 * triple_bonds + aromatic_bonds +
                                             lower_h_count + h_count;
                    if (valence_from_bonds != *at.atom.explicit_valence) {
                      if (settings.debug_level == DebugLevel::Verbose) {
                        std::cout
                            << "  Skipping due to explicit valence constraint: "
                            << *at.atom.explicit_valence << std::endl;
                      }
                      continue;
                    }
                  }

                  if (!is_valence_consistent(at.atom, h_count, charge,
                                             single_bonds, double_bonds,
                                             triple_bonds, aromatic_bonds,
                                             aromflag, settings.debug_level)) {
                    continue;
                  }
                  const int serialized_degree = degree;
                  if (!degree_values.empty() &&
                      std::find(degree_values.begin(), degree_values.end(),
                                serialized_degree) == degree_values.end()) {
                    if (settings.debug_level == DebugLevel::Verbose) {
                      std::cout << "  Skipping due to degree constraint: "
                                << serialized_degree << std::endl;
                    }
                    continue;
                  }
                  auto key = std::make_tuple(
                      degree, h_count, charge, single_bonds, double_bonds,
                      triple_bonds, aromatic_bonds, aromflag ? 1 : 0);
                  if (seen_candidates.insert(key).second) {
                    candidates.push_back({degree, h_count, charge, single_bonds,
                                          double_bonds, triple_bonds,
                                          aromatic_bonds, aromflag});
                  }
                }
              }
            }
          }
        }
      }
    }
    if (settings.debug_level == DebugLevel::Verbose) {
      std::cout << "Generated " << candidates.size()
                << " candidates before sorting." << std::endl;
      for (const auto &c : candidates) {
        std::cout << "Generated Candidate: degree=" << c.degree << " h_count=" << c.h_count
                  << " charge=" << c.charge << " single_bonds=" << c.single_bonds
                  << " double_bonds=" << c.double_bonds
                  << " triple_bonds=" << c.triple_bonds
                  << " aromatic_bonds=" << c.aromatic_bonds
                  << " arom=" << c.arom
                  << std::endl;
      }
    }
    if (candidates.empty()) {
      // if (settings.debug_level == DebugLevel::Verbose) {
      std::stringstream err;
      err << "No valid DoF atom alternatives could be generated for atom_idx="
          << at.atom.atom_idx << " (atomic_number=" << at.atom.atomic_number
          << ", formal_charge=" << at.atom.formal_charge
          << ", num_single_bonds=" << at.atom.num_single_bonds
          << ", num_double_bonds=" << at.atom.num_double_bonds
          << ", num_triple_bonds=" << at.atom.num_triple_bonds
          << ", num_aromatic_bonds=" << at.atom.num_aromatic_bonds
          << ", explicit_H=";

      if (at.atom.explicit_H.has_value()) {
        err << *at.atom.explicit_H;
      } else {
        err << "None";
      }
      err << ", explicit_charge=";
      if (at.atom.explicit_charge.has_value()) {
        err << *at.atom.explicit_charge;
      } else {
        err << "None";
      }
      err << ", aromatic=" << (at.atom.is_aromatic ? "true" : "false")
          << ", aliphatic=" << (at.atom.is_aliphatic ? "true" : "false")
          << ", source_atom_smarts='" << at.atom.source_atom_smarts
          << " source_smarts='" << at.atom.smarts_pattern
          << "'). Consider relaxing SMARTS constraints or enumeration bounds.";
      throw std::runtime_error(err.str());
      //   std::cout << err.str() << std::endl;
      //   return at.atom.source_atom_smarts;
      // }
    }
    const auto charge_sort_key = [](int charge) {
      if (charge >= 0) {
        return std::string("+") + std::to_string(charge);
      }
      return std::to_string(charge);
    };

    std::sort(candidates.begin(), candidates.end(),
              [charge_sort_key](const Candidate &a, const Candidate &b) {
                if (a.degree != b.degree) {
                  return a.degree < b.degree;
                }
                if (a.h_count != b.h_count) {
                  return a.h_count < b.h_count;
                }
                if (a.single_bonds != b.single_bonds) {
                  return a.single_bonds < b.single_bonds;
                }
                if (a.double_bonds != b.double_bonds) {
                  return a.double_bonds < b.double_bonds;
                }
                if (a.triple_bonds != b.triple_bonds) {
                  return a.triple_bonds < b.triple_bonds;
                }
                if (a.aromatic_bonds != b.aromatic_bonds) {
                  return a.aromatic_bonds < b.aromatic_bonds;
                }
                if (a.arom != b.arom) {
                  return a.arom < b.arom;
                }
                return charge_sort_key(a.charge) < charge_sort_key(b.charge);
              });

    const auto serialized_degree_value = [&](const Candidate &c) {
      return c.degree;
    };
    const auto all_same_degree = std::all_of(
        candidates.begin(), candidates.end(), [&](const Candidate &c) {
          return serialized_degree_value(c) ==
                 serialized_degree_value(candidates.front());
        });
    const auto all_same_h = std::all_of(
        candidates.begin(), candidates.end(), [&](const Candidate &c) {
          return c.h_count == candidates.front().h_count;
        });
    const auto all_same_charge = std::all_of(
        candidates.begin(), candidates.end(), [&](const Candidate &c) {
          return c.charge == candidates.front().charge;
        });
    const auto candidate_x_value = [](const Candidate &c) {
      return c.degree + c.h_count;
    };
    const auto all_same_x = std::all_of(
        candidates.begin(), candidates.end(), [&](const Candidate &c) {
          return candidate_x_value(c) == candidate_x_value(candidates.front());
        });
    // v (valence) = single + 2*double + 3*triple + aromatic + h_count
    const auto candidate_v_value = [](const Candidate &c) {
      return c.single_bonds + 2 * c.double_bonds + 3 * c.triple_bonds +
             c.aromatic_bonds + c.h_count;
    };
    const auto all_same_v = std::all_of(
        candidates.begin(), candidates.end(), [&](const Candidate &c) {
          return candidate_v_value(c) == candidate_v_value(candidates.front());
        });
    const auto all_same_arom = std::all_of(
        candidates.begin(), candidates.end(),
        [&](const Candidate &c) { return c.arom == candidates.front().arom; });
      // for (const auto &c : candidates) {
      //     std::cout << "Aromaticity varies among candidates: " << c.arom
      //               << " vs " << candidates.front().arom << std::endl;
      // }
    const auto all_same_single_bonds = std::all_of(
        candidates.begin(), candidates.end(), [&](const Candidate &c) {
          return c.single_bonds == candidates.front().single_bonds;
        });
    const auto all_same_double_bonds = std::all_of(
        candidates.begin(), candidates.end(), [&](const Candidate &c) {
          return c.double_bonds == candidates.front().double_bonds;
        });
    const auto all_same_triple_bonds = std::all_of(
        candidates.begin(), candidates.end(), [&](const Candidate &c) {
          return c.triple_bonds == candidates.front().triple_bonds;
        });
    const auto all_same_aromatic_bonds = std::all_of(
        candidates.begin(), candidates.end(), [&](const Candidate &c) {
          return c.aromatic_bonds == candidates.front().aromatic_bonds;
        });

    // Preserve semantic fidelity: when H is explicit in the input and D was
    // not explicitly requested, inferring and emitting D can over-constrain
    // isotope-bearing matches (e.g., [2H]X patterns).
    const bool include_degree_terms =
      at.atom.explicit_D.has_value() || !at.atom.explicit_H.has_value();
    const bool include_x_terms = use_Xs;
    const bool include_h_terms = true;
    const bool include_charge_terms = true;
    const bool include_arom_terms = true;
    const bool include_valence_terms = force_valence_terms;
    const bool include_ring_count_terms = force_ring_count_terms;
    const bool include_ring_size_terms = force_ring_size_terms;
    const bool include_ring_bond_count_terms = force_ring_bond_count_terms;

    // Negation shortcut: when a dimension varies only because of exclusion
    // constraints (no positive value was specified in the original SMARTS),
    // emit the negation(s) as shared prefix constraints instead of
    // enumerating all valid alternatives in OR branches.  This avoids
    // combinatorial explosion for patterns like [#16&!H0].
    const bool h_via_negation = !force_h_terms && !all_same_h &&
        !at.atom.explicit_H.has_value() &&
        !at.atom.excluded_h_counts.empty();
    const bool charge_via_negation = !force_charge_terms && !all_same_charge &&
        !at.atom.explicit_charge.has_value() &&
        !at.atom.excluded_charges.empty();
    const bool degree_via_negation = !force_degree_terms && !all_same_degree &&
        !at.atom.explicit_D.has_value() &&
        !at.atom.excluded_D_values.empty();

    // Unconstrained shortcut: when a dimension was never specified or
    // excluded in the original SMARTS, it adds no discriminating power
    // beyond what the element + other constraints already provide.
    // Omit it from the output to avoid redundant enumeration.
    const bool h_unconstrained = !force_h_terms && !all_same_h &&
        !at.atom.explicit_H.has_value() &&
        at.atom.excluded_h_counts.empty();
    const bool charge_unconstrained = !force_charge_terms && !all_same_charge &&
        !at.atom.explicit_charge.has_value() &&
        at.atom.excluded_charges.empty();
    const bool degree_unconstrained = !force_degree_terms && !all_same_degree &&
        !at.atom.explicit_D.has_value() &&
        at.atom.excluded_D_values.empty();

    // When candidates include both aromatic and aliphatic forms and the
    // atomic number is already in the prefix, the A/a terms are redundant
    // — #N already covers both.  Drop them to avoid pointless OR branches.
    // const bool arom_both_present = !all_same_arom;

    const bool no_varying_primitive_terms =
        (!include_degree_terms || all_same_degree || degree_via_negation || degree_unconstrained) &&
        (!include_x_terms || all_same_x) &&
        (!include_h_terms || all_same_h || h_via_negation || h_unconstrained) &&
        (!include_charge_terms || all_same_charge || charge_via_negation || charge_unconstrained) &&
        (!include_arom_terms || all_same_arom ) &&
        (!include_valence_terms || all_same_v);

    const bool include_bond_constraints = settings.enumerate_bond_order;
    const bool no_varying_bonds =
      !include_bond_constraints ||
      (all_same_single_bonds && all_same_double_bonds &&
       all_same_triple_bonds && all_same_aromatic_bonds);
    const bool no_varying_terms =
      no_varying_primitive_terms && no_varying_bonds;
    std::stringstream out;

    const auto append_added_bond_constraints =
        [&](std::stringstream &ss, int add_single, int add_double,
            int add_triple, int add_aromatic, int &next_atom_map_idx,
            bool map_new_atoms) {
          const auto emit = [&](const std::string &bond_token,
                                const std::string &atom_token, int count) {
            for (int n = 0; n < count; ++n) {
              const int map_idx = next_atom_map_idx++;
              const std::string token = bond_token + "[" + atom_token + ":" +
                                        std::to_string(map_idx) + "]";
              ss << "(" << token << ")";
            }
          };

          emit("-", "*", add_single);
          emit("=", "*", add_double);
          emit("#", "*", add_triple);
          emit(":", "a", add_aromatic);
        };

    const auto ring_constraint_token = [&]() -> std::string {
      if (at.atom.explicit_in_ring.has_value() &&
          !(*at.atom.explicit_in_ring)) {
        return "!R";
      }
      if (at.atom.explicit_ring_count.has_value() &&
          *at.atom.explicit_ring_count >= 0) {
        return "R" + std::to_string(*at.atom.explicit_ring_count);
      }
      if (at.atom.explicit_min_ring_size.has_value() &&
          *at.atom.explicit_min_ring_size > 0) {
        return "R" + std::to_string(*at.atom.explicit_min_ring_size);
      }
      if (at.atom.explicit_in_ring.has_value() && *at.atom.explicit_in_ring) {
        return "R";
      }
      return "";
    }();

    // If the atom is attached to a double or triple bond it must be
    // aliphatic, so emit the element symbol instead of the generic #N
    // token.  This ensures [#6]=[#8] and [C]=[O] produce the same output.
    const bool has_multi_bond =
        at.atom.num_double_bonds > 0 || at.atom.num_triple_bonds > 0;

    static const std::unordered_map<int, std::string> aliphatic_symbol = {
        {5, "B"},   {6, "C"},   {7, "N"},   {8, "O"},   {9, "F"},
        {14, "Si"}, {15, "P"},  {16, "S"},  {17, "Cl"}, {33, "As"},
        {34, "Se"}, {35, "Br"}, {53, "I"},
    };

    std::string atom_prefix;
    if (has_multi_bond && !at.atom.is_aromatic) {
      auto sym_it = aliphatic_symbol.find(at.atom.atomic_number);
      if (sym_it != aliphatic_symbol.end()) {
        atom_prefix = sym_it->second;
      } else {
        atom_prefix = "#" + std::to_string(at.atom.atomic_number);
      }
    } else {
      atom_prefix = "#" + std::to_string(at.atom.atomic_number);
    }

    // When the prefix is an element symbol (e.g. "N", "C") rather than
    // the generic "#N" form, it already implies aliphatic (uppercase) or
    // aromatic (lowercase).  Track this so we can suppress a redundant ;A / ;a.
    const bool prefix_implies_aliphatic =
        !atom_prefix.empty() && atom_prefix[0] != '#' &&
        std::isupper(static_cast<unsigned char>(atom_prefix[0]));
    const bool prefix_implies_aromatic =
        !atom_prefix.empty() && atom_prefix[0] != '#' &&
        std::islower(static_cast<unsigned char>(atom_prefix[0]));

    out << "[" << atom_prefix;
    if (!ring_constraint_token.empty()) {
      out << ";" << ring_constraint_token;
    }

    // if (settings.debug_level == DebugLevel::Verbose) {
    //   std::cout << "[enumerate_atom_alternatives] ;#alt=" << candidates.size()
    //             << std::endl;
    //   for (const auto &c : candidates) {
    //     std::cout << "[enumerate_atom_alternatives] Candidate: D=" << c.degree
    //               << " H=" << c.h_count << " Charge=" << c.charge
    //               << " Single=" << c.single_bonds
    //               << " Double=" << c.double_bonds
    //               << " Triple=" << c.triple_bonds
    //               << " Aromatic=" << c.aromatic_bonds << " Arom=" << c.arom
    //               << std::endl;
    //   }
    // }

    if (h_via_negation) {
      for (int excl : at.atom.excluded_h_counts) {
        out << ";!H" << excl;
      }
    } else if (include_h_terms && all_same_h) {
      out << ";H" << candidates.front().h_count;
    }
    if (degree_via_negation) {
      for (int excl : at.atom.excluded_D_values) {
        out << ";!D" << excl;
      }
    } else if (include_degree_terms && all_same_degree) {
      out << ";D" << serialized_degree_value(candidates.front());
    }
    if (include_x_terms && all_same_x) {
      out << ";X" << candidate_x_value(candidates.front());
    }
    if (charge_via_negation) {
      for (int excl : at.atom.excluded_charges) {
        out << ";!" << charge_token(excl);
      }
    } else if (include_charge_terms && all_same_charge) {
      out << ";" << charge_token(candidates.front().charge);
    }

    // Valence (v): calculated from enumerated bonds + Hs
    if (include_valence_terms && all_same_v) {
      out << ";v" << candidate_v_value(candidates.front());
    }

    // Ring count (R): from input atom type
    if (include_ring_count_terms && ring_constraint_token.empty()) {
      if (at.atom.explicit_ring_count.has_value()) {
        out << ";R" << *at.atom.explicit_ring_count;
      } else if (at.atom.explicit_in_ring.has_value()) {
        if (*at.atom.explicit_in_ring) {
          out << ";R";
        } else {
          out << ";R0";
        }
      }
    }

    // Ring size (r): smallest ring size from input
    if (include_ring_size_terms) {
      if (at.atom.explicit_min_ring_size.has_value() &&
          *at.atom.explicit_min_ring_size > 0) {
        out << ";r" << *at.atom.explicit_min_ring_size;
      } else if (at.atom.ring_size > 0) {
        out << ";r" << at.atom.ring_size;
      }
    }

    // Ring bond count (x): from input atom type
    if (include_ring_bond_count_terms) {
      out << ";x" << at.atom.num_ring_bonds;
    }

    if (lower_h_count > 0) {
      out << ";h" << lower_h_count;
    }

    if (include_arom_terms && all_same_arom) {
      if (candidates.front().arom) {
        if (!prefix_implies_aromatic) {
          out << ";a";
        }
      } else {
        if (!prefix_implies_aliphatic) {
          out << ";A";
        }
      }
    }

    if (include_bond_constraints && no_varying_bonds) {
      const int add_single =
          candidates.front().single_bonds - at.atom.num_single_bonds;
      const int add_double =
          candidates.front().double_bonds - at.atom.num_double_bonds;
      const int add_triple =
          candidates.front().triple_bonds - at.atom.num_triple_bonds;
      const int add_aromatic =
          candidates.front().aromatic_bonds - at.atom.num_aromatic_bonds;
      if (add_single > 0 || add_double > 0 || add_triple > 0 ||
          add_aromatic > 0) {
        out << ";$(" << "[#" << std::to_string(at.atom.atomic_number) << "]";
        append_added_bond_constraints(out, add_single, add_double, add_triple,
                                      add_aromatic, max_amap, map_new_atoms);
        out << ")";
      }
    }

    if (settings.debug_level == DebugLevel::Verbose) {
      std::cout << "[enumerate_atom_alternatives] join flags "
                << no_varying_bonds << " " << no_varying_terms << " "
                << " " << all_same_aromatic_bonds << " "
                << all_same_single_bonds << " " << all_same_double_bonds << " "
                << all_same_triple_bonds << " " << all_same_arom << " " << std::endl;
    }

    if (no_varying_terms) {
      out << "]";
      return out.str();
    }

    out << ";";
    int next_atom_map_idx_hold = max_amap;

    for (size_t i = 0; i < candidates.size(); ++i) {
      bool wrote_any = false;
      if (include_degree_terms && !all_same_degree && !degree_via_negation && !degree_unconstrained) {
        out << "D" << serialized_degree_value(candidates[i]);
        wrote_any = true;
      }
      if (include_x_terms && !all_same_x) {
        if (wrote_any) {
          out << "&";
        }
        out << "X" << candidate_x_value(candidates[i]);
        wrote_any = true;
      }
      if (include_h_terms && !all_same_h && !h_via_negation && !h_unconstrained) {
        if (wrote_any) {
          out << "&";
        }
        out << "H" << candidates[i].h_count;
        wrote_any = true;
      }
      if (include_charge_terms && !all_same_charge && !charge_via_negation && !charge_unconstrained) {
        if (wrote_any) {
          out << "&";
        }
        out << charge_token(candidates[i].charge);
        wrote_any = true;
      }

      if (include_valence_terms && !all_same_v) {
        if (wrote_any) {
          out << "&";
        }
        out << "v" << candidate_v_value(candidates[i]);
        wrote_any = true;
      }

      if (include_arom_terms && !all_same_arom) {
        const bool skip_arom_token =
            (candidates[i].arom && prefix_implies_aromatic) ||
            (!candidates[i].arom && prefix_implies_aliphatic);
        if (!skip_arom_token) {
          if (wrote_any) {
            out << "&";
          }
          out << (candidates[i].arom ? "a" : "A");
          wrote_any = true;
        }
      }

      if (include_bond_constraints && !no_varying_bonds) {
        const int add_single =
            candidates[i].single_bonds - at.atom.num_single_bonds;
        const int add_double =
            candidates[i].double_bonds - at.atom.num_double_bonds;
        const int add_triple =
            candidates[i].triple_bonds - at.atom.num_triple_bonds;
        const int add_aromatic =
            candidates[i].aromatic_bonds - at.atom.num_aromatic_bonds;
        if (add_single > 0 || add_double > 0 || add_triple > 0 ||
            add_aromatic > 0) {
          if (wrote_any) {
            out << "&";
          }
          out << "$(" << "[#" << std::to_string(at.atom.atomic_number) << "]";
          append_added_bond_constraints(out, add_single, add_double, add_triple,
                                        add_aromatic, max_amap, map_new_atoms);
          out << ")";
          wrote_any = true;
        }
        // max_amap = next_atom_map_idx_hold;
      }
      if (!wrote_any) {
        if (settings.debug_level == DebugLevel::Trace) {
          std::cout
              << "Warning: candidate " << i
              << " has no varying terms but was not caught by all_same checks."
              << std::endl;
        }
        // out << "D" << candidates[i].degree << "&H" << candidates[i].h_count
        //     << "&" << charge_token(candidates[i].charge);
      }

      if (i + 1 < candidates.size()) {
        if (wrote_any) {
          out << ",";
        }
      }
    }
    out << "]";

    return out.str();
  }

  /**
   * Process a molecule and extract all atom types
   */
  std::vector<AtomType> process_molecule(RDKit::ROMol *mol,
                                         bool use_query_constraints = false) {
    if (!mol) {
      throw std::runtime_error("Invalid molecule");
    }

    // Assign CIP (R/S) stereochemistry codes to atoms
    // Parameters: cleanIt=true, force=true, flagPossibleStereoCenters=false
    RDKit::MolOps::assignStereochemistry(*mol, true, true, true);

    std::vector<AtomType> atom_types;
    atom_types.reserve(mol->getNumAtoms());

    for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
      auto *atom = mol->getAtomWithIdx(i);
      if (use_query_constraints) {
        atom_types.push_back(extract_smarts_atom_type(atom, mol, i));
      } else {
        atom_types.push_back(extract_atom_type(atom, mol, i));
      }
    }

    return atom_types;
  }
};

// AtomTyper implementation
AtomTyper::AtomTyper() : pimpl(std::make_unique<Impl>()) {}

AtomTyper::~AtomTyper() = default;

std::vector<PatternItem> AtomTyper::type_atoms_from_smiles(
  const std::string &smiles) {
  if (smiles.empty()) {
    throw std::runtime_error("Empty SMILES string");
  }

  // Parse SMILES
  std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(smiles));
  if (!mol) {
    throw std::runtime_error("Failed to parse SMILES: " + smiles);
  }

  // Apply canonicalization if enabled
  if (pimpl->use_canonical) {
    // Use RDKit's canonical SMILES
    std::string canonical = RDKit::MolToSmiles(*mol);
    mol.reset(RDKit::SmilesToMol(canonical));
  }

  // Process molecule and extract atom/bond pattern items
  auto atom_types = pimpl->process_molecule(mol.get(), false);
  std::vector<PatternItem> items;
  if (atom_types.empty()) {
    return items;
  }

  items.reserve(atom_types.size() * 2 - 1);
  for (size_t i = 0; i < atom_types.size(); ++i) {
    PatternItem atom_item;
    atom_item.kind = PatternItemKind::Atom;
    atom_item.atom = atom_types[i];
    items.push_back(atom_item);

    if (i + 1 < atom_types.size()) {
      if (auto *bond = mol->getBondBetweenAtoms(static_cast<unsigned int>(i),
                                                static_cast<unsigned int>(i + 1))) {
        PatternItem bond_item;
        bond_item.kind = PatternItemKind::Bond;
        bond_item.bond = pimpl->make_bond_type(bond);
        items.push_back(bond_item);
      }
    }
  }

  return items;
}


std::string AtomTyper::type_atoms_from_smarts(const std::string &smarts) {
  int max_amap = 0;
  return type_atoms_from_smarts(smarts, false, max_amap, false, false);
}

std::string AtomTyper::type_atoms_from_smarts(const std::string &smarts,
                                              DebugLevel debug_level) {
  int max_amap = 0;
  const bool verbose = debug_level != DebugLevel::Off;
  return type_atoms_from_smarts(smarts, 0, 4, -1, 1, verbose, debug_level,
                                false, max_amap, false, true);
}

std::string AtomTyper::type_atoms_from_smarts(
    const std::string &smarts, bool map_new_atoms, int &max_amap, bool verbose,
    bool include_x_in_reserialization, bool enumerate_bond_order,
    unsigned int forced_primitives_mask) {
  return type_atoms_from_smarts(smarts, 0, 4, -1, 1, verbose,
                                pimpl->default_debug_level, map_new_atoms,
                                max_amap, include_x_in_reserialization,
                                enumerate_bond_order,
                                forced_primitives_mask);
}

std::string AtomTyper::type_atoms_from_smarts(
    const std::string &smarts, int h_min, int h_max, int charge_min,
    int charge_max, bool map_new_atoms, int &max_amap,
    bool include_x_in_reserialization, bool enumerate_bond_order,
    unsigned int forced_primitives_mask) {
  return type_atoms_from_smarts(smarts, h_min, h_max, charge_min, charge_max,
                                false, pimpl->default_debug_level,
                                map_new_atoms, max_amap,
                                include_x_in_reserialization,
                                enumerate_bond_order,
                                forced_primitives_mask);
}

std::string AtomTyper::type_atoms_from_smarts(
    const std::string &smarts, int h_min, int h_max, int charge_min,
    int charge_max, bool verbose, DebugLevel debug_level, bool map_new_atoms,
    int &max_amap, bool include_x_in_reserialization,
  bool enumerate_bond_order, unsigned int forced_primitives_mask) {
  ZoneScopedN("AtomTyper::enumerate_dof_smarts");
  if (smarts.empty()) {
    std::cout << "[type_atoms_from_smarts] Error: Empty SMARTS string provided."
              << std::endl;
    throw std::runtime_error("Empty SMARTS string");
  }

  std::unique_ptr<RDKit::ROMol> mol(RDKit::SmartsToMol(smarts));
  if (!mol) {
    std::cout << "[type_atoms_from_smarts] Failed to parse SMARTS: " << smarts
              << std::endl;
    throw std::runtime_error("Failed to parse SMARTS: " + smarts);
  }

  // Compute implicit valence / property cache BEFORE process_molecule so that
  // assignStereochemistry and extract_atom_type can safely call
  // getTotalNumHs() on query atoms (e.g. [#7&h1] has df_noImplicit=false).
  try {
    mol->updatePropertyCache(false);
  } catch (...) {
    // Best-effort: SMARTS molecules may not have enough info for all atoms.
  }

  AtomTyper::Impl::EnumerationSettings settings;
  settings.h_min = std::min(h_min, h_max);
  settings.h_max = std::max(h_min, h_max);
  settings.charge_min = std::min(charge_min, charge_max);
  settings.charge_max = std::max(charge_min, charge_max);
  settings.debug_level = verbose ? debug_level : DebugLevel::Off;
  settings.include_x_in_reserialization = include_x_in_reserialization;
  settings.enumerate_bond_order = enumerate_bond_order;
  settings.forced_primitives_mask = forced_primitives_mask;

  if (debug_level == DebugLevel::Verbose) {
    verbose = true;
    settings.debug_level = DebugLevel::Verbose;
  }
  if (verbose) {
    settings.debug_level = DebugLevel::Verbose;
  }

  auto atom_types = pimpl->process_molecule(mol.get(), true);
  RDKit::RWMol rewritten(*mol);
  std::vector<std::unique_ptr<RDKit::ROMol>> query_mol_owners;
  query_mol_owners.reserve(rewritten.getNumAtoms());

  try {
    rewritten.updatePropertyCache(false);
    RDKit::MolOps::findSSSR(rewritten);
  } catch (...) {
    // Ring perception may fail for partial query graphs; keep best-effort
    // behavior and fall back to atom-level ring flags.
  }

  if (verbose) {
    std::cout << "[type_atoms_from_smarts] Original SMARTS: " << smarts
              << std::endl;
    for (size_t i = 0; i < atom_types.size(); ++i) {
      const auto &at = atom_types[i];
      std::cout << "[type_atoms_from_smarts] Atom idx " << i
                << ": Z=" << at.atomic_number << " H=" << at.num_hydrogens
                << " q=" << at.formal_charge << " s=" << at.num_single_bonds
                << " d=" << at.num_double_bonds << " t=" << at.num_triple_bonds
                << " a=" << at.num_aromatic_bonds
                << " arom=" << (at.is_aromatic ? 1 : 0)
                << " alip=" << (at.is_aliphatic ? 1 : 0) << std::endl;
    }
  }

  max_amap =
      max_amap + 1;  // Start new maps from max existing + 1 to avoid conflicts

  const auto is_single_primitive_atom_pattern = [](std::string token) {
    if (token.size() < 2 || token.front() != '[' || token.back() != ']') {
      return false;
    }

    token = token.substr(1, token.size() - 2);

    const auto colon = token.rfind(':');
    if (colon != std::string::npos && colon + 1 < token.size() &&
        std::all_of(token.begin() + static_cast<std::ptrdiff_t>(colon + 1),
                    token.end(), [](char c) {
                      return std::isdigit(static_cast<unsigned char>(c));
                    })) {
      token = token.substr(0, colon);
    }

    if (token == "*" || token == "a" || token == "A") {
      return true;
    }
    if (token == "!#1") {
      return true;
    }

    if (token.size() >= 2 && token[0] == '#') {
      return std::all_of(token.begin() + 1, token.end(), [](char c) {
        return std::isdigit(static_cast<unsigned char>(c));
      });
    }

    const auto is_letter = [](char c) {
      return std::isalpha(static_cast<unsigned char>(c)) != 0;
    };

    if (std::all_of(token.begin(), token.end(), is_letter)) {
      if (token.size() == 1) {
        return true;
      }
      if (token.size() == 2) {
        // Aliphatic element symbols like Cl, Br, Si ...
        if (std::isupper(static_cast<unsigned char>(token[0])) &&
            std::islower(static_cast<unsigned char>(token[1]))) {
          return true;
        }
        // Aromatic two-letter symbols in SMARTS: as, se
        if ((token == "as") || (token == "se")) {
          return true;
        }
      }
    }

    return false;
  };

  const auto is_atomic_anchor_with_recursive_only =
      [](const std::string &token, const AtomType &at) {
        if (token.find("$(") == std::string::npos) {
          return false;
        }
        if (!at.explicit_atomic_num.has_value()) {
          return false;
        }

        // Keep this narrow: only skip DoF typing when the atom is an
        // atomic-number anchor with recursive constraints and no explicit
        // local DoF/typing primitives to preserve.
        return !at.explicit_charge.has_value() && !at.explicit_H.has_value() &&
               !at.explicit_lower_h.has_value() && !at.explicit_D.has_value() &&
               !at.explicit_X.has_value() && !at.explicit_valence.has_value() &&
               !at.explicit_hybridization.has_value() &&
               !at.explicit_aromatic.has_value() &&
               !at.explicit_aliphatic.has_value() &&
               !at.explicit_in_ring.has_value() &&
               !at.explicit_ring_count.has_value() &&
               !at.explicit_min_ring_size.has_value() &&
               at.excluded_charges.empty() && at.excluded_h_counts.empty() &&
               at.excluded_D_values.empty() && at.excluded_X_values.empty() &&
               !at.excluded_aromatic && !at.excluded_aliphatic;
      };
  if (verbose) {
    std::cout << "[type_atoms_from_smarts] before atom enumeration, SMARTS: "
              << RDKit::MolToSmarts(rewritten) << std::endl;
  }

  const bool force_extract_any = settings.forced_primitives_mask != 0u;

  for (unsigned int i = 0; i < rewritten.getNumAtoms() && i < atom_types.size();
       ++i) {
    const std::string source_atom_smarts =
        RDKit::SmartsWrite::GetAtomSmarts(rewritten.getAtomWithIdx(i));

    if (!force_extract_any &&
      is_atomic_anchor_with_recursive_only(source_atom_smarts,
                         atom_types[i])) {
      if (verbose) {
        std::cout << "[type_atoms_from_smarts] Skipping DoF enumeration for atom idx "
                  << i
                  << " (atomic-number anchor with recursive constraint): "
                  << source_atom_smarts << std::endl;
      }
      continue;
    }
    
    if (is_single_primitive_atom_pattern(source_atom_smarts)) {
      // If the atom uses #N notation and is attached to a double or triple
      // bond, we may convert it to the aliphatic element symbol for
      // normalization. Keep #6 as atomic-number form to avoid forcing carbon
      // to aliphatic-only in mixed/aromatic-capable contexts.
      static const std::unordered_map<int, std::string> aliphatic_sym = {
          {5, "B"},   {6, "C"},   {7, "N"},   {8, "O"},   {9, "F"},
          {14, "Si"}, {15, "P"},  {16, "S"},  {17, "Cl"}, {33, "As"},
          {34, "Se"}, {35, "Br"}, {53, "I"},
      };

      // Strip brackets to get the bare token (e.g. "#6" from "[#6]")
      std::string token = source_atom_smarts;
      if (token.size() >= 2 && token.front() == '[' && token.back() == ']') {
        token = token.substr(1, token.size() - 2);
      }
      // Strip trailing atom map (e.g. ":1")
      auto colon = token.rfind(':');
      if (colon != std::string::npos) {
        token = token.substr(0, colon);
      }

      if (token.size() >= 2 && token[0] == '#') {
        bool has_multi_bond = false;
        const auto *rw_atom = rewritten.getAtomWithIdx(i);
        for (const auto &nbr : rewritten.atomNeighbors(rw_atom)) {
          const auto *bond =
              rewritten.getBondBetweenAtoms(i, nbr->getIdx());
          if (!bond) continue;
          const std::string bs = RDKit::SmartsWrite::GetBondSmarts(
              bond, static_cast<int>(bond->getBeginAtomIdx()));
          if (bs.find('=') != std::string::npos ||
              bs.find('#') != std::string::npos) {
            has_multi_bond = true;
            break;
          }
          auto bt = bond->getBondType();
          if (bt == RDKit::Bond::DOUBLE || bt == RDKit::Bond::TRIPLE) {
            has_multi_bond = true;
            break;
          }
        }

        if (has_multi_bond) {
          int anum = 0;
          try {
            anum = std::stoi(token.substr(1));
          } catch (...) {}
          if (anum == 6) {
            // Preserve [#6] on multiple bonds instead of forcing [C].
            // This avoids narrowing matches for carbonyl-like motifs.
            if (verbose) {
              std::cout << "[type_atoms_from_smarts] Preserving "
                        << source_atom_smarts
                        << " (multi-bond carbon kept as atomic-number form)"
                        << std::endl;
            }
          } else {
            auto sym_it = aliphatic_sym.find(anum);
            if (sym_it != aliphatic_sym.end()) {
              // Build replacement: symbol form with original atom map
              std::string new_smarts = "[" + sym_it->second;
              // Preserve the atom map number from the original token
              auto orig_colon = source_atom_smarts.rfind(':');
              if (orig_colon != std::string::npos) {
                new_smarts += source_atom_smarts.substr(orig_colon,
                    source_atom_smarts.size() - orig_colon - 1);  // skip trailing ]
              }
              new_smarts += "]";

              // Parse and replace the atom query on the rewritten mol
              try {
                std::unique_ptr<RDKit::RWMol> frag(
                    RDKit::SmartsToMol(new_smarts));
                if (frag && frag->getNumAtoms() == 1) {
                  auto *new_atom = frag->getAtomWithIdx(0);
                  if (new_atom->hasQuery()) {
                    rewritten.getAtomWithIdx(i)->setQuery(
                        new_atom->getQuery()->copy());
                  }
                  if (verbose) {
                    std::cout << "[type_atoms_from_smarts] Converted "
                              << source_atom_smarts << " -> " << new_smarts
                              << " (multi-bond forces aliphatic)" << std::endl;
                  }
                }
              } catch (...) {}
            }
          }
        }
      }

      if (!force_extract_any) {
        if (verbose) {
          std::cout
              << "[type_atoms_from_smarts] Skipping DoF enumeration for atom idx "
              << i << " (single primitive atom pattern): "
              << source_atom_smarts << std::endl;
        }
        continue;
      }
    }

    PatternItem atom_item;
    atom_item.kind = PatternItemKind::Atom;
    atom_item.atom = atom_types[i];
    std::string atom_smarts;
    if (verbose) {
      std::cout
          << "[type_atoms_from_smarts] Enumerating DoF alternatives for atom idx "
          << i << ": " << source_atom_smarts << std::endl;
    }

    try {
      atom_smarts = pimpl->enumerate_atom_alternatives(atom_item, settings,
                                                       map_new_atoms, max_amap);
      if (verbose) {
        //   std::cout << i << ": " << std::endl;
        std::cout << "[type_atoms_from_smarts] Enumerated SMARTS for atom idx "
                  << i << ": " << source_atom_smarts << " -> " << atom_smarts
                  << std::endl;
      }
    } catch (const std::exception &e) {
      // debate whether to discard whole subpattern or just handle this
      // subsubinstance
      if (verbose) {
        std::stringstream err;
        err << "[type_atoms_from_smarts] Failed to enumerate DoF alternatives for atom idx="
            << i << " " << source_atom_smarts << ":";
        std::cout << err.str() << " Error: " << e.what() << std::endl;
      }
      continue;
      //   atom_smarts = source_atom_smarts;
      //   throw std::runtime_error(err.str());
    }
    // if (verbose) {
    //   //   std::cout << i << ": " << std::endl;
    //   std::cout << "[type_atoms_from_smarts] Enumerated SMARTS for atom idx "
    //             << i << ": " << source_atom_smarts << " -> " << atom_smarts
    //             << std::endl;
    // }
    const bool parser_risky_atom_alt =
        atom_smarts.size() > 256 &&
        atom_smarts.find("$(") != std::string::npos &&
        atom_smarts.find(',') != std::string::npos;
    if (parser_risky_atom_alt) {
      if (verbose) {
        std::cout << "[type_atoms_from_smarts] Skipping atom idx " << i
                  << " due to parser-risky generated alternative (keeping "
                     "original query): "
                  << atom_smarts << std::endl;
      }
      continue;
    }

    std::unique_ptr<RDKit::ROMol> one_atom_mol(RDKit::SmartsToMol(atom_smarts));
    if (!one_atom_mol || one_atom_mol->getNumAtoms() != 1) {
      if (verbose) {
        std::cout << "[type_atoms_from_smarts] Skipping atom idx " << i
                  << " because generated SMARTS could not be parsed as a "
                     "single atom query: "
                  << atom_smarts << std::endl;
      }
      continue;
    }
    RDKit::ROMol *one_atom_owner = one_atom_mol.get();
    auto *src_atom = one_atom_owner->getAtomWithIdx(0);
    auto *src_qatom = dynamic_cast<RDKit::QueryAtom *>(src_atom);
    if (!src_qatom || !src_qatom->getQuery()) {
      if (verbose) {
        std::cout << "[type_atoms_from_smarts] Invalid query atom for idx " << i
                  << std::endl;
      }
      continue;
    }

    // ---- Carry through RecursiveStructure queries from the original atom ----
    // The generated atom_smarts encodes only flat primitives (D, H, charge,
    // etc.).  Any $(…) recursive expressions on the original atom would be
    // lost.  Walk the original query tree, collect RecursiveStructure children,
    // and graft them onto the replacement query.
    {
      using QATOM_QUERY =
          Queries::Query<int, const RDKit::Atom *, true>;
      std::vector<QATOM_QUERY *> rec_queries;

      // Recursive collector: walks AND/OR sub-trees of the original query
      // looking for RecursiveStructure leaves.
      std::function<void(const QATOM_QUERY *)> collect_recursive;
      collect_recursive = [&](const QATOM_QUERY *node) {
        if (!node) {
          return;
        }
        if (node->getDescription() == "RecursiveStructure") {
          rec_queries.push_back(node->copy());
          return;
        }
        for (auto it = node->beginChildren(); it != node->endChildren(); ++it) {
          collect_recursive(it->get());
        }
      };

      const auto *orig_atom =
          dynamic_cast<const RDKit::QueryAtom *>(mol->getAtomWithIdx(i));
      if (orig_atom && orig_atom->getQuery()) {
        collect_recursive(orig_atom->getQuery());
      }

      if (!rec_queries.empty()) {
        // Wrap the existing replacement query + recursive queries in an AND.
        auto *combined = new RDKit::ATOM_AND_QUERY;
        combined->setDescription("AtomAnd");
        combined->addChild(
            QATOM_QUERY::CHILD_TYPE(src_qatom->getQuery()->copy()));
        for (auto *rq : rec_queries) {
          combined->addChild(QATOM_QUERY::CHILD_TYPE(rq));
        }
        src_qatom->setQuery(combined);
        if (verbose) {
          std::cout
              << "[type_atoms_from_smarts] Re-attached "
              << rec_queries.size()
              << " recursive expression(s) to atom idx " << i << ": "
              << RDKit::SmartsWrite::GetAtomSmarts(src_qatom) << std::endl;
        }
      }
    }

    try {
      rewritten.replaceAtom(i, src_atom, false, true);
      query_mol_owners.push_back(std::move(one_atom_mol));
    } catch (const std::exception &e) {
      if (verbose) {
        std::cout
            << "[type_atoms_from_smarts] Error replacing atom query for idx "
            << i << ": " << e.what() << std::endl;
      }
    }
  }

  const auto atom_is_aromatic_only = [&](unsigned int atom_idx) -> bool {
    if (atom_idx >= atom_types.size()) {
      return false;
    }
    const auto &at = atom_types[atom_idx];
    if (at.explicit_aromatic.has_value()) {
      return *at.explicit_aromatic;
    }
    if (at.explicit_aliphatic.has_value()) {
      return !(*at.explicit_aliphatic);
    }
    return at.is_aromatic && !at.is_aliphatic;
  };

  const auto atom_is_aliphatic_only = [&](unsigned int atom_idx) -> bool {
    if (atom_idx >= atom_types.size()) {
      return false;
    }
    const auto &at = atom_types[atom_idx];
    if (at.explicit_aliphatic.has_value()) {
      return *at.explicit_aliphatic;
    }
    if (at.explicit_aromatic.has_value()) {
      return !(*at.explicit_aromatic);
    }
    return at.is_aliphatic && !at.is_aromatic;
  };

  const auto set_bond_from_token = [&](RDKit::Bond *dst_bond,
                                       const std::string &bond_token) {
    if (!dst_bond || bond_token.empty()) {
      return false;
    }
    std::unique_ptr<RDKit::ROMol> probe(
        RDKit::SmartsToMol("[*]" + bond_token + "[*]"));
    if (!probe || probe->getNumBonds() != 1) {
      return false;
    }
    const auto *src_bond = probe->getBondWithIdx(0);
    if (!src_bond || !src_bond->getQuery()) {
      return false;
    }
    dst_bond->setQuery(src_bond->getQuery()->copy());
    return true;
  };

  for (auto *bond : rewritten.bonds()) {
    if (!bond) {
      continue;
    }

    const unsigned int begin_idx = bond->getBeginAtomIdx();
    const unsigned int end_idx = bond->getEndAtomIdx();
    const bool begin_aromatic_only = atom_is_aromatic_only(begin_idx);
    const bool end_aromatic_only = atom_is_aromatic_only(end_idx);
    const bool begin_aliphatic_only = atom_is_aliphatic_only(begin_idx);
    const bool end_aliphatic_only = atom_is_aliphatic_only(end_idx);
    const bool bond_in_ring = rewritten.getRingInfo() &&
                  rewritten.getRingInfo()->numBondRings(
                    bond->getIdx()) > 0;
    const bool endpoints_in_ring =
      begin_idx < atom_types.size() && end_idx < atom_types.size() &&
      atom_types[begin_idx].is_in_ring && atom_types[end_idx].is_in_ring;

    const std::string bond_smarts =
        RDKit::SmartsWrite::GetBondSmarts(bond, begin_idx);
    // Preserve explicit any-bond (~) from input SMARTS.
    // Only truly unspecified bonds (empty token) are eligible for
    // normalization to '-' or ':'.
    const bool unspecified_bond = bond_smarts.empty();

    if (unspecified_bond) {
      if (begin_aliphatic_only && end_aliphatic_only) {
        // Do not force '-' for ring connections when the input bond was
        // unspecified; preserve ring ambiguity unless caller explicitly used
        // '-'.
        if (!(bond_in_ring || endpoints_in_ring)) {
          set_bond_from_token(bond, "-");
        }
      } else if (begin_aromatic_only && end_aromatic_only) {
        set_bond_from_token(bond, ":");
      }
      continue;
    }

    const bool aromatic_bond = bond_smarts.find(':') != std::string::npos;
    const bool aliphatic_bond = bond_smarts.find('-') != std::string::npos ||
                                bond_smarts.find('=') != std::string::npos ||
                                bond_smarts.find('#') != std::string::npos ||
                                bond_smarts.find('/') != std::string::npos ||
                                bond_smarts.find('\\') != std::string::npos;

    // Aromatic bonds are invalid only when an endpoint is explicitly forced
    // aliphatic. If aromaticity is ambiguous (e.g. [#6]), keep the variant.
    if (aromatic_bond && (begin_aliphatic_only || end_aliphatic_only)) {
      throw std::runtime_error(
          "Discarded DoF variant: aromatic bond cannot connect non-aromatic "
          "atom pair at bond idx=" +
          std::to_string(bond->getIdx()));
    }

    // if (aliphatic_bond && begin_aromatic_only && end_aromatic_only) {
    //   throw std::runtime_error(
    //       "Discarded DoF variant: aliphatic bond cannot connect two aromatic
    //       " "atoms at bond idx=" + std::to_string(bond->getIdx()));
    // }
  }

  try {
    std::string rewritten_smarts = RDKit::MolToSmarts(rewritten);

    const auto is_structurally_balanced_smarts =
        [](const std::string &s) -> bool {
      int bracket_depth = 0;
      int paren_depth = 0;
      for (const char ch : s) {
        // In SMARTS, '\\' is a directional bond token (cis/trans), not a
        // generic escape character. So we must still parse the following atom
        // token normally (e.g. ...\\[C]).
        if (ch == '[') {
          ++bracket_depth;
        } else if (ch == ']') {
          --bracket_depth;
          if (bracket_depth < 0) {
            return false;
          }
        } else if (ch == '(') {
          ++paren_depth;
        } else if (ch == ')') {
          --paren_depth;
          if (paren_depth < 0) {
            return false;
          }
        }
      }
      return bracket_depth == 0 && paren_depth == 0;
    };

    if (!is_structurally_balanced_smarts(rewritten_smarts)) {
      throw std::runtime_error(
          "Generated malformed SMARTS during enumerate_dof_smarts: " +
          rewritten_smarts);
    }

    return rewritten_smarts;

  } catch (const std::exception &e) {
    throw std::runtime_error("Error converting rewritten molecule to SMARTS: " +
                             std::string(e.what()));
    // return smarts;  // Fallback to original SMARTS on error
  }
}

std::string AtomTyper::reorder_query_tree_by_embedding(
    const std::string &smarts, const std::map<std::string, double> &embedding) {
  return query_reorder::reorder_query_tree_by_embedding_smarts(smarts,
                                                               embedding);
}

std::map<std::string, double> AtomTyper::get_default_query_embedding() const {
  // Lower score => rarer/more selective primitive.
  // Default embedding seeded from observed primitive frequencies.
  return {
      {"O", 49042},   {"!O", 252807}, {"D1", 68386},   {"!D1", 233463},
      {"H0", 126406}, {"!H0", 175443}, {"C", 117448},   {"!C", 184401},
      {"D3", 91399},  {"!D3", 210450}, {"+0", 298536},  {"!+0", 3313},
      {"c", 92059},   {"!c", 209790}, {"N", 19923},    {"!N", 281926},
      {"H1", 105310}, {"!H1", 196539}, {"D2", 133464},  {"!D2", 168385},
      {"H2", 51180},  {"!H2", 250669}, {"S", 2707},     {"!S", 299142},
      {"#11", 133},   {"!#11", 301716}, {"Cl", 1998},   {"!Cl", 299851},
      {"#7", 30837},  {"!#7", 271012}, {"a", 104226},   {"!a", 197623},
      {"n", 10914},   {"!n", 290935}, {"Br", 385},      {"!Br", 301464},
      {"H3", 18932},  {"!H3", 282917}, {"I", 318},      {"!I", 301531},
      {"-1", 1680},   {"!-1", 300169}, {"P", 1366},     {"!P", 300483},
      {"#16", 3371},  {"!#16", 298478}, {"F", 3413},    {"!F", 298436},
      {"+1", 1158},   {"!+1", 300691}, {"#19", 51},     {"!#19", 301798},
      {"D4", 7450},   {"!D4", 294399}, {"#12", 61},     {"!#12", 301788},
      {"#14", 51},    {"!#14", 301798}, {"D0", 1126},   {"!D0", 300723},
      {"#8", 49623},  {"!#8", 252226}, {"H4", 12},      {"!H4", 301837},
      {"#55", 1},     {"!#55", 301848}, {"#3", 10},     {"!#3", 301839},
      {"#34", 26},    {"!#34", 301823}, {"#50", 6},     {"!#50", 301843},
      {"B", 84},      {"!B", 301765}, {"#5", 84},      {"!#5", 301765},
      {"#30", 34},    {"!#30", 301815}, {"+3", 140},    {"!+3", 301709},
      {"s", 664},     {"!s", 301185}, {"o", 581},      {"!o", 301268},
      {"#15", 1367},  {"!#15", 300482}, {"-2", 106},    {"!-2", 301743},
      {"#25", 9},     {"!#25", 301840}, {"#24", 14},    {"!#24", 301835},
      {"#13", 55},    {"!#13", 301794}, {"D5", 9},      {"!D5", 301840},
      {"#31", 10},    {"!#31", 301839}, {"#33", 23},    {"!#33", 301826},
      {"#83", 15},    {"!#83", 301834}, {"#52", 1},     {"!#52", 301848},
      {"#81", 4},     {"!#81", 301845}, {"*", 301849},  {"#51", 11},
      {"!#51", 301838}, {"+2", 184},   {"!+2", 301665}, {"#22", 2},
      {"!#22", 301847}, {"#79", 9},    {"!#79", 301840}, {"#26", 90},
      {"!#26", 301759}, {"#47", 9},    {"!#47", 301840}, {"#29", 23},
      {"!#29", 301826}, {"!#32", 301849}, {"p", 1},     {"!p", 301848},
      {"#44", 3},     {"!#44", 301846}, {"#1", 70},     {"!#1", 301779},
      {"#89", 1},     {"!#89", 301848}, {"#80", 20},    {"!#80", 301829},
      {"#82", 1},     {"!#82", 301848}, {"D6", 14},     {"!D6", 301835},
      {"#58", 2},     {"!#58", 301847}, {"#74", 4},     {"!#74", 301845},
      {"#20", 34},    {"!#20", 301815}, {"#28", 2},     {"!#28", 301847},
      {"#49", 7},     {"!#49", 301842}, {"#40", 4},     {"!#40", 301845},
      {"#54", 4},     {"!#54", 301845}, {"#73", 7},     {"!#73", 301842},
      {"+5", 5},      {"!+5", 301844}, {"#56", 2},      {"!#56", 301847},
      {"#78", 11},    {"!#78", 301838}, {"#48", 1},     {"!#48", 301848},
      {"#46", 2},     {"!#46", 301847}, {"#42", 24},    {"!#42", 301825},
      {"#23", 20},    {"!#23", 301829}, {"!#77", 301849}, {"#53", 318},
      {"!#53", 301531}, {"#27", 16},   {"!#27", 301833}, {"#57", 3},
      {"!#57", 301846}, {"!D8", 301849}, {"se", 7},     {"!se", 301842},
      {"!#45", 301849}, {"#62", 3},    {"!#62", 301846}, {"+4", 22},
      {"!+4", 301827}, {"#39", 2},     {"!#39", 301847}, {"!#72", 301849},
      {"-3", 5},      {"!-3", 301844}, {"-4", 4},       {"!-4", 301845},
      {"+6", 8},      {"!+6", 301841}, {"#17", 1998},   {"!#17", 299851},
      {"#35", 385},   {"!#35", 301464}, {"#9", 3413},   {"!#9", 298436},
      {"#76", 1},     {"!#76", 301848}, {"#75", 2},     {"!#75", 301847},
      {"!#93", 301849}, {"!#92", 301849}, {"#41", 1},   {"!#41", 301848},
      {"#64", 13},    {"!#64", 301836}, {"!#95", 301849}, {"!#94", 301849},
      {"!H7", 301849}, {"!H6", 301849}, {"!#63", 301849}, {"!#65", 301849},
      {"#38", 6},     {"!#38", 301843}, {"#43", 14},    {"!#43", 301835},
      {"!H5", 301849},
  };
}

bool AtomTyper::inspect_tautomer(const std::string &smarts,
                                 bool verbose) const {
  if (smarts.empty()) {
    return false;
  }

  std::unique_ptr<RDKit::ROMol> mol(RDKit::SmartsToMol(smarts));
  if (!mol) {
    return false;
  }

  struct TautomerPattern {
    const char *name;
    const char *query;
  };

  static const std::vector<TautomerPattern> patterns = {
      // Nitro-like resonance/tautomeric query forms
      {"nitro_neutral", "[#7](=[#8])=[#8]"},
      {"nitro_charged", "[#7+](=[#8])-[#8-]"},

      // Common prototropic motifs (keto/enol, amide/imidic, amidine-like)
      {"keto_enol", "[#6X3](=[#8X1])-[#6]"},
      {"amide_imidic", "[#7]-[#6X3](=[#8X1])"},
      {"amidine", "[#7]-[#6X3]=[#7]"},

      // Aromatic [nH] tautomerizable systems
      {"aromatic_nh", "[nH]"},
  };

  for (const auto &pat : patterns) {
    std::unique_ptr<RDKit::ROMol> q(RDKit::SmartsToMol(pat.query));
    if (!q) {
      continue;
    }
    std::vector<RDKit::MatchVectType> matches;
    if (RDKit::SubstructMatch(*mol, *q, matches) && !matches.empty()) {
      if (verbose) {
        std::cout << "inspect_tautomer: matched pattern '" << pat.name
                  << "' | smarts='" << smarts << "'" << std::endl;
      }
      return true;
    }
  }

  return false;
}

bool AtomTyper::is_valid_valence_smarts(const std::string &smarts,
                                        bool verbose) const {
  ZoneScopedN("AtomTyper::is_valid_valence_smarts");
  const auto fail = [&](const std::string &reason) {
    if (inspect_tautomer(smarts, verbose)) {
      if (verbose) {
        std::cout << "is_valid_valence_smarts: tautomer override | reason='"
                  << reason << "' | smarts='" << smarts << "'" << std::endl;
      }
      return true;
    }
    if (verbose) {
      std::cout << "is_valid_valence_smarts: " << reason << " | smarts='"
                << smarts << "'" << std::endl;
    }
    return false;
  };

  if (smarts.empty()) {
    return fail("empty SMARTS");
  }

  std::unique_ptr<RDKit::ROMol> mol(RDKit::SmartsToMol(smarts));
  if (!mol) {
    return fail("SMARTS parse failed");
  }

  mol->updatePropertyCache(false);
  RDKit::MolOps::findSSSR(*mol);

  const auto atom_types = pimpl->process_molecule(mol.get(), true);
  if (atom_types.size() != mol->getNumAtoms()) {
    return fail("typed atom count does not match molecule atom count");
  }

  using BondQuery = Queries::Query<int, const RDKit::Bond *, true>;

  const auto atom_can_be_aromatic = [&](unsigned int atom_idx) -> bool {
    if (atom_idx >= atom_types.size()) {
      return false;
    }
    const auto &at = atom_types[atom_idx];
    return at.is_aromatic;
  };

  const auto atom_can_be_aliphatic = [&](unsigned int atom_idx) -> bool {
    if (atom_idx >= atom_types.size()) {
      return false;
    }
    const auto &at = atom_types[atom_idx];
    return at.is_aliphatic;
  };

  for (const auto *bond : mol->bonds()) {
    if (!bond) {
      continue;
    }

    const unsigned int begin_idx = bond->getBeginAtomIdx();
    const unsigned int end_idx = bond->getEndAtomIdx();
    bool begin_can_aromatic = atom_can_be_aromatic(begin_idx);
    bool end_can_aromatic = atom_can_be_aromatic(end_idx);
    bool begin_can_aliphatic = atom_can_be_aliphatic(begin_idx);
    bool end_can_aliphatic = atom_can_be_aliphatic(end_idx);

    const auto is_primitive_delim = [](char ch) {
      return ch == '[' || ch == ']' || ch == ';' || ch == '&' || ch == ',' ||
             ch == ':' || ch == ')' || ch == '(';
    };

    const auto contains_primitive = [&](const std::string &token, char prim,
                                        bool negated) {
      for (size_t i = 0; i < token.size(); ++i) {
        if (!negated) {
          if (token[i] != prim) {
            continue;
          }
          const bool prev_ok = (i == 0) || is_primitive_delim(token[i - 1]) ||
                               token[i - 1] == '!';
          const bool next_ok =
              (i + 1 >= token.size()) || is_primitive_delim(token[i + 1]);
          if (prev_ok && next_ok) {
            return true;
          }
        } else {
          if (token[i] != '!') {
            continue;
          }
          if (i + 1 < token.size() && token[i + 1] == prim) {
            const bool prev_ok = (i == 0) || is_primitive_delim(token[i - 1]);
            const bool next_ok =
                (i + 2 >= token.size()) || is_primitive_delim(token[i + 2]);
            if (prev_ok && next_ok) {
              return true;
            }
          }
        }
      }
      return false;
    };

    const auto refine_from_atom_smarts =
        [&](unsigned int atom_idx, bool &can_aromatic, bool &can_aliphatic) {
          if (atom_idx >= mol->getNumAtoms()) {
            return;
          }
          const std::string atom_smarts =
              RDKit::SmartsWrite::GetAtomSmarts(mol->getAtomWithIdx(atom_idx));
          const bool has_a = contains_primitive(atom_smarts, 'a', false);
          const bool has_A = contains_primitive(atom_smarts, 'A', false);
          const bool has_not_a = contains_primitive(atom_smarts, 'a', true);
          const bool has_not_A = contains_primitive(atom_smarts, 'A', true);

          if (has_a) {
            can_aliphatic = false;
            can_aromatic = true;
          }
          if (has_A) {
            can_aromatic = false;
            can_aliphatic = true;
          }
          if (has_not_a) {
            // !a means NOT aromatic → must be aliphatic
            can_aromatic = false;
            can_aliphatic = true;
          }
          if (has_not_A) {
            // !A means NOT aliphatic → must be aromatic
            can_aliphatic = false;
            can_aromatic = true;
          }
          // If no a/A primitives at all, atom could be either.
          if (!has_a && !has_A && !has_not_a && !has_not_A) {
            can_aromatic = true;
            can_aliphatic = true;
          }
        };

    refine_from_atom_smarts(begin_idx, begin_can_aromatic, begin_can_aliphatic);
    refine_from_atom_smarts(end_idx, end_can_aromatic, end_can_aliphatic);

    bool has_aromatic_arm = false;
    bool has_aliphatic_arm = false;

    const std::string bond_smarts =
        RDKit::SmartsWrite::GetBondSmarts(bond, begin_idx);
    has_aromatic_arm =
        has_aromatic_arm || (bond_smarts.find(':') != std::string::npos);
    has_aliphatic_arm =
        has_aliphatic_arm || (bond_smarts.find('-') != std::string::npos ||
                              bond_smarts.find('=') != std::string::npos ||
                              bond_smarts.find('#') != std::string::npos ||
                              bond_smarts.find('/') != std::string::npos ||
                              bond_smarts.find('\\') != std::string::npos);

    if (bond->getBondType() == RDKit::Bond::BondType::AROMATIC ||
        bond->getIsAromatic()) {
      has_aromatic_arm = true;
    }
    if (bond->getBondType() == RDKit::Bond::BondType::SINGLE ||
        bond->getBondType() == RDKit::Bond::BondType::DOUBLE ||
        bond->getBondType() == RDKit::Bond::BondType::TRIPLE) {
      has_aliphatic_arm = true;
    }

    std::function<void(const BondQuery *)> inspect_bond_query;
    inspect_bond_query = [&](const BondQuery *query) {
      if (!query) {
        return;
      }

      const auto &desc = query->getDescription();
      if (desc == "BondIsAromatic") {
        has_aromatic_arm = true;
      }

      const auto *eq = dynamic_cast<const RDKit::BOND_EQUALS_QUERY *>(query);
      if (eq) {
        const auto val = static_cast<RDKit::Bond::BondType>(eq->getVal());
        if (val == RDKit::Bond::BondType::AROMATIC) {
          has_aromatic_arm = true;
        }
        if (val == RDKit::Bond::BondType::SINGLE ||
            val == RDKit::Bond::BondType::DOUBLE ||
            val == RDKit::Bond::BondType::TRIPLE) {
          has_aliphatic_arm = true;
        }
      }

      for (auto it = query->beginChildren(); it != query->endChildren(); ++it) {
        inspect_bond_query(it->get());
      }
    };

    inspect_bond_query(bond->getQuery());

    bool bond_has_compatible_assignment = false;
    if (has_aromatic_arm && begin_can_aromatic && end_can_aromatic) {
      bond_has_compatible_assignment = true;
    }
    // Non-aromatic bonds (single/double/triple/etc.) can connect aromatic and
    // aliphatic atoms (e.g. exocyclic substituents on aromatic rings). So for
    // aliphatic bond arms, only require that each endpoint is chemically
    // feasible in at least one state.
    if (has_aliphatic_arm && (begin_can_aromatic || begin_can_aliphatic) &&
        (end_can_aromatic || end_can_aliphatic)) {
      bond_has_compatible_assignment = true;
    }

    // If a bond carries explicit aromatic/aliphatic semantics, require at
    // least one atom-compatibility assignment to exist.
    if ((has_aromatic_arm || has_aliphatic_arm) &&
        !bond_has_compatible_assignment) {
      return fail("bond compatibility failure at bond idx=" + std::to_string(bond->getIdx()) );
          // " (begin_idx=" + std::to_string(begin_idx) + "\n" +
          // (has_aliphatic_arm ? "has aliphatic arm" : "no aliphatic arm") +
          // ", end_idx=" + std::to_string(end_idx) + "\n" +
          // (has_aromatic_arm ? "has aromatic arm" : "no aromatic arm") + "\n" +
          // (end_can_aliphatic ? "end can aliphatic" : "no end can aliphatic") +
          // "\n" +
          // (end_can_aromatic ? "end can aromatic" : "no end can aromatic") +
          // "\n" +
          // (begin_can_aliphatic ? "begin can aliphatic"
          //                      : "no begin can aliphatic") +
          // "\n" +
          // (begin_can_aromatic ? "begin can aromatic"
          //                     : "no begin can aromatic") +
          // "\n" + ") ");
    }
  }

  for (const auto &at : atom_types) {
    const auto is_skippable_atom_pattern = [](std::string token) {
      if (token.size() < 2 || token.front() != '[' || token.back() != ']') {
        return false;
      }
      token = token.substr(1, token.size() - 2);

      const auto colon = token.rfind(':');
      if (colon != std::string::npos && colon + 1 < token.size() &&
          std::all_of(token.begin() + static_cast<std::ptrdiff_t>(colon + 1),
                      token.end(), [](char c) {
                        return std::isdigit(static_cast<unsigned char>(c));
                      })) {
        token = token.substr(0, colon);
      }

      return token == "*" || token == "a";
    };

    if (static_cast<unsigned int>(at.atom_idx) < mol->getNumAtoms()) {
      const std::string atom_smarts =
          RDKit::SmartsWrite::GetAtomSmarts(mol->getAtomWithIdx(at.atom_idx));
      if (is_skippable_atom_pattern(atom_smarts)) {
        continue;
      }
    }

    std::vector<int> h_candidates;
    if (at.explicit_H.has_value()) {
      h_candidates.push_back(*at.explicit_H);
    } else {
      for (int h = 0; h <= 4; ++h) {
        h_candidates.push_back(h);
      }
    }
    h_candidates.erase(
        std::remove_if(h_candidates.begin(), h_candidates.end(),
                       [&](int h) {
                         return std::find(at.excluded_h_counts.begin(),
                                          at.excluded_h_counts.end(),
                                          h) != at.excluded_h_counts.end();
                       }),
        h_candidates.end());

    std::vector<int> charge_candidates;
    if (at.explicit_charge.has_value()) {
      charge_candidates.push_back(*at.explicit_charge);
    } else {
      for (int c = -2; c <= 2; ++c) {
        charge_candidates.push_back(c);
      }
    }
    charge_candidates.erase(
        std::remove_if(charge_candidates.begin(), charge_candidates.end(),
                       [&](int c) {
                         return std::find(at.excluded_charges.begin(),
                                          at.excluded_charges.end(),
                                          c) != at.excluded_charges.end();
                       }),
        charge_candidates.end());

    if (at.explicit_D.has_value() &&
        std::find(at.excluded_D_values.begin(), at.excluded_D_values.end(),
                  *at.explicit_D) != at.excluded_D_values.end()) {
      return fail("atom idx=" + std::to_string(at.atom_idx) +
                  " has explicit D excluded by negation");
    }
    if (at.explicit_X.has_value() &&
        std::find(at.excluded_X_values.begin(), at.excluded_X_values.end(),
                  *at.explicit_X) != at.excluded_X_values.end()) {
      return fail("atom idx=" + std::to_string(at.atom_idx) +
                  " has explicit X excluded by negation");
    }

    std::vector<bool> aromatic_options;
    if (at.explicit_aromatic.has_value()) {
      aromatic_options.push_back(*at.explicit_aromatic);
    } else if (at.explicit_aliphatic.has_value()) {
      aromatic_options.push_back(!(*at.explicit_aliphatic));
    } else if (at.is_aromatic && at.is_aliphatic) {
      aromatic_options.push_back(true);
      aromatic_options.push_back(false);
    } else {
      aromatic_options.push_back(at.is_aromatic);
    }

    aromatic_options.erase(
        std::remove_if(aromatic_options.begin(), aromatic_options.end(),
                       [&](bool a) {
                         return (a && at.excluded_aromatic) ||
                                (!a && at.excluded_aliphatic);
                       }),
        aromatic_options.end());
    std::sort(aromatic_options.begin(), aromatic_options.end());
    aromatic_options.erase(
        std::unique(aromatic_options.begin(), aromatic_options.end()),
        aromatic_options.end());

    if (h_candidates.empty() || charge_candidates.empty() ||
        aromatic_options.empty()) {
      return fail("atom idx=" + std::to_string(at.atom_idx) +
                  " has no feasible H/charge/aromatic assignment candidates");
    }

    bool any_valid = false;

    // Bonds with unspecified/query order in SMARTS (including ':' aromatic
    // query bonds in some parsed forms) may not be reflected in
    // num_single/num_double/num_triple/num_aromatic counts, but they still
    // contribute to topological degree constraints such as D/X.

    int lower_h = at.explicit_lower_h.value_or(0);

    const int known_bond_degree =
        std::max(at.num_single_bonds + at.num_double_bonds +
                     at.num_triple_bonds + at.num_aromatic_bonds,
                 at.min_bonds);
    const int known_degree = known_bond_degree + lower_h;
    const int target_degree =
        at.explicit_D.value_or(at.explicit_X.value_or(known_degree));
    const int missing_degree = std::max(0, target_degree - known_degree);

    const int test_single_bonds = at.num_single_bonds + missing_degree;
    const int test_double_bonds = at.num_double_bonds;
    const int test_triple_bonds = at.num_triple_bonds;
    const int test_aromatic_bonds = at.num_aromatic_bonds;

    for (int h_count : h_candidates) {
      for (int charge : charge_candidates) {
        for (bool aromflag : aromatic_options) {
          // if (verbose) {
          //   std::cout << "Trying aromatic flag=" << aromflag << std::endl;
          // }


          if (pimpl->is_valence_consistent(
                  at, h_count, charge, test_single_bonds, test_double_bonds,
                  test_triple_bonds, test_aromatic_bonds, aromflag,
                  verbose ? DebugLevel::Verbose : DebugLevel::Off)) {
            any_valid = true;
            break;
          }
        }
        if (any_valid) {
          break;
        }
      }
      if (any_valid) {
        break;
      }
    }

    // Fallback: allow protonated ring-hetero assignment when the baseline
    // candidate grid has no valence-consistent solution.
    // This helps tautomeric/prototropic ring systems where a heteroatom can
    // gain one proton with formal +1 while remaining a valid query target.
    if (!any_valid) {


      const bool ring_member = at.explicit_in_ring.value_or(at.is_in_ring);
      const bool hetero_atom = at.atomic_number > 1 && at.atomic_number != 6;
      const bool allow_charge_plus_one =
          !at.explicit_charge.has_value() &&
          std::find(at.excluded_charges.begin(), at.excluded_charges.end(),
                    +1) == at.excluded_charges.end();
      const bool allow_h_increment = !at.explicit_H.has_value();

      if (ring_member && hetero_atom && allow_charge_plus_one &&
          allow_h_increment) {
        std::set<int> boosted_h_candidates;
        for (int h_count : h_candidates) {
          const int boosted_h = h_count + 1;
          if (boosted_h > 4) {
            continue;
          }
          if (std::find(at.excluded_h_counts.begin(), at.excluded_h_counts.end(),
                        boosted_h) != at.excluded_h_counts.end()) {
            continue;
          }
          boosted_h_candidates.insert(boosted_h);
        }

        for (int h_count : boosted_h_candidates) {
          for (bool aromflag : aromatic_options) {
            if (pimpl->is_valence_consistent(
                    at, h_count, +1, test_single_bonds, test_double_bonds,
                    test_triple_bonds, test_aromatic_bonds, aromflag,
                    verbose ? DebugLevel::Verbose : DebugLevel::Off)) {
              any_valid = true;
              if (verbose) {
                std::cout
                    << "is_valid_valence_smarts: accepted protonated ring-hetero "
                    << "fallback at atom idx=" << at.atom_idx
                    << " (Z=" << at.atomic_number << ", H=" << h_count
                    << ", charge=+1)" << std::endl;
              }
              break;
            }
          }
          if (any_valid) {
            break;
          }
        }
      }
    }

    if (!any_valid) {
      return fail("no valence-consistent assignment for atom idx=" +
                  std::to_string(at.atom_idx) +
                  " (Z=" + std::to_string(at.atomic_number) + ") " + std::to_string(test_single_bonds) + " " + std::to_string(test_double_bonds) + " " + std::to_string(test_triple_bonds) + " " + std::to_string(test_aromatic_bonds) + " " + std::to_string(at.is_aromatic) + " " + std::to_string(at.is_aliphatic) + "  " + std::to_string(at.explicit_H.value_or(-1)));
    }
  }

  return true;
}

std::vector<std::string> AtomTyper::filter_invalid_valence_smarts(
    const std::vector<std::string> &patterns) const {
  std::vector<std::string> out;
  out.reserve(patterns.size());
  for (const auto &pat : patterns) {
    if (is_valid_valence_smarts(pat)) {
      out.push_back(pat);
    }
  }
  return out;
}

std::vector<std::string> AtomTyper::consolidate_smarts_by_atom_maps(
    const std::vector<std::string> &patterns) const {
  std::vector<std::string> out;
  if (patterns.empty()) {
    return out;
  }

  const auto strip_brackets = [](const std::string &token) {
    if (token.size() >= 2 && token.front() == '[' && token.back() == ']') {
      return token.substr(1, token.size() - 2);
    }
    return token;
  };

  const auto parse_map_from_inner = [](const std::string &inner) {
    const auto colon = inner.rfind(':');
    if (colon == std::string::npos || colon + 1 >= inner.size()) {
      return 0U;
    }
    for (size_t i = colon + 1; i < inner.size(); ++i) {
      if (!std::isdigit(static_cast<unsigned char>(inner[i]))) {
        return 0U;
      }
    }
    return static_cast<unsigned int>(std::stoul(inner.substr(colon + 1)));
  };

  const auto strip_map_suffix = [&](const std::string &token) {
    std::string inner = strip_brackets(token);
    const auto colon = inner.rfind(':');
    if (colon == std::string::npos) {
      return inner;
    }
    bool all_digits = (colon + 1 < inner.size());
    for (size_t i = colon + 1; i < inner.size() && all_digits; ++i) {
      all_digits = std::isdigit(static_cast<unsigned char>(inner[i]));
    }
    if (all_digits) {
      return inner.substr(0, colon);
    }
    return inner;
  };

  const std::regex atom_re("\\[[^\\]]+\\]");

  const auto split_terms = [](const std::string &inner) {
    std::vector<std::string> terms;
    std::string cur;
    int paren_depth = 0;
    for (char ch : inner) {
      if (ch == '(') {
        ++paren_depth;
      } else if (ch == ')') {
        --paren_depth;
      }
      if (ch == '&' && paren_depth == 0) {
        if (!cur.empty()) {
          terms.push_back(cur);
        }
        cur.clear();
      } else {
        cur.push_back(ch);
      }
    }
    if (!cur.empty()) {
      terms.push_back(cur);
    }
    return terms;
  };

  const auto join_terms = [](const std::vector<std::string> &terms,
                             const std::string &sep) {
    std::stringstream ss;
    for (size_t i = 0; i < terms.size(); ++i) {
      if (i) {
        ss << sep;
      }
      ss << terms[i];
    }
    return ss.str();
  };

  struct ClusterData {
    std::string template_smarts;
    std::map<unsigned int, std::set<std::string>> map_to_tokens;
  };

  // length bucket -> signature bucket
  std::map<size_t, std::map<std::string, ClusterData>> groups;

  for (const auto &pat : patterns) {
    if (pat.empty()) {
      continue;
    }
    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmartsToMol(pat));
    if (!mol) {
      continue;
    }

    size_t mapped_count = 0;
    std::map<unsigned int, std::set<std::string>> local_tokens;
    for (const auto *atom : mol->atoms()) {
      if (!atom) {
        continue;
      }
      const std::string atom_token = RDKit::SmartsWrite::GetAtomSmarts(atom);
      const std::string inner = strip_brackets(atom_token);
      const unsigned int amap = parse_map_from_inner(inner);
      if (amap == 0U) {
        continue;
      }
      ++mapped_count;
      local_tokens[amap].insert(strip_map_suffix(atom_token));
    }

    // Build topology signature by replacing each mapped atom token with
    // a neutral [*:<map>] placeholder while preserving all non-atom text.
    std::string signature;
    size_t cursor = 0;
    for (std::sregex_iterator it(pat.begin(), pat.end(), atom_re), end;
         it != end; ++it) {
      const auto &m = *it;
      const size_t start = static_cast<size_t>(m.position());
      const size_t len = static_cast<size_t>(m.length());
      signature += pat.substr(cursor, start - cursor);

      const std::string tok = m.str();
      const unsigned int amap = parse_map_from_inner(strip_brackets(tok));
      if (amap == 0U) {
        signature += tok;
      } else {
        signature += "[*:" + std::to_string(amap) + "]";
      }
      cursor = start + len;
    }
    signature += pat.substr(cursor);

    auto &g = groups[mapped_count][signature];
    if (g.template_smarts.empty()) {
      g.template_smarts = pat;
    }
    for (const auto &kv : local_tokens) {
      auto &dst = g.map_to_tokens[kv.first];
      dst.insert(kv.second.begin(), kv.second.end());
    }
  }

  for (const auto &len_bucket : groups) {
    for (const auto &sig_bucket : len_bucket.second) {
      const auto &g = sig_bucket.second;
      if (g.template_smarts.empty()) {
        continue;
      }

      std::string consolidated;
      size_t cursor = 0;

      for (std::sregex_iterator
               it(g.template_smarts.begin(), g.template_smarts.end(), atom_re),
           end;
           it != end; ++it) {
        const auto &m = *it;
        const size_t start = static_cast<size_t>(m.position());
        const size_t len = static_cast<size_t>(m.length());
        consolidated += g.template_smarts.substr(cursor, start - cursor);

        const std::string token = m.str();
        const std::string inner = strip_brackets(token);
        const unsigned int amap = parse_map_from_inner(inner);
        if (amap == 0U) {
          consolidated += token;
          cursor = start + len;
          continue;
        }

        auto found = g.map_to_tokens.find(amap);
        if (found == g.map_to_tokens.end() || found->second.empty()) {
          consolidated += token;
          cursor = start + len;
          continue;
        }

        std::vector<std::string> alts(found->second.begin(),
                                      found->second.end());
        if (alts.size() == 1) {
          consolidated += "[" + alts.front() + ":" + std::to_string(amap) + "]";
          cursor = start + len;
          continue;
        }

        // Recursive expressions (especially negated recursive constraints)
        // are sensitive to OR/AND precedence. The common-term factoring below
        // can rewrite
        //   [common;alt1,alt2]
        // in a way that changes semantics for !$()/$(...) terms.
        // Preserve exact alternatives when recursion is present.
        bool has_recursive_alt = false;
        for (const auto &a : alts) {
          if (a.find("$(") != std::string::npos ||
              a.find("!$(") != std::string::npos) {
            has_recursive_alt = true;
            break;
          }
        }
        if (has_recursive_alt) {
          consolidated += "[" + join_terms(alts, ",") + ":" +
                          std::to_string(amap) + "]";
          cursor = start + len;
          continue;
        }

        // Factor common '&'-terms to produce compact forms like:
        // [#6;A;D2;H0&+,H0&-,H1&+0:3]
        std::vector<std::vector<std::string>> alt_terms;
        alt_terms.reserve(alts.size());
        for (const auto &a : alts) {
          alt_terms.push_back(split_terms(a));
        }

        std::vector<std::string> common = alt_terms.front();
        for (size_t ai = 1; ai < alt_terms.size(); ++ai) {
          std::vector<std::string> next_common;
          for (const auto &term : common) {
            if (std::find(alt_terms[ai].begin(), alt_terms[ai].end(), term) !=
                alt_terms[ai].end()) {
              next_common.push_back(term);
            }
          }
          common = std::move(next_common);
        }

        // Preserve original order from first alternative.
        std::vector<std::string> ordered_common;
        for (const auto &term : alt_terms.front()) {
          if (std::find(common.begin(), common.end(), term) != common.end()) {
            ordered_common.push_back(term);
          }
        }

        std::set<std::string> residual_set;
        bool has_empty_residual = false;
        for (const auto &terms : alt_terms) {
          std::vector<std::string> residual;
          for (const auto &t : terms) {
            if (std::find(ordered_common.begin(), ordered_common.end(), t) ==
                ordered_common.end()) {
              residual.push_back(t);
            }
          }
          if (residual.empty()) {
            has_empty_residual = true;
          } else {
            residual_set.insert(join_terms(residual, "&"));
          }
        }

        std::string rebuilt_inner;
        if (ordered_common.empty() || has_empty_residual) {
          // Fallback to simple OR of full alternatives when factoring would be
          // ambiguous.
          rebuilt_inner = join_terms(alts, ",");
        } else {
          std::vector<std::string> residuals(residual_set.begin(),
                                             residual_set.end());
          rebuilt_inner = join_terms(ordered_common, ";") + ";" +
                          join_terms(residuals, ",");
        }

        consolidated += "[" + rebuilt_inner + ":" + std::to_string(amap) + "]";
        cursor = start + len;
      }
      consolidated += g.template_smarts.substr(cursor);

      // Keep only syntactically valid consolidated SMARTS.
      std::unique_ptr<RDKit::ROMol> check(RDKit::SmartsToMol(consolidated));
      if (check) {
        out.push_back(consolidated);
      }
    }
  }

  return out;
}

std::vector<std::string> AtomTyper::consolidate_smarts_by_atom_maps_recanon(
    const std::vector<std::string> &patterns,
    const std::map<std::string, double> &embedding) const {
  if (patterns.empty()) {
    return {};
  }

  const bool trace_enabled =
      static_cast<int>(pimpl->default_debug_level) >=
      static_cast<int>(DebugLevel::Trace);

  const std::map<std::string, double> &active_embedding =
      embedding.empty() ? get_default_query_embedding() : embedding;

  query_reorder::set_comparison_trace_enabled(pimpl->recanon_comparison_trace);

  if (trace_enabled) {
    std::cout << "[consolidate_smarts_by_atom_maps_recanon] input_count="
              << patterns.size() << " embedding="
              << (embedding.empty() ? "default" : "custom")
              << " embedding_size=" << active_embedding.size() << std::endl;
  }

  std::vector<std::string> reordered_patterns;
  reordered_patterns.reserve(patterns.size());

  for (size_t i = 0; i < patterns.size(); ++i) {
    const auto &pat = patterns[i];
    if (pat.empty()) {
      if (trace_enabled) {
        std::cout
            << "[consolidate_smarts_by_atom_maps_recanon] pattern[" << i
            << "] skipped (empty input)" << std::endl;
      }
      continue;
    }

    std::string reordered;
    bool reorder_failed = false;
    try {
      reordered =
          query_reorder::reorder_query_tree_by_embedding_smarts(pat,
                                                                active_embedding);
    } catch (const std::exception &e) {
      reordered = pat;
      reorder_failed = true;
      if (trace_enabled) {
        std::cout << "[consolidate_smarts_by_atom_maps_recanon] pattern["
                  << i << "] reorder_failed=yes (" << e.what()
                  << "), using original" << std::endl;
      }
    }

    reordered_patterns.push_back(reordered);

    if (trace_enabled) {
      std::cout << "[consolidate_smarts_by_atom_maps_recanon] pattern[" << i
                << "] reordered_changed=" << (reordered != pat ? "yes" : "no")
                << " reorder_failed=" << (reorder_failed ? "yes" : "no")
                << " map_normalized=no"
                << "\n  in : " << pat << "\n  out: " << reordered
                << "\n  norm: " << reordered
                << std::endl;
    }
  }

  const auto consolidated = consolidate_smarts_by_atom_maps(reordered_patterns);

  if (trace_enabled) {
    std::cout << "[consolidate_smarts_by_atom_maps_recanon] reordered_count="
              << reordered_patterns.size()
              << " consolidated_count=" << consolidated.size() << std::endl;
    for (size_t i = 0; i < consolidated.size(); ++i) {
      std::cout << "[consolidate_smarts_by_atom_maps_recanon] consolidated["
                << i << "] " << consolidated[i] << std::endl;
    }
  }

  return consolidated;
}

std::string AtomTyper::consolidate_smarts_with_recursive_paths(
    const std::vector<std::string> &patterns) const {
  if (patterns.empty()) {
    return std::string();
  }

  struct ParsedPattern {
    std::string smarts;
    std::unique_ptr<RDKit::ROMol> mol;
    std::map<unsigned int, unsigned int> map_to_atom_idx;
  };

  const auto strip_brackets = [](const std::string &token) {
    if (token.size() >= 2 && token.front() == '[' && token.back() == ']') {
      return token.substr(1, token.size() - 2);
    }
    return token;
  };

  const auto parse_map_from_inner = [](const std::string &inner) {
    const auto colon = inner.rfind(':');
    if (colon == std::string::npos || colon + 1 >= inner.size()) {
      return 0U;
    }
    for (size_t i = colon + 1; i < inner.size(); ++i) {
      if (!std::isdigit(static_cast<unsigned char>(inner[i]))) {
        return 0U;
      }
    }
    return static_cast<unsigned int>(std::stoul(inner.substr(colon + 1)));
  };

  const auto strip_map_suffix = [&](const std::string &token) {
    std::string inner = strip_brackets(token);
    const auto colon = inner.rfind(':');
    if (colon == std::string::npos) {
      return inner;
    }
    bool all_digits = (colon + 1 < inner.size());
    for (size_t i = colon + 1; i < inner.size() && all_digits; ++i) {
      all_digits = std::isdigit(static_cast<unsigned char>(inner[i]));
    }
    if (all_digits) {
      return inner.substr(0, colon);
    }
    return inner;
  };

  std::vector<ParsedPattern> parsed;
  parsed.reserve(patterns.size());
  for (const auto &pat : patterns) {
    if (pat.empty()) {
      continue;
    }
    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmartsToMol(pat));
    if (!mol) {
      continue;
    }

    ParsedPattern entry;
    entry.smarts = pat;
    entry.mol = std::move(mol);
    for (const auto *atom : entry.mol->atoms()) {
      if (!atom) {
        continue;
      }
      const std::string atom_token = RDKit::SmartsWrite::GetAtomSmarts(atom);
      // std::cout << "atom_token: " << atom_token << std::endl;
      const unsigned int amap =
          parse_map_from_inner(strip_brackets(atom_token));
      if (amap != 0U) {
        entry.map_to_atom_idx[amap] = atom->getIdx();
      }
    }
    parsed.push_back(std::move(entry));
  }

  if (parsed.empty()) {
    return std::string();
  }
  if (parsed.size() == 1) {
    return parsed.front().smarts;
  }

  // Use smallest mapped pattern as base scaffold.
  size_t base_idx = 0;
  for (size_t i = 1; i < parsed.size(); ++i) {
    if (parsed[i].map_to_atom_idx.size() <
        parsed[base_idx].map_to_atom_idx.size()) {
      base_idx = i;
    }
  }

  // Common map set across all patterns.
  std::set<unsigned int> common_maps;
  for (const auto &kv : parsed[base_idx].map_to_atom_idx) {
    common_maps.insert(kv.first);
  }
  for (const auto &p : parsed) {
    std::set<unsigned int> keep;
    for (const auto &m : common_maps) {
      if (p.map_to_atom_idx.count(m)) {
        keep.insert(m);
      }
    }
    common_maps = std::move(keep);
  }

  if (common_maps.empty()) {
    return parsed[base_idx].smarts;
  }

  std::function<std::string(const RDKit::ROMol &, unsigned int, int,
                            std::set<unsigned int> &)>
      serialize_subtree;
  serialize_subtree = [&](const RDKit::ROMol &m, unsigned int root, int parent,
                          std::set<unsigned int> &path) -> std::string {
    if (root >= m.getNumAtoms()) {
      return std::string();
    }
    if (path.count(root)) {
      return RDKit::SmartsWrite::GetAtomSmarts(m.getAtomWithIdx(root));
    }

    path.insert(root);
    std::string out = RDKit::SmartsWrite::GetAtomSmarts(m.getAtomWithIdx(root));
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
          serialize_subtree(m, nidx, static_cast<int>(root), path);
      if (!child.empty()) {
        out +=
            "(" + RDKit::SmartsWrite::GetBondSmarts(bond, root) + child + ")";
      }
    }
    path.erase(root);
    return out;
  };

  const auto map_local_signature = [&](const ParsedPattern &p,
                                       unsigned int amap) {
    auto found = p.map_to_atom_idx.find(amap);
    if (found == p.map_to_atom_idx.end()) {
      return std::string("<missing>");
    }
    const auto *atom = p.mol->getAtomWithIdx(found->second);
    std::string sig = strip_map_suffix(RDKit::SmartsWrite::GetAtomSmarts(atom));

    std::vector<std::string> extras;
    for (const auto *nbr : p.mol->atomNeighbors(atom)) {
      const unsigned int nidx = nbr->getIdx();
      const unsigned int nbr_map = parse_map_from_inner(
          strip_brackets(RDKit::SmartsWrite::GetAtomSmarts(nbr)));
      if (nbr_map != 0U && common_maps.count(nbr_map)) {
        continue;
      }
      auto *bond = p.mol->getBondBetweenAtoms(atom->getIdx(), nidx);
      if (!bond) {
        continue;
      }
      std::set<unsigned int> path;
      path.insert(atom->getIdx());
      extras.push_back(RDKit::SmartsWrite::GetBondSmarts(bond, atom->getIdx()) +
                       serialize_subtree(*p.mol, nidx,
                                         static_cast<int>(atom->getIdx()),
                                         path));
    }
    std::sort(extras.begin(), extras.end());
    for (const auto &e : extras) {
      sig += "|" + e;
    }
    return sig;
  };

  std::vector<unsigned int> divergent_maps;
  for (const auto &amap : common_maps) {
    std::set<std::string> local;
    for (const auto &p : parsed) {
      local.insert(map_local_signature(p, amap));
    }
    if (local.size() > 1) {
      divergent_maps.push_back(amap);
    }
  }

  if (divergent_maps.empty()) {
    return parsed[base_idx].smarts;
  }

  const auto canonicalize_without_atom_maps = [&](const std::string &smarts) {
    const std::regex atom_re_local("\\[[^\\]]+\\]");
    std::string out;
    size_t cursor_local = 0;
    for (std::sregex_iterator it(smarts.begin(), smarts.end(), atom_re_local),
         end;
         it != end; ++it) {
      const auto &m = *it;
      const size_t start = static_cast<size_t>(m.position());
      const size_t len = static_cast<size_t>(m.length());
      out += smarts.substr(cursor_local, start - cursor_local);

      const std::string tok = m.str();
      out += "[" + strip_map_suffix(tok) + "]";

      cursor_local = start + len;
    }
    out += smarts.substr(cursor_local);
    return out;
  };

  const auto atom_map_signature = [&](const std::string &smarts) {
    const std::regex atom_re_local("\\[[^\\]]+\\]");
    std::vector<unsigned int> sig;
    for (std::sregex_iterator it(smarts.begin(), smarts.end(), atom_re_local),
         end;
         it != end; ++it) {
      const std::string tok = it->str();
      const unsigned int amap = parse_map_from_inner(strip_brackets(tok));
      if (amap != 0U) {
        sig.push_back(amap);
      }
    }
    return sig;
  };

  const auto recursive_option_to_mapped_single_atom =
      [&](const std::string &recursive_option,
          unsigned int amap) -> std::string {
    if (recursive_option.size() < 4 || recursive_option.rfind("$(", 0) != 0 ||
        recursive_option.back() != ')') {
      return std::string();
    }
    const std::string inner =
        recursive_option.substr(2, recursive_option.size() - 3);
    std::unique_ptr<RDKit::ROMol> inner_mol(RDKit::SmartsToMol(inner));
    if (!inner_mol || inner_mol->getNumAtoms() != 1) {
      return std::string();
    }
    auto *atom = inner_mol->getAtomWithIdx(0);
    if (!atom) {
      return std::string();
    }
    atom->setAtomMapNum(amap);
    return RDKit::MolToSmarts(*inner_mol);
  };

  const auto remap_single_atom_smarts = [&](const std::string &smarts,
                                            unsigned int amap) {
    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmartsToMol(smarts));
    if (!mol || mol->getNumAtoms() != 1) {
      return std::string();
    }
    auto *atom = mol->getAtomWithIdx(0);
    if (!atom) {
      return std::string();
    }
    atom->setAtomMapNum(amap);
    return RDKit::MolToSmarts(*mol);
  };

  const auto has_top_level_comma_in_atom_inner =
      [&](const std::string &single_atom_smarts) {
    if (single_atom_smarts.size() < 2 || single_atom_smarts.front() != '[' ||
        single_atom_smarts.back() != ']') {
      return false;
    }

    const std::string inner =
        single_atom_smarts.substr(1, single_atom_smarts.size() - 2);
    int paren_depth = 0;
    for (const char ch : inner) {
      if (ch == '(') {
        ++paren_depth;
        continue;
      }
      if (ch == ')') {
        --paren_depth;
        continue;
      }
      if (ch == ',' && paren_depth == 0) {
        return true;
      }
    }
    return false;
  };

  std::map<unsigned int, std::string> replacements;
  for (const auto divergent_map : divergent_maps) {
    std::map<std::string, std::string> canonical_to_best_option;
    std::map<std::string, std::vector<unsigned int>> canonical_to_best_sig;
    for (const auto &p : parsed) {
      auto found = p.map_to_atom_idx.find(divergent_map);
      if (found == p.map_to_atom_idx.end()) {
        continue;
      }
      const auto *atom = p.mol->getAtomWithIdx(found->second);
      std::string local =
          "[" + strip_map_suffix(RDKit::SmartsWrite::GetAtomSmarts(atom)) + "]";
      // std::cout << "local: " << local << std::endl;
      std::vector<std::string> extras;
      for (const auto *nbr : p.mol->atomNeighbors(atom)) {
        const unsigned int nidx = nbr->getIdx();
        const unsigned int nbr_map = parse_map_from_inner(
            strip_brackets(RDKit::SmartsWrite::GetAtomSmarts(nbr)));
        if (nbr_map != 0U && common_maps.count(nbr_map)) {
          continue;
        }
        auto *bond = p.mol->getBondBetweenAtoms(atom->getIdx(), nidx);
        if (!bond) {
          continue;
        }
        std::set<unsigned int> path;
        path.insert(atom->getIdx());
        extras.push_back(
            "(" + RDKit::SmartsWrite::GetBondSmarts(bond, atom->getIdx()) +
            serialize_subtree(*p.mol, nidx, static_cast<int>(atom->getIdx()),
                              path) +
            ")");
      }
      std::sort(extras.begin(), extras.end());
      for (const auto &e : extras) {
        local += e;
      }

      const std::string option = "$(" + local + ")";
      // std::cout << "option: " << option << std::endl;
      const std::string canonical = canonicalize_without_atom_maps(option);
      const std::vector<unsigned int> sig = atom_map_signature(option);
      // std::cout << "canonical: " << canonical << std::endl;
      // std::cout << "sig: " << sig.size() << " elements" << std::endl;
      // for (const auto &s : sig) {
        // std::cout << "  " << s << std::endl;
      // }

      auto best_it = canonical_to_best_option.find(canonical);
      if (best_it == canonical_to_best_option.end()) {
        canonical_to_best_option[canonical] = option;
        canonical_to_best_sig[canonical] = sig;
        continue;
      }

      const auto &best_sig = canonical_to_best_sig[canonical];
      if (sig < best_sig) {
        canonical_to_best_option[canonical] = option;
        canonical_to_best_sig[canonical] = sig;
      }
    }

    if (canonical_to_best_option.empty()) {
      // std::cout << "No options found for divergent map " << divergent_map
                // << ", skipping consolidation for this map." << std::endl;
      continue;
    }

    std::vector<std::string> recursive_options;
    recursive_options.reserve(canonical_to_best_option.size());
    for (const auto &kv : canonical_to_best_option) {
      recursive_options.push_back(kv.second);
      // std::cout << "recursive option for map " << divergent_map << ": "
                // << kv.second << std::endl;
    }
    std::sort(recursive_options.begin(), recursive_options.end());

    // Second-stage consolidation for recursive OR alternatives:
    // convert $(...) options into mapped single-atom SMARTS, consolidate them
    // with map-aware factoring, then use the compact one-atom result when
    // possible.
    {
      std::vector<std::string> mapped_recursive_options;
      mapped_recursive_options.reserve(recursive_options.size());
      bool all_convertible = true;
      bool has_nested_or_in_input = false;
      for (const auto &opt : recursive_options) {
        const std::string mapped_one_atom =
            recursive_option_to_mapped_single_atom(opt, divergent_map);
        if (mapped_one_atom.empty()) {
          all_convertible = false;
          break;
        }
        if (has_top_level_comma_in_atom_inner(mapped_one_atom)) {
          has_nested_or_in_input = true;
          break;
        }
        mapped_recursive_options.push_back(mapped_one_atom);
      }

      if (all_convertible && !has_nested_or_in_input &&
          !mapped_recursive_options.empty()) {
        const auto compact_options =
            consolidate_smarts_by_atom_maps_recanon(mapped_recursive_options);
        if (compact_options.size() == 1) {
          const std::string remapped =
              remap_single_atom_smarts(compact_options.front(), divergent_map);
          if (!remapped.empty()) {
            replacements[divergent_map] = remapped;
            continue;
          }
        }
      }
    }

    std::stringstream replacement_ss;
    replacement_ss << "[";
    bool first = true;
    for (const auto &opt : recursive_options) {
      if (!first) {
        replacement_ss << ",";
      }
      replacement_ss << opt;
      first = false;
    }
    replacement_ss << ":" << divergent_map << "]";
    replacements[divergent_map] = replacement_ss.str();
  }

  if (replacements.empty()) {
    return parsed[base_idx].smarts;
  }

  const std::string &base_smarts = parsed[base_idx].smarts;
  std::unique_ptr<RDKit::ROMol> base_mol(RDKit::SmartsToMol(base_smarts));
  if (!base_mol) {
    return parsed[base_idx].smarts;
  }

  RDKit::RWMol rewritten(*base_mol);
  std::vector<std::unique_ptr<RDKit::ROMol>> replacement_owners;
  replacement_owners.reserve(replacements.size());
  std::vector<unsigned int> recursive_replaced_atom_indices;

  for (unsigned int idx = 0; idx < rewritten.getNumAtoms(); ++idx) {
    const auto *atom = rewritten.getAtomWithIdx(idx);
    if (!atom) {
      continue;
    }

    const std::string atom_token = RDKit::SmartsWrite::GetAtomSmarts(atom);
    const unsigned int amap = parse_map_from_inner(strip_brackets(atom_token));
    if (amap == 0U) {
      continue;
    }

    const auto rep_it = replacements.find(amap);
    if (rep_it == replacements.end()) {
      continue;
    }
    const bool replacement_has_recursive =
        rep_it->second.find("$(") != std::string::npos;

    std::unique_ptr<RDKit::ROMol> one_atom_mol(
        RDKit::SmartsToMol(rep_it->second));
    if (!one_atom_mol || one_atom_mol->getNumAtoms() != 1) {
      continue;
    }

    auto *src_atom = one_atom_mol->getAtomWithIdx(0);
    rewritten.replaceAtom(idx, src_atom, false, true);
    replacement_owners.push_back(std::move(one_atom_mol));
    if (replacement_has_recursive) {
      recursive_replaced_atom_indices.push_back(idx);
    }
  }

  if (!recursive_replaced_atom_indices.empty()) {
    std::set<unsigned int> atoms_to_remove;

    const auto should_keep_component =
        [&](unsigned int start_idx,
            unsigned int blocked_idx,
            std::vector<unsigned int> &component) {
          component.clear();
          std::set<unsigned int> visited;
          std::vector<unsigned int> stack;
          stack.push_back(start_idx);
          bool has_common_mapped_atom = false;

          while (!stack.empty()) {
            const unsigned int current = stack.back();
            stack.pop_back();
            if (!visited.insert(current).second) {
              continue;
            }
            component.push_back(current);

            const auto *cur_atom = rewritten.getAtomWithIdx(current);
            if (cur_atom && common_maps.count(cur_atom->getAtomMapNum())) {
              has_common_mapped_atom = true;
            }

            for (const auto *nbr : rewritten.atomNeighbors(cur_atom)) {
              const unsigned int nidx = nbr->getIdx();
              if (nidx == blocked_idx) {
                continue;
              }
              if (!visited.count(nidx)) {
                stack.push_back(nidx);
              }
            }
          }
          return has_common_mapped_atom;
        };

    for (const auto anchor_idx : recursive_replaced_atom_indices) {
      if (anchor_idx >= rewritten.getNumAtoms()) {
        continue;
      }
      const auto *anchor_atom = rewritten.getAtomWithIdx(anchor_idx);
      if (!anchor_atom) {
        continue;
      }

      std::vector<unsigned int> neighbor_indices;
      neighbor_indices.reserve(anchor_atom->getDegree());
      for (const auto *nbr : rewritten.atomNeighbors(anchor_atom)) {
        neighbor_indices.push_back(nbr->getIdx());
      }

      for (const auto nidx : neighbor_indices) {
        const auto *nbr_atom = rewritten.getAtomWithIdx(nidx);
        if (!nbr_atom) {
          continue;
        }
        const unsigned int nbr_map = nbr_atom->getAtomMapNum();
        if (nbr_map != 0U && common_maps.count(nbr_map)) {
          continue;
        }

        std::vector<unsigned int> component;
        const bool keep_component =
            should_keep_component(nidx, anchor_idx, component);
        if (!keep_component) {
          atoms_to_remove.insert(component.begin(), component.end());
        }
      }
    }

    std::vector<unsigned int> remove_order(atoms_to_remove.begin(),
                                           atoms_to_remove.end());
    std::sort(remove_order.begin(), remove_order.end(), std::greater<>());
    for (const auto atom_idx : remove_order) {
      if (atom_idx < rewritten.getNumAtoms()) {
        rewritten.removeAtom(atom_idx);
      }
    }
  }

  std::string consolidated = RDKit::MolToSmarts(rewritten);
  std::unique_ptr<RDKit::ROMol> check(RDKit::SmartsToMol(consolidated));
  if (!check) {
    return parsed[base_idx].smarts;
  }
  return consolidated;
}

std::vector<PatternItem> AtomTyper::type_pattern_from_smarts(
    const std::string &smarts) {
  if (smarts.empty()) {
    throw std::runtime_error("Empty SMARTS string");
  }

  std::unique_ptr<RDKit::ROMol> mol(RDKit::SmartsToMol(smarts));
  if (!mol) {
    throw std::runtime_error("Failed to parse SMARTS: " + smarts);
  }

  mol->updatePropertyCache();
  RDKit::MolOps::findSSSR(*mol);

  auto atom_types = pimpl->process_molecule(mol.get(), true);
  std::vector<PatternItem> items;
  if (atom_types.empty()) {
    return items;
  }

  items.reserve(atom_types.size() * 2 - 1);
  for (size_t i = 0; i < atom_types.size(); ++i) {
    PatternItem atom_item;
    atom_item.kind = PatternItemKind::Atom;
    atom_item.atom = atom_types[i];
    items.push_back(atom_item);

    if (i + 1 < atom_types.size()) {
      if (auto *bond = mol->getBondBetweenAtoms(
              static_cast<unsigned int>(i), static_cast<unsigned int>(i + 1))) {
        PatternItem bond_item;
        bond_item.kind = PatternItemKind::Bond;
        bond_item.bond = pimpl->make_bond_type(bond);
        items.push_back(bond_item);
      }
    }
  }

  return items;
}

std::string AtomTyper::get_atom_types_string(
    const std::vector<AtomType> &atom_types) {
  std::stringstream ss;
  const auto optional_int_to_string = [](const std::optional<int> &v) {
    return v.has_value() ? std::to_string(*v) : std::string("None");
  };
  const auto optional_bool_to_string = [](const std::optional<bool> &v) {
    if (!v.has_value()) {
      return std::string("None");
    }
    return *v ? std::string("true") : std::string("false");
  };
  ss << "Atom Types:\n";
  ss << "===========\n";

  for (const auto &at : atom_types) {
    ss << "Atom " << at.atom_idx << ": "
       << "Element=" << at.atomic_number << ", MinBonds=" << at.min_bonds
       << ", MaxValence=" << at.max_valence
       << ", Aromatic=" << (at.is_aromatic ? "Yes" : "No")
       << ", Aliphatic=" << (at.is_aliphatic ? "Yes" : "No")
       << ", InRing=" << (at.is_in_ring ? "Yes" : "No")
       << ", RingSize=" << at.ring_size << ", Hybrid=" << at.hybridization
       << ", AtomTypeEnum=" << at.atom_type_enumeration
       << ", SingleBonds=" << at.num_single_bonds
       << ", DoubleBonds=" << at.num_double_bonds
       << ", TripleBonds=" << at.num_triple_bonds
       << ", AromaticBonds=" << at.num_aromatic_bonds
       << ", RemainingValence=" << at.remaining_valence
       << ", SMARTS=" << at.smarts_pattern;
    ss << ", explicit_atomic_num="
       << optional_int_to_string(at.explicit_atomic_num)
       << ", explicit_charge=" << optional_int_to_string(at.explicit_charge)
       << ", explicit_H=" << optional_int_to_string(at.explicit_H)
       << ", explicit_lower_h=" << optional_int_to_string(at.explicit_lower_h)
       << ", explicit_D=" << optional_int_to_string(at.explicit_D)
       << ", explicit_X=" << optional_int_to_string(at.explicit_X)
       << ", explicit_valence=" << optional_int_to_string(at.explicit_valence)
       << ", explicit_hybridization="
       << optional_int_to_string(at.explicit_hybridization)
       << ", explicit_aromatic="
       << optional_bool_to_string(at.explicit_aromatic)
       << ", explicit_aliphatic="
       << optional_bool_to_string(at.explicit_aliphatic)
       << ", explicit_in_ring=" << optional_bool_to_string(at.explicit_in_ring)
       << ", explicit_min_ring_size="
       << optional_int_to_string(at.explicit_min_ring_size);

    if (!at.bond_details.empty()) {
      ss << ", BondTypes=";
      for (size_t i = 0; i < at.bond_details.size(); ++i) {
        const auto &bt = at.bond_details[i];
        ss << "(" << bt.begin_atom_idx << "-" << bt.end_atom_idx << ":"
           << bt.bond_type_name << ",SMARTS=" << bt.smarts_pattern << ")";
        if (i + 1 < at.bond_details.size()) {
          ss << ";";
        }
      }
    }

    ss << "\n";
  }

  return ss.str();
}

std::string AtomTyper::get_smarts_from_pattern_types(
    const std::vector<PatternItem> &items) {
  std::stringstream ss;
  for (const auto &item : items) {
    if (item.kind == PatternItemKind::Atom) {
      ss << item.atom.smarts_pattern;
    } else {
      ss << item.bond.smarts_pattern;
    }
  }
  return ss.str();
}

std::string AtomTyper::get_pattern_types_string(
    const std::vector<PatternItem> &items) {
  std::stringstream ss;
  const auto optional_int_to_string = [](const std::optional<int> &v) {
    return v.has_value() ? std::to_string(*v) : std::string("None");
  };
  const auto optional_bool_to_string = [](const std::optional<bool> &v) {
    if (!v.has_value()) {
      return std::string("None");
    }
    return *v ? std::string("true") : std::string("false");
  };
  ss << "Pattern Types:\n";
  ss << "==============\n";

  for (const auto &item : items) {
    if (item.kind == PatternItemKind::Atom) {
      const auto &at = item.atom;
      ss << "Atom " << at.atom_idx << ": "
         << "Element=" << at.atomic_number << ", MinBonds=" << at.min_bonds
         << ", MaxValence=" << at.max_valence
         << ", Aromatic=" << (at.is_aromatic ? "Yes" : "No")
         << ", Aliphatic=" << (at.is_aliphatic ? "Yes" : "No")
         << ", InRing=" << (at.is_in_ring ? "Yes" : "No")
         << ", RingSize=" << at.ring_size << ", Hybrid=" << at.hybridization
         << ", SingleBonds=" << at.num_single_bonds
         << ", DoubleBonds=" << at.num_double_bonds
         << ", TripleBonds=" << at.num_triple_bonds
         << ", AromaticBonds=" << at.num_aromatic_bonds
         << ", RemainingValence=" << at.remaining_valence
         << ", SMARTS=" << at.smarts_pattern << ", explicit_atomic_num="
         << optional_int_to_string(at.explicit_atomic_num)
         << ", explicit_charge=" << optional_int_to_string(at.explicit_charge)
         << ", explicit_H=" << optional_int_to_string(at.explicit_H)
         << ", explicit_lower_h=" << optional_int_to_string(at.explicit_lower_h)
         << ", explicit_D=" << optional_int_to_string(at.explicit_D)
         << ", explicit_X=" << optional_int_to_string(at.explicit_X)
         << ", explicit_valence=" << optional_int_to_string(at.explicit_valence)
         << ", explicit_hybridization="
         << optional_int_to_string(at.explicit_hybridization)
         << ", explicit_aromatic="
         << optional_bool_to_string(at.explicit_aromatic)
         << ", explicit_aliphatic="
         << optional_bool_to_string(at.explicit_aliphatic)
         << ", explicit_in_ring="
         << optional_bool_to_string(at.explicit_in_ring)
         << ", explicit_min_ring_size="
         << optional_int_to_string(at.explicit_min_ring_size);
      ss << "\n";
    } else {
      const auto &bt = item.bond;
      ss << "Bond " << bt.bond_idx << ": " << bt.begin_atom_idx << "-"
         << bt.end_atom_idx << ", Type=" << bt.bond_type_name
         << ", SMARTS=" << bt.smarts_pattern
         << ", Aromatic=" << (bt.is_aromatic ? "Yes" : "No")
         << ", InRing=" << (bt.is_in_ring ? "Yes" : "No") << "\n";
    }
  }

  return ss.str();
}

void AtomTyper::set_use_canonical(bool use_canonical) {
  pimpl->use_canonical = use_canonical;
}

void AtomTyper::set_debug_level(DebugLevel level) {
  pimpl->default_debug_level = level;
}

void AtomTyper::set_recanon_comparison_trace(bool enabled) {
  pimpl->recanon_comparison_trace = enabled;
  query_reorder::set_comparison_trace_enabled(enabled);
}

}  // namespace atom_typer
