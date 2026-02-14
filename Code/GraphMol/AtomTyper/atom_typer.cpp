#include "atom_typer.hpp"
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/RWMol.h>
#include <Query/EqualityQuery.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace atom_typer {

/**
 * Private implementation class (PIMPL pattern)
 */
class AtomTyper::Impl {
 public:
  bool use_canonical = false;
  DebugLevel default_debug_level = DebugLevel::Off;

  struct EnumerationSettings {
    int h_min = 0;
    int h_max = 4;
    int charge_min = -1;
    int charge_max = 1;
    DebugLevel debug_level = DebugLevel::Off;
  };

  static int debug_level_to_int(DebugLevel level) {
    return static_cast<int>(level);
  }

  struct QueryConstraints {
    std::optional<int> atomic_number;
    std::optional<int> formal_charge;
    std::optional<int> h_count;
    std::optional<int> explicit_degree_query;  // SMARTS D primitive
    std::optional<int> total_degree;
    std::optional<int> valence;
    std::optional<int> hybridization;
    std::optional<bool> aromatic;
    std::optional<bool> aliphatic;
    std::optional<bool> in_ring;
    std::optional<int> min_ring_size;
    int min_bonds = 0;  // explicit topology neighbors in SMARTS graph
  };

  void collect_query_constraints(const RDKit::Atom::QUERYATOM_QUERY *query,
                                 QueryConstraints &constraints) {
    if (!query) {
      return;
    }

    const auto &descr = query->getDescription();
    const bool negated = query->getNegation();
    if (!negated) {
      auto eq = dynamic_cast<const RDKit::ATOM_EQUALS_QUERY *>(query);
      if ((descr == "AtomAtomicNum" || descr == "AtomNum") && eq) {
        constraints.atomic_number = eq->getVal();
      } else if (descr == "AtomFormalCharge" && eq) {
        constraints.formal_charge = eq->getVal();
      } else if ((descr == "AtomHCount" || descr == "AtomImplicitHCount") &&
                 eq) {
        constraints.h_count = eq->getVal();
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
      } else if (descr == "AtomMinRingSize" && eq) {
        constraints.min_ring_size = eq->getVal();
        constraints.in_ring = true;
      }
    }

    for (auto it = query->beginChildren(); it != query->endChildren(); ++it) {
      collect_query_constraints(it->get(), constraints);
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
          default:
            break;
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
    if (!atom->hasQuery()) {
      return at;
    }

    QueryConstraints constraints;
    collect_query_constraints(atom->getQuery(), constraints);
    at.explicit_atomic_num = constraints.atomic_number;
    at.explicit_charge = constraints.formal_charge;
    at.explicit_H = constraints.h_count;
    at.explicit_D = constraints.explicit_degree_query;
    at.explicit_X = constraints.total_degree;
    at.explicit_valence = constraints.valence;
    at.explicit_hybridization = constraints.hybridization;
    at.explicit_aromatic = constraints.aromatic;
    at.explicit_aliphatic = constraints.aliphatic;
    at.explicit_in_ring = constraints.in_ring;
    at.explicit_min_ring_size = constraints.min_ring_size;

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
    const bool has_aromatic_constraint = constraints.aromatic.has_value();
    const bool has_aliphatic_constraint = constraints.aliphatic.has_value();
    if (has_aromatic_constraint) {
      at.is_aromatic = *constraints.aromatic;
      at.is_aliphatic = false;
    }
    if (has_aliphatic_constraint) {
      at.is_aliphatic = *constraints.aliphatic;
      at.is_aromatic = false;
    }

    // SMARTS like [#6] often carries only atomic-number constraints and should
    // be treated as potentially matching either aliphatic (A) or aromatic (a).
    if (constraints.atomic_number.has_value() && !has_aromatic_constraint &&
        !has_aliphatic_constraint) {
      at.is_aromatic = true;
      at.is_aliphatic = true;
    }
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

    // if (constraints.valence.has_value()) {
    //   const double used = used_connectivity_valence(atom, mol);
    //   at.remaining_valence =
    //       std::max(0, at.max_valence - static_cast<int>(std::floor(used)));
    // } else {
    //   at.remaining_valence =
    //       static_cast<int>(remaining_valence_from_connectivity(
    //           at.atomic_number, at.formal_charge, atom, mol));
    // }
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
    // temp_atom->setNumExplicitHs(h_count);

    const bool want_aromatic = arom;
    temp_atom->setIsAromatic(want_aromatic);

    const unsigned int aidx = temp_mol.addAtom(temp_atom, false, true);

    // --- if aromatic, embed into a minimal aromatic ring (benzene-like) ---
    // This guarantees sanitization/aromatic valence bookkeeping is possible.
    if (want_aromatic) {
      // In an aromatic ring, the atom will have 2 aromatic bonds in-ring.
      // If your counting expects something else, you'll need a fused-ring
      // builder. For now, we treat ring bonds as satisfying aromaticity
      // requirements.
      if (aromatic_bonds != 0 && aromatic_bonds != 2) {
        // Conservative: refuse cases we can't model without fused systems.
        // (You can relax this if you implement fused-ring construction.)
        return false;
      }

      // Build 5 additional aromatic atoms to make a 6-member ring
      unsigned int ring[6];
      ring[0] = aidx;

      for (int i = 1; i < 6; ++i) {
        auto *c = new RDKit::Atom(6);  // aromatic carbon
        c->setIsAromatic(true);
        ring[i] = temp_mol.addAtom(c, false, true);
      }

      // Connect ring with aromatic bonds: 0-1-2-3-4-5-0
      for (int i = 0; i < 6; ++i) {
        int j = (i + 1) % 6;
        temp_mol.addBond(ring[i], ring[j], RDKit::Bond::BondType::AROMATIC);
        auto *b = temp_mol.getBondBetweenAtoms(ring[i], ring[j]);
        b->setIsAromatic(true);
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

    // IMPORTANT: do NOT add aromatic bonds to standalone neighbors here.
    // Aromatic bond types generally require both atoms to be aromatic *and*
    // part of an aromatic system; otherwise sanitization/property cache will
    // complain. We treat aromaticity as handled by the ring embedding above.

    // --- force valence computation / cache update ---
    try {
      // Minimal sanitation that updates properties/valences without trying to
      // perceive aromaticity from scratch (we already set aromatic
      // flags/bonds).
      unsigned int san;
      //   RDKit::MolOps::sanitizeMol(temp_mol, san,
      //                              RDKit::MolOps::SANITIZE_PROPERTIES |
      //                                  RDKit::MolOps::SANITIZE_SYMMRINGS |
      //                                  RDKit::MolOps::SANITIZE_KEKULIZE |
      //                                  RDKit::MolOps::SANITIZE_SETAROMATICITY);
      RDKit::MolOps::sanitizeMol(temp_mol);

      //   std::cerr << "Sanitization flags: " << san << std::endl;
      //   temp_mol.updatePropertyCache(true);
    } catch (...) {
      if (debug_level_to_int(debug_level) >=
          debug_level_to_int(DebugLevel::Trace)) {
        std::cerr << "Sanitization failed for candidate: Z=" << at.atomic_number
                  << " H=" << h_count << " q=" << charge
                  << " s=" << single_bonds << " d=" << double_bonds
                  << " t=" << triple_bonds << " a=" << aromatic_bonds
                  << " arom=" << (arom ? 1 : 0) << std::endl;
      }
      return false;
    }

    if (temp_mol.getAtomWithIdx(aidx)->getTotalNumHs() != h_count) {
      return false;
    }

    if (debug_level_to_int(debug_level) >=
        debug_level_to_int(DebugLevel::Verbose) &&
      !temp_mol.getAtomWithIdx(aidx)->hasValenceViolation()) {
      int deg = single_bonds + double_bonds + triple_bonds + aromatic_bonds;
      std::cerr << RDKit::MolToSmiles(temp_mol) << " " << h_count << " "
                << charge << " " << single_bonds << " " << double_bonds << " "
                << triple_bonds << " " << aromatic_bonds << " "
                << temp_mol.getAtomWithIdx(0)->getNumRadicalElectrons() << " "
                << deg << std::endl;
    }

    return !temp_mol.getAtomWithIdx(aidx)->hasValenceViolation();
  }
  std::string enumerate_atom_alternatives(const PatternItem &at,
                                          const EnumerationSettings &settings) {
    struct Candidate {
      int degree;
      int h_count;
      int charge;
      //   int single_bonds;
      //   int double_bonds;
      //   int triple_bonds;
      //   int aromatic_bonds;
      int arom;
    };

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

    std::vector<int> charges;
    if (at.atom.explicit_charge.has_value()) {
      charges.push_back(*at.atom.explicit_charge);
    } else {
      const int charge_min = std::min(settings.charge_min, settings.charge_max);
      const int charge_max = std::max(settings.charge_min, settings.charge_max);
      for (int c = charge_min; c <= charge_max; ++c) {
        charges.push_back(c);
      }
    }

    std::vector<int> degree_values;
    if (at.atom.explicit_D.has_value()) {
      degree_values.push_back(*at.atom.explicit_D);
    } else {
      for (int d = at.atom.min_bonds; d <= 5; ++d) {
        degree_values.push_back(d);
      }
    }

    std::vector<int> single_bond_values;
    for (int d = at.atom.num_single_bonds; d <= 4; ++d) {
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
    if (at.atom.is_aromatic) {
      arom.push_back(true);
      for (int d = at.atom.num_aromatic_bonds; d <= 3; ++d) {
        aromatic_bond_values.push_back(d);
      }
    } else {
      aromatic_bond_values.push_back(0);
    }
    if (at.atom.is_aliphatic) {
      arom.push_back(false);
    }

    // std::cerr << arom.size() << std::endl;
    // for (int a : arom) {
    //   std::cerr << "Arom: " << a << std::endl;
    // }

    std::vector<Candidate> candidates;
    std::set<std::tuple<int, int, int, int>> seen_candidates;
    candidates.reserve(64);
    for (int h_count : hs) {
      for (int charge : charges) {
        for (int single_bonds : single_bond_values) {
          for (int double_bonds : double_bond_values) {
            for (int triple_bonds : triple_bond_values) {
              for (int aromatic_bonds : aromatic_bond_values) {
                for (bool aromflag : arom) {
                  if (!aromflag) {
                    aromatic_bonds = 0;
                  }
                  if (!is_valence_consistent(
                          at.atom, h_count, charge, single_bonds, double_bonds,
                      triple_bonds, aromatic_bonds, aromflag,
                      settings.debug_level)) {
                    continue;
                  }
                  int degree = single_bonds + double_bonds + triple_bonds +
                               aromatic_bonds;
                  //   candidates.push_back({degree, h_count, charge,
                  //   single_bonds,
                  //                         double_bonds, triple_bonds,
                  //                         aromatic_bonds, aromflag});
                  auto key = std::make_tuple(degree, h_count, charge,
                                             aromflag ? 1 : 0);
                  if (seen_candidates.insert(key).second) {
                    candidates.push_back(
                        {degree, h_count, charge, aromflag ? 1 : 0});
                  }
                }
              }
            }
          }
        }
      }
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
                // if (a.single_bonds != b.single_bonds) {
                //   return a.single_bonds < b.single_bonds;
                // }
                // if (a.double_bonds != b.double_bonds) {
                //   return a.double_bonds < b.double_bonds;
                // }
                // if (a.triple_bonds != b.triple_bonds) {
                //   return a.triple_bonds < b.triple_bonds;
                // }
                // if (a.aromatic_bonds != b.aromatic_bonds) {
                //   return a.aromatic_bonds < b.aromatic_bonds;
                // }
                if (a.arom != b.arom) {
                  return a.arom < b.arom;
                }
                return charge_sort_key(a.charge) < charge_sort_key(b.charge);
              });
    // for (const auto &c : candidates) {
    //   std::cerr << "Candidate: D=" << c.degree << " H=" << c.h_count
    //             << " Charge=" << c.charge << " Single=" << c.single_bonds
    //             << " Double=" << c.double_bonds << " Triple=" <<
    //             c.triple_bonds
    //             << " Aromatic=" << c.aromatic_bonds << " Arom=" << c.arom
    //             << std::endl;
    // }
    if (candidates.empty()) {
      std::stringstream empty;
      empty << "[" << atom_symbol_for_smarts(at.atom) << "]";
      return empty.str();
    }

    const auto all_same_degree = std::all_of(
        candidates.begin(), candidates.end(), [&](const Candidate &c) {
          return c.degree == candidates.front().degree;
        });
    const auto all_same_h = std::all_of(
        candidates.begin(), candidates.end(), [&](const Candidate &c) {
          return c.h_count == candidates.front().h_count;
        });
    const auto all_same_charge = std::all_of(
        candidates.begin(), candidates.end(), [&](const Candidate &c) {
          return c.charge == candidates.front().charge;
        });
    const auto all_same_arom = std::all_of(
        candidates.begin(), candidates.end(),
        [&](const Candidate &c) { return c.arom == candidates.front().arom; });

    const bool no_varying_terms =
        all_same_degree && all_same_h && all_same_charge && all_same_arom;

    std::stringstream out;
    // const auto symbols = atom_symbol_for_smarts(at.atom);
    out << "[" << "#" << std::to_string(at.atom.atomic_number);;
    // for (size_t i = 0; i < symbols.size(); ++i) {
    //   out << symbols[i];
    //   if (i + 1 < symbols.size()) {
    //     out << ",";
    //   }
    // }

    if (all_same_h) {
      out << ";H" << candidates.front().h_count;
    }
    if (all_same_degree) {
      out << ";D" << candidates.front().degree;
    }
    if (all_same_charge) {
      out << ";" << charge_token(candidates.front().charge);
    }

    if (all_same_arom) {
      if (candidates.front().arom) {
        out << ";a";
      } else {
        out << ";A";
      }
    }

    if (no_varying_terms) {
      out << "]";
      return out.str();
    }

    out << ";";
    for (size_t i = 0; i < candidates.size(); ++i) {
      bool wrote_any = false;
      if (!all_same_degree) {
        out << "D" << candidates[i].degree;
        wrote_any = true;
      }
      if (!all_same_h) {
        if (wrote_any) {
          out << "&";
        }
        out << "H" << candidates[i].h_count;
        wrote_any = true;
      }
      if (!all_same_charge) {
        if (wrote_any) {
          out << "&";
        }
        out << charge_token(candidates[i].charge);
        wrote_any = true;
      }

      if (!all_same_arom) {
        if (wrote_any) {
          out << "&";
        }
        out << (candidates[i].arom ? "a" : "A");
        wrote_any = true;
      }

      if (!wrote_any) {
        out << "D" << candidates[i].degree << "&H" << candidates[i].h_count
            << "&" << charge_token(candidates[i].charge);
      }

      if (i + 1 < candidates.size()) {
        out << ",";
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

std::vector<AtomType> AtomTyper::type_atoms_from_smiles(
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

  // Process molecule and extract atom types
  auto atom_types = pimpl->process_molecule(mol.get(), false);

  return atom_types;
}

std::vector<AtomType> AtomTyper::type_atoms_from_smarts(
    const std::string &smarts) {
  if (smarts.empty()) {
    throw std::runtime_error("Empty SMARTS string");
  }

  // Parse SMARTS
  std::unique_ptr<RDKit::ROMol> mol(RDKit::SmartsToMol(smarts));
  if (!mol) {
    throw std::runtime_error("Failed to parse SMARTS: " + smarts);
  }

  // Update property cache and find rings for SMARTS molecules
  mol->updatePropertyCache();
  RDKit::MolOps::findSSSR(*mol);

  return pimpl->process_molecule(mol.get(), true);
}

std::string AtomTyper::enumerate_dof_smarts(const std::string &smarts) {
  return enumerate_dof_smarts(smarts, 0, 4, -1, 1, true,
                              pimpl->default_debug_level);
}

std::string AtomTyper::enumerate_dof_smarts(const std::string &smarts,
                                            int h_min, int h_max,
                                            int charge_min,
                                            int charge_max) {
  return enumerate_dof_smarts(smarts, h_min, h_max, charge_min, charge_max,
                              true, pimpl->default_debug_level);
}

std::string AtomTyper::enumerate_dof_smarts(const std::string &smarts,
                                            int h_min, int h_max,
                                            int charge_min, int charge_max,
                                            bool verbose,
                                            DebugLevel debug_level) {
  if (smarts.empty()) {
    throw std::runtime_error("Empty SMARTS string");
  }

  std::unique_ptr<RDKit::ROMol> mol(RDKit::SmartsToMol(smarts));
  if (!mol) {
    throw std::runtime_error("Failed to parse SMARTS: " + smarts);
  }

  //   mol->updatePropertyCache();
  //   RDKit::MolOps::findSSSR(*mol);

  AtomTyper::Impl::EnumerationSettings settings;
  settings.h_min = std::min(h_min, h_max);
  settings.h_max = std::max(h_min, h_max);
  settings.charge_min = std::min(charge_min, charge_max);
  settings.charge_max = std::max(charge_min, charge_max);
  settings.debug_level = verbose ? debug_level : DebugLevel::Off;

  const auto atom_types = type_pattern_from_smarts(smarts);

  std::stringstream out;
  for (const auto &atom_type : atom_types) {
    if (atom_type.kind == PatternItemKind::Atom) {
      out << pimpl->enumerate_atom_alternatives(atom_type, settings);
    } else {
      out << atom_type.bond.smarts_pattern;
    }
  }

  return out.str();
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

}  // namespace atom_typer
