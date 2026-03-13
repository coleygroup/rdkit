#include "reaction_extractor.hpp"

#include "expression_builder.hpp"

#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include <GraphMol/Subgraphs/Subgraphs.h>

#include <boost/shared_ptr.hpp>

#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace atom_typer {
namespace {

using MolPtr = boost::shared_ptr<RDKit::ROMol>;
using BondType = RDKit::Bond::BondType;
using PathType = std::vector<int>;

static unsigned add_unmapped_atom_and_bond(
    const std::string &to_add, int idx_to_add_to, BondType bond_type_to_add,
    std::unique_ptr<RDKit::RWMol> &merged_template,
    std::unordered_map<const RDKit::Atom *, int> &unmapped_idx_map,
    const RDKit::Atom *source_atom_ptr) {
  unsigned ni;
  auto it = unmapped_idx_map.find(source_atom_ptr);
  if (it == unmapped_idx_map.end()) {
    std::unique_ptr<RDKit::ROMol> one(RDKit::SmartsToMol(to_add));
    if (!one || one->getNumAtoms() != 1) {
      throw std::invalid_argument("Failed to parse atom token: " + to_add);
    }
    auto *qnode = one->getAtomWithIdx(0)->getQuery();
    auto *src = one->getAtomWithIdx(0);
    ni = merged_template->addAtom(src->copy(), true, true);
    if (qnode) {
      merged_template->getAtomWithIdx(ni)->setQuery(qnode->copy());
    }
    unmapped_idx_map.emplace(source_atom_ptr, static_cast<int>(ni));
  } else {
    ni = static_cast<unsigned>(it->second);
  }

  if (idx_to_add_to >= 0 &&
      merged_template->getBondBetweenAtoms(idx_to_add_to, ni) == nullptr) {
    merged_template->addBond(idx_to_add_to, ni, bond_type_to_add);
  }
  return ni;
}

static void add_atom_to_merged_template(
    const std::string &to_add, int idx_to_add_to, BondType bond_type_to_add,
    std::unique_ptr<RDKit::RWMol> &merged_template,
    std::unordered_map<int, int> &merged_map_to_id) {
  std::unique_ptr<RDKit::ROMol> one(RDKit::SmartsToMol(to_add));
  if (!one || one->getNumAtoms() != 1) {
    throw std::invalid_argument("Failed to parse atom token: " + to_add);
  }
  auto *qnode = one->getAtomWithIdx(0)->getQuery();
  auto *src = one->getAtomWithIdx(0);

  unsigned ni;
  if (merged_map_to_id.count(src->getAtomMapNum()) == 0) {
    ni = merged_template->addAtom(src->copy(), true, true);
    if (qnode) {
      merged_template->getAtomWithIdx(ni)->setQuery(qnode->copy());
    }
    merged_map_to_id[src->getAtomMapNum()] = static_cast<int>(ni);
  } else {
    ni = static_cast<unsigned>(merged_map_to_id[src->getAtomMapNum()]);
  }

  if (idx_to_add_to == -1) {
    return;
  }
  if (merged_template->getBondBetweenAtoms(idx_to_add_to, ni) != nullptr) {
    return;
  }
  merged_template->addBond(idx_to_add_to, ni, bond_type_to_add);
}

static std::string dump_template(
    std::vector<std::unique_ptr<RDKit::RWMol>> &merged_reactant_templates,
    std::unique_ptr<RDKit::RWMol> &merged_product_template) {
  std::ostringstream oss;
  bool added = false;
  for (size_t i = 0; i < merged_reactant_templates.size(); ++i) {
    std::string smarts = RDKit::MolToSmarts(*merged_reactant_templates[i]);
    if (smarts.empty()) {
      continue;
    }

    if (i > 0 && added) {
      oss << ".";
    }
    if (smarts.find('.') != std::string::npos) {
      oss << "(" << smarts << ")";
    } else {
      oss << smarts;
    }
    added = true;
  }
  oss << ">>" << RDKit::MolToSmarts(*merged_product_template);
  return oss.str();
}

template <typename MolPtrVec>
static std::vector<RadiusTemplate> fragments_for_mols(
    const MolPtrVec &reactants, const MolPtrVec &products, unsigned maxRadius,
    bool verbose, const std::vector<std::string> &primitiveList) {
  std::unordered_map<int, std::pair<RDKit::Atom *, MolPtr>> product_keys;
  for (const MolPtr &kv : products) {
    if (!kv) {
      continue;
    }
    for (auto *atom : kv->atoms()) {
      if (atom->getAtomMapNum() > 0) {
        product_keys[atom->getAtomMapNum()] = {atom, kv};
      }
    }
  }

  std::unordered_map<int, std::pair<RDKit::Atom *, MolPtr>> reactant_keys;

  std::vector<std::unique_ptr<RDKit::RWMol>> merged_reactant_templates(
      reactants.size());
  for (size_t i = 0; i < reactants.size(); ++i) {
    merged_reactant_templates[i] = std::unique_ptr<RDKit::RWMol>(new RDKit::RWMol());
  }
  std::unique_ptr<RDKit::RWMol> merged_product_template(new RDKit::RWMol());

  std::unordered_map<int, int> merged_reactant_map_to_id;
  std::unordered_map<int, int> merged_product_map_to_id;
  std::unordered_map<int, int> amap_to_reactant_idx;
  std::unordered_map<const RDKit::Atom *, int> unmapped_reactant_idx;
  std::unordered_map<const RDKit::Atom *, int> unmapped_product_idx;

  int reactant_idx = 0;
  for (const MolPtr &kv : reactants) {
    if (!kv) {
      reactant_idx++;
      continue;
    }
    for (auto *atom : kv->atoms()) {
      const int amap = atom->getAtomMapNum();
      if (amap <= 0) {
        continue;
      }
      amap_to_reactant_idx[amap] = reactant_idx;
      reactant_keys[amap] = {atom, kv};
    }
    reactant_idx++;
  }

  std::unordered_set<int> changed_amaps;
  auto build_neighbor_signature = [](const RDKit::ROMol &mol,
                                     const RDKit::Atom *a) {
    std::unordered_set<std::string> sig;
    RDKit::ROMol::ADJ_ITER it, end;
    boost::tie(it, end) = mol.getAtomNeighbors(a);
    for (; it != end; ++it) {
      const RDKit::Atom *nb = mol.getAtomWithIdx(*it);
      const int nb_map = nb->getAtomMapNum();
      if (nb_map <= 0) {
        continue;
      }
      const RDKit::Bond *b = mol.getBondBetweenAtoms(a->getIdx(), nb->getIdx());
      int bt = static_cast<int>(b ? b->getBondType() : RDKit::Bond::UNSPECIFIED);
      sig.insert(std::to_string(nb_map) + "#" + std::to_string(bt));
    }
    return sig;
  };

  for (const auto &kvp : product_keys) {
    const int amap = kvp.first;
    if (reactant_keys.count(amap) == 0) {
      continue;
    }
    const RDKit::Atom *p_atom = kvp.second.first;
    const MolPtr p_mol = kvp.second.second;
    const RDKit::Atom *r_atom = reactant_keys[amap].first;
    const MolPtr r_mol = reactant_keys[amap].second;

    bool changed = false;
    if (p_atom->getFormalCharge() != r_atom->getFormalCharge()) {
      changed = true;
    }
    if (!changed && p_atom->getIsAromatic() != r_atom->getIsAromatic()) {
      changed = true;
    }
    if (!changed && p_atom->getTotalNumHs() != r_atom->getTotalNumHs()) {
      changed = true;
    }
    if (!changed) {
      auto p_sig = build_neighbor_signature(*p_mol, p_atom);
      auto r_sig = build_neighbor_signature(*r_mol, r_atom);
      if (p_sig != r_sig) {
        changed = true;
      }
    }
    if (changed) {
      changed_amaps.insert(amap);
    }
  }

  if (changed_amaps.empty()) {
    for (const auto &kvp : product_keys) {
      const int amap = kvp.first;
      if (reactant_keys.count(amap) == 0) {
        continue;
      }
      const RDKit::Atom *p_atom = kvp.second.first;
      const MolPtr p_mol = kvp.second.second;
      const RDKit::Atom *r_atom = reactant_keys[amap].first;
      const MolPtr r_mol = reactant_keys[amap].second;

      PathType reactant_env =
          RDKit::findAtomEnvironmentOfRadiusN(*r_mol, 1, r_atom->getIdx());
      PathType product_env =
          RDKit::findAtomEnvironmentOfRadiusN(*p_mol, 1, p_atom->getIdx());
      std::unique_ptr<RDKit::ROMol> reactant_env_mol(
          RDKit::Subgraphs::pathToSubmol(*r_mol, reactant_env, false));
      std::unique_ptr<RDKit::ROMol> product_env_mol(
          RDKit::Subgraphs::pathToSubmol(*p_mol, product_env, false));
      std::string r_sm = RDKit::MolToSmiles(*reactant_env_mol, true, false, -1, true);
      std::string p_sm = RDKit::MolToSmiles(*product_env_mol, true, false, -1, true);
      if (r_sm != p_sm) {
        changed_amaps.insert(amap);
      }
    }
  }

  for (const int amap : changed_amaps) {
    const RDKit::Atom *r_atom = reactant_keys[amap].first;
    const RDKit::Atom *p_atom = product_keys[amap].first;
    const int ridx = amap_to_reactant_idx[amap];
    const std::string parsed_r_atom =
        build_atom_primitive_token(r_atom, &r_atom->getOwningMol(),
                                   primitiveList, true);
    const std::string parsed_p_atom =
        build_atom_primitive_token(p_atom, &p_atom->getOwningMol(),
                                   primitiveList, true);
    add_atom_to_merged_template(parsed_r_atom, -1, BondType::SINGLE,
                                merged_reactant_templates[ridx],
                                merged_reactant_map_to_id);
    add_atom_to_merged_template(parsed_p_atom, -1, BondType::SINGLE,
                                merged_product_template,
                                merged_product_map_to_id);
  }

  std::vector<RadiusTemplate> results;
  results.emplace_back(0, dump_template(merged_reactant_templates,
                                        merged_product_template));

  if (verbose) {
    std::cout << std::get<1>(results.back()) << std::endl;
  }

  std::unordered_set<int> reactant_frontier;
  std::unordered_set<int> reactant_processed;
  for (const int amap : changed_amaps) {
    reactant_frontier.insert(amap);
  }

  int radius = 1;
  std::unordered_set<int> processed_atoms;
  while (radius <= static_cast<int>(maxRadius)) {
    std::unordered_set<int> next_reactant_frontier;
    for (const int amap_seed : reactant_frontier) {
      if (reactant_processed.count(amap_seed) > 0) {
        continue;
      }
      if (reactant_keys.count(amap_seed) == 0) {
        continue;
      }

      const RDKit::Atom *r_atom_seed = reactant_keys[amap_seed].first;
      const MolPtr r_mol_seed = reactant_keys[amap_seed].second;
      const int r_template_idx_seed = amap_to_reactant_idx[amap_seed];

      RDKit::ROMol::ADJ_ITER r_it, r_end;
      boost::tie(r_it, r_end) = r_mol_seed->getAtomNeighbors(r_atom_seed);
      for (; r_it != r_end; ++r_it) {
        const unsigned rv = *r_it;
        RDKit::Atom *r_neighbor_atom = r_mol_seed->getAtomWithIdx(rv);
        const int neighbor_amap = r_neighbor_atom->getAtomMapNum();
        const RDKit::Bond *r_b = r_mol_seed->getBondBetweenAtoms(
            r_neighbor_atom->getIdx(), r_atom_seed->getIdx());
        const std::string parsed_nr_atom =
            build_atom_primitive_token(r_neighbor_atom, r_mol_seed.get(),
                                       primitiveList, true);

        auto it_anchor_react = merged_reactant_map_to_id.find(amap_seed);
        if (neighbor_amap <= 0) {
          int anchor_idx = (it_anchor_react != merged_reactant_map_to_id.end())
                               ? it_anchor_react->second
                               : -1;
          add_unmapped_atom_and_bond(
              parsed_nr_atom, anchor_idx,
              r_b ? r_b->getBondType() : BondType::SINGLE,
              merged_reactant_templates[r_template_idx_seed],
              unmapped_reactant_idx, r_neighbor_atom);
          continue;
        }
        if (reactant_keys.count(neighbor_amap) == 0) {
          continue;
        }

        const int neighbor_template_idx = amap_to_reactant_idx[neighbor_amap];

        if (neighbor_template_idx == r_template_idx_seed) {
          if (it_anchor_react != merged_reactant_map_to_id.end() && r_b) {
            add_atom_to_merged_template(
                parsed_nr_atom, it_anchor_react->second, r_b->getBondType(),
                merged_reactant_templates[r_template_idx_seed],
                merged_reactant_map_to_id);
          } else {
            add_atom_to_merged_template(
                parsed_nr_atom, -1, BondType::SINGLE,
                merged_reactant_templates[r_template_idx_seed],
                merged_reactant_map_to_id);
          }
        } else {
          add_atom_to_merged_template(
              parsed_nr_atom, -1, BondType::SINGLE,
              merged_reactant_templates[neighbor_template_idx],
              merged_reactant_map_to_id);
        }

        if (reactant_processed.count(neighbor_amap) == 0) {
          next_reactant_frontier.insert(neighbor_amap);
        }
      }
    }

    int num_atoms = static_cast<int>(merged_product_template->getNumAtoms());
    for (int i = 0; i < num_atoms; ++i) {
      const RDKit::Atom *p_template_atom = merged_product_template->getAtomWithIdx(i);
      const int amap = p_template_atom->getAtomMapNum();
      if (processed_atoms.count(amap) > 0) {
        continue;
      }
      processed_atoms.insert(amap);

      const RDKit::Atom *p_atom = product_keys[amap].first;
      const MolPtr p_mol = product_keys[amap].second;
      const RDKit::Atom *r_atom = reactant_keys[amap].first;
      const MolPtr r_mol = reactant_keys[amap].second;
      const int r_atom_idx = r_atom->getIdx();
      const int r_template_idx = amap_to_reactant_idx[amap];

      RDKit::ROMol::ADJ_ITER it, end;
      boost::tie(it, end) = p_mol->getAtomNeighbors(p_atom);
      for (; it != end; ++it) {
        const unsigned u = *it;
        RDKit::Atom *neighbor_atom = p_mol->getAtomWithIdx(u);
        int neighbor_amap = neighbor_atom->getAtomMapNum();
        if (neighbor_amap <= 0) {
          continue;
        }
        if (reactant_keys.count(neighbor_amap) == 0) {
          continue;
        }
        if (reactant_processed.count(amap) == 0 &&
            reactant_frontier.count(amap) == 0 &&
            reactant_frontier.count(neighbor_amap) == 0) {
          continue;
        }
        RDKit::Atom *r_neighbor_atom = reactant_keys[neighbor_amap].first;

        std::string parsed_np_atom =
            build_atom_primitive_token(neighbor_atom, p_mol.get(), primitiveList,
                                       true);
        std::string parsed_nr_atom =
            build_atom_primitive_token(r_neighbor_atom, r_mol.get(),
                                       primitiveList, true);

        const RDKit::Bond *b = p_mol->getBondBetweenAtoms(u, p_atom->getIdx());
        auto it_anchor_prod = merged_product_map_to_id.find(amap);
        if (it_anchor_prod != merged_product_map_to_id.end() && b) {
          if (neighbor_amap > 0) {
            add_atom_to_merged_template(parsed_np_atom, it_anchor_prod->second,
                                        b->getBondType(), merged_product_template,
                                        merged_product_map_to_id);
          } else {
            add_unmapped_atom_and_bond(parsed_np_atom, it_anchor_prod->second,
                                       b->getBondType(), merged_product_template,
                                       unmapped_product_idx, neighbor_atom);
          }
        }

        if (neighbor_amap > 0) {
          const int neighbor_template_idx = amap_to_reactant_idx[neighbor_amap];
          if (neighbor_template_idx == r_template_idx) {
            const RDKit::Bond *r_b =
                r_mol->getBondBetweenAtoms(r_neighbor_atom->getIdx(), r_atom_idx);
            auto it_anchor_react = merged_reactant_map_to_id.find(amap);
            if (it_anchor_react != merged_reactant_map_to_id.end() && r_b) {
              add_atom_to_merged_template(
                  parsed_nr_atom, it_anchor_react->second, r_b->getBondType(),
                  merged_reactant_templates[r_template_idx],
                  merged_reactant_map_to_id);
            } else {
              add_atom_to_merged_template(
                  parsed_nr_atom, -1, BondType::SINGLE,
                  merged_reactant_templates[r_template_idx],
                  merged_reactant_map_to_id);
            }
          } else {
            add_atom_to_merged_template(
                parsed_nr_atom, -1, BondType::SINGLE,
                merged_reactant_templates[neighbor_template_idx],
                merged_reactant_map_to_id);
          }
        } else {
          auto it_anchor_react = merged_reactant_map_to_id.find(amap);
          const RDKit::Bond *r_b =
              r_mol->getBondBetweenAtoms(r_neighbor_atom->getIdx(), r_atom_idx);
          int anchor_idx = (it_anchor_react != merged_reactant_map_to_id.end())
                               ? it_anchor_react->second
                               : -1;
          add_unmapped_atom_and_bond(
              parsed_nr_atom, anchor_idx,
              r_b ? r_b->getBondType() : BondType::SINGLE,
              merged_reactant_templates[r_template_idx], unmapped_reactant_idx,
              r_neighbor_atom);
        }
      }
    }

    for (const int a_done : reactant_frontier) {
      reactant_processed.insert(a_done);
    }
    reactant_frontier = std::move(next_reactant_frontier);

    results.emplace_back(radius, dump_template(merged_reactant_templates,
                                               merged_product_template));
    if (verbose) {
      std::cout << std::get<1>(results.back()) << std::endl;
    }

    radius++;
  }

  return results;
}

}  // namespace

std::vector<RadiusTemplate> extract_single_root_template(
    RDKit::ChemicalReaction *rxn, const std::vector<std::string> &primitiveList,
    unsigned maxRadius, bool verbose) {
  if (!rxn) {
    throw std::invalid_argument("Reaction pointer is null");
  }
  if (normalize_primitive_list(primitiveList).empty()) {
    throw std::invalid_argument("Primitive list is empty");
  }

  const auto &reactants = rxn->getReactants();
  const auto &products = rxn->getProducts();
  for (auto i : reactants) {
    if (i) {
      i->updatePropertyCache(false);
    }
  }
  for (auto i : products) {
    if (i) {
      i->updatePropertyCache(false);
    }
  }

  return fragments_for_mols(reactants, products, maxRadius, verbose,
                            primitiveList);
}

std::vector<RadiusTemplate> extract_single_root_template(
    const std::string &reactionText,
    const std::vector<std::string> &primitiveList,
    unsigned maxRadius, bool verbose, bool useSmiles) {
  if (reactionText.empty()) {
    throw std::invalid_argument("Reaction text is empty");
  }

  std::unique_ptr<RDKit::ChemicalReaction> rxn(
      RDKit::RxnSmartsToChemicalReaction(reactionText, nullptr, useSmiles));
  if (!rxn) {
    throw std::invalid_argument("Invalid reaction text");
  }

  return extract_single_root_template(rxn.get(), primitiveList, maxRadius,
                                      verbose);
}

}  // namespace atom_typer
