
// _rxn_backend_min.cpp
//
// Build (example):
//   cl /O2 /std:c++20 /MD /EHsc /DBOOST_PYTHON_DYN_LINK ^
//     /I"%CONDA_PREFIX%\Include" /I"%CONDA_PREFIX%\Library\include" ^
//     /I"%RDBASE%\Code" ^
//     /c _rxn_backend_min.cpp && ^
//   link /DLL /OUT:_rxn_backend_min.pyd _rxn_backend_min.obj ^
//     rdkitchemreactions.lib rdkitgraphmol.lib rdkitrdgeneral.lib ^
//     rdkitfingerprint.lib rdkitSubstruct.lib
//     boost_python311-vc143-mt-x64-1_85.lib
//
// Python usage:
//   from _rxn_backend_min import ReactionEngine
//   eng = ReactionEngine()
//   eng.load_reactions(reaction_list)
//   eng.load_molecules_into_pool(mol_list)   # sanitized, deduped
//   eng.run_pending(); eng.add_products_to_pool()

#include <algorithm>
#include <atomic>
#include <chrono>
#include <csignal>
#include <cstdio> // setvbuf
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <boost/noncopyable.hpp>
#include <boost/python.hpp>

// Simple interrupt flag for Ctrl+C handling
std::atomic<bool> g_interrupt{false};

void signal_handler(int) { g_interrupt = true; }

#include <DataStructs/BitVects.h>
#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemReactions/ReactionRunner.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <RDGeneral/types.h>

namespace bp = boost::python;
using namespace RDKit;

using VectMatchVectType = std::vector<MatchVectType>;

class ReactionEngine {
public:
  using MolPtr = ROMOL_SPTR;

  ReactionEngine(const ReactionEngine &) = delete;
  ReactionEngine &operator=(const ReactionEngine &) = delete;
  ReactionEngine(ReactionEngine &&) = default;
  ReactionEngine &operator=(ReactionEngine &&) = default;

  ReactionEngine() { std::signal(SIGINT, signal_handler); }

  // Keep signature; ignore BVs (accepted for compatibility)
  void load_reactions(const bp::list &reaction_list,
                      const bp::list & /*per_reaction_per_role_bvs*/) {
    reactions_.clear();
    reactant_templates_.clear();
    reaction_to_template_indices_.clear();

    const int n_rxns = bp::len(reaction_list);
    reactions_.reserve(n_rxns);

    for (int i = 0; i < n_rxns; ++i) {
      bp::extract<RDKit::ChemicalReaction *> rxn_ex(reaction_list[i]);
      if (!rxn_ex.check() || !rxn_ex()) {
        throw std::runtime_error(
            "load_reactions: element is not a ChemicalReaction");
      }
      auto *rxn = rxn_ex();
      rxn->initReactantMatchers(false);
      reactions_.push_back(rxn);
    }

    std::unordered_map<std::string, int> idx_of;
    idx_of.reserve(64);
    for (int ridx = 0; ridx < static_cast<int>(reactions_.size()); ++ridx) {
      auto *rxn = reactions_[ridx];
      const MOL_SPTR_VECT &rtemps = rxn->getReactants();

      int ta = -1, tb = -1;
      for (size_t pos = 0; pos < rtemps.size(); ++pos) {
        const auto &tmpl = rtemps[pos];
        const std::string smarts = RDKit::MolToSmarts(*tmpl);
        int t_idx;
        auto it = idx_of.find(smarts);
        if (it == idx_of.end()) {
          t_idx = static_cast<int>(reactant_templates_.size());
          idx_of.emplace(smarts, t_idx);
          reactant_templates_.emplace_back(tmpl, smarts);
        } else {
          t_idx = it->second;
        }
        reactant_templates_[t_idx].reaction_positions.emplace_back(
            ridx, static_cast<int>(pos));
        if (pos == 0)
          ta = t_idx;
        else if (pos == 1)
          tb = t_idx;
      }
      reaction_to_template_indices_[ridx] = {ta, tb};
    }
  }

  // Keep signature; ignore BVs
  void load_molecules_into_pool(const bp::list &smiles_list) {
    const int N = bp::len(smiles_list);
    molecule_pool_.clear();
    seen_reaction_combos_.clear();
    lineage_map_.clear();
    seen_smiles_.clear();
    molecule_pool_.reserve(N);

    for (int i = 0; i < N; ++i) {
      boost::shared_ptr<RDKit::ROMol> m;
      try {
        m = bp::extract<boost::shared_ptr<RDKit::ROMol>>(smiles_list[i]);
      } catch (...) {
        continue;
      }
      // if (sm.empty())
      //   continue;

      try {
        // Parse with default sanitize=true.
        // ROMOL_SPTR m(SmilesToMol(sm));
        // if (!m)
        // continue;

        // Some RDKit builds only expose sanitizeMol(RWMol&): sanitize via a
        // copy.
        // {
        //   RWMol rw(*m);
        //   RDKit::MolOps::sanitizeMol(rw);
        //   m.reset(new ROMol(rw));
        // }

        // RDKit::MolOps::assignStereochemistry(*m, true, true);
        m->updatePropertyCache(
            true); // strict: ensure implicit valence is computed

        const std::string can = MolToSmiles(*m, true, false, -1, true);
        if (seen_smiles_.insert(can).second)
          molecule_pool_.emplace_back(std::move(m));
      } catch (const std::exception &e) {
        std::cerr << "[load_molecules_into_pool] skipping \"" << "i"
                  << "\": " << e.what() << "\n";
      } catch (...) {
        std::cerr << "[load_molecules_into_pool] skipping \"" << "i"
                  << "\": unknown RDKit error\n";
      }
    }
    int status = perform_template_matching();
    if (status == -1) {
      std::cerr << "[load_molecules_into_pool] template matching ended early\n";
      return;
    }
    generate_new_combinations();
  }

  void run_pending(bool verbose = false) { run_pending_serial(verbose); }

  void run_one_atom_level_reaction(int reaction_index, bool verbose = false) {
    if (reaction_index < 0 ||
        reaction_index >=
            static_cast<int>(pending_atom_level_reactions_.size())) {
      std::cout << "Reaction " << reaction_index << " is not in pending\n";
      return;
    }
    const AtomLevelReaction &rxn =
        pending_atom_level_reactions_[reaction_index];
    auto res = process_one_atom_level_reaction(rxn);
    if (verbose) {
      std::cout << "Running reaction " << reaction_index << " of "
                << pending_atom_level_reactions_.size() << "\n";
      std::cout << "SMARTS:"
                << ChemicalReactionToRxnSmarts(
                       *reactions_[rxn.reaction_template_index])
                << "\n";
      std::cout << "Reactants:"
                << MolToSmiles(*molecule_pool_[rxn.combo_key.a_index]) << " "
                << (rxn.combo_key.b_index >= 0
                        ? MolToSmiles(*molecule_pool_[rxn.combo_key.b_index])
                        : "")
                << "\n";
      std::cout << "Atom indices:" << rxn.atom_indices[0].size() << " "
                << (rxn.atom_indices.size() > 1 ? rxn.atom_indices[1].size()
                                                : 0)
                << "\n";
      if (!std::get<1>(res).empty()) {
        std::string res_smiles = MolToSmiles(*std::get<1>(res)[0]);
        std::cout << "Result count: " << std::get<1>(res).size() << "\n";
        std::cout << "Result smiles: " << res_smiles << "\n";
        std::cout << "Reaction index: " << std::get<2>(res) << "\n";
      }
      std::cout << "--------------------------------\n";
    }
  }

  void run_pending_serial(bool verbose = false) {
    if (!timer_started_) {
      start_time_ = std::chrono::steady_clock::now();
      timer_started_ = true;
    }
    pending_results_.clear();
    int total = static_cast<int>(pending_atom_level_reactions_.size());

    if (pending_atom_level_reactions_.size() > 30000) {
      std::cerr << "Warning: Pending results size ("
                << pending_atom_level_reactions_.size()
                << ") too large. Limiting.\n";
      total = 30000;
      // return;
    }

    for (int i = 0; i < total; ++i) {
      // Check for interrupt signal
      if (g_interrupt.load()) {
        std::cerr << "\n[C++] Interrupted after " << i << " reactions\n";
        break;
      }

      // stop on max runtime
      current_time_ = std::chrono::duration_cast<std::chrono::duration<double>>(
                          std::chrono::steady_clock::now() - start_time_)
                          .count();
      if (current_time_ > max_runtime_) {
        std::cerr << "\n[C++] run_pending_serial interrupted (max runtime)\n";
        break;
      }

      const AtomLevelReaction &rxn = pending_atom_level_reactions_[i];
      auto res = process_one_atom_level_reaction(rxn);
      if (!std::get<1>(res).empty())
        pending_results_.push_back(res);
      if (i % 1000 == 0) {
        if (verbose) {
          std::cout << "> Running reaction " << i << " of " << total << " ("
                    << current_time_ << "s)\n";
          std::cout << "> Pending results size: " << pending_results_.size()
                    << " (" << current_time_ << "s)\n";
        }
      }
    }
  }

  void add_products_to_pool(bool parallelize = false, int num_threads = 0,
                            bool verbose = false) {
    if (pending_results_.size() > 50001) {
      std::cout << "Warning: Pending results size (" << pending_results_.size()
                << ") too large. Limiting.\n";
      pending_results_.resize(50001);
      return;
    }

    if (verbose) {
      std::cout << "pending_results_.size() = " << pending_results_.size()
                << std::endl;
    }

    int result_idx = 0;
    int pend_idx = 0;
    for (auto &result_pair : pending_results_) {
      const ComboKey combo_key = std::get<0>(result_pair);
      std::vector<MolPtr> &product_set = std::get<1>(result_pair);
      const int reaction_idx = std::get<2>(result_pair);

      for (MolPtr &product : product_set) {
        try {
          // Skip sanitization for now to avoid ROMol constructor issues
          // TODO: Re-enable sanitization when linking issues are resolved
          product->updatePropertyCache(false);
          const std::string sm = MolToSmiles(*product, true, false, -1, true);
          if (seen_smiles_.insert(sm).second) {
            int new_idx = static_cast<int>(molecule_pool_.size());
            molecule_pool_.push_back(std::move(product));
            lineage_map_.emplace(new_idx, LineageNode(combo_key, reaction_idx));
          }
        } catch (const std::exception &e) {
          std::cerr << "[add_products_to_pool] dropping invalid product: "
                    << e.what() << "\n";
        } catch (...) {
          std::cerr << "[add_products_to_pool] dropping invalid product: "
                       "unknown RDKit error\n";
        }
      }
      if (verbose) {
        if (result_idx % 10000 == 0) {
          std::cout << ">> completed " << result_idx << " product sets of "
                    << pending_results_.size() << std::endl;
        }
      }

      result_idx++;
    }
    pending_results_.clear();
    if (verbose) {
      std::cout << "..molecule_pool_.size() = " << molecule_pool_.size()
                << std::endl;
    }
    perform_template_matching(parallelize, num_threads, verbose);
    if (verbose) {
      std::cout << "...template_match_matrix_.size() = "
                << template_match_matrix_.size() << std::endl;
    }
    generate_new_combinations(verbose);
  }

  // ---- Introspection / getters ----
  bp::dict get_molecule_lineage(int mol_index) {
    bp::dict lineage_info;
    auto it = lineage_map_.find(mol_index);
    if (it != lineage_map_.end()) {
      const LineageNode &node = it->second;
      bp::list reactants;
      if (node.reactant_combo.a_index != -1)
        reactants.append(node.reactant_combo.a_index);
      if (node.reactant_combo.b_index != -1)
        reactants.append(node.reactant_combo.b_index);
      lineage_info["reactants"] = reactants;
      lineage_info["reaction"] = node.reaction_index;
    } else {
      lineage_info["reactants"] = bp::list();
      lineage_info["reaction"] = -1;
    }
    return lineage_info;
  }

  bp::list get_full_lineage() {
    bp::list full_tree;
    for (const auto &pair : lineage_map_) {
      int product_index = pair.first;
      const LineageNode &node = pair.second;
      bp::dict node_info;
      node_info["product"] = product_index;
      bp::list reactants;
      if (node.reactant_combo.a_index != -1)
        reactants.append(node.reactant_combo.a_index);
      if (node.reactant_combo.b_index != -1)
        reactants.append(node.reactant_combo.b_index);
      node_info["reactants"] = reactants;
      node_info["reaction"] = node.reaction_index;
      full_tree.append(node_info);
    }
    return full_tree;
  }

  int get_pool_size() const { return static_cast<int>(molecule_pool_.size()); }

  bp::list get_all_molecules_as_smiles() const {
    bp::list result;
    for (const auto &mol : molecule_pool_) {
      try {
        result.append(MolToSmiles(*mol));
      } catch (...) {
        result.append("");
      }
    }
    return result;
  }

  int get_num_reactant_templates() const {
    return static_cast<int>(reactant_templates_.size());
  }

  int get_num_pending_atom_level_reactions() const {
    return static_cast<int>(pending_atom_level_reactions_.size());
  }

  bp::str get_molecule_smiles_from_pool(int mol_index) const {
    if (mol_index >= 0 && mol_index < static_cast<int>(molecule_pool_.size()))
      return bp::str(MolToSmiles(*molecule_pool_[mol_index]));
    return bp::str("");
  }

  bp::dict get_reactant_template_info(int template_index) const {
    bp::dict info;
    if (template_index >= 0 &&
        template_index < static_cast<int>(reactant_templates_.size())) {
      const auto &templ = reactant_templates_[template_index];
      info["smarts"] = templ.template_smarts;
      bp::list positions;
      for (const auto &pos : templ.reaction_positions) {
        bp::dict p;
        p["reaction_index"] = pos.first;
        p["reactant_position"] = pos.second;
        positions.append(p);
      }
      try {
        info["smiles"] = MolToSmiles(*templ.template_mol);
      } catch (...) {
        info["smiles"] = "";
      }
    }
    return info;
  }

  bp::list get_template_matches(int template_index, int molecule_index) const {
    bp::list matches;
    if (template_index >= 0 &&
        template_index < static_cast<int>(template_match_matrix_.size()) &&
        molecule_index >= 0 &&
        molecule_index <
            static_cast<int>(template_match_matrix_[template_index].size())) {
      const VectMatchVectType &mol_matches =
          template_match_matrix_[template_index][molecule_index];
      bp::list atom_indices;
      for (const auto &match_vec : mol_matches) {
        for (const auto &pair : match_vec) {
          bp::dict ai;
          ai["mol_atom_idx"] = pair.first;
          ai["pattern_atom_idx"] = pair.second;
          atom_indices.append(ai);
        }
      }
      matches.append(atom_indices);
    }
    return matches;
  }

  bp::dict get_all_template_matches() const {
    bp::dict all_matches;
    for (size_t t = 0; t < template_match_matrix_.size(); ++t) {
      bp::dict template_matches;
      for (size_t m = 0; m < template_match_matrix_[t].size(); ++m) {
        bp::list matches;
        const VectMatchVectType &mol_matches = template_match_matrix_[t][m];
        for (const auto &match_vec : mol_matches)
          for (const auto &pair : match_vec) {
            bp::dict ai;
            ai["mol_atom_idx"] = pair.first;
            ai["pattern_atom_idx"] = pair.second;
            matches.append(ai);
          }
        template_matches[static_cast<int>(m)] = matches;
      }
      all_matches[static_cast<int>(t)] = template_matches;
    }
    return all_matches;
  }

  bp::list get_pending_combinations_info() const {
    bp::list combinations_info;
    for (size_t i = 0; i < pending_atom_level_reactions_.size(); ++i) {
      const auto &par = pending_atom_level_reactions_[i];
      bp::dict combo_info;
      bp::list reactant_indices, reactant_smiles, atom_indices_list;

      reactant_indices.append(par.combo_key.a_index);
      reactant_smiles.append(
          MolToSmiles(*molecule_pool_[par.combo_key.a_index]));
      bp::list this_atom_indices_a;
      for (const auto &pr : par.atom_indices[0]) {
        bp::dict ai;
        ai["mol_atom_idx"] = pr.first;
        ai["pattern_atom_idx"] = pr.second;
        this_atom_indices_a.append(ai);
      }
      atom_indices_list.append(this_atom_indices_a);

      if (par.atom_indices.size() > 1) {
        reactant_indices.append(par.combo_key.b_index);
        reactant_smiles.append(
            MolToSmiles(*molecule_pool_[par.combo_key.b_index]));
        bp::list this_atom_indices_b;
        for (const auto &pr : par.atom_indices[1]) {
          bp::dict ai;
          ai["mol_atom_idx"] = pr.first;
          ai["pattern_atom_idx"] = pr.second;
          this_atom_indices_b.append(ai);
        }
        atom_indices_list.append(this_atom_indices_b);
      }

      combo_info["reactant_indices"] = reactant_indices;
      combo_info["reactant_smiles"] = reactant_smiles;
      combo_info["reaction_smarts_template"] =
          ChemicalReactionToRxnSmarts(*reactions_[par.reaction_template_index]);
      combo_info["reaction_template_index"] = par.reaction_template_index;
      combo_info["atom_indices"] = atom_indices_list;

      combinations_info.append(combo_info);
    }
    return combinations_info;
  }

  bp::list get_atom_level_reactions() const {
    bp::list atom_reactions;
    for (const auto &atom_rxn : pending_atom_level_reactions_) {
      bp::dict reaction_info;
      bp::list mol_indices;
      for (int mol_idx :
           {atom_rxn.combo_key.a_index, atom_rxn.combo_key.b_index})
        mol_indices.append(mol_idx);
      reaction_info["molecule_indices"] = mol_indices;

      bp::list atom_indices_list;
      for (const auto &atom_indices : atom_rxn.atom_indices) {
        bp::list atom_list;
        for (const auto &pair : atom_indices) {
          bp::dict ai;
          ai["mol_atom_idx"] = pair.first;
          ai["pattern_atom_idx"] = pair.second;
          atom_list.append(ai);
        }
        atom_indices_list.append(atom_list);
      }
      reaction_info["atom_indices"] = atom_indices_list;
      reaction_info["reaction_template"] = atom_rxn.reaction_template_index;
      atom_reactions.append(reaction_info);
    }
    return atom_reactions;
  }

  bp::list get_reaction_forming_product(bp::str product_smiles) {
    std::string product_smiles_str = bp::extract<std::string>(product_smiles);
    for (int i = 0; i < static_cast<int>(molecule_pool_.size()); ++i) {
      if (MolToSmiles(*molecule_pool_[i]) == product_smiles_str) {
        auto it = lineage_map_.find(i);
        if (it == lineage_map_.end())
          break;
        const auto &lineage = it->second;
        bp::list out;
        if (lineage.reactant_combo.b_index < 0) {
          out.append(
              MolToSmiles(*molecule_pool_[lineage.reactant_combo.a_index]));
        } else {
          out.append(
              MolToSmiles(*molecule_pool_[lineage.reactant_combo.a_index]) +
              "." +
              MolToSmiles(*molecule_pool_[lineage.reactant_combo.b_index]));
        }
        out.append(lineage.reaction_index);
        return out;
      }
    }
    return bp::list();
  }

  bp::str get_reaction_smarts_from_index(int reaction_index) const {
    return bp::str(ChemicalReactionToRxnSmarts(*reactions_[reaction_index]));
  }

private:
  // (simplified) no external limit configuration retained
  // runtime timer state
  int max_runtime_ = 10;
  bool timer_started_ = false;
  std::chrono::steady_clock::time_point start_time_;
  double current_time_ = 0.0;

  struct ComboKey {
    int a_index, b_index;
    ComboKey(int a = -1, int b = -1) : a_index(a), b_index(b) {}
    bool operator==(const ComboKey &o) const noexcept {
      return a_index == o.a_index && b_index == o.b_index;
    }
  };
  struct ComboKeyHash {
    std::size_t operator()(const ComboKey &k) const noexcept {
      auto h1 = std::hash<int>()(k.a_index);
      auto h2 = std::hash<int>()(k.b_index);
      return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
    }
  };

  struct LineageNode {
    ComboKey reactant_combo;
    int reaction_index;
    LineageNode(const ComboKey &combo, int rxn_idx)
        : reactant_combo(combo), reaction_index(rxn_idx) {}
  };

  struct AtomLevelReaction {
    ComboKey combo_key;
    std::vector<MatchVectType> atom_indices;
    int reaction_template_index;
    AtomLevelReaction(const ComboKey &combo,
                      const std::vector<MatchVectType> &atoms, int rxn_idx)
        : combo_key(combo), atom_indices(atoms),
          reaction_template_index(rxn_idx) {}
  };

  struct ReactantTemplate {
    MolPtr template_mol;
    std::string template_smarts;
    std::vector<std::pair<int, int>> reaction_positions; // (rxn_idx, role)
    std::unique_ptr<ExplicitBitVect> fp; // Morgan r=2, nBits=2048
    ReactantTemplate(MolPtr mol, const std::string &smarts)
        : template_mol(mol), template_smarts(smarts), fp(nullptr) {}
  };

  struct ReactionComboKey {
    int reaction_index;
    ComboKey combo_key;
    bool operator==(const ReactionComboKey &o) const noexcept {
      return reaction_index == o.reaction_index &&
             combo_key.a_index == o.combo_key.a_index &&
             combo_key.b_index == o.combo_key.b_index;
    }
  };
  struct ReactionComboKeyHash {
    std::size_t operator()(const ReactionComboKey &k) const noexcept {
      std::size_t h1 = std::hash<int>()(k.reaction_index);
      std::size_t h2 = std::hash<int>()(k.combo_key.a_index);
      std::size_t h3 = std::hash<int>()(k.combo_key.b_index);
      std::size_t h12 =
          h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
      return h12 ^ (h3 + 0x9e3779b97f4a7c15ULL + (h12 << 6) + (h12 >> 2));
    }
  };

  static std::unique_ptr<ExplicitBitVect> morgan_fp(const ROMol &mol,
                                                    unsigned int nBits = 2048,
                                                    unsigned int radius = 2) {
    return std::unique_ptr<ExplicitBitVect>(
        MorganFingerprints::getFingerprintAsBitVect(mol, radius, nBits));
  }
  static bool is_subset_fp(const ExplicitBitVect &pattern,
                           const ExplicitBitVect &mol) {
    IntVect onbits;
    pattern.getOnBits(onbits);
    for (auto b : onbits)
      if (!mol.getBit(b))
        return false;
    return true;
  }

  std::vector<ChemicalReaction *> reactions_;
  std::vector<MolPtr> molecule_pool_;
  std::unordered_set<std::string> seen_smiles_;
  std::unordered_map<int, LineageNode> lineage_map_; // product_index -> lineage
  std::vector<std::tuple<ComboKey, std::vector<MolPtr>, int>> pending_results_;

  std::unordered_set<ReactionComboKey, ReactionComboKeyHash>
      seen_reaction_combos_;
  std::vector<ReactantTemplate> reactant_templates_;
  std::vector<std::vector<VectMatchVectType>>
      template_match_matrix_; // [template][mol]
  std::unordered_map<int, std::pair<int, int>>
      reaction_to_template_indices_; // rxn -> (tA,tB)
  std::vector<AtomLevelReaction> pending_atom_level_reactions_;

  std::vector<MolPtr> prepare_reactants(const ComboKey &combo_key) {
    std::vector<MolPtr> combo_mols;
    for (int idx : {combo_key.a_index, combo_key.b_index}) {
      if (idx >= 0 && idx < static_cast<int>(molecule_pool_.size()) &&
          molecule_pool_[idx])
        combo_mols.push_back(molecule_pool_[idx]);
    }
    return combo_mols;
  }

  std::tuple<ComboKey, std::vector<MolPtr>, int>
  process_one_atom_level_reaction(const AtomLevelReaction &atom_rxn) {
    const ComboKey combo_key = atom_rxn.combo_key;
    std::vector<MolPtr> combo_mols = prepare_reactants(combo_key);
    ChemicalReaction *rxn = reactions_[atom_rxn.reaction_template_index];

    std::vector<MatchVectType> reactantsMatch =
        (combo_mols.size() == 1)
            ? std::vector<MatchVectType>{atom_rxn.atom_indices[0]}
            : std::vector<MatchVectType>{atom_rxn.atom_indices[0],
                                         atom_rxn.atom_indices[1]};

    MOL_SPTR_VECT prods = ReactionRunnerUtils::generateOneProductSet(
        *rxn, combo_mols, reactantsMatch);
    return std::make_tuple(combo_key, prods, atom_rxn.reaction_template_index);
  }

  void generate_new_combinations(bool verbose = false) {
    pending_atom_level_reactions_.clear();
    const int R = static_cast<int>(reactions_.size());
    for (int rxn_idx = 0; rxn_idx < R; ++rxn_idx) {
      auto [tA, tB] = reaction_to_template_indices_[rxn_idx];
      const auto &A = template_match_matrix_[tA];
      if (tB < 0) {
        for (int i = 0, M = static_cast<int>(A.size()); i < M; ++i) {
          if (A[i].empty())
            continue;
          ReactionComboKey rck{rxn_idx, ComboKey(i, -1)};
          if (!seen_reaction_combos_.insert(rck).second)
            continue;
          for (const auto &ma : A[i])
            pending_atom_level_reactions_.emplace_back(
                ComboKey(i, -1), std::vector<MatchVectType>{ma}, rxn_idx);
        }
      } else {
        const auto &B = template_match_matrix_[tB];
        for (int i = 0, Mi = static_cast<int>(A.size()); i < Mi; ++i) {
          if (A[i].empty())
            continue;
          for (int j = 0, Mj = static_cast<int>(B.size()); j < Mj; ++j) {
            if (B[j].empty())
              continue;
            ReactionComboKey rck{rxn_idx, ComboKey(i, j)};
            if (!seen_reaction_combos_.insert(rck).second)
              continue;
            for (const auto &ma : A[i])
              for (const auto &mb : B[j])
                pending_atom_level_reactions_.emplace_back(
                    ComboKey(i, j), std::vector<MatchVectType>{ma, mb},
                    rxn_idx);
          }
        }
      }
      // if (verbose) {
      //   if (rxn_idx % (R / 10) == 0) {
      //     std::cout << ">>>> completed combinatorials for rxn " << rxn_idx
      //               << " of " << R << " (" << current_time_ << "s)\n";
      //   }
      // }
    }
  }

  int perform_template_matching(bool parallelize = false, int num_threads = 0,
                                bool verbose = false) {
    if (!timer_started_) {
      start_time_ = std::chrono::steady_clock::now();
      timer_started_ = true;
    }
    // if (verbose) {
    //   std::cout << "....reactant_templates_.size() = "
    //             << reactant_templates_.size() << std::endl;

    //   current_time_ =
    //   std::chrono::duration_cast<std::chrono::duration<double>>(
    //                       std::chrono::steady_clock::now() - start_time_)
    //                       .count();
    //   std::cout << "curr_time = " << current_time_ << std::endl;
    // }
    const size_t T = reactant_templates_.size();
    const size_t M = molecule_pool_.size();
    template_match_matrix_.assign(T, std::vector<VectMatchVectType>(M));

    RDKit::SubstructMatchParameters params;
    params.maxMatches = 1000000;
    params.uniquify = true;

    if (!parallelize || T == 0 || M == 0) {
      // single-threaded path
      if (verbose) {
        std::cout << "....single-threaded path"
                  << " T = " << T << " M = " << M << std::endl;
      }
      for (size_t t = 0; t < T; ++t) {
        // Check for interrupt signal
        if (g_interrupt.load()) {
          std::cout << "\n[C++] Template matching interrupted\n";
          return -1;
        }

        // if (verbose) {
        //   std::cout << ".....template " << t << " of " << T << std::endl;
        // }

        for (size_t m = 0; m < M; ++m) {
          current_time_ =
              std::chrono::duration_cast<std::chrono::duration<double>>(
                  std::chrono::steady_clock::now() - start_time_)
                  .count();
          // if (verbose) {
          // if (m % (M / 10) == 0) {
          //   std::cout << ".....molecule " << m << " of " << M << " ("
          //             << current_time_ << "s)\n";
          // }
          // }

          if (current_time_ > max_runtime_) {
            std::cout << "\n[C++] Template matching interrupted [max runtime] ("
                      << current_time_ << "s)\n";
            return -1;
          }

          template_match_matrix_[t][m] = RDKit::SubstructMatch(
              *molecule_pool_[m], *reactant_templates_[t].template_mol, params);
        }

        if (verbose) {
          if (t % 10 == 0) {
            std::cout << ">>> Completed template matching for template " << t;
            std::cout << " of " << T << " (" << current_time_ << "s)\n";
          }
        }
      }
      return 0;
    }
    // multi-threaded path
    const size_t total = T * M;
    std::atomic<size_t> next_idx{0};

    size_t nthreads = (num_threads > 0) ? static_cast<size_t>(num_threads)
                                        : std::thread::hardware_concurrency();
    if (nthreads == 0)
      nthreads = 1; // fallback

    std::vector<std::thread> workers;
    workers.reserve(nthreads);

    for (size_t tid = 0; tid < nthreads; ++tid) {
      workers.emplace_back([&, T, M, total]() {
        while (true) {
          const size_t idx = next_idx.fetch_add(1, std::memory_order_relaxed);
          if (idx >= total)
            break;
          const size_t t = idx / M;
          const size_t m = idx % M;

          // safe: read-only molecules; each thread writes a unique cell
          template_match_matrix_[t][m] = RDKit::SubstructMatch(
              *molecule_pool_[m], *reactant_templates_[t].template_mol, params);
        }
      });
    }

    for (auto &th : workers) {
      th.join();
    }
    return 0;
  }
}; // end of ReactionEngine class

BOOST_PYTHON_MODULE(_rxn_backend_min) {
  // Make C/C++ prints show immediately when called from Python
  std::cout.setf(std::ios::unitbuf);
  std::cerr.setf(std::ios::unitbuf);
  setvbuf(stdout, nullptr, _IONBF, 0);
  setvbuf(stderr, nullptr, _IONBF, 0);

  bp::class_<ReactionEngine, boost::noncopyable>("ReactionEngine")
      .def("load_reactions", &ReactionEngine::load_reactions)
      .def("load_molecules_into_pool",
           &ReactionEngine::load_molecules_into_pool)
      .def("run_pending", &ReactionEngine::run_pending,
           (bp::arg("verbose") = false))
      .def("run_one_atom_level_reaction",
           &ReactionEngine::run_one_atom_level_reaction,
           (bp::arg("reaction_index"), bp::arg("verbose") = false))
      .def("add_products_to_pool", &ReactionEngine::add_products_to_pool,
           (bp::arg("parallelize") = false, bp::arg("num_threads") = 0,
            bp::arg("verbose") = false))
      .def("get_molecule_lineage", &ReactionEngine::get_molecule_lineage)
      .def("get_full_lineage", &ReactionEngine::get_full_lineage)
      .def("get_pool_size", &ReactionEngine::get_pool_size)
      .def("get_all_molecules_as_smiles",
           &ReactionEngine::get_all_molecules_as_smiles)
      .def("get_num_reactant_templates",
           &ReactionEngine::get_num_reactant_templates)
      .def("get_reactant_template_info",
           &ReactionEngine::get_reactant_template_info)
      .def("get_template_matches", &ReactionEngine::get_template_matches)
      .def("get_all_template_matches",
           &ReactionEngine::get_all_template_matches)
      .def("get_pending_combinations_info",
           &ReactionEngine::get_pending_combinations_info)
      .def("get_num_pending_atom_level_reactions",
           &ReactionEngine::get_num_pending_atom_level_reactions)
      .def("get_molecule_smiles_from_pool",
           &ReactionEngine::get_molecule_smiles_from_pool)
      .def("get_atom_level_reactions",
           &ReactionEngine::get_atom_level_reactions)
      .def("get_reaction_forming_product",
           &ReactionEngine::get_reaction_forming_product)
      .def("get_reaction_smarts_from_index",
           &ReactionEngine::get_reaction_smarts_from_index);
}
