#include <RDBoost/python.h>
#include <RDBoost/Wrap.h>

#include <boost/python/stl_iterator.hpp>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

#include <GraphMol/AtomTyper/atom_typer.hpp>
#include <GraphMol/AtomTyper/smarts_analyzer.hpp>
#include <GraphMol/AtomTyper/expression_builder.hpp>

namespace python = boost::python;

namespace {
python::list vectorIntToList(const std::vector<int> &vals) {
  python::list res;
  for (const auto v : vals) {
    res.append(v);
  }
  return res;
}

python::dict mapIntIntToDict(const std::map<int, int> &vals) {
  python::dict res;
  for (const auto &kv : vals) {
    res[kv.first] = kv.second;
  }
  return res;
}

std::vector<atom_typer::AtomType> extractAtomTypes(
    const std::vector<atom_typer::PatternItem> &items) {
  std::vector<atom_typer::AtomType> atoms;
  atoms.reserve(items.size());
  for (const auto &item : items) {
    if (item.kind == atom_typer::PatternItemKind::Atom) {
      atoms.push_back(item.atom);
    }
  }
  return atoms;
}

python::list typeAtomsFromSmiles(atom_typer::AtomTyper &self,
                                 const std::string &smiles) {
  auto atomTypes = extractAtomTypes(self.type_atoms_from_smiles(smiles));
  python::list res;
  for (const auto &at : atomTypes) {
    res.append(at);
  }
  return res;
}

std::string typeAtomsFromSmarts(atom_typer::AtomTyper &self,
                                const std::string &smarts) {
  return self.type_atoms_from_smarts(smarts);
}

python::tuple typeAtomsFromSmartsMapNewAtoms(
    atom_typer::AtomTyper &self, const std::string &smarts, bool map_new_atoms,
    int max_amap, bool verbose, bool include_x_in_reserialization) {
  auto res = self.type_atoms_from_smarts(smarts, map_new_atoms, max_amap,
                                         verbose, include_x_in_reserialization);
  return python::make_tuple(res, max_amap);
}

python::tuple typeAtomsFromSmartsWithRanges(
    atom_typer::AtomTyper &self, const std::string &smarts, int h_min,
    int h_max, int charge_min, int charge_max, bool map_new_atoms,
    int max_amap, bool include_x_in_reserialization) {
  auto res =
      self.type_atoms_from_smarts(smarts, h_min, h_max, charge_min,
                                  charge_max, map_new_atoms, max_amap,
                                  include_x_in_reserialization);
  return python::make_tuple(res, max_amap);
}

python::tuple typeAtomsFromSmartsWithRangesDebug(
    atom_typer::AtomTyper &self, const std::string &smarts, int h_min,
    int h_max, int charge_min, int charge_max, bool verbose,
    atom_typer::DebugLevel debug_level, bool map_new_atoms, int max_amap,
    bool include_x_in_reserialization) {
  auto res = self.type_atoms_from_smarts(
      smarts, h_min, h_max, charge_min, charge_max, verbose, debug_level,
      map_new_atoms, max_amap, include_x_in_reserialization);
  return python::make_tuple(res, max_amap);
}

std::string typeAtomsFromSmartsWithDebugLevel(atom_typer::AtomTyper &self,
                                              const std::string &smarts,
                                              atom_typer::DebugLevel level) {
  return self.type_atoms_from_smarts(smarts, level);
}

std::string getAtomTypesString(atom_typer::AtomTyper &self,
                               const python::object &pyAtomTypes) {
  std::vector<atom_typer::AtomType> atomTypes;
  python::stl_input_iterator<atom_typer::AtomType> begin(pyAtomTypes), end;
  for (auto it = begin; it != end; ++it) {
    atomTypes.push_back(*it);
  }
  return self.get_atom_types_string(atomTypes);
}

python::list enumerateVariants(atom_typer::SmartsAnalyzer &self,
                               const std::string &smarts, int maxVariants) {
  auto variants = self.enumerate_variants(smarts, maxVariants);
  python::list res;
  for (const auto &v : variants) {
    res.append(v);
  }
  return res;
}

python::list ringMembershipList(const atom_typer::AtomType &self) {
  return vectorIntToList(self.ring_membership_list);
}

python::list neighborsList(const atom_typer::AtomType &self) {
  return vectorIntToList(self.neighbors);
}

python::dict bondTypesDict(const atom_typer::AtomType &self) {
  return mapIntIntToDict(self.bond_types);
}

python::list vectorStringToList(const std::vector<std::string> &vals) {
  python::list res;
  for (const auto &v : vals) {
    res.append(v);
  }
  return res;
}

python::list atomTypesToList(const std::vector<atom_typer::AtomType> &vals) {
  python::list res;
  for (const auto &v : vals) {
    res.append(v);
  }
  return res;
}

python::object typeSmiles(const python::object &smiles_or_list,
                          bool use_canonical = true,
                          bool reserialize = true) {
  std::vector<std::string> inputs;

  if (PyUnicode_Check(smiles_or_list.ptr())) {
    inputs.push_back(python::extract<std::string>(smiles_or_list));
  } else {
    python::stl_input_iterator<std::string> begin(smiles_or_list), end;
    for (auto it = begin; it != end; ++it) {
      inputs.push_back(*it);
    }
  }

  if (inputs.empty()) {
    PyErr_SetString(
        PyExc_ValueError,
        "type_smiles requires a SMILES string or non-empty iterable of SMILES strings");
    python::throw_error_already_set();
  }

  atom_typer::AtomTyper typer;
  typer.set_use_canonical(use_canonical);

  if (PyUnicode_Check(smiles_or_list.ptr())) {
    auto pattern_items = typer.type_atoms_from_smiles(inputs.front());
    if (reserialize) {
      return python::object(typer.get_smarts_from_pattern_types(pattern_items));
    }
    auto atom_types = extractAtomTypes(pattern_items);
    return atomTypesToList(atom_types);
  }

  python::list out;
  for (const auto &smi : inputs) {
    auto pattern_items = typer.type_atoms_from_smiles(smi);
    if (reserialize) {
      out.append(typer.get_smarts_from_pattern_types(pattern_items));
    } else {
      auto atom_types = extractAtomTypes(pattern_items);
      out.append(atomTypesToList(atom_types));
    }
  }
  return out;
}

python::object typeSmarts(const python::object &smarts_or_list,
                          bool verbose = false,
                          bool include_x_in_reserialization = false,
                          bool ignoreValence = false,
                          bool catchErrors = true,
                          bool enumerate_bond_order = true,
                          bool log_enabled = false,
                          unsigned int log_flags = atom_typer::SmartsAnalyzer::LogAll,
                          unsigned int extracted_primitives_mask = atom_typer::SmartsAnalyzer::ExtractPrimitiveNone,
                          const python::object &factoring_priority = python::object(),
                            bool remove_aa_wildcard = true,
                          atom_typer::SmartsAnalyzer::SymbolForm symbol_form =
                              atom_typer::SmartsAnalyzer::SymbolForm::Unchanged,
                          bool fold_singleton_or = false,
                            bool explicit_charge_values = false,
                            const python::object &or_primitive_rewrite_atomic_nums =
                              python::object()) {
  std::vector<std::string> inputs;

  if (PyUnicode_Check(smarts_or_list.ptr())) {
    inputs.push_back(python::extract<std::string>(smarts_or_list));
  } else {
    python::stl_input_iterator<std::string> begin(smarts_or_list), end;
    for (auto it = begin; it != end; ++it) {
      inputs.push_back(*it);
    }
  }

  if (inputs.empty()) {
    PyErr_SetString(PyExc_ValueError,
                    "type_smarts requires a SMARTS string or non-empty iterable of SMARTS strings");
    python::throw_error_already_set();
  }

  std::vector<std::string> factoring_priority_vec;
  if (factoring_priority.ptr() != Py_None) {
    python::stl_input_iterator<std::string> begin(factoring_priority), end;
    for (auto it = begin; it != end; ++it) {
      factoring_priority_vec.push_back(*it);
    }
  }

  std::vector<int> rewrite_atomic_nums_vec;
  if (or_primitive_rewrite_atomic_nums.ptr() != Py_None) {
    python::stl_input_iterator<int> begin(or_primitive_rewrite_atomic_nums),
        end;
    for (auto it = begin; it != end; ++it) {
      rewrite_atomic_nums_vec.push_back(*it);
    }
  }

  atom_typer::SmartsAnalyzer analyzer;
  atom_typer::SmartsAnalyzer::StandardSmartsWorkflowOptions workflow_options(
      include_x_in_reserialization, enumerate_bond_order,
      extracted_primitives_mask, std::move(factoring_priority_vec));
  workflow_options.remove_aa_wildcard = remove_aa_wildcard;
  workflow_options.symbol_form = symbol_form;
  workflow_options.fold_singleton_or = fold_singleton_or;
  workflow_options.explicit_charge_values = explicit_charge_values;
    workflow_options.rewrite_or_primitives_to_negated =
      !rewrite_atomic_nums_vec.empty();
    workflow_options.or_primitive_rewrite_atomic_nums =
      std::move(rewrite_atomic_nums_vec);
  atom_typer::SmartsAnalyzer::StandardSmartsLogOptions log_options(
      log_enabled, log_flags);
  auto standardized = analyzer.standard_smarts(
      inputs, verbose, ignoreValence, catchErrors,
      workflow_options, log_options);

  if (PyUnicode_Check(smarts_or_list.ptr())) {
    if (standardized.empty()) {
      return python::object(std::string(""));
    }
    return python::object(standardized.front());
  }
  return vectorStringToList(standardized);
}

}  // namespace

BOOST_PYTHON_MODULE(rdAtomTyper) {
  python::scope().attr("__doc__") =
      "AtomTyper: type atoms in SMILES/SMARTS strings based on atomic primitives "
      "and local connectivity.";

  python::enum_<atom_typer::Level>("Level")
      .value("MINIMAL", atom_typer::Level::MINIMAL)
      .value("STANDARD", atom_typer::Level::STANDARD)
      .value("DETAILED", atom_typer::Level::DETAILED)
      .value("COMPLETE", atom_typer::Level::COMPLETE);

  python::enum_<atom_typer::SmartsAnalyzer::SymbolForm>("SymbolForm")
      .value("Unchanged", atom_typer::SmartsAnalyzer::SymbolForm::Unchanged)
      .value("Expanded",  atom_typer::SmartsAnalyzer::SymbolForm::Expanded)
      .value("Condensed", atom_typer::SmartsAnalyzer::SymbolForm::Condensed);

    python::enum_<atom_typer::DebugLevel>("DebugLevel")
      .value("Off", atom_typer::DebugLevel::Off)
      .value("Basic", atom_typer::DebugLevel::Basic)
      .value("Verbose", atom_typer::DebugLevel::Verbose)
      .value("Trace", atom_typer::DebugLevel::Trace);

  python::class_<atom_typer::AtomType>("AtomType")
      .def(python::init<>())
      .def_readwrite("atom_idx", &atom_typer::AtomType::atom_idx)
      .def_readwrite("atomic_number", &atom_typer::AtomType::atomic_number)
      .def_readwrite("formal_charge", &atom_typer::AtomType::formal_charge)
      .def_readwrite("num_hydrogens", &atom_typer::AtomType::num_hydrogens)
      .def_readwrite("min_bonds", &atom_typer::AtomType::min_bonds)
      .def_readwrite("max_valence", &atom_typer::AtomType::max_valence)
      .def_readwrite("is_aromatic", &atom_typer::AtomType::is_aromatic)
      .def_readwrite("is_aliphatic", &atom_typer::AtomType::is_aliphatic)
      .def_readwrite("is_in_ring", &atom_typer::AtomType::is_in_ring)
      .def_readwrite("ring_size", &atom_typer::AtomType::ring_size)
      .def_readwrite("hybridization", &atom_typer::AtomType::hybridization)
      .def_readwrite("smarts_pattern", &atom_typer::AtomType::smarts_pattern)
      .def_readwrite("num_ring_bonds", &atom_typer::AtomType::num_ring_bonds)
      .def_readwrite("num_aliphatic_rings",
                     &atom_typer::AtomType::num_aliphatic_rings)
      .def_readwrite("num_aromatic_rings",
                     &atom_typer::AtomType::num_aromatic_rings)
      .def_readwrite("chirality", &atom_typer::AtomType::chirality)
      .def_readwrite("ring_connectivity",
                     &atom_typer::AtomType::ring_connectivity)
      .def_readwrite("atom_type_enumeration",
             &atom_typer::AtomType::atom_type_enumeration)
      .def_readwrite("num_single_bonds",
             &atom_typer::AtomType::num_single_bonds)
      .def_readwrite("num_double_bonds",
             &atom_typer::AtomType::num_double_bonds)
      .def_readwrite("num_triple_bonds",
             &atom_typer::AtomType::num_triple_bonds)
      .def_readwrite("num_aromatic_bonds",
             &atom_typer::AtomType::num_aromatic_bonds)
      .def_readwrite("remaining_valence",
             &atom_typer::AtomType::remaining_valence)
      .def_readwrite("source_atom_smarts",
             &atom_typer::AtomType::source_atom_smarts)
      // STL containers: expose read-only views for consistency and to avoid
      // requiring Boost.Python container converters.
      .add_property("neighbors", &neighborsList)
      .add_property("ring_membership_list", &ringMembershipList)
      .add_property("bond_types", &bondTypesDict);

  python::class_<atom_typer::AtomTyper,
                 boost::shared_ptr<atom_typer::AtomTyper>, boost::noncopyable>(
      "AtomTyper")
      .def(python::init<>())
      .def("type_atoms_from_smiles", &typeAtomsFromSmiles,
           (python::arg("smiles")))
      .def("type_atoms_from_smarts", &typeAtomsFromSmarts,
           (python::arg("smarts")))
       .def("type_atoms_from_smarts", &typeAtomsFromSmartsMapNewAtoms,
         (python::arg("smarts"), python::arg("map_new_atoms"),
          python::arg("max_amap"), python::arg("verbose") = false,
          python::arg("include_x_in_reserialization") = false))
       .def("type_atoms_from_smarts", &typeAtomsFromSmartsWithRanges,
         (python::arg("smarts"), python::arg("h_min"),
          python::arg("h_max"), python::arg("charge_min"),
          python::arg("charge_max"), python::arg("map_new_atoms"),
          python::arg("max_amap"),
          python::arg("include_x_in_reserialization") = false))
       .def("type_atoms_from_smarts", &typeAtomsFromSmartsWithRangesDebug,
         (python::arg("smarts"), python::arg("h_min"),
          python::arg("h_max"), python::arg("charge_min"),
          python::arg("charge_max"), python::arg("verbose"),
          python::arg("debug_level"), python::arg("map_new_atoms"),
          python::arg("max_amap"),
          python::arg("include_x_in_reserialization") = false))
       .def("type_atoms_from_smarts", &typeAtomsFromSmartsWithDebugLevel,
         (python::arg("smarts"), python::arg("debug_level")))
      .def("get_atom_types_string", &getAtomTypesString,
           (python::arg("atom_types")))
      .def("set_use_canonical", &atom_typer::AtomTyper::set_use_canonical,
         (python::arg("use_canonical")))
       .def("set_debug_level", &atom_typer::AtomTyper::set_debug_level,
         (python::arg("debug_level")));

  python::class_<atom_typer::SmartsAnalyzer,
                 boost::shared_ptr<atom_typer::SmartsAnalyzer>,
                 boost::noncopyable>("SmartsAnalyzer")
      .def(python::init<>())
      .def("calculate_dof", &atom_typer::SmartsAnalyzer::calculate_dof,
           (python::arg("smarts")))
      .def("enumerate_variants", &enumerateVariants,
           (python::arg("smarts"), python::arg("max")));

  python::def("smiles_to_smarts", &atom_typer::smiles_to_smarts,
              (python::arg("smiles"), python::arg("level")));

  python::scope().attr("LOG_ALL") = atom_typer::SmartsAnalyzer::LogAll;

  python::scope().attr("EXTRACT_PRIMITIVE_NONE") =
      atom_typer::SmartsAnalyzer::ExtractPrimitiveNone;
  python::scope().attr("EXTRACT_PRIMITIVE_H") =
      atom_typer::SmartsAnalyzer::ExtractPrimitiveH;
  python::scope().attr("EXTRACT_PRIMITIVE_D") =
      atom_typer::SmartsAnalyzer::ExtractPrimitiveD;
  python::scope().attr("EXTRACT_PRIMITIVE_X") =
      atom_typer::SmartsAnalyzer::ExtractPrimitiveX;
  python::scope().attr("EXTRACT_PRIMITIVE_CHARGE") =
      atom_typer::SmartsAnalyzer::ExtractPrimitiveCharge;
  python::scope().attr("EXTRACT_PRIMITIVE_V") =
      atom_typer::SmartsAnalyzer::ExtractPrimitiveV;
  python::scope().attr("EXTRACT_PRIMITIVE_RING_COUNT") =
      atom_typer::SmartsAnalyzer::ExtractPrimitiveRingCount;
  python::scope().attr("EXTRACT_PRIMITIVE_RING_SIZE") =
      atom_typer::SmartsAnalyzer::ExtractPrimitiveRingSize;
  python::scope().attr("EXTRACT_PRIMITIVE_RING_BOND_COUNT") =
      atom_typer::SmartsAnalyzer::ExtractPrimitiveRingBondCount;

  python::def(
      "type_smarts", &typeSmarts,
      (python::arg("smarts_or_list"), python::arg("verbose") = false,
       python::arg("include_x_in_reserialization") = false,
       python::arg("ignoreValence") = false,
       python::arg("catchErrors") = true,
       python::arg("enumerate_bond_order") = false,
       python::arg("log_enabled") = false,
       python::arg("log_flags") = atom_typer::SmartsAnalyzer::LogAll,
       python::arg("extracted_primitives_mask") =
           atom_typer::SmartsAnalyzer::ExtractPrimitiveNone,
       python::arg("factoring_priority") = python::object(),
         python::arg("remove_aa_wildcard") = true,
       python::arg("symbol_form") =
           atom_typer::SmartsAnalyzer::SymbolForm::Unchanged,
       python::arg("fold_singleton_or") = false,
         python::arg("explicit_charge_values") = false,
         python::arg("or_primitive_rewrite_atomic_nums") = python::object()),
      "Standardize SMARTS via SmartsAnalyzer::standard_smarts.\n"
      "Accepts either a single SMARTS string or an iterable of SMARTS strings.\n\n"
      "Parameters:\n"
      "  smarts_or_list: str or iterable of str\n"
      "  verbose: print debug info (default False)\n"
      "  include_x_in_reserialization: include [#0] in output (default False)\n"
      "  ignoreValence: skip valence checks (default False)\n"
      "  catchErrors: catch and skip errors (default True)\n"
      "  enumerate_bond_order: enumerate bond-order variants (default False)\n"
      "  log_enabled: enable logging (default False)\n"
      "  log_flags: bitmask controlling log categories (default LOG_ALL)\n"
      "  extracted_primitives_mask: bitmask of primitives to extract\n"
      "      (default EXTRACT_PRIMITIVE_NONE, combine with |)\n"
      "  factoring_priority: list of SMARTS primitives for factoring order\n"
      "      (default None)\n"
      "  remove_aa_wildcard: remove OR(A,a) aromaticity tautologies (default True)\n"
      "  symbol_form: SymbolForm.Unchanged / Expanded (C->[#6&A]) /\n"
      "      Condensed ([#6&A]->C) (default Unchanged)\n"
      "  fold_singleton_or: collapse single-arm OR nodes (default False)\n"
      "  explicit_charge_values: expand +/-/++/-- to +1/-1/+2/-2 (default False)\n"
      "  or_primitive_rewrite_atomic_nums: iterable of atomic numbers where\n"
      "      OR->negated primitive rewrites are enabled\n"
      "      (e.g. [6] for carbon only; default None/empty = disabled)\n");

    python::def(
      "type_smiles", &typeSmiles,
      (python::arg("smiles_or_list"), python::arg("use_canonical") = true,
       python::arg("reserialize") = true),
      "Type SMILES via AtomTyper::type_atoms_from_smiles. Accepts a single SMILES string or an iterable of SMILES strings. If reserialize=True, returns a typed SMARTS-like string reconstructed from PatternItem atom/bond tokens instead of AtomType objects.");
}
