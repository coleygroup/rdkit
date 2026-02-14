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

python::list typeAtomsFromSmiles(atom_typer::AtomTyper &self,
                                 const std::string &smiles) {
  auto atomTypes = self.type_atoms_from_smiles(smiles);
  python::list res;
  for (const auto &at : atomTypes) {
    res.append(at);
  }
  return res;
}

python::list typeAtomsFromSmarts(atom_typer::AtomTyper &self,
                                 const std::string &smarts) {
  auto atomTypes = self.type_atoms_from_smarts(smarts);
  python::list res;
  for (const auto &at : atomTypes) {
    res.append(at);
  }
  return res;
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

  python::class_<atom_typer::AtomType>("AtomType")
      .def(python::init<>())
      .def_readwrite("atom_idx", &atom_typer::AtomType::atom_idx)
      .def_readwrite("atomic_number", &atom_typer::AtomType::atomic_number)
      .def_readwrite("formal_charge", &atom_typer::AtomType::formal_charge)
      .def_readwrite("num_hydrogens", &atom_typer::AtomType::num_hydrogens)
      .def_readwrite("min_bonds", &atom_typer::AtomType::min_bonds)
      .def_readwrite("valence", &atom_typer::AtomType::valence)
      .def_readwrite("is_aromatic", &atom_typer::AtomType::is_aromatic)
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
       .def("enumerate_dof_smarts", &atom_typer::AtomTyper::enumerate_dof_smarts,
         (python::arg("smarts")))
      .def("get_atom_types_string", &getAtomTypesString,
           (python::arg("atom_types")))
      .def("set_use_canonical", &atom_typer::AtomTyper::set_use_canonical,
           (python::arg("use_canonical")));

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
}
