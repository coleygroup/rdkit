#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "atom_typer.hpp"
#include "smarts_analyzer.hpp"
#include "expression_builder.hpp"

namespace py = pybind11;

PYBIND11_MODULE(atom_typer, m) {
    m.doc() = "Atom Typer - Type atoms in SMILES/SMARTS strings based on atomic primitives and local connectivity";
    
    // Bind AtomType struct
    py::class_<atom_typer::AtomType>(m, "AtomType")
        .def(py::init<>())
        .def_readwrite("atom_idx", &atom_typer::AtomType::atom_idx)
        .def_readwrite("atomic_number", &atom_typer::AtomType::atomic_number)
        .def_readwrite("formal_charge", &atom_typer::AtomType::formal_charge)
        .def_readwrite("num_hydrogens", &atom_typer::AtomType::num_hydrogens)
        .def_readwrite("degree", &atom_typer::AtomType::degree)
        .def_readwrite("valence", &atom_typer::AtomType::valence)
        .def_readwrite("is_aromatic", &atom_typer::AtomType::is_aromatic)
        .def_readwrite("is_in_ring", &atom_typer::AtomType::is_in_ring)
        .def_readwrite("ring_size", &atom_typer::AtomType::ring_size)
        .def_readwrite("hybridization", &atom_typer::AtomType::hybridization)
        .def_readwrite("neighbors", &atom_typer::AtomType::neighbors)
        .def_readwrite("smarts_pattern", &atom_typer::AtomType::smarts_pattern)
        .def("__repr__", [](const atom_typer::AtomType& at) {
            return "<AtomType idx=" + std::to_string(at.atom_idx) + 
                   " element=" + std::to_string(at.atomic_number) + 
                   " smarts=" + at.smarts_pattern + ">";
        });
    
    // Bind AtomTyper class
    py::class_<atom_typer::AtomTyper>(m, "AtomTyper")
        .def(py::init<>())
        .def("type_atoms_from_smiles", &atom_typer::AtomTyper::type_atoms_from_smiles,
             py::arg("smiles"),
             "Type all atoms in a SMILES string")
        .def("type_atoms_from_smarts", &atom_typer::AtomTyper::type_atoms_from_smarts,
             py::arg("smarts"),
             "Type all atoms in a SMARTS string")
        .def("get_atom_types_string", &atom_typer::AtomTyper::get_atom_types_string,
             py::arg("atom_types"),
             "Get atom types as a formatted string")
        .def("set_use_canonical", &atom_typer::AtomTyper::set_use_canonical,
             py::arg("use_canonical"),
             "Set whether to use canonical SMILES");
    
    // Bind SmartsAnalyzer class
    py::class_<atom_typer::SmartsAnalyzer>(m, "SmartsAnalyzer")
        .def(py::init<>())
        .def("calculate_dof", &atom_typer::SmartsAnalyzer::calculate_dof,
             py::arg("smarts"),
             "Calculate the degrees of freedom (number of possible variants) for a SMARTS pattern")
        .def("enumerate_variants", &atom_typer::SmartsAnalyzer::enumerate_variants,
             py::arg("smarts"), py::arg("max"),
             "Enumerate all possible SMARTS variants up to a maximum number");
    
    // Bind Level enum
    py::enum_<atom_typer::Level>(m, "Level")
        .value("MINIMAL", atom_typer::Level::MINIMAL, "Basic element types only")
        .value("STANDARD", atom_typer::Level::STANDARD, "Elements + bond types + aromaticity")
        .value("DETAILED", atom_typer::Level::DETAILED, "Add hydrogen counts, formal charges")
        .value("COMPLETE", atom_typer::Level::COMPLETE, "All features: valence, connectivity, ring membership")
        .export_values();
    
    // Bind smiles_to_smarts function
    m.def("smiles_to_smarts", &atom_typer::smiles_to_smarts,
          py::arg("smiles"), py::arg("level"),
          "Convert a SMILES string to a SMARTS pattern with specified detail level");
}
