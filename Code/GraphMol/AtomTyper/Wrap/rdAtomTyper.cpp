#include <RDBoost/python.h>
#include <RDBoost/Wrap.h>

#include <boost/python/stl_iterator.hpp>

#include <GraphMol/AtomTyper/atom_typer.hpp>
#include <GraphMol/AtomTyper/expression_builder.hpp>
#include <GraphMol/AtomTyper/reaction_extractor.hpp>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/GraphMol.h>
// #include <GraphMol/Atom.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/QueryOps.h>

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

python::list vectorVectorStringToList(
    const std::vector<std::vector<std::string>> &vals) {
  python::list res;
  for (const auto &inner : vals) {
    res.append(vectorStringToList(inner));
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

// ---------------------------------------------------------------------------
// Query-tree serialisation: walk RDKit's atom/bond query trees and return
// them as plain Python dicts so the visualiser frontend can render them.
// ---------------------------------------------------------------------------

using ATOM_QUERY = Queries::Query<int, const RDKit::Atom *, true>;
using BOND_QUERY = Queries::Query<int, const RDKit::Bond *, true>;

// Forward declaration
python::dict walkAtomQuery(const ATOM_QUERY *q);

python::dict walkAtomQuery(const ATOM_QUERY *q) {
  python::dict d;
  if (!q) {
    d["description"] = "AtomNull";
    d["op"] = "leaf";
    return d;
  }
  const auto &desc = q->getDescription();
  d["description"] = desc;
  d["negated"] = q->getNegation();

  // Classify the node
  if (desc == "AtomAnd" || desc == "And") {
    d["op"] = "and";
  } else if (desc == "AtomOr" || desc == "Or") {
    d["op"] = "or";
  } else if (desc == "AtomNot" || desc == "Not") {
    d["op"] = "not";
  } else if (desc == "RecursiveSmarts" || desc == "RecursiveStructure") {
    d["op"] = "recursive";
  } else {
    d["op"] = "leaf";
    // Try to extract the integer value via ATOM_EQUALS_QUERY
    try {
      const auto *eq =
          dynamic_cast<const Queries::EqualityQuery<int, const RDKit::Atom *, true> *>(q);
      if (eq) {
        if (desc == "AtomHybridization") {
              switch (eq->getVal()) {
                case RDKit::QueryAtom::S:
                  d["value"] = 0;
                  break;
                case RDKit::QueryAtom::SP:
                  d["value"] = 1;
                  break;
                case RDKit::QueryAtom::SP2:
                  d["value"] = 2;
                  break;
                case RDKit::QueryAtom::SP3:
                  d["value"] = 3;
                  break;
                case RDKit::QueryAtom::SP3D:
                  d["value"] = 4;
                  break;
                case RDKit::QueryAtom::SP3D2:
                  d["value"] = 5;
                  break;
              }
        } else{
          d["value"] = eq->getVal();
        }
      }
    } catch (...) {
      // not an equality query – leave value unset
    }
  }

  // Recurse into children
  python::list children;
  for (auto it = q->beginChildren(); it != q->endChildren(); ++it) {
    children.append(walkAtomQuery(it->get()));
  }
  if (python::len(children) > 0) {
    d["children"] = children;
  }
  return d;
}

python::dict walkBondQuery(const BOND_QUERY *q) {
  python::dict d;
  if (!q) {
    d["description"] = "BondNull";
    d["op"] = "leaf";
    return d;
  }
  const auto &desc = q->getDescription();
  d["description"] = desc;
  d["negated"] = q->getNegation();

  if (desc == "BondAnd" || desc == "And") {
    d["op"] = "and";
  } else if (desc == "BondOr" || desc == "Or") {
    d["op"] = "or";
  } else if (desc == "BondNot" || desc == "Not") {
    d["op"] = "not";
  } else {
    d["op"] = "leaf";
    try {
      const auto *eq =
          dynamic_cast<const Queries::EqualityQuery<int, const RDKit::Bond *, true> *>(q);
      if (eq) {
        d["value"] = eq->getVal();
      }
    } catch (...) {}
  }

  python::list children;
  for (auto it = q->beginChildren(); it != q->endChildren(); ++it) {
    children.append(walkBondQuery(it->get()));
  }
  if (python::len(children) > 0) {
    d["children"] = children;
  }
  return d;
}

/**
 * Parse a SMARTS string and return a dict containing all atom and bond
 * query trees.  This is the main entry point for the visualiser frontend.
 *
 * Returns: {
 *   "smarts": <canonical SMARTS>,
 *   "atoms": [
 *     { "idx": 0, "atom_map": 0, "smarts_token": "[#6]",
 *       "query_tree": { ... }, "neighbors": [1, 2] },
 *     ...
 *   ],
 *   "bonds": [
 *     { "idx": 0, "begin_atom_idx": 0, "end_atom_idx": 1,
 *       "bond_query": { ... } },
 *     ...
 *   ]
 * }
 */
python::dict querymolFromSmarts(const std::string &smarts) {
  python::dict result;
  auto *mol = RDKit::SmartsToMol(smarts);
  if (!mol) {
    throw std::runtime_error("Invalid SMARTS: " + smarts);
  }
  std::unique_ptr<RDKit::ROMol> mol_owner(mol);

  // Property cache and ring info
  try { mol->updatePropertyCache(false); } catch (...) {}
  RDKit::MolOps::fastFindRings(*mol);

  result["smarts"] = RDKit::MolToSmarts(*mol);

  // --- Atoms ---
  python::list atoms;
  for (const auto *atom : mol->atoms()) {
    python::dict ad;
    ad["idx"] = static_cast<int>(atom->getIdx());
    ad["atom_map"] = atom->getAtomMapNum();

    // Get this atom's SMARTS token
    if (atom->hasQuery()) {
      ad["smarts_token"] = RDKit::SmartsWrite::GetAtomSmarts(
          static_cast<const RDKit::QueryAtom *>(atom));
    } else {
      ad["smarts_token"] = std::string("?");
    }

    // Walk query tree
    const auto *query = atom->getQuery();
    ad["query_tree"] = walkAtomQuery(
        static_cast<const ATOM_QUERY *>(query));

    // Neighbors
    python::list nbrs;
    for (const auto *nbr : mol->atomNeighbors(atom)) {
      nbrs.append(static_cast<int>(nbr->getIdx()));
    }
    ad["neighbors"] = nbrs;

    atoms.append(ad);
  }
  result["atoms"] = atoms;

  // --- Bonds ---
  python::list bonds;
  for (const auto *bond : mol->bonds()) {
    python::dict bd;
    bd["idx"] = static_cast<int>(bond->getIdx());
    bd["begin_atom_idx"] = static_cast<int>(bond->getBeginAtomIdx());
    bd["end_atom_idx"] = static_cast<int>(bond->getEndAtomIdx());

    if (bond->hasQuery()) {
      bd["bond_query"] = walkBondQuery(
          static_cast<const BOND_QUERY *>(bond->getQuery()));
    } else {
      python::dict bq;
      bq["description"] = "BondOrder";
      bq["op"] = "leaf";
      bq["value"] = static_cast<int>(bond->getBondType());
      bq["negated"] = false;
      bd["bond_query"] = bq;
    }
    bonds.append(bd);
  }
  result["bonds"] = bonds;

  return result;
}

std::vector<std::string> smilesInputToVector(const python::object &smiles_or_list,
                                             bool *was_single = nullptr) {
  std::vector<std::string> smiles_list;
  const bool single = PyUnicode_Check(smiles_or_list.ptr());
  if (was_single) {
    *was_single = single;
  }

  if (single) {
    smiles_list.push_back(python::extract<std::string>(smiles_or_list));
    return smiles_list;
  }

  python::stl_input_iterator<std::string> begin(smiles_or_list), end;
  for (auto it = begin; it != end; ++it) {
    smiles_list.push_back(*it);
  }
  return smiles_list;
}

std::vector<std::string> primitivesToVector(const python::object &primitives);

python::object smilesToSmartsByPrimitives(const python::object &smiles_or_list,
                                          const python::object &primitives) {
  std::vector<std::string> primitive_list;
  if (PyUnicode_Check(primitives.ptr())) {
    primitive_list.push_back(python::extract<std::string>(primitives));
  } else {
    python::stl_input_iterator<std::string> begin(primitives), end;
    for (auto it = begin; it != end; ++it) {
      primitive_list.push_back(*it);
    }
  }

  if (primitive_list.empty()) {
    PyErr_SetString(PyExc_ValueError,
                    "primitives must be a non-empty string or iterable of strings");
    python::throw_error_already_set();
  }

  bool was_single = false;
  const auto smiles_list = smilesInputToVector(smiles_or_list, &was_single);
  const auto out = atom_typer::smiles_to_smarts(smiles_list, primitive_list);
  if (was_single) {
    return python::object(out.front());
  }
  return vectorStringToList(out);
}

python::object smilesToAtomCenteredSmartsByPrimitives(
    const python::object &smiles_or_list, const python::object &primitives,
  unsigned int radius = 0, bool wildcardNeighbors = false,
  bool includePrimitiveSubsets = false, bool deduplicate = false) {
  const auto primitive_list = primitivesToVector(primitives);

  bool was_single = false;
  const auto smiles_list = smilesInputToVector(smiles_or_list, &was_single);
  const auto out = atom_typer::smiles_to_atom_centered_smarts(
      smiles_list, primitive_list, radius, wildcardNeighbors,
      includePrimitiveSubsets, deduplicate);

  if (was_single) {
    if (out.empty()) {
      return python::list();
    }
    return vectorStringToList(out.front());
  }

  if (deduplicate) {
    std::vector<std::string> flattened;
    size_t totalSize = 0;
    for (const auto &perMol : out) {
      totalSize += perMol.size();
    }
    flattened.reserve(totalSize);
    for (const auto &perMol : out) {
      flattened.insert(flattened.end(), perMol.begin(), perMol.end());
    }
    return vectorStringToList(flattened);
  }

  return vectorVectorStringToList(out);
}

std::vector<std::string> primitivesToVector(const python::object &primitives) {
  std::vector<std::string> primitive_list;
  if (PyUnicode_Check(primitives.ptr())) {
    primitive_list.push_back(python::extract<std::string>(primitives));
  } else {
    python::stl_input_iterator<std::string> begin(primitives), end;
    for (auto it = begin; it != end; ++it) {
      primitive_list.push_back(*it);
    }
  }
  if (primitive_list.empty()) {
    PyErr_SetString(PyExc_ValueError,
                    "primitives must be a non-empty string or iterable of strings");
    python::throw_error_already_set();
  }
  return primitive_list;
}

python::list radiusTemplatesToPy(
    const std::vector<atom_typer::RadiusTemplate> &results) {
  python::list out;
  for (const auto &entry : results) {
    out.append(python::make_tuple(std::get<0>(entry), std::get<1>(entry)));
  }
  return out;
}

python::list extractSingleRootTemplateFromReaction(
    RDKit::ChemicalReaction &rxn, const python::object &primitives,
    unsigned int max_radius = 3, bool verbose = false) {
  const auto primitive_list = primitivesToVector(primitives);
  const auto results = atom_typer::extract_single_root_template(
      &rxn, primitive_list, max_radius, verbose);
  return radiusTemplatesToPy(results);
}

python::list extractSingleRootTemplateFromText(
    const std::string &reaction_text, const python::object &primitives,
    unsigned int max_radius = 3, bool verbose = false,
    bool use_smiles = true) {
  const auto primitive_list = primitivesToVector(primitives);
  const auto results = atom_typer::extract_single_root_template(
      reaction_text, primitive_list, max_radius, verbose, use_smiles);
  return radiusTemplatesToPy(results);
}

}  // namespace

BOOST_PYTHON_MODULE(rdAtomTyper) {
  python::scope().attr("__doc__") =
      "AtomTyper: type atoms in SMILES/SMARTS strings based on atomic primitives "
      "and local connectivity.";

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

  python::def("smiles_to_smarts", &smilesToSmartsByPrimitives,
              (python::arg("smiles"), python::arg("primitives")),
              "Convert SMILES to SMARTS using custom primitives, e.g. "
              "['X','D','R'] or '[charge, R, D]'.");

  python::def("smiles_to_atom_centered_smarts",
              &smilesToAtomCenteredSmartsByPrimitives,
              (python::arg("smiles"), python::arg("primitives"),
               python::arg("radius") = 0,
               python::arg("wildcardNeighbors") = false,
               python::arg("includePrimitiveSubsets") = false,
               python::arg("deduplicate") = false),
              "Generate per-atom SMARTS neighborhoods rooted at each atom. "
              "Accepts a SMILES string or iterable of SMILES strings. "
              "If wildcardNeighbors=True, all non-center atoms are emitted as [*]. "
              "If includePrimitiveSubsets=True, emits all non-empty subsets of rooted atom primitives. "
              "If deduplicate=True, duplicate SMARTS are removed while preserving first occurrence order. "
              "For iterable input with deduplicate=True, deduplication is global across the whole batch "
              "and a single flattened list is returned. "
              "Otherwise returns one SMARTS fragment list per input molecule.");

  python::def("extract_single_root_template", &extractSingleRootTemplateFromReaction,
              (python::arg("reaction"), python::arg("primitives"),
               python::arg("max_radius") = 3, python::arg("verbose") = false),
              "Extract radius-indexed reaction core templates from a ChemicalReaction.");

  python::def("extract_single_root_template", &extractSingleRootTemplateFromText,
              (python::arg("reaction_text"), python::arg("primitives"),
               python::arg("max_radius") = 3, python::arg("verbose") = false,
               python::arg("use_smiles") = true),
              "Extract radius-indexed reaction core templates from reaction text.");

    python::def(
      "type_smiles", &typeSmiles,
      (python::arg("smiles_or_list"), python::arg("use_canonical") = true,
       python::arg("reserialize") = true),
      "Type SMILES via AtomTyper::type_atoms_from_smiles. Accepts a single SMILES string or an iterable of SMILES strings. If reserialize=True, returns a typed SMARTS-like string reconstructed from PatternItem atom/bond tokens instead of AtomType objects.");

    python::def(
      "querymol_from_smarts", &querymolFromSmarts,
      (python::arg("smarts")),
      "Parse a SMARTS string and return a dict containing atom and bond\n"
      "query trees.  Each atom entry includes the full query tree as a\n"
      "nested dict with keys: description, op (and/or/not/leaf/recursive),\n"
      "negated, value (for leaf nodes), and children.\n\n"
      "Returns a dict with keys: smarts, atoms, bonds.\n");
}
