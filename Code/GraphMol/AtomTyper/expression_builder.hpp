#ifndef EXPRESSION_BUILDER_HPP
#define EXPRESSION_BUILDER_HPP

#include <RDGeneral/export.h>
#include <string>
#include <vector>

// Forward declarations
namespace RDKit {
    class Atom;
    class ROMol;
    class RWMol;
}

namespace atom_typer {

/**
 * Detail levels for SMARTS pattern generation
 */
enum Level { 
    MINIMAL,   // Basic element types only
    STANDARD,  // Elements + bond types + aromaticity
    DETAILED,  // Add hydrogen counts, formal charges
    COMPLETE   // All features: valence, connectivity, ring membership
};

/**
 * Convert a SMILES string to a SMARTS pattern
 * 
 * @param smiles Input SMILES string representing a molecule
 * @param level Detail level for the generated SMARTS pattern
 * @return SMARTS string with query features based on the specified level
 */
RDKIT_ATOMTYPER_EXPORT std::string smiles_to_smarts(const std::string& smiles, Level level);

/**
 * Batch convert SMILES strings to SMARTS patterns using a predefined detail level.
 *
 * @param smilesList Input SMILES strings
 * @param level Detail level for generated SMARTS patterns
 * @return SMARTS strings in the same order as smilesList
 */
RDKIT_ATOMTYPER_EXPORT std::vector<std::string> smiles_to_smarts(
    const std::vector<std::string>& smilesList, Level level);

/**
 * Convert a SMILES string to a SMARTS pattern using a custom primitive list.
 *
 * Primitive names may be passed either as individual
 * entries (e.g. {"X", "D", "R"}) or bracketed/comma-delimited groups
 * (e.g. {"[X,D,R]"} or {"[charge, R, D]"}).
 *
 * Supported primitives:
 * - a / AtomIsAromatic: aromaticity (a or !a)
 * - A / AtomIsAliphatic: aliphaticity (A or !A)
 * - D / AtomExplicitDegree: explicit degree
 * - H / AtomHCount: total hydrogen count
 * - h / AtomImplicitHCount: implicit hydrogen count
 * - R / AtomInNRings / atom_ring_count: number of rings containing atom
 * - r / AtomMinRingSize: minimum ring size containing atom
 * - charge / formalcharge: formal charge (+1, -1, +0)
 * - v / AtomTotalValence: total valence
 * - X / AtomTotalDegree: total degree
 * - x / AtomRingBondCount: number of ring bonds at atom
 * - atomic_number / AtomAtomicNum: atomic number (#6, #7, ...)
 *
 * Extra convenience aliases are also supported (for example "element" and
 * "atom") to generate element-symbol-based atom typing.
 */
RDKIT_ATOMTYPER_EXPORT std::string smiles_to_smarts(
    const std::string& smiles, const std::vector<std::string>& primitiveList);

/**
 * Batch convert SMILES strings to SMARTS patterns using a shared primitive list.
 *
 * @param smilesList Input SMILES strings
 * @param primitiveList Primitive list using same API as smiles_to_smarts()
 * @return SMARTS strings in the same order as smilesList
 */
RDKIT_ATOMTYPER_EXPORT std::vector<std::string> smiles_to_smarts(
    const std::vector<std::string>& smilesList,
    const std::vector<std::string>& primitiveList);

/**
 * Generate per-atom SMARTS neighborhoods from a SMILES string.
 *
 * For each atom in the molecule, returns a SMARTS fragment rooted at that atom
 * and including all atoms/bonds within the specified graph radius.
 * Atom typing is controlled by primitiveList (same API as smiles_to_smarts()).
 *
 * @param smiles Input SMILES string
 * @param primitiveList Primitive list using same API as smiles_to_smarts()
 * @param radius Graph radius around each atom center (0 = center atom only)
 * @return SMARTS fragment per atom in atom index order
 */
RDKIT_ATOMTYPER_EXPORT std::vector<std::string> smiles_to_atom_centered_smarts(
    const std::string& smiles, const std::vector<std::string>& primitiveList,
    unsigned int radius);

/**
 * Batch generate per-atom SMARTS neighborhoods for multiple SMILES strings.
 *
 * @param smilesList Input SMILES strings
 * @param primitiveList Primitive list using same API as smiles_to_smarts()
 * @param radius Graph radius around each atom center (0 = center atom only)
 * @return Per-molecule list of per-atom SMARTS fragments
 */
RDKIT_ATOMTYPER_EXPORT std::vector<std::vector<std::string>>
smiles_to_atom_centered_smarts(
    const std::vector<std::string>& smilesList,
    const std::vector<std::string>& primitiveList, unsigned int radius);

/**
 * Normalize primitive list input into ordered primitive tokens.
 *
 * Input entries can be individual tokens or bracketed/comma-delimited groups.
 */
RDKIT_ATOMTYPER_EXPORT std::vector<std::string> normalize_primitive_list(
    const std::vector<std::string>& primitiveList);

/**
 * Build a single atom SMARTS token using expression_builder primitive logic.
 *
 * @param atom Atom to featurize
 * @param mol Parent molecule used for topology/ring primitives
 * @param primitiveList Primitive list using same API as smiles_to_smarts()
 * @param includeAtomMap Include :map in output when atom has map number
 * @return Bracket atom SMARTS token, e.g. [#6;D3;H1:4]
 */
RDKIT_ATOMTYPER_EXPORT std::string build_atom_primitive_token(
    const RDKit::Atom* atom, const RDKit::ROMol* mol,
    const std::vector<std::string>& primitiveList,
    bool includeAtomMap = false);

/**
 * Convert a query molecule to a SMARTS string
 * 
 * @param mol Original molecule (for accessing atom properties)
 * @param queryMol Query molecule built with appropriate query features
 * @param level Detail level for formatting the SMARTS string
 * @return Formatted SMARTS string
 */
RDKIT_ATOMTYPER_EXPORT std::string queryMoleculeToSmarts(const RDKit::ROMol* mol, 
                                   const RDKit::RWMol* queryMol,
                                   Level level);

} // namespace atom_typer

#endif // EXPRESSION_BUILDER_HPP
