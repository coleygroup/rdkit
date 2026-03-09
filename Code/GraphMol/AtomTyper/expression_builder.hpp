#ifndef EXPRESSION_BUILDER_HPP
#define EXPRESSION_BUILDER_HPP

#include <RDGeneral/export.h>
#include <string>
#include <vector>

// Forward declarations
namespace RDKit {
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
