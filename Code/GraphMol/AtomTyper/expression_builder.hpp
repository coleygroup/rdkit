#ifndef EXPRESSION_BUILDER_HPP
#define EXPRESSION_BUILDER_HPP

#include <RDGeneral/export.h>
#include <string>

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
