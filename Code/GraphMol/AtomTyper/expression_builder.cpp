#include "expression_builder.hpp"
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/PeriodicTable.h>
#include <stdexcept>
#include <string>
#include <sstream>

namespace atom_typer {

// Helper function to get element symbol from atomic number
std::string getElementSymbol(int atomicNum) {
    return RDKit::PeriodicTable::getTable()->getElementSymbol(atomicNum);
}

/**
 * Convert a query molecule to a SMARTS string based on detail level
 * 
 * @param mol Original molecule (for accessing atom properties)
 * @param queryMol Query molecule built with appropriate query features
 * @param level Detail level for formatting the SMARTS string
 * @return Formatted SMARTS string
 */
std::string queryMoleculeToSmarts(const RDKit::ROMol* mol, 
                                   const RDKit::RWMol* queryMol,
                                   Level level) {
    std::stringstream ss;
    
    if (level == Level::MINIMAL) {
        // Format: [C][C][O]
        for (size_t i = 0; i < mol->getNumAtoms(); ++i) {
            const auto atom = mol->getAtomWithIdx(i);
            ss << "[" << getElementSymbol(atom->getAtomicNum()) << "]";
        }
    }
    else if (level == Level::STANDARD) {
        // Format: [C;D1][C;D2][O;D1]
        for (size_t i = 0; i < mol->getNumAtoms(); ++i) {
            const auto atom = mol->getAtomWithIdx(i);
            ss << "[" << getElementSymbol(atom->getAtomicNum()) 
               << ";D" << atom->getDegree() << "]";
        }
    }
    else if (level == Level::DETAILED) {
        // Format: [C;D1;H3][C;D2;H2][O;D1;H1]
        for (size_t i = 0; i < mol->getNumAtoms(); ++i) {
            const auto atom = mol->getAtomWithIdx(i);
            ss << "[" << getElementSymbol(atom->getAtomicNum())
               << ";D" << atom->getDegree()
               << ";H" << atom->getTotalNumHs(true) << "]";
        }
    }
    else if (level == Level::COMPLETE) {
        // Use RDKit's MolToSmarts for complete level
        return RDKit::MolToSmarts(*queryMol);
    }
    
    return ss.str();
}

/**
 * Convert a SMILES string to a SMARTS pattern with specified detail level
 * 
 * @param smiles Input SMILES string
 * @param level Detail level for the SMARTS pattern
 * @return SMARTS string representation
 */
std::string smiles_to_smarts(const std::string& smiles, Level level) {
    /*
    Parse SMILES string to RDKit molecule
    Based on level, add appropriate query features:
        - MINIMAL: Just element types
        - STANDARD: Elements + bond types + aromaticity
        - DETAILED: Add hydrogen counts, charges
        - COMPLETE: Add all features (valence, connectivity, ring membership)
    Convert to SMARTS string
    */ 

    // For an empty SMILES string, throw an exception
    if (smiles.empty()) {
        throw std::invalid_argument("Input SMILES string is empty");
    }

    // Parse SMILES to RDKit molecule
    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(smiles));

    // if not valid molecule, throw an exception
    if (!mol) {
        throw std::invalid_argument("Invalid SMILES string");
    }

    // Create a new RWMol (editable molecule) to build our query molecule -> SMARTS
    std::unique_ptr<RDKit::RWMol> queryMol(new RDKit::RWMol());

    // For MINIMAL level, smiles "CCO" -> "[C][C][O]" 
    // Each element is represented only by its atomic symbol
    if (level == Level::MINIMAL) {
        /*
        SMARTS string in the form of only atomic symbols

        Create new molecule with only atomic numbers preserved
        Use QueryAtom objects for atoms
        Include bonds with UNSPECIFIED type
        */

        for (const auto atom : mol->atoms()) {

            // Get the atomic number and add to query molecule
            int atomicNum = atom->getAtomicNum();
            RDKit::QueryAtom* queryAtom = new RDKit::QueryAtom(atom->getAtomicNum());
            queryMol->addAtom(queryAtom, false, true);
        }

        // Add bonds and connectivity, but no specific bond types
        for (auto bond : mol->bonds()) {
            queryMol->addBond(bond->getBeginAtomIdx(), bond->getEndAtomIdx(), 
            RDKit::Bond::UNSPECIFIED);
        }
    }

    else if (level == Level::STANDARD) {
        /*
        SMARTS string with elements, aromaticity, and bond types

        Creates new molecule, get atomic number
        Get aromaticity from original molecule
        Add bond types (single, double, aromatic, etc.)
        */
        for (const auto atom : mol->atoms()) {

            // Get the atomic number and add to query molecule
            int atomicNum = atom->getAtomicNum();
            RDKit::QueryAtom* queryAtom = new RDKit::QueryAtom(atom->getAtomicNum());

            //Check if it is aromatic
            if (atom->getIsAromatic()) {
                queryAtom->setIsAromatic(true);
            }

            queryMol->addAtom(queryAtom, false, true);
        }

        for (auto bond : mol->bonds()) {

            // Get bond type and add to query molecule
            RDKit::Bond::BondType bondType = bond->getBondType();
            queryMol->addBond(bond->getBeginAtomIdx(), bond->getEndAtomIdx(), bondType);

            // Check bond aromaticity
            if (bond->getIsAromatic()) {
                queryMol->getBondWithIdx(queryMol->getNumBonds() - 1)->setIsAromatic(true);
            }
        }
    }

    else if (level == Level::DETAILED) {
        /*
        Create SMARTS string with elements, aromaticity, bond types,
        explicit hydrogen counts, and formal charges

        Add specific bond types as in STANDARD level
        */

        for (const auto atom :  mol->atoms()) {

            // Get the atomic number and add to query molecule
            int atomicNum = atom->getAtomicNum();
            RDKit::QueryAtom* queryAtom = new RDKit::QueryAtom(atom->getAtomicNum());

            // Check for aromaticity
            if (atom->getIsAromatic()) {
                queryAtom->setIsAromatic(true);
            }

            // Add explicit hydrogen count query
            int explicitHCount = atom->getTotalNumHs(true);

            queryAtom->expandQuery(
                RDKit::makeAtomHCountQuery(explicitHCount),
                Queries::COMPOSITE_AND);

            // Add formal charge
            int formalCharge = atom->getFormalCharge();
            // Check for formal charge
            if (formalCharge != 0) {
                queryAtom->expandQuery(
                    RDKit::makeAtomFormalChargeQuery(formalCharge),
                    Queries::COMPOSITE_AND);
            }
            // Get the atomic number and add to query molecule
            queryMol->addAtom(queryAtom, false, true);
        }

        for (auto bond : mol->bonds()) {

            // Get bond type and add to query molecule same as STANDARD level
            RDKit::Bond::BondType bondType = bond->getBondType();
            queryMol->addBond(bond->getBeginAtomIdx(), bond->getEndAtomIdx(), bondType);
            if (bond->getIsAromatic()) {
                queryMol->getBondWithIdx(queryMol->getNumBonds() - 1)->setIsAromatic(true);
            }
        }
    }

    else if (level == Level::COMPLETE) {
        /*
        Create SMARTS string with atomic number, aromaticity, bond types,
        explicit hydrogen counts, formal charges, valence, connectivity,
        ring membership, hybridization

        Add bonds with specific types like STANDARD and DETAILED levels
        */

        for (const auto atom : mol->atoms()) {

            // Get the atomic number and add to query molecule
            int atomicNum = atom->getAtomicNum();
            RDKit::QueryAtom* queryAtom = new RDKit::QueryAtom(atom->getAtomicNum());

            // Check for aromaticity
            if (atom->getIsAromatic()) {
                queryAtom->setIsAromatic(true);
            }

            // Add explicit hydrogen count query
            int explicitHCount = atom->getTotalNumHs(true);
            queryAtom->expandQuery(
                RDKit::makeAtomHCountQuery(explicitHCount),
                Queries::COMPOSITE_AND);

            // Add formal charge
            int formalCharge = atom->getFormalCharge();
            if (formalCharge != 0) {
                queryAtom->expandQuery(
                    RDKit::makeAtomFormalChargeQuery(formalCharge),
                    Queries::COMPOSITE_AND);
            }

            // Add total valence query
            int totalValence = atom->getTotalValence();
            queryAtom->expandQuery(
                RDKit::makeAtomTotalValenceQuery(totalValence),
                Queries::COMPOSITE_AND);

            // Add ring membership query
            RDKit::RingInfo* ringInfo = mol->getRingInfo();
            bool isInRing = ringInfo->numAtomRings(atom->getIdx()) > 0;
            queryAtom->expandQuery(
                RDKit::makeAtomInRingQuery(),
                Queries::COMPOSITE_AND,
                !isInRing);

            // Add hybridization query
            RDKit::Atom::HybridizationType hybridization = atom->getHybridization();
            if (hybridization != RDKit::Atom::UNSPECIFIED) {
                queryAtom->expandQuery(
                RDKit::makeAtomHybridizationQuery(hybridization),
                Queries::COMPOSITE_AND);
            }
            
            // Get the atomic number and add to query molecule
            queryMol->addAtom(queryAtom, false, true);
        }

        for (auto bond : mol->bonds()) {

            // Get bond type and add to query molecule same as STANDARD and DETAILED levels
            RDKit::Bond::BondType bondType = bond->getBondType();
            queryMol->addBond(bond->getBeginAtomIdx(), bond->getEndAtomIdx(), bondType);
            if (bond->getIsAromatic()) {
                queryMol->getBondWithIdx(queryMol->getNumBonds() - 1)->setIsAromatic(true);
            }
        }
    }

    // Convert query molecule to SMARTS string
    return queryMoleculeToSmarts(mol.get(), queryMol.get(), level);
}

} // namespace atom_typer

