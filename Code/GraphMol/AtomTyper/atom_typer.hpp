#ifndef ATOM_TYPER_HPP
#define ATOM_TYPER_HPP

#include <RDGeneral/export.h>
#include <string>
#include <vector>
#include <map>
#include <memory>

namespace atom_typer {

/**
 * Structure to hold atom type information
 */
struct AtomType {
    int atom_idx;               // Index of the atom in the molecule
    int atomic_number;          // Atomic number (element)
    int formal_charge;          // Formal charge
    int num_hydrogens;          // Number of attached hydrogens
    int degree;                 // Degree (number of bonds)
    int valence;                // Valence
    bool is_aromatic;           // Aromaticity flag
    bool is_in_ring;            // Ring membership
    int ring_size;              // Size of smallest ring (0 if not in ring)
    std::string hybridization;  // Hybridization (SP, SP2, SP3, etc.)
    std::vector<int> neighbors; // Indices of neighboring atoms
    std::string smarts_pattern; // SMARTS pattern for this atom type
    std::vector<int> ring_membership_list;      // List of members in ring
    int num_ring_bonds;                         // number of bonds in ring
    int num_aliphatic_rings;                    // number of aliphatic rings
    int num_aromatic_rings;                     // number of aromatic rings
    std::string chirality;                      // Chirality information
    int ring_connectivity;                      // Ring connectivity
    std::map<int, int> bond_types;              // Map of bond types to counts  
};

/**
 * Main AtomTyper class for typing atoms in SMILES/SMARTS strings
 */
class RDKIT_ATOMTYPER_EXPORT AtomTyper {
public:
    /**
     * Constructor
     */
    AtomTyper();
    
    /**
     * Destructor
     */
    ~AtomTyper();
    
    /**
     * Type all atoms in a SMILES string
     * @param smiles The input SMILES string
     * @return Vector of AtomType structures for each atom
     */
    std::vector<AtomType> type_atoms_from_smiles(const std::string& smiles);
    
    /**
     * Type all atoms in a SMARTS string
     * @param smarts The input SMARTS string
     * @return Vector of AtomType structures for each atom
     */
    std::vector<AtomType> type_atoms_from_smarts(const std::string& smarts);
    
    /**
     * Get atom types as a formatted string
     * @param atom_types Vector of atom types
     * @return Formatted string representation
     */
    std::string get_atom_types_string(const std::vector<AtomType>& atom_types);
    
    /**
     * Set whether to use canonical SMILES
     * @param use_canonical Flag to enable/disable canonicalization
     */
    void set_use_canonical(bool use_canonical);

private:
    class Impl;
    std::unique_ptr<Impl> pimpl;
};

} // namespace atom_typer

#endif // ATOM_TYPER_HPP
