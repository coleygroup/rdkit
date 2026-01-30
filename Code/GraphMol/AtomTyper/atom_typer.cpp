#include "atom_typer.hpp"
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <sstream>
#include <stdexcept>

namespace atom_typer {

/**
 * Private implementation class (PIMPL pattern)
 */
class AtomTyper::Impl {
public:
    bool use_canonical = false;
    
    /**
     * Extract atom type information from an RDKit atom
     */
    AtomType extract_atom_type(RDKit::Atom* atom, RDKit::ROMol* mol, int idx) {
        AtomType at;
        at.atom_idx = idx;
        at.atomic_number = atom->getAtomicNum();
        at.formal_charge = atom->getFormalCharge();
        at.num_hydrogens = atom->getTotalNumHs();
        at.degree = atom->getDegree();
        at.valence = atom->getTotalValence();
        at.is_aromatic = atom->getIsAromatic();
        at.is_in_ring = mol->getRingInfo()->numAtomRings(idx) > 0;
        
        // Initialize ring-related fields
        at.ring_size = 0;
        at.num_ring_bonds = 0;
        at.num_aliphatic_rings = 0;
        at.num_aromatic_rings = 0;
        at.ring_connectivity = 0;
        // ring_connectivity counts neighbors that are also in rings
        for (const auto& nbr : mol->atomNeighbors(atom)) {
            if (mol->getRingInfo()->numAtomRings(nbr->getIdx()) > 0) {
                at.ring_connectivity++;
            }
        }
        
        //find number of aliphatic/aromatic rings
        //start both counts at 0
        if (!atom -> getIsAromatic()) {
            //the atom is aliphatic, so all rings aliphatic
            at.num_aliphatic_rings = mol -> getRingInfo() -> numAtomRings(idx);
            at.num_aromatic_rings = 0;
        }
        else {
            at.num_aliphatic_rings = 0;
            at.num_aromatic_rings = 0;

            for (const auto& ring : mol -> getRingInfo() -> atomRings()) {
                if (std::find(ring.begin(), ring.end(), idx) != ring.end()) {
                    //check if ring is aromatic
                    bool is_aromatic_ring = true;
                    for (int atom_id : ring) {
                        RDKit::Atom* ring_atom = mol -> getAtomWithIdx(atom_id);
                        if (!ring_atom -> getIsAromatic()) {
                            is_aromatic_ring = false;
                            break;
                        }
                    }

                    //add if aromatic or aliphatic
                    if (is_aromatic_ring) {
                        at.num_aromatic_rings += 1;
                    }
                    else {
                        at.num_aliphatic_rings += 1;
                    }
                }
            }
        }

        // Get smallest ring size
        // Note: ring_size is stored as int to match RDKit's AtomType convention
        // and to maintain API consistency. Max ring size in practice is << INT_MAX.
        at.ring_size = 0;

        //declare smallest ring
        std::vector<int> smallest_ring;
        if (at.is_in_ring) {
            const auto& ring_info = mol->getRingInfo();
            for (const auto& ring : ring_info->atomRings()) {
                if (std::find(ring.begin(), ring.end(), idx) != ring.end()) {

                    //find the smallest ring/size
                    int current_ring_size = static_cast<int>(ring.size());
                    if (at.ring_size == 0 || current_ring_size < at.ring_size) {
                        at.ring_size = current_ring_size;
                        smallest_ring = ring;
                    }
                }
            }
        }
        //use smallest_ring to get the member atoms' indices
        for (int ind : smallest_ring) {
            at.ring_membership_list.push_back(ind);
        }

        //find number of ring bonds (bonds that are in ANY ring, not just smallest)
        at.num_ring_bonds = 0;
        for (const auto& bond : mol->atomBonds(atom)) {
            if (mol->getRingInfo()->numBondRings(bond->getIdx()) > 0) {
                at.num_ring_bonds++;
            }
        }       
        
        // Get hybridization
        switch (atom->getHybridization()) {
            case RDKit::Atom::SP:
                at.hybridization = "SP";
                break;
            case RDKit::Atom::SP2:
                at.hybridization = "SP2";
                break;
            case RDKit::Atom::SP3:
                at.hybridization = "SP3";
                break;
            case RDKit::Atom::SP3D:
                at.hybridization = "SP3D";
                break;
            case RDKit::Atom::SP3D2:
                at.hybridization = "SP3D2";
                break;
            default:
                at.hybridization = "UNSPECIFIED";
        }

        // Get Chirality using R/S notation
        // Use RDKit's CIP assignment to get R/S labels
        std::string cip_code;
        if (atom->hasProp("_CIPCode")) {
            atom->getProp("_CIPCode", cip_code);
            at.chirality = cip_code;  // Will be "R" or "S"
        } else {
            // No CIP code assigned
            switch (atom->getChiralTag()) {
                case RDKit::Atom::CHI_TETRAHEDRAL_CW:
                case RDKit::Atom::CHI_TETRAHEDRAL_CCW:
                    at.chirality = "UNSPECIFIED";  // Has chirality but no R/S assignment
                    break;
                default:
                    at.chirality = "UNSPECIFIED";
            }
        }
        
        // Get neighbors & bond types 
        std::map<int, int> bond_type_counts;
        for (const auto& nbr : mol->atomNeighbors(atom)) {
            at.neighbors.push_back(nbr->getIdx());

            RDKit::Bond* bond = mol->getBondBetweenAtoms(idx, nbr->getIdx());
            if (bond) {
                int bond_type_int = static_cast<int>(bond->getBondType());
                bond_type_counts[bond_type_int]++;
            }
        }
        at.bond_types = bond_type_counts;
        
        // Generate SMARTS pattern for this atom type
        at.smarts_pattern = generate_smarts_pattern(at);
        
        return at;
    }
    
    /**
     * Generate a SMARTS pattern based on atom type
     */
    std::string generate_smarts_pattern(const AtomType& at) {
        std::stringstream ss;
        ss << "[";
        
        // Add atomic number
        if (at.atomic_number == 6) ss << "C";
        else if (at.atomic_number == 7) ss << "N";
        else if (at.atomic_number == 8) ss << "O";
        else if (at.atomic_number == 16) ss << "S";
        else if (at.atomic_number == 1) ss << "H";
        else ss << "#" << at.atomic_number;
        
        // Add degree
        if (at.degree > 0) {
            ss << "D" << at.degree;
        }
        
        // Add hydrogen count
        if (at.num_hydrogens > 0) {
            ss << "H" << at.num_hydrogens;
        }
        
        // Add charge
        if (at.formal_charge != 0) {
            if (at.formal_charge > 0) {
                ss << "+" << at.formal_charge;
            } else {
                ss << at.formal_charge;
            }
        }
        
        // Add aromaticity
        if (at.is_aromatic) {
            ss << "a";
        }
        
        // Add ring membership
        if (at.is_in_ring) {
            ss << "R";
            if (at.ring_size > 0) {
                ss << at.ring_size;
            }
        }
        
        ss << "]";
        return ss.str();
    }
    
    /**
     * Process a molecule and extract all atom types
     */
    std::vector<AtomType> process_molecule(RDKit::ROMol* mol) {
        if (!mol) {
            throw std::runtime_error("Invalid molecule");
        }
        
        // Assign CIP (R/S) stereochemistry codes to atoms
        // Parameters: cleanIt=true, force=true, flagPossibleStereoCenters=false
        RDKit::MolOps::assignStereochemistry(*mol, true, true, true);
        
        std::vector<AtomType> atom_types;
        atom_types.reserve(mol->getNumAtoms());
        
        for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
            auto* atom = mol->getAtomWithIdx(i);
            atom_types.push_back(extract_atom_type(atom, mol, i));
        }
        
        return atom_types;
    }
};

// AtomTyper implementation
AtomTyper::AtomTyper() : pimpl(std::make_unique<Impl>()) {}

AtomTyper::~AtomTyper() = default;

std::vector<AtomType> AtomTyper::type_atoms_from_smiles(const std::string& smiles) {
    if (smiles.empty()) {
        throw std::runtime_error("Empty SMILES string");
    }
    
    // Parse SMILES
    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(smiles));
    if (!mol) {
        throw std::runtime_error("Failed to parse SMILES: " + smiles);
    }
    
    // Apply canonicalization if enabled
    if (pimpl->use_canonical) {
        // Use RDKit's canonical SMILES
        std::string canonical = RDKit::MolToSmiles(*mol);
        mol.reset(RDKit::SmilesToMol(canonical));
    }
    
    // Process molecule and extract atom types
    auto atom_types = pimpl->process_molecule(mol.get());
    
    return atom_types;
}

std::vector<AtomType> AtomTyper::type_atoms_from_smarts(const std::string& smarts) {
    if (smarts.empty()) {
        throw std::runtime_error("Empty SMARTS string");
    }
    
    // Parse SMARTS
    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmartsToMol(smarts));
    if (!mol) {
        throw std::runtime_error("Failed to parse SMARTS: " + smarts);
    }
    
    // Update property cache and find rings for SMARTS molecules
    mol->updatePropertyCache();
    RDKit::MolOps::findSSSR(*mol);
    
    return pimpl->process_molecule(mol.get());
}

std::string AtomTyper::get_atom_types_string(const std::vector<AtomType>& atom_types) {
    std::stringstream ss;
    ss << "Atom Types:\n";
    ss << "===========\n";
    
    for (const auto& at : atom_types) {
        ss << "Atom " << at.atom_idx << ": "
           << "Element=" << at.atomic_number
           << ", Charge=" << at.formal_charge
           << ", H=" << at.num_hydrogens
           << ", Degree=" << at.degree
           << ", Valence=" << at.valence
           << ", Aromatic=" << (at.is_aromatic ? "Yes" : "No")
           << ", InRing=" << (at.is_in_ring ? "Yes" : "No")
           << ", RingSize=" << at.ring_size
           << ", Hybrid=" << at.hybridization
           << ", SMARTS=" << at.smarts_pattern
           << "\n";
    }
    
    return ss.str();
}

void AtomTyper::set_use_canonical(bool use_canonical) {
    pimpl->use_canonical = use_canonical;
}

} // namespace atom_typer
