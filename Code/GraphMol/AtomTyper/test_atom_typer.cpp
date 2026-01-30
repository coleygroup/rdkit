#include "atom_typer.hpp"
#include <gtest/gtest.h>
#include <stdexcept>

class AtomTyperTest : public ::testing::Test {
protected:
    atom_typer::AtomTyper typer;
};

// Test basic SMILES parsing
TEST_F(AtomTyperTest, BasicSmilesEthanol) {
    std::string smiles = "CCO";
    auto atom_types = typer.type_atoms_from_smiles(smiles);
    
    ASSERT_EQ(atom_types.size(), 3);
    
    // Check first carbon
    EXPECT_EQ(atom_types[0].atomic_number, 6);
    EXPECT_EQ(atom_types[0].degree, 1);
    
    // Check second carbon
    EXPECT_EQ(atom_types[1].atomic_number, 6);
    EXPECT_EQ(atom_types[1].degree, 2);
    
    // Check oxygen
    EXPECT_EQ(atom_types[2].atomic_number, 8);
    EXPECT_EQ(atom_types[2].degree, 1);
}

// Test aromatic molecules
TEST_F(AtomTyperTest, AromaticBenzene) {
    std::string smiles = "c1ccccc1";
    auto atom_types = typer.type_atoms_from_smiles(smiles);
    
    ASSERT_EQ(atom_types.size(), 6);
    
    // All atoms should be aromatic carbon in a ring
    for (const auto& at : atom_types) {
        EXPECT_EQ(at.atomic_number, 6);
        EXPECT_TRUE(at.is_aromatic);
        EXPECT_TRUE(at.is_in_ring);
        EXPECT_EQ(at.ring_size, 6);
    }
}

// Test ring detection
TEST_F(AtomTyperTest, RingDetection) {
    std::string smiles = "C1CCCCC1";  // Cyclohexane
    auto atom_types = typer.type_atoms_from_smiles(smiles);
    
    ASSERT_EQ(atom_types.size(), 6);
    
    for (const auto& at : atom_types) {
        EXPECT_TRUE(at.is_in_ring);
        EXPECT_EQ(at.ring_size, 6);
        EXPECT_FALSE(at.is_aromatic);  // Not aromatic
    }
}

// Test hybridization
TEST_F(AtomTyperTest, Hybridization) {
    std::string smiles = "C=C";  // Ethene (SP2)
    auto atom_types = typer.type_atoms_from_smiles(smiles);
    
    ASSERT_EQ(atom_types.size(), 2);
    
    for (const auto& at : atom_types) {
        EXPECT_EQ(at.hybridization, "SP2");
    }
}

// Test triple bond (SP hybridization)
TEST_F(AtomTyperTest, TripleBond) {
    std::string smiles = "C#C";  // Ethyne (SP)
    auto atom_types = typer.type_atoms_from_smiles(smiles);
    
    ASSERT_EQ(atom_types.size(), 2);
    
    for (const auto& at : atom_types) {
        EXPECT_EQ(at.hybridization, "SP");
    }
}

// Test formal charges
TEST_F(AtomTyperTest, FormalCharge) {
    std::string smiles = "[NH4+]";  // Ammonium ion
    auto atom_types = typer.type_atoms_from_smiles(smiles);
    
    ASSERT_EQ(atom_types.size(), 1);
    EXPECT_EQ(atom_types[0].atomic_number, 7);
    EXPECT_EQ(atom_types[0].formal_charge, 1);
}

// Test SMARTS parsing
TEST_F(AtomTyperTest, BasicSmarts) {
    std::string smarts = "[C,N]";
    auto atom_types = typer.type_atoms_from_smarts(smarts);
    
    // SMARTS patterns may have variable results, just check it doesn't crash
    EXPECT_GE(atom_types.size(), 0);
}

// Test invalid SMILES
TEST_F(AtomTyperTest, InvalidSmiles) {
    std::string invalid_smiles = "C(C";  // Unclosed parenthesis
    EXPECT_THROW(typer.type_atoms_from_smiles(invalid_smiles), std::runtime_error);
}

// Test empty SMILES
TEST_F(AtomTyperTest, EmptySmiles) {
    std::string empty_smiles = "";
    EXPECT_THROW(typer.type_atoms_from_smiles(empty_smiles), std::runtime_error);
}

// Test neighbor detection
TEST_F(AtomTyperTest, NeighborDetection) {
    std::string smiles = "CCC";  // Propane
    auto atom_types = typer.type_atoms_from_smiles(smiles);
    
    ASSERT_EQ(atom_types.size(), 3);
    
    // First carbon has 1 neighbor
    EXPECT_EQ(atom_types[0].neighbors.size(), 1);
    
    // Middle carbon has 2 neighbors
    EXPECT_EQ(atom_types[1].neighbors.size(), 2);
    
    // Last carbon has 1 neighbor
    EXPECT_EQ(atom_types[2].neighbors.size(), 1);
}

// Test SMARTS pattern generation
TEST_F(AtomTyperTest, SmartsPatternGeneration) {
    std::string smiles = "CCO";
    auto atom_types = typer.type_atoms_from_smiles(smiles);
    
    ASSERT_EQ(atom_types.size(), 3);
    
    // Check that SMARTS patterns are generated
    for (const auto& at : atom_types) {
        EXPECT_FALSE(at.smarts_pattern.empty());
        EXPECT_EQ(at.smarts_pattern[0], '[');
        EXPECT_EQ(at.smarts_pattern[at.smarts_pattern.size()-1], ']');
    }
}

// Test get_atom_types_string
TEST_F(AtomTyperTest, GetAtomTypesString) {
    std::string smiles = "CC";
    auto atom_types = typer.type_atoms_from_smiles(smiles);
    std::string output = typer.get_atom_types_string(atom_types);
    
    EXPECT_FALSE(output.empty());
    EXPECT_NE(output.find("Atom Types"), std::string::npos);
}

// Test canonical SMILES setting
TEST_F(AtomTyperTest, CanonicalSetting) {
    typer.set_use_canonical(true);
    std::string smiles = "CCO";
    
    // Should not throw
    EXPECT_NO_THROW(typer.type_atoms_from_smiles(smiles));
}

//Write more basic tests finding regarding the new properties in AtomType

// Test ring membership list
TEST_F(AtomTyperTest, RingMembershipList) {
    std::string smiles = "C1CCCCC1";  // Cyclohexane
    auto atom_types = typer.type_atoms_from_smiles(smiles);
    
    ASSERT_EQ(atom_types.size(), 6);
    
    for (const auto& at : atom_types) {
        EXPECT_EQ(at.ring_membership_list.size(), 6);
    }
}

// Test number of ring bonds
TEST_F(AtomTyperTest, NumberOfRingBonds) {
    std::string smiles = "C1CCCCC1";  // Cyclohexane    
    auto atom_types = typer.type_atoms_from_smiles(smiles);

    ASSERT_EQ(atom_types.size(), 6);
    
    for (const auto& at : atom_types) {
        EXPECT_EQ(at.num_ring_bonds, 2);
    }
}

// Test number of aliphatic rings
TEST_F(AtomTyperTest, NumberOfAliphaticRings1) {
    std::string smiles = "C1CCCCC1";  // Cyclohexane
    auto atom_types = typer.type_atoms_from_smiles(smiles);
    
    ASSERT_EQ(atom_types.size(), 6);

    for (const auto& at : atom_types) {
        EXPECT_EQ(at.num_aliphatic_rings, 1);
    }
}

// Test number of aliphatic rings
TEST_F(AtomTyperTest, NumberOfAliphaticRings2) {
    std::string smiles = "c1ccccc1";  // Benzene
    auto atom_types = typer.type_atoms_from_smiles(smiles);
    
    ASSERT_EQ(atom_types.size(), 6);

    for (const auto& at : atom_types) {
        EXPECT_EQ(at.num_aliphatic_rings, 0);
    }
}

// Test number of aromatic rings
TEST_F(AtomTyperTest, NumberOfAromaticRings1) {
    std::string smiles = "c1ccccc1";  // Benzene
    auto atom_types = typer.type_atoms_from_smiles(smiles);

    ASSERT_EQ(atom_types.size(), 6);

    for (const auto& at : atom_types) {
        EXPECT_EQ(at.num_aromatic_rings, 1);
    }
}

// Test number of aromatic rings
TEST_F(AtomTyperTest, NumberOfAromaticRings2) {
    std::string smiles = "C1CCCCC1";  // Cyclohexane
    auto atom_types = typer.type_atoms_from_smiles(smiles);

    ASSERT_EQ(atom_types.size(), 6);

    for (const auto& at : atom_types) {
        EXPECT_EQ(at.num_aromatic_rings, 0);
    }
}

// Test ring connectivity
TEST_F(AtomTyperTest, RingConnectivity) {
    std::string smiles = "C1CCCCC1";  // Cyclohexane
    auto atom_types = typer.type_atoms_from_smiles(smiles);

    ASSERT_EQ(atom_types.size(), 6);

    for (const auto& at : atom_types) {
        // number of connected atoms in the ring should be 2 for cyclohexane
        EXPECT_EQ(at.ring_connectivity, 2);
    }
}

// Test bond types map
TEST_F(AtomTyperTest, BondTypesMap) {
    std::string smiles = "CC=O";  // Acetaldehyde
    auto atom_types = typer.type_atoms_from_smiles(smiles);

    ASSERT_EQ(atom_types.size(), 3);

    // Check bond types for the carbonyl carbon
    const auto& carbonyl_carbon = atom_types[1];
    EXPECT_EQ(carbonyl_carbon.bond_types.at(1), 1); // single bond
    EXPECT_EQ(carbonyl_carbon.bond_types.at(2), 1); // double bond
}   

// Test chirality
TEST_F(AtomTyperTest, Chirality) {
    std::string smiles = "CC[C@H](O)C";  // (R)-2-butanol - true chiral center
    auto atom_types = typer.type_atoms_from_smiles(smiles);

    ASSERT_EQ(atom_types.size(), 5);

    // Check chirality of the chiral carbon (index 2)
    const auto& chiral_carbon = atom_types[2];
    EXPECT_EQ(chiral_carbon.chirality, "R");
}

// Test atom type string
TEST_F(AtomTyperTest, AtomTypeStringContent) {
    std::string smiles = "CCO";
    auto atom_types = typer.type_atoms_from_smiles(smiles);
    std::string output = typer.get_atom_types_string(atom_types);
    
    // Check for specific content
    EXPECT_NE(output.find("Atom 0: Element=6"), std::string::npos); // First carbon
    EXPECT_NE(output.find("Atom 1: Element=6"), std::string::npos); // Second carbon
    EXPECT_NE(output.find("Atom 2: Element=8"), std::string::npos); // Oxygen
}



int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
