#include "atom_typer.hpp"
#include <catch2/catch_all.hpp>
#include <stdexcept>

struct AtomTyperFixture {
    atom_typer::AtomTyper typer;
};

// Test basic SMILES parsing
TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: Basic SMILES (ethanol)", "[AtomTyper]") {
    std::string smiles = "CCO";
    auto atom_types = typer.type_atoms_from_smiles(smiles);
    
    REQUIRE(atom_types.size() == 3);
    
    // Check first carbon
    CHECK(atom_types[0].atomic_number == 6);
    CHECK(atom_types[0].degree == 1);
    
    // Check second carbon
    CHECK(atom_types[1].atomic_number == 6);
    CHECK(atom_types[1].degree == 2);
    
    // Check oxygen
    CHECK(atom_types[2].atomic_number == 8);
    CHECK(atom_types[2].degree == 1);
}

// Test aromatic molecules
TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: Aromatic benzene", "[AtomTyper]") {
    std::string smiles = "c1ccccc1";
    auto atom_types = typer.type_atoms_from_smiles(smiles);
    
    REQUIRE(atom_types.size() == 6);
    
    // All atoms should be aromatic carbon in a ring
    for (const auto& at : atom_types) {
        CHECK(at.atomic_number == 6);
        CHECK(at.is_aromatic);
        CHECK(at.is_in_ring);
        CHECK(at.ring_size == 6);
    }
}

// Test ring detection
TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: Ring detection", "[AtomTyper]") {
    std::string smiles = "C1CCCCC1";  // Cyclohexane
    auto atom_types = typer.type_atoms_from_smiles(smiles);
    
    REQUIRE(atom_types.size() == 6);
    
    for (const auto& at : atom_types) {
        CHECK(at.is_in_ring);
        CHECK(at.ring_size == 6);
        CHECK_FALSE(at.is_aromatic);  // Not aromatic
    }
}

// Test hybridization
TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: Hybridization (SP2)", "[AtomTyper]") {
    std::string smiles = "C=C";  // Ethene (SP2)
    auto atom_types = typer.type_atoms_from_smiles(smiles);
    
    REQUIRE(atom_types.size() == 2);
    
    for (const auto& at : atom_types) {
        CHECK(at.hybridization == "SP2");
    }
}

// Test triple bond (SP hybridization)
TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: Hybridization (SP)", "[AtomTyper]") {
    std::string smiles = "C#C";  // Ethyne (SP)
    auto atom_types = typer.type_atoms_from_smiles(smiles);
    
    REQUIRE(atom_types.size() == 2);
    
    for (const auto& at : atom_types) {
        CHECK(at.hybridization == "SP");
    }
}

// Test formal charges
TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: Formal charge", "[AtomTyper]") {
    std::string smiles = "[NH4+]";  // Ammonium ion
    auto atom_types = typer.type_atoms_from_smiles(smiles);
    
    REQUIRE(atom_types.size() == 1);
    CHECK(atom_types[0].atomic_number == 7);
    CHECK(atom_types[0].formal_charge == 1);
}

// Test SMARTS parsing
TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: Basic SMARTS", "[AtomTyper]") {
    std::string smarts = "[C,N]";
    auto atom_types = typer.type_atoms_from_smarts(smarts);
    
    // SMARTS patterns may have variable results, just check it doesn't crash
    CHECK_NOTHROW(typer.type_atoms_from_smarts(smarts));
}

// Test invalid SMILES
TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: Invalid SMILES throws", "[AtomTyper]") {
    std::string invalid_smiles = "C(C";  // Unclosed parenthesis
    CHECK_THROWS_AS(typer.type_atoms_from_smiles(invalid_smiles), std::runtime_error);
}

// Test empty SMILES
TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: Empty SMILES throws", "[AtomTyper]") {
    std::string empty_smiles = "";
    CHECK_THROWS_AS(typer.type_atoms_from_smiles(empty_smiles), std::runtime_error);
}

// Test neighbor detection
TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: Neighbor detection", "[AtomTyper]") {
    std::string smiles = "CCC";  // Propane
    auto atom_types = typer.type_atoms_from_smiles(smiles);
    
    REQUIRE(atom_types.size() == 3);
    
    // First carbon has 1 neighbor
    CHECK(atom_types[0].neighbors.size() == 1);
    
    // Middle carbon has 2 neighbors
    CHECK(atom_types[1].neighbors.size() == 2);
    
    // Last carbon has 1 neighbor
    CHECK(atom_types[2].neighbors.size() == 1);
}

// Test SMARTS pattern generation
TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: SMARTS pattern generation", "[AtomTyper]") {
    std::string smiles = "CCO";
    auto atom_types = typer.type_atoms_from_smiles(smiles);
    
    REQUIRE(atom_types.size() == 3);
    
    // Check that SMARTS patterns are generated
    for (const auto& at : atom_types) {
        CHECK_FALSE(at.smarts_pattern.empty());
        CHECK(at.smarts_pattern.front() == '[');
        CHECK(at.smarts_pattern.back() == ']');
    }
}

// Test get_atom_types_string
TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: get_atom_types_string", "[AtomTyper]") {
    std::string smiles = "CC";
    auto atom_types = typer.type_atoms_from_smiles(smiles);
    std::string output = typer.get_atom_types_string(atom_types);
    
    CHECK_FALSE(output.empty());
    CHECK(output.find("Atom Types") != std::string::npos);
}

// Test canonical SMILES setting
TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: canonical setting", "[AtomTyper]") {
    typer.set_use_canonical(true);
    std::string smiles = "CCO";
    
    // Should not throw
    CHECK_NOTHROW(typer.type_atoms_from_smiles(smiles));
}

//Write more basic tests finding regarding the new properties in AtomType

// Test ring membership list
TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: ring membership list", "[AtomTyper]") {
    std::string smiles = "C1CCCCC1";  // Cyclohexane
    auto atom_types = typer.type_atoms_from_smiles(smiles);
    
    REQUIRE(atom_types.size() == 6);
    
    for (const auto& at : atom_types) {
        CHECK(at.ring_membership_list.size() == 6);
    }
}

// Test number of ring bonds
TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: number of ring bonds", "[AtomTyper]") {
    std::string smiles = "C1CCCCC1";  // Cyclohexane    
    auto atom_types = typer.type_atoms_from_smiles(smiles);

    REQUIRE(atom_types.size() == 6);
    
    for (const auto& at : atom_types) {
        CHECK(at.num_ring_bonds == 2);
    }
}

// Test number of aliphatic rings
TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: number of aliphatic rings (cyclohexane)", "[AtomTyper]") {
    std::string smiles = "C1CCCCC1";  // Cyclohexane
    auto atom_types = typer.type_atoms_from_smiles(smiles);
    
    REQUIRE(atom_types.size() == 6);

    for (const auto& at : atom_types) {
        CHECK(at.num_aliphatic_rings == 1);
    }
}

// Test number of aliphatic rings
TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: number of aliphatic rings (benzene)", "[AtomTyper]") {
    std::string smiles = "c1ccccc1";  // Benzene
    auto atom_types = typer.type_atoms_from_smiles(smiles);
    
    REQUIRE(atom_types.size() == 6);

    for (const auto& at : atom_types) {
        CHECK(at.num_aliphatic_rings == 0);
    }
}

// Test number of aromatic rings
TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: number of aromatic rings (benzene)", "[AtomTyper]") {
    std::string smiles = "c1ccccc1";  // Benzene
    auto atom_types = typer.type_atoms_from_smiles(smiles);

    REQUIRE(atom_types.size() == 6);

    for (const auto& at : atom_types) {
        CHECK(at.num_aromatic_rings == 1);
    }
}

// Test number of aromatic rings
TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: number of aromatic rings (cyclohexane)", "[AtomTyper]") {
    std::string smiles = "C1CCCCC1";  // Cyclohexane
    auto atom_types = typer.type_atoms_from_smiles(smiles);

    REQUIRE(atom_types.size() == 6);

    for (const auto& at : atom_types) {
        CHECK(at.num_aromatic_rings == 0);
    }
}

// Test ring connectivity
TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: ring connectivity", "[AtomTyper]") {
    std::string smiles = "C1CCCCC1";  // Cyclohexane
    auto atom_types = typer.type_atoms_from_smiles(smiles);

    REQUIRE(atom_types.size() == 6);

    for (const auto& at : atom_types) {
        // number of connected atoms in the ring should be 2 for cyclohexane
        CHECK(at.ring_connectivity == 2);
    }
}

// Test bond types map
TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: bond types map", "[AtomTyper]") {
    std::string smiles = "CC=O";  // Acetaldehyde
    auto atom_types = typer.type_atoms_from_smiles(smiles);

    REQUIRE(atom_types.size() == 3);

    // Check bond types for the carbonyl carbon
    const auto& carbonyl_carbon = atom_types[1];
    CHECK(carbonyl_carbon.bond_types.at(1) == 1); // single bond
    CHECK(carbonyl_carbon.bond_types.at(2) == 1); // double bond
}   

// Test chirality
TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: chirality", "[AtomTyper]") {
    std::string smiles = "CC[C@H](O)C";  // (R)-2-butanol - true chiral center
    auto atom_types = typer.type_atoms_from_smiles(smiles);

    REQUIRE(atom_types.size() == 5);

    // Check chirality of the chiral carbon (index 2)
    const auto& chiral_carbon = atom_types[2];
    CHECK(chiral_carbon.chirality == "R");
}

// Test atom type string
TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: atom type string content", "[AtomTyper]") {
    std::string smiles = "CCO";
    auto atom_types = typer.type_atoms_from_smiles(smiles);
    std::string output = typer.get_atom_types_string(atom_types);
    
    // Check for specific content
    CHECK(output.find("Atom 0: Element=6") != std::string::npos); // First carbon
    CHECK(output.find("Atom 1: Element=6") != std::string::npos); // Second carbon
    CHECK(output.find("Atom 2: Element=8") != std::string::npos); // Oxygen
}
