#include <GraphMol/AtomTyper/atom_typer.hpp>
#include <GraphMol/AtomTyper/smarts_analyzer.hpp>
#include <catch2/catch_all.hpp>
#include <map>
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
    CHECK(atom_types[0].min_bonds == 1);
    
    // Check second carbon
    CHECK(atom_types[1].atomic_number == 6);
    CHECK(atom_types[1].min_bonds == 2);
    
    // Check oxygen
    CHECK(atom_types[2].atomic_number == 8);
    CHECK(atom_types[2].min_bonds == 1);
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

TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: SMARTS explicit H defaults to 0", "[AtomTyper]") {
    std::string smarts = "[C]";
    auto atom_types = typer.type_atoms_from_smarts(smarts);

    REQUIRE(atom_types.size() == 1);
    CHECK(atom_types[0].num_hydrogens == 0);
}

TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: SMARTS #6 supports aromatic and aliphatic", "[AtomTyper]") {
    std::string smarts = "[#6]";
    auto atom_types = typer.type_atoms_from_smarts(smarts);

    REQUIRE(atom_types.size() == 1);
    CHECK(atom_types[0].is_aromatic);
    CHECK(atom_types[0].is_aliphatic);
    CHECK(atom_types[0].smarts_pattern.find("A,a") != std::string::npos);
}

TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: SMARTS explicit H count preserved", "[AtomTyper]") {
    std::string smarts = "[CH3]";
    auto atom_types = typer.type_atoms_from_smarts(smarts);

    REQUIRE(atom_types.size() == 1);
    CHECK(atom_types[0].num_hydrogens == 3);
}

TEST_CASE_METHOD(AtomTyperFixture,
                 "AtomTyper: enumerate_dof_smarts emits expected flags",
                 "[AtomTyper]") {
    const std::string smarts = "[CH2][CX3]=[O]";
    const std::string enumerated = typer.enumerate_dof_smarts(smarts);

    CHECK_FALSE(enumerated.empty());
    CHECK(enumerated.find("D") != std::string::npos);
    CHECK(enumerated.find("H") != std::string::npos);
    const bool has_charge_token =
        enumerated.find("+") != std::string::npos ||
        enumerated.find("-") != std::string::npos;
    CHECK(has_charge_token);
}

TEST_CASE_METHOD(AtomTyperFixture,
                 "AtomTyper: enumerate_dof_smarts preserves branching",
                 "[AtomTyper]") {
    const std::string smarts = "[#6:1](O)[$([c:3][#6:5]O)]";
    const std::string enumerated = typer.enumerate_dof_smarts(smarts);

    CHECK_FALSE(enumerated.empty());
    CHECK(enumerated.find("[") != std::string::npos);
    CHECK(enumerated.find("]") != std::string::npos);
    CHECK(enumerated.find("([") != std::string::npos);
}

TEST_CASE_METHOD(AtomTyperFixture,
                 "AtomTyper: enumerate_dof_smarts maps newly introduced atoms",
                 "[AtomTyper][enumerate_dof_map_new_atoms]") {
    const std::string smarts = "[C:1]=[C:2]";
    const std::string enumerated = typer.enumerate_dof_smarts(smarts, true);

    CHECK_FALSE(enumerated.empty());
    CHECK(enumerated.find(":1") != std::string::npos);
    CHECK(enumerated.find(":2") != std::string::npos);
}

TEST_CASE_METHOD(AtomTyperFixture,
                 "AtomTyper: enumerate_dof_smarts handles complex bond OR variants",
                 "[AtomTyper][enumerate_dof_complex_bonds]") {
    atom_typer::SmartsAnalyzer analyzer;
    const auto variants = analyzer.enumerate_variants("[#6:1]=,-,:[#6:2]", 1000);

    REQUIRE(variants.size() == 3);
    for (const auto& v : variants) {
        INFO("Variant: " << v);
        const auto enumerated = typer.enumerate_dof_smarts(v);
        CHECK_FALSE(enumerated.empty());
    }
}

TEST_CASE_METHOD(AtomTyperFixture,
                 "AtomTyper: aromatic constrained atom is not collapsed",
                 "[AtomTyper]") {
    const std::string smarts = "[#6:1]:[c&H2&+0]";
    CHECK_THROWS_WITH(
        typer.enumerate_dof_smarts(smarts),
        Catch::Matchers::ContainsSubstring(
            "No valid DoF atom alternatives could be generated") &&
            Catch::Matchers::ContainsSubstring("source_smarts='"));
}

TEST_CASE_METHOD(AtomTyperFixture,
                 "AtomTyper: reorder_query_tree_by_embedding orders AtomAnd",
                 "[AtomTyper]") {
    const std::string smarts = "[#6&H1&+0]";
    const std::map<std::string, double> embedding = {
        {"#6", 10.0}, {"H1", 1.0}, {"+0", 5.0}};

    const std::string reordered =
        typer.reorder_query_tree_by_embedding(smarts, embedding);
    CHECK_FALSE(reordered.empty());
    CHECK(reordered.find("H1") != std::string::npos);
    CHECK(reordered.find("+0") != std::string::npos);
    CHECK(reordered.find("#6") != std::string::npos);
}

TEST_CASE_METHOD(AtomTyperFixture,
                 "AtomTyper: reorder_query_tree_by_embedding orders AtomOr",
                 "[AtomTyper]") {
    const std::string smarts = "[#6,#7,#8]";
    const std::map<std::string, double> embedding = {
        {"#6", 10.0}, {"#7", 1.0}, {"#8", 5.0}};

    const std::string reordered =
        typer.reorder_query_tree_by_embedding(smarts, embedding);
    CHECK_FALSE(reordered.empty());
    CHECK(reordered.find("#6") != std::string::npos);
    CHECK(reordered.find("#7") != std::string::npos);
    CHECK(reordered.find("#8") != std::string::npos);
}

TEST_CASE_METHOD(AtomTyperFixture,
                 "AtomTyper: default embedding is usable",
                 "[AtomTyper]") {
    const auto embedding = typer.get_default_query_embedding();
    CHECK_FALSE(embedding.empty());
    CHECK(embedding.count("#6") > 0);
    CHECK(embedding.count("!H0") > 0);

    const std::string reordered =
        typer.reorder_query_tree_by_embedding("[#6&H1&+0]", embedding);
    CHECK_FALSE(reordered.empty());
}

TEST_CASE_METHOD(AtomTyperFixture,
                 "AtomTyper: is_valid_valence_smarts for fully typed patterns",
                 "[AtomTyper]") {
    const std::string valid =
        "[#6&D4&H3&+0&A:1]-[#6&D4&H3&+0&A:2]";
    const std::string invalid = "[#6&D5&H0&+0&A]";

    CHECK(typer.is_valid_valence_smarts(valid));
    CHECK_FALSE(typer.is_valid_valence_smarts(invalid));
}

TEST_CASE_METHOD(AtomTyperFixture,
                 "AtomTyper: is_valid_valence_smarts rejects aromatic bond between aliphatic-only atoms",
                 "[AtomTyper]") {
    const std::string invalid = "[#6&D0&H3&+&A:1]:[#6&D0&H3&+&A:2]";
    CHECK_FALSE(typer.is_valid_valence_smarts(invalid));
}

TEST_CASE_METHOD(AtomTyperFixture,
                 "AtomTyper: is_valid_valence_smarts rejects aliphatic bond between aromatic-only atoms",
                 "[AtomTyper]") {
    const std::string invalid = "[#6&D3&H0&+0&a:1]-[#6&D3&H0&+0&a:2]";
    CHECK_FALSE(typer.is_valid_valence_smarts(invalid));
}

TEST_CASE_METHOD(AtomTyperFixture,
                 "AtomTyper: is_valid_valence_smarts handles negated primitives",
                 "[AtomTyper]") {
    CHECK(typer.is_valid_valence_smarts("[#6;!H1]"));
    CHECK_FALSE(typer.is_valid_valence_smarts("[#6;H1;!H1]"));
    CHECK_FALSE(typer.is_valid_valence_smarts("[#6;+1;!+1]"));
}

TEST_CASE_METHOD(AtomTyperFixture,
                 "AtomTyper: filter_invalid_valence_smarts discards invalid",
                 "[AtomTyper]") {
    const std::vector<std::string> input = {
    "[#6&D4&H3&+0&A:1]-[#6&D4&H3&+0&A:2]",
        "[#6&D0&H3&+&A:1]:[#6&D0&H3&+&A:2]",
    "[#6;H1;!H1]",
        "[#6&D5&H0&+0&A]"};

    const auto filtered = typer.filter_invalid_valence_smarts(input);
    REQUIRE(filtered.size() == 1);
    CHECK(filtered[0] == input[0]);
}

TEST_CASE_METHOD(AtomTyperFixture,
                 "AtomTyper: enumerate_dof_smarts preserves negated atomic-number atoms",
                 "[AtomTyper]") {
    const std::string smarts = "[!#6:1]=[!#6:2]";
    const std::string enumerated = typer.enumerate_dof_smarts(smarts);

    CHECK_FALSE(enumerated.empty());
    CHECK(enumerated.find("!#6") != std::string::npos);
    CHECK(enumerated.find("#0") == std::string::npos);
}

TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: SMARTS bond counts and hybridization", "[AtomTyper]") {
    std::string smarts = "[C]=[C]";
    auto atom_types = typer.type_atoms_from_smarts(smarts);

    REQUIRE(atom_types.size() == 2);
    for (const auto& at : atom_types) {
        CHECK(at.num_single_bonds == 0);
        CHECK(at.num_double_bonds == 1);
        CHECK(at.num_triple_bonds == 0);
        CHECK(at.num_aromatic_bonds == 0);
        CHECK(at.hybridization == "SP2");
    }

    const auto pattern_items = typer.type_pattern_from_smarts(smarts);
    const auto typed_pattern = typer.get_pattern_types_string(pattern_items);
    CHECK(typed_pattern.find("SingleBonds=0") != std::string::npos);
    CHECK(typed_pattern.find("DoubleBonds=1") != std::string::npos);
}

TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: remaining valence uses RDKit valence API", "[AtomTyper]") {
    std::string smiles = "CCO";
    auto atom_types = typer.type_atoms_from_smiles(smiles);

    REQUIRE(atom_types.size() == 3);
    CHECK(atom_types[0].remaining_valence == 3);
    CHECK(atom_types[1].remaining_valence == 2);
    CHECK(atom_types[2].remaining_valence == 1);
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
