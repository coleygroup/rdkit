#include "expression_builder.hpp"
#include <catch2/catch_all.hpp>
#include <stdexcept>
#include <algorithm>
#include <set>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/QueryAtom.h>

TEST_CASE("ExpressionBuilder: SmilesToSmarts minimal", "[ExpressionBuilder]") {
    std::string smarts = atom_typer::smiles_to_smarts("CCO", atom_typer::Level::MINIMAL);
    CHECK(smarts == "[C][C][O]");
}

TEST_CASE("ExpressionBuilder: SmilesToSmarts standard", "[ExpressionBuilder]") {
    std::string smarts = atom_typer::smiles_to_smarts("CCO", atom_typer::Level::STANDARD);
    CHECK(smarts == "[C;D1][C;D2][O;D1]");
}

TEST_CASE("ExpressionBuilder: SmilesToSmarts detailed", "[ExpressionBuilder]") {
    std::string smarts = atom_typer::smiles_to_smarts("CCO", atom_typer::Level::DETAILED);
    CHECK(smarts == "[C;D1;H3][C;D2;H2][O;D1;H1]");
}

// Test aromatic molecule (benzene)
TEST_CASE("ExpressionBuilder: SmilesToSmarts aromatic minimal", "[ExpressionBuilder]") {
    std::string smarts = atom_typer::smiles_to_smarts("c1ccccc1", atom_typer::Level::MINIMAL);
    CHECK(smarts == "[C][C][C][C][C][C]");
}

TEST_CASE("ExpressionBuilder: SmilesToSmarts aromatic standard", "[ExpressionBuilder]") {
    std::string smarts = atom_typer::smiles_to_smarts("c1ccccc1", atom_typer::Level::STANDARD);
    // All carbons in benzene have degree 2
    CHECK(smarts == "[C;D2][C;D2][C;D2][C;D2][C;D2][C;D2]");
}

// Test branched molecule
TEST_CASE("ExpressionBuilder: SmilesToSmarts branched minimal", "[ExpressionBuilder]") {
    std::string smarts = atom_typer::smiles_to_smarts("CC(C)C", atom_typer::Level::MINIMAL);
    CHECK(smarts == "[C][C][C][C]");
}

TEST_CASE("ExpressionBuilder: SmilesToSmarts branched standard", "[ExpressionBuilder]") {
    std::string smarts = atom_typer::smiles_to_smarts("CC(C)C", atom_typer::Level::STANDARD);
    // First C: D1, middle C: D3, branch C: D1, last C: D1
    CHECK(smarts == "[C;D1][C;D3][C;D1][C;D1]");
}

// Test charged molecule
TEST_CASE("ExpressionBuilder: SmilesToSmarts charged detailed", "[ExpressionBuilder]") {
    std::string smarts = atom_typer::smiles_to_smarts("[NH4+]", atom_typer::Level::DETAILED);
    // Ammonium ion: N with 4 hydrogens and +1 charge
    CHECK(smarts == "[N;D0;H4]");
}

// Test single atom
TEST_CASE("ExpressionBuilder: SmilesToSmarts single atom minimal", "[ExpressionBuilder]") {
    std::string smarts = atom_typer::smiles_to_smarts("C", atom_typer::Level::MINIMAL);
    CHECK(smarts == "[C]");
}

TEST_CASE("ExpressionBuilder: SmilesToSmarts single atom standard", "[ExpressionBuilder]") {
    std::string smarts = atom_typer::smiles_to_smarts("C", atom_typer::Level::STANDARD);
    CHECK(smarts == "[C;D0]");
}

TEST_CASE("ExpressionBuilder: SmilesToSmarts single atom detailed", "[ExpressionBuilder]") {
    std::string smarts = atom_typer::smiles_to_smarts("C", atom_typer::Level::DETAILED);
    CHECK(smarts == "[C;D0;H4]");
}

// Test nitrogen-containing molecule
TEST_CASE("ExpressionBuilder: SmilesToSmarts nitrogen minimal", "[ExpressionBuilder]") {
    std::string smarts = atom_typer::smiles_to_smarts("CCN", atom_typer::Level::MINIMAL);
    CHECK(smarts == "[C][C][N]");
}

TEST_CASE("ExpressionBuilder: SmilesToSmarts nitrogen detailed", "[ExpressionBuilder]") {
    std::string smarts = atom_typer::smiles_to_smarts("CCN", atom_typer::Level::DETAILED);
    CHECK(smarts == "[C;D1;H3][C;D2;H2][N;D1;H2]");
}

// Test double bond molecule
TEST_CASE("ExpressionBuilder: SmilesToSmarts double bond minimal", "[ExpressionBuilder]") {
    std::string smarts = atom_typer::smiles_to_smarts("CC=O", atom_typer::Level::MINIMAL);
    CHECK(smarts == "[C][C][O]");
}

TEST_CASE("ExpressionBuilder: SmilesToSmarts double bond standard", "[ExpressionBuilder]") {
    std::string smarts = atom_typer::smiles_to_smarts("CC=O", atom_typer::Level::STANDARD);
    CHECK(smarts == "[C;D1][C;D2][O;D1]");
}

// Test cyclic molecule (cyclopropane)
TEST_CASE("ExpressionBuilder: SmilesToSmarts cyclic standard", "[ExpressionBuilder]") {
    std::string smarts = atom_typer::smiles_to_smarts("C1CC1", atom_typer::Level::STANDARD);
    // All carbons in cyclopropane have degree 2
    CHECK(smarts == "[C;D2][C;D2][C;D2]");
}

// Test error handling - empty SMILES
TEST_CASE("ExpressionBuilder: SmilesToSmarts empty string throws", "[ExpressionBuilder]") {
    CHECK_THROWS_AS(atom_typer::smiles_to_smarts("", atom_typer::Level::MINIMAL), std::invalid_argument);
}

// Test error handling - invalid SMILES
TEST_CASE("ExpressionBuilder: SmilesToSmarts invalid SMILES throws", "[ExpressionBuilder]") {
    CHECK_THROWS_AS(atom_typer::smiles_to_smarts("CCCZ", atom_typer::Level::MINIMAL), std::invalid_argument);
}

// Test complete level returns non-empty result
TEST_CASE("ExpressionBuilder: SmilesToSmarts complete not empty", "[ExpressionBuilder]") {
    std::string smarts = atom_typer::smiles_to_smarts("CCO", atom_typer::Level::COMPLETE);
    CHECK_FALSE(smarts.empty());
    CHECK(smarts.length() > 0);
}

// Test sulfur-containing molecule
TEST_CASE("ExpressionBuilder: SmilesToSmarts sulfur minimal", "[ExpressionBuilder]") {
    std::string smarts = atom_typer::smiles_to_smarts("CCS", atom_typer::Level::MINIMAL);
    CHECK(smarts == "[C][C][S]");
}

// Test complex molecule with multiple heteroatoms
TEST_CASE("ExpressionBuilder: SmilesToSmarts multiple heteroatoms minimal", "[ExpressionBuilder]") {
    std::string smarts = atom_typer::smiles_to_smarts("C(N)O", atom_typer::Level::MINIMAL);
    CHECK(smarts == "[C][N][O]");
}

TEST_CASE("ExpressionBuilder: SmilesToSmarts multiple heteroatoms standard", "[ExpressionBuilder]") {
    std::string smarts = atom_typer::smiles_to_smarts("C(N)O", atom_typer::Level::STANDARD);
    // C: D2, N: D1, O: D1
    CHECK(smarts == "[C;D2][N;D1][O;D1]");
}

// Test queryMoleculeToSmarts function directly
TEST_CASE("ExpressionBuilder: queryMoleculeToSmarts minimal", "[ExpressionBuilder]") {
    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol("CCO"));
    REQUIRE(mol != nullptr);
    
    std::unique_ptr<RDKit::RWMol> queryMol(new RDKit::RWMol());
    for (const auto atom : mol->atoms()) {
        RDKit::QueryAtom* queryAtom = new RDKit::QueryAtom(atom->getAtomicNum());
        queryMol->addAtom(queryAtom, false, true);
    }
    
    std::string smarts = atom_typer::queryMoleculeToSmarts(mol.get(), queryMol.get(), 
                                                            atom_typer::Level::MINIMAL);
    CHECK(smarts == "[C][C][O]");
}

TEST_CASE("ExpressionBuilder: queryMoleculeToSmarts standard", "[ExpressionBuilder]") {
    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol("CCO"));
    REQUIRE(mol != nullptr);
    
    std::unique_ptr<RDKit::RWMol> queryMol(new RDKit::RWMol());
    for (const auto atom : mol->atoms()) {
        RDKit::QueryAtom* queryAtom = new RDKit::QueryAtom(atom->getAtomicNum());
        queryMol->addAtom(queryAtom, false, true);
    }
    
    std::string smarts = atom_typer::queryMoleculeToSmarts(mol.get(), queryMol.get(), 
                                                            atom_typer::Level::STANDARD);
    CHECK(smarts == "[C;D1][C;D2][O;D1]");
}

TEST_CASE("ExpressionBuilder: queryMoleculeToSmarts detailed", "[ExpressionBuilder]") {
    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol("CCO"));
    REQUIRE(mol != nullptr);
    
    std::unique_ptr<RDKit::RWMol> queryMol(new RDKit::RWMol());
    for (const auto atom : mol->atoms()) {
        RDKit::QueryAtom* queryAtom = new RDKit::QueryAtom(atom->getAtomicNum());
        queryMol->addAtom(queryAtom, false, true);
    }
    
    std::string smarts = atom_typer::queryMoleculeToSmarts(mol.get(), queryMol.get(), 
                                                            atom_typer::Level::DETAILED);
    CHECK(smarts == "[C;D1;H3][C;D2;H2][O;D1;H1]");
}

TEST_CASE("ExpressionBuilder: queryMoleculeToSmarts simple carbon", "[ExpressionBuilder]") {
    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol("C"));
    REQUIRE(mol != nullptr);
    
    std::unique_ptr<RDKit::RWMol> queryMol(new RDKit::RWMol());
    for (const auto atom : mol->atoms()) {
        RDKit::QueryAtom* queryAtom = new RDKit::QueryAtom(atom->getAtomicNum());
        queryMol->addAtom(queryAtom, false, true);
    }
    
    std::string smartsMinimal = atom_typer::queryMoleculeToSmarts(mol.get(), queryMol.get(), 
                                                                   atom_typer::Level::MINIMAL);
    CHECK(smartsMinimal == "[C]");
    
    std::string smartsStandard = atom_typer::queryMoleculeToSmarts(mol.get(), queryMol.get(), 
                                                                    atom_typer::Level::STANDARD);
    CHECK(smartsStandard == "[C;D0]");
    
    std::string smartsDetailed = atom_typer::queryMoleculeToSmarts(mol.get(), queryMol.get(), 
                                                                    atom_typer::Level::DETAILED);
    CHECK(smartsDetailed == "[C;D0;H4]");
}