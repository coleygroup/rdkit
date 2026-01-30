#include "expression_builder.hpp"
#include <gtest/gtest.h>
#include <stdexcept>
#include <algorithm>
#include <set>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/QueryAtom.h>

class ExpressionBuilderTest : public ::testing::Test {
protected:
    // No member variables needed 
};

TEST_F(ExpressionBuilderTest, SmilesToSmarts_Minimal) {
    std::string smarts = atom_typer::smiles_to_smarts("CCO", atom_typer::Level::MINIMAL);
    EXPECT_EQ(smarts, "[C][C][O]");
}

TEST_F(ExpressionBuilderTest, SmilesToSmarts_Standard) {
    std::string smarts = atom_typer::smiles_to_smarts("CCO", atom_typer::Level::STANDARD);
    EXPECT_EQ(smarts, "[C;D1][C;D2][O;D1]");
}

TEST_F(ExpressionBuilderTest, SmilesToSmarts_Detailed) {
    std::string smarts = atom_typer::smiles_to_smarts("CCO", atom_typer::Level::DETAILED);
    EXPECT_EQ(smarts, "[C;D1;H3][C;D2;H2][O;D1;H1]");
}

// Test aromatic molecule (benzene)
TEST_F(ExpressionBuilderTest, SmilesToSmarts_Aromatic_Minimal) {
    std::string smarts = atom_typer::smiles_to_smarts("c1ccccc1", atom_typer::Level::MINIMAL);
    EXPECT_EQ(smarts, "[C][C][C][C][C][C]");
}

TEST_F(ExpressionBuilderTest, SmilesToSmarts_Aromatic_Standard) {
    std::string smarts = atom_typer::smiles_to_smarts("c1ccccc1", atom_typer::Level::STANDARD);
    // All carbons in benzene have degree 2
    EXPECT_EQ(smarts, "[C;D2][C;D2][C;D2][C;D2][C;D2][C;D2]");
}

// Test branched molecule
TEST_F(ExpressionBuilderTest, SmilesToSmarts_Branched_Minimal) {
    std::string smarts = atom_typer::smiles_to_smarts("CC(C)C", atom_typer::Level::MINIMAL);
    EXPECT_EQ(smarts, "[C][C][C][C]");
}

TEST_F(ExpressionBuilderTest, SmilesToSmarts_Branched_Standard) {
    std::string smarts = atom_typer::smiles_to_smarts("CC(C)C", atom_typer::Level::STANDARD);
    // First C: D1, middle C: D3, branch C: D1, last C: D1
    EXPECT_EQ(smarts, "[C;D1][C;D3][C;D1][C;D1]");
}

// Test charged molecule
TEST_F(ExpressionBuilderTest, SmilesToSmarts_Charged_Detailed) {
    std::string smarts = atom_typer::smiles_to_smarts("[NH4+]", atom_typer::Level::DETAILED);
    // Ammonium ion: N with 4 hydrogens and +1 charge
    EXPECT_EQ(smarts, "[N;D0;H4]");
}

// Test single atom
TEST_F(ExpressionBuilderTest, SmilesToSmarts_SingleAtom_Minimal) {
    std::string smarts = atom_typer::smiles_to_smarts("C", atom_typer::Level::MINIMAL);
    EXPECT_EQ(smarts, "[C]");
}

TEST_F(ExpressionBuilderTest, SmilesToSmarts_SingleAtom_Standard) {
    std::string smarts = atom_typer::smiles_to_smarts("C", atom_typer::Level::STANDARD);
    EXPECT_EQ(smarts, "[C;D0]");
}

TEST_F(ExpressionBuilderTest, SmilesToSmarts_SingleAtom_Detailed) {
    std::string smarts = atom_typer::smiles_to_smarts("C", atom_typer::Level::DETAILED);
    EXPECT_EQ(smarts, "[C;D0;H4]");
}

// Test nitrogen-containing molecule
TEST_F(ExpressionBuilderTest, SmilesToSmarts_Nitrogen_Minimal) {
    std::string smarts = atom_typer::smiles_to_smarts("CCN", atom_typer::Level::MINIMAL);
    EXPECT_EQ(smarts, "[C][C][N]");
}

TEST_F(ExpressionBuilderTest, SmilesToSmarts_Nitrogen_Detailed) {
    std::string smarts = atom_typer::smiles_to_smarts("CCN", atom_typer::Level::DETAILED);
    EXPECT_EQ(smarts, "[C;D1;H3][C;D2;H2][N;D1;H2]");
}

// Test double bond molecule
TEST_F(ExpressionBuilderTest, SmilesToSmarts_DoubleBond_Minimal) {
    std::string smarts = atom_typer::smiles_to_smarts("CC=O", atom_typer::Level::MINIMAL);
    EXPECT_EQ(smarts, "[C][C][O]");
}

TEST_F(ExpressionBuilderTest, SmilesToSmarts_DoubleBond_Standard) {
    std::string smarts = atom_typer::smiles_to_smarts("CC=O", atom_typer::Level::STANDARD);
    EXPECT_EQ(smarts, "[C;D1][C;D2][O;D1]");
}

// Test cyclic molecule (cyclopropane)
TEST_F(ExpressionBuilderTest, SmilesToSmarts_Cyclic_Standard) {
    std::string smarts = atom_typer::smiles_to_smarts("C1CC1", atom_typer::Level::STANDARD);
    // All carbons in cyclopropane have degree 2
    EXPECT_EQ(smarts, "[C;D2][C;D2][C;D2]");
}

// Test error handling - empty SMILES
TEST_F(ExpressionBuilderTest, SmilesToSmarts_EmptyString_ThrowsException) {
    EXPECT_THROW({
        atom_typer::smiles_to_smarts("", atom_typer::Level::MINIMAL);
    }, std::invalid_argument);
}

// Test error handling - invalid SMILES
TEST_F(ExpressionBuilderTest, SmilesToSmarts_InvalidSmiles_ThrowsException) {
    EXPECT_THROW({
        atom_typer::smiles_to_smarts("CCCZ", atom_typer::Level::MINIMAL);
    }, std::invalid_argument);
}

// Test complete level returns non-empty result
TEST_F(ExpressionBuilderTest, SmilesToSmarts_Complete_NotEmpty) {
    std::string smarts = atom_typer::smiles_to_smarts("CCO", atom_typer::Level::COMPLETE);
    EXPECT_FALSE(smarts.empty());
    EXPECT_GT(smarts.length(), 0);
}

// Test sulfur-containing molecule
TEST_F(ExpressionBuilderTest, SmilesToSmarts_Sulfur_Minimal) {
    std::string smarts = atom_typer::smiles_to_smarts("CCS", atom_typer::Level::MINIMAL);
    EXPECT_EQ(smarts, "[C][C][S]");
}

// Test complex molecule with multiple heteroatoms
TEST_F(ExpressionBuilderTest, SmilesToSmarts_MultipleHeteroatoms_Minimal) {
    std::string smarts = atom_typer::smiles_to_smarts("C(N)O", atom_typer::Level::MINIMAL);
    EXPECT_EQ(smarts, "[C][N][O]");
}

TEST_F(ExpressionBuilderTest, SmilesToSmarts_MultipleHeteroatoms_Standard) {
    std::string smarts = atom_typer::smiles_to_smarts("C(N)O", atom_typer::Level::STANDARD);
    // C: D2, N: D1, O: D1
    EXPECT_EQ(smarts, "[C;D2][N;D1][O;D1]");
}

// Test queryMoleculeToSmarts function directly
TEST_F(ExpressionBuilderTest, QueryMoleculeToSmarts_Minimal) {
    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol("CCO"));
    ASSERT_TRUE(mol != nullptr);
    
    std::unique_ptr<RDKit::RWMol> queryMol(new RDKit::RWMol());
    for (const auto atom : mol->atoms()) {
        RDKit::QueryAtom* queryAtom = new RDKit::QueryAtom(atom->getAtomicNum());
        queryMol->addAtom(queryAtom, false, true);
    }
    
    std::string smarts = atom_typer::queryMoleculeToSmarts(mol.get(), queryMol.get(), 
                                                            atom_typer::Level::MINIMAL);
    EXPECT_EQ(smarts, "[C][C][O]");
}

TEST_F(ExpressionBuilderTest, QueryMoleculeToSmarts_Standard) {
    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol("CCO"));
    ASSERT_TRUE(mol != nullptr);
    
    std::unique_ptr<RDKit::RWMol> queryMol(new RDKit::RWMol());
    for (const auto atom : mol->atoms()) {
        RDKit::QueryAtom* queryAtom = new RDKit::QueryAtom(atom->getAtomicNum());
        queryMol->addAtom(queryAtom, false, true);
    }
    
    std::string smarts = atom_typer::queryMoleculeToSmarts(mol.get(), queryMol.get(), 
                                                            atom_typer::Level::STANDARD);
    EXPECT_EQ(smarts, "[C;D1][C;D2][O;D1]");
}

TEST_F(ExpressionBuilderTest, QueryMoleculeToSmarts_Detailed) {
    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol("CCO"));
    ASSERT_TRUE(mol != nullptr);
    
    std::unique_ptr<RDKit::RWMol> queryMol(new RDKit::RWMol());
    for (const auto atom : mol->atoms()) {
        RDKit::QueryAtom* queryAtom = new RDKit::QueryAtom(atom->getAtomicNum());
        queryMol->addAtom(queryAtom, false, true);
    }
    
    std::string smarts = atom_typer::queryMoleculeToSmarts(mol.get(), queryMol.get(), 
                                                            atom_typer::Level::DETAILED);
    EXPECT_EQ(smarts, "[C;D1;H3][C;D2;H2][O;D1;H1]");
}

TEST_F(ExpressionBuilderTest, QueryMoleculeToSmarts_SimpleCarbon) {
    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol("C"));
    ASSERT_TRUE(mol != nullptr);
    
    std::unique_ptr<RDKit::RWMol> queryMol(new RDKit::RWMol());
    for (const auto atom : mol->atoms()) {
        RDKit::QueryAtom* queryAtom = new RDKit::QueryAtom(atom->getAtomicNum());
        queryMol->addAtom(queryAtom, false, true);
    }
    
    std::string smartsMinimal = atom_typer::queryMoleculeToSmarts(mol.get(), queryMol.get(), 
                                                                   atom_typer::Level::MINIMAL);
    EXPECT_EQ(smartsMinimal, "[C]");
    
    std::string smartsStandard = atom_typer::queryMoleculeToSmarts(mol.get(), queryMol.get(), 
                                                                    atom_typer::Level::STANDARD);
    EXPECT_EQ(smartsStandard, "[C;D0]");
    
    std::string smartsDetailed = atom_typer::queryMoleculeToSmarts(mol.get(), queryMol.get(), 
                                                                    atom_typer::Level::DETAILED);
    EXPECT_EQ(smartsDetailed, "[C;D0;H4]");
}