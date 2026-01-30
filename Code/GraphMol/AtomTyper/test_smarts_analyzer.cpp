#include "smarts_analyzer.hpp"
#include <gtest/gtest.h>
#include <stdexcept>
#include <algorithm>
#include <set>

class SmartsAnalyzerTest : public ::testing::Test {
protected:
    atom_typer::SmartsAnalyzer analyzer;
};


// calculate_dof Tests

// Test basic SMARTS with no OR queries (DOF = 1)
TEST_F(SmartsAnalyzerTest, CalculateDOF_Simple) {
    EXPECT_EQ(analyzer.calculate_dof("[C]"), 1);
    EXPECT_EQ(analyzer.calculate_dof("[#6]"), 1);
    EXPECT_EQ(analyzer.calculate_dof("[CH3]"), 1);
}

// Test SMARTS with single OR query
TEST_F(SmartsAnalyzerTest, CalculateDOF_SingleOR) {
    // [C,N] = carbon OR nitrogen = 2 possibilities
    EXPECT_EQ(analyzer.calculate_dof("[C,N]"), 2);
    
    // [C,N,O] = 3 possibilities
    EXPECT_EQ(analyzer.calculate_dof("[C,N,O]"), 3);
    
    // [#6,#7,#8,#15] = 4 elements
    EXPECT_EQ(analyzer.calculate_dof("[#6,#7,#8,#15]"), 4);
}

// Test SMARTS with multiple independent OR queries (multiply DOFs)
TEST_F(SmartsAnalyzerTest, CalculateDOF_MultipleOR) {
    // [C,N][C,O] = 2 * 2 = 4
    EXPECT_EQ(analyzer.calculate_dof("[C,N][C,O]"), 4);
    
    // [C,N,O][H1,H2] = 3 * 2 = 6
    EXPECT_EQ(analyzer.calculate_dof("[C,N,O][H1,H2]"), 6);
    
    // [C,N][O,S][H1,H2,H3] = 2 * 2 * 3 = 12
    EXPECT_EQ(analyzer.calculate_dof("[C,N][O,S][H1,H2,H3]"), 12);
}

// Test SMARTS with AND queries (constraints reduce possibilities)
TEST_F(SmartsAnalyzerTest, CalculateDOF_AND) {
    // [C&H3] = carbon AND 3 hydrogens = 1 (specific constraint)
    EXPECT_EQ(analyzer.calculate_dof("[C&H3]"), 1);
    
    // [C&!H0] = carbon AND not zero hydrogens
    EXPECT_GE(analyzer.calculate_dof("[C&!H0]"), 1);
}

// Test SMARTS with NOT queries
TEST_F(SmartsAnalyzerTest, CalculateDOF_NOT) {
    // [!C] = not carbon (but still finite due to exclusion)
    EXPECT_GE(analyzer.calculate_dof("[!C]"), 1);
    
    // [!#6] = not atomic number 6
    EXPECT_GE(analyzer.calculate_dof("[!#6]"), 1);
}

// Test complex nested queries
TEST_F(SmartsAnalyzerTest, CalculateDOF_Complex) {
    // Aromatic carbon or nitrogen: [c,n]
    EXPECT_EQ(analyzer.calculate_dof("[c,n]"), 2);
    
    // Carbon with 2 or 3 hydrogens: [CH2,CH3]
    EXPECT_EQ(analyzer.calculate_dof("[CH2,CH3]"), 2);
}

// Test bond queries
TEST_F(SmartsAnalyzerTest, CalculateDOF_Bonds) {
    // Single bonds don't add DOF
    EXPECT_EQ(analyzer.calculate_dof("C-C"), 1);
    
    // Double bonds don't add DOF
    EXPECT_EQ(analyzer.calculate_dof("C=C"), 1);
    
    // OR in bonds: C-,=C would be 2 if bond types counted
    EXPECT_GE(analyzer.calculate_dof("C-C"), 1);
}

// Test empty or invalid SMARTS
TEST_F(SmartsAnalyzerTest, CalculateDOF_Invalid) {
    // Empty string should return 0 or throw
    EXPECT_EQ(analyzer.calculate_dof(""), 0);
}

// Test ring queries
TEST_F(SmartsAnalyzerTest, CalculateDOF_Rings) {
    // Ring membership doesn't change DOF
    EXPECT_EQ(analyzer.calculate_dof("[C;R]"), 1);
    
    // Ring size specification
    EXPECT_EQ(analyzer.calculate_dof("[C;r6]"), 1);
}


// enumerate_variants Tests


// Test simple SMARTS with no variations
TEST_F(SmartsAnalyzerTest, EnumerateVariants_Simple) {
    auto variants = analyzer.enumerate_variants("[C]", 10);
    ASSERT_EQ(variants.size(), 1);
    EXPECT_EQ(variants[0], "[C]");
}

// Test SMARTS with single OR query
TEST_F(SmartsAnalyzerTest, EnumerateVariants_SingleOR) {
    auto variants = analyzer.enumerate_variants("[C,N]", 10);
    ASSERT_EQ(variants.size(), 2);
    
    // Should contain both [C] and [N]
    std::set<std::string> variant_set(variants.begin(), variants.end());
    EXPECT_TRUE(variant_set.count("[C]") > 0 || variant_set.count("[#6]") > 0);
    EXPECT_TRUE(variant_set.count("[N]") > 0 || variant_set.count("[#7]") > 0);
}

// Test SMARTS with three-way OR
TEST_F(SmartsAnalyzerTest, EnumerateVariants_ThreeWayOR) {
    auto variants = analyzer.enumerate_variants("[C,N,O]", 10);
    ASSERT_EQ(variants.size(), 3);
    
    // Should be unique
    std::set<std::string> unique_variants(variants.begin(), variants.end());
    EXPECT_EQ(unique_variants.size(), 3);
}

// Test SMARTS with multiple OR queries (Cartesian product)
TEST_F(SmartsAnalyzerTest, EnumerateVariants_MultipleOR) {
    // [C,N][O,S] should give 4 variants: CO, CS, NO, NS
    auto variants = analyzer.enumerate_variants("[C,N][O,S]", 10);
    EXPECT_EQ(variants.size(), 4);
    
    // All should be unique
    std::set<std::string> unique_variants(variants.begin(), variants.end());
    EXPECT_EQ(unique_variants.size(), 4);
}

// Test max limit on variants
TEST_F(SmartsAnalyzerTest, EnumerateVariants_MaxLimit) {
    // [C,N,O][H1,H2,H3] = 9 possible, but limit to 5
    auto variants = analyzer.enumerate_variants("[C,N,O][H1,H2,H3]", 5);
    EXPECT_LE(variants.size(), 5);
    
    // Should still be unique
    std::set<std::string> unique_variants(variants.begin(), variants.end());
    EXPECT_EQ(unique_variants.size(), variants.size());
}

// Test that max=0 returns empty
TEST_F(SmartsAnalyzerTest, EnumerateVariants_MaxZero) {
    auto variants = analyzer.enumerate_variants("[C,N]", 0);
    EXPECT_EQ(variants.size(), 0);
}

// Test deduplication
TEST_F(SmartsAnalyzerTest, EnumerateVariants_Deduplication) {
    // Any enumeration should not have duplicates
    auto variants = analyzer.enumerate_variants("[C,N][O,S]", 100);
    std::set<std::string> unique_variants(variants.begin(), variants.end());
    EXPECT_EQ(unique_variants.size(), variants.size());
}

// Test with AND queries
TEST_F(SmartsAnalyzerTest, EnumerateVariants_AND) {
    auto variants = analyzer.enumerate_variants("[C&H3]", 10);
    ASSERT_GE(variants.size(), 1);
    
    // Should generate valid SMARTS
    for (const auto& v : variants) {
        EXPECT_FALSE(v.empty());
    }
}

// Test with NOT queries
TEST_F(SmartsAnalyzerTest, EnumerateVariants_NOT) {
    auto variants = analyzer.enumerate_variants("[!H0]", 10);
    
    // Should generate some variants
    EXPECT_GE(variants.size(), 1);
    
    // Should all be valid non-empty strings
    for (const auto& v : variants) {
        EXPECT_FALSE(v.empty());
    }
}

// Test aromatic queries
TEST_F(SmartsAnalyzerTest, EnumerateVariants_Aromatic) {
    auto variants = analyzer.enumerate_variants("[c,n]", 10);
    EXPECT_EQ(variants.size(), 2);
}

// Test with bonds
TEST_F(SmartsAnalyzerTest, EnumerateVariants_WithBonds) {
    auto variants = analyzer.enumerate_variants("[C,N]-[O,S]", 10);
    EXPECT_EQ(variants.size(), 4);
    
    // All should contain bonds
    for (const auto& v : variants) {
        EXPECT_TRUE(v.find("-") != std::string::npos || 
                    v.size() > 3); // May have implicit bonds
    }
}

// Test that variants are all valid SMARTS
TEST_F(SmartsAnalyzerTest, EnumerateVariants_ValidOutput) {
    auto variants = analyzer.enumerate_variants("[C,N,O]", 10);
    
    for (const auto& v : variants) {
        // Should not be empty
        EXPECT_FALSE(v.empty());
        
        // Should contain brackets for atom primitives
        EXPECT_TRUE(v.find("[") != std::string::npos);
        EXPECT_TRUE(v.find("]") != std::string::npos);
    }
}

// Test edge case: very large DOF with small max
TEST_F(SmartsAnalyzerTest, EnumerateVariants_LargeDOFSmallMax) {
    // [C,N,O,S,P][H0,H1,H2,H3][+0,+1,-1] = 5*4*3 = 60 possibilities
    auto variants = analyzer.enumerate_variants("[C,N,O,S,P][H0,H1,H2,H3][+0,+1,-1]", 10);
    
    // Should respect max limit
    EXPECT_LE(variants.size(), 10);
    
    // Should be unique
    std::set<std::string> unique_variants(variants.begin(), variants.end());
    EXPECT_EQ(unique_variants.size(), variants.size());
}

// Test charge queries
TEST_F(SmartsAnalyzerTest, EnumerateVariants_Charges) {
    auto variants = analyzer.enumerate_variants("[C+1,C+0]", 10);
    EXPECT_EQ(variants.size(), 2);
}

// Test hydrogen count queries
TEST_F(SmartsAnalyzerTest, EnumerateVariants_HydrogenCounts) {
    auto variants = analyzer.enumerate_variants("[CH0,CH1,CH2]", 10);
    EXPECT_EQ(variants.size(), 3);
}

// Test empty SMARTS
TEST_F(SmartsAnalyzerTest, EnumerateVariants_Empty) {
    auto variants = analyzer.enumerate_variants("", 10);
    EXPECT_EQ(variants.size(), 0);
}

// Test realistic chemistry example
TEST_F(SmartsAnalyzerTest, EnumerateVariants_RealisticChemistry) {
    // Carbonyl carbon: [C;$(C=O)]
    auto variants = analyzer.enumerate_variants("[C,N]=O", 10);
    EXPECT_GE(variants.size(), 1);
    EXPECT_LE(variants.size(), 10);
}

// Test that enumerate is consistent with calculate_dof
TEST_F(SmartsAnalyzerTest, ConsistencyCheck_DOFvsEnumerate) {
    std::string smarts = "[C,N][O,S]";
    int dof = analyzer.calculate_dof(smarts);
    auto variants = analyzer.enumerate_variants(smarts, 100);
    
    // Number of variants should match DOF (if max is high enough)
    EXPECT_EQ(variants.size(), dof);
}

// Test another consistency check
TEST_F(SmartsAnalyzerTest, ConsistencyCheck_ThreeElements) {
    std::string smarts = "[C,N,O]";
    int dof = analyzer.calculate_dof(smarts);
    auto variants = analyzer.enumerate_variants(smarts, 100);
    
    EXPECT_EQ(dof, 3);
    EXPECT_EQ(variants.size(), 3);
}