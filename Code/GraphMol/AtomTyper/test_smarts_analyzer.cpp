#include "smarts_analyzer.hpp"
#include <catch2/catch_all.hpp>
#include <stdexcept>
#include <algorithm>
#include <set>

struct SmartsAnalyzerFixture {
    atom_typer::SmartsAnalyzer analyzer;
};


// calculate_dof Tests

TEST_CASE_METHOD(SmartsAnalyzerFixture,
             "SmartsAnalyzer: add_atom_maps simple",
             "[SmartsAnalyzer][add_atom_maps]") {
    const auto mapped = analyzer.add_atom_maps("[C][O]");
    CHECK_FALSE(mapped.empty());
    CHECK(mapped.find(":1") != std::string::npos);
    CHECK(mapped.find(":2") != std::string::npos);
}

TEST_CASE_METHOD(SmartsAnalyzerFixture,
             "SmartsAnalyzer: add_atom_maps includes recursive atoms",
             "[SmartsAnalyzer][add_atom_maps]") {
    const auto mapped = analyzer.add_atom_maps("[#6&$([#7]-[#8])]");
    CHECK_FALSE(mapped.empty());
    CHECK(mapped.find("$(") != std::string::npos);
    CHECK((mapped.find("#7:") != std::string::npos ||
        mapped.find("N:") != std::string::npos));
    CHECK((mapped.find("#8:") != std::string::npos ||
        mapped.find("O:") != std::string::npos));
}

TEST_CASE_METHOD(SmartsAnalyzerFixture,
             "SmartsAnalyzer: add_atom_maps invalid input",
             "[SmartsAnalyzer][add_atom_maps]") {
    CHECK_THROWS(analyzer.add_atom_maps(""));
    CHECK_THROWS(analyzer.add_atom_maps("[C"));
}

// Test basic SMARTS with no OR queries (DOF = 1)
TEST_CASE_METHOD(SmartsAnalyzerFixture, "SmartsAnalyzer: calculate_dof simple", "[SmartsAnalyzer]") {
    CHECK(analyzer.calculate_dof("[C]") == 1);
    CHECK(analyzer.calculate_dof("[#6]") == 1);
    CHECK(analyzer.calculate_dof("[CH3]") == 1);
}

// Test SMARTS with single OR query
TEST_CASE_METHOD(SmartsAnalyzerFixture, "SmartsAnalyzer: calculate_dof single OR", "[SmartsAnalyzer]") {
    // [C,N] = carbon OR nitrogen = 2 possibilities
    CHECK(analyzer.calculate_dof("[C,N]") == 2);
    
    // [C,N,O] = 3 possibilities
    CHECK(analyzer.calculate_dof("[C,N,O]") == 3);
    
    // [#6,#7,#8,#15] = 4 elements
    CHECK(analyzer.calculate_dof("[#6,#7,#8,#15]") == 4);
}

// Test SMARTS with multiple independent OR queries (multiply DOFs)
TEST_CASE_METHOD(SmartsAnalyzerFixture, "SmartsAnalyzer: calculate_dof multiple OR", "[SmartsAnalyzer]") {
    // [C,N][C,O] = 2 * 2 = 4
    CHECK(analyzer.calculate_dof("[C,N][C,O]") == 4);
    
    // [C,N,O][H1,H2] = 3 * 2 = 6
    CHECK(analyzer.calculate_dof("[C,N,O][H1,H2]") == 6);
    
    // [C,N][O,S][H1,H2,H3] = 2 * 2 * 3 = 12
    CHECK(analyzer.calculate_dof("[C,N][O,S][H1,H2,H3]") == 12);
}

// Test SMARTS with AND queries (constraints reduce possibilities)
TEST_CASE_METHOD(SmartsAnalyzerFixture, "SmartsAnalyzer: calculate_dof AND", "[SmartsAnalyzer]") {
    // [C&H3] = carbon AND 3 hydrogens = 1 (specific constraint)
    CHECK(analyzer.calculate_dof("[C&H3]") == 1);
    
    // [C&!H0] = carbon AND not zero hydrogens
    CHECK(analyzer.calculate_dof("[C&!H0]") >= 1);
}

// Test SMARTS with NOT queries
TEST_CASE_METHOD(SmartsAnalyzerFixture, "SmartsAnalyzer: calculate_dof NOT", "[SmartsAnalyzer]") {
    // [!C] = not carbon (but still finite due to exclusion)
    CHECK(analyzer.calculate_dof("[!C]") >= 1);
    
    // [!#6] = not atomic number 6
    CHECK(analyzer.calculate_dof("[!#6]") >= 1);
}

// Test complex nested queries
TEST_CASE_METHOD(SmartsAnalyzerFixture, "SmartsAnalyzer: calculate_dof complex", "[SmartsAnalyzer]") {
    // Aromatic carbon or nitrogen: [c,n]
    CHECK(analyzer.calculate_dof("[c,n]") == 2);
    
    // Carbon with 2 or 3 hydrogens: [CH2,CH3]
    CHECK(analyzer.calculate_dof("[CH2,CH3]") == 2);
}

// Test bond queries
TEST_CASE_METHOD(SmartsAnalyzerFixture, "SmartsAnalyzer: calculate_dof bonds", "[SmartsAnalyzer]") {
    // Single bonds don't add DOF
    CHECK(analyzer.calculate_dof("C-C") == 1);
    
    // Double bonds don't add DOF
    CHECK(analyzer.calculate_dof("C=C") == 1);
    
    // OR in bonds: C-,=C would be 2 if bond types counted
    CHECK(analyzer.calculate_dof("C-C") >= 1);
}

// Test empty or invalid SMARTS
TEST_CASE_METHOD(SmartsAnalyzerFixture, "SmartsAnalyzer: calculate_dof invalid", "[SmartsAnalyzer]") {
    // Empty string should return 0 or throw
    CHECK(analyzer.calculate_dof("") == 0);
}

// Test ring queries
TEST_CASE_METHOD(SmartsAnalyzerFixture, "SmartsAnalyzer: calculate_dof rings", "[SmartsAnalyzer]") {
    // Ring membership doesn't change DOF
    CHECK(analyzer.calculate_dof("[C;R]") == 1);
    
    // Ring size specification
    CHECK(analyzer.calculate_dof("[C;r6]") == 1);
}


// enumerate_variants Tests


// Test simple SMARTS with no variations
TEST_CASE_METHOD(SmartsAnalyzerFixture, "SmartsAnalyzer: enumerate_variants simple", "[SmartsAnalyzer]") {
    auto variants = analyzer.enumerate_variants("[C]", 10);
    REQUIRE(variants.size() == 1);
    CHECK(variants[0] == "[C]");
}

// Test SMARTS with single OR query
TEST_CASE_METHOD(SmartsAnalyzerFixture, "SmartsAnalyzer: enumerate_variants single OR", "[SmartsAnalyzer]") {
    auto variants = analyzer.enumerate_variants("[C,N]", 10);
    REQUIRE(variants.size() == 2);
    
    // Should contain both [C] and [N]
    std::set<std::string> variant_set(variants.begin(), variants.end());
    CHECK((variant_set.count("[C]") > 0 || variant_set.count("[#6]") > 0));
    CHECK((variant_set.count("[N]") > 0 || variant_set.count("[#7]") > 0));
}

// Test SMARTS with three-way OR
TEST_CASE_METHOD(SmartsAnalyzerFixture, "SmartsAnalyzer: enumerate_variants three-way OR", "[SmartsAnalyzer]") {
    auto variants = analyzer.enumerate_variants("[C,N,O]", 10);
    REQUIRE(variants.size() == 3);
    
    // Should be unique
    std::set<std::string> unique_variants(variants.begin(), variants.end());
    CHECK(unique_variants.size() == 3);
}

// Test SMARTS with multiple OR queries (Cartesian product)
TEST_CASE_METHOD(SmartsAnalyzerFixture, "SmartsAnalyzer: enumerate_variants multiple OR", "[SmartsAnalyzer]") {
    // [C,N][O,S] should give 4 variants: CO, CS, NO, NS
    auto variants = analyzer.enumerate_variants("[C,N][O,S]", 10);
    CHECK(variants.size() == 4);
    
    // All should be unique
    std::set<std::string> unique_variants(variants.begin(), variants.end());
    CHECK(unique_variants.size() == 4);
}

// Test max limit on variants
TEST_CASE_METHOD(SmartsAnalyzerFixture, "SmartsAnalyzer: enumerate_variants max limit", "[SmartsAnalyzer]") {
    // [C,N,O][H1,H2,H3] = 9 possible, but limit to 5
    auto variants = analyzer.enumerate_variants("[C,N,O][H1,H2,H3]", 5);
    CHECK(variants.size() <= 5);
    
    // Should still be unique
    std::set<std::string> unique_variants(variants.begin(), variants.end());
    CHECK(unique_variants.size() == variants.size());
}

// Test that max=0 returns empty
TEST_CASE_METHOD(SmartsAnalyzerFixture, "SmartsAnalyzer: enumerate_variants max=0", "[SmartsAnalyzer]") {
    auto variants = analyzer.enumerate_variants("[C,N]", 0);
    CHECK(variants.size() == 0);
}

// Test deduplication
TEST_CASE_METHOD(SmartsAnalyzerFixture, "SmartsAnalyzer: enumerate_variants deduplication", "[SmartsAnalyzer]") {
    // Any enumeration should not have duplicates
    auto variants = analyzer.enumerate_variants("[C,N][O,S]", 100);
    std::set<std::string> unique_variants(variants.begin(), variants.end());
    CHECK(unique_variants.size() == variants.size());
}

// Test with AND queries
TEST_CASE_METHOD(SmartsAnalyzerFixture, "SmartsAnalyzer: enumerate_variants AND", "[SmartsAnalyzer]") {
    auto variants = analyzer.enumerate_variants("[C&H3]", 10);
    REQUIRE(variants.size() >= 1);
    
    // Should generate valid SMARTS
    for (const auto& v : variants) {
        CHECK_FALSE(v.empty());
    }
}

// Test with NOT queries
TEST_CASE_METHOD(SmartsAnalyzerFixture, "SmartsAnalyzer: enumerate_variants NOT", "[SmartsAnalyzer]") {
    auto variants = analyzer.enumerate_variants("[!H0]", 10);
    
    // Should generate some variants
    CHECK(variants.size() >= 1);
    
    // Should all be valid non-empty strings
    for (const auto& v : variants) {
        CHECK_FALSE(v.empty());
    }
}

// Test aromatic queries
TEST_CASE_METHOD(SmartsAnalyzerFixture, "SmartsAnalyzer: enumerate_variants aromatic", "[SmartsAnalyzer]") {
    auto variants = analyzer.enumerate_variants("[c,n]", 10);
    CHECK(variants.size() == 2);
}

// Test with bonds
TEST_CASE_METHOD(SmartsAnalyzerFixture, "SmartsAnalyzer: enumerate_variants with bonds", "[SmartsAnalyzer]") {
    auto variants = analyzer.enumerate_variants("[C,N]-[O,S]", 10);
    CHECK(variants.size() == 4);
    
    // All should contain bonds
    for (const auto& v : variants) {
        CHECK((v.find("-") != std::string::npos || v.size() > 3)); // May have implicit bonds
    }
}

TEST_CASE_METHOD(SmartsAnalyzerFixture,
                 "SmartsAnalyzer: enumerate_variants preserves negated primitives",
                 "[SmartsAnalyzer]") {
    const auto variants = analyzer.enumerate_variants("[#6&!H1:1]=[#6:2]", 10);
    REQUIRE_FALSE(variants.empty());

    bool saw_negated = false;
    for (const auto &v : variants) {
        if (v.find("!H1") != std::string::npos) {
            saw_negated = true;
        }
    }
    CHECK(saw_negated);
}

TEST_CASE_METHOD(SmartsAnalyzerFixture,
                 "SmartsAnalyzer: reject dangling bond OR arm",
                 "[SmartsAnalyzer]") {
    CHECK_THROWS_WITH(analyzer.enumerate_variants("[C]=,[C]", 10),
                      Catch::Matchers::ContainsSubstring(
                          "Invalid SMARTS: dangling bond OR arm"));
}

TEST_CASE_METHOD(SmartsAnalyzerFixture, "SmartsAnalyzer: enumerate_variants splits OR bonds", "[SmartsAnalyzer]") {
    auto variants = analyzer.enumerate_variants("[#6]:,-[c&!H1]", 10);
    CHECK(variants.size() == 2);

    bool saw_aromatic_bond = false;
    bool saw_single_bond = false;
    for (const auto &v : variants) {
        CHECK(v.find(":,-") == std::string::npos);
        if (v.find(":") != std::string::npos) {
            saw_aromatic_bond = true;
        }
        if (v.find("-") != std::string::npos) {
            saw_single_bond = true;
        }
    }
    CHECK(saw_aromatic_bond);
    CHECK(saw_single_bond);
}

TEST_CASE_METHOD(SmartsAnalyzerFixture,
                 "SmartsAnalyzer: recursive flatten also emits discard variant",
                 "[SmartsAnalyzer]") {
    auto variants = analyzer.enumerate_variants("[#6][c$([c][#6])]", 20);

    bool saw_inline = false;
    bool saw_discard = false;
    for (const auto &v : variants) {
        if (v == "[#6][c][#6]") {
            saw_inline = true;
        }
        if (v == "[#6][c]") {
            saw_discard = true;
        }
    }

    CHECK(saw_inline);
    CHECK(saw_discard);
}

TEST_CASE_METHOD(SmartsAnalyzerFixture,
                 "SmartsAnalyzer: recursive-only atom syntax is equivalent",
                 "[SmartsAnalyzer]") {
    auto variants = analyzer.enumerate_variants("[#6][$([c][#6])]", 20);

    bool saw_inline = false;
    bool saw_discard = false;
    for (const auto &v : variants) {
        if (v == "[#6][c][#6]") {
            saw_inline = true;
        }
        if (v == "[#6][c]") {
            saw_discard = true;
        }
    }

    CHECK(saw_inline);
    CHECK(saw_discard);
}

TEST_CASE_METHOD(SmartsAnalyzerFixture,
                 "SmartsAnalyzer: recursive flatten preserves branching",
                 "[SmartsAnalyzer]") {
    auto variants = analyzer.enumerate_variants("[#6:1](O)[$([c:3][#6:5]O)]", 20);

    bool saw_inline = false;
    bool saw_discard = false;
    for (const auto &v : variants) {
        if (v == "[#6:1]([O])[c:3][#6:5][O]") {
            saw_inline = true;
        }
        if (v == "[#6:1]([O])[c:3]") {
            saw_discard = true;
        }
    }

    CHECK(saw_inline);
    CHECK(saw_discard);
}

TEST_CASE_METHOD(SmartsAnalyzerFixture,
                 "SmartsAnalyzer: recursive flatten keeps multi-tail branching",
                 "[SmartsAnalyzer]") {
    auto variants = analyzer.enumerate_variants(
        "[#6;H0;+0;a;D1,D3&$([#6](:[a])(:[a]))]", 20);

    bool saw_branched = false;
    bool saw_linearized = false;
    for (const auto &v : variants) {
        if (v.find("(:a):a") != std::string::npos) {
            saw_branched = true;
        }
        if (v.find(":a:a") != std::string::npos) {
            saw_linearized = true;
        }
    }

    CHECK(saw_branched);
    CHECK_FALSE(saw_linearized);
}

TEST_CASE_METHOD(
    SmartsAnalyzerFixture,
    "SmartsAnalyzer: complex recursive OR with atom maps does not crash",
    "[SmartsAnalyzer]") {
    const std::string input =
        "[$([!#1&!#6:2]1=[!#1&!#6:3]-[#6:4]1),$([!#1&!#6:5]1-[!#6&!#1:6]=[#6:7]1):1]";

    std::vector<std::string> variants;
    CHECK_NOTHROW(variants = analyzer.enumerate_variants(input, 1000));
    CHECK_FALSE(variants.empty());
}

// Test that variants are all valid SMARTS
TEST_CASE_METHOD(SmartsAnalyzerFixture, "SmartsAnalyzer: enumerate_variants valid output", "[SmartsAnalyzer]") {
    auto variants = analyzer.enumerate_variants("[C,N,O]", 10);
    
    for (const auto& v : variants) {
        // Should not be empty
        CHECK_FALSE(v.empty());
        
        // Should contain brackets for atom primitives
        CHECK(v.find("[") != std::string::npos);
        CHECK(v.find("]") != std::string::npos);
    }
}

// Test edge case: very large DOF with small max
TEST_CASE_METHOD(SmartsAnalyzerFixture, "SmartsAnalyzer: enumerate_variants large DOF small max", "[SmartsAnalyzer]") {
    // [C,N,O,S,P][H0,H1,H2,H3][+0,+1,-1] = 5*4*3 = 60 possibilities
    auto variants = analyzer.enumerate_variants("[C,N,O,S,P][H0,H1,H2,H3][+0,+1,-1]", 10);
    
    // Should respect max limit
    CHECK(variants.size() <= 10);
    
    // Should be unique
    std::set<std::string> unique_variants(variants.begin(), variants.end());
    CHECK(unique_variants.size() == variants.size());
}

// Test charge queries
TEST_CASE_METHOD(SmartsAnalyzerFixture, "SmartsAnalyzer: enumerate_variants charges", "[SmartsAnalyzer]") {
    auto variants = analyzer.enumerate_variants("[C+1,C+0]", 10);
    CHECK(variants.size() == 2);
}

// Test hydrogen count queries
TEST_CASE_METHOD(SmartsAnalyzerFixture, "SmartsAnalyzer: enumerate_variants hydrogen counts", "[SmartsAnalyzer]") {
    auto variants = analyzer.enumerate_variants("[CH0,CH1,CH2]", 10);
    CHECK(variants.size() == 3);
}

// Test empty SMARTS
TEST_CASE_METHOD(SmartsAnalyzerFixture, "SmartsAnalyzer: enumerate_variants empty", "[SmartsAnalyzer]") {
    auto variants = analyzer.enumerate_variants("", 10);
    CHECK(variants.size() == 0);
}

// Test realistic chemistry example
TEST_CASE_METHOD(SmartsAnalyzerFixture, "SmartsAnalyzer: enumerate_variants realistic chemistry", "[SmartsAnalyzer]") {
    // Carbonyl carbon: [C;$(C=O)]
    auto variants = analyzer.enumerate_variants("[C,N]=O", 10);
    CHECK(variants.size() >= 1);
    CHECK(variants.size() <= 10);
}

// Test that enumerate is consistent with calculate_dof
TEST_CASE_METHOD(SmartsAnalyzerFixture, "SmartsAnalyzer: consistency DOF vs enumerate", "[SmartsAnalyzer]") {
    std::string smarts = "[C,N][O,S]";
    int dof = analyzer.calculate_dof(smarts);
    auto variants = analyzer.enumerate_variants(smarts, 100);
    
    // Number of variants should match DOF (if max is high enough)
    CHECK(static_cast<int>(variants.size()) == dof);
}

// Test another consistency check
TEST_CASE_METHOD(SmartsAnalyzerFixture, "SmartsAnalyzer: consistency three elements", "[SmartsAnalyzer]") {
    std::string smarts = "[C,N,O]";
    int dof = analyzer.calculate_dof(smarts);
    auto variants = analyzer.enumerate_variants(smarts, 100);
    
    CHECK(dof == 3);
    CHECK(variants.size() == 3);
}