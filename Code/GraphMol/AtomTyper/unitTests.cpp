#include <GraphMol/AtomTyper/atom_typer.hpp>
#include <GraphMol/AtomTyper/expression_builder.hpp>
#include <GraphMol/AtomTyper/reaction_extractor.hpp>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <catch2/catch_all.hpp>

#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

  TEST_CASE("ExpressionBuilder: SmilesToSmarts minimal", "[ExpressionBuilder]") {
    const std::string smarts = atom_typer::smiles_to_smarts("CCO", {"element"});
    std::unique_ptr<RDKit::ROMol> query(RDKit::SmartsToMol(smarts));
    std::unique_ptr<RDKit::ROMol> target(RDKit::SmilesToMol("CCO"));
    REQUIRE(query != nullptr);
    REQUIRE(target != nullptr);
    CHECK(query->getNumAtoms() == target->getNumAtoms());
    CHECK(query->getNumBonds() == target->getNumBonds());
    CHECK(RDKit::SubstructMatch(*target, *query).size() > 0);
  }

  TEST_CASE("ExpressionBuilder: SmilesToSmarts standard", "[ExpressionBuilder]") {
    const std::string smarts = atom_typer::smiles_to_smarts("CCO", {"element", "D"});
    std::unique_ptr<RDKit::ROMol> query(RDKit::SmartsToMol(smarts));
    std::unique_ptr<RDKit::ROMol> target(RDKit::SmilesToMol("CCO"));
    REQUIRE(query != nullptr);
    REQUIRE(target != nullptr);
    CHECK(query->getNumAtoms() == target->getNumAtoms());
    CHECK(query->getNumBonds() == target->getNumBonds());
    CHECK(RDKit::SubstructMatch(*target, *query).size() > 0);
  }

  TEST_CASE("ExpressionBuilder: SmilesToSmarts detailed", "[ExpressionBuilder]") {
    const std::string smarts = atom_typer::smiles_to_smarts("CCO", {"element", "D", "H"});
    std::unique_ptr<RDKit::ROMol> query(RDKit::SmartsToMol(smarts));
    std::unique_ptr<RDKit::ROMol> target(RDKit::SmilesToMol("CCO"));
    REQUIRE(query != nullptr);
    REQUIRE(target != nullptr);
    CHECK(query->getNumAtoms() == target->getNumAtoms());
    CHECK(query->getNumBonds() == target->getNumBonds());
    CHECK(RDKit::SubstructMatch(*target, *query).size() > 0);
  }

  TEST_CASE("ExpressionBuilder: SmilesToSmarts aromatic minimal", "[ExpressionBuilder]") {
    const std::string smarts = atom_typer::smiles_to_smarts("c1ccccc1", {"element", "a", "D"});
    std::unique_ptr<RDKit::ROMol> query(RDKit::SmartsToMol(smarts));
    std::unique_ptr<RDKit::ROMol> target(RDKit::SmilesToMol("c1ccccc1"));
    REQUIRE(query != nullptr);
    REQUIRE(target != nullptr);
    CHECK(query->getNumAtoms() == target->getNumAtoms());
    CHECK(query->getNumBonds() == target->getNumBonds());
    CHECK(RDKit::SubstructMatch(*target, *query).size() > 0);
  }

  TEST_CASE("ExpressionBuilder: SmilesToSmarts branched standard", "[ExpressionBuilder]") {
    const std::string smarts = atom_typer::smiles_to_smarts("CC(C)C", {"element", "D"});
    std::unique_ptr<RDKit::ROMol> query(RDKit::SmartsToMol(smarts));
    std::unique_ptr<RDKit::ROMol> target(RDKit::SmilesToMol("CC(C)C"));
    REQUIRE(query != nullptr);
    REQUIRE(target != nullptr);
    CHECK(query->getNumAtoms() == target->getNumAtoms());
    CHECK(query->getNumBonds() == target->getNumBonds());
    CHECK(RDKit::SubstructMatch(*target, *query).size() > 0);
  }

  TEST_CASE("ExpressionBuilder: SmilesToSmarts charged detailed", "[ExpressionBuilder]") {
    const std::string smarts = atom_typer::smiles_to_smarts("[NH4+]", {"element", "D", "H", "charge"});
    std::unique_ptr<RDKit::ROMol> query(RDKit::SmartsToMol(smarts));
    std::unique_ptr<RDKit::ROMol> target(RDKit::SmilesToMol("[NH4+]"));
    REQUIRE(query != nullptr);
    REQUIRE(target != nullptr);
    CHECK(query->getNumAtoms() == target->getNumAtoms());
    CHECK(query->getNumBonds() == target->getNumBonds());
    CHECK(RDKit::SubstructMatch(*target, *query).size() > 0);
  }

  TEST_CASE("ExpressionBuilder: SmilesToSmarts nitrogen detailed", "[ExpressionBuilder]") {
    const std::string smarts = atom_typer::smiles_to_smarts("CCN", {"element", "D", "H"});
    std::unique_ptr<RDKit::ROMol> query(RDKit::SmartsToMol(smarts));
    std::unique_ptr<RDKit::ROMol> target(RDKit::SmilesToMol("CCN"));
    REQUIRE(query != nullptr);
    REQUIRE(target != nullptr);
    CHECK(query->getNumAtoms() == target->getNumAtoms());
    CHECK(query->getNumBonds() == target->getNumBonds());
    CHECK(RDKit::SubstructMatch(*target, *query).size() > 0);
  }

  TEST_CASE("ExpressionBuilder: SmilesToSmarts empty string throws", "[ExpressionBuilder]") {
    CHECK_THROWS_AS(atom_typer::smiles_to_smarts("", {"element"}), std::invalid_argument);
  }

  TEST_CASE("ExpressionBuilder: SmilesToSmarts invalid SMILES throws", "[ExpressionBuilder]") {
    CHECK_THROWS_AS(atom_typer::smiles_to_smarts("CCCZ", {"element"}), std::invalid_argument);
  }

  TEST_CASE("ExpressionBuilder: SmilesToSmarts empty primitive list throws", "[ExpressionBuilder]") {
    CHECK_THROWS_AS(atom_typer::smiles_to_smarts("CCO", {}), std::invalid_argument);
  }

  TEST_CASE("ExpressionBuilder: SmilesToSmarts sulfur minimal", "[ExpressionBuilder]") {
    const std::string smarts = atom_typer::smiles_to_smarts("CCS", {"element"});
    std::unique_ptr<RDKit::ROMol> query(RDKit::SmartsToMol(smarts));
    std::unique_ptr<RDKit::ROMol> target(RDKit::SmilesToMol("CCS"));
    REQUIRE(query != nullptr);
    REQUIRE(target != nullptr);
    CHECK(query->getNumAtoms() == target->getNumAtoms());
    CHECK(query->getNumBonds() == target->getNumBonds());
    CHECK(RDKit::SubstructMatch(*target, *query).size() > 0);
  }

  TEST_CASE("ExpressionBuilder: SmilesToSmarts multiple heteroatoms standard", "[ExpressionBuilder]") {
    const std::string smarts = atom_typer::smiles_to_smarts("C(N)O", {"element", "D"});
    std::unique_ptr<RDKit::ROMol> query(RDKit::SmartsToMol(smarts));
    std::unique_ptr<RDKit::ROMol> target(RDKit::SmilesToMol("C(N)O"));
    REQUIRE(query != nullptr);
    REQUIRE(target != nullptr);
    CHECK(query->getNumAtoms() == target->getNumAtoms());
    CHECK(query->getNumBonds() == target->getNumBonds());
    CHECK(RDKit::SubstructMatch(*target, *query).size() > 0);
  }

  TEST_CASE("ExpressionBuilder: primitive list supports bracket syntax", "[ExpressionBuilder]") {
    const std::string smarts = atom_typer::smiles_to_smarts("CCO", {"[atomic_number,D,R]"});
    std::unique_ptr<RDKit::ROMol> query(RDKit::SmartsToMol(smarts));
    std::unique_ptr<RDKit::ROMol> target(RDKit::SmilesToMol("CCO"));
    REQUIRE(query != nullptr);
    REQUIRE(target != nullptr);
    CHECK(query->getNumBonds() == target->getNumBonds());
    CHECK(RDKit::SubstructMatch(*target, *query).size() > 0);
  }

  TEST_CASE("ExpressionBuilder: primitive list without element", "[ExpressionBuilder]") {
    const std::string smarts = atom_typer::smiles_to_smarts("C1CC1", {"[R,D]"});
    std::unique_ptr<RDKit::ROMol> query(RDKit::SmartsToMol(smarts));
    std::unique_ptr<RDKit::ROMol> target(RDKit::SmilesToMol("C1CC1"));
    REQUIRE(query != nullptr);
    REQUIRE(target != nullptr);
    CHECK(query->getNumBonds() == target->getNumBonds());
    CHECK(RDKit::SubstructMatch(*target, *query).size() > 0);
  }

  TEST_CASE("ExpressionBuilder: primitive list with charge", "[ExpressionBuilder]") {
    const std::string smarts = atom_typer::smiles_to_smarts("[NH4+]", {"[charge, R, D]"});
    std::unique_ptr<RDKit::ROMol> query(RDKit::SmartsToMol(smarts));
    std::unique_ptr<RDKit::ROMol> target(RDKit::SmilesToMol("[NH4+]"));
    REQUIRE(query != nullptr);
    REQUIRE(target != nullptr);
    CHECK(query->getNumBonds() == target->getNumBonds());
    CHECK(RDKit::SubstructMatch(*target, *query).size() > 0);
  }

  TEST_CASE("ExpressionBuilder: primitive list includes total degree and ring bond count", "[ExpressionBuilder]") {
    const std::string smarts = atom_typer::smiles_to_smarts("C1CC1", {"[X,x]"});
    std::unique_ptr<RDKit::ROMol> query(RDKit::SmartsToMol(smarts));
    std::unique_ptr<RDKit::ROMol> target(RDKit::SmilesToMol("C1CC1"));
    REQUIRE(query != nullptr);
    REQUIRE(target != nullptr);
    CHECK(query->getNumBonds() == target->getNumBonds());
    CHECK(RDKit::SubstructMatch(*target, *query).size() > 0);
  }

  TEST_CASE("ExpressionBuilder: primitive list rejects unknown primitive", "[ExpressionBuilder]") {
    CHECK_THROWS_AS(atom_typer::smiles_to_smarts("CCO", {"[X,not_a_primitive]"}),
            std::invalid_argument);
  }

  TEST_CASE("ExpressionBuilder: atomic number is always first when present", "[ExpressionBuilder]") {
    const std::string smarts = atom_typer::smiles_to_smarts("CCO", {"[D,R,atomic_number]"});
    CHECK(smarts.find("[#6") != std::string::npos);
    CHECK(smarts.find("-") != std::string::npos);
  }

  TEST_CASE("ExpressionBuilder: element symbol is moved ahead of other primitives", "[ExpressionBuilder]") {
    const std::string smarts = atom_typer::smiles_to_smarts("CCO", {"[D,R,element]"});
    CHECK(smarts.find("[C") != std::string::npos);
    CHECK(smarts.find("-") != std::string::npos);
  }

  TEST_CASE("ExpressionBuilder: batch smiles_to_smarts by primitive list", "[ExpressionBuilder]") {
    const std::vector<std::string> smarts = atom_typer::smiles_to_smarts(
      std::vector<std::string>{"CCO", "C1CC1"}, {"[atomic_number,D]"});
    REQUIRE(smarts.size() == 2);

    std::unique_ptr<RDKit::ROMol> query1(RDKit::SmartsToMol(smarts[0]));
    std::unique_ptr<RDKit::ROMol> target1(RDKit::SmilesToMol("CCO"));
    REQUIRE(query1 != nullptr);
    REQUIRE(target1 != nullptr);
    CHECK(query1->getNumBonds() == target1->getNumBonds());
    CHECK(RDKit::SubstructMatch(*target1, *query1).size() > 0);

    std::unique_ptr<RDKit::ROMol> query2(RDKit::SmartsToMol(smarts[1]));
    std::unique_ptr<RDKit::ROMol> target2(RDKit::SmilesToMol("C1CC1"));
    REQUIRE(query2 != nullptr);
    REQUIRE(target2 != nullptr);
    CHECK(query2->getNumBonds() == target2->getNumBonds());
    CHECK(RDKit::SubstructMatch(*target2, *query2).size() > 0);
  }

  TEST_CASE("ExpressionBuilder: batch smiles_to_smarts handles empty list", "[ExpressionBuilder]") {
    const auto smarts = atom_typer::smiles_to_smarts(std::vector<std::string>{},
                             std::vector<std::string>{"element"});
    CHECK(smarts.empty());
  }

  TEST_CASE("ExpressionBuilder: atom-centered deduplicate suppresses duplicates", "[ExpressionBuilder]") {
    const auto smarts = atom_typer::smiles_to_atom_centered_smarts(
      "c1ccccc1", {"[atomic_number,H,D]"}, 1, true, false, true);

    REQUIRE(smarts.size() == 1);
    CHECK(smarts.front().find("[*]") != std::string::npos);
  }

  TEST_CASE("ExpressionBuilder: batch atom-centered deduplicate is global across batch", "[ExpressionBuilder]") {
    const auto smarts = atom_typer::smiles_to_atom_centered_smarts(
      std::vector<std::string>{"C", "C"}, {"[atomic_number,H,D]"}, 1,
      true, false, true);

    REQUIRE(smarts.size() == 2);
    CHECK(smarts[0].size() == 1);
    CHECK(smarts[1].empty());
  }

  TEST_CASE("ExpressionBuilder: atom-centered wildcard neighbors preserve chirality", "[ExpressionBuilder]") {
    const auto smarts = atom_typer::smiles_to_atom_centered_smarts(
      "CC[C@H](O)C",
      {"atom", "A", "D", "H", "charge", "v", "R", "r", "X", "x", "^"},
      1, true, false, false);

    REQUIRE(smarts.size() == 5);
    CHECK(smarts[2].find("[#6@") != std::string::npos);

    std::unique_ptr<RDKit::ROMol> query(RDKit::SmartsToMol(smarts[2]));
    REQUIRE(query != nullptr);
  }

  TEST_CASE("ExpressionBuilder: atom-centered preserves lexical chirality from minimal smiles", "[ExpressionBuilder]") {
    const auto smarts = atom_typer::smiles_to_atom_centered_smarts(
      "[C@H][C]",
      {"atom", "A", "D", "H", "charge", "v", "R", "r", "X", "x", "^"},
      1, true, false, false);

    REQUIRE(smarts.size() == 2);
    CHECK(smarts[0].find("[#6@") != std::string::npos);

    std::unique_ptr<RDKit::ROMol> query(RDKit::SmartsToMol(smarts[0]));
    REQUIRE(query != nullptr);
  }

  TEST_CASE("ReactionExtractor: reaction smiles extraction", "[ReactionExtractor]") {
    const std::string rxn = "[CH3:1][Br:2].[OH:3]>>[CH3:1][OH:3].[Br-:2]";
    const auto results = atom_typer::extract_single_root_template(
      rxn, {"[atomic_number,D,H,charge]"}, 2, false, true);

    REQUIRE(results.size() == 3);
    CHECK(std::get<0>(results[0]) == 0);
    CHECK(std::get<0>(results[1]) == 1);
    CHECK(std::get<0>(results[2]) == 2);
    CHECK(std::get<1>(results[0]).find(">>") != std::string::npos);
    CHECK(std::get<1>(results[0]).find(":1") != std::string::npos);
    CHECK(std::get<1>(results[0]).find("#") != std::string::npos);
  }

  TEST_CASE("ReactionExtractor: primitive list required", "[ReactionExtractor]") {
    const std::string rxn = "[CH3:1][Br:2].[OH:3]>>[CH3:1][OH:3].[Br-:2]";
    CHECK_THROWS_AS(atom_typer::extract_single_root_template(
              rxn, {}, 1, false, true),
            std::invalid_argument);
  }

  TEST_CASE("ReactionExtractor: invalid reaction text throws", "[ReactionExtractor]") {
    CHECK_THROWS_AS(atom_typer::extract_single_root_template(
              "not-a-reaction", {"element"}, 1, false, true),
            std::invalid_argument);
  }

  struct AtomTyperFixture {
    atom_typer::AtomTyper typer;
  };

  static std::vector<atom_typer::AtomType> extract_atom_types(
    const std::vector<atom_typer::PatternItem> &pattern_items) {
    std::vector<atom_typer::AtomType> atom_types;
    atom_types.reserve(pattern_items.size());
    for (const auto &item : pattern_items) {
      if (item.kind == atom_typer::PatternItemKind::Atom) {
        atom_types.push_back(item.atom);
      }
    }
    return atom_types;
  }

  TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: Basic SMILES (ethanol)", "[AtomTyper]") {
    std::string smiles = "CCO";
    auto atom_types = extract_atom_types(typer.type_atoms_from_smiles(smiles));

    REQUIRE(atom_types.size() == 3);
    CHECK(atom_types[0].atomic_number == 6);
    CHECK(atom_types[0].min_bonds == 1);
    CHECK(atom_types[1].atomic_number == 6);
    CHECK(atom_types[1].min_bonds == 2);
    CHECK(atom_types[2].atomic_number == 8);
    CHECK(atom_types[2].min_bonds == 1);
  }

  TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: Aromatic benzene", "[AtomTyper]") {
    std::string smiles = "c1ccccc1";
    auto atom_types = extract_atom_types(typer.type_atoms_from_smiles(smiles));

    REQUIRE(atom_types.size() == 6);
    for (const auto& at : atom_types) {
      CHECK(at.atomic_number == 6);
      CHECK(at.is_aromatic);
      CHECK(at.is_in_ring);
      CHECK(at.ring_size == 6);
    }
  }

  TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: Ring detection", "[AtomTyper]") {
    std::string smiles = "C1CCCCC1";
    auto atom_types = extract_atom_types(typer.type_atoms_from_smiles(smiles));

    REQUIRE(atom_types.size() == 6);
    for (const auto& at : atom_types) {
      CHECK(at.is_in_ring);
      CHECK(at.ring_size == 6);
      CHECK_FALSE(at.is_aromatic);
    }
  }

  TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: Hybridization (SP2)", "[AtomTyper]") {
    std::string smiles = "C=C";
    auto atom_types = extract_atom_types(typer.type_atoms_from_smiles(smiles));

    REQUIRE(atom_types.size() == 2);
    for (const auto& at : atom_types) {
      CHECK(at.hybridization == "SP2");
    }
  }

  TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: Hybridization (SP)", "[AtomTyper]") {
    std::string smiles = "C#C";
    auto atom_types = extract_atom_types(typer.type_atoms_from_smiles(smiles));

    REQUIRE(atom_types.size() == 2);
    for (const auto& at : atom_types) {
      CHECK(at.hybridization == "SP");
    }
  }

  TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: Formal charge", "[AtomTyper]") {
    std::string smiles = "[NH4+]";
    auto atom_types = extract_atom_types(typer.type_atoms_from_smiles(smiles));

    REQUIRE(atom_types.size() == 1);
    CHECK(atom_types[0].atomic_number == 7);
    CHECK(atom_types[0].formal_charge == 1);
  }

  TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: Basic SMARTS", "[AtomTyper]") {
    std::string smarts = "[C,N]";
    CHECK_NOTHROW(typer.type_atoms_from_smarts(smarts));
  }

  TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: SMARTS explicit H defaults to 0", "[AtomTyper]") {
    std::string smarts = "[C]";
    auto atom_types = extract_atom_types(typer.type_pattern_from_smarts(smarts));

    REQUIRE(atom_types.size() == 1);
    CHECK(atom_types[0].num_hydrogens == 0);
  }

  TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: SMARTS #6 supports aromatic and aliphatic", "[AtomTyper]") {
    std::string smarts = "[#6]";
    auto atom_types = extract_atom_types(typer.type_pattern_from_smarts(smarts));

    REQUIRE(atom_types.size() == 1);
    CHECK(atom_types[0].is_aromatic);
    CHECK(atom_types[0].is_aliphatic);
    CHECK(atom_types[0].smarts_pattern.find("A,a") != std::string::npos);
  }

  TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: SMARTS explicit H count preserved", "[AtomTyper]") {
    std::string smarts = "[CH3]";
    auto atom_types = extract_atom_types(typer.type_pattern_from_smarts(smarts));

    REQUIRE(atom_types.size() == 1);
    CHECK(atom_types[0].num_hydrogens == 3);
  }

  TEST_CASE_METHOD(AtomTyperFixture,
           "AtomTyper: type_atoms_from_smarts emits expected flags",
           "[AtomTyper]") {
    const std::string smarts = "[CH2][CX3]=[O]";
    const std::string enumerated = typer.type_atoms_from_smarts(smarts);

    CHECK_FALSE(enumerated.empty());
    CHECK(enumerated.find("D") != std::string::npos);
    CHECK(enumerated.find("H") != std::string::npos);
    const bool has_charge_token =
      enumerated.find("+") != std::string::npos ||
      enumerated.find("-") != std::string::npos;
    CHECK(has_charge_token);
  }

  TEST_CASE_METHOD(AtomTyperFixture,
        "AtomTyper: type_atoms_from_smarts preserves branching",
           "[AtomTyper]") {
    const std::string smarts = "[#6:1](O)[$([c:3][#6:5]O)]";
      const std::string enumerated = typer.type_atoms_from_smarts(smarts);

    CHECK_FALSE(enumerated.empty());
    CHECK(enumerated.find("[") != std::string::npos);
    CHECK(enumerated.find("]") != std::string::npos);
    CHECK(enumerated.find("([") != std::string::npos);
  }

  TEST_CASE_METHOD(AtomTyperFixture,
        "AtomTyper: type_atoms_from_smarts maps newly introduced atoms",
           "[AtomTyper][enumerate_dof_map_new_atoms]") {
    const std::string smarts = "[C:1]=[C:2]";
      int max_amap = 2;
      const std::string enumerated = typer.type_atoms_from_smarts(smarts, true, max_amap);

    CHECK_FALSE(enumerated.empty());
    CHECK(enumerated.find(":1") != std::string::npos);
    CHECK(enumerated.find(":2") != std::string::npos);
  }

  TEST_CASE_METHOD(AtomTyperFixture,
        "AtomTyper: type_atoms_from_smarts handles complex bond OR variants",
           "[AtomTyper][enumerate_dof_complex_bonds]") {
    const std::vector<std::string> variants = {
      "[#6:1]=[#6:2]",
      "[#6:1]-[#6:2]",
      "[#6:1]:[#6:2]",
    };

    REQUIRE(variants.size() == 3);
    for (const auto& v : variants) {
      INFO("Variant: " << v);
      const auto enumerated = typer.type_atoms_from_smarts(v);
      CHECK_FALSE(enumerated.empty());
    }
  }

  TEST_CASE_METHOD(AtomTyperFixture,
           "AtomTyper: aromatic constrained atom is not collapsed",
           "[AtomTyper]") {
    const std::string smarts = "[#6:1]:[c&H2&+0]";
    CHECK_THROWS_WITH(
      typer.type_atoms_from_smarts(smarts),
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
           "AtomTyper: inspect_tautomer recognizes nitro-like motif",
           "[AtomTyper]") {
    CHECK(typer.inspect_tautomer("[#7](=[#8])=[#8]"));
    CHECK_FALSE(typer.inspect_tautomer("[#6]-[#6]"));
  }

  TEST_CASE_METHOD(AtomTyperFixture,
           "AtomTyper: is_valid_valence_smarts allows tautomer override",
           "[AtomTyper]") {
    const std::string tautomeric = "[#7](=[#8])=[#8]";
    CHECK(typer.is_valid_valence_smarts(tautomeric));
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
           "AtomTyper: type_atoms_from_smarts preserves negated atomic-number atoms",
           "[AtomTyper]") {
    const std::string smarts = "[!#6:1]=[!#6:2]";
    const std::string enumerated = typer.type_atoms_from_smarts(smarts);

    CHECK_FALSE(enumerated.empty());
    CHECK(enumerated.find("!#6") != std::string::npos);
    CHECK(enumerated.find("#0") == std::string::npos);
  }

  TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: SMARTS bond counts and hybridization", "[AtomTyper]") {
    std::string smarts = "[C]=[C]";
    auto atom_types = extract_atom_types(typer.type_pattern_from_smarts(smarts));

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
    auto atom_types = extract_atom_types(typer.type_atoms_from_smiles(smiles));

    REQUIRE(atom_types.size() == 3);
    CHECK(atom_types[0].remaining_valence == 3);
    CHECK(atom_types[1].remaining_valence == 2);
    CHECK(atom_types[2].remaining_valence == 1);
  }

  TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: Invalid SMILES throws", "[AtomTyper]") {
    std::string invalid_smiles = "C(C";
    CHECK_THROWS_AS(typer.type_atoms_from_smiles(invalid_smiles), std::runtime_error);
  }

  TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: Empty SMILES throws", "[AtomTyper]") {
    std::string empty_smiles = "";
    CHECK_THROWS_AS(typer.type_atoms_from_smiles(empty_smiles), std::runtime_error);
  }

  TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: Neighbor detection", "[AtomTyper]") {
    std::string smiles = "CCC";
    auto atom_types = extract_atom_types(typer.type_atoms_from_smiles(smiles));

    REQUIRE(atom_types.size() == 3);
    CHECK(atom_types[0].neighbors.size() == 1);
    CHECK(atom_types[1].neighbors.size() == 2);
    CHECK(atom_types[2].neighbors.size() == 1);
  }

  TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: SMARTS pattern generation", "[AtomTyper]") {
    std::string smiles = "CCO";
    auto atom_types = extract_atom_types(typer.type_atoms_from_smiles(smiles));

    REQUIRE(atom_types.size() == 3);
    for (const auto& at : atom_types) {
      CHECK_FALSE(at.smarts_pattern.empty());
      CHECK(at.smarts_pattern.front() == '[');
      CHECK(at.smarts_pattern.back() == ']');
    }
  }

  TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: get_atom_types_string", "[AtomTyper]") {
    std::string smiles = "CC";
    auto atom_types = extract_atom_types(typer.type_atoms_from_smiles(smiles));
    std::string output = typer.get_atom_types_string(atom_types);

    CHECK_FALSE(output.empty());
    CHECK(output.find("Atom Types") != std::string::npos);
  }

  TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: canonical setting", "[AtomTyper]") {
    typer.set_use_canonical(true);
    std::string smiles = "CCO";
    CHECK_NOTHROW(typer.type_atoms_from_smiles(smiles));
  }

  TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: ring membership list", "[AtomTyper]") {
    std::string smiles = "C1CCCCC1";
    auto atom_types = extract_atom_types(typer.type_atoms_from_smiles(smiles));

    REQUIRE(atom_types.size() == 6);
    for (const auto& at : atom_types) {
      CHECK(at.ring_membership_list.size() == 6);
    }
  }

  TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: number of ring bonds", "[AtomTyper]") {
    std::string smiles = "C1CCCCC1";
    auto atom_types = extract_atom_types(typer.type_atoms_from_smiles(smiles));

    REQUIRE(atom_types.size() == 6);
    for (const auto& at : atom_types) {
      CHECK(at.num_ring_bonds == 2);
    }
  }

  TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: number of aliphatic rings (cyclohexane)", "[AtomTyper]") {
    std::string smiles = "C1CCCCC1";
    auto atom_types = extract_atom_types(typer.type_atoms_from_smiles(smiles));

    REQUIRE(atom_types.size() == 6);
    for (const auto& at : atom_types) {
      CHECK(at.num_aliphatic_rings == 1);
    }
  }

  TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: number of aliphatic rings (benzene)", "[AtomTyper]") {
    std::string smiles = "c1ccccc1";
    auto atom_types = extract_atom_types(typer.type_atoms_from_smiles(smiles));

    REQUIRE(atom_types.size() == 6);
    for (const auto& at : atom_types) {
      CHECK(at.num_aliphatic_rings == 0);
    }
  }

  TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: number of aromatic rings (benzene)", "[AtomTyper]") {
    std::string smiles = "c1ccccc1";
    auto atom_types = extract_atom_types(typer.type_atoms_from_smiles(smiles));

    REQUIRE(atom_types.size() == 6);
    for (const auto& at : atom_types) {
      CHECK(at.num_aromatic_rings == 1);
    }
  }

  TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: number of aromatic rings (cyclohexane)", "[AtomTyper]") {
    std::string smiles = "C1CCCCC1";
    auto atom_types = extract_atom_types(typer.type_atoms_from_smiles(smiles));

    REQUIRE(atom_types.size() == 6);
    for (const auto& at : atom_types) {
      CHECK(at.num_aromatic_rings == 0);
    }
  }

  TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: ring connectivity", "[AtomTyper]") {
    std::string smiles = "C1CCCCC1";
    auto atom_types = extract_atom_types(typer.type_atoms_from_smiles(smiles));

    REQUIRE(atom_types.size() == 6);
    for (const auto& at : atom_types) {
      CHECK(at.ring_connectivity == 2);
    }
  }

  TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: bond types map", "[AtomTyper]") {
    std::string smiles = "CC=O";
    auto atom_types = extract_atom_types(typer.type_atoms_from_smiles(smiles));

    REQUIRE(atom_types.size() == 3);
    const auto& carbonyl_carbon = atom_types[1];
    CHECK(carbonyl_carbon.bond_types.at(1) == 1);
    CHECK(carbonyl_carbon.bond_types.at(2) == 1);
  }

  TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: chirality", "[AtomTyper]") {
    std::string smiles = "CC[C@H](O)C";
    auto atom_types = extract_atom_types(typer.type_atoms_from_smiles(smiles));

    REQUIRE(atom_types.size() == 5);
    const auto& chiral_carbon = atom_types[2];
    CHECK(chiral_carbon.chirality == "R");
  }

  TEST_CASE_METHOD(AtomTyperFixture, "AtomTyper: atom type string content", "[AtomTyper]") {
    std::string smiles = "CCO";
    auto atom_types = extract_atom_types(typer.type_atoms_from_smiles(smiles));
    std::string output = typer.get_atom_types_string(atom_types);

    CHECK(output.find("Atom 0: Element=6") != std::string::npos);
    CHECK(output.find("Atom 1: Element=6") != std::string::npos);
    CHECK(output.find("Atom 2: Element=8") != std::string::npos);
  }
