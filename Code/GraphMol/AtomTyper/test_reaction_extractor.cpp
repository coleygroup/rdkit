#include "reaction_extractor.hpp"

#include <catch2/catch_all.hpp>

#include <string>
#include <tuple>
#include <vector>

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
