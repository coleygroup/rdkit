#include <GraphMol/AtomTyper/smarts_analyzer.hpp>
#include <GraphMol/AtomTyper/atom_typer.hpp>
#include <iostream>
#include <string>
#include <vector>

int main() {
  atom_typer::SmartsAnalyzer sa;
  atom_typer::AtomTyper at;
  // const std::string test1 = "[#6]:[c]";

  // const auto pattern_types = at.type_pattern_from_smarts(test1);
  // std::cout << "Pattern item count: " << pattern_types.size() << std::endl;
  // std::cout << at.get_pattern_types_string(pattern_types) << std::endl;

  // const auto smarts_from_pattern =
  //     at.get_smarts_from_pattern_types(pattern_types);

  // std::cout << smarts_from_pattern << std::endl;

  // std::cout << sa.calculate_dof(test1) << std::endl;

  // const std::vector<std::string> dof_example = {"[CH2][CX3]=[O]", "[c][C]=[O;H0;+0]", "[#6][c][C]=[O;H0;+0]"};


    const std::vector<std::string> dof_example = { "[#6][c]"};


  for (const auto &dof_example : dof_example) {
    const auto pattern_types = at.type_pattern_from_smarts(dof_example);
      std::cout << "Pattern item count: " << pattern_types.size() << std::endl;
      std::cout << at.get_pattern_types_string(pattern_types) << std::endl;
  
      // std::cout << "SMARTS from pattern: " << at.get_smarts_from_pattern_types(pattern_types) << std::endl;
    std::cout << "Enumerated DOF SMARTS for " << dof_example << ":\n"
              << at.enumerate_dof_smarts(dof_example) << std::endl;
  }
  // const auto pattern_types2 = at.type_pattern_from_smarts(dof_example);

  // const auto variants = sa.enumerate_variants(test1, 1000);
  // for (const auto &variant : variants) {
  //   std::cout << variant << std::endl;
  // }

  //   if(atom_type.kind == atom_typer::PatternItemKind::Atom) {
  //     std::cout << "Atom: " << atom_type.atom.smarts_pattern << std::endl;
  //   } else {
  //     std::cout << "Bond: " << atom_type.bond.smarts_pattern << std::endl;
  //   }
  //   // std::cout << atom_type.atom.smarts_pattern << std::endl;
  // }

  return 0;
}
