#include <GraphMol/AtomTyper/smarts_analyzer.hpp>
#include <GraphMol/AtomTyper/atom_typer.hpp>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <iostream>
#include <string>
#include <vector>

int main() {
  const std::vector<std::string> smarts_list = {
      "[CH2][CX3]=[O]", "[c][C]", "[#6][c][C]=[O;H0;+0]",
      "[CX4]([F,Cl,Br,I])([F,Cl,Br,I])", "[#6X3](=[OX1])[#8X2][#6X3](=[OX1])"};

  // const std::vector<std::string> smarts_list = {"[#6][c][C]=[O;H0;+0]"};

  atom_typer::SmartsAnalyzer sa;
  atom_typer::AtomTyper at;
  std::vector<std::string> results = sa.standard_smarts(smarts_list, true);
  for (int i = 0; i < results.size(); ++i) {
    const auto &s = results[i];
    std::cout << "input: " << smarts_list[i] << std::endl;
    std::cout << "output: " << s << std::endl;
  }
  return 0;
}
