#include <GraphMol/AtomTyper/smarts_analyzer.hpp>
#include <GraphMol/AtomTyper/atom_typer.hpp>
#include <GraphMol/AtomTyper/smarts_analyzer_smoke_benchmark.hpp>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <iostream>
#include <string>
#include <vector>

int main() {
  // const std::vector<std::string> smarts_list = {
  //     "[CH2][CX3]=[O]",
  //     "[c][C]",
  //     "[#6][c][C]=[O;H0;+0]",
  //     "[CX4]([F,Cl,Br,I])([F,Cl,Br,I])",
  //     "[#6X3](=[OX1])[#8X2][#6X3](=[OX1])",
  //     "[#6R][CX3R](=[SX1])[NX3R]([#1,#6&!$([CX3]=[OX1,SX1,NX2])])[#6R&$([CX3R]=[OX1,SX1,NX2])]",
  //     "[CH2,$([CH1][#6]),$([CX4]([#6])[#6])]([NX3]([#1,*])[#1,*])[NX3]([#1,*])[#1,*]",
  //     "[#1,#6][CX3]([#1,#6])=[NX2][NX3]([#1,#6&!$([CX3]=[OX1,SX1,NX2])])[CX3](=[SX1])[NX3]([#1,#6&!$([CX3]=[OX1,SX1,NX2])])[#1,#6&!$([CX3]=[OX1,SX1,NX2])]"};

  // const std::vector<std::string> smarts_list = {
  // "[$([NX2]=[NX2+]=[NX1-]),$([NX2]=[NX2+]=N),$([NX2-]-[NX2+]#[NX1]),$([NX3](=[OX1])=[OX1]),$([NX3+](=[OX1])[O-]),$([NX2+]#[CX1-]),$([NX2]#[CX1]),$([OX1-,OH1][#7X4+]([*])([*])([*])),$([OX1]=[#7X4v5]([*])([*])([*])),$([OX1-,OH1][#7X3+R](~[R])(~[R])),$([OX1]=[#7v5X3R](~[R])(~[R])),$([*+]~[*-])]"};
  // const std::vector<std::string> smarts_list = {
  //     "[$([N&X2]=[N&X2&+:2]=[N&X1&-:3])]",
  //     "[$([NX2]=[NX2+]=[NX1-]),$([NX2]=[NX2+]=N)]",
  //     "[$([NX2]=[NX2+]=[NX1-]),$([NX2]=[NX2+]=N),$([NX2-]-[NX2+]#[NX1]),$([NX3](=[OX1])=[OX1])]"};
  // const std::vector<std::string> smarts_list = {
  // "[$([NX2]=[NX2+]=[NX1-]),$([NX2]=[NX2+]=N),$([NX2-]-[NX2+]#[NX1]),$([NX3](=[OX1])-[OX1]),$([NX3+](=[OX1])[O-]),$([NX2+]#[CX1-]),$([NX2]#[CX1]),$([OX1-,OH1][#7X4+]([*])([*])([*])),$([OX1]=[#7X4v5]([*])([*])([*])),$([OX1-,OH1][#7X3+R](~[R])(~[R])),$([OX1]=[#7v5X3R](~[R])(~[R])),$([*+]~[*-])]"
  // ,
  // "$([OX1-,OH1][#7X4+]([*])([*])([*])),$([OX1]=[#7X4v5]([*])([*])([*])),$([OX1-,OH1][#7X3+R](~[R])(~[R])),$([OX1]=[#7v5X3R](~[R])(~[R])),$([*+]~[*-])]"
  // , $([NX2]=[NX2+]=[NX1-]),$([NX2]=[NX2+]=N),
  // "[$([NX2]=[NX2+]=[NX1-]),$([NX2]=[NX2+]=N),$([NX2-]-[NX2+]#[NX1]),$([NX3](=[OX1])=[OX1]),$([NX3+](=[OX1])[O-]),$([NX2+]#[CX1-]),$([NX2]#[CX1]),$([OX1-,OH1][#7X4+]([*])([*])([*])),$([OX1]=[#7X4v5]([*])([*])([*])),$([OX1-,OH1][#7X3+R](~[R])(~[R])),$([OX1]=[#7v5X3R](~[R])(~[R])),$([*+]~[*-])]"
  // };

  const std::vector<std::string> smarts_list = {
      // "[#6&!$(C#N)][CX3](=[OX1])",
      // "[H][CX3]([#1,#6])=[OX1]",
      // "aF",
      // "[#6][Sv6X4](=[O])(=[O])[OH1,O-]",
      // // "[$([NX3](=[OX1])=[OX1])]",
      // "[c;H0;D3;+0]-[c;H0;D3;+0]",
      // "[#16]-[N&H2&D1&+0]",
      // "[$([!#1!#6]1~[#6]=[#6]~1),$([!#1!#6]1=[#6]~[#6]~1)]",
      // // "[#6;$([#6]1=[#6]~[#6]~1)]",
      // // "[c]:[c&H1&D2&+0]:[c]"
      // "[c]:[n&H1&D2&+0]:[c]",
      // "[$([!#1!#6]1~[#6]=[#6]~1)]",
      // "[!#1!#6]1=[#6]~[#6]~1",
      // "[!#1!#6]1~[#6]=[#6]~1",
      "[NR1]1=[CR1][CR1][CR1][CR1]1"
      // "[SR1]1[SR1][CR1]=[CR1][CR1]1",
      // "[C]/[C&H1&D2&+0]=[C&H1&D2&+0]\\[C]",
      // "O/[N&H0&D2&+0]=[C&H0&D3&+0](\[N&H2&D1&+0])-[c&H0&D3&+0](:[c]):[c]"
      // "[#6&!$([CX4]([SX2])([#7,O,S,F,Cl,Br,I,P]))&!$([CX3]([SX2])=[OX1,SX1,NX2,C])&!$([#6][SX2]C#N)][SX2][#6&!$([CX4]([SX2])([#7,O,S,F,Cl,Br,I,P]))&!$([CX3]([SX2])=[OX1,SX1,NX2,C])&!$([#6][SX2]C#N)]",
      // "c[OX2&!$([OX2r3])]c",
      // "c[OX2&!$([OX2r3])]c",
      // "c[OX2&!$([OX2r3])]c",
      // "c[OX2&!$([OX2r3])]c",
      // "c[OX2&!$([OX2r3])]c",
      // "c[OX2&!$([OX2r3])]c",
      // "c[OX2&!$([OX2r3])]c",
      // "c[OX2&!$([OX2r3])]c",
      // "c[OX2&!$([OX2r3])]c",
      // "c[OX2&!$([OX2r3])]c",
      // "c[OX2&!$([OX2r3])]c",
      // "c[OX2&!$([OX2r3])]c",
      // "c[OX2&!$([OX2r3])]c",
      // "c[OX2&!$([OX2r3])]c",
      // "c[OX2&!$([OX2r3])]c",
      // "c[OX2&!$([OX2r3])]c",
      // "c[OX2&!$([OX2r3])]c",
      // "c[OX2&!$([OX2r3])]c",
      // "c[OX2&!$([OX2r3])]c",
      // "c[OX2&!$([OX2r3])]c",
      // "c[OX2&!$([OX2r3])]c",
      // "c[OX2&!$([OX2r3])]c",

      // "[#1,#6][CX3]([#1,#6])=[NX2][NX3]([#1,#6&!$([CX3]=[OX1,SX1,NX2])])[CX3](=[OX1])[NX3]([#1,#6&!$([CX3]=[OX1,SX1,NX2])])[#1,#6&!$([CX3]=[OX1,SX1,NX2])]"

  };
  atom_typer::SmartsAnalyzer sa;
  atom_typer::AtomTyper at;
  std::vector<std::string> results;
  try {
    results = sa.standard_smarts(smarts_list, false, false, false);
  } catch (const std::exception &e) {
    std::cerr << "!!! Error during SMARTS analysis: " << e.what() << std::endl;
    return 1;
  }
  std::cout << "Generated " << results.size()
            << " standardized SMARTS variants." << std::endl;
  for (int i = 0; i < results.size(); ++i) {
    const auto &s = results[i];
    const std::string mapped_smarts = sa.add_atom_maps(smarts_list[i]);

    std::cout << "input: " << mapped_smarts << std::endl;
    std::cout << "output: " << s << std::endl << std::endl;

    try {
      const auto bench =
          atom_typer::smoke_benchmark::benchmark_smarts_pair(mapped_smarts, s);
      std::cout << "benchmark dataset: " << bench.dataset_path.string()
                << std::endl;
      std::cout << "target_count=" << bench.target_count
                << " before=" << bench.initial_matches.size()
                << " after=" << bench.final_matches.size()
                << " intersection=" << bench.intersection.size()
                << " before_only=" << bench.initial_only.size()
                << " after_only=" << bench.final_only.size() << std::endl;
      std::cout << "jaccard=" << bench.jaccard
                << " precision=" << bench.precision
                << " recall=" << bench.recall << " f1=" << bench.f1 << std::endl
                << std::endl;
      std::cout << "before-only rows" << std::endl;
      for (const auto &row_smiles : bench.initial_only_rows) {
        std::cout << row_smiles.second << std::endl;
      }
    } catch (const std::exception &e) {
      std::cerr << "benchmark error: " << e.what() << std::endl;
    }
  }
  return 0;
}
