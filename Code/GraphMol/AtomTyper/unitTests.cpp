#include <GraphMol/AtomTyper/smarts_analyzer.hpp>
#include <GraphMol/AtomTyper/atom_typer.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <tuple>

int main() {
  // Each entry: { description, smarts_A, smarts_B }
  // Both should produce the same standardized output.
  const std::vector<std::tuple<std::string, std::string, std::string>> test_cases = {

      // ── 1. Atomic number vs element symbol ──
      {"Carbon: #6 vs C",
       "[#6;X4;H3]",
       "[C;X4;H3]"},

      {"Nitrogen: #7 vs N",
       "[#7;H2;D1]",
       "[N;H2;D1]"},

      {"Aromatic carbon: #6&a vs c",
       "[#6&a]",
       "[c]"},

      {"Aromatic nitrogen: #7&a vs n",
       "[#7&a]",
       "[n]"},

      {"Sulfur: #16 vs S",
       "[#16;X2;H1]",
       "[S;X2;H1]"},

      // ── 2. Operator reordering (AND high vs AND low) ──
      {"AND operator order: H before D vs D before H",
       "[C;H0;D3]",
       "[C;D3;H0]"},

      {"Nested AND: charge before degree",
       "[N;+1;D3;H0]",
       "[N;H0;D3;+1]"},

      // ── 3. High-precedence & vs low-precedence ; ──
      {"High-AND vs low-AND mixing",
       "[C&H0&D3;+0]",
       "[C;H0;D3;+0]"},

      {"All high-AND",
       "[N&H1&D2&+0]",
       "[N;H1;D2;+0]"},

      // ── 4. Redundant aromaticity primitives ──
      {"Aromatic carbon: c&a vs c",
       "[c&a]",
       "[c]"},

      {"Aliphatic carbon: C&A vs C",
       "[C&A]",
       "[C]"},

      {"Atomic num + aromatic flag: #6&a vs c",
       "[#6&a;H1;D2]",
       "[c;H1;D2]"},

      {"Atomic num + aliphatic flag: #6&A vs C",
       "[#6&A;H3;D1]",
       "[C;H3;D1]"},

      // ── 5. OR branches in different orders ──
      {"OR branch order: F,Cl,Br vs Br,Cl,F",
       "[C]([F,Cl,Br])",
       "[C]([Br,Cl,F])"},

      {"OR charge order",
       "[N;+0,+1]",
       "[N;+1,+0]"},

      // ── 6. Negation equivalences ──
      {"!H0 vs H1,H2,H3,H4 on carbon",
       "[C;!H0]",
       "[C;H1,H2,H3,H4]"},

      {"Ring constraint: R vs !R0 (in ring)",
       "[C;R]",
       "[C;!R0]"},

      // ── 7. Recursive SMARTS with equivalent anchors ──
      {"Recursive anchor: $(C=O) vs $(#6=O)",
       "[N;$([N]C=O)]",
       "[N;$([N][#6]=O)]"},

      {"Recursive with aromatic anchor",
       "[C;$(C-c)]",
       "[C;$(C-[#6&a])]"},

      // ── 8. Ring membership ──
      {"Not in ring: !R vs R0",
       "[C;!R;D3]",
       "[C;R0;D3]"},

      // ── 9. Connectivity X vs D equivalence (when no implicit H ambiguity) ──
      {"Explicit degree on saturated carbon",
       "[C;X4;H0;D4]",
       "[C;D4;H0;X4]"},

      // ── 10. Multi-atom patterns with bond type variations ──
      {"Bond type: explicit single vs default",
       "[C]-[C]",
       "[C][C]"},

      {"Aromatic bond: colon vs lowercase",
       "[c]:[c]",
       "c:c"},

      {"Ring bond between aromatic atoms",
       "[c]:[n]:[c]",
       "c:n:c"},

      // ── 11. Complex mixed patterns ──
      {"Carbonyl: different notations",
       "[CX3](=[OX1])[#6]",
       "[#6;X3](=[#8;X1])[#6]"},

      {"Amide nitrogen",
       "[NX3;H1][CX3](=[OX1])[#6]",
       "[N;X3;H1][C;X3](=[O;X1])[#6]"},

      {"Sulfonamide",
       "[#7][SX4](=[OX1])(=[OX1])[#6]",
       "[N][S;X4](=[O;X1])(=[O;X1])[C,c]"},

      // ── 12. Charge notation equivalences ──
      {"Positive charge: +1 vs +",
       "[N;+1;H0;D3]",
       "[N;+;H0;D3]"},

      {"Negative charge: -1 vs -",
       "[O;-1;H0;D1]",
       "[O;-;H0;D1]"},

      {"Neutral: +0 explicit",
       "[C;+0;H3;D1]",
       "[C;H3;D1;+0]"},

      // ── 13. Wildcard / any-atom equivalences ──
      {"Any atom notations",
       "[*;H0]",
       "[!#0;H0]"},

      // ── 14. Degree constraints with implicit coverage ──
      {"Degree on aromatic nitrogen",
       "[n;H0;D2]",
       "[n;D2;H0]"},

      {"Degree on aromatic sulfur",
       "[s;D2;H0]",
       "[#16&a;D2;H0]"},

      // ── 15. Patterns from real force fields ──
      {"GAFF-style sp3 carbon",
       "[C;X4;!D4]-[C;D3]",
       "[#6;A;X4;!D4]-[#6;A;D3]"},

      {"Thiol",
       "[S;X2;H1][C,c]",
       "[#16;X2;H1][#6]"},

      {"Primary amine",
       "[N;X3;H2;D1;+0]",
       "[#7;A;H2;D1;+0;X3]"},

      {"Aromatic ring nitrogen in 6-ring",
       "[n;H0;D2;$(a1aaaaa1)]",
       "[#7&a;H0;D2;$(a1aaaaa1)]"},

      {"Ether oxygen",
       "[O;X2;H0;!R][C;!R]",
       "[#8;A;X2;H0;R0][#6;A;R0]"},

      // ── 16. Hydrogen count: total vs explicit ──
      {"Total H count",
       "[C;H2;D2]",
       "[#6;A;H2;D2]"},

      // ── 17. Complex recursive with hoist ──
      {"Recursive with hoistable anchor on aromatic S",
       "[#6]:[#16;$(a1aaaa1)]",
       "[#6]:[s;$(a1aaaa1)]"},

      {"Negated recursive preserved",
       "[C;!$(C=O)]",
       "[#6;A;!$(C=O)]"},

      // ── 18. Multi-bond patterns ──
      {"Double bond carbonyl",
       "[C]=[O]",
       "[#6]=[#8]"},

      {"Triple bond nitrile",
       "[C]#[N]",
       "[#6]#[#7]"},

      // ── 19. Patterns that should NOT match each other ──
      // (Negative tests — outputs should differ)
      // These are separated so we can test both positive and negative cases.
  };

  // Negative test cases: these should produce DIFFERENT standardized output
  const std::vector<std::tuple<std::string, std::string, std::string>> negative_cases = {
      {"Aromatic vs aliphatic carbon",
       "[c]",
       "[C]"},

      {"Different charges",
       "[N;+1]",
       "[N;-1]"},

      {"Different H counts",
       "[C;H0]",
       "[C;H3]"},

      {"In ring vs not in ring",
       "[C;R]",
       "[C;!R]"},

      {"Different degree",
       "[C;D2]",
       "[C;D4]"},
  };

  atom_typer::SmartsAnalyzer sa;
  atom_typer::SmartsAnalyzer::StandardSmartsLogOptions log_options;
  atom_typer::SmartsAnalyzer::StandardSmartsWorkflowOptions workflow_options;
  workflow_options.include_x_in_reserialization = false;
  workflow_options.enumerate_bond_order = false;
  log_options.enabled = false;

  int pass = 0, fail = 0, error = 0;

  std::cout << "=== EQUIVALENCE TESTS (should produce same output) ===" << std::endl;
  for (const auto &[desc, smarts_a, smarts_b] : test_cases) {
    try {
      auto result_a = sa.standard_smarts({smarts_a}, false, false, false,
                                          workflow_options, log_options);
      auto result_b = sa.standard_smarts({smarts_b}, false, false, false,
                                          workflow_options, log_options);

      const std::string out_a = result_a.empty() ? "(empty)" : result_a[0];
      const std::string out_b = result_b.empty() ? "(empty)" : result_b[0];

      if (out_a == out_b) {
        std::cout << "  PASS: " << desc << std::endl;
        std::cout << "    A: " << smarts_a << " -> " << out_a << std::endl;
        std::cout << "    B: " << smarts_b << " -> " << out_b << std::endl;
        ++pass;
      } else {
        std::cout << "  FAIL: " << desc << std::endl;
        std::cout << "    A: " << smarts_a << " -> " << out_a << std::endl;
        std::cout << "    B: " << smarts_b << " -> " << out_b << std::endl;
        ++fail;
      }
    } catch (const std::exception &e) {
      std::cout << "  ERROR: " << desc << " — " << e.what() << std::endl;
      std::cout << "    A: " << smarts_a << std::endl;
      std::cout << "    B: " << smarts_b << std::endl;
      ++error;
    }
  }

  std::cout << "\n=== NEGATIVE TESTS (should produce different output) ===" << std::endl;
  for (const auto &[desc, smarts_a, smarts_b] : negative_cases) {
    try {
      auto result_a = sa.standard_smarts({smarts_a}, false, false, false,
                                          workflow_options, log_options);
      auto result_b = sa.standard_smarts({smarts_b}, false, false, false,
                                          workflow_options, log_options);

      const std::string out_a = result_a.empty() ? "(empty)" : result_a[0];
      const std::string out_b = result_b.empty() ? "(empty)" : result_b[0];

      if (out_a != out_b) {
        std::cout << "  PASS: " << desc << std::endl;
        std::cout << "    A: " << smarts_a << " -> " << out_a << std::endl;
        std::cout << "    B: " << smarts_b << " -> " << out_b << std::endl;
        ++pass;
      } else {
        std::cout << "  FAIL: " << desc << " (should differ but got same output)" << std::endl;
        std::cout << "    A: " << smarts_a << " -> " << out_a << std::endl;
        std::cout << "    B: " << smarts_b << " -> " << out_b << std::endl;
        ++fail;
      }
    } catch (const std::exception &e) {
      std::cout << "  ERROR: " << desc << " — " << e.what() << std::endl;
      ++error;
    }
  }

  std::cout << "\n=== SUMMARY ===" << std::endl;
  std::cout << "  Total: " << (pass + fail + error) << std::endl;
  std::cout << "  Pass:  " << pass << std::endl;
  std::cout << "  Fail:  " << fail << std::endl;
  std::cout << "  Error: " << error << std::endl;

  return (fail + error) > 0 ? 1 : 0;
}