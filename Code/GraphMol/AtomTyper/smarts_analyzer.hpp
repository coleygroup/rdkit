#ifndef SMARTS_ANALYZER_HPP
#define SMARTS_ANALYZER_HPP

#if __has_include(<RDGeneral/export.h>)
#include <RDGeneral/export.h>
#else
#include <../../RDGeneral/export.h>
#endif
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <utility>
#include <cstdint>

namespace atom_typer {

/**
 * Main SmartsAnalyzer class for operations/DOF calcs
 */
class RDKIT_ATOMTYPER_EXPORT SmartsAnalyzer {
 public:
    static constexpr unsigned int LogSummary = 1u << 0;
    static constexpr unsigned int LogMapping = 1u << 1;
    static constexpr unsigned int LogVariants = 1u << 2;
    static constexpr unsigned int LogAtomTyping = 1u << 3;
    static constexpr unsigned int LogValidation = 1u << 4;
    static constexpr unsigned int LogFinal = 1u << 5;
    static constexpr unsigned int LogErrors = 1u << 6;
    static constexpr unsigned int LogRecanon = 1u << 7;
    static constexpr unsigned int LogRecanonComparisons = 1u << 8;
    static constexpr unsigned int LogAll = 0xFFFFFFFFu;

    static constexpr unsigned int ExtractPrimitiveNone = 0u;
    static constexpr unsigned int ExtractPrimitiveH = 1u << 0;
    static constexpr unsigned int ExtractPrimitiveD = 1u << 1;
    static constexpr unsigned int ExtractPrimitiveX = 1u << 2;
    static constexpr unsigned int ExtractPrimitiveCharge = 1u << 3;
    static constexpr unsigned int ExtractPrimitiveV = 1u << 4;      // v (valence, calculated from bonds+Hs)
    static constexpr unsigned int ExtractPrimitiveRingCount = 1u << 5;  // R (ring count)
    static constexpr unsigned int ExtractPrimitiveRingSize = 1u << 6;   // r (smallest ring size)
    static constexpr unsigned int ExtractPrimitiveRingBondCount = 1u << 7; // x (ring bond count)

    struct StandardSmartsLogOptions {
        constexpr StandardSmartsLogOptions(bool enabled_in = false,
                                           unsigned int flags_in = LogAll)
            : enabled(enabled_in), flags(flags_in) {}

        bool enabled;
        unsigned int flags;
    };

    /// Controls how element symbols are encoded in the final SMARTS output.
    enum class SymbolForm {
      Unchanged = 0,  ///< no change (default)
      Expanded  = 1,  ///< C -> [#6&A], c -> [#6&a]
      Condensed = -1, ///< [#6&A] -> C, [#6&a] -> c
    };

    struct StandardSmartsWorkflowOptions {
        StandardSmartsWorkflowOptions(
            bool include_x_in_reserialization_in = false,
            bool enumerate_bond_order_in = true,
            unsigned int extracted_primitives_mask_in =
                ExtractPrimitiveNone,
            std::vector<std::string> factoring_priority_in = {})
            : include_x_in_reserialization(include_x_in_reserialization_in),
              enumerate_bond_order(enumerate_bond_order_in),
              extracted_primitives_mask(extracted_primitives_mask_in),
              factoring_priority(std::move(factoring_priority_in)) {}

        bool include_x_in_reserialization;
        bool enumerate_bond_order;
        unsigned int extracted_primitives_mask;
        std::vector<std::string> factoring_priority;

        /// Remove [A,a] OR sub-trees that cover all aromaticity states
        /// (always-true within an AND context).
        bool remove_aa_wildcard = true;

        /// Expand element symbols to explicit atomic-num + aromaticity form,
        /// or condense back to the SMARTS organic-subset shorthand.
        SymbolForm symbol_form = SymbolForm::Unchanged;

        /// Collapse OR nodes that contain only a single arm into that arm.
        bool fold_singleton_or = false;

        /// Convert abbreviated charge notation to explicit integer form:
        /// + -> +1, - -> -1, ++ -> +2, -- -> -2, etc.
        bool explicit_charge_values = false;

        /// Rewrite certain OR-joined primitive sets to equivalent negated
        /// primitives on selected atom types during normalize_smarts_encoding.
        /// Current rule set includes: H1,H2,H3,H4 -> !H0.
        bool rewrite_or_primitives_to_negated = false;

        /// Atomic numbers for which OR->negated primitive rewrites are enabled
        /// (e.g. {6} enables carbon-only behavior).
        std::vector<int> or_primitive_rewrite_atomic_nums;
    };

  SmartsAnalyzer();   // Constructor
  ~SmartsAnalyzer();  // Destructor

  /**
   * Calculate the degrees of freedom of a SMARTS string
   * @param smarts The input SMARTS string
   * @return Degrees of freedom (the number of possible variants based on
   * SMARTS)
   */
  int calculate_dof(const std::string &smarts);

  /**
   * Find the vector of all possible SMARTS string variants, with a max
   * @param smarts The input SMARTS string
   * @param max Max number of SMARTS string variants in vector return
   * @return Vector of SMARTS string variants
   */
  std::vector<std::string> enumerate_variants(const std::string &,
                                              int max = 1000,
                                              bool verbose = false,
                                              bool carry_atom_maps = false,
                                              bool enumerate_bond_order = true);

  /**
   * Add atom map numbers to each atom in the SMARTS, including atoms inside
   * recursive SMARTS expressions ($(...)).
   * @param smarts Input SMARTS
   * @param start_map First map index to assign (default 1)
   * @return SMARTS with atom maps assigned
   */
  std::string add_atom_maps(const std::string &smarts,
                            unsigned int start_map = 1);

  std::vector<std::vector<std::string>> generate_all_combos(
      std::vector<std::string> smarts_list, bool verbose,
      bool ignoreValence = false, bool catchErrors = true,
      const StandardSmartsWorkflowOptions &workflow_options =
          StandardSmartsWorkflowOptions(),
      const StandardSmartsLogOptions &log_options = StandardSmartsLogOptions());

  std::vector<std::vector<std::string>> generate_all_combos(
      std::vector<std::string> smarts_list, bool verbose,
      bool include_x_in_reserialization, bool ignoreValence, bool catchErrors,
      const StandardSmartsLogOptions &log_options = StandardSmartsLogOptions());

  std::vector<std::string> standard_smarts(
      const std::vector<std::string> &smarts_list, bool verbose,
      bool ignoreValence = false, bool catchErrors = true,
      const StandardSmartsWorkflowOptions &workflow_options =
          StandardSmartsWorkflowOptions(),
      const StandardSmartsLogOptions &log_options = StandardSmartsLogOptions());

  std::vector<std::string> standard_smarts(
      const std::vector<std::string> &smarts_list, bool verbose,
      bool include_x_in_reserialization, bool ignoreValence, bool catchErrors,
      const StandardSmartsLogOptions &log_options = StandardSmartsLogOptions());

 private:
  class Impl;
  std::unique_ptr<Impl> pimpl;
};

}  // namespace atom_typer

#endif  // SMARTS_ANALYZER_HPP