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

    struct StandardSmartsLogOptions {
        constexpr StandardSmartsLogOptions(bool enabled_in = false,
                                           unsigned int flags_in = LogAll)
            : enabled(enabled_in), flags(flags_in) {}

        bool enabled;
        unsigned int flags;
    };

    struct StandardSmartsWorkflowOptions {
        constexpr StandardSmartsWorkflowOptions(
            bool include_x_in_reserialization_in = false,
            bool enumerate_bond_order_in = true)
            : include_x_in_reserialization(include_x_in_reserialization_in),
              enumerate_bond_order(enumerate_bond_order_in) {}

        bool include_x_in_reserialization;
        bool enumerate_bond_order;
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