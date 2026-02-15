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
                                              bool carry_atom_maps = false);

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
      std::vector<std::string> smarts_list, bool verbose);

  std::vector<std::string> standard_smarts(
      const std::vector<std::string> &smarts_list, bool verbose);

 private:
  class Impl;
  std::unique_ptr<Impl> pimpl;
};

}  // namespace atom_typer

#endif  // SMARTS_ANALYZER_HPP