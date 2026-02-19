#ifndef ATOM_TYPER_HPP
#define ATOM_TYPER_HPP

#if __has_include(<RDGeneral/export.h>)
#include <RDGeneral/export.h>
#else
#include <../../RDGeneral/export.h>
#endif
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <optional>

namespace atom_typer {

/**
 * Structure to hold bond type information
 */
struct BondType {
  int bond_idx;                // Bond index in molecule
  int begin_atom_idx;          // Begin atom index
  int end_atom_idx;            // End atom index
  int bond_type_enum;          // RDKit Bond::BondType enum as int
  std::string bond_type_name;  // SINGLE/DOUBLE/TRIPLE/AROMATIC/etc
  std::string smarts_pattern;  // SMARTS token for this bond
  bool is_aromatic;            // Aromatic bond flag
  bool is_in_ring;             // Ring bond flag
};

/**
 * Structure to hold atom type information
 */
struct AtomType {
  int atom_idx;                // Index of the atom in the molecule
  int atomic_number;           // Atomic number (element)
  int formal_charge;           // Formal charge
  int num_hydrogens;           // Number of attached hydrogens
  int min_bonds;               // Minimum explicit bonds (SMARTS D token)
  int max_valence;             // Total valence allowed for this atom type
  bool is_aromatic;            // Aromaticity flag
  bool is_aliphatic;           // Aliphatic flag
  bool is_in_ring;             // Ring membership
  int ring_size;               // Size of smallest ring (0 if not in ring)
  std::string hybridization;   // Hybridization (SP, SP2, SP3, etc.)
  std::vector<int> neighbors;  // Indices of neighboring atoms
  std::string smarts_pattern;  // SMARTS pattern for this atom type
  std::vector<int> ring_membership_list;   // List of members in ring
  int num_ring_bonds;                      // number of bonds in ring
  int num_aliphatic_rings;                 // number of aliphatic rings
  int num_aromatic_rings;                  // number of aromatic rings
  std::string chirality;                   // Chirality information
  int ring_connectivity;                   // Ring connectivity
  int atom_type_enumeration;               // RDKit atom type encoding
  int num_single_bonds;                    // Number of single bonds
  int num_double_bonds;                    // Number of double bonds
  int num_triple_bonds;                    // Number of triple bonds
  int num_aromatic_bonds;                  // Number of aromatic bonds
  int remaining_valence;                   // Total valence - explicit valence
  std::optional<int> explicit_atomic_num;  // explicit atomic-number query
  std::optional<int> explicit_charge;      // explicit charge query
  std::optional<int> explicit_H;           // explicit H-count query
  std::optional<int> explicit_lower_h;     // explicit implicit-H (h) query
  std::optional<int> explicit_D;           // explicit degree (D) query
  std::optional<int> explicit_X;           // explicit total degree (X) query
  std::optional<int> explicit_valence;     // explicit valence query
  std::optional<int> explicit_hybridization;  // explicit hybridization query
  std::optional<bool> explicit_aromatic;      // explicit aromatic query
  std::optional<bool> explicit_aliphatic;     // explicit aliphatic query
  std::optional<bool> inferred_aromatic;      // inferred aromaticity flag
  std::optional<bool> inferred_aliphatic;     // inferred aliphaticity flag
  std::optional<bool> explicit_in_ring;       // explicit ring-membership query
  std::optional<int> explicit_ring_count;      // explicit ring-count query (R)
  std::optional<int> explicit_min_ring_size;  // explicit ring-size query
  std::vector<int> excluded_charges;          // negated charge primitives
  std::vector<int> excluded_h_counts;         // negated H-count primitives
  std::vector<int> excluded_D_values;         // negated explicit-degree values
  std::vector<int> excluded_X_values;         // negated total-degree values
  bool excluded_aromatic = false;             // negated aromatic primitive (!a)
  bool excluded_aliphatic = false;     // negated aliphatic primitive (!A)
  std::string source_atom_smarts;      // original atom SMARTS token
  std::vector<BondType> bond_details;  // Detailed typed bonds for atom
  std::map<int, int> bond_types;       // Map of bond types to counts
};

enum class PatternItemKind {
  Atom,
  Bond
};

struct PatternItem {
  PatternItemKind kind;
  AtomType atom;
  BondType bond;
};

enum class DebugLevel {
  Off = 0,
  Basic = 1,
  Verbose = 2,
  Trace = 3
};

/**
 * Main AtomTyper class for typing atoms in SMILES/SMARTS strings
 */
class RDKIT_ATOMTYPER_EXPORT AtomTyper {
 public:
  /**
   * Constructor
   */
  AtomTyper();

  /**
   * Destructor
   */
  ~AtomTyper();

  /**
   * Type a SMILES string into an ordered sequence of atoms and bonds
   * (e.g. C=C => Atom0, Bond0, Atom1).
   * @param smiles The input SMILES string
   * @return Vector of PatternItem structures for atoms/bonds in index order
   */
  std::vector<PatternItem> type_atoms_from_smiles(const std::string &smiles);

  // /**
  //  * Type all atoms in a SMARTS string
  //  * @param smarts The input SMARTS string
  //  * @return Vector of AtomType structures for each atom
  //  */
  // std::vector<AtomType> type_atoms_from_smarts(const std::string &smarts);

  /**
   * Enumerate SMARTS atom local degrees-of-freedom (H, charge, D, X)
   * while preserving the original atom/bond order.
   */
  std::string type_atoms_from_smarts(const std::string &smarts);

  /**
   * Enumerate DOF SMARTS and optionally map newly introduced atoms.
   * When enabled, only currently unmapped atoms in the rewritten SMARTS
   * (including recursive SMARTS expressions) receive map numbers starting
   * at the next map index after the maximum existing map.
   */
  std::string type_atoms_from_smarts(const std::string &smarts,
                                   bool map_new_atoms, int &max_amap,
                                   bool verbose=false,
                                   bool include_x_in_reserialization=false,
                                   bool enumerate_bond_order=true,
                                   unsigned int forced_primitives_mask=0u);

  /**
   * Overload: enumerate DOF SMARTS using custom default ranges for H and
   * formal charge when those constraints are not explicit in the query.
   */
  std::string type_atoms_from_smarts(const std::string &smarts, int h_min,
                                   int h_max, int charge_min, int charge_max,
                                   bool map_new_atoms, int &max_amap,
                                   bool include_x_in_reserialization=false,
                                   bool enumerate_bond_order=true,
                                   unsigned int forced_primitives_mask=0u);

  /**
   * Overload: enumerate with custom ranges plus runtime verbosity controls.
   * If verbose=false, debug logging is disabled regardless of debug_level.
   */
  std::string type_atoms_from_smarts(const std::string &smarts, int h_min,
                                   int h_max, int charge_min, int charge_max,
                                   bool verbose, DebugLevel debug_level,
                                   bool map_new_atoms, int &max_amap,
                                   bool include_x_in_reserialization=false,
                                   bool enumerate_bond_order=true,
                                   unsigned int forced_primitives_mask=0u);
  // Overload: enumerate with custom debug level but default H/charge ranges and
  // verbosity=true
  std::string type_atoms_from_smarts(const std::string &smarts,
                                   DebugLevel debug_level);

//   std::string enumerate_dof_smarts(const std::string &smarts,
//                                    bool map_new_atoms, int &max_amap);

  /**
   * Reorder each atom-query tree in a SMARTS expression using an embedding
   * of primitive-token scores.
   *
   * Ordering rule:
   * - `AtomAnd`: rarer branches first (ascending score)
   * - `AtomOr`: rarer branches last (descending score)
   *
   * The embedding key should be the primitive SMARTS token without brackets
   * (e.g. `#6`, `H1`, `+0`, `!H0`).
   */
  std::string reorder_query_tree_by_embedding(
      const std::string &smarts,
      const std::map<std::string, double> &embedding);

  /**
   * Return a built-in embedding for query-tree reordering.
   * Lower values represent rarer/more selective primitives.
   */
  std::map<std::string, double> get_default_query_embedding() const;

  /**
   * Validate a (fully typed) SMARTS pattern using the same valence logic as
   * internal DOF enumeration.
   */
  bool is_valid_valence_smarts(const std::string &smarts,
                               bool verbose = false) const;

  /**
   * Inspect whether a SMARTS contains a known tautomeric motif that should
   * be tolerated by valence screening.
   */
  bool inspect_tautomer(const std::string &smarts,
                        bool verbose = false) const;

  /**
   * Filter input SMARTS patterns, keeping only valence-valid patterns.
   */
  std::vector<std::string> filter_invalid_valence_smarts(
      const std::vector<std::string> &patterns) const;

    /**
     * Consolidate a set of mapped SMARTS into one expression per token-length
     * bucket by merging atom alternatives at matching atom-map numbers.
     *
     * Grouping key: number of mapped atom tokens in each input SMARTS.
     */
    std::vector<std::string> consolidate_smarts_by_atom_maps(
      const std::vector<std::string> &patterns) const;

      /**
       * Re-canonicalize each input SMARTS by reordering atom-query trees using an
       * embedding, then consolidate by atom-map index.
       *
       * If `embedding` is empty, the built-in default embedding is used.
       */
      std::vector<std::string> consolidate_smarts_by_atom_maps_recanon(
        const std::vector<std::string> &patterns,
        const std::map<std::string, double> &embedding = {}) const;

      /**
       * Consolidate mapped SMARTS of different lengths into a single SMARTS by
       * embedding divergent mapped-atom pathways as recursive alternatives at the
       * divergent atom-map position.
       */
      std::string consolidate_smarts_with_recursive_paths(
        const std::vector<std::string> &patterns) const;

  /**
   * Type SMARTS into an ordered sequence of atoms and bonds as they appear
   * in the pattern path (e.g. C=C => Atom0, Bond0, Atom1).
   */
  std::vector<PatternItem> type_pattern_from_smarts(const std::string &smarts);

  /**
   * Get atom types as a formatted string
   * @param atom_types Vector of atom types
   * @return Formatted string representation
   */
  std::string get_atom_types_string(const std::vector<AtomType> &atom_types);

  std::string get_smarts_from_pattern_types(
      const std::vector<PatternItem> &items);

  /**
   * Get mixed atom/bond sequence as a formatted string
   */
  std::string get_pattern_types_string(const std::vector<PatternItem> &items);

  /**
   * Set whether to use canonical SMILES
   * @param use_canonical Flag to enable/disable canonicalization
   */
  void set_use_canonical(bool use_canonical);

  /**
   * Set default debug level used by enumeration workflows.
   */
  void set_debug_level(DebugLevel level);

  /**
   * Enable/disable detailed query-reorder comparison tracing used by
   * recanonicalization workflows.
   */
  void set_recanon_comparison_trace(bool enabled);

 private:
  class Impl;
  std::unique_ptr<Impl> pimpl;
};

}  // namespace atom_typer

#endif  // ATOM_TYPER_HPP
