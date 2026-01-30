#include "smarts_analyzer.hpp" //header file for smarts_analyzer.cpp
#include <GraphMol/GraphMol.h> // Core RDKit molecule functionality
#include <GraphMol/MolOps.h> // For molecule operations from RDKit molecule
#include <GraphMol/SmilesParse/SmilesParse.h> // For parsing SMILES strings
#include <GraphMol/SmilesParse/SmartsWrite.h> // For writing SMARTS strings
#include <GraphMol/Descriptors/MolDescriptors.h> // For molecular descriptors like the AND/OR/NOT queries
#include <GraphMol/QueryAtom.h>  // For QueryAtom and Query functionality
#include <GraphMol/QueryOps.h>   // For ATOM_EQUALS_QUERY and other query types
#include <Query/QueryObjects.h> // For the base Query class
#include <Query/EqualityQuery.h> // For Equality Query functionality
#include <sstream> // For string stream operations, work with data types better
#include <stdexcept> // For std exceptions
#include <algorithm> // For std algorithms like std::max, std::min, std::sort, std::unique

//stop variable name clashing with other libraries
namespace atom_typer {

/**
 * Private implementation class PIMPL, from the header file
 * This class contains the actual implementation details for SmartsAnalyzer
 */
class SmartsAnalyzer::Impl{
public:

    // Returns number of logical alternatives encoded by this query (for single operator)
    // So, OR(x, y, z) returns 3, AND(x, y) returns 1, NOT(x) returns 1
    int count_query_dof(const Queries::Query<int, const RDKit::Atom*, true>* query) {
        
        // no query to consider, only one logical option
        if (!query) return 1;

        // retrieve name of query, such as "AtomAnd", "AtomOr", "AtomNot", or specific atom property
        const std::string& desc = query->getDescription();

        // OR node, go through number of branches from this node
        // RDKit uses "AtomOr" for atom-level OR queries
        if (desc == "AtomOr" || desc == "Or") {
            int sum = 0;

            // each option under OR node adds combinatorially, sums up via iteration/recursion
            for (auto it = query->beginChildren(); it != query->endChildren(); ++it) {

                // calls count_query_dof recursively on each child, accumulating total
                sum += count_query_dof(it->get());
            }

            return sum;
        }

        // AND node: all constraints must be satisfied (no branching at this level)
        // However, children may still have OR nodes that create branches
        // RDKit uses "AtomAnd" for atom-level AND queries
        if (desc == "AtomAnd" || desc == "And") {

            // start at 1, since if no children, AND is trivially satisfied
            // multiplies if has OR nodes in children
            int product = 1;
            
            // Recursively count DOF in children and multiply
            // Example: C&(H1,H2) means Carbon AND (1H OR 2H) = 2 alternatives
            for (auto it = query->beginChildren(); it != query->endChildren(); ++it) {
                product *= count_query_dof(it->get());
            }
            
            return product;
        }
        
        // NOT node: negation doesn't branch (simplified - just counts as 1)
        // RDKit uses "AtomNot" for atom-level NOT queries
        if (desc == "AtomNot" || desc == "Not") {
            return 1;
        }

        // If it's not a logical operator, it's a leaf node, one option
        return 1;
    }

    //class of atom constraints
    //TODO: expand constraints in a case by case basis for different atom properties
    struct AtomConstraints {

        // number of valence electrons (can calculate bond connections from this)
        int min_valence = 0;   
        int max_valence = 8;  // octet rule cap (temp)

        int min_h = 0;  // no hydrogens
        int max_h = 4;  // max hydrogens (temp)

        // formal charge (case by case for organic molecules)
        int min_charge = -8; 
        int max_charge = 8;

        // flag for contradiction, if true, no valid atom can satisfy constraints
        bool contradiction = false;
        
        // Track excluded values for NOT queries
        std::vector<int> excluded_atomic_nums;
        std::vector<int> excluded_charges;
        std::vector<int> excluded_h_counts;
        
        // Check if constraints are satisfied
        bool is_satisfied() const {
            return !contradiction;
        }
        
        // Deduplicate exclusion lists
        // makes sure there are no duplicates in the exclusion vectors and sorts them
        void deduplicate_exclusions() {
            std::sort(excluded_atomic_nums.begin(), excluded_atomic_nums.end());
            excluded_atomic_nums.erase(std::unique(excluded_atomic_nums.begin(), 
                                                   excluded_atomic_nums.end()),
                                      excluded_atomic_nums.end());
            
            std::sort(excluded_charges.begin(), excluded_charges.end());
            excluded_charges.erase(std::unique(excluded_charges.begin(),
                                              excluded_charges.end()),
                                  excluded_charges.end());
            
            std::sort(excluded_h_counts.begin(), excluded_h_counts.end());
            excluded_h_counts.erase(std::unique(excluded_h_counts.begin(),
                                               excluded_h_counts.end()),
                                   excluded_h_counts.end());
        }
    };

    // Check rules and contradictions for AND
    // returns True if there is a contradiction, and vice versa
    // merges list of child constraints into one and sets limits accordingly
    AtomConstraints merge_and(const std::vector<AtomConstraints>& children) {
        AtomConstraints result;

        // loop through the list of child constraints
        for (const auto& c : children) {

            // Make sure constraints are tightened
            result.min_valence = std::max(result.min_valence, c.min_valence);
            result.max_valence = std::min(result.max_valence, c.max_valence);

            result.min_h = std::max(result.min_h, c.min_h);
            result.max_h = std::min(result.max_h, c.max_h);

            result.min_charge = std::max(result.min_charge, c.min_charge);
            result.max_charge = std::min(result.max_charge, c.max_charge);
            
            // Merge exclusion lists (union of all exclusions)
            result.excluded_atomic_nums.insert(result.excluded_atomic_nums.end(),
                                               c.excluded_atomic_nums.begin(),
                                               c.excluded_atomic_nums.end());
            result.excluded_charges.insert(result.excluded_charges.end(),
                                          c.excluded_charges.begin(),
                                          c.excluded_charges.end());
            result.excluded_h_counts.insert(result.excluded_h_counts.end(),
                                           c.excluded_h_counts.begin(),
                                           c.excluded_h_counts.end());
        }

        // Contradicting constraints
        if (result.min_valence > result.max_valence ||
            result.min_h > result.max_h ||
            result.min_charge > result.max_charge) {
            result.contradiction = true;
            return result;
        }

        // Note: Element-specific valence limits are already set by apply_element_valence()
        // when the atomic number is known from leaf queries. No additional enforcement needed here.
        
        // Deduplicate exclusion lists
        result.deduplicate_exclusions();
        
        return result;
    }

    // Check rules for OR
    AtomConstraints merge_or(const std::vector<AtomConstraints>& children) {

        // OR takes the union (most permissive) of valid branches
        AtomConstraints result;
        result.contradiction = true;  // Start assuming failure
        
        // Loop through children to find valid branches
        for (const auto& c : children) {

            // Only consider non-contradictory branches
            if (!c.contradiction) {

                // First valid branch sets initial constraints
                if (result.contradiction) {
                    result = c;  // First valid branch
                    result.contradiction = false;
                } else {
                    // Expand constraints to union
                    result.min_valence = std::min(result.min_valence, c.min_valence);
                    result.max_valence = std::max(result.max_valence, c.max_valence);
                    result.min_h = std::min(result.min_h, c.min_h);
                    result.max_h = std::max(result.max_h, c.max_h);
                    result.min_charge = std::min(result.min_charge, c.min_charge);
                    result.max_charge = std::max(result.max_charge, c.max_charge);
                }
            }
        }
        
        return result;
    }

    // Check rules for NOT - tracks excluded values
    AtomConstraints apply_not(const AtomConstraints& negated, const std::string& negated_desc) {
        AtomConstraints result;
        
        // Try to cast the negated query to extract the specific value being excluded
        // For example, "NOT AtomicNum=6" means exclude carbon
        
        // If specific values were constrained in the negated query, we exclude them
        if (negated_desc == "AtomicNum") {

            // If negated specified a specific element, exclude it
            // This would need the actual value, which we'd get from analyzing the child
            result.excluded_atomic_nums = negated.excluded_atomic_nums;
        }

        else if (negated_desc == "AtomFormalCharge") {
            // Track excluded charges
            if (negated.min_charge == negated.max_charge) {
                result.excluded_charges.push_back(negated.min_charge);
            }
        }

        else if (negated_desc == "AtomHCount" || negated_desc == "AtomTotalHCount") {
            // Track excluded H counts
            if (negated.min_h == negated.max_h) {
                result.excluded_h_counts.push_back(negated.min_h);
            }
        }
        
        // Copy any existing exclusions from the negated constraints
        // merges existing exclusions from negated constraints
        result.excluded_atomic_nums.insert(result.excluded_atomic_nums.end(),
                                           negated.excluded_atomic_nums.begin(),
                                           negated.excluded_atomic_nums.end());
        result.excluded_charges.insert(result.excluded_charges.end(),
                                       negated.excluded_charges.begin(),
                                       negated.excluded_charges.end());
        result.excluded_h_counts.insert(result.excluded_h_counts.end(),
                                        negated.excluded_h_counts.begin(),
                                        negated.excluded_h_counts.end());
        
        return result;
    }

    // analyze valence query
    // extract constraints from a query node
    // recursively traverse query tree to extract constraints 
    AtomConstraints analyze_query(const Queries::Query<int, const RDKit::Atom*, true>* q) {

        // No query means any atom, no constraints
        if (!q) return {};

        // get description of query, such as "AtomAnd", "AtomOr", "AtomNot", or specific atom property
        const auto& desc = q->getDescription();

        // recursively handle logical operators
        // RDKit uses "AtomAnd" for atom-level AND queries
        if (desc == "AtomAnd" || desc == "And") {
            std::vector<AtomConstraints> children;

            // gather constraints from each child
            for (auto it = q->beginChildren(); it != q->endChildren(); ++it) {
                children.push_back(analyze_query(it->get()));
            }

            // for each member of children, merge constraints
            return merge_and(children); 
        }

        // OR operator
        // RDKit uses "AtomOr" for atom-level OR queries
        if (desc == "AtomOr" || desc == "Or") {
            std::vector<AtomConstraints> children;

            for (auto it = q->beginChildren(); it != q->endChildren(); ++it) {
                children.push_back(analyze_query(it->get()));
            }

            // merge constraints
            return merge_or(children);
        }

        // NOT operator
        // RDKit uses "AtomNot" for atom-level NOT queries
        if (desc == "AtomNot" || desc == "Not") {
            
            // analyze negated child
            auto it = q->beginChildren();
            if (it != q->endChildren()) {
                AtomConstraints negated = analyze_query(it->get());
                return apply_not(negated, it->get()->getDescription());
            }
            return AtomConstraints{};
        }

        // Leaf queries - extract values from specific query types
        AtomConstraints c;

        // Try to cast to EqualityQuery to extract value
        const auto* eq = dynamic_cast<const Queries::EqualityQuery<int, const RDKit::Atom*, true>*>(q);
        
        // If it's an equality query, extract the value and set constraints accordingly
        if (eq) {
            int val = eq->getVal();
            
            if (desc == "AtomicNum") {

                // Set valence limits based on element
                apply_element_valence(c, val);

                // Store the atomic number so NOT queries can exclude it
                c.excluded_atomic_nums.clear();

                // Note: for positive match, we don't add to exclusions
                // Exclusions are populated by NOT operations
            }

            else if (desc == "AtomHCount" || desc == "AtomTotalHCount" || desc == "AtomImplicitHCount") {
                apply_hydrogen(c, val);
            }

            else if (desc == "AtomFormalCharge") {
                apply_charge(c, val);
            }
            else if (desc == "AtomExplicitValence" || desc == "AtomTotalValence") {
                apply_valence(c, val);
            }
            else if (desc == "AtomExplicitDegree" || desc == "AtomTotalDegree") {
                // Degree is closely related to valence
                c.min_valence = val;
                c.max_valence = val;
            }
        }
        
        return c;
    }

    // Helper to set valence based on element
    void apply_element_valence(AtomConstraints& c, int atomic_num) {

        // Set realistic valence limits based on common elements
        switch (atomic_num) {
            case 1:  // H
                c.min_valence = 0;
                c.max_valence = 1;
                break;
            case 6:  // C
                c.min_valence = 0;
                c.max_valence = 4;
                break;
            case 7:  // N
                c.min_valence = 0;
                c.max_valence = 4;  // 3 + lone pair or positively charged
                break;
            case 8:  // O
                c.min_valence = 0;
                c.max_valence = 2;
                break;
            case 9:  // F
                c.min_valence = 0;
                c.max_valence = 1;
                break;
            case 15: // P
                c.min_valence = 0;
                c.max_valence = 5;
                break;
            case 16: // S
                c.min_valence = 0;
                c.max_valence = 6;
                break;
            case 17: // Cl
                c.min_valence = 0;
                c.max_valence = 7;  // Can be hypervalent (e.g., ClO4-)
                break;
            case 35: // Br
                c.min_valence = 0;
                c.max_valence = 7;  // Can be hypervalent
                break;
            case 53: // I
                c.min_valence = 0;
                c.max_valence = 7;  // Can be hypervalent
                break;
            default:
                // For other elements, use periodic table groups as guide
                c.min_valence = 0;
                c.max_valence = 8;
                break;
        }
    }

    // Apply valence constraints
    void apply_valence(AtomConstraints& c, int val) {
        c.min_valence = val;
        c.max_valence = val;
    }

    void apply_hydrogen(AtomConstraints& c, int val) {
        c.min_h = val;
        c.max_h = val;
    }

    void apply_charge(AtomConstraints& c, int val) {
        c.min_charge = val;
        c.max_charge = val;
    }

    // Extract SMARTS string representation for a single query node
    std::string query_to_smarts(const Queries::Query<int, const RDKit::Atom*, true>* q) {
        if (!q) return "*";
        
        const std::string& desc = q->getDescription();
        
        // Handle AND nodes by combining all child properties
        if (desc == "AtomAnd" || desc == "And") {
            std::string element = "";
            std::string properties = "";
            bool is_aromatic = false;
            
            // Collect all properties from AND children
            for (auto it = q->beginChildren(); it != q->endChildren(); ++it) {
                const auto* child = it->get();
                const std::string& child_desc = child->getDescription();
                const auto* eq = dynamic_cast<const Queries::EqualityQuery<int, const RDKit::Atom*, true>*>(child);
                
                if (eq) {
                    int val = eq->getVal();
                    
                    if (child_desc == "AtomicNum") {
                        // Element symbols
                        if (val == 6) element = is_aromatic ? "c" : "C";
                        else if (val == 7) element = is_aromatic ? "n" : "N";
                        else if (val == 8) element = is_aromatic ? "o" : "O";
                        else if (val == 16) element = is_aromatic ? "s" : "S";
                        else if (val == 15) element = is_aromatic ? "p" : "P";
                        else element = "#" + std::to_string(val);
                    } else if (child_desc == "AtomType") {
                        int atomic_num = val & 0xFF;
                        is_aromatic = (val > 255);
                        if (atomic_num == 6) element = is_aromatic ? "c" : "C";
                        else if (atomic_num == 7) element = is_aromatic ? "n" : "N";
                        else if (atomic_num == 8) element = is_aromatic ? "o" : "O";
                        else if (atomic_num == 16) element = is_aromatic ? "s" : "S";
                        else if (atomic_num == 15) element = is_aromatic ? "p" : "P";
                        else element = "#" + std::to_string(atomic_num);
                    } else if (child_desc == "AtomHCount" || child_desc == "AtomTotalHCount") {
                        properties += "H" + std::to_string(val);
                    } else if (child_desc == "AtomFormalCharge") {
                        if (val > 0) properties += "+" + std::to_string(val);
                        else if (val < 0) properties += std::to_string(val);
                        else properties += "+0";
                    }
                }
            }
            
            if (element.empty()) element = "*";
            return "[" + element + properties + "]";
        }
        
        // Try to extract value from equality queries
        const auto* eq = dynamic_cast<const Queries::EqualityQuery<int, const RDKit::Atom*, true>*>(q);
        
        if (eq) {
            int val = eq->getVal();
            
            if (desc == "AtomicNum") {
                // Return element symbol in brackets
                return "[#" + std::to_string(val) + "]";
            } else if (desc == "AtomType") {
                // AtomType encodes both atomic number and aromaticity
                int atomic_num = val & 0xFF;
                bool aromatic = (val > 255);
                
                // Common elements
                if (atomic_num == 6) return aromatic ? "[c]" : "[C]";
                if (atomic_num == 7) return aromatic ? "[n]" : "[N]";
                if (atomic_num == 8) return aromatic ? "[o]" : "[O]";
                if (atomic_num == 16) return aromatic ? "[s]" : "[S]";
                if (atomic_num == 15) return aromatic ? "[p]" : "[P]";
                
                return "[#" + std::to_string(atomic_num) + "]";
            } else if (desc == "AtomHCount" || desc == "AtomTotalHCount") {
                return "[H" + std::to_string(val) + "]";
            } else if (desc == "AtomFormalCharge") {
                if (val > 0) return "[+" + std::to_string(val) + "]";
                if (val < 0) return "[" + std::to_string(val) + "]";
                return "[+0]";
            }
        }
        
        // For complex queries (NOT, etc.), return generic wildcard
        return "[*]";
    }
    
    // Enumerate SMARTS strings from a query node
    std::vector<std::string> enumerate_query_smarts(const Queries::Query<int, const RDKit::Atom*, true>* q) {
        std::vector<std::string> results;
        
        if (!q) {
            results.push_back("[*]");
            return results;
        }
        
        const std::string& desc = q->getDescription();
        
        // OR node: recursively enumerate all branches
        if (desc == "AtomOr" || desc == "Or") {
            for (auto it = q->beginChildren(); it != q->endChildren(); ++it) {
                auto child_variants = enumerate_query_smarts(it->get());
                results.insert(results.end(), child_variants.begin(), child_variants.end());
            }
            return results;
        }
        
        // For non-OR nodes, convert to SMARTS string
        results.push_back(query_to_smarts(q));
        return results;
    }

    // Helper to check if a query value is in the exclusion list
    bool is_excluded_value(const Queries::Query<int, const RDKit::Atom*, true>* q, 
                           const AtomConstraints& constraints) {
        const auto* eq = dynamic_cast<const Queries::EqualityQuery<int, const RDKit::Atom*, true>*>(q);
        if (!eq) {
            // Not an equality query - could be complex (AND/OR/NOT)
            // For complex queries, check recursively through children
            const std::string& desc = q->getDescription();
            if (desc == "AtomAnd" || desc == "And" || desc == "AtomOr" || desc == "Or") {
                for (auto it = q->beginChildren(); it != q->endChildren(); ++it) {
                    if (is_excluded_value(it->get(), constraints)) {
                        return true;
                    }
                }
            }
            return false;
        }
        
        const std::string& desc = q->getDescription();
        int val = eq->getVal();
        
        // Check if this specific value is excluded
        if (desc == "AtomicNum") {
            return std::find(constraints.excluded_atomic_nums.begin(),
                           constraints.excluded_atomic_nums.end(), val) 
                   != constraints.excluded_atomic_nums.end();
        }
        else if (desc == "AtomFormalCharge") {
            return std::find(constraints.excluded_charges.begin(),
                           constraints.excluded_charges.end(), val)
                   != constraints.excluded_charges.end();
        }
        else if (desc == "AtomHCount" || desc == "AtomTotalHCount" || desc == "AtomImplicitHCount") {
            return std::find(constraints.excluded_h_counts.begin(),
                           constraints.excluded_h_counts.end(), val)
                   != constraints.excluded_h_counts.end();
        }
        
        return false;
    }

    void enumerate_cartesian(
        const std::vector<std::vector<const Queries::Query<int, const RDKit::Atom*, true>*>>& atom_variants,
        size_t idx,
        RDKit::ROMol* mol,
        std::vector<std::string>& out,
        int max
    ) {
        // Base case: all atoms processed
        if (out.size() >= max) return;

        // If we've assigned queries for all atoms, generate SMARTS
        if (idx == atom_variants.size()) {
            out.push_back(RDKit::MolToSmarts(*mol));
            return;
        }

        // Process current atom
        RDKit::Atom* atom = mol->getAtomWithIdx(idx);
        RDKit::QueryAtom* qatom = dynamic_cast<RDKit::QueryAtom*>(atom);
        if (!qatom) {
            // If not a query atom, skip to next
            enumerate_cartesian(atom_variants, idx + 1, mol, out, max);
            return;
        }

        // Save original query to restore later
        auto* original_query = qatom->getQuery();
        
        // If no variants for this atom, just keep the original and continue
        if (atom_variants[idx].empty()) {
            enumerate_cartesian(atom_variants, idx + 1, mol, out, max);
            return;
        }
        
        // Track if we found any valid variants
        bool found_valid = false;
        bool modified_query = false;

        // Iterate through all variants for this atom           
        for (auto* q : atom_variants[idx]) {

            // Analyze constraints for this query variant
            AtomConstraints constraints = analyze_query(q);
            
            // Skip if constraints are contradictory (chemically impossible)
            if (!constraints.is_satisfied()) {
                continue;
            }
            
            // Check if this specific query value is in the exclusion list
            if (is_excluded_value(q, constraints)) {
                continue;
            }
            
            found_valid = true;
            
            // Only modify the query if it's different from the original
            if (q != original_query) {
                // Clone the query to avoid ownership issues with sub-nodes
                auto* cloned_query = q->copy();
                qatom->setQuery(cloned_query);
                modified_query = true;
            }
            enumerate_cartesian(atom_variants, idx + 1, mol, out, max);
            if (out.size() >= max) break;
        }
        
        // If no valid variants were found, all were filtered out
        // Continue with original query to at least try
        if (!found_valid) {
            qatom->setQuery(original_query);
            enumerate_cartesian(atom_variants, idx + 1, mol, out, max);
            modified_query = false;  // We set it back to original
        }
        
        // Restore original query only if we modified it
        if (modified_query) {
            qatom->setQuery(original_query);
        }
    }
};

// SmartsAnalyzer implementation
SmartsAnalyzer::SmartsAnalyzer() : pimpl(std::make_unique<Impl>()) {}

SmartsAnalyzer::~SmartsAnalyzer() = default;

/**
 * Find the vector of all strings of SMARTS variants w/ max
 */
std::vector<std::string> SmartsAnalyzer::enumerate_variants(const std::string& smarts, int max) {
    std::vector<std::string> results;
    
    if (smarts.empty()) {
        return results;
    }
    
    if (max <= 0) {
        return results;
    }
    
    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmartsToMol(smarts));
    if (!mol) throw std::runtime_error("Invalid SMARTS");
    
    // Enumerate SMARTS alternatives for each atom
    std::vector<std::vector<std::string>> atom_smarts;
    
    for (const auto atom : mol->atoms()) {
        const RDKit::QueryAtom* qatom = dynamic_cast<const RDKit::QueryAtom*>(atom);
        
        if (qatom && qatom->getQuery()) {
            // Get SMARTS alternatives for this atom
            auto variants = pimpl->enumerate_query_smarts(qatom->getQuery());
            atom_smarts.push_back(variants);
        } else {
            // Non-query atom - shouldn't happen in valid SMARTS, but handle it
            atom_smarts.push_back({smarts});
        }
    }
    
    // Generate cartesian product of all atom SMARTS combinations
    std::function<void(size_t, std::string)> generate_combinations;
    generate_combinations = [&](size_t idx, std::string current) {
        if (results.size() >= static_cast<size_t>(max)) return;
        
        if (idx == atom_smarts.size()) {
            // Validate the generated SMARTS before adding
            try {
                std::unique_ptr<RDKit::ROMol> test_mol(RDKit::SmartsToMol(current));
                if (test_mol && test_mol->getNumAtoms() == mol->getNumAtoms()) {
                    results.push_back(current);
                }
            } catch (...) {
                // Invalid SMARTS, skip it
            }
            return;
        }
        
        // Try each alternative for this atom
        for (const auto& atom_smarts_str : atom_smarts[idx]) {
            std::string next = current + atom_smarts_str;
            
            // Add bonds if not last atom
            if (idx + 1 < atom_smarts.size()) {
                // For simplicity, assume implicit single bonds
                // In a full implementation, we'd need to handle explicit bonds from the original SMARTS
            }
            
            generate_combinations(idx + 1, next);
            
            if (results.size() >= static_cast<size_t>(max)) break;
        }
    };
    
    generate_combinations(0, "");
    
    // Deduplicate
    std::sort(results.begin(), results.end());
    results.erase(std::unique(results.begin(), results.end()), results.end());
    
    return results;
}

int SmartsAnalyzer::calculate_dof(const std::string& smarts) {
    
    // nothing in SMARTS, return 0
    if (smarts.empty()) {
        return 0;
    }

    else {
        // Parse SMARTS into RDKit molecule
        std::unique_ptr<RDKit::ROMol> mol(RDKit::SmartsToMol(smarts));
        if (!mol) {
            throw std::runtime_error("Failed to parse SMARTS: " + smarts);
        }

        //DOF calculation logic goes here
        int dof = 1; // combinatorial logic, starts at 1
        const int MAX_DOF = 1000000; // Prevent overflow

        //iterate over all the atoms in the RDKit molecule
        for (const auto atom: mol->atoms()) {
            
            const auto* qatom = dynamic_cast<const RDKit::QueryAtom*>(atom);
            if (!qatom) continue;

            int atom_dof = pimpl->count_query_dof(qatom->getQuery());

            // Check for overflow before multiplying
            if (dof > 0 && atom_dof > MAX_DOF / dof) {
                // Would overflow - return max value
                return MAX_DOF;
            }
            
            // Multiply into total DOF
            dof *= atom_dof;
        }

        return dof;
    }
}
} // namespace atom_typer