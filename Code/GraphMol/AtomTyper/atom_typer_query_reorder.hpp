#ifndef ATOM_TYPER_QUERY_REORDER_HPP
#define ATOM_TYPER_QUERY_REORDER_HPP

#include <map>
#include <string>
#include <vector>

namespace atom_typer::query_reorder {

void set_comparison_trace_enabled(bool enabled);
bool comparison_trace_enabled();

std::string reorder_query_tree_by_embedding_smarts(
    const std::string &smarts,
    const std::map<std::string, double> &embedding);

std::vector<std::string> reorder_patterns_by_embedding(
    const std::vector<std::string> &patterns,
    const std::map<std::string, double> &embedding);

}  // namespace atom_typer::query_reorder

#endif  // ATOM_TYPER_QUERY_REORDER_HPP
