#ifndef ATOM_TYPER_QUERY_REORDER_HPP
#define ATOM_TYPER_QUERY_REORDER_HPP

#include <map>
#include <string>

namespace atom_typer::query_reorder {

void set_comparison_trace_enabled(bool enabled);

std::string reorder_query_tree_by_embedding_smarts(
    const std::string &smarts,
    const std::map<std::string, double> &embedding);

}  // namespace atom_typer::query_reorder

#endif  // ATOM_TYPER_QUERY_REORDER_HPP
