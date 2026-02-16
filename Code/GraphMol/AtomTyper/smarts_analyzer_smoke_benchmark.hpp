#pragma once

#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/MolOps.h>
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace atom_typer::smoke_benchmark {

struct TargetMol {
  size_t row_idx = 0;
  std::string smiles;
  std::unique_ptr<RDKit::ROMol> mol;
};

struct SubstructBenchmarkResult {
  std::filesystem::path dataset_path;
  size_t target_count = 0;

  std::set<size_t> initial_matches;
  std::set<size_t> final_matches;
  std::set<size_t> intersection;
  std::set<size_t> initial_only;
  std::set<size_t> final_only;

  std::vector<std::pair<size_t, std::string>> initial_only_rows;
  std::vector<std::pair<size_t, std::string>> final_only_rows;

  double jaccard = 0.0;
  double precision = 0.0;
  double recall = 0.0;
  double f1 = 0.0;
};

inline std::string trim(std::string s) {
  auto not_space = [](unsigned char c) { return !std::isspace(c); };
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
  s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
  return s;
}

inline std::string to_lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return s;
}

inline std::vector<std::string> split_csv_line(const std::string &line) {
  std::vector<std::string> fields;
  std::string cur;
  bool in_quotes = false;

  for (size_t i = 0; i < line.size(); ++i) {
    const char c = line[i];
    if (c == '"') {
      if (in_quotes && i + 1 < line.size() && line[i + 1] == '"') {
        cur.push_back('"');
        ++i;
      } else {
        in_quotes = !in_quotes;
      }
      continue;
    }

    if (c == ',' && !in_quotes) {
      fields.push_back(cur);
      cur.clear();
      continue;
    }

    cur.push_back(c);
  }

  fields.push_back(cur);
  return fields;
}

inline std::filesystem::path resolve_default_dataset_path() {
  const std::vector<std::filesystem::path> candidates = {
      std::filesystem::path(
          "C:/Projects/atype/rdkit/Code/GraphMol/AtomTyper/all_substruct_data.csv"),
      std::filesystem::path("Code/GraphMol/AtomTyper/all_substruct_data.csv"),
      std::filesystem::path(
          "rdkit/Code/GraphMol/AtomTyper/all_substruct_data.csv")};

  for (const auto &cand : candidates) {
    if (std::filesystem::exists(cand)) {
      return cand;
    }
  }

  return candidates.front();
}

inline std::vector<TargetMol> load_targets_from_csv(
    const std::filesystem::path &path) {
  std::ifstream in(path);
  if (!in.good()) {
    throw std::runtime_error("Failed to open CSV file: " + path.string());
  }

  std::vector<TargetMol> out;
  std::string line;
  size_t line_no = 0;
  int smiles_col = 0;
  bool header_checked = false;

  while (std::getline(in, line)) {
    ++line_no;
    if (line.empty()) {
      continue;
    }

    auto fields = split_csv_line(line);
    for (auto &f : fields) {
      f = trim(f);
    }

    if (!header_checked) {
      header_checked = true;
      bool has_header = false;
      for (size_t i = 0; i < fields.size(); ++i) {
        const std::string low = to_lower(fields[i]);
        if (low.find("random_substrate_smiles") != std::string::npos) {
          smiles_col = static_cast<int>(i);
          has_header = true;
          break;
        }
      }
      if (has_header) {
        continue;
      }
    }

    if (smiles_col < 0 || static_cast<size_t>(smiles_col) >= fields.size()) {
      continue;
    }

    const std::string smiles = trim(fields[smiles_col]);
    if (smiles.empty()) {
      continue;
    }

    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(smiles, 0, false));
    try{
        mol->updatePropertyCache(false);
        RDKit::MolOps::findSSSR(*mol);
    } catch (...) {
        continue;
    }
    if (!mol) {
      continue;
    }

    TargetMol tm;
    tm.row_idx = line_no;
    tm.smiles = smiles;
    tm.mol = std::move(mol);
    out.push_back(std::move(tm));
  }

  return out;
}

inline std::set<size_t> collect_match_rows(
    const std::string &query_smarts, const std::vector<TargetMol> &targets) {
  std::set<size_t> rows;
  std::unique_ptr<RDKit::ROMol> query(RDKit::SmartsToMol(query_smarts));
  if (!query) {
    throw std::runtime_error("Invalid SMARTS query: " + query_smarts);
  }

  for (const auto &t : targets) {
    if (t.mol && !RDKit::SubstructMatch(*t.mol, *query).empty()) {
      rows.insert(t.row_idx);
    }
  }

  return rows;
}

inline double ratio(size_t num, size_t den) {
  return den == 0 ? 0.0 : static_cast<double>(num) / static_cast<double>(den);
}

inline SubstructBenchmarkResult benchmark_smarts_pair(
    const std::string &initial_smarts, const std::string &final_smarts,
    const std::filesystem::path &dataset_path = std::filesystem::path()) {
  if (initial_smarts.empty() || final_smarts.empty()) {
    throw std::runtime_error("Both initial and final SMARTS must be non-empty");
  }

  const std::filesystem::path resolved_dataset =
      dataset_path.empty() ? resolve_default_dataset_path() : dataset_path;

  const auto targets = load_targets_from_csv(resolved_dataset);

  SubstructBenchmarkResult out;
  out.dataset_path = resolved_dataset;
  out.target_count = targets.size();

  out.initial_matches = collect_match_rows(initial_smarts, targets);
  out.final_matches = collect_match_rows(final_smarts, targets);

  std::set_intersection(
      out.initial_matches.begin(), out.initial_matches.end(),
      out.final_matches.begin(), out.final_matches.end(),
      std::inserter(out.intersection, out.intersection.begin()));

  std::set_difference(
      out.initial_matches.begin(), out.initial_matches.end(),
      out.final_matches.begin(), out.final_matches.end(),
      std::inserter(out.initial_only, out.initial_only.begin()));

  std::set_difference(out.final_matches.begin(), out.final_matches.end(),
                      out.initial_matches.begin(), out.initial_matches.end(),
                      std::inserter(out.final_only, out.final_only.begin()));

  std::set<size_t> union_set;
  std::set_union(out.initial_matches.begin(), out.initial_matches.end(),
                 out.final_matches.begin(), out.final_matches.end(),
                 std::inserter(union_set, union_set.begin()));

  out.jaccard = ratio(out.intersection.size(), union_set.size());
  out.precision = ratio(out.intersection.size(), out.final_matches.size());
  out.recall = ratio(out.intersection.size(), out.initial_matches.size());
  out.f1 =
      (out.precision + out.recall == 0.0)
          ? 0.0
          : (2.0 * out.precision * out.recall) / (out.precision + out.recall);

  std::unordered_map<size_t, std::string> row_to_smiles;
  row_to_smiles.reserve(targets.size());
  for (const auto &t : targets) {
    row_to_smiles.emplace(t.row_idx, t.smiles);
  }

  for (const auto row : out.initial_only) {
    const auto it = row_to_smiles.find(row);
    out.initial_only_rows.emplace_back(row, (it == row_to_smiles.end())
                                                ? std::string("<unknown>")
                                                : it->second);
  }
  for (const auto row : out.final_only) {
    const auto it = row_to_smiles.find(row);
    out.final_only_rows.emplace_back(row, (it == row_to_smiles.end())
                                              ? std::string("<unknown>")
                                              : it->second);
  }

  return out;
}

}  // namespace atom_typer::smoke_benchmark
