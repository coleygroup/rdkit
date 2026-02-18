#include <GraphMol/AtomTyper/atom_typer.hpp>
#include <GraphMol/AtomTyper/smarts_analyzer.hpp>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/MolOps.h>

#include <algorithm>
#include <chrono>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <optional>
#include <random>
#include <set>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#if __has_include(<tracy/Tracy.hpp>)
#include <tracy/Tracy.hpp>
#else
#define ZoneScoped
#define ZoneScopedN(x)
#endif

namespace {

struct TargetMol {
  size_t row_idx = 0;
  std::string smiles;
  std::unique_ptr<RDKit::ROMol> mol;
};

std::string to_lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return s;
}

std::string trim(std::string s) {
  auto not_space = [](unsigned char c) { return !std::isspace(c); };
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
  s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
  return s;
}

std::vector<std::string> split_csv_line(const std::string &line) {
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

std::filesystem::path resolve_dataset_path(int argc, char **argv) {
  if (argc > 1 && argv[1] && std::string(argv[1]).size() > 0) {
    return std::filesystem::path(argv[1]);
  }

  const std::vector<std::filesystem::path> candidates = {
      std::filesystem::path(
          "C:\\Projects\\atype\\rdkit\\Code\\GraphMol\\AtomTyper\\all_substruct_data.csv"),
      std::filesystem::path("Code/GraphMol/AtomTyper/all_substruct_data.csv"),
      std::filesystem::path("Code/GraphMol/AtomTyper/sigma_smiles.csv"),
      std::filesystem::path(
          "rdkit/Code/GraphMol/AtomTyper/all_substruct_data.csv")};

  for (const auto &cand : candidates) {
    if (std::filesystem::exists(cand)) {
      std::cout << "Using dataset file: " << cand.string() << std::endl;
      return cand;
    }
  }

  return candidates.front();
}

std::filesystem::path resolve_simple_templates_path(int argc, char **argv) {
  if (argc > 1 && argv[1] && std::string(argv[1]).size() > 0) {
    return std::filesystem::path(argv[1]);
  }

  const std::vector<std::filesystem::path> candidates = {
      std::filesystem::path(
          "C:\\Projects\\atype\\rdkit\\Code\\GraphMol\\AtomTyper\\simple_smarts_list.csv"),
      std::filesystem::path("Code/GraphMol/AtomTyper/simple_smarts_list.csv"),
      std::filesystem::path(
          "rdkit/Code/GraphMol/AtomTyper/simple_smarts_list.csv")};

  for (const auto &cand : candidates) {
    if (std::filesystem::exists(cand)) {
      std::cout << "Using simple templates file: " << cand.string()
                << std::endl;
      return cand;
    }
  }

  return candidates.front();
}

std::filesystem::path resolve_templates_path(int argc, char **argv) {
  if (argc > 3 && argv[3] && std::string(argv[3]).size() > 0) {
    return std::filesystem::path(argv[3]);
  }

  const std::vector<std::filesystem::path> candidates = {
      std::filesystem::path(
          "C:\\Projects\\atype\\rdkit\\Code\\GraphMol\\AtomTyper\\noncanon_efg_templates_20240108.csv"),
      std::filesystem::path(
          "Code/GraphMol/AtomTyper/noncanon_efg_templates_20240108.csv"),
      std::filesystem::path(
          "rdkit/Code/GraphMol/AtomTyper/noncanon_efg_templates_20240108.csv")};

  for (const auto &cand : candidates) {
    if (std::filesystem::exists(cand)) {
      std::cout << "Using template file: " << cand.string() << std::endl;
      return cand;
    }
  }

  return candidates.front();
}


std::filesystem::path resolve_large_templates_path(int argc, char **argv) {
  if (argc > 3 && argv[3] && std::string(argv[3]).size() > 0) {
    return std::filesystem::path(argv[3]);
  }

  const std::vector<std::filesystem::path> candidates = {
      std::filesystem::path(
          "C:\\Projects\\atype\\rdkit\\Code\\GraphMol\\AtomTyper\\ValidationPatterns.csv"),
      std::filesystem::path(
          "Code/GraphMol/AtomTyper/ValidationPatterns.csv"),
      std::filesystem::path(
          "rdkit/Code/GraphMol/AtomTyper/ValidationPatterns.csv")};

  for (const auto &cand : candidates) {
    if (std::filesystem::exists(cand)) {
      std::cout << "Using template file: " << cand.string() << std::endl;
      return cand;
    }
  }

  return candidates.front();
}


std::vector<TargetMol> load_targets_from_csv(
    const std::filesystem::path &path, size_t target_limit = 0,
    bool random_selection = false,
    std::optional<unsigned int> random_seed = std::nullopt) {
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
          // if (low.find("smarts") != std::string::npos) {
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
    if (!mol) {
      continue;
    }
    try {
      if (mol) {
        mol->updatePropertyCache(false);
        RDKit::MolOps::findSSSR(*mol);
        for (auto &a : mol->atoms()) {
          a->setIsotope(0);
        }
      }
    } catch (...) {
      continue;
    }

    TargetMol tm;
    tm.row_idx = line_no;
    tm.smiles = smiles;
    tm.mol = std::move(mol);
    out.push_back(std::move(tm));
  }

  if (target_limit > 0 && out.size() > target_limit) {
    if (random_selection) {
      std::mt19937 rng(
          random_seed.has_value() ? *random_seed : std::random_device{}());
      std::shuffle(out.begin(), out.end(), rng);
    }
    out.resize(target_limit);
  }

  return out;
}

std::vector<std::string> load_templates_from_csv(
    const std::filesystem::path &path, int max = -1) {
  std::ifstream in(path);
  if (!in.good()) {
    throw std::runtime_error("Failed to open template CSV file: " +
                             path.string());
  }

  std::vector<std::string> out;
  std::set<std::string> seen;
  std::string line;
  bool header_checked = false;
  int template_col = 0;

  while (std::getline(in, line)) {
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
        if (low.find("noncanon_efg_templates") != std::string::npos ||
            low.find("template") != std::string::npos ||
            low.find("smarts") != std::string::npos) {
          template_col = static_cast<int>(i);
          has_header = true;
          break;
        }
      }
      if (has_header) {
        continue;
      }
    }

    if (template_col < 0 ||
        static_cast<size_t>(template_col) >= fields.size()) {
      continue;
    }

    const std::string tpl = trim(fields[template_col]);
    if (tpl.empty()) {
      continue;
    }

    if (seen.insert(tpl).second) {
      out.push_back(tpl);
    }
  }

  if (max >= 0 && static_cast<size_t>(max) < out.size()) {
    out.resize(max);
  }

  return out;
}

std::vector<std::string> filter_templates_by_max_atoms(
    const std::vector<std::string> &templates, size_t min_atoms,
    size_t &skipped_too_small, size_t &skipped_invalid) {
  std::vector<std::string> out;
  out.reserve(templates.size());
  skipped_too_small = 0;
  skipped_invalid = 0;

  for (const auto &tpl : templates) {
    std::unique_ptr<RDKit::ROMol> query(RDKit::SmartsToMol(tpl));
    if (!query) {
      ++skipped_invalid;
      continue;
    }
    if (query->getNumAtoms() > min_atoms) {
      ++skipped_too_small;
      continue;
    }
    out.push_back(tpl);
  }

  return out;
}

std::set<size_t> collect_match_rows(const std::string &query_smarts,
                                    const std::vector<TargetMol> &targets,
                                    int num_workers = 6) {
  std::set<size_t> rows;
  if (num_workers < 1) {
    num_workers = 1;
  }

  const size_t n_targets = targets.size();
  if (n_targets == 0) {
    return rows;
  }

  size_t workers = static_cast<size_t>(num_workers);
  workers = std::min(workers, n_targets);

  std::vector<std::unique_ptr<RDKit::ROMol>> queries(workers);
  for (size_t i = 0; i < workers; ++i) {
    queries[i].reset(RDKit::SmartsToMol(query_smarts));
    if (!queries[i]) {
      throw std::runtime_error("Invalid SMARTS query: " + query_smarts);
    }
  }

  if (workers == 1) {
    for (const auto &t : targets) {
      if (t.mol && !RDKit::SubstructMatch(*t.mol, *queries[0]).empty()) {
        rows.insert(t.row_idx);
      }
    }
    return rows;
  }

  std::vector<std::set<size_t>> partial_rows(workers);
  std::vector<std::thread> threads;
  threads.reserve(workers);

  const size_t chunk = (n_targets + workers - 1) / workers;
  for (size_t w = 0; w < workers; ++w) {
    const size_t begin = w * chunk;
    const size_t end = std::min(begin + chunk, n_targets);
    if (begin >= end) {
      continue;
    }

    threads.emplace_back([&targets, &queries, &partial_rows, w, begin, end]() {
      auto &local_rows = partial_rows[w];
      const auto &query = queries[w];
      for (size_t i = begin; i < end; ++i) {
        const auto &t = targets[i];
        if (t.mol && !RDKit::SubstructMatch(*t.mol, *query).empty()) {
          local_rows.insert(t.row_idx);
        }
      }
    });
  }

  for (auto &th : threads) {
    th.join();
  }

  for (const auto &local : partial_rows) {
    rows.insert(local.begin(), local.end());
  }

  return rows;
}

int resolve_num_workers(int argc, char **argv) {
  constexpr int kDefaultWorkers = 6;
  if (argc > 4 && argv[4] && std::string(argv[4]).size() > 0) {
    try {
      int parsed = std::stoi(argv[4]);
      return parsed > 0 ? parsed : kDefaultWorkers;
    } catch (...) {
      return kDefaultWorkers;
    }
  }
  return kDefaultWorkers;
}

size_t resolve_target_limit(int argc, char **argv) {
  if (argc > 5 && argv[5] && std::string(argv[5]).size() > 0) {
    try {
      return static_cast<size_t>(std::stoull(argv[5]));
    } catch (...) {
      return 0;
    }
  }
  return 0;
}

bool resolve_random_selection(int argc, char **argv) {
  if (argc > 6 && argv[6] && std::string(argv[6]).size() > 0) {
    const std::string v = to_lower(trim(argv[6]));
    return v == "1" || v == "true" || v == "yes" || v == "y";
  }
  return false;
}

std::optional<unsigned int> resolve_random_seed(int argc, char **argv) {
  if (argc > 7 && argv[7] && std::string(argv[7]).size() > 0) {
    try {
      return static_cast<unsigned int>(std::stoul(argv[7]));
    } catch (...) {
      return std::nullopt;
    }
  }
  return std::nullopt;
}

std::set<size_t> set_intersection_of(const std::set<size_t> &a,
                                     const std::set<size_t> &b) {
  std::set<size_t> out;
  std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                        std::inserter(out, out.begin()));
  return out;
}

std::set<size_t> set_diff_of(const std::set<size_t> &a,
                             const std::set<size_t> &b) {
  std::set<size_t> out;
  std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
                      std::inserter(out, out.begin()));
  return out;
}

std::filesystem::path resolve_output_path(int argc, char **argv) {
  if (argc > 2 && argv[2] && std::string(argv[2]).size() > 0) {
    return std::filesystem::path(argv[2]);
  }

  const std::vector<std::filesystem::path> candidates = {
    std::filesystem::path("/home/bmahjour/test_rdpp/manuscript/results/smarts_analyzer_matchset_report.txt"),
      std::filesystem::path(
          "C:\\Projects\\atype\\results\\smarts_analyzer_matchset_report.txt"),
      std::filesystem::path("results\\smarts_analyzer_matchset_report.txt")};

  for (const auto &cand : candidates) {
    if (cand.is_absolute()) {
      const auto parent = cand.parent_path();
      if (parent.empty() || std::filesystem::exists(parent)) {
        return cand;
      }
    } else {
      return cand;
    }
  }

  return candidates.front();
}

double ratio(size_t num, size_t den) {
  return den == 0 ? 0.0 : static_cast<double>(num) / static_cast<double>(den);
}

double elapsed_seconds(
    const std::chrono::steady_clock::time_point &program_start,
    const std::chrono::steady_clock::time_point &t) {
  return std::chrono::duration<double>(t - program_start).count();
}

std::string csv_escape(const std::string &s) {
  bool needs_quotes = false;
  std::string out;
  out.reserve(s.size() + 8);
  for (char c : s) {
    if (c == '"') {
      out.push_back('"');
      out.push_back('"');
      needs_quotes = true;
    } else {
      if (c == ',' || c == '\n' || c == '\r') {
        needs_quotes = true;
      }
      out.push_back(c);
    }
  }
  if (needs_quotes) {
    return "\"" + out + "\"";
  }
  return out;
}

void write_set_details(
    std::ofstream &out, const std::string &label, const std::set<size_t> &rows,
    const std::unordered_map<size_t, std::string> &row_to_smiles) {
  out << "  " << label << " (" << rows.size() << ")" << std::endl;
  for (const auto row : rows) {
    const auto it = row_to_smiles.find(row);
    const std::string smiles =
        (it == row_to_smiles.end()) ? "<unknown>" : it->second;
    out << "    row=" << row << " smiles=" << smiles << std::endl;
  }
}

}  // namespace

int main(int argc, char **argv) {
  ZoneScopedN("smartsAnalyzerMatchsetSmoke::main");
  try {
    constexpr size_t kMinTemplateAtoms = 10;
    const auto program_start = std::chrono::steady_clock::now();
    const std::filesystem::path csv_path = resolve_dataset_path(argc, argv);
    // const std::filesystem::path template_path =
        // resolve_simple_templates_path(argc, argv);
    // const std::filesystem::path template_path =
        // resolve_templates_path(argc, argv);
    const std::filesystem::path template_path =
        resolve_large_templates_path(argc, argv);
    const std::filesystem::path report_path = resolve_output_path(argc, argv);
    const int num_workers = resolve_num_workers(argc, argv);
    const size_t target_limit = resolve_target_limit(argc, argv);
    const bool random_selection = resolve_random_selection(argc, argv);
    const std::optional<unsigned int> random_seed =
      resolve_random_seed(argc, argv);
    const std::filesystem::path csv_report_path =
        report_path.parent_path() /
        (report_path.stem().string() + "_stats.csv");
    const auto targets =
      load_targets_from_csv(csv_path, 5000, random_selection, 1);
    const auto raw_smarts_list = load_templates_from_csv(template_path, 100000);
    size_t skipped_too_small = 0;
    size_t skipped_invalid_templates = 0;
    const auto smarts_list = filter_templates_by_max_atoms(
      raw_smarts_list, kMinTemplateAtoms, skipped_too_small,
      skipped_invalid_templates);
    std::cout << "Dataset loaded: " << targets.size() << " molecules"
              << std::endl;
    std::cout << "Templates loaded: " << raw_smarts_list.size() << std::endl;
    std::cout << "Templates kept (>= " << kMinTemplateAtoms
          << " atoms): " << smarts_list.size() << std::endl;
    std::cout << "Templates skipped (< " << kMinTemplateAtoms
          << " atoms): " << skipped_too_small << std::endl;
    std::cout << "Templates skipped (invalid SMARTS): "
          << skipped_invalid_templates << std::endl;
    std::ofstream report(report_path);
    std::ofstream csv_report(csv_report_path);
    if (!report.good()) {
      throw std::runtime_error("Failed to open report file for writing: " +
                               report_path.string());
    }
    if (!csv_report.good()) {
      throw std::runtime_error("Failed to open CSV report file for writing: " +
                               csv_report_path.string());
    }
    if (smarts_list.empty()) {
      throw std::runtime_error(
          "No SMARTS templates remaining after filtering in: " +
          template_path.string());
    }

    std::unordered_map<size_t, std::string> row_to_smiles;
    row_to_smiles.reserve(targets.size());
    for (const auto &t : targets) {
      row_to_smiles.emplace(t.row_idx, t.smiles);
    }

    std::cout << "Loaded target molecules: " << targets.size() << std::endl;
    std::cout << "CSV path: " << csv_path.string() << std::endl << std::endl;
    std::cout << "Templates loaded: " << raw_smarts_list.size() << std::endl;
    std::cout << "Templates kept (>= " << kMinTemplateAtoms
          << " atoms): " << smarts_list.size() << std::endl;
    std::cout << "Templates skipped (< " << kMinTemplateAtoms
          << " atoms): " << skipped_too_small << std::endl;
    std::cout << "Templates skipped (invalid SMARTS): "
          << skipped_invalid_templates << std::endl;
    std::cout << "Template CSV path: " << template_path.string() << std::endl
              << std::endl;
    std::cout << "Worker threads: " << num_workers << std::endl;
    std::cout << "Target limit: " << target_limit << std::endl;
    std::cout << "Random selection: " << (random_selection ? "true" : "false")
              << std::endl;
    if (random_seed.has_value()) {
      std::cout << "Random seed: " << *random_seed << std::endl;
    }
    std::cout << "Writing report to: " << report_path.string() << std::endl;

    report << "SMARTS Matchset Comparison Report" << std::endl;
    report << "CSV path: " << csv_path.string() << std::endl;
    report << "Template CSV path: " << template_path.string() << std::endl;
    report << "Target molecules loaded: " << targets.size() << std::endl;
        report << "Templates loaded: " << raw_smarts_list.size() << std::endl;
        report << "Template min atoms filter: " << kMinTemplateAtoms
          << std::endl;
        report << "Templates kept (>= min atoms): " << smarts_list.size()
          << std::endl;
        report << "Templates skipped (< min atoms): " << skipped_too_small
          << std::endl;
        report << "Templates skipped (invalid SMARTS): "
          << skipped_invalid_templates << std::endl;
    report << "Worker threads: " << num_workers << std::endl;
    report << "Target limit: " << target_limit << std::endl;
    report << "Random selection: " << (random_selection ? "true" : "false")
           << std::endl;
    if (random_seed.has_value()) {
      report << "Random seed: " << *random_seed << std::endl;
    }
    report << std::endl;

    csv_report << "query_idx,input_smarts,standardized_smarts,before_matches,"
                  "after_matches,intersection_size,union_size,before_only_size,"
                  "after_only_size,jaccard,precision,recall,f1,exact_match,"
                  "time_start,time_end"
               << std::endl;

    atom_typer::SmartsAnalyzer::StandardSmartsWorkflowOptions workflow_options;
    workflow_options.include_x_in_reserialization = false;
    workflow_options.enumerate_bond_order = false;
    

    atom_typer::SmartsAnalyzer sa;

    size_t total_before = 0;
    size_t total_after = 0;
    size_t total_intersection = 0;
    size_t total_union = 0;
    size_t total_before_only = 0;
    size_t total_after_only = 0;
    size_t exact_match_queries = 0;

    size_t query_idx = 0;
    size_t input_failed = 0;
    for (const auto &input_smarts : smarts_list) {
      ZoneScopedN("matchset_query");
      const auto entry_start = std::chrono::steady_clock::now();
      ++query_idx;
      std::vector<std::string> standardized_vec;
      try {
        standardized_vec = sa.standard_smarts({input_smarts}, false, false, false, workflow_options);
      } catch (const std::exception &e) {
        std::cerr << "Failed while standardizing input SMARTS '" << input_smarts
                  << "': " << e.what() << std::endl;
        continue;
      }
      if (standardized_vec.empty()) {
        std::cerr << "Failed to standardize SMARTS: " << input_smarts
                  << std::endl;
        continue;
      }

      const std::string standardized = standardized_vec.front();

      const auto before_set =
          collect_match_rows(input_smarts, targets, num_workers);
      const auto after_set =
          collect_match_rows(standardized, targets, num_workers);

      const auto intersection_set = set_intersection_of(before_set, after_set);
      const std::set<size_t> union_set = [&]() {
        std::set<size_t> out;
        std::set_union(before_set.begin(), before_set.end(), after_set.begin(),
                       after_set.end(), std::inserter(out, out.begin()));
        return out;
      }();
      const auto before_only = set_diff_of(before_set, after_set);
      const auto after_only = set_diff_of(after_set, before_set);

      const double jaccard = ratio(intersection_set.size(), union_set.size());
      const double precision = ratio(intersection_set.size(), after_set.size());
      const double recall = ratio(intersection_set.size(), before_set.size());
      const double f1 = (precision + recall == 0.0)
                            ? 0.0
                            : (2.0 * precision * recall) / (precision + recall);
      const auto entry_end = std::chrono::steady_clock::now();
      const double time_start = elapsed_seconds(program_start, entry_start);
      const double time_end = elapsed_seconds(program_start, entry_end);

      total_before += before_set.size();
      total_after += after_set.size();
      total_intersection += intersection_set.size();
      total_union += union_set.size();
      total_before_only += before_only.size();
      total_after_only += after_only.size();
      if (before_set == after_set) {
        ++exact_match_queries;
      }

      // std::cout << "input:      " << input_smarts << std::endl;
      // std::cout << "standard:   " << standardized << std::endl;
      // std::cout << "  before matches: " << before_set.size() << std::endl;
      // std::cout << "  after matches:  " << after_set.size() << std::endl;
      // std::cout << "  union matches:  " << union_set.size() << std::endl;
      // std::cout << "  jaccard: " << jaccard << std::endl;
      // std::cout << std::endl;

      csv_report << query_idx << "," << csv_escape(input_smarts) << ","
                 << csv_escape(standardized) << "," << before_set.size() << ","
                 << after_set.size() << "," << intersection_set.size() << ","
                 << union_set.size() << "," << before_only.size() << ","
                 << after_only.size() << "," << jaccard << "," << precision
                 << "," << recall << "," << f1 << ","
                 << ((before_set == after_set) ? 1 : 0) << "," << time_start
                 << "," << time_end << std::endl;

      if (union_set.size() == 0) {
        continue;
      }
      report << "============================================================"
             << std::endl;
      report << "input SMARTS:        " << input_smarts << std::endl;
      report << "standardized SMARTS: " << standardized << std::endl;
      report << "before matches:      " << before_set.size() << std::endl;
      report << "after matches:       " << after_set.size() << std::endl;
      report << "intersection size:   " << intersection_set.size() << std::endl;
      report << "union size:          " << union_set.size() << std::endl;
      report << "before-only size:    " << before_only.size() << std::endl;
      report << "after-only size:     " << after_only.size() << std::endl;
      report << "jaccard:             " << jaccard << std::endl;
      report << "precision(after|before): " << precision << std::endl;
      report << "recall(before|after):    " << recall << std::endl;
      report << "f1:                  " << f1 << std::endl;

      write_set_details(report, "before-only", before_only, row_to_smiles);
      write_set_details(report, "after-only", after_only, row_to_smiles);
      report << std::endl;
    }

    const auto program_end = std::chrono::steady_clock::now();
    const double total_elapsed = elapsed_seconds(program_start, program_end);

    const size_t query_count = smarts_list.size();
    report << "############################################################"
           << std::endl;
    report << "GLOBAL SUMMARY" << std::endl;
    report << "queries:                 " << query_count << std::endl;
    report << "exact-match queries:     " << exact_match_queries << std::endl;
    report << "sum(before):             " << total_before << std::endl;
    report << "sum(after):              " << total_after << std::endl;
    report << "sum(intersection):       " << total_intersection << std::endl;
    report << "sum(union):              " << total_union << std::endl;
    report << "sum(before-only):        " << total_before_only << std::endl;
    report << "sum(after-only):         " << total_after_only << std::endl;
    report << "global jaccard (sum):    "
           << ratio(total_intersection, total_union) << std::endl;
    report << "global precision (sum):  "
           << ratio(total_intersection, total_after) << std::endl;
    report << "global recall (sum):     "
           << ratio(total_intersection, total_before) << std::endl;

    csv_report << "GLOBAL,,," << total_before << "," << total_after << ","
               << total_intersection << "," << total_union << ","
               << total_before_only << "," << total_after_only << ","
               << ratio(total_intersection, total_union) << ","
               << ratio(total_intersection, total_after) << ","
               << ratio(total_intersection, total_before) << ",," << 0.0 << ","
               << total_elapsed << std::endl;

    std::cout << "Report complete: " << report_path.string() << std::endl;
    std::cout << "CSV stats complete: " << csv_report_path.string()
              << std::endl;
  } catch (const std::exception &e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
