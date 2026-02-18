#include "atom_typer_query_reorder.hpp"

#include <GraphMol/GraphMol.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <Query/EqualityQuery.h>

#include <algorithm>
#include <cctype>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>

namespace atom_typer::query_reorder {
namespace {

bool g_comparison_trace_enabled = false;

bool is_atom_and_query(const std::string &desc) {
  return desc == "AtomAnd" || desc == "And";
}

bool is_atom_or_query(const std::string &desc) {
  return desc == "AtomOr" || desc == "Or";
}

bool is_recursive_query(const std::string &desc) {
  return desc.find("Recursive") != std::string::npos ||
         desc.find("recursive") != std::string::npos;
}

bool is_atom_identity_primitive(const RDKit::Atom::QUERYATOM_QUERY *query) {
  if (!query || query->getNegation()) {
    return false;
  }
  const std::string desc = query->getDescription();
  if (desc == "AtomAtomicNum" || desc == "AtomNum" || desc == "AtomType") {
    return true;
  }
  return false;
}

bool is_atom_identity_token(const std::string &token) {
  if (token.empty()) {
    return false;
  }
  if (token[0] == '#') {
    if (token.size() < 2) {
      return false;
    }
    return std::all_of(token.begin() + 1, token.end(), [](char c) {
      return std::isdigit(static_cast<unsigned char>(c)) != 0;
    });
  }

  // Element/aromatic symbols and wildcard primitive.
  if (token == "*" || token == "a" || token == "A") {
    return true;
  }
  if (token.size() == 1 &&
      std::isalpha(static_cast<unsigned char>(token[0])) != 0) {
    return true;
  }
  if (token.size() == 2 &&
      std::isalpha(static_cast<unsigned char>(token[0])) != 0 &&
      std::isalpha(static_cast<unsigned char>(token[1])) != 0) {
    return true;
  }
  return false;
}

std::string move_identity_term_first_in_atom_inner(const std::string &inner) {
  if (inner.empty()) {
    return inner;
  }

  int paren_depth = 0;
  int bracket_depth = 0;
  std::vector<size_t> top_level_ampersands;
  size_t top_level_map_colon = std::string::npos;

  for (size_t i = 0; i < inner.size(); ++i) {
    const char ch = inner[i];
    if (ch == '[') {
      ++bracket_depth;
      continue;
    }
    if (ch == ']') {
      --bracket_depth;
      continue;
    }
    if (ch == '(') {
      ++paren_depth;
      continue;
    }
    if (ch == ')') {
      --paren_depth;
      continue;
    }
    if (paren_depth == 0 && bracket_depth == 0) {
      if (ch == '&') {
        top_level_ampersands.push_back(i);
      } else if (ch == ':') {
        top_level_map_colon = i;
      }
    }
  }

  if (top_level_ampersands.empty()) {
    return inner;
  }

  std::string body = inner;
  std::string map_suffix;
  if (top_level_map_colon != std::string::npos &&
      top_level_map_colon + 1 < inner.size()) {
    bool digits_only = true;
    for (size_t i = top_level_map_colon + 1; i < inner.size(); ++i) {
      if (!std::isdigit(static_cast<unsigned char>(inner[i]))) {
        digits_only = false;
        break;
      }
    }
    if (digits_only) {
      body = inner.substr(0, top_level_map_colon);
      map_suffix = inner.substr(top_level_map_colon);
    }
  }

  std::vector<std::string> terms;
  std::string cur;
  paren_depth = 0;
  bracket_depth = 0;
  for (size_t i = 0; i < body.size(); ++i) {
    const char ch = body[i];
    if (ch == '[') {
      ++bracket_depth;
      cur.push_back(ch);
      continue;
    }
    if (ch == ']') {
      --bracket_depth;
      cur.push_back(ch);
      continue;
    }
    if (ch == '(') {
      ++paren_depth;
      cur.push_back(ch);
      continue;
    }
    if (ch == ')') {
      --paren_depth;
      cur.push_back(ch);
      continue;
    }

    if (ch == '&' && paren_depth == 0 && bracket_depth == 0) {
      terms.push_back(cur);
      cur.clear();
    } else {
      cur.push_back(ch);
    }
  }
  terms.push_back(cur);

  size_t id_idx = terms.size();
  for (size_t i = 0; i < terms.size(); ++i) {
    if (is_atom_identity_token(terms[i])) {
      id_idx = i;
      break;
    }
  }
  if (id_idx == terms.size() || id_idx == 0) {
    return inner;
  }

  std::string rebuilt = terms[id_idx];
  for (size_t i = 0; i < terms.size(); ++i) {
    if (i == id_idx) {
      continue;
    }
    rebuilt += "&" + terms[i];
  }
  rebuilt += map_suffix;
  return rebuilt;
}

std::string force_atom_identity_first_in_smarts(const std::string &smarts) {
  if (smarts.empty()) {
    return smarts;
  }

  std::string out;
  out.reserve(smarts.size());

  size_t cursor = 0;
  size_t token_start = std::string::npos;
  int bracket_depth = 0;
  for (size_t i = 0; i < smarts.size(); ++i) {
    const char ch = smarts[i];
    if (ch == '[') {
      if (bracket_depth == 0) {
        token_start = i;
      }
      ++bracket_depth;
      continue;
    }
    if (ch == ']') {
      --bracket_depth;
      if (bracket_depth == 0 && token_start != std::string::npos) {
        const size_t token_end = i;
        out += smarts.substr(cursor, token_start - cursor);

        const std::string token = smarts.substr(token_start, token_end - token_start + 1);
        if (token.size() >= 2) {
          const std::string inner = token.substr(1, token.size() - 2);
          out += "[" + move_identity_term_first_in_atom_inner(inner) + "]";
        } else {
          out += token;
        }
        cursor = token_end + 1;
        token_start = std::string::npos;
      }
    }
  }

  out += smarts.substr(cursor);
  return out;
}

size_t matching_paren(const std::string &txt, size_t open_pos) {
  int depth = 0;
  for (size_t i = open_pos; i < txt.size(); ++i) {
    if (txt[i] == '(') {
      ++depth;
    } else if (txt[i] == ')') {
      --depth;
      if (depth == 0) {
        return i;
      }
    }
  }
  return std::string::npos;
}

std::string normalize_embedding_token(std::string token) {
  if (token.size() >= 2 && token.front() == '[' && token.back() == ']') {
    token = token.substr(1, token.size() - 2);
  }
  return token;
}

std::string primitive_token_from_query(
    const RDKit::Atom::QUERYATOM_QUERY *query) {
  if (!query) {
    return "*";
  }

  const auto *eq = dynamic_cast<const RDKit::ATOM_EQUALS_QUERY *>(query);
  const std::string desc = query->getDescription();
  std::string token;

  if (eq) {
    const int val = eq->getVal();
    if (desc == "AtomFormalCharge") {
      token = (val >= 0) ? ("+" + std::to_string(val)) : std::to_string(val);
    } else if (desc == "AtomHCount") {
      token = "H" + std::to_string(val);
    } else if (desc == "AtomImplicitHCount") {
      token = "h" + std::to_string(val);
    } else if (desc == "AtomExplicitDegree") {
      token = "D" + std::to_string(val);
    } else if (desc == "AtomTotalDegree") {
      token = "X" + std::to_string(val);
    } else if (desc == "AtomAtomicNum" || desc == "AtomNum") {
      token = "#" + std::to_string(val);
    }
  } else if (desc == "AtomAromatic") {
    token = "a";
  } else if (desc == "AtomAliphatic") {
    token = "A";
  }

  if (!token.empty()) {
    if (query->getNegation()) {
      return "!" + token;
    }
    return token;
  }

  RDKit::QueryAtom tmp(0);
  tmp.setQuery(query->copy());
  return normalize_embedding_token(RDKit::SmartsWrite::GetAtomSmarts(&tmp));
}

double default_embedding_score(const std::map<std::string, double> &embedding) {
  if (embedding.empty()) {
    return 1.0;
  }
  double max_val = std::numeric_limits<double>::lowest();
  for (const auto &kv : embedding) {
    max_val = std::max(max_val, kv.second);
  }
  return max_val + 1.0;
}

double score_query_node(const RDKit::Atom::QUERYATOM_QUERY *query,
                        const std::map<std::string, double> &embedding,
                        double default_score) {
  if (!query) {
    return default_score;
  }

  const std::string desc = query->getDescription();
  if (is_atom_and_query(desc) || is_atom_or_query(desc)) {
    double score = default_score;
    bool has_child = false;
    for (auto it = query->beginChildren(); it != query->endChildren(); ++it) {
      score = std::min(score,
                       score_query_node(it->get(), embedding, default_score));
      has_child = true;
    }
    return has_child ? score : default_score;
  }

  const std::string token = primitive_token_from_query(query);
  const auto found = embedding.find(token);
  if (found != embedding.end()) {
    return found->second;
  }
  return default_score;
}

std::unique_ptr<RDKit::Atom::QUERYATOM_QUERY> query_from_atom_token(
    const std::string &token) {
  const std::string atom_smarts = "[" + token + "]";
  std::unique_ptr<RDKit::ROMol> one_atom_mol(RDKit::SmartsToMol(atom_smarts));
  if (!one_atom_mol || one_atom_mol->getNumAtoms() != 1) {
    return nullptr;
  }

  auto *src_atom = one_atom_mol->getAtomWithIdx(0);
  auto *src_qatom = dynamic_cast<RDKit::QueryAtom *>(src_atom);
  if (!src_qatom || !src_qatom->getQuery()) {
    return nullptr;
  }

  return std::unique_ptr<RDKit::Atom::QUERYATOM_QUERY>(
      src_qatom->getQuery()->copy());
}

std::string reorder_recursive_smarts_keep_anchor(
    const std::string &recursive_smarts,
    const std::map<std::string, double> &embedding, double default_score);

std::unique_ptr<RDKit::Atom::QUERYATOM_QUERY> reorder_atom_query_tree(
    const RDKit::Atom::QUERYATOM_QUERY *query,
    const std::map<std::string, double> &embedding, double default_score,
    bool lock_first_child = false) {
  if (!query) {
    return nullptr;
  }

  const std::string desc = query->getDescription();
  const bool is_and = is_atom_and_query(desc);
  const bool is_or = is_atom_or_query(desc);

  if (!is_and && !is_or && is_recursive_query(desc)) {
    const std::string token = primitive_token_from_query(query);
    const size_t start = token.find("$(");
    if (start != std::string::npos) {
      const size_t open_pos = start + 1;
      const size_t close_pos = matching_paren(token, open_pos);
      if (close_pos != std::string::npos && close_pos > open_pos + 1) {
        const std::string prefix = token.substr(0, start + 2);
        const std::string inner =
            token.substr(open_pos + 1, close_pos - open_pos - 1);
        const std::string suffix = token.substr(close_pos);

        const std::string rewritten_inner =
            reorder_recursive_smarts_keep_anchor(inner, embedding,
                                                 default_score);

        const std::string rewritten_token =
            prefix + rewritten_inner + suffix;
        auto rebuilt = query_from_atom_token(rewritten_token);
        if (rebuilt) {
          return rebuilt;
        }
      }
    }
  }

  if (!is_and && !is_or) {
    return std::unique_ptr<RDKit::Atom::QUERYATOM_QUERY>(query->copy());
  }

  struct RankedChild {
    std::unique_ptr<RDKit::Atom::QUERYATOM_QUERY> query;
    double score = 0.0;
    size_t original_idx = 0;
    std::string token;
  };

  std::vector<RankedChild> children;
  size_t idx = 0;
  for (auto it = query->beginChildren(); it != query->endChildren();
       ++it, ++idx) {
    auto child =
        reorder_atom_query_tree(it->get(), embedding, default_score, false);
    if (!child) {
      continue;
    }
    const double score = score_query_node(child.get(), embedding, default_score);
    children.push_back(
        {std::move(child), score, idx, primitive_token_from_query(child.get())});
  }

  if (g_comparison_trace_enabled && is_and) {
    std::cout << "[query_reorder] AND node before sort" << std::endl;
    for (const auto &c : children) {
      std::cout << "  idx=" << c.original_idx << " token='" << c.token
                << "' score=" << c.score << std::endl;
    }
  }

  RankedChild locked_child;
  bool has_locked_child = false;

  // Keep a stable atom identity primitive first in AtomAnd nodes.
  if (is_and) {
    auto lock_it = std::find_if(children.begin(), children.end(),
                                [](const RankedChild &c) {
                                  return is_atom_identity_primitive(c.query.get());
                                });
    if (lock_it != children.end()) {
      locked_child = std::move(*lock_it);
      children.erase(lock_it);
      has_locked_child = true;
      if (g_comparison_trace_enabled) {
        std::cout << "[query_reorder] AND lock atom-identity token='"
                  << locked_child.token << "'" << std::endl;
      }
    }
  }

  // Recursive-anchor lock has higher precedence when explicitly requested.
  if (lock_first_child && !children.empty()) {
    locked_child = std::move(children.front());
    children.erase(children.begin());
    has_locked_child = true;
    if (g_comparison_trace_enabled) {
      std::cout << "[query_reorder] lock_first_child override token='"
                << locked_child.token << "'" << std::endl;
    }
  }

  std::stable_sort(children.begin(), children.end(), [is_and](
                                                      const RankedChild &a,
                                                      const RankedChild &b) {
    if (g_comparison_trace_enabled && is_and) {
      std::cout << "[query_reorder] compare AND lhs=('" << a.token
                << "', score=" << a.score << ") rhs=('" << b.token
                << "', score=" << b.score << ")" << std::endl;
    }
                     if (a.score == b.score) {
                       return a.original_idx < b.original_idx;
                     }
                     return is_and ? (a.score < b.score) : (a.score > b.score);
                   });

  if (g_comparison_trace_enabled && is_and) {
    std::cout << "[query_reorder] AND node after sort" << std::endl;
    for (const auto &c : children) {
      std::cout << "  idx=" << c.original_idx << " token='" << c.token
                << "' score=" << c.score << std::endl;
    }
  }

  std::unique_ptr<RDKit::Atom::QUERYATOM_QUERY> out(
      is_and ? static_cast<RDKit::Atom::QUERYATOM_QUERY *>(
                   new RDKit::ATOM_AND_QUERY)
             : static_cast<RDKit::Atom::QUERYATOM_QUERY *>(
                   new RDKit::ATOM_OR_QUERY));
  out->setDescription(is_and ? "AtomAnd" : "AtomOr");
  out->setNegation(query->getNegation());
  out->setTypeLabel(query->getTypeLabel());

    if (has_locked_child) {
    out->addChild(
      RDKit::Atom::QUERYATOM_QUERY::CHILD_TYPE(locked_child.query.release()));
  }

  for (auto &child : children) {
    out->addChild(RDKit::Atom::QUERYATOM_QUERY::CHILD_TYPE(child.query.release()));
  }

  return out;
}

std::string reorder_recursive_smarts_keep_anchor(
    const std::string &recursive_smarts,
    const std::map<std::string, double> &embedding, double default_score) {
  std::unique_ptr<RDKit::ROMol> recursive_mol(RDKit::SmartsToMol(recursive_smarts));
  if (!recursive_mol) {
    return recursive_smarts;
  }

  RDKit::RWMol rewritten_recursive(*recursive_mol);
  for (auto *atom : rewritten_recursive.atoms()) {
    auto *qatom = dynamic_cast<RDKit::QueryAtom *>(atom);
    if (!qatom || !qatom->hasQuery() || !qatom->getQuery()) {
      continue;
    }

    const bool lock_first = (atom->getIdx() == 0u);
    auto reordered = reorder_atom_query_tree(qatom->getQuery(), embedding,
                                             default_score, lock_first);
    if (!reordered) {
      continue;
    }
    qatom->setQuery(reordered.release());
  }

  return RDKit::MolToSmarts(rewritten_recursive);
}

}  // namespace

void set_comparison_trace_enabled(bool enabled) {
  g_comparison_trace_enabled = enabled;
}

bool comparison_trace_enabled() {
  return g_comparison_trace_enabled;
}

std::string reorder_query_tree_by_embedding_smarts(
    const std::string &smarts, const std::map<std::string, double> &embedding) {
  if (smarts.empty()) {
    throw std::runtime_error("Empty SMARTS string");
  }

  std::unique_ptr<RDKit::ROMol> mol(RDKit::SmartsToMol(smarts));
  if (!mol) {
    throw std::runtime_error("Failed to parse SMARTS: " + smarts);
  }

  RDKit::RWMol rewritten(*mol);
  const double default_score = default_embedding_score(embedding);

  for (auto *atom : rewritten.atoms()) {
    auto *qatom = dynamic_cast<RDKit::QueryAtom *>(atom);
    if (!qatom || !qatom->hasQuery() || !qatom->getQuery()) {
      continue;
    }

    auto reordered =
      reorder_atom_query_tree(qatom->getQuery(), embedding, default_score,
                  false);
    if (!reordered) {
      continue;
    }
    qatom->setQuery(reordered.release());
  }

  return force_atom_identity_first_in_smarts(RDKit::MolToSmarts(rewritten));
}

std::vector<std::string> reorder_patterns_by_embedding(
    const std::vector<std::string> &patterns,
    const std::map<std::string, double> &embedding) {
  std::vector<std::string> out;
  out.reserve(patterns.size());
  for (const auto &p : patterns) {
    if (p.empty()) {
      continue;
    }
    out.push_back(reorder_query_tree_by_embedding_smarts(p, embedding));
  }
  return out;
}

}  // namespace atom_typer::query_reorder
