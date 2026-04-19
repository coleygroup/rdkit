#include "expression_builder.hpp"
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include <GraphMol/QueryAtom.h>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>
#include <tuple>
#include <cctype>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <map>

namespace atom_typer {

namespace {

std::string trim(const std::string &input) {
    size_t first = 0;
    while (first < input.size() &&
           std::isspace(static_cast<unsigned char>(input[first]))) {
        ++first;
    }
    size_t last = input.size();
    while (last > first &&
           std::isspace(static_cast<unsigned char>(input[last - 1]))) {
        --last;
    }
    return input.substr(first, last - first);
}

std::string toLowerCopy(std::string value) {
    for (auto &ch : value) {
        ch = static_cast<char>(
            std::tolower(static_cast<unsigned char>(ch)));
    }
    return value;
}

std::vector<std::string> normalizePrimitiveList(
    const std::vector<std::string> &primitiveList) {
    std::vector<std::string> tokens;
    for (const auto &entry : primitiveList) {
        std::string work = trim(entry);
        if (work.empty()) {
            continue;
        }
        if (!work.empty() && work.front() == '[') {
            work.erase(work.begin());
        }
        if (!work.empty() && work.back() == ']') {
            work.pop_back();
        }

        size_t start = 0;
        while (start < work.size()) {
            const size_t comma = work.find(',', start);
            const size_t end = (comma == std::string::npos) ? work.size() : comma;
            std::string token = trim(work.substr(start, end - start));
            if (!token.empty()) {
                tokens.push_back(token);
            }
            if (comma == std::string::npos) {
                break;
            }
            start = comma + 1;
        }
    }
    auto primitiveOrderPriority = [](const std::string &token) {
        const std::string lower = toLowerCopy(token);
        if (lower == "atomic_number" || lower == "atomicnum" ||
            lower == "atomatomicnum") {
            return 0;
        }
        if (lower == "element" || lower == "atom" || lower == "symbol") {
            return 1;
        }
        return 2;
    };

    // Keep user order for most primitives, but force atomic identity keys first.
    std::stable_sort(tokens.begin(), tokens.end(),
                     [&](const std::string &lhs, const std::string &rhs) {
                         return primitiveOrderPriority(lhs) <
                                primitiveOrderPriority(rhs);
                     });

    return tokens;
}

int bondSerializationPriority(const RDKit::Bond *bond) {
    if (!bond) {
        return 5;
    }
    if (bond->getIsAromatic()) {
        return 2;
    }
    switch (bond->getBondType()) {
        case RDKit::Bond::BondType::TRIPLE:
            return 0;
        case RDKit::Bond::BondType::DOUBLE:
            return 1;
        case RDKit::Bond::BondType::SINGLE:
            return 3;
        default:
            return 4;
    }
}

std::string bondToken(const RDKit::Bond *bond) {
    if (!bond) {
        return "~";
    }
    if (bond->getIsAromatic()) {
        return ":";
    }
    switch (bond->getBondType()) {
        case RDKit::Bond::BondType::TRIPLE:
            return "#";
        case RDKit::Bond::BondType::DOUBLE:
            return "=";
        case RDKit::Bond::BondType::SINGLE:
            return "-";
        default:
            return "~";
    }
}

std::string bracketUnbracketedWildcardAtoms(const std::string &smarts) {
    std::string out;
    out.reserve(smarts.size() + 8);

    int bracketDepth = 0;
    for (char ch : smarts) {
        if (ch == '[') {
            ++bracketDepth;
            out.push_back(ch);
            continue;
        }
        if (ch == ']') {
            if (bracketDepth > 0) {
                --bracketDepth;
            }
            out.push_back(ch);
            continue;
        }
        if (ch == '*' && bracketDepth == 0) {
            out += "[*]";
            continue;
        }
        out.push_back(ch);
    }

    return out;
}

std::string buildWildcardNeighborSmarts(
    const RDKit::ROMol &mol, unsigned int centerIdx,
    const RDKit::PATH_TYPE &environmentBonds,
    const std::string &rootAtomToken) {
    std::vector<const RDKit::Bond *> centerBonds;
    centerBonds.reserve(environmentBonds.size());

    for (const auto bondIdx : environmentBonds) {
        const auto *bond = mol.getBondWithIdx(static_cast<unsigned int>(bondIdx));
        if (!bond) {
            continue;
        }
        if (bond->getBeginAtomIdx() == centerIdx ||
            bond->getEndAtomIdx() == centerIdx) {
            centerBonds.push_back(bond);
        }
    }

    if (centerBonds.empty()) {
        return rootAtomToken;
    }

    std::stable_sort(centerBonds.begin(), centerBonds.end(),
                     [&](const RDKit::Bond *lhs, const RDKit::Bond *rhs) {
                         const int lhsPriority = bondSerializationPriority(lhs);
                         const int rhsPriority = bondSerializationPriority(rhs);
                         if (lhsPriority != rhsPriority) {
                             return lhsPriority < rhsPriority;
                         }

                         const unsigned int lhsOther =
                             lhs->getBeginAtomIdx() == centerIdx
                                 ? lhs->getEndAtomIdx()
                                 : lhs->getBeginAtomIdx();
                         const unsigned int rhsOther =
                             rhs->getBeginAtomIdx() == centerIdx
                                 ? rhs->getEndAtomIdx()
                                 : rhs->getBeginAtomIdx();
                         if (lhsOther != rhsOther) {
                             return lhsOther < rhsOther;
                         }

                         return lhs->getIdx() < rhs->getIdx();
                     });

    std::string out = rootAtomToken;
    if (centerBonds.size() == 1) {
        out += bondToken(centerBonds.front()) + "[*]";
        return out;
    }

    for (size_t i = 0; i < centerBonds.size(); ++i) {
        const std::string segment = bondToken(centerBonds[i]) + "[*]";
        if (i + 1 == centerBonds.size()) {
            out += segment;
        } else {
            out += "(" + segment + ")";
        }
    }
    return out;
}

size_t rootedAtomTokenEnd(const std::string &smarts) {
    if (smarts.empty() || smarts.front() != '[') {
        return std::string::npos;
    }
    int bracketDepth = 0;
    for (size_t i = 0; i < smarts.size(); ++i) {
        const char ch = smarts[i];
        if (ch == '[') {
            ++bracketDepth;
        } else if (ch == ']') {
            --bracketDepth;
            if (bracketDepth == 0) {
                return i;
            }
        }
    }
    return std::string::npos;
}

std::string buildAtomTokenFromPrimitiveValues(
    const std::vector<std::string> &primitiveValues,
    std::uint64_t subsetMask) {
    std::stringstream ss;
    ss << "[";
    bool first = true;
    for (size_t i = 0; i < primitiveValues.size(); ++i) {
        if ((subsetMask & (std::uint64_t{1} << i)) == 0) {
            continue;
        }
        if (!first) {
            ss << "&";
        }
        ss << primitiveValues[i];
        first = false;
    }
    ss << "]";
    return ss.str();
}

std::unique_ptr<RDKit::ROMol> reorderSubmolBondsDeterministically(
    const RDKit::ROMol &envSubmol, const std::map<int, int> &atomIdxMap,
    unsigned int parentCenterIdx) {
    std::vector<int> subToParent(envSubmol.getNumAtoms(), -1);
    for (const auto &entry : atomIdxMap) {
        const int parentIdx = entry.first;
        const int subIdx = entry.second;
        if (subIdx >= 0 && static_cast<size_t>(subIdx) < subToParent.size()) {
            subToParent[static_cast<size_t>(subIdx)] = parentIdx;
        }
    }

    int subCenterIdx = -1;
    auto centerIt = atomIdxMap.find(static_cast<int>(parentCenterIdx));
    if (centerIt != atomIdxMap.end()) {
        subCenterIdx = centerIt->second;
    }

    std::vector<const RDKit::Bond *> bonds;
    bonds.reserve(envSubmol.getNumBonds());
    for (const auto bond : envSubmol.bonds()) {
        bonds.push_back(bond);
    }

    std::stable_sort(bonds.begin(), bonds.end(),
                     [&](const RDKit::Bond *lhs, const RDKit::Bond *rhs) {
                         const unsigned int lhsBegin = lhs->getBeginAtomIdx();
                         const unsigned int lhsEnd = lhs->getEndAtomIdx();
                         const unsigned int rhsBegin = rhs->getBeginAtomIdx();
                         const unsigned int rhsEnd = rhs->getEndAtomIdx();

                         const bool lhsTouchesCenter =
                             (static_cast<int>(lhsBegin) == subCenterIdx ||
                              static_cast<int>(lhsEnd) == subCenterIdx);
                         const bool rhsTouchesCenter =
                             (static_cast<int>(rhsBegin) == subCenterIdx ||
                              static_cast<int>(rhsEnd) == subCenterIdx);

                         if (lhsTouchesCenter != rhsTouchesCenter) {
                             return lhsTouchesCenter > rhsTouchesCenter;
                         }

                         const int lhsPriority = bondSerializationPriority(lhs);
                         const int rhsPriority = bondSerializationPriority(rhs);
                         if (lhsPriority != rhsPriority) {
                             return lhsPriority < rhsPriority;
                         }

                         const int lhsP0 = subToParent[lhsBegin];
                         const int lhsP1 = subToParent[lhsEnd];
                         const int rhsP0 = subToParent[rhsBegin];
                         const int rhsP1 = subToParent[rhsEnd];

                         const auto lhsPair = std::minmax(lhsP0, lhsP1);
                         const auto rhsPair = std::minmax(rhsP0, rhsP1);
                         if (lhsPair != rhsPair) {
                             return lhsPair < rhsPair;
                         }

                         return lhs->getIdx() < rhs->getIdx();
                     });

    auto *ordered = new RDKit::RWMol();
    for (size_t i = 0; i < envSubmol.getNumAtoms(); ++i) {
        const auto *atom = envSubmol.getAtomWithIdx(i);
        auto *queryAtom = new RDKit::QueryAtom(atom->getAtomicNum());
        if (atom->hasQuery()) {
            queryAtom->setQuery(atom->getQuery()->copy());
        }
        queryAtom->setFormalCharge(atom->getFormalCharge());
        queryAtom->setIsAromatic(atom->getIsAromatic());
        ordered->addAtom(queryAtom, false, true);
    }

    for (const auto *bond : bonds) {
        ordered->addBond(bond->getBeginAtomIdx(), bond->getEndAtomIdx(),
                         bond->getBondType());
        auto *newBond =
            ordered->getBondWithIdx(ordered->getNumBonds() - 1);
        newBond->setIsAromatic(bond->getIsAromatic());
    }

    return std::unique_ptr<RDKit::ROMol>(ordered);
}

std::string primitiveValue(const RDKit::Atom *atom,
                           const RDKit::ROMol *mol,
                           const std::string &primitive) {
    const std::string primitiveLower = toLowerCopy(primitive);
    if (primitiveLower == "element" || primitiveLower == "atom" ||
        primitiveLower == "symbol") {
        return "#" + std::to_string(atom->getAtomicNum());
    }
    if (primitive == "D" || primitiveLower == "degree" ||
        primitiveLower == "atomexplicitdegree") {
        return "D" + std::to_string(atom->getDegree());
    }
    if (primitive == "H" || primitiveLower == "hcount" ||
        primitiveLower == "atomhcount") {
        return "H" + std::to_string(atom->getTotalNumHs(true));
    }
    if (primitive == "h" || primitiveLower == "atomimplicithcount") {
        return "h" + std::to_string(atom->getNumImplicitHs());
    }
    if (primitiveLower == "charge" || primitiveLower == "formalcharge" ||
        primitiveLower == "atomformalcharge") {
        const int charge = atom->getFormalCharge();
        if (charge == 1) {
            return "+";
        }
        if (charge == -1) {
            return "-";
        }
        return (charge >= 0 ? "+" : "") + std::to_string(charge);
    }
    if (primitive == "v" || primitiveLower == "valence" ||
        primitiveLower == "atomtotalvalence") {
        return "v" + std::to_string(atom->getTotalValence());
    }
    if (primitive == "R" || primitiveLower == "ring" ||
        primitiveLower == "atominnrings" ||
        primitiveLower == "atomringcount" ||
        primitiveLower == "atom_ring_count") {
        return "R" + std::to_string(mol->getRingInfo()->numAtomRings(atom->getIdx()));
    }
    if (primitive == "r" || primitiveLower == "atomminringsize") {
        return "r" + std::to_string(mol->getRingInfo()->minAtomRingSize(atom->getIdx()));
    }
    if (primitive == "a" || primitiveLower == "aromatic" ||
        primitiveLower == "atomisaromatic") {
        return atom->getIsAromatic() ? "a" : "!a";
    }
    if (primitive == "A" || primitiveLower == "aliphatic" ||
        primitiveLower == "atomisaliphatic") {
        return atom->getIsAromatic() ? "!A" : "A";
    }
    if (primitive == "X" || primitiveLower == "atomtotaldegree") {
        return "X" + std::to_string(atom->getTotalDegree());
    }
    if (primitive == "x" || primitiveLower == "atomringbondcount") {
        int ringBondCount = 0;
        for (const auto bond : mol->atomBonds(atom)) {
            if (mol->getRingInfo()->numBondRings(bond->getIdx()) > 0) {
                ++ringBondCount;
            }
        }
        return "x" + std::to_string(ringBondCount);
    }
    if (primitiveLower == "atomic_number" || primitiveLower == "atomicnum" ||
        primitiveLower == "atomatomicnum") {
        return "#" + std::to_string(atom->getAtomicNum());
    }
    if (primitiveLower == "^" || primitiveLower == "hybridization" ||
        primitiveLower == "AtomHybridization") {
              switch (atom->getHybridization()) {
                case RDKit::QueryAtom::S:
                  return "^0";
                case RDKit::QueryAtom::SP:
                  return "^1";
                case RDKit::QueryAtom::SP2:
                  return "^2";
                case RDKit::QueryAtom::SP3:
                  return "^3";
                case RDKit::QueryAtom::SP3D:
                  return "^4";
                case RDKit::QueryAtom::SP3D2:
                  return "^5";
                default:
                  return "";
              }
    }

    throw std::invalid_argument("Unsupported primitive: " + primitive);
}

std::string buildAtomPrimitiveTokenFromTokens(
    const RDKit::Atom *atom, const RDKit::ROMol *mol,
    const std::vector<std::string> &tokens, bool includeAtomMap) {
    std::stringstream ss;
    ss << "[";
    for (size_t i = 0; i < tokens.size(); ++i) {
        if (i > 0) {
            ss << "&";
        }
        ss << primitiveValue(atom, mol, tokens[i]);
    }
    if (includeAtomMap && atom->getAtomMapNum() > 0) {
        ss << ":" << atom->getAtomMapNum();
    }
    ss << "]";
    return ss.str();
}

std::string smilesToSmartsFromNormalizedPrimitives(
    const std::string &smiles, const std::vector<std::string> &tokens) {
    if (smiles.empty()) {
        throw std::invalid_argument("Input SMILES string is empty");
    }
    if (tokens.empty()) {
        throw std::invalid_argument("Primitive list is empty");
    }

    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(smiles));
    if (!mol) {
        throw std::invalid_argument("Invalid SMILES string");
    }

    RDKit::RWMol queryMol;
    for (size_t i = 0; i < mol->getNumAtoms(); ++i) {
        const auto atom = mol->getAtomWithIdx(i);
        const std::string atomToken =
            buildAtomPrimitiveTokenFromTokens(atom, mol.get(), tokens, false);

        std::unique_ptr<RDKit::ROMol> atomQueryMol(RDKit::SmartsToMol(atomToken));
        if (!atomQueryMol || atomQueryMol->getNumAtoms() != 1) {
            throw std::invalid_argument("Failed to build atom query from primitives: " +
                                        atomToken);
        }

        const RDKit::Atom *parsedAtom = atomQueryMol->getAtomWithIdx(0);
        auto *queryAtom = new RDKit::QueryAtom(parsedAtom->getAtomicNum());
        if (parsedAtom->hasQuery()) {
            queryAtom->setQuery(parsedAtom->getQuery()->copy());
        }
        queryMol.addAtom(queryAtom, false, true);
    }

    for (const auto bond : mol->bonds()) {
        queryMol.addBond(bond->getBeginAtomIdx(), bond->getEndAtomIdx(),
                         bond->getBondType());
        if (bond->getIsAromatic()) {
            queryMol.getBondWithIdx(queryMol.getNumBonds() - 1)->setIsAromatic(true);
        }
    }

    return RDKit::MolToSmarts(queryMol);
}

std::vector<std::string> smilesToAtomCenteredSmartsFromNormalizedPrimitives(
    const std::string &smiles, const std::vector<std::string> &tokens,
    unsigned int radius, bool wildcardNeighbors,
    bool includePrimitiveSubsets, bool deduplicate,
    std::unordered_set<std::string> *seen) {
    if (smiles.empty()) {
        throw std::invalid_argument("Input SMILES string is empty");
    }
    if (tokens.empty()) {
        throw std::invalid_argument("Primitive list is empty");
    }

    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(smiles));
    if (!mol) {
        throw std::invalid_argument("Invalid SMILES string");
    }

    std::vector<std::string> out;
    if (!includePrimitiveSubsets) {
        out.reserve(mol->getNumAtoms());
    }
    std::unordered_set<std::string> localSeen;
    if (deduplicate && !seen) {
        localSeen.reserve(mol->getNumAtoms());
        seen = &localSeen;
    }

    auto appendOut = [&](const std::string &value) {
        if (!deduplicate || seen->insert(value).second) {
            out.push_back(value);
        }
    };

    for (size_t centerIdx = 0; centerIdx < mol->getNumAtoms(); ++centerIdx) {
        std::unordered_map<unsigned int, unsigned int> atomDistanceMap;
        const RDKit::PATH_TYPE environmentBonds =
            RDKit::findAtomEnvironmentOfRadiusN(*mol, radius,
                                                static_cast<unsigned int>(centerIdx),
                                                false, false, &atomDistanceMap);

        std::map<int, int> atomIdxMap;
        std::unique_ptr<RDKit::ROMol> envSubmol;
        if (!environmentBonds.empty()) {
            std::vector<int> sortedEnvironmentBonds(environmentBonds.begin(),
                                                    environmentBonds.end());
            std::sort(sortedEnvironmentBonds.begin(), sortedEnvironmentBonds.end());
            envSubmol.reset(
                RDKit::Subgraphs::pathToSubmol(*mol, sortedEnvironmentBonds, true,
                                               atomIdxMap));
        } else {
            auto *single = new RDKit::RWMol();
            const auto centerAtom = mol->getAtomWithIdx(centerIdx);
            auto *queryAtom = new RDKit::QueryAtom(centerAtom->getAtomicNum());
            const std::string centerToken =
                buildAtomPrimitiveTokenFromTokens(centerAtom, mol.get(), tokens,
                                                  false);
            std::unique_ptr<RDKit::ROMol> centerQueryMol(
                RDKit::SmartsToMol(centerToken));
            if (!centerQueryMol || centerQueryMol->getNumAtoms() != 1) {
                throw std::invalid_argument(
                    "Failed to build atom query from primitives: " + centerToken);
            }
            const RDKit::Atom *parsedCenter = centerQueryMol->getAtomWithIdx(0);
            if (parsedCenter->hasQuery()) {
                queryAtom->setQuery(parsedCenter->getQuery()->copy());
            }
            single->addAtom(queryAtom, false, true);
            envSubmol.reset(single);
            atomIdxMap[static_cast<int>(centerIdx)] = 0;
        }

        for (const auto &mapping : atomIdxMap) {
            const int originalIdx = mapping.first;
            const int submolIdx = mapping.second;
            const auto parentAtom = mol->getAtomWithIdx(static_cast<unsigned int>(originalIdx));
            const bool isNeighbor =
                static_cast<unsigned int>(originalIdx) != centerIdx;
            const std::string token =
                (wildcardNeighbors && isNeighbor)
                    ? "[*]"
                    : buildAtomPrimitiveTokenFromTokens(parentAtom, mol.get(), tokens,
                                                        false);
            std::unique_ptr<RDKit::ROMol> atomQueryMol(RDKit::SmartsToMol(token));
            if (!atomQueryMol || atomQueryMol->getNumAtoms() != 1) {
                throw std::invalid_argument(
                    "Failed to build atom query from primitives: " + token);
            }
            const RDKit::Atom *parsedAtom = atomQueryMol->getAtomWithIdx(0);
            auto *queryAtom = static_cast<RDKit::QueryAtom *>(
                envSubmol->getAtomWithIdx(static_cast<unsigned int>(submolIdx)));
            if (parsedAtom->hasQuery()) {
                queryAtom->setQuery(parsedAtom->getQuery()->copy());
            }
        }

        auto centerIt = atomIdxMap.find(static_cast<int>(centerIdx));
        if (centerIt == atomIdxMap.end()) {
            throw std::runtime_error("Center atom not present in atom environment");
        }

        std::vector<std::string> centerPrimitiveValues;
        if (includePrimitiveSubsets) {
            centerPrimitiveValues.reserve(tokens.size());
            const auto centerAtom =
                mol->getAtomWithIdx(static_cast<unsigned int>(centerIdx));
            for (const auto &token : tokens) {
                centerPrimitiveValues.push_back(
                    primitiveValue(centerAtom, mol.get(), token));
            }
            if (centerPrimitiveValues.size() >= 64) {
                throw std::invalid_argument(
                    "Too many primitives for subset expansion (max 63)");
            }
        }

        std::string smartsOut;
        if (wildcardNeighbors) {
            const auto centerAtom = mol->getAtomWithIdx(
                static_cast<unsigned int>(centerIdx));
            const std::string rootToken =
                buildAtomPrimitiveTokenFromTokens(centerAtom, mol.get(), tokens,
                                                  false);
            smartsOut = buildWildcardNeighborSmarts(
                *mol, static_cast<unsigned int>(centerIdx), environmentBonds,
                rootToken);
        } else {
            std::unique_ptr<RDKit::ROMol> orderedEnvSubmol =
                reorderSubmolBondsDeterministically(*envSubmol, atomIdxMap,
                                                    static_cast<unsigned int>(centerIdx));

            RDKit::SmilesWriteParams params;
            params.doIsomericSmiles = true;
            params.allBondsExplicit = true;
            params.canonical = false;
            params.rootedAtAtom = centerIt->second;
            smartsOut = bracketUnbracketedWildcardAtoms(
                RDKit::MolToSmarts(*orderedEnvSubmol, params));
        }

        if (!includePrimitiveSubsets) {
            appendOut(smartsOut);
            continue;
        }

        const size_t rootEnd = rootedAtomTokenEnd(smartsOut);
        if (rootEnd == std::string::npos) {
            appendOut(smartsOut);
            continue;
        }

        const std::string suffix = smartsOut.substr(rootEnd + 1);
        const std::uint64_t fullMask =
            (std::uint64_t{1} << centerPrimitiveValues.size()) - 1;

        // Emit full primitive token first for backwards-friendly ordering,
        // then all remaining non-empty subsets.
        appendOut(buildAtomTokenFromPrimitiveValues(centerPrimitiveValues,
                                                    fullMask) +
                  suffix);
        for (std::uint64_t mask = fullMask - 1; mask > 0; --mask) {
            appendOut(buildAtomTokenFromPrimitiveValues(centerPrimitiveValues,
                                                        mask) +
                      suffix);
        }
    }

    return out;
}

}  // namespace

std::vector<std::string> normalize_primitive_list(
    const std::vector<std::string> &primitiveList) {
    return normalizePrimitiveList(primitiveList);
}

std::string build_atom_primitive_token(
    const RDKit::Atom *atom, const RDKit::ROMol *mol,
    const std::vector<std::string> &primitiveList, bool includeAtomMap) {
    if (!atom || !mol) {
        throw std::invalid_argument("Atom or molecule pointer is null");
    }
    const std::vector<std::string> tokens = normalizePrimitiveList(primitiveList);
    if (tokens.empty()) {
        throw std::invalid_argument("Primitive list is empty");
    }
    return buildAtomPrimitiveTokenFromTokens(atom, mol, tokens, includeAtomMap);
}

std::string smiles_to_smarts(const std::string &smiles,
                             const std::vector<std::string> &primitiveList) {
    const std::vector<std::string> tokens = normalize_primitive_list(primitiveList);
    return smilesToSmartsFromNormalizedPrimitives(smiles, tokens);
}

std::vector<std::string> smiles_to_smarts(
    const std::vector<std::string> &smilesList,
    const std::vector<std::string> &primitiveList) {
    const std::vector<std::string> tokens = normalize_primitive_list(primitiveList);

    std::vector<std::string> out;
    out.reserve(smilesList.size());
    for (const auto &smiles : smilesList) {
        out.push_back(smilesToSmartsFromNormalizedPrimitives(smiles, tokens));
    }
    return out;
}

std::vector<std::string> smiles_to_atom_centered_smarts(
    const std::string &smiles, const std::vector<std::string> &primitiveList,
    unsigned int radius) {
    return smiles_to_atom_centered_smarts(smiles, primitiveList, radius,
                                          false, false, false);
}

std::vector<std::string> smiles_to_atom_centered_smarts(
    const std::string &smiles, const std::vector<std::string> &primitiveList,
    unsigned int radius, bool wildcardNeighbors) {
    return smiles_to_atom_centered_smarts(smiles, primitiveList, radius,
                                          wildcardNeighbors, false, false);
}

std::vector<std::string> smiles_to_atom_centered_smarts(
    const std::string &smiles, const std::vector<std::string> &primitiveList,
    unsigned int radius, bool wildcardNeighbors,
    bool includePrimitiveSubsets, bool deduplicate) {
    const std::vector<std::string> tokens = normalize_primitive_list(primitiveList);
    return smilesToAtomCenteredSmartsFromNormalizedPrimitives(smiles, tokens,
                                                              radius,
                                                              wildcardNeighbors,
                                                              includePrimitiveSubsets,
                                                              deduplicate,
                                                              nullptr);
}

std::vector<std::vector<std::string>> smiles_to_atom_centered_smarts(
    const std::vector<std::string> &smilesList,
    const std::vector<std::string> &primitiveList, unsigned int radius) {
    return smiles_to_atom_centered_smarts(smilesList, primitiveList, radius,
                                          false, false, false);
}

std::vector<std::vector<std::string>> smiles_to_atom_centered_smarts(
    const std::vector<std::string> &smilesList,
    const std::vector<std::string> &primitiveList, unsigned int radius,
    bool wildcardNeighbors) {
    return smiles_to_atom_centered_smarts(smilesList, primitiveList, radius,
                                          wildcardNeighbors, false, false);
}

std::vector<std::vector<std::string>> smiles_to_atom_centered_smarts(
    const std::vector<std::string> &smilesList,
    const std::vector<std::string> &primitiveList, unsigned int radius,
    bool wildcardNeighbors, bool includePrimitiveSubsets,
    bool deduplicate) {
    const std::vector<std::string> tokens = normalize_primitive_list(primitiveList);

    std::vector<std::vector<std::string>> out;
    out.reserve(smilesList.size());
    std::unordered_set<std::string> seen;
    if (deduplicate) {
        seen.reserve(smilesList.size() * 4);
    }
    for (const auto &smiles : smilesList) {
        try {
        out.push_back(smilesToAtomCenteredSmartsFromNormalizedPrimitives(smiles,
                                                                          tokens,
                                                                          radius,
                                                                          wildcardNeighbors,
                                                                          includePrimitiveSubsets,
                                                                          deduplicate,
                                                                          deduplicate ? &seen : nullptr));
        } catch (const std::exception &e) {
            // If there's an error processing a SMILES string, log it and continue with the next one
            std::cerr << "Error processing SMILES '" << smiles << "': " << e.what() << std::endl;
            out.push_back({}); // Add an empty vector for this entry to maintain alignment with input list
        }   
    }
    return out;
}

} // namespace atom_typer

