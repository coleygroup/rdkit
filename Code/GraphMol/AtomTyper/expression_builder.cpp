#include "expression_builder.hpp"
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/PeriodicTable.h>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>
#include <cctype>
#include <algorithm>

namespace atom_typer {

std::string getElementSymbol(int atomicNum);

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

std::vector<std::string> primitivesForLevel(Level level) {
    if (level == Level::MINIMAL) {
        return {"element"};
    }
    if (level == Level::STANDARD) {
        return {"element", "D"};
    }
    if (level == Level::DETAILED) {
        return {"element", "D", "H"};
    }
    return {"element", "a", "D", "H", "charge", "v", "R", "r", "X", "x"};
}

std::string primitiveValue(const RDKit::Atom *atom,
                           const RDKit::ROMol *mol,
                           const std::string &primitive) {
    const std::string primitiveLower = toLowerCopy(primitive);
    if (primitiveLower == "element" || primitiveLower == "atom" ||
        primitiveLower == "symbol") {
        return getElementSymbol(atom->getAtomicNum());
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

    throw std::invalid_argument("Unsupported primitive: " + primitive);
}

std::string buildAtomPrimitiveToken(const RDKit::Atom *atom,
                                    const RDKit::ROMol *mol,
                                    const std::vector<std::string> &tokens) {
    std::stringstream ss;
    ss << "[";
    for (size_t i = 0; i < tokens.size(); ++i) {
        if (i > 0) {
            ss << ";";
        }
        ss << primitiveValue(atom, mol, tokens[i]);
    }
    ss << "]";
    return ss.str();
}

}  // namespace

// Helper function to get element symbol from atomic number
std::string getElementSymbol(int atomicNum) {
    return RDKit::PeriodicTable::getTable()->getElementSymbol(atomicNum);
}

/**
 * Convert a query molecule to a SMARTS string based on detail level
 * 
 * @param mol Original molecule (for accessing atom properties)
 * @param queryMol Query molecule built with appropriate query features
 * @param level Detail level for formatting the SMARTS string
 * @return Formatted SMARTS string
 */
std::string queryMoleculeToSmarts(const RDKit::ROMol* mol, 
                                   const RDKit::RWMol* queryMol,
                                   Level level) {
    std::stringstream ss;
    
    if (level == Level::MINIMAL) {
        // Format: [C][C][O]
        for (size_t i = 0; i < mol->getNumAtoms(); ++i) {
            const auto atom = mol->getAtomWithIdx(i);
            ss << "[" << getElementSymbol(atom->getAtomicNum()) << "]";
        }
    }
    else if (level == Level::STANDARD) {
        // Format: [C;D1][C;D2][O;D1]
        for (size_t i = 0; i < mol->getNumAtoms(); ++i) {
            const auto atom = mol->getAtomWithIdx(i);
            ss << "[" << getElementSymbol(atom->getAtomicNum()) 
               << ";D" << atom->getDegree() << "]";
        }
    }
    else if (level == Level::DETAILED) {
        // Format: [C;D1;H3][C;D2;H2][O;D1;H1]
        for (size_t i = 0; i < mol->getNumAtoms(); ++i) {
            const auto atom = mol->getAtomWithIdx(i);
            ss << "[" << getElementSymbol(atom->getAtomicNum())
               << ";D" << atom->getDegree()
               << ";H" << atom->getTotalNumHs(true) << "]";
        }
    }
    else if (level == Level::COMPLETE) {
        // Use RDKit's MolToSmarts for complete level
        return RDKit::MolToSmarts(*queryMol);
    }
    
    return ss.str();
}

/**
 * Convert a SMILES string to a SMARTS pattern with specified detail level
 * 
 * @param smiles Input SMILES string
 * @param level Detail level for the SMARTS pattern
 * @return SMARTS string representation
 */
std::string smiles_to_smarts(const std::string& smiles, Level level) {
    if (smiles.empty()) {
        throw std::invalid_argument("Input SMILES string is empty");
    }

    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(smiles));
    if (!mol) {
        throw std::invalid_argument("Invalid SMILES string");
    }

    std::stringstream ss;
    if (level == Level::MINIMAL) {
        for (size_t i = 0; i < mol->getNumAtoms(); ++i) {
            const auto atom = mol->getAtomWithIdx(i);
            ss << "[" << getElementSymbol(atom->getAtomicNum()) << "]";
        }
        return ss.str();
    }
    if (level == Level::STANDARD) {
        for (size_t i = 0; i < mol->getNumAtoms(); ++i) {
            const auto atom = mol->getAtomWithIdx(i);
            ss << "[" << getElementSymbol(atom->getAtomicNum())
               << ";D" << atom->getDegree() << "]";
        }
        return ss.str();
    }
    if (level == Level::DETAILED) {
        for (size_t i = 0; i < mol->getNumAtoms(); ++i) {
            const auto atom = mol->getAtomWithIdx(i);
            ss << "[" << getElementSymbol(atom->getAtomicNum())
               << ";D" << atom->getDegree()
               << ";H" << atom->getTotalNumHs(true) << "]";
        }
        return ss.str();
    }

    return smiles_to_smarts(smiles, primitivesForLevel(level));
}

std::string smiles_to_smarts(const std::string &smiles,
                             const std::vector<std::string> &primitiveList) {
    if (smiles.empty()) {
        throw std::invalid_argument("Input SMILES string is empty");
    }

    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(smiles));
    if (!mol) {
        throw std::invalid_argument("Invalid SMILES string");
    }

    const std::vector<std::string> tokens = normalizePrimitiveList(primitiveList);
    if (tokens.empty()) {
        throw std::invalid_argument("Primitive list is empty");
    }

    RDKit::RWMol queryMol;
    for (size_t i = 0; i < mol->getNumAtoms(); ++i) {
        const auto atom = mol->getAtomWithIdx(i);
        const std::string atomToken = buildAtomPrimitiveToken(atom, mol.get(), tokens);

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

} // namespace atom_typer

