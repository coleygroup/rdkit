#ifndef REACTION_EXTRACTOR_HPP
#define REACTION_EXTRACTOR_HPP

#include <RDGeneral/export.h>

#include <string>
#include <tuple>
#include <vector>

namespace RDKit {
class ChemicalReaction;
}

namespace atom_typer {

using RadiusTemplate = std::tuple<int, std::string>;

RDKIT_ATOMTYPER_EXPORT std::vector<RadiusTemplate> extract_single_root_template(
    RDKit::ChemicalReaction *rxn, const std::vector<std::string> &primitiveList,
    unsigned maxRadius = 3, bool verbose = false);

RDKIT_ATOMTYPER_EXPORT std::vector<RadiusTemplate> extract_single_root_template(
    const std::string &reactionText,
    const std::vector<std::string> &primitiveList,
    unsigned maxRadius = 3, bool verbose = false, bool useSmiles = true);

}  // namespace atom_typer

#endif  // REACTION_EXTRACTOR_HPP
