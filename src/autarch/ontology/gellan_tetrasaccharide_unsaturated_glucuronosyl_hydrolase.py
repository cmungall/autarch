"""gellan tetrasaccharide unsaturated glucuronosyl hydrolase.

Catalysis of hydrolytic cleavage of supported unsaturated glucuronosyl
oligosaccharides to a glycan product and 5-dehydro-4-deoxy-glucuronate.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H2O
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass


class GellanTetrasaccharideUnsaturatedGlucuronosylHydrolase(ReactionClass):
    """gellan tetrasaccharide unsaturated glucuronosyl hydrolase.

    Catalysis of hydrolytic cleavage of supported unsaturated glucuronosyl
    oligosaccharides to a glycan product and 5-dehydro-4-deoxy-glucuronate.
    """

    EC_NUMBER_PREFIX = "3.2.1.179"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    SUBSTRATE_TO_PRODUCTS = {
        "CHEBI:134390": {"CHEBI:134389", "CHEBI:17117"},
        "CHEBI:63263": {"CHEBI:17117", "CHEBI:28497"},
        "CHEBI:63274": {"CHEBI:63278", "CHEBI:28497"},
        "CHEBI:63280": {"CHEBI:17117", "CHEBI:28497"},
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = [participant.chebi_id for participant in reaction.left_participants if participant.chebi_id]
        right_ids = [participant.chebi_id for participant in reaction.right_participants if participant.chebi_id]

        if CHEBI_H2O not in left_ids:
            return ClassificationResult(is_member=False, explanation="Requires water substrate")

        left_core = [chebi_id for chebi_id in left_ids if chebi_id != CHEBI_H2O]
        if len(left_core) != 1 or len(right_ids) != 2:
            return ClassificationResult(is_member=False, explanation="Expected one glycan substrate and two products")

        substrate = left_core[0]
        if set(right_ids) != self.SUBSTRATE_TO_PRODUCTS.get(substrate, set()):
            return ClassificationResult(is_member=False, explanation="Products do not match a supported unsaturated glucuronosyl hydrolase branch")

        return ClassificationResult(is_member=True, explanation="Gellan tetrasaccharide unsaturated glucuronosyl hydrolase: hydrolytic cleavage of a supported unsaturated glucuronosyl glycan")
