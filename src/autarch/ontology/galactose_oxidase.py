"""galactose oxidase activity.

Catalysis of the reaction: D-galactose + O2 = D-galacto-hexodialdose + H2O2.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H2O2, CHEBI_O2
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_GALACTOSE = "CHEBI:4139"
CHEBI_GALACTO_HEXODIALDOSE = "CHEBI:16222"


class GalactoseOxidase(ReactionClass):
    """galactose oxidase activity.

    Catalysis of the reaction: D-galactose + O2 = D-galacto-hexodialdose + H2O2.
    """

    GO_ID = "GO:0045480"
    EC_NUMBER_PREFIX = "1.1.3.9"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_chebis = {participant.chebi_id for participant in reaction.left_participants if participant.chebi_id}
        right_chebis = {participant.chebi_id for participant in reaction.right_participants if participant.chebi_id}

        if left_chebis != {CHEBI_GALACTOSE, CHEBI_O2}:
            return ClassificationResult(
                is_member=False,
                explanation="Requires D-galactose and dioxygen as the complete substrate set",
            )
        if right_chebis != {CHEBI_GALACTO_HEXODIALDOSE, CHEBI_H2O2}:
            return ClassificationResult(
                is_member=False,
                explanation="Requires D-galacto-hexodialdose and hydrogen peroxide as the complete product set",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Galactose oxidase: oxygen-dependent oxidation of D-galactose to the corresponding dialdose with hydrogen peroxide formation",
        )
