"""limonene 1,2-monooxygenase.

Catalysis of the reaction: (4R)- or (4S)-limonene + NAD(P)H + O2 + H(+) =
limonene 1,2-epoxide + NAD(P)(+) + H2O.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import (
    CHEBI_H2O,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
)
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_R_LIMONENE = "CHEBI:15382"
CHEBI_S_LIMONENE = "CHEBI:15383"
CHEBI_LIMONENE_12_EPOXIDE = "CHEBI:16431"


class Limonene12Monooxygenase(ReactionClass):
    """limonene 1,2-monooxygenase.

    Catalysis of the reaction: (4R)- or (4S)-limonene + NAD(P)H + O2 + H(+) =
    limonene 1,2-epoxide + NAD(P)(+) + H2O.
    """

    EC_NUMBER_PREFIX = "1.14.13.107"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    SUBSTRATES = {CHEBI_R_LIMONENE, CHEBI_S_LIMONENE}
    REDOX_PAIRS = {
        CHEBI_NADPH: CHEBI_NADP_PLUS,
        CHEBI_NADH: CHEBI_NAD_PLUS,
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = [participant.chebi_id for participant in reaction.left_participants if participant.chebi_id]
        right_ids = [participant.chebi_id for participant in reaction.right_participants if participant.chebi_id]

        if CHEBI_O2 not in left_ids or CHEBI_H2O not in right_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Requires dioxygen substrate and water product",
            )
        if CHEBI_LIMONENE_12_EPOXIDE not in right_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Requires limonene 1,2-epoxide product",
            )

        substrate = next((chebi_id for chebi_id in left_ids if chebi_id in self.SUBSTRATES), None)
        if substrate is None:
            return ClassificationResult(
                is_member=False,
                explanation="Requires (4R)- or (4S)-limonene substrate",
            )

        for reduced, oxidized in self.REDOX_PAIRS.items():
            if reduced in left_ids and oxidized in right_ids:
                return ClassificationResult(
                    is_member=True,
                    explanation="Limonene 1,2-monooxygenase: NAD(P)H-dependent epoxidation of limonene",
                )

        return ClassificationResult(
            is_member=False,
            explanation="Missing the NAD(P)H to NAD(P)(+) cofactor conversion",
        )
