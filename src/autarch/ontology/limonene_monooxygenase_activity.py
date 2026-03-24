"""limonene monooxygenase activity.

Catalysis of limonene oxidation by incorporation of one atom of molecular
oxygen into a limonene-derived product.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_NADH, CHEBI_NADP_PLUS, CHEBI_NAD_PLUS, CHEBI_O2
from autarch.ontology.lipid_utils import carbon_count, oxygen_count
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_R_LIMONENE = "CHEBI:15382"
CHEBI_S_LIMONENE = "CHEBI:15383"
CHEBI_REDUCED_HEMOPROTEIN_REDUCTASE = "CHEBI:57618"
CHEBI_OXIDIZED_HEMOPROTEIN_REDUCTASE = "CHEBI:58210"
CHEBI_REDUCED_FMN = "CHEBI:57783"
CHEBI_OXIDIZED_FMN = "CHEBI:58349"
SUPPORTED_OXYGENATED_PRODUCTS = {
    "CHEBI:10782",
    "CHEBI:15388",
    "CHEBI:15389",
    "CHEBI:15406",
    "CHEBI:16431",
}


class LimoneneMonooxygenaseActivity(ReactionClass):
    """limonene monooxygenase activity.

    Catalysis of limonene oxidation by incorporation of one atom of molecular
    oxygen into a limonene-derived product.
    """

    GO_ID = "GO:0019113"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        forward = self._check_direction(reaction.left_participants, reaction.right_participants)
        if forward.is_member:
            return forward

        reverse = self._check_direction(reaction.right_participants, reaction.left_participants)
        if reverse.is_member:
            return ClassificationResult(
                is_member=True,
                explanation=f"{reverse.explanation} (reverse reaction orientation)",
            )
        return forward

    def _check_direction(
        self,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> ClassificationResult:
        left_ids = {participant.chebi_id for participant in left_participants if participant.chebi_id}
        right_ids = {participant.chebi_id for participant in right_participants if participant.chebi_id}
        if not ({CHEBI_R_LIMONENE, CHEBI_S_LIMONENE} & left_ids):
            return ClassificationResult(
                is_member=False,
                explanation="Requires limonene substrate",
            )
        if CHEBI_O2 not in left_ids or CHEBI_H2O not in right_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Requires dioxygen consumption with water production",
            )

        if not self._has_supported_cofactor_pair(left_ids, right_ids):
            return ClassificationResult(
                is_member=False,
                explanation="Missing the supported monooxygenase cofactor conversion",
            )

        products = [
            participant
            for participant in right_participants
            if participant.chebi_id not in {
                CHEBI_H2O,
                CHEBI_OXIDIZED_HEMOPROTEIN_REDUCTASE,
                CHEBI_OXIDIZED_FMN,
                CHEBI_NAD_PLUS,
                CHEBI_NADP_PLUS,
            }
        ]
        if any(participant.chebi_id in SUPPORTED_OXYGENATED_PRODUCTS for participant in products):
            return ClassificationResult(
                is_member=True,
                explanation="Limonene monooxygenase activity: limonene monooxygenation with the supported redox cofactor conversion",
            )

        oxygenated_products = [
            participant
            for participant in products
            if participant.chebi_id not in {CHEBI_R_LIMONENE, CHEBI_S_LIMONENE}
            and carbon_count(participant) == 10
            and oxygen_count(participant) >= 1
        ]
        if oxygenated_products:
            return ClassificationResult(
                is_member=True,
                explanation="Limonene monooxygenase activity: limonene monooxygenation to an oxygenated monoterpene product",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No supported oxygenated limonene product detected",
        )

    @staticmethod
    def _has_supported_cofactor_pair(left_ids: set[str], right_ids: set[str]) -> bool:
        return any(
            reduced in left_ids and oxidized in right_ids
            for reduced, oxidized in {
                CHEBI_REDUCED_HEMOPROTEIN_REDUCTASE: CHEBI_OXIDIZED_HEMOPROTEIN_REDUCTASE,
                CHEBI_REDUCED_FMN: CHEBI_OXIDIZED_FMN,
                CHEBI_NADH: CHEBI_NAD_PLUS,
            }.items()
        )
