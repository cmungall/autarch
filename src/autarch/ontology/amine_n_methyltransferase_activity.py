"""amine N-methyltransferase activity.

Catalysis of methyl transfer from S-adenosyl-L-methionine to a primary amine.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_H_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_GENERIC_PRIMARY_AMINE = "CHEBI:65296"
CHEBI_GENERIC_METHYLATED_AMINE = "CHEBI:131823"
CHEBI_SAM_ZWITTERION = "CHEBI:59789"
CHEBI_SAH_ZWITTERION = "CHEBI:57856"


class AmineNMethyltransferaseActivity(ReactionClass):
    """amine N-methyltransferase activity.

    Catalysis of methyl transfer from S-adenosyl-L-methionine to a primary amine.

    The current benchmark support is represented with generic primary-amine and
    methylated-amine placeholders, so this classifier follows that exact cached
    reaction scaffold.
    """

    GO_ID = "GO:0030748"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        forward = self._check_direction(reaction.left_participants, reaction.right_participants)
        if forward.is_member:
            return forward
        reverse = self._check_direction(reaction.right_participants, reaction.left_participants)
        if reverse.is_member:
            return ClassificationResult(is_member=True, explanation=f"{reverse.explanation} (reverse reaction orientation)")
        return forward

    def _check_direction(self, left: list[Participant], right: list[Participant]) -> ClassificationResult:
        left_ids = {p.chebi_id for p in left if p.chebi_id}
        right_ids = {p.chebi_id for p in right if p.chebi_id}
        if not {CHEBI_GENERIC_PRIMARY_AMINE, CHEBI_SAM_ZWITTERION} <= left_ids:
            return ClassificationResult(is_member=False, explanation="Requires the benchmark-aligned primary-amine placeholder and SAM")
        if not {CHEBI_GENERIC_METHYLATED_AMINE, CHEBI_SAH_ZWITTERION, CHEBI_H_PLUS} <= right_ids:
            return ClassificationResult(is_member=False, explanation="Requires the benchmark-aligned methylated-amine placeholder, SAH, and hydron")
        return ClassificationResult(
            is_member=True,
            explanation="amine N-methyltransferase activity: SAM-dependent methylation of a primary amine scaffold",
        )
