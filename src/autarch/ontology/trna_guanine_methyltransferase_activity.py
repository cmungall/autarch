"""tRNA (guanine) methyltransferase activity.

Catalysis of the transfer of a methyl group from S-adenosyl-L-methionine to guanine in tRNA.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H_PLUS, CHEBI_SAH, CHEBI_SAM
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_GUANOSINE_IN_TRNA = "CHEBI:74269"
METHYLATED_GUANOSINE_IN_TRNA = {
    "CHEBI:73542",
    "CHEBI:74445",
    "CHEBI:74480",
    "CHEBI:74481",
    "CHEBI:74513",
}


class TRNAGuanineMethyltransferaseActivity(ReactionClass):
    """tRNA (guanine) methyltransferase activity.

    Catalysis of the transfer of a methyl group from S-adenosyl-L-methionine to guanine in tRNA.
    """

    GO_ID = "GO:0016423"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        forward = self._check_direction(
            left_ids=[
                participant.chebi_id
                for participant in reaction.left_participants
                if participant.chebi_id
            ],
            right_ids=[
                participant.chebi_id
                for participant in reaction.right_participants
                if participant.chebi_id
            ],
            donor=CHEBI_SAM,
            coproduct=CHEBI_SAH,
        )
        if forward:
            return ClassificationResult(
                is_member=True,
                explanation="tRNA (guanine) methyltransferase activity: SAM-dependent methylation of guanosine in tRNA",
            )

        reverse = self._check_direction(
            left_ids=[
                participant.chebi_id
                for participant in reaction.right_participants
                if participant.chebi_id
            ],
            right_ids=[
                participant.chebi_id
                for participant in reaction.left_participants
                if participant.chebi_id
            ],
            donor=CHEBI_SAM,
            coproduct=CHEBI_SAH,
        )
        if reverse:
            return ClassificationResult(
                is_member=True,
                explanation="tRNA (guanine) methyltransferase activity: SAM-dependent methylation of guanosine in tRNA (reverse reaction orientation)",
            )

        return ClassificationResult(
            is_member=False,
            explanation="Requires SAM/SAH coupling and the supported guanosine-in-tRNA to methylated-guanosine-in-tRNA branch",
        )

    @staticmethod
    def _check_direction(
        left_ids: list[str],
        right_ids: list[str],
        donor: str,
        coproduct: str,
    ) -> bool:
        if donor not in left_ids or coproduct not in right_ids:
            return False

        left_core = {chebi_id for chebi_id in left_ids if chebi_id not in {donor, CHEBI_H_PLUS}}
        right_core = {chebi_id for chebi_id in right_ids if chebi_id not in {coproduct, CHEBI_H_PLUS}}
        return left_core == {CHEBI_GUANOSINE_IN_TRNA} and bool(
            right_core & METHYLATED_GUANOSINE_IN_TRNA
        )
