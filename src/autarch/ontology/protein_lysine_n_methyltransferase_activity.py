"""protein-lysine N-methyltransferase activity.

Catalysis of the transfer of a methyl group from S-adenosyl-L-methionine to the epsilon-amino group of a lysine residue in a protein substrate.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H_PLUS, CHEBI_SAH, CHEBI_SAM
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_PROTEIN_LYSINE = "CHEBI:29969"
CHEBI_METHYL_PROTEIN_LYSINE = "CHEBI:61929"


class _ProteinLysineMethyltransferaseBase(ReactionClass):
    """Shared SAM-dependent protein-lysine methyltransferase scaffold."""

    EXCLUDE_FROM_DISCOVERY = True
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE
    CONCEPT_PHRASE = "protein-lysine N-methyltransferase activity"

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
            substrate=CHEBI_PROTEIN_LYSINE,
            product=CHEBI_METHYL_PROTEIN_LYSINE,
        )
        if forward:
            return ClassificationResult(
                is_member=True,
                explanation=f"{self.CONCEPT_PHRASE}: SAM-dependent methyl transfer to a protein lysine residue",
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
            substrate=CHEBI_PROTEIN_LYSINE,
            product=CHEBI_METHYL_PROTEIN_LYSINE,
        )
        if reverse:
            return ClassificationResult(
                is_member=True,
                explanation=f"{self.CONCEPT_PHRASE}: SAM-dependent methyl transfer to a protein lysine residue (reverse reaction orientation)",
            )

        return ClassificationResult(
            is_member=False,
            explanation="Requires SAM/SAH coupling and the supported protein lysine to methyl-lysine scaffold",
        )

    @staticmethod
    def _check_direction(
        left_ids: list[str],
        right_ids: list[str],
        donor: str,
        coproduct: str,
        substrate: str,
        product: str,
    ) -> bool:
        if donor not in left_ids or coproduct not in right_ids:
            return False
        left_core = [chebi_id for chebi_id in left_ids if chebi_id not in {donor, CHEBI_H_PLUS}]
        right_core = [chebi_id for chebi_id in right_ids if chebi_id not in {coproduct, CHEBI_H_PLUS}]
        return left_core == [substrate] and right_core == [product]


class ProteinLysineNMethyltransferaseActivity(_ProteinLysineMethyltransferaseBase):
    """protein-lysine N-methyltransferase activity.

    Catalysis of the transfer of a methyl group from S-adenosyl-L-methionine to the epsilon-amino group of a lysine residue in a protein substrate.
    """

    GO_ID = "GO:0016279"
    CONCEPT_PHRASE = "protein-lysine N-methyltransferase activity"
