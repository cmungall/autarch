"""protein methyltransferase activity.

Catalysis of the transfer of a methyl group to a protein substrate.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_SAH, CHEBI_SAM
from autarch.ontology.biochemical_context_utils import is_protein_like
from autarch.ontology.reaction import IDENTIFIER_AND_POLYMER_FRIENDLY_EVIDENCE, ReactionClass


class ProteinMethyltransferaseActivity(ReactionClass):
    """protein methyltransferase activity.

    Catalysis of the transfer of a methyl group to a protein substrate.
    """

    GO_ID = "GO:0008276"
    EVALUATION_EVIDENCE = IDENTIFIER_AND_POLYMER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = {participant.chebi_id for participant in reaction.left_participants if participant.chebi_id}
        right_ids = {participant.chebi_id for participant in reaction.right_participants if participant.chebi_id}
        if CHEBI_SAM not in left_ids or CHEBI_SAH not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires SAM consumption with SAH production")

        left_proteins = [participant for participant in reaction.left_participants if is_protein_like(participant)]
        right_proteins = [participant for participant in reaction.right_participants if is_protein_like(participant)]
        if not left_proteins or not right_proteins:
            return ClassificationResult(is_member=False, explanation="Requires protein-like participants on both sides of the reaction")

        return ClassificationResult(
            is_member=True,
            explanation="Protein methyltransferase activity: SAM-dependent methyl transfer on a protein-like substrate",
        )
