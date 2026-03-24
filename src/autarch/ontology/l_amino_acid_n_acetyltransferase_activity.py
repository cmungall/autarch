"""L-amino-acid N-acetyltransferase activity."""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.moiety import Moiety
from autarch.molecules import CHEBI_ACETYL_COA, CHEBI_COA, CHEBI_PHOSPHATE
from autarch.ontology.biochemical_context_utils import is_amino_acid_like
from autarch.ontology.lipid_utils import carbon_count
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_ACETYL_PHOSPHATE = "CHEBI:15351"


class LAminoAcidNAcetyltransferaseActivity(ReactionClass):
    """L-amino-acid N-acetyltransferase activity."""

    GO_ID = "GO:0140085"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        forward = self._check_direction(reaction.left_participants, reaction.right_participants)
        if forward.is_member:
            return forward
        reverse = self._check_direction(reaction.right_participants, reaction.left_participants)
        if reverse.is_member:
            return ClassificationResult(is_member=True, explanation=f"{reverse.explanation} (reverse reaction orientation)")
        return forward

    @staticmethod
    def _check_direction(left_participants: list[Participant], right_participants: list[Participant]) -> ClassificationResult:
        left_ids = {participant.chebi_id for participant in left_participants if participant.chebi_id}
        right_ids = {participant.chebi_id for participant in right_participants if participant.chebi_id}
        has_acetyl_donor = CHEBI_ACETYL_COA in left_ids or CHEBI_ACETYL_PHOSPHATE in left_ids
        has_leaving_group = CHEBI_COA in right_ids or CHEBI_PHOSPHATE in right_ids
        if not (has_acetyl_donor and has_leaving_group):
            return ClassificationResult(is_member=False, explanation="Requires an acetyl donor with CoA or phosphate leaving-group release")

        substrates = [participant for participant in left_participants if is_amino_acid_like(participant)]
        products = [participant for participant in right_participants if participant.has_moiety(Moiety.CARBOXYL) and carbon_count(participant) >= 4]
        if not substrates or not products:
            return ClassificationResult(is_member=False, explanation="Requires amino-acid substrate and acetylated organic product")

        for substrate in substrates:
            for product in products:
                if carbon_count(product) in {carbon_count(substrate) + 2, carbon_count(substrate)}:
                    return ClassificationResult(is_member=True, explanation="L-amino-acid N-acetyltransferase activity: acetyl transfer to an amino-acid substrate")

        return ClassificationResult(is_member=False, explanation="No amino-acid acetylation pattern detected")
