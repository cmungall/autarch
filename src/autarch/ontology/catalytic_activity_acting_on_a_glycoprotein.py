"""catalytic activity, acting on a glycoprotein."""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.biochemical_context_utils import is_protein_like
from autarch.ontology.glycosyltransferase import Glycosyltransferase
from autarch.ontology.reaction import IDENTIFIER_AND_POLYMER_FRIENDLY_EVIDENCE, ReactionClass


class CatalyticActivityActingOnAGlycoprotein(ReactionClass):
    """catalytic activity, acting on a glycoprotein."""

    GO_ID = "GO:0140103"
    EVALUATION_EVIDENCE = IDENTIFIER_AND_POLYMER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_proteins = [participant for participant in reaction.left_participants if is_protein_like(participant)]
        right_proteins = [participant for participant in reaction.right_participants if is_protein_like(participant)]
        if not left_proteins or not right_proteins:
            return ClassificationResult(is_member=False, explanation="Requires protein-like participants on both sides of the reaction")

        glyco = Glycosyltransferase().check_membership_impl(reaction)
        if not glyco.is_member:
            return ClassificationResult(is_member=False, explanation="No glycoprotein-directed glycosyl transfer pattern detected")

        return ClassificationResult(is_member=True, explanation="Catalytic activity acting on a glycoprotein: glycosyl transfer on a protein-like substrate")
