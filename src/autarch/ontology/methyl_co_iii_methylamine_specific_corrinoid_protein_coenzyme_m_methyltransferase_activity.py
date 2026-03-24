"""[methyl-Co(III) methylamine-specific corrinoid protein]:coenzyme M methyltransferase activity.

Catalysis of methyl transfer from methyl-Co(III) corrinoid protein to coenzyme M.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_H_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_METHYL_CORRINOID = "CHEBI:85035"
CHEBI_CORRINOID_CO_I = "CHEBI:85033"
CHEBI_COENZYME_M = "CHEBI:58319"
CHEBI_METHYL_COM = "CHEBI:58286"


class MethylCoIIIMethylamineSpecificCorrinoidProteinCoenzymeMMethyltransferaseActivity(ReactionClass):
    """[methyl-Co(III) methylamine-specific corrinoid protein]:coenzyme M methyltransferase activity.

    Catalysis of methyl transfer from methyl-Co(III) corrinoid protein to coenzyme M.

    The current structured RHEA cache materializes the methylamine-, dimethylamine-,
    and trimethylamine-specific corrinoid carriers with the same ChEBI identifiers,
    so this classifier follows the benchmark-aligned generic corrinoid pattern.
    """

    GO_ID = "GO:0043833"
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
        if not {CHEBI_METHYL_CORRINOID, CHEBI_COENZYME_M} <= left_ids:
            return ClassificationResult(is_member=False, explanation="Requires methyl-Co(III) corrinoid protein and coenzyme M as substrates")
        if not {CHEBI_CORRINOID_CO_I, CHEBI_METHYL_COM, CHEBI_H_PLUS} <= right_ids:
            return ClassificationResult(is_member=False, explanation="Requires Co(I) corrinoid protein, methyl-coenzyme M, and hydron as products")
        return ClassificationResult(
            is_member=True,
            explanation="corrinoid-protein:coenzyme M methyltransferase activity: methyl transfer from corrinoid protein to coenzyme M",
        )
