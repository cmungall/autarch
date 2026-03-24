"""aureusidin synthase activity.

Catalysis of oxidative cyclization of chalcone glucosides to aureusidin or
bracteatin glucosides.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS, CHEBI_O2
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_PENTAHYDROXYCHALCONE_GLUCOSIDE = "CHEBI:77622"
CHEBI_TETRAHYDROXYCHALCONE_GLUCOSIDE = "CHEBI:66906"
CHEBI_AUREUSIDIN_GLUCOSIDE = "CHEBI:66905"
CHEBI_BRACTEATIN_GLUCOSIDE = "CHEBI:66907"


class AureusidinSynthaseActivity(ReactionClass):
    """aureusidin synthase activity.

    Catalysis of oxidative cyclization of chalcone glucosides to aureusidin or
    bracteatin glucosides.
    """

    GO_ID = "GO:0033793"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE
    SUBSTRATE_IDS = {CHEBI_PENTAHYDROXYCHALCONE_GLUCOSIDE, CHEBI_TETRAHYDROXYCHALCONE_GLUCOSIDE}
    PRODUCT_IDS = {CHEBI_AUREUSIDIN_GLUCOSIDE, CHEBI_BRACTEATIN_GLUCOSIDE}

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        forward = self._check_direction(reaction.left_participants, reaction.right_participants)
        if forward.is_member:
            return forward
        reverse = self._check_direction(reaction.right_participants, reaction.left_participants)
        if reverse.is_member:
            return ClassificationResult(is_member=True, explanation=f"{reverse.explanation} (reverse reaction orientation)")
        return forward

    def _check_direction(self, left: list[Participant], right: list[Participant]) -> ClassificationResult:
        left_ids = [p.chebi_id for p in left if p.chebi_id]
        right_ids = [p.chebi_id for p in right if p.chebi_id]
        left_set = set(left_ids)
        right_set = set(right_ids)
        if not (self.SUBSTRATE_IDS & left_set):
            return ClassificationResult(is_member=False, explanation="Requires a supported chalcone glucoside substrate")
        if not (self.PRODUCT_IDS & right_set):
            return ClassificationResult(is_member=False, explanation="Requires aureusidin or bracteatin glucoside product")
        if CHEBI_O2 not in left_set or CHEBI_H2O not in right_set:
            return ClassificationResult(is_member=False, explanation="Requires dioxygen consumption and water formation")
        hydrons = left_ids.count(CHEBI_H_PLUS)
        if hydrons not in {0, 1, 2}:
            return ClassificationResult(is_member=False, explanation="Unexpected hydron stoichiometry for aureusidin synthase chemistry")
        return ClassificationResult(
            is_member=True,
            explanation="aureusidin synthase activity: oxidative cyclization of chalcone glucosides",
        )
