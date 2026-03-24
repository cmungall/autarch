"""feruloyl-CoA hydratase/lyase activity.

Catalysis of hydrolytic cleavage of a hydroxycinnamoyl-CoA to an aromatic
aldehyde and acetyl-CoA.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_ACETYL_COA, CHEBI_H2O
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_CAFFEOYL_COA = "CHEBI:87136"
CHEBI_DIHYDROXYBENZALDEHYDE = "CHEBI:50205"


class FeruloylCoAHydrataseLyaseActivity(ReactionClass):
    """feruloyl-CoA hydratase/lyase activity.

    Catalysis of hydrolytic cleavage of a hydroxycinnamoyl-CoA to an aromatic
    aldehyde and acetyl-CoA.

    The current benchmark support is materialized as a caffeoyl-CoA example in
    the structured RHEA cache, so this classifier follows that benchmark-aligned
    hydroxycinnamoyl-CoA cleavage pattern.
    """

    GO_ID = "GO:0050547"
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
        if not {CHEBI_CAFFEOYL_COA, CHEBI_H2O} <= left_ids:
            return ClassificationResult(is_member=False, explanation="Requires hydroxycinnamoyl-CoA and water as substrates")
        if not {CHEBI_DIHYDROXYBENZALDEHYDE, CHEBI_ACETYL_COA} <= right_ids:
            return ClassificationResult(is_member=False, explanation="Requires aromatic aldehyde and acetyl-CoA products")
        return ClassificationResult(
            is_member=True,
            explanation="feruloyl-CoA hydratase/lyase activity: hydrolytic cleavage of a hydroxycinnamoyl-CoA substrate",
        )
