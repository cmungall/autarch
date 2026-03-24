"""acyl-CoA desaturase activity.

Catalysis of oxygen-dependent dehydrogenation of a fatty acyl-CoA to a more
unsaturated fatty acyl-CoA.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS, CHEBI_O2
from autarch.ontology.lipid_utils import (
    carbon_count,
    double_bond_count,
    is_fatty_acyl_coa_like,
)
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass


class AcylCoADesaturaseActivity(ReactionClass):
    """acyl-CoA desaturase activity.

    Catalysis of oxygen-dependent dehydrogenation of a fatty acyl-CoA to a
    more unsaturated fatty acyl-CoA.
    """

    GO_ID = "GO:0016215"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        forward = self._check_direction(
            reaction.left_participants,
            reaction.right_participants,
        )
        if forward.is_member:
            return forward
        reverse = self._check_direction(
            reaction.right_participants,
            reaction.left_participants,
        )
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
        if not any(participant.chebi_id == CHEBI_O2 for participant in left_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Requires dioxygen",
            )
        if not any(participant.chebi_id == CHEBI_H2O for participant in right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Requires water production",
            )

        left_core = [
            participant
            for participant in left_participants
            if participant.chebi_id not in {CHEBI_O2, CHEBI_H_PLUS}
        ]
        right_core = [
            participant
            for participant in right_participants
            if participant.chebi_id not in {CHEBI_H2O, CHEBI_H_PLUS}
        ]

        substrates = [participant for participant in left_core if is_fatty_acyl_coa_like(participant)]
        products = [participant for participant in right_core if is_fatty_acyl_coa_like(participant)]
        if len(substrates) != 1 or len(products) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Expected one fatty acyl-CoA substrate and one desaturated fatty acyl-CoA product",
            )

        if carbon_count(substrates[0]) != carbon_count(products[0]):
            return ClassificationResult(
                is_member=False,
                explanation="Acyl-CoA desaturation preserves chain length",
            )
        if double_bond_count(products[0]) != double_bond_count(substrates[0]) + 1:
            return ClassificationResult(
                is_member=False,
                explanation="Desaturase product must introduce exactly one additional C=C bond",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Acyl-CoA desaturase activity: oxygen-dependent introduction of one double bond into a fatty acyl-CoA",
        )
