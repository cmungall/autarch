"""Ribulose-bisphosphate carboxylase reaction classifier."""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_CO2
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_RIBULOSE_15_BISPHOSPHATE = "CHEBI:57870"
CHEBI_3_PHOSPHOGLYCERATE = "CHEBI:58272"


def _count_chebi(participants, chebi_id: str) -> int:
    """Count total stoichiometry for a ChEBI ID."""
    return sum(p.count for p in participants if p.chebi_id == chebi_id)


class RibuloseBisphosphateCarboxylase(ReactionClass):
    """ribulose-bisphosphate carboxylase"""

    GO_ID = "GO:0016984"  # ribulose-bisphosphate carboxylase activity
    EC_NUMBER_PREFIX = "4.1.1.39"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for the RuBisCO carboxylation/decarboxylation signature."""
        left_ribulose_bp = _count_chebi(
            reaction.left_participants, CHEBI_RIBULOSE_15_BISPHOSPHATE
        )
        right_ribulose_bp = _count_chebi(
            reaction.right_participants, CHEBI_RIBULOSE_15_BISPHOSPHATE
        )
        left_co2 = _count_chebi(reaction.left_participants, CHEBI_CO2)
        right_co2 = _count_chebi(reaction.right_participants, CHEBI_CO2)
        left_pga = _count_chebi(reaction.left_participants, CHEBI_3_PHOSPHOGLYCERATE)
        right_pga = _count_chebi(reaction.right_participants, CHEBI_3_PHOSPHOGLYCERATE)

        forward = left_ribulose_bp >= 1 and left_co2 >= 1 and right_pga >= 2
        reverse = right_ribulose_bp >= 1 and right_co2 >= 1 and left_pga >= 2

        if forward:
            return ClassificationResult(
                is_member=True,
                explanation="RuBisCO carboxylation: ribulose-1,5-bisphosphate + CO2 -> 2 3-phosphoglycerate",
            )
        if reverse:
            return ClassificationResult(
                is_member=True,
                explanation="RuBisCO reverse direction detected: 2 3-phosphoglycerate -> ribulose-1,5-bisphosphate + CO2",
            )

        return ClassificationResult(
            is_member=False,
            explanation="Missing ribulose-bisphosphate carboxylase substrate/product signature",
        )
