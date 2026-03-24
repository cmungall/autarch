"""carnosine synthase activity.

Catalysis of the reaction: beta-alanine + L-histidine + ATP = carnosine + ADP + phosphate + H+.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_ADP, CHEBI_ATP, CHEBI_H_PLUS, CHEBI_PHOSPHATE
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_BETA_ALANINE = "CHEBI:57966"
CHEBI_HISTIDINE = "CHEBI:57595"
CHEBI_CARNOSINE = "CHEBI:57485"
PHOSPHATE_CHEBIS = {CHEBI_PHOSPHATE, "CHEBI:16838"}


class CarnosineSynthase(ReactionClass):
    """carnosine synthase activity.

    Catalysis of the reaction: beta-alanine + L-histidine + ATP = carnosine + ADP + phosphate + H+.
    """

    GO_ID = "GO:0047730"
    EC_NUMBER_PREFIX = "6.3.2.11"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        if reaction.is_transport_reaction():
            return ClassificationResult(is_member=False, explanation="Transport reaction - not carnosine synthase")

        left_chebis = {p.chebi_id for p in reaction.left_participants if p.chebi_id}
        right_chebis = {p.chebi_id for p in reaction.right_participants if p.chebi_id}

        required_left = {CHEBI_BETA_ALANINE, CHEBI_HISTIDINE, CHEBI_ATP}
        required_right = {CHEBI_CARNOSINE, CHEBI_ADP}
        if not required_left.issubset(left_chebis):
            return ClassificationResult(is_member=False, explanation="Missing beta-alanine, histidine, or ATP substrate")
        if not required_right.issubset(right_chebis) or not (right_chebis & PHOSPHATE_CHEBIS):
            return ClassificationResult(is_member=False, explanation="Missing carnosine, ADP, or phosphate product")

        substantive_left = {chebi for chebi in left_chebis if chebi not in {CHEBI_H_PLUS}}
        substantive_right = {chebi for chebi in right_chebis if chebi not in {CHEBI_H_PLUS}}
        if substantive_left - required_left or substantive_right - (required_right | PHOSPHATE_CHEBIS):
            return ClassificationResult(
                is_member=False,
                explanation="Contains additional substantive participants outside carnosine synthase chemistry",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Carnosine synthase: ATP-dependent ligation of beta-alanine and histidine to form carnosine",
        )
