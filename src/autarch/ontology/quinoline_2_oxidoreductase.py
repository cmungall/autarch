"""quinoline 2-oxidoreductase.

Catalysis of quinoline 2-oxidation for supported quinoline substrates using the
characteristic quinoline oxidoreductase carrier placeholders present in the
cached RHEA reactions.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_GENERIC_OXIDANT = "CHEBI:15377"
CHEBI_GENERIC_CARRIER_LEFT = "CHEBI:13193"
CHEBI_GENERIC_CARRIER_RIGHT = "CHEBI:17499"


class Quinoline2Oxidoreductase(ReactionClass):
    """quinoline 2-oxidoreductase.

    Catalysis of quinoline 2-oxidation for supported quinoline substrates using the
    characteristic quinoline oxidoreductase carrier placeholders present in the
    cached RHEA reactions.
    """

    EC_NUMBER_PREFIX = "1.3.99.17"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    SUBSTRATE_TO_PRODUCT = {
        "CHEBI:17362": "CHEBI:18289",
        "CHEBI:48993": "CHEBI:48995",
        "CHEBI:48994": "CHEBI:48996",
        "CHEBI:48980": "CHEBI:48987",
        "CHEBI:48981": "CHEBI:48988",
        "CHEBI:20140": "CHEBI:20114",
        "CHEBI:48983": "CHEBI:48986",
        "CHEBI:48984": "CHEBI:48989",
        "CHEBI:48985": "CHEBI:48990",
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = [
            participant.chebi_id
            for participant in reaction.left_participants
            if participant.chebi_id
            for _ in range(max(1, participant.count))
        ]
        right_ids = [
            participant.chebi_id
            for participant in reaction.right_participants
            if participant.chebi_id
            for _ in range(max(1, participant.count))
        ]

        if CHEBI_GENERIC_CARRIER_LEFT not in left_ids or CHEBI_GENERIC_CARRIER_RIGHT not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires quinoline oxidoreductase carrier placeholders")
        if left_ids.count(CHEBI_GENERIC_OXIDANT) < 2:
            return ClassificationResult(is_member=False, explanation="Requires the duplicated oxidant placeholder pattern from the cached quinoline reactions")

        left_core = [chebi_id for chebi_id in left_ids if chebi_id not in {CHEBI_GENERIC_CARRIER_LEFT, CHEBI_GENERIC_OXIDANT}]
        right_core = [chebi_id for chebi_id in right_ids if chebi_id != CHEBI_GENERIC_CARRIER_RIGHT]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(is_member=False, explanation="Expected one quinoline substrate and one quinolinone product")

        if self.SUBSTRATE_TO_PRODUCT.get(left_core[0]) != right_core[0]:
            return ClassificationResult(is_member=False, explanation="Quinoline branch does not match a supported quinoline 2-oxidoreductase pair")

        return ClassificationResult(is_member=True, explanation="Quinoline 2-oxidoreductase: supported quinoline 2-oxidation branch with cached carrier placeholders")
