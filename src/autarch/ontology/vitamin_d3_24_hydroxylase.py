"""vitamin D3 24-hydroxylase.

Catalysis of side-chain hydroxylation of vitamin D3 derivatives using reduced
adrenodoxin as electron donor.
"""

from typing import Sequence

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS, CHEBI_O2
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_ADRENOXIN_REDUCED = "CHEBI:33738"
CHEBI_ADRENOXIN_OXIDIZED = "CHEBI:33737"


class VitaminD324Hydroxylase(ReactionClass):
    """vitamin D3 24-hydroxylase.

    Catalysis of side-chain hydroxylation of vitamin D3 derivatives using reduced
    adrenodoxin as electron donor.
    """

    EC_NUMBER_PREFIX = "1.14.15.16"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    SUBSTRATE_TO_PRODUCT = {
        "CHEBI:17823": "CHEBI:47799",
        "CHEBI:17933": "CHEBI:28818",
        "CHEBI:47799": "CHEBI:47812",
        "CHEBI:47812": "CHEBI:47813",
        "CHEBI:47818": "CHEBI:47820",
        "CHEBI:47820": "CHEBI:58715",
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = _expand_chebi_ids(reaction.left_participants)
        right_ids = _expand_chebi_ids(reaction.right_participants)

        if left_ids.count(CHEBI_ADRENOXIN_REDUCED) != 2 or right_ids.count(CHEBI_ADRENOXIN_OXIDIZED) != 2:
            return ClassificationResult(is_member=False, explanation="Requires two reduced adrenodoxin donors converted to two oxidized adrenodoxin products")
        if left_ids.count(CHEBI_O2) != 1 or left_ids.count(CHEBI_H_PLUS) not in {1, 2}:
            return ClassificationResult(is_member=False, explanation="Requires dioxygen and hydrons for the hydroxylation step")
        if CHEBI_H2O not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires water product")

        left_core = [chebi_id for chebi_id in left_ids if chebi_id not in {CHEBI_ADRENOXIN_REDUCED, CHEBI_O2, CHEBI_H_PLUS}]
        right_core = [chebi_id for chebi_id in right_ids if chebi_id not in {CHEBI_ADRENOXIN_OXIDIZED, CHEBI_H2O}]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(is_member=False, explanation="Expected one vitamin D substrate and one hydroxylated vitamin D product")

        substrate = left_core[0]
        product = right_core[0]
        if self.SUBSTRATE_TO_PRODUCT.get(substrate) != product:
            return ClassificationResult(is_member=False, explanation="Sterol branch does not match a supported vitamin D3 24-hydroxylase pair")

        return ClassificationResult(is_member=True, explanation="Vitamin D3 24-hydroxylase: adrenodoxin-dependent hydroxylation of a vitamin D3 derivative")


def _expand_chebi_ids(participants: Sequence[Participant]) -> list[str]:
    ids: list[str] = []
    for participant in participants:
        chebi_id = getattr(participant, "chebi_id", None)
        if not chebi_id:
            continue
        count = max(int(getattr(participant, "count", 1) or 1), 1)
        ids.extend([chebi_id] * count)
    return ids
