"""homomethionine N-monooxygenase.

Catalysis of oxidation of homomethionine derivatives to the corresponding oximes
with reduced flavoprotein and dioxygen.
"""

from typing import Sequence

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_CO2, CHEBI_H2O, CHEBI_H_PLUS, CHEBI_O2
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_REDUCED_HEMOPROTEIN_REDUCTASE = "CHEBI:57618"
CHEBI_OXIDIZED_HEMOPROTEIN_REDUCTASE = "CHEBI:58210"


class HomomethionineNMonooxygenase(ReactionClass):
    """homomethionine N-monooxygenase.

    Catalysis of oxidation of homomethionine derivatives to the corresponding oximes
    with reduced flavoprotein and dioxygen.
    """

    EC_NUMBER_PREFIX = "1.14.14.42"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    SUBSTRATE_TO_PRODUCT = {
        "CHEBI:134632": "CHEBI:134682",
        "CHEBI:134633": "CHEBI:134681",
        "CHEBI:134634": "CHEBI:134683",
        "CHEBI:134635": "CHEBI:134684",
        "CHEBI:134636": "CHEBI:134685",
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = _expand_chebi_ids(reaction.left_participants)
        right_ids = _expand_chebi_ids(reaction.right_participants)

        if left_ids.count(CHEBI_REDUCED_HEMOPROTEIN_REDUCTASE) != 2 or right_ids.count(CHEBI_OXIDIZED_HEMOPROTEIN_REDUCTASE) != 2:
            return ClassificationResult(is_member=False, explanation="Requires two reduced flavoprotein cofactors converted to two oxidized products")
        if left_ids.count(CHEBI_O2) != 2 or right_ids.count(CHEBI_H2O) != 3 or CHEBI_CO2 not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires dioxygen consumption with carbon dioxide and three water products")
        if right_ids.count(CHEBI_H_PLUS) != 2:
            return ClassificationResult(is_member=False, explanation="Requires two hydron products")

        left_core = [chebi_id for chebi_id in left_ids if chebi_id not in {CHEBI_REDUCED_HEMOPROTEIN_REDUCTASE, CHEBI_O2}]
        right_core = [chebi_id for chebi_id in right_ids if chebi_id not in {CHEBI_OXIDIZED_HEMOPROTEIN_REDUCTASE, CHEBI_CO2, CHEBI_H2O, CHEBI_H_PLUS}]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(is_member=False, explanation="Expected one homomethionine substrate and one oxime product")

        substrate = left_core[0]
        product = right_core[0]
        if self.SUBSTRATE_TO_PRODUCT.get(substrate) != product:
            return ClassificationResult(is_member=False, explanation="Organic branch does not match a supported homomethionine N-monooxygenase pair")

        return ClassificationResult(is_member=True, explanation="Homomethionine N-monooxygenase: flavoprotein-dependent oxidation of a homomethionine derivative to the corresponding oxime")


def _expand_chebi_ids(participants: Sequence[Participant]) -> list[str]:
    ids: list[str] = []
    for participant in participants:
        chebi_id = getattr(participant, "chebi_id", None)
        if not chebi_id:
            continue
        count = max(int(getattr(participant, "count", 1) or 1), 1)
        ids.extend([chebi_id] * count)
    return ids
