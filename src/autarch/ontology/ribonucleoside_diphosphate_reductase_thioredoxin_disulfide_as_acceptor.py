"""ribonucleoside-diphosphate reductase activity, thioredoxin disulfide as acceptor.

Catalysis of the reduction of ribonucleoside diphosphates to the corresponding
2'-deoxyribonucleoside diphosphates with thioredoxin as reductant.
"""

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_H2O
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_THIOREDOXIN_DISULFIDE = "CHEBI:50058"
CHEBI_THIOREDOXIN_DITHIOL = "CHEBI:57930"
CHEBI_GENERIC_DEOXY_NDP = "CHEBI:73316"
CHEBI_GENERIC_NDP = "CHEBI:29950"


class RibonucleosideDiphosphateReductaseThioredoxinDisulfideAsAcceptor(ReactionClass):
    """ribonucleoside-diphosphate reductase activity, thioredoxin disulfide as acceptor.

    Catalysis of the reduction of ribonucleoside diphosphates to the corresponding
    2'-deoxyribonucleoside diphosphates with thioredoxin as reductant.
    """

    GO_ID = "GO:0004748"
    EC_NUMBER_PREFIX = "1.17.4.1"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    SUBSTRATE_TO_PRODUCT = {
        CHEBI_GENERIC_DEOXY_NDP: CHEBI_GENERIC_NDP,
        "CHEBI:60471": "CHEBI:58223",
        "CHEBI:58595": "CHEBI:58189",
        "CHEBI:57667": "CHEBI:456216",
        "CHEBI:58593": "CHEBI:58069",
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        forward = self._check_direction(reaction.left_participants, reaction.right_participants)
        if forward.is_member:
            return forward

        reverse = self._check_direction(reaction.right_participants, reaction.left_participants)
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
        left_ids = [participant.chebi_id for participant in left_participants if participant.chebi_id]
        right_ids = _expand_chebi_ids(right_participants)
        if CHEBI_THIOREDOXIN_DISULFIDE not in left_ids or CHEBI_H2O not in left_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Requires thioredoxin disulfide and water substrates",
            )
        thioredoxin_products = _extract_thioredoxin_products(right_ids)
        if thioredoxin_products is None:
            return ClassificationResult(
                is_member=False,
                explanation="Requires reduced thioredoxin products",
            )

        left_core = [
            chebi_id
            for chebi_id in left_ids
            if chebi_id not in {CHEBI_THIOREDOXIN_DISULFIDE, CHEBI_H2O}
        ]
        right_core = right_ids.copy()
        for chebi_id in thioredoxin_products:
            right_core.remove(chebi_id)
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Expected one nucleotide diphosphate substrate and one product",
            )

        substrate = left_core[0]
        product = right_core[0]
        if self.SUBSTRATE_TO_PRODUCT.get(substrate) != product:
            return ClassificationResult(
                is_member=False,
                explanation="Nucleotide branch does not match a supported ribonucleoside-diphosphate reductase pair",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Ribonucleoside-diphosphate reductase: thioredoxin-dependent reduction of a nucleotide diphosphate",
        )


def _expand_chebi_ids(participants: list[Participant]) -> list[str]:
    ids: list[str] = []
    for participant in participants:
        if participant.chebi_id:
            ids.extend([participant.chebi_id] * max(participant.count, 1))
    return ids


def _extract_thioredoxin_products(right_ids: list[str]) -> list[str] | None:
    if CHEBI_THIOREDOXIN_DITHIOL in right_ids and CHEBI_GENERIC_NDP in right_ids:
        return [CHEBI_THIOREDOXIN_DITHIOL, CHEBI_GENERIC_NDP]
    if right_ids.count(CHEBI_GENERIC_NDP) >= 2:
        return [CHEBI_GENERIC_NDP, CHEBI_GENERIC_NDP]
    return None
