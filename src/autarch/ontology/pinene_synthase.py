"""pinene synthase activity.

Catalysis of the reaction: geranyl diphosphate = pinene + diphosphate. This reaction can produce (1R,5R)-alpha-pinene, (1S,5S)-alpha-pinene, (1R,5R)-beta-pinene and (1S,5S)-beta-pinene.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_DIPHOSPHATE
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE
from autarch.ontology.reaction import ReactionClass

CHEBI_GPP = "CHEBI:58057"
PINENE_PRODUCTS = {
    "CHEBI:28660",
    "CHEBI:28261",
    "CHEBI:28359",
    "CHEBI:50026",
}


class PineneSynthase(ReactionClass):
    """pinene synthase activity.

    Catalysis of the reaction: geranyl diphosphate = pinene + diphosphate. This reaction can produce (1R,5R)-alpha-pinene, (1S,5S)-alpha-pinene, (1R,5R)-beta-pinene and (1S,5S)-beta-pinene.
    """

    GO_ID = "GO:0050550"
    EC_NUMBER_PREFIX = "4.2.3.-"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        if reaction.is_transport_reaction():
            return ClassificationResult(is_member=False, explanation="Transport reaction - not pinene synthase")

        left_chebis = {p.chebi_id for p in reaction.left_participants if p.chebi_id}
        right_chebis = {p.chebi_id for p in reaction.right_participants if p.chebi_id}

        if left_chebis != {CHEBI_GPP}:
            return ClassificationResult(
                is_member=False,
                explanation="Pinene synthase requires geranyl diphosphate as the sole substantive substrate",
            )

        pinene_products = right_chebis & PINENE_PRODUCTS
        if len(pinene_products) != 1 or CHEBI_DIPHOSPHATE not in right_chebis:
            return ClassificationResult(
                is_member=False,
                explanation="Missing a pinene product or diphosphate release",
            )

        if right_chebis != pinene_products | {CHEBI_DIPHOSPHATE}:
            return ClassificationResult(
                is_member=False,
                explanation="Contains additional substantive products outside pinene synthase chemistry",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Pinene synthase: cyclization of geranyl diphosphate to a pinene with diphosphate release",
        )
