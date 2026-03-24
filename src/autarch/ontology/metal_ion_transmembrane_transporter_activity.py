"""metal ion transmembrane transporter activity."""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass
from autarch.ontology.transport_utils import coupled_transported_pairs, is_metal_ion


class MetalIonTransmembraneTransporterActivity(ReactionClass):
    """metal ion transmembrane transporter activity."""

    GO_ID = "GO:0046873"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        if not reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="No transported substrate between locations",
            )

        transported = coupled_transported_pairs(reaction)
        if not transported:
            return ClassificationResult(
                is_member=False,
                explanation="No transported substrate beyond coupling spectators",
            )

        if not any(is_metal_ion(left_participant) for left_participant, _ in transported):
            return ClassificationResult(
                is_member=False,
                explanation="No transported metal ion detected",
            )

        if not all(
            is_metal_ion(left_participant) or left_participant.chebi_id == "CHEBI:15378"
            for left_participant, _ in transported
        ):
            return ClassificationResult(
                is_member=False,
                explanation="Non-metal co-transported substrate detected",
            )
        return ClassificationResult(
            is_member=True,
            explanation="Metal ion transmembrane transporter activity",
        )
