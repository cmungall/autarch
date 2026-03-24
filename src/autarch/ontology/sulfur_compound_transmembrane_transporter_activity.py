"""sulfur compound transmembrane transporter activity."""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass
from autarch.ontology.transport_utils import coupled_transported_pairs, is_sulfur_compound


class SulfurCompoundTransmembraneTransporterActivity(ReactionClass):
    """sulfur compound transmembrane transporter activity."""

    GO_ID = "GO:1901682"
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

        if not any(is_sulfur_compound(left_participant) for left_participant, _ in transported):
            return ClassificationResult(
                is_member=False,
                explanation="No transported sulfur compound detected",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Sulfur compound transmembrane transporter activity",
        )
