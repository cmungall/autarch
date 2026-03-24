"""acting on the ch-ch group of donors with other acceptors."""

from autarch.datamodel import Reaction
from autarch.ontology.ec_prefix_aggregate import (
    aggregate_explicit_child_membership,
    aggregate_explicit_child_supports_evaluation,
)
from autarch.ontology.quinoline_2_oxidoreductase import Quinoline2Oxidoreductase
from autarch.ontology.reaction import ReactionClass


class ActingOnTheCHCHGroupOfDonorsWithOtherAcceptors(ReactionClass):
    """acting on the ch-ch group of donors with other acceptors."""

    EC_NUMBER_PREFIX = "1.3.99.-"
    CHILD_CLASSES = (
        Quinoline2Oxidoreductase,
    )

    @classmethod
    def supports_evaluation(cls, reaction: Reaction) -> bool:
        return aggregate_explicit_child_supports_evaluation(reaction, cls.CHILD_CLASSES)

    def check_membership_impl(self, reaction: Reaction):
        return aggregate_explicit_child_membership(
            reaction,
            self.CHILD_CLASSES,
            "Acting on the CH-CH group of donors with other acceptors",
        )
