"""transketolases and transaldolases."""

from autarch.datamodel import Reaction
from autarch.ontology.ec_prefix_aggregate import (
    aggregate_explicit_child_membership,
    aggregate_explicit_child_supports_evaluation,
)
from autarch.ontology.aldehyde_ketone_transferase import AldehydeKetoneTransferase
from autarch.ontology.reaction import ReactionClass


class TransketolasesAndTransaldolases(ReactionClass):
    """transketolases and transaldolases."""

    EC_NUMBER_PREFIX = "2.2.1.-"
    CHILD_CLASSES = (
        AldehydeKetoneTransferase,
    )

    @classmethod
    def supports_evaluation(cls, reaction: Reaction) -> bool:
        return aggregate_explicit_child_supports_evaluation(reaction, cls.CHILD_CLASSES)

    def check_membership_impl(self, reaction: Reaction):
        return aggregate_explicit_child_membership(
            reaction,
            self.CHILD_CLASSES,
            "Transketolases and transaldolases",
        )
