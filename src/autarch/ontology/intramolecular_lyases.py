"""intramolecular lyases."""

from autarch.datamodel import Reaction
from autarch.ontology.ec_prefix_aggregate import (
    aggregate_explicit_child_membership,
    aggregate_explicit_child_supports_evaluation,
)
from autarch.ontology.intramolecular_lyase import IntramolecularLyase
from autarch.ontology.reaction import ReactionClass


class IntramolecularLyases(ReactionClass):
    """intramolecular lyases."""

    EC_NUMBER_PREFIX = "5.5.1.-"
    CHILD_CLASSES = (
        IntramolecularLyase,
    )

    @classmethod
    def supports_evaluation(cls, reaction: Reaction) -> bool:
        return aggregate_explicit_child_supports_evaluation(reaction, cls.CHILD_CLASSES)

    def check_membership_impl(self, reaction: Reaction):
        return aggregate_explicit_child_membership(
            reaction,
            self.CHILD_CLASSES,
            "Intramolecular lyases",
        )
