"""carbon-sulfur lyases."""

from autarch.datamodel import Reaction
from autarch.ontology.ec_prefix_aggregate import (
    aggregate_explicit_child_membership,
    aggregate_explicit_child_supports_evaluation,
)
from autarch.ontology.carbon_sulfur_lyase import CarbonSulfurLyase
from autarch.ontology.reaction import ReactionClass


class CarbonSulfurLyases(ReactionClass):
    """carbon-sulfur lyases."""

    EC_NUMBER_PREFIX = "4.4.1.-"
    CHILD_CLASSES = (
        CarbonSulfurLyase,
    )

    @classmethod
    def supports_evaluation(cls, reaction: Reaction) -> bool:
        return aggregate_explicit_child_supports_evaluation(reaction, cls.CHILD_CLASSES)

    def check_membership_impl(self, reaction: Reaction):
        return aggregate_explicit_child_membership(
            reaction,
            self.CHILD_CLASSES,
            "Carbon-sulfur lyases",
        )
