"""phosphotransferases (phosphomutases)."""

from autarch.datamodel import Reaction
from autarch.ontology.ec_prefix_aggregate import (
    aggregate_explicit_child_membership,
    aggregate_explicit_child_supports_evaluation,
)
from autarch.ontology.intramolecular_phosphotransferase import IntramolecularPhosphotransferase
from autarch.ontology.reaction import ReactionClass


class PhosphotransferasesPhosphomutases(ReactionClass):
    """phosphotransferases (phosphomutases)."""

    EC_NUMBER_PREFIX = "5.4.2.-"
    CHILD_CLASSES = (
        IntramolecularPhosphotransferase,
    )

    @classmethod
    def supports_evaluation(cls, reaction: Reaction) -> bool:
        return aggregate_explicit_child_supports_evaluation(reaction, cls.CHILD_CLASSES)

    def check_membership_impl(self, reaction: Reaction):
        return aggregate_explicit_child_membership(
            reaction,
            self.CHILD_CLASSES,
            "Phosphotransferases (phosphomutases)",
        )
