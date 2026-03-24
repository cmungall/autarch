"""adenylyltransferase activity."""

from autarch.datamodel import Reaction
from autarch.ontology.ampylase import AMPylase
from autarch.ontology.ec_prefix_aggregate import (
    aggregate_explicit_child_membership,
    aggregate_explicit_child_supports_evaluation,
)
from autarch.ontology.reaction import ReactionClass


class Adenylyltransferase(ReactionClass):
    """adenylyltransferase activity."""

    GO_ID = "GO:0070566"
    CHILD_CLASSES = (AMPylase,)

    @classmethod
    def supports_evaluation(cls, reaction: Reaction) -> bool:
        return aggregate_explicit_child_supports_evaluation(reaction, cls.CHILD_CLASSES)

    def check_membership_impl(self, reaction: Reaction):
        return aggregate_explicit_child_membership(
            reaction,
            self.CHILD_CLASSES,
            "Adenylyltransferase",
        )
