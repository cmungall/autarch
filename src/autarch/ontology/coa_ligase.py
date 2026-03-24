"""CoA-ligase activity."""

from autarch.datamodel import Reaction
from autarch.ontology.cholate_coa_ligase import CholateCoALigase
from autarch.ontology.ec_prefix_aggregate import (
    aggregate_explicit_child_membership,
    aggregate_explicit_child_supports_evaluation,
)
from autarch.ontology.reaction import ReactionClass


class CoALigase(ReactionClass):
    """CoA-ligase activity."""

    GO_ID = "GO:0016405"
    CHILD_CLASSES = (CholateCoALigase,)

    @classmethod
    def supports_evaluation(cls, reaction: Reaction) -> bool:
        return aggregate_explicit_child_supports_evaluation(reaction, cls.CHILD_CLASSES)

    def check_membership_impl(self, reaction: Reaction):
        return aggregate_explicit_child_membership(
            reaction,
            self.CHILD_CLASSES,
            "CoA-ligase",
        )
