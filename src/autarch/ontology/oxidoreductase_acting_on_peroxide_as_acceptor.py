"""oxidoreductase acting on peroxide as acceptor.

Catalysis of an oxidation-reduction (redox) reaction in which peroxide acts as
a hydrogen or electron acceptor.
"""

from autarch.datamodel import Reaction
from autarch.ontology.ec_prefix_aggregate import (
    aggregate_explicit_child_membership,
    aggregate_explicit_child_supports_evaluation,
)
from autarch.ontology.peroxidase import Peroxidase
from autarch.ontology.reaction import ReactionClass


class OxidoreductaseActingOnPeroxideAsAcceptor(ReactionClass):
    """oxidoreductase acting on peroxide as acceptor.

    Catalysis of an oxidation-reduction (redox) reaction in which peroxide acts
    as a hydrogen or electron acceptor.
    """

    GO_ID = "GO:0016684"
    EC_NUMBER_PREFIX = "1.11.-.-"
    CHILD_CLASSES = (Peroxidase,)

    @classmethod
    def supports_evaluation(cls, reaction: Reaction) -> bool:
        return aggregate_explicit_child_supports_evaluation(reaction, cls.CHILD_CLASSES)

    def check_membership_impl(self, reaction: Reaction):
        return aggregate_explicit_child_membership(
            reaction,
            self.CHILD_CLASSES,
            "Oxidoreductase acting on peroxide as acceptor",
        )
