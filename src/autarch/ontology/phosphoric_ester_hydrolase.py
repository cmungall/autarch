"""phosphoric ester hydrolase activity."""

from autarch.datamodel import Reaction
from autarch.ontology.diphosphoric_monoester_hydrolase import DiphosphoricMonoesterHydrolase
from autarch.ontology.ec_prefix_aggregate import (
    aggregate_explicit_child_membership,
    aggregate_explicit_child_supports_evaluation,
)
from autarch.ontology.phosphatase import Phosphatase
from autarch.ontology.phosphoric_diester_hydrolase import PhosphoricDiesterHydrolase
from autarch.ontology.reaction import ReactionClass


class PhosphoricEsterHydrolase(ReactionClass):
    """phosphoric ester hydrolase activity."""

    GO_ID = "GO:0042578"
    CHILD_CLASSES = (
        Phosphatase,
        PhosphoricDiesterHydrolase,
        DiphosphoricMonoesterHydrolase,
    )

    @classmethod
    def supports_evaluation(cls, reaction: Reaction) -> bool:
        return aggregate_explicit_child_supports_evaluation(reaction, cls.CHILD_CLASSES)

    def check_membership_impl(self, reaction: Reaction):
        return aggregate_explicit_child_membership(
            reaction,
            self.CHILD_CLASSES,
            "Phosphoric ester hydrolase",
        )
