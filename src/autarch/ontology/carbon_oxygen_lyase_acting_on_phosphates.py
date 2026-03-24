"""carbon-oxygen lyase activity, acting on phosphates."""

from autarch.datamodel import Reaction
from autarch.ontology.ec_prefix_aggregate import (
    aggregate_explicit_child_membership,
    aggregate_explicit_child_supports_evaluation,
)
from autarch.ontology.gamma_humulene_synthase import GammaHumuleneSynthase
from autarch.ontology.pinene_synthase import PineneSynthase
from autarch.ontology.reaction import ReactionClass
from autarch.ontology.terpene_synthase import TerpeneSynthase


class CarbonOxygenLyaseActingOnPhosphates(ReactionClass):
    """carbon-oxygen lyase activity, acting on phosphates."""

    GO_ID = "GO:0016838"
    EC_NUMBER_PREFIX = "4.2.3.-"
    CHILD_CLASSES = (
        TerpeneSynthase,
        PineneSynthase,
        GammaHumuleneSynthase,
    )

    @classmethod
    def supports_evaluation(cls, reaction: Reaction) -> bool:
        return aggregate_explicit_child_supports_evaluation(reaction, cls.CHILD_CLASSES)

    def check_membership_impl(self, reaction: Reaction):
        return aggregate_explicit_child_membership(
            reaction,
            self.CHILD_CLASSES,
            "Carbon-oxygen lyase acting on phosphates",
        )
