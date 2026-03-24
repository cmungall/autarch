"""transposing c=c bonds."""

from autarch.datamodel import Reaction
from autarch.ontology.ec_prefix_aggregate import (
    aggregate_explicit_child_membership,
    aggregate_explicit_child_supports_evaluation,
)
from autarch.ontology.delta3_delta2_enoyl_coa_isomerase import Delta3Delta2EnoylCoAIsomerase
from autarch.ontology.reaction import ReactionClass


class TransposingCCBonds(ReactionClass):
    """transposing c=c bonds."""

    EC_NUMBER_PREFIX = "5.3.3.-"
    CHILD_CLASSES = (
        Delta3Delta2EnoylCoAIsomerase,
    )

    @classmethod
    def supports_evaluation(cls, reaction: Reaction) -> bool:
        return aggregate_explicit_child_supports_evaluation(reaction, cls.CHILD_CLASSES)

    def check_membership_impl(self, reaction: Reaction):
        return aggregate_explicit_child_membership(
            reaction,
            self.CHILD_CLASSES,
            "Transposing C=C bonds",
        )
