"""O-methyltransferase activity."""

from autarch.datamodel import Reaction
from autarch.ontology.ec_prefix_aggregate import (
    aggregate_explicit_child_membership,
    aggregate_explicit_child_supports_evaluation,
)
from autarch.ontology.gibberellin_a4_carboxyl_methyltransferase import (
    GibberellinA4CarboxylMethyltransferase,
)
from autarch.ontology.gibberellin_a9_o_methyltransferase import (
    GibberellinA9OMethyltransferase,
)
from autarch.ontology.reaction import ReactionClass
from autarch.ontology.ten_hydroxydihydrosanguinarine_10_o_methyltransferase import (
    TenHydroxydihydrosanguinarine10OMethyltransferase,
)
from autarch.ontology.tocopherol_c_methyltransferase import TocopherolCMethyltransferase


class OMethyltransferase(ReactionClass):
    """O-methyltransferase activity."""

    GO_ID = "GO:0008171"
    CHILD_CLASSES = (
        TocopherolCMethyltransferase,
        TenHydroxydihydrosanguinarine10OMethyltransferase,
        GibberellinA9OMethyltransferase,
        GibberellinA4CarboxylMethyltransferase,
    )

    @classmethod
    def supports_evaluation(cls, reaction: Reaction) -> bool:
        return aggregate_explicit_child_supports_evaluation(reaction, cls.CHILD_CLASSES)

    def check_membership_impl(self, reaction: Reaction):
        return aggregate_explicit_child_membership(
            reaction,
            self.CHILD_CLASSES,
            "O-methyltransferase",
        )
