"""RNA methyltransferase activity."""

from autarch.datamodel import Reaction
from autarch.ontology.ec_prefix_aggregate import (
    aggregate_explicit_child_membership,
    aggregate_explicit_child_supports_evaluation,
)
from autarch.ontology.reaction import ReactionClass
from autarch.ontology.rrna_small_subunit_pseudouridine_methyltransferase_nep1 import (
    RRNASmallSubunitPseudouridineMethyltransferaseNep1,
)
from autarch.ontology.trna_cytidine_5_methyltransferase import TRNACytidine5Methyltransferase


class RNAMethyltransferase(ReactionClass):
    """RNA methyltransferase activity."""

    GO_ID = "GO:0008173"
    CHILD_CLASSES = (
        TRNACytidine5Methyltransferase,
        RRNASmallSubunitPseudouridineMethyltransferaseNep1,
    )

    @classmethod
    def supports_evaluation(cls, reaction: Reaction) -> bool:
        return aggregate_explicit_child_supports_evaluation(reaction, cls.CHILD_CLASSES)

    def check_membership_impl(self, reaction: Reaction):
        return aggregate_explicit_child_membership(
            reaction,
            self.CHILD_CLASSES,
            "RNA methyltransferase",
        )
