"""nucleobase-containing compound kinase activity."""

from autarch.datamodel import Reaction
from autarch.ontology.deoxynucleoside_kinase import DeoxynucleosideKinase
from autarch.ontology.ec_prefix_aggregate import (
    aggregate_explicit_child_membership,
    aggregate_explicit_child_supports_evaluation,
)
from autarch.ontology.nucleoside_diphosphate_kinase import NucleosideDiphosphateKinase
from autarch.ontology.reaction import ReactionClass
from autarch.ontology.uridine_cytidine_kinase import UridineCytidineKinase


class NucleobaseContainingCompoundKinase(ReactionClass):
    """nucleobase-containing compound kinase activity."""

    GO_ID = "GO:0019205"
    CHILD_CLASSES = (
        DeoxynucleosideKinase,
        UridineCytidineKinase,
        NucleosideDiphosphateKinase,
    )

    @classmethod
    def supports_evaluation(cls, reaction: Reaction) -> bool:
        return aggregate_explicit_child_supports_evaluation(reaction, cls.CHILD_CLASSES)

    def check_membership_impl(self, reaction: Reaction):
        return aggregate_explicit_child_membership(
            reaction,
            self.CHILD_CLASSES,
            "Nucleobase-containing compound kinase",
        )
