"""cis-trans isomerases."""

from autarch.datamodel import Reaction
from autarch.ontology.ec_prefix_aggregate import (
    aggregate_explicit_child_membership,
    aggregate_explicit_child_supports_evaluation,
)
from autarch.ontology.peptidyl_prolyl_cis_trans_isomerase import PeptidylProlylCisTransIsomerase
from autarch.ontology.reaction import ReactionClass


class CisTransIsomerases(ReactionClass):
    """cis-trans isomerases."""

    EC_NUMBER_PREFIX = "5.2.1.-"
    CHILD_CLASSES = (
        PeptidylProlylCisTransIsomerase,
    )

    @classmethod
    def supports_evaluation(cls, reaction: Reaction) -> bool:
        return aggregate_explicit_child_supports_evaluation(reaction, cls.CHILD_CLASSES)

    def check_membership_impl(self, reaction: Reaction):
        return aggregate_explicit_child_membership(
            reaction,
            self.CHILD_CLASSES,
            "Cis-trans isomerases",
        )
