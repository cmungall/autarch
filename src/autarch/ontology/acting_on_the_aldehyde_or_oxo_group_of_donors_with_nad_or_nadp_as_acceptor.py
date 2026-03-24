"""acting on the aldehyde or oxo group of donors with nad(+) or nadp(+) as acceptor."""

from autarch.datamodel import Reaction
from autarch.ontology.ec_prefix_aggregate import (
    aggregate_explicit_child_membership,
    aggregate_explicit_child_supports_evaluation,
)
from autarch.ontology.oxidoreductase_acting_on_the_aldehyde_or_oxo_group_of_donors_nad_or_nadp_as_acceptor import OxidoreductaseActingOnTheAldehydeOrOxoGroupOfDonorsNADOrNADPAsAcceptor
from autarch.ontology.reaction import ReactionClass


class ActingOnTheAldehydeOrOxoGroupOfDonorsWithNADOrNADPAsAcceptor(ReactionClass):
    """acting on the aldehyde or oxo group of donors with nad(+) or nadp(+) as acceptor."""

    EC_NUMBER_PREFIX = "1.2.1.-"
    CHILD_CLASSES = (
        OxidoreductaseActingOnTheAldehydeOrOxoGroupOfDonorsNADOrNADPAsAcceptor,
    )

    @classmethod
    def supports_evaluation(cls, reaction: Reaction) -> bool:
        return aggregate_explicit_child_supports_evaluation(reaction, cls.CHILD_CLASSES)

    def check_membership_impl(self, reaction: Reaction):
        return aggregate_explicit_child_membership(
            reaction,
            self.CHILD_CLASSES,
            "Acting on the aldehyde or oxo group of donors with NAD(+) or NADP(+) as acceptor",
        )
