"""oxidoreductase acting on the aldehyde or oxo group of donors.

Catalysis of an oxidation-reduction (redox) reaction in which an aldehyde or
ketone (oxo) group acts as a hydrogen or electron donor and reduces a hydrogen
or electron acceptor.
"""

from autarch.datamodel import Reaction
from autarch.ontology.ec_prefix_aggregate import (
    aggregate_explicit_child_membership,
    aggregate_explicit_child_supports_evaluation,
)
from autarch.ontology.oxidoreductase_acting_on_the_aldehyde_or_oxo_group_of_donors_iron_sulfur_protein_as_acceptor import (
    OxidoreductaseActingOnTheAldehydeOrOxoGroupOfDonorsIronSulfurProteinAsAcceptor,
)
from autarch.ontology.oxidoreductase_acting_on_the_aldehyde_or_oxo_group_of_donors_nad_or_nadp_as_acceptor import (
    OxidoreductaseActingOnTheAldehydeOrOxoGroupOfDonorsNADOrNADPAsAcceptor,
)
from autarch.ontology.oxidoreductase_acting_on_the_aldehyde_or_oxo_group_of_donors_oxygen_as_acceptor import (
    OxidoreductaseActingOnTheAldehydeOrOxoGroupOfDonorsOxygenAsAcceptor,
)
from autarch.ontology.reaction import ReactionClass


class OxidoreductaseActingOnTheAldehydeOrOxoGroupOfDonors(ReactionClass):
    """oxidoreductase acting on the aldehyde or oxo group of donors.

    Catalysis of an oxidation-reduction (redox) reaction in which an aldehyde
    or ketone (oxo) group acts as a hydrogen or electron donor and reduces a
    hydrogen or electron acceptor.
    """

    GO_ID = "GO:0016903"
    EC_NUMBER_PREFIX = "1.2.-.-"
    CHILD_CLASSES = (
        OxidoreductaseActingOnTheAldehydeOrOxoGroupOfDonorsNADOrNADPAsAcceptor,
        OxidoreductaseActingOnTheAldehydeOrOxoGroupOfDonorsOxygenAsAcceptor,
        OxidoreductaseActingOnTheAldehydeOrOxoGroupOfDonorsIronSulfurProteinAsAcceptor,
    )

    @classmethod
    def supports_evaluation(cls, reaction: Reaction) -> bool:
        return aggregate_explicit_child_supports_evaluation(reaction, cls.CHILD_CLASSES)

    def check_membership_impl(self, reaction: Reaction):
        return aggregate_explicit_child_membership(
            reaction,
            self.CHILD_CLASSES,
            "Oxidoreductase acting on the aldehyde or oxo group of donors",
        )
