"""aldehyde dehydrogenase [NAD(P)+] activity."""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.acting_on_the_aldehyde_or_oxo_group_of_donors_with_nad_or_nadp_as_acceptor import ActingOnTheAldehydeOrOxoGroupOfDonorsWithNADOrNADPAsAcceptor


class AldehydeDehydrogenaseNADOrNADP(ExplicitGoAggregate):
    """aldehyde dehydrogenase [NAD(P)+] activity."""

    GO_ID = "GO:0004030"
    CONCEPT_PHRASE = "aldehyde dehydrogenase [NAD(P)+] activity"
    CHILD_CLASSES = (ActingOnTheAldehydeOrOxoGroupOfDonorsWithNADOrNADPAsAcceptor,)
