"""steroid dehydrogenase activity."""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.steroid_dehydrogenase_activity_acting_on_the_ch_oh_group_of_donors_nad_or_nadp_as_acceptor import (
    SteroidDehydrogenaseActivityActingOnTheCHOHGroupOfDonorsNADOrNADPAsAcceptor,
)


class SteroidDehydrogenaseActivity(ExplicitGoAggregate):
    """steroid dehydrogenase activity."""

    GO_ID = "GO:0016229"
    EC_BROAD_XREFS = ["1.-.-.-"]
    CONCEPT_PHRASE = "steroid dehydrogenase activity"
    CHILD_CLASSES = (SteroidDehydrogenaseActivityActingOnTheCHOHGroupOfDonorsNADOrNADPAsAcceptor,)
