"""amino-acid dehydrogenase [NAD(P)+] activity."""

from typing import ClassVar, Optional

from autarch.ontology.oxidoreductase_acting_on_the_ch_nh2_group_of_donors_nad_or_nadp_as_acceptor import (
    OxidoreductaseActingOnTheCHNH2GroupOfDonorsNADOrNADPAsAcceptor,
)


class AminoAcidDehydrogenaseNADOrNADP(
    OxidoreductaseActingOnTheCHNH2GroupOfDonorsNADOrNADPAsAcceptor
):
    """amino-acid dehydrogenase [NAD(P)+] activity."""

    GO_ID = "GO:0050018"
    EC_NUMBER_PREFIX: ClassVar[Optional[str]] = None
