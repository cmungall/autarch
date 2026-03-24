"""UDP-galactosyltransferase activity."""

from autarch.molecules import CHEBI_UDP
from autarch.ontology.specific_hexosyltransferase_activity import (
    SpecificHexosyltransferaseActivity,
)


class UDPGalactosyltransferaseActivity(SpecificHexosyltransferaseActivity):
    """UDP-galactosyltransferase activity."""

    GO_ID = "GO:0035250"
    CONCEPT_PHRASE = "UDP-galactosyltransferase activity"
    DONOR_IDS = frozenset({"CHEBI:18307", "CHEBI:66914"})
    PRODUCT_IDS = frozenset({CHEBI_UDP})
