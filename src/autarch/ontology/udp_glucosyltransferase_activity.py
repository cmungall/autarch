"""UDP-glucosyltransferase activity."""

from autarch.molecules import CHEBI_UDP
from autarch.ontology.specific_hexosyltransferase_activity import (
    SpecificHexosyltransferaseActivity,
)


class UDPGlucosyltransferaseActivity(SpecificHexosyltransferaseActivity):
    """UDP-glucosyltransferase activity."""

    GO_ID = "GO:0035251"
    CONCEPT_PHRASE = "UDP-glucosyltransferase activity"
    DONOR_IDS = frozenset({"CHEBI:18066", "CHEBI:58885"})
    PRODUCT_IDS = frozenset({CHEBI_UDP})
