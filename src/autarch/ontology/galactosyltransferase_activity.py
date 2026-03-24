"""galactosyltransferase activity."""

from autarch.molecules import CHEBI_UDP
from autarch.ontology.specific_hexosyltransferase_activity import (
    SpecificHexosyltransferaseActivity,
)


class GalactosyltransferaseActivity(SpecificHexosyltransferaseActivity):
    """galactosyltransferase activity."""

    GO_ID = "GO:0008378"
    CONCEPT_PHRASE = "galactosyltransferase activity"
    DONOR_IDS = frozenset({"CHEBI:18307", "CHEBI:66914"})
    PRODUCT_IDS = frozenset({CHEBI_UDP})
