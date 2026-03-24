"""mannosyltransferase activity."""

from autarch.molecules import CHEBI_GDP
from autarch.ontology.specific_hexosyltransferase_activity import (
    SpecificHexosyltransferaseActivity,
)


class MannosyltransferaseActivity(SpecificHexosyltransferaseActivity):
    """mannosyltransferase activity."""

    GO_ID = "GO:0000030"
    CONCEPT_PHRASE = "mannosyltransferase activity"
    DONOR_IDS = frozenset({"CHEBI:57527"})
    PRODUCT_IDS = frozenset({CHEBI_GDP})
