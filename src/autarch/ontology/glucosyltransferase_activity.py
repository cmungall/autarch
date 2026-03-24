"""glucosyltransferase activity."""

from autarch.molecules import CHEBI_ADP, CHEBI_GDP, CHEBI_UDP
from autarch.ontology.specific_hexosyltransferase_activity import (
    SpecificHexosyltransferaseActivity,
)


class GlucosyltransferaseActivity(SpecificHexosyltransferaseActivity):
    """glucosyltransferase activity."""

    GO_ID = "GO:0046527"
    CONCEPT_PHRASE = "glucosyltransferase activity"
    DONOR_IDS = frozenset(
        {
            "CHEBI:18066",
            "CHEBI:58885",
            "CHEBI:57498",
            "CHEBI:76533",
            "CHEBI:62230",
        }
    )
    PRODUCT_IDS = frozenset({CHEBI_UDP, CHEBI_ADP, CHEBI_GDP})
