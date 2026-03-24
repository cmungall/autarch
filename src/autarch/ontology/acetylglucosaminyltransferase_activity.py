"""acetylglucosaminyltransferase activity."""

from autarch.molecules import CHEBI_UDP
from autarch.ontology.specific_hexosyltransferase_activity import (
    SpecificHexosyltransferaseActivity,
)


class AcetylglucosaminyltransferaseActivity(SpecificHexosyltransferaseActivity):
    """acetylglucosaminyltransferase activity."""

    GO_ID = "GO:0008375"
    CONCEPT_PHRASE = "acetylglucosaminyltransferase activity"
    DONOR_IDS = frozenset({"CHEBI:57705"})
    PRODUCT_IDS = frozenset({CHEBI_UDP})
