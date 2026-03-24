"""rRNA methyltransferase activity."""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.rrna_small_subunit_pseudouridine_methyltransferase_nep1 import RRNASmallSubunitPseudouridineMethyltransferaseNep1


class RRNAMethyltransferase(ExplicitGoAggregate):
    """rRNA methyltransferase activity."""

    GO_ID = "GO:0008649"
    CONCEPT_PHRASE = "rRNA methyltransferase activity"
    CHILD_CLASSES = (RRNASmallSubunitPseudouridineMethyltransferaseNep1,)
