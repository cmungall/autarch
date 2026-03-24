"""catalytic activity, acting on a rRNA."""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.rrna_methyltransferase import RRNAMethyltransferase


class CatalyticActivityActingOnARRNA(ExplicitGoAggregate):
    """catalytic activity, acting on a rRNA."""

    GO_ID = "GO:0140102"
    CONCEPT_PHRASE = "catalytic activity, acting on a rRNA"
    CHILD_CLASSES = (RRNAMethyltransferase,)
