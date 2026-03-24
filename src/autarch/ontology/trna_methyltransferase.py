"""tRNA methyltransferase activity."""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.trna_cytidine_5_methyltransferase import TRNACytidine5Methyltransferase


class TRNAMethyltransferase(ExplicitGoAggregate):
    """tRNA methyltransferase activity."""

    GO_ID = "GO:0008175"
    CONCEPT_PHRASE = "tRNA methyltransferase activity"
    CHILD_CLASSES = (TRNACytidine5Methyltransferase,)
