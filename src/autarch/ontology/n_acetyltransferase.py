"""N-acetyltransferase activity."""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.acetyltransferase import Acetyltransferase


class NAcetyltransferase(ExplicitGoAggregate):
    """N-acetyltransferase activity."""

    GO_ID = "GO:0008080"
    CONCEPT_PHRASE = "N-acetyltransferase activity"
    CHILD_CLASSES = (Acetyltransferase,)
