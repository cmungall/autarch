"""N-acyltransferase activity."""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.acetyltransferase import Acetyltransferase


class NAcyltransferase(ExplicitGoAggregate):
    """N-acyltransferase activity."""

    GO_ID = "GO:0016410"
    CONCEPT_PHRASE = "N-acyltransferase activity"
    CHILD_CLASSES = (Acetyltransferase,)
