"""O-acyltransferase activity."""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.polysialic_acid_o_acetyltransferase import (
    PolysialicAcidOAcetyltransferase,
)


class OAcyltransferase(ExplicitGoAggregate):
    """O-acyltransferase activity."""

    GO_ID = "GO:0008374"
    CONCEPT_PHRASE = "O-acyltransferase activity"
    CHILD_CLASSES = (PolysialicAcidOAcetyltransferase,)
