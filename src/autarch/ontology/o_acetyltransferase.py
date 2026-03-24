"""O-acetyltransferase activity."""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.polysialic_acid_o_acetyltransferase import (
    PolysialicAcidOAcetyltransferase,
)


class OAcetyltransferase(ExplicitGoAggregate):
    """O-acetyltransferase activity."""

    GO_ID = "GO:0016413"
    CONCEPT_PHRASE = "O-acetyltransferase activity"
    CHILD_CLASSES = (PolysialicAcidOAcetyltransferase,)
