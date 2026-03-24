"""antioxidant activity."""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.peroxidase import Peroxidase
from autarch.ontology.superoxide_dismutase import SuperoxideDismutase


class AntioxidantActivity(ExplicitGoAggregate):
    """antioxidant activity."""

    GO_ID = "GO:0016209"
    CONCEPT_PHRASE = "antioxidant activity"
    CHILD_CLASSES = (Peroxidase, SuperoxideDismutase,)
