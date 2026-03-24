"""carbohydrate kinase activity."""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.hexokinase import Hexokinase


class CarbohydrateKinaseActivity(ExplicitGoAggregate):
    """carbohydrate kinase activity."""

    GO_ID = "GO:0019200"
    CONCEPT_PHRASE = "carbohydrate kinase activity"
    CHILD_CLASSES = (Hexokinase,)
