"""lipase activity."""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.triacylglycerol_lipase import TriacylglycerolLipase


class Lipase(ExplicitGoAggregate):
    """lipase activity."""

    GO_ID = "GO:0016298"
    CONCEPT_PHRASE = "lipase activity"
    CHILD_CLASSES = (TriacylglycerolLipase,)
