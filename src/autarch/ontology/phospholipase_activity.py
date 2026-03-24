"""phospholipase activity."""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.glycerophospholipase_activity import GlycerophospholipaseActivity


class PhospholipaseActivity(ExplicitGoAggregate):
    """phospholipase activity."""

    GO_ID = "GO:0120569"
    CONCEPT_PHRASE = "phospholipase activity"
    CHILD_CLASSES = (GlycerophospholipaseActivity,)
