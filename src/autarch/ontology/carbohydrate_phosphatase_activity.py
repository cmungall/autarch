"""carbohydrate phosphatase activity."""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.inositol_phosphate_phosphatase_activity import (
    InositolPhosphatePhosphataseActivity,
)
from autarch.ontology.phosphatidylinositol_phosphate_phosphatase_activity import (
    PhosphatidylinositolPhosphatePhosphataseActivity,
)
from autarch.ontology.sugar_phosphatase_activity import SugarPhosphataseActivity


class CarbohydratePhosphataseActivity(ExplicitGoAggregate):
    """carbohydrate phosphatase activity."""

    GO_ID = "GO:0019203"
    CONCEPT_PHRASE = "carbohydrate phosphatase activity"
    CHILD_CLASSES = (
        SugarPhosphataseActivity,
        InositolPhosphatePhosphataseActivity,
        PhosphatidylinositolPhosphatePhosphataseActivity,
    )
