"""carbohydrate transmembrane transporter activity."""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.abc_type_carbohydrate_transporter import ABCTypeCarbohydrateTransporter


class CarbohydrateTransmembraneTransporterActivity(ExplicitGoAggregate):
    """carbohydrate transmembrane transporter activity."""

    GO_ID = "GO:0015144"
    CONCEPT_PHRASE = "carbohydrate transmembrane transporter activity"
    CHILD_CLASSES = (ABCTypeCarbohydrateTransporter,)
