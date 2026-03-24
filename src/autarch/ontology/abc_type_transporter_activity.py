"""ABC-type transporter activity."""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.abc_type_carbohydrate_transporter import ABCTypeCarbohydrateTransporter
from autarch.ontology.abc_type_polar_amino_acid_transporter import ABCTypePolarAminoAcidTransporter


class ABCTypeTransporterActivity(ExplicitGoAggregate):
    """ABC-type transporter activity."""

    GO_ID = "GO:0140359"
    CONCEPT_PHRASE = "ABC-type transporter activity"
    CHILD_CLASSES = (ABCTypeCarbohydrateTransporter, ABCTypePolarAminoAcidTransporter,)
