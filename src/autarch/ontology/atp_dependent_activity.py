"""ATP-dependent activity."""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.abc_type_carbohydrate_transporter import ABCTypeCarbohydrateTransporter
from autarch.ontology.abc_type_polar_amino_acid_transporter import ABCTypePolarAminoAcidTransporter
from autarch.ontology.atpase_coupled_transmembrane_transporter import ATPaseCoupledTransmembraneTransporter
from autarch.ontology.helicase import Helicase
from autarch.ontology.translocase_linked_to_hydrolysis import TranslocaseLinkedToHydrolysis


class ATPDependentActivity(ExplicitGoAggregate):
    """ATP-dependent activity."""

    GO_ID = "GO:0140657"
    CONCEPT_PHRASE = "ATP-dependent activity"
    CHILD_CLASSES = (ATPaseCoupledTransmembraneTransporter, ABCTypeCarbohydrateTransporter, ABCTypePolarAminoAcidTransporter, Helicase, TranslocaseLinkedToHydrolysis,)
