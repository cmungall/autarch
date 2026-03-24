"""monoatomic cation transmembrane transporter activity."""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.cytochrome_c_oxidase import CytochromeCOxidase
from autarch.ontology.translocation_of_hydrons import TranslocationOfHydrons
from autarch.ontology.translocation_of_inorganic_cations_linked_to_the_hydrolysis_of_a_nucleoside_triphosphate import TranslocationOfInorganicCationsLinkedToTheHydrolysisOfANucleosideTriphosphate


class MonoatomicCationTransmembraneTransporterActivity(ExplicitGoAggregate):
    """monoatomic cation transmembrane transporter activity."""

    GO_ID = "GO:0008324"
    CONCEPT_PHRASE = "monoatomic cation transmembrane transporter activity"
    CHILD_CLASSES = (CytochromeCOxidase, TranslocationOfHydrons, TranslocationOfInorganicCationsLinkedToTheHydrolysisOfANucleosideTriphosphate,)
