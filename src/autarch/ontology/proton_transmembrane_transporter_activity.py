"""proton transmembrane transporter activity."""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.cytochrome_c_oxidase import CytochromeCOxidase
from autarch.ontology.translocation_of_hydrons import TranslocationOfHydrons


class ProtonTransmembraneTransporterActivity(ExplicitGoAggregate):
    """proton transmembrane transporter activity."""

    GO_ID = "GO:0015078"
    CONCEPT_PHRASE = "proton transmembrane transporter activity"
    CHILD_CLASSES = (CytochromeCOxidase, TranslocationOfHydrons,)
