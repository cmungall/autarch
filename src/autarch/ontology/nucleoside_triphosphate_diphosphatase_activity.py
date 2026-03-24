"""nucleoside triphosphate diphosphatase activity."""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.atp_diphosphatase import ATPDiphosphatase
from autarch.ontology.nucleotide_diphosphatase import NucleotideDiphosphatase


class NucleosideTriphosphateDiphosphataseActivity(ExplicitGoAggregate):
    """nucleoside triphosphate diphosphatase activity."""

    GO_ID = "GO:0047429"
    CONCEPT_PHRASE = "nucleoside triphosphate diphosphatase activity"
    CHILD_CLASSES = (ATPDiphosphatase, NucleotideDiphosphatase,)
