"""ligase activity, forming carbon-oxygen bonds."""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.aminoacyl_trna_ligase import AminoacylTRNALigase


class LigaseActivityFormingCarbonOxygenBonds(ExplicitGoAggregate):
    """ligase activity, forming carbon-oxygen bonds."""

    GO_ID = "GO:0016875"
    CONCEPT_PHRASE = "ligase activity, forming carbon-oxygen bonds"
    CHILD_CLASSES = (AminoacylTRNALigase,)
