"""RNA dihydrouridine synthase activity."""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.trna_dihydrouridine_1617_synthase_nad_p import TRNADihydrouridine1617SynthaseNADP
from autarch.ontology.trna_dihydrouridine_20a20b_synthase_nad_p import TRNADihydrouridine20A20BSynthaseNADP


class RNADihydrouridineSynthase(ExplicitGoAggregate):
    """RNA dihydrouridine synthase activity."""

    GO_ID = "GO:0106413"
    CONCEPT_PHRASE = "RNA dihydrouridine synthase activity"
    CHILD_CLASSES = (TRNADihydrouridine1617SynthaseNADP, TRNADihydrouridine20A20BSynthaseNADP,)
