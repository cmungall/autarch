"""catalytic activity, acting on a tRNA."""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.aminoacyl_trna_ligase import AminoacylTRNALigase
from autarch.ontology.d_aminoacyl_trna_deacylase import DAminoacylTRNADeacylase
from autarch.ontology.trna_cytidine_5_methyltransferase import TRNACytidine5Methyltransferase
from autarch.ontology.trna_dihydrouridine_1617_synthase_nad_p import TRNADihydrouridine1617SynthaseNADP
from autarch.ontology.trna_dihydrouridine_20a20b_synthase_nad_p import TRNADihydrouridine20A20BSynthaseNADP


class CatalyticActivityActingOnATRNA(ExplicitGoAggregate):
    """catalytic activity, acting on a tRNA."""

    GO_ID = "GO:0140101"
    CONCEPT_PHRASE = "catalytic activity, acting on a tRNA"
    CHILD_CLASSES = (AminoacylTRNALigase, DAminoacylTRNADeacylase, TRNACytidine5Methyltransferase, TRNADihydrouridine1617SynthaseNADP, TRNADihydrouridine20A20BSynthaseNADP,)
