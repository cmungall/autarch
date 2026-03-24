"""catalytic activity, acting on RNA."""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.aminoacyl_trna_ligase import AminoacylTRNALigase
from autarch.ontology.d_aminoacyl_trna_deacylase import DAminoacylTRNADeacylase
from autarch.ontology.helicase import Helicase
from autarch.ontology.rna_methyltransferase import RNAMethyltransferase
from autarch.ontology.rna_nuclease import RNANuclease
from autarch.ontology.rna_polymerase import RNAPolymerase
from autarch.ontology.rrna_small_subunit_pseudouridine_methyltransferase_nep1 import RRNASmallSubunitPseudouridineMethyltransferaseNep1
from autarch.ontology.trna_cytidine_5_methyltransferase import TRNACytidine5Methyltransferase
from autarch.ontology.trna_dihydrouridine_1617_synthase_nad_p import TRNADihydrouridine1617SynthaseNADP
from autarch.ontology.trna_dihydrouridine_20a20b_synthase_nad_p import TRNADihydrouridine20A20BSynthaseNADP


class CatalyticActivityActingOnRNA(ExplicitGoAggregate):
    """catalytic activity, acting on RNA."""

    GO_ID = "GO:0140098"
    CONCEPT_PHRASE = "catalytic activity, acting on RNA"
    CHILD_CLASSES = (AminoacylTRNALigase, DAminoacylTRNADeacylase, Helicase, RNAMethyltransferase, RNANuclease, RNAPolymerase, TRNACytidine5Methyltransferase, RRNASmallSubunitPseudouridineMethyltransferaseNep1, TRNADihydrouridine1617SynthaseNADP, TRNADihydrouridine20A20BSynthaseNADP,)
