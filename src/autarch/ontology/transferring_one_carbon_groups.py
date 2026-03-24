"""transferring one-carbon groups."""

from autarch.ontology.ec_prefix_aggregate import ExplicitEcAggregate
from autarch.ontology.caffeine_synthase import CaffeineSynthase
from autarch.ontology.carboxyl_or_carbamoyltransferase import CarboxylOrCarbamoyltransferase
from autarch.ontology.demethylmenaquinone_methyltransferase import DemethylmenaquinoneMethyltransferase
from autarch.ontology.five_methyltetrahydropteroyltriglutamate_homocysteine_s_methyltransferase import FiveMethyltetrahydropteroyltriglutamateHomocysteineSMethyltransferase
from autarch.ontology.gibberellin_a4_carboxyl_methyltransferase import GibberellinA4CarboxylMethyltransferase
from autarch.ontology.gibberellin_a9_o_methyltransferase import GibberellinA9OMethyltransferase
from autarch.ontology.hydroxymethyl_formyl_and_related_transferase import HydroxymethylFormylAndRelatedTransferase
from autarch.ontology.l_histidine_n_alpha_methyltransferase import LHistidineNAlphaMethyltransferase
from autarch.ontology.methionine_s_methyltransferase import MethionineSMethyltransferase
from autarch.ontology.methoxylated_aromatic_compound_corrinoid_protein_co_methyltransferase import MethoxylatedAromaticCompoundCorrinoidProteinCoMethyltransferase
from autarch.ontology.methyltransferase import Methyltransferase
from autarch.ontology.nicotinamide_n_methyltransferase import NicotinamideNMethyltransferase
from autarch.ontology.phosphoethanolamine_n_methyltransferase import PhosphoethanolamineNMethyltransferase
from autarch.ontology.precorrin_3b_c17_methyltransferase import Precorrin3BC17Methyltransferase
from autarch.ontology.rrna_small_subunit_pseudouridine_methyltransferase_nep1 import RRNASmallSubunitPseudouridineMethyltransferaseNep1
from autarch.ontology.rs1_benzyl_1234_tetrahydroisoquinoline_n_methyltransferase import RS1Benzyl1234TetrahydroisoquinolineNMethyltransferase
from autarch.ontology.ten_hydroxydihydrosanguinarine_10_o_methyltransferase import TenHydroxydihydrosanguinarine10OMethyltransferase
from autarch.ontology.tocopherol_c_methyltransferase import TocopherolCMethyltransferase
from autarch.ontology.trimethylsulfonium_tetrahydrofolate_n_methyltransferase import TrimethylsulfoniumTetrahydrofolateNMethyltransferase
from autarch.ontology.trna_cytidine_5_methyltransferase import TRNACytidine5Methyltransferase
from autarch.ontology.tyramine_n_methyltransferase import TyramineNMethyltransferase


class TransferringOneCarbonGroups(ExplicitEcAggregate):
    """transferring one-carbon groups."""

    EC_NUMBER_PREFIX = '2.1.-.-'
    CHILD_CLASSES = (
        CaffeineSynthase,
        DemethylmenaquinoneMethyltransferase,
        FiveMethyltetrahydropteroyltriglutamateHomocysteineSMethyltransferase,
        GibberellinA4CarboxylMethyltransferase,
        GibberellinA9OMethyltransferase,
        LHistidineNAlphaMethyltransferase,
        MethionineSMethyltransferase,
        MethoxylatedAromaticCompoundCorrinoidProteinCoMethyltransferase,
        NicotinamideNMethyltransferase,
        PhosphoethanolamineNMethyltransferase,
        Precorrin3BC17Methyltransferase,
        RRNASmallSubunitPseudouridineMethyltransferaseNep1,
        RS1Benzyl1234TetrahydroisoquinolineNMethyltransferase,
        TRNACytidine5Methyltransferase,
        TenHydroxydihydrosanguinarine10OMethyltransferase,
        TocopherolCMethyltransferase,
        TrimethylsulfoniumTetrahydrofolateNMethyltransferase,
        TyramineNMethyltransferase,
        CarboxylOrCarbamoyltransferase,
        HydroxymethylFormylAndRelatedTransferase,
        Methyltransferase,
    )
