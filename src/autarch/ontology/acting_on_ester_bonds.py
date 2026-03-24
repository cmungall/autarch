"""acting on ester bonds."""

from autarch.ontology.ec_prefix_aggregate import ExplicitEcAggregate
from autarch.ontology.acyl_coa_hydrolase import AcylCoAHydrolase
from autarch.ontology.carboxylic_ester_hydrolase import CarboxylicEsterHydrolase
from autarch.ontology.d_aminoacyl_trna_deacylase import DAminoacylTRNADeacylase
from autarch.ontology.diphosphoric_monoester_hydrolase import DiphosphoricMonoesterHydrolase
from autarch.ontology.monoterpene_epsilon_lactone_hydrolase import MonoterpeneEpsilonLactoneHydrolase
from autarch.ontology.phosphatase import Phosphatase
from autarch.ontology.phosphoric_diester_hydrolase import PhosphoricDiesterHydrolase
from autarch.ontology.rna_nuclease import RNANuclease
from autarch.ontology.sulfuric_ester_hydrolase import SulfuricEsterHydrolase
from autarch.ontology.thiolester_hydrolase import ThiolesterHydrolase
from autarch.ontology.triacylglycerol_lipase import TriacylglycerolLipase
from autarch.ontology.two_three_cyclic_nucleotide_two_phosphodiesterase import TwoThreeCyclicNucleotideTwoPhosphodiesterase


class ActingOnEsterBonds(ExplicitEcAggregate):
    """acting on ester bonds."""

    EC_NUMBER_PREFIX = '3.1.-.-'
    CHILD_CLASSES = (
        DAminoacylTRNADeacylase,
        MonoterpeneEpsilonLactoneHydrolase,
        TriacylglycerolLipase,
        TwoThreeCyclicNucleotideTwoPhosphodiesterase,
        AcylCoAHydrolase,
        CarboxylicEsterHydrolase,
        DiphosphoricMonoesterHydrolase,
        Phosphatase,
        PhosphoricDiesterHydrolase,
        RNANuclease,
        SulfuricEsterHydrolase,
        ThiolesterHydrolase,
    )
