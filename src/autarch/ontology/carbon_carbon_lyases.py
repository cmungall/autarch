"""carbon-carbon lyases."""

from autarch.ontology.ec_prefix_aggregate import ExplicitEcAggregate
from autarch.ontology.aldehyde_lyase import AldehydeLyase
from autarch.ontology.carboxy_lyase import CarboxyLyase
from autarch.ontology.methyl_ethyl_malonyl_coa_decarboxylase import MethylEthylMalonylCoADecarboxylase
from autarch.ontology.other_carbon_carbon_lyases import OtherCarbonCarbonLyases
from autarch.ontology.oxo_acid_lyase import OxoAcidLyase
from autarch.ontology.ribulose_bisphosphate_carboxylase import RibuloseBisphosphateCarboxylase


class CarbonCarbonLyases(ExplicitEcAggregate):
    """carbon-carbon lyases."""

    EC_NUMBER_PREFIX = '4.1.-.-'
    CHILD_CLASSES = (
        MethylEthylMalonylCoADecarboxylase,
        RibuloseBisphosphateCarboxylase,
        AldehydeLyase,
        CarboxyLyase,
        OtherCarbonCarbonLyases,
        OxoAcidLyase,
    )
