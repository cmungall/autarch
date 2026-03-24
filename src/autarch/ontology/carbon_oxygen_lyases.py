"""carbon-oxygen lyases."""

from autarch.ontology.ec_prefix_aggregate import ExplicitEcAggregate
from autarch.ontology.adp_dependent_nadh_or_nadph_hydrate_dehydratase import ADPDependentNADHOrNADPHHydrateDehydratase
from autarch.ontology.atp_dependent_nadh_or_nadph_hydrate_dehydratase import ATPDependentNADHOrNADPHHydrateDehydratase
from autarch.ontology.carbon_oxygen_lyase_acting_on_phosphates import CarbonOxygenLyaseActingOnPhosphates
from autarch.ontology.carbonate_dehydratase import CarbonateDehydratase
from autarch.ontology.gamma_humulene_synthase import GammaHumuleneSynthase
from autarch.ontology.glycerol_dehydratase import GlycerolDehydratase
from autarch.ontology.hydro_lyase import HydroLyase
from autarch.ontology.pectate_lyase import PectateLyase
from autarch.ontology.pinene_synthase import PineneSynthase
from autarch.ontology.prephenate_dehydratase import PrephenateDehydratase
from autarch.ontology.terpene_synthase import TerpeneSynthase
from autarch.ontology.three_hydroxyacyl_coa_dehydratase import ThreeHydroxyacylCoADehydratase


class CarbonOxygenLyases(ExplicitEcAggregate):
    """carbon-oxygen lyases."""

    EC_NUMBER_PREFIX = '4.2.-.-'
    CHILD_CLASSES = (
        ADPDependentNADHOrNADPHHydrateDehydratase,
        ATPDependentNADHOrNADPHHydrateDehydratase,
        CarbonateDehydratase,
        GammaHumuleneSynthase,
        GlycerolDehydratase,
        PectateLyase,
        PrephenateDehydratase,
        CarbonOxygenLyaseActingOnPhosphates,
        HydroLyase,
        PineneSynthase,
        TerpeneSynthase,
        ThreeHydroxyacylCoADehydratase,
    )
