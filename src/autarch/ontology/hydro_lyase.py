"""hydro-lyase.

Catalysis of the cleavage of a carbon-oxygen bond by elimination of water.
"""

from autarch.datamodel import Reaction
from autarch.ontology.adp_dependent_nadh_or_nadph_hydrate_dehydratase import (
    ADPDependentNADHOrNADPHHydrateDehydratase,
)
from autarch.ontology.atp_dependent_nadh_or_nadph_hydrate_dehydratase import (
    ATPDependentNADHOrNADPHHydrateDehydratase,
)
from autarch.ontology.carbonate_dehydratase import CarbonateDehydratase
from autarch.ontology.ec_prefix_aggregate import (
    aggregate_explicit_child_membership,
    aggregate_explicit_child_supports_evaluation,
)
from autarch.ontology.glycerol_dehydratase import GlycerolDehydratase
from autarch.ontology.prephenate_dehydratase import PrephenateDehydratase
from autarch.ontology.reaction import ReactionClass
from autarch.ontology.three_hydroxyacyl_coa_dehydratase import (
    ThreeHydroxyacylCoADehydratase,
)


class HydroLyase(ReactionClass):
    """hydro-lyase.

    Catalysis of the cleavage of a carbon-oxygen bond by elimination of water.
    """

    GO_ID = "GO:0016836"
    EC_NUMBER_PREFIX = "4.2.1.-"
    CHILD_CLASSES = (
        ThreeHydroxyacylCoADehydratase,
        GlycerolDehydratase,
        PrephenateDehydratase,
        ATPDependentNADHOrNADPHHydrateDehydratase,
        ADPDependentNADHOrNADPHHydrateDehydratase,
        CarbonateDehydratase,
    )

    @classmethod
    def supports_evaluation(cls, reaction: Reaction) -> bool:
        return aggregate_explicit_child_supports_evaluation(reaction, cls.CHILD_CLASSES)

    def check_membership_impl(self, reaction: Reaction):
        return aggregate_explicit_child_membership(
            reaction,
            self.CHILD_CLASSES,
            "Hydro-lyase",
        )
