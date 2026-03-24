"""hydrolase acting on carbon-carbon bonds in ketonic substances.

Catalysis of the hydrolysis of carbon-carbon bonds in ketonic substances.
"""

from autarch.ontology.hydrolase_acting_on_acid_carbon_carbon_bonds_in_ketonic_substances import (
    HydrolaseActingOnAcidCarbonCarbonBondsInKetonicSubstances,
)


class HydrolaseActingOnCarbonCarbonBondsInKetonicSubstances(
    HydrolaseActingOnAcidCarbonCarbonBondsInKetonicSubstances
):
    """hydrolase acting on carbon-carbon bonds in ketonic substances.

    Catalysis of the hydrolysis of carbon-carbon bonds in ketonic substances.
    """

    GO_ID = "GO:0016823"
    EC_NUMBER_PREFIX = "3.7.1.-"

    def _implementation_only(self) -> None:
        """Concrete benchmark wrapper."""
