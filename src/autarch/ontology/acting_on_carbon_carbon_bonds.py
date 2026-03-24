"""acting on carbon-carbon bonds."""

from autarch.ontology.ec_prefix_aggregate import ExplicitEcAggregate
from autarch.ontology.hydrolase_acting_on_carbon_carbon_bonds_in_ketonic_substances import HydrolaseActingOnCarbonCarbonBondsInKetonicSubstances


class ActingOnCarbonCarbonBonds(ExplicitEcAggregate):
    """acting on carbon-carbon bonds."""

    EC_NUMBER_PREFIX = '3.7.-.-'
    CHILD_CLASSES = (
        HydrolaseActingOnCarbonCarbonBondsInKetonicSubstances,
    )
