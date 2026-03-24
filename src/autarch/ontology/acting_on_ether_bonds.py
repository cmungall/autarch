"""acting on ether bonds."""

from autarch.ontology.ec_prefix_aggregate import ExplicitEcAggregate
from autarch.ontology.epoxide_hydrolase import EpoxideHydrolase
from autarch.ontology.ether_hydrolase import EtherHydrolase


class ActingOnEtherBonds(ExplicitEcAggregate):
    """acting on ether bonds."""

    EC_NUMBER_PREFIX = '3.3.-.-'
    CHILD_CLASSES = (
        EpoxideHydrolase,
        EtherHydrolase,
    )
