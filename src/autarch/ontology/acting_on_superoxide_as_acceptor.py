"""acting on superoxide as acceptor."""

from autarch.ontology.ec_prefix_aggregate import ExplicitEcAggregate
from autarch.ontology.superoxide_dismutase import SuperoxideDismutase


class ActingOnSuperoxideAsAcceptor(ExplicitEcAggregate):
    """acting on superoxide as acceptor."""

    EC_NUMBER_PREFIX = '1.15.-.-'
    CHILD_CLASSES = (
        SuperoxideDismutase,
    )
