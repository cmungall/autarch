"""transferring sulfur-containing groups."""

from autarch.ontology.ec_prefix_aggregate import ExplicitEcAggregate
from autarch.ontology.coa_transferase import CoATransferase
from autarch.ontology.sulfotransferase import Sulfotransferase
from autarch.ontology.sulfurtransferase import Sulfurtransferase


class TransferringSulfurContainingGroups(ExplicitEcAggregate):
    """transferring sulfur-containing groups."""

    EC_NUMBER_PREFIX = '2.8.-.-'
    CHILD_CLASSES = (
        CoATransferase,
        Sulfotransferase,
        Sulfurtransferase,
    )
