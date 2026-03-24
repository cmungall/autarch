"""translocation of hydrons."""

from autarch.ontology.ec_prefix_aggregate import ExplicitEcAggregate
from autarch.ontology.cytochrome_c_oxidase import CytochromeCOxidase


class TranslocationOfHydrons(ExplicitEcAggregate):
    """translocation of hydrons."""

    EC_NUMBER_PREFIX = '7.1.-.-'
    CHILD_CLASSES = (
        CytochromeCOxidase,
    )
