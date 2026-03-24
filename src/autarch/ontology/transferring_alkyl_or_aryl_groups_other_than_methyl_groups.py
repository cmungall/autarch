"""transferring alkyl or aryl groups other than methyl groups."""

from autarch.ontology.ec_prefix_aggregate import ExplicitEcAggregate
from autarch.ontology.corrinoid_adenosyltransferase import CorrinoidAdenosyltransferase
from autarch.ontology.prenyltransferase import Prenyltransferase


class TransferringAlkylOrArylGroupsOtherThanMethylGroups(ExplicitEcAggregate):
    """transferring alkyl or aryl groups other than methyl groups."""

    EC_NUMBER_PREFIX = '2.5.-.-'
    CHILD_CLASSES = (
        CorrinoidAdenosyltransferase,
        Prenyltransferase,
    )
