"""forming carbon-oxygen bonds."""

from autarch.ontology.ec_prefix_aggregate import ExplicitEcAggregate
from autarch.ontology.aminoacyl_trna_ligase import AminoacylTRNALigase


class FormingCarbonOxygenBonds(ExplicitEcAggregate):
    """forming carbon-oxygen bonds."""

    EC_NUMBER_PREFIX = '6.1.-.-'
    CHILD_CLASSES = (
        AminoacylTRNALigase,
    )
