"""carbon-nitrogen lyases."""

from autarch.ontology.ec_prefix_aggregate import ExplicitEcAggregate
from autarch.ontology.ammonia_lyase import AmmoniaLyase
from autarch.ontology.lyases_acting_on_amides_amidines_etc import LyasesActingOnAmidesAmidinesEtc


class CarbonNitrogenLyases(ExplicitEcAggregate):
    """carbon-nitrogen lyases."""

    EC_NUMBER_PREFIX = '4.3.-.-'
    CHILD_CLASSES = (
        AmmoniaLyase,
        LyasesActingOnAmidesAmidinesEtc,
    )
