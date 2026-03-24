try:
    from autarch._version import __version__, __version_tuple__
except ImportError:  # pragma: no cover
    __version__ = "0.0.0"
    __version_tuple__ = (0, 0, 0)

# Import data models
from autarch.datamodel import (
    GoTerm,
    RheaTerm,
    Reaction,
    Participant,
    ClassificationResult,
)

# Import ETL functions
from autarch.etl.chebi_etl import (
    fetch_go_enzyme_mappings,
    mine_go_enzyme_classifications,
    serialize_go_terms,
)
from autarch.etl.rhea_etl import (
    fetch_rhea_reactions,
    serialize_rhea_reactions,
)

__all__ = [
    "__version__",
    "__version_tuple__",
    # Data models
    "GoTerm",
    "RheaTerm",
    "Reaction",
    "Participant",
    "ClassificationResult",
    # ETL functions
    "fetch_go_enzyme_mappings",
    "mine_go_enzyme_classifications",
    "serialize_go_terms",
    "fetch_rhea_reactions",
    "serialize_rhea_reactions",
]
