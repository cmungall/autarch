"""ETL (Extract, Transform, Load) utilities for importing chemical data.

This package contains modules for importing and processing chemical data
from various sources like ChEBI and Rhea databases.
"""

from autarch.etl.chebi_etl import (
    fetch_go_enzyme_mappings,
    mine_go_enzyme_classifications,
    serialize_go_terms,
)

from autarch.etl.rhea_etl import (
    fetch_rhea_reactions,
    serialize_rhea_reactions,
    parse_location_from_label,
)

from autarch.etl.chebi_smiles import (
    fetch_chebi_smiles_from_rhea,
)

from autarch.etl.modelseed_etl import (
    cache_modelseed_dataset,
    summarize_modelseed_cache,
)

__all__ = [
    # ChEBI/GO ETL
    "fetch_go_enzyme_mappings",
    "mine_go_enzyme_classifications",
    "serialize_go_terms",
    # Rhea ETL
    "fetch_rhea_reactions",
    "serialize_rhea_reactions",
    "parse_location_from_label",
    # ChEBI SMILES
    "fetch_chebi_smiles_from_rhea",
    # ModelSEED bridge
    "cache_modelseed_dataset",
    "summarize_modelseed_cache",
]
