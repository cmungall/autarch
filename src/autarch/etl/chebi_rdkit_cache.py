"""Cache RDKit-derived molecular data by ChEBI ID.

Caches: InChI, canonical SMILES (with/without stereo), chiral centers.
Populated at ETL time, reused during classification.

The cache avoids expensive RDKit computations during classification by
pre-computing and storing molecular properties indexed by ChEBI ID.
"""

import json
from pathlib import Path
from typing import Optional

from rdkit import Chem
from rdkit.Chem.inchi import MolToInchi

DEFAULT_CACHE_FILE = "cache/chebi_rdkit.json"

# Global cache instance (loaded lazily)
_cache: Optional[dict[str, dict]] = None
_cache_modified: bool = False


def compute_rdkit_data(smiles: str) -> Optional[dict]:
    """Compute RDKit-derived data from SMILES.

    Args:
        smiles: SMILES string for the molecule

    Returns:
        Dict with inchi, smiles_canonical, smiles_no_stereo, chiral_centers
        or None if SMILES cannot be parsed.

    Example:
        >>> data = compute_rdkit_data("C[C@H](N)C(=O)O")  # L-alanine
        >>> data["inchi"][:20]
        'InChI=1S/C3H7NO2/c1-'
        >>> len(data["chiral_centers"])
        1
    """
    if not smiles:
        return None

    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None

    # Assign stereochemistry from wedge/hash notation
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

    # Find chiral centers
    chiral_centers = []
    for idx, config in Chem.FindMolChiralCenters(mol, includeUnassigned=True):
        chiral_centers.append({"atom_idx": idx, "config": config})

    return {
        "inchi": MolToInchi(mol),
        "smiles_canonical": Chem.MolToSmiles(mol, isomericSmiles=True),
        "smiles_no_stereo": Chem.MolToSmiles(mol, isomericSmiles=False),
        "chiral_centers": chiral_centers,
    }


def load_cache(cache_file: str = DEFAULT_CACHE_FILE) -> dict[str, dict]:
    """Load cache from disk.

    Args:
        cache_file: Path to the cache JSON file

    Returns:
        Dict mapping ChEBI ID to RDKit data
    """
    global _cache

    if _cache is not None:
        return _cache

    cache_path = Path(cache_file)
    if cache_path.exists():
        with open(cache_path) as f:
            _cache = json.load(f)
    else:
        _cache = {}

    return _cache


def save_cache(cache_file: str = DEFAULT_CACHE_FILE) -> None:
    """Save cache to disk.

    Args:
        cache_file: Path to the cache JSON file
    """
    global _cache, _cache_modified

    if _cache is None or not _cache_modified:
        return

    cache_path = Path(cache_file)
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    with open(cache_path, "w") as f:
        json.dump(_cache, f, indent=2)

    _cache_modified = False


def get_rdkit_data(
    chebi_id: str,
    smiles: str,
    cache_file: str = DEFAULT_CACHE_FILE
) -> Optional[dict]:
    """Get RDKit data for a molecule, computing and caching if needed.

    Args:
        chebi_id: ChEBI identifier (e.g., "CHEBI:16977")
        smiles: SMILES string for the molecule
        cache_file: Path to the cache JSON file

    Returns:
        Dict with inchi, smiles_canonical, smiles_no_stereo, chiral_centers
        or None if data cannot be computed.

    Example:
        >>> data = get_rdkit_data("CHEBI:16977", "C[C@H](N)C(=O)O")
        >>> "inchi" in data
        True
    """
    global _cache, _cache_modified

    if not chebi_id or not smiles:
        return None

    cache = load_cache(cache_file)

    # Check if already cached
    if chebi_id in cache:
        return cache[chebi_id]

    # Compute and cache
    data = compute_rdkit_data(smiles)
    if data:
        cache[chebi_id] = data
        _cache_modified = True

    return data


def get_inchi(
    chebi_id: str,
    smiles: str,
    cache_file: str = DEFAULT_CACHE_FILE
) -> Optional[str]:
    """Get InChI for a molecule (convenience function).

    Args:
        chebi_id: ChEBI identifier
        smiles: SMILES string

    Returns:
        InChI string or None
    """
    data = get_rdkit_data(chebi_id, smiles, cache_file)
    return data["inchi"] if data else None


def clear_cache() -> None:
    """Clear the in-memory cache (for testing)."""
    global _cache, _cache_modified
    _cache = None
    _cache_modified = False


def populate_cache_from_smiles_map(
    chebi_smiles: dict[str, str],
    cache_file: str = DEFAULT_CACHE_FILE,
    progress: bool = True
) -> int:
    """Bulk populate cache from a ChEBI → SMILES mapping.

    Args:
        chebi_smiles: Dict mapping ChEBI ID to SMILES
        cache_file: Path to the cache JSON file
        progress: Whether to print progress

    Returns:
        Number of new entries added
    """
    cache = load_cache(cache_file)
    added = 0

    for chebi_id, smiles in chebi_smiles.items():
        if chebi_id not in cache:
            data = compute_rdkit_data(smiles)
            if data:
                cache[chebi_id] = data
                added += 1

    if added > 0:
        global _cache_modified
        _cache_modified = True
        save_cache(cache_file)

    if progress:
        print(f"Added {added} new entries to RDKit cache")

    return added
