"""CHEBI to SMILES conversion utilities.

Provides bidirectional mapping between ChEBI IDs and SMILES strings:
- ChEBI → SMILES: For looking up structures from identifiers
- SMILES → ChEBI: For reverse lookup when parsing reaction SMILES

The reverse lookup uses RDKit canonicalization to handle different
SMILES representations of the same molecule.
"""

from typing import Dict, Optional, Set, Tuple
import json
import sqlite3
from pathlib import Path

from rdkit import Chem


def canonicalize_smiles(smiles: str) -> Optional[str]:
    """Canonicalize a SMILES string using RDKit.

    Args:
        smiles: Input SMILES string

    Returns:
        Canonical SMILES or None if parsing fails

    Examples:
        >>> canonicalize_smiles("O")
        'O'
        >>> canonicalize_smiles("[H]O[H]")
        'O'
        >>> canonicalize_smiles("C(C)O")
        'CCO'
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True)


def fetch_chebi_smiles_from_rhea(cache_dir: str = "cache") -> Dict[str, str]:
    """Fetch SMILES for CHEBI IDs referenced in cached RHEA reactions.

    Args:
        cache_dir: Directory containing cached RHEA reactions

    Returns:
        Dictionary mapping CHEBI ID to SMILES string
    """
    cache_path = Path(cache_dir)
    rhea_cache = cache_path / "rhea_reactions.jsonl"

    if not rhea_cache.exists():
        raise FileNotFoundError(f"RHEA cache not found at {rhea_cache}")

    # Collect all unique CHEBI IDs from reactions
    all_chebi_ids: Set[str] = set()
    with open(rhea_cache) as f:
        for line in f:
            reaction_data = json.loads(line)
            # Handle new format with reaction.left_participants/right_participants
            if "reaction" in reaction_data and reaction_data["reaction"]:
                reaction = reaction_data["reaction"]
                if "left_participants" in reaction:
                    for p in reaction["left_participants"]:
                        if p and p.get("chebi_id"):
                            all_chebi_ids.add(p["chebi_id"])
                if "right_participants" in reaction:
                    for p in reaction["right_participants"]:
                        if p and p.get("chebi_id"):
                            all_chebi_ids.add(p["chebi_id"])
            # Also handle old format for backward compatibility
            all_chebi_ids.update(reaction_data.get("inputs", []))
            all_chebi_ids.update(reaction_data.get("outputs", []))

    print(f"Found {len(all_chebi_ids)} unique CHEBI IDs in RHEA reactions")

    # Connect to CHEBI database
    chebi_db = Path.home() / ".data" / "oaklib" / "chebi.db"
    if not chebi_db.exists():
        raise FileNotFoundError(f"CHEBI database not found at {chebi_db}")

    conn = sqlite3.connect(chebi_db)
    cursor = conn.cursor()

    chebi_to_smiles = {}
    failed_ids = []

    for i, chebi_id in enumerate(all_chebi_ids, 1):
        if i % 100 == 0:
            print(f"Progress: {i}/{len(all_chebi_ids)} CHEBI IDs processed...")

        try:
            # Query for SMILES in statements table
            # The statements table has subject, predicate, and value columns
            # SMILES is typically stored with a specific predicate
            cursor.execute(
                """
                SELECT value 
                FROM statements 
                WHERE subject = ? 
                AND predicate LIKE '%SMILES%'
                LIMIT 1
            """,
                (chebi_id,),
            )

            result = cursor.fetchone()
            if result and result[0]:
                chebi_to_smiles[chebi_id] = result[0]
            else:
                # Try alternate query - check in entity annotations
                cursor.execute(
                    """
                    SELECT value
                    FROM statements
                    WHERE subject = ?
                    AND predicate = 'chebi:smiles'
                    LIMIT 1
                """,
                    (chebi_id,),
                )

                result = cursor.fetchone()
                if result and result[0]:
                    chebi_to_smiles[chebi_id] = result[0]
                else:
                    failed_ids.append(chebi_id)
        except Exception:
            failed_ids.append(chebi_id)
            continue

    conn.close()

    print(f"Successfully fetched SMILES for {len(chebi_to_smiles)} CHEBI IDs")
    print(f"Failed to fetch SMILES for {len(failed_ids)} CHEBI IDs")

    return chebi_to_smiles


def load_chebi_smiles_cache(cache_dir: str = "cache") -> Dict[str, str]:
    """Load CHEBI to SMILES mapping from cache.

    Args:
        cache_dir: Directory containing cached data

    Returns:
        Dictionary mapping CHEBI ID to SMILES string

    Raises:
        FileNotFoundError: If cache file doesn't exist
    """
    cache_path = Path(cache_dir) / "chebi_smiles.json"
    if not cache_path.exists():
        raise FileNotFoundError(f"CHEBI SMILES cache not found at {cache_path}")

    with open(cache_path) as f:
        return json.load(f)


# =============================================================================
# SMILES → ChEBI reverse lookup
# =============================================================================


def build_smiles_to_chebi_lookup(
    chebi_db_path: Optional[Path] = None,
    progress_interval: int = 10000,
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Build SMILES→ChEBI reverse lookup from ChEBI database.

    Uses RDKit canonicalization to ensure consistent SMILES matching.

    Args:
        chebi_db_path: Path to ChEBI SQLite database (default: ~/.data/oaklib/chebi.db)
        progress_interval: Print progress every N entries

    Returns:
        Tuple of (canonical_smiles_to_chebi, chebi_to_name) dictionaries

    Examples:
        >>> lookup, names = build_smiles_to_chebi_lookup()  # doctest: +SKIP
        >>> lookup.get("O")  # doctest: +SKIP
        'CHEBI:15377'
    """
    if chebi_db_path is None:
        chebi_db_path = Path.home() / ".data" / "oaklib" / "chebi.db"

    if not chebi_db_path.exists():
        raise FileNotFoundError(f"ChEBI database not found at {chebi_db_path}")

    conn = sqlite3.connect(chebi_db_path)
    cursor = conn.cursor()

    # Get all SMILES from ChEBI using the correct predicate
    cursor.execute(
        """
        SELECT subject, value
        FROM statements
        WHERE predicate = 'chemrof:smiles_string'
        """
    )
    smiles_rows = cursor.fetchall()
    print(f"Found {len(smiles_rows)} SMILES entries in ChEBI database")

    # Get names for all ChEBI IDs
    cursor.execute(
        """
        SELECT subject, value
        FROM statements
        WHERE predicate = 'rdfs:label'
        AND subject LIKE 'CHEBI:%'
        """
    )
    name_rows = cursor.fetchall()
    chebi_to_name = {row[0]: row[1] for row in name_rows}

    conn.close()

    # Build canonical SMILES lookup
    smiles_to_chebi: Dict[str, str] = {}
    skipped = 0

    for i, (chebi_id, smiles) in enumerate(smiles_rows):
        if i > 0 and i % progress_interval == 0:
            print(f"Progress: {i}/{len(smiles_rows)} entries processed...")

        # Handle SMILES: prefix if present
        if smiles.startswith("SMILES:"):
            smiles = smiles[7:]

        # Canonicalize SMILES
        canonical = canonicalize_smiles(smiles)
        if canonical is None:
            skipped += 1
            continue

        # Store mapping (first ChEBI ID wins for duplicates)
        if canonical not in smiles_to_chebi:
            smiles_to_chebi[canonical] = chebi_id

    print(f"Built lookup with {len(smiles_to_chebi)} canonical SMILES")
    print(f"Skipped {skipped} invalid SMILES")
    print(f"Loaded {len(chebi_to_name)} ChEBI names")

    return smiles_to_chebi, chebi_to_name


def cache_smiles_to_chebi(cache_dir: str = "cache") -> None:
    """Build and cache SMILES→ChEBI lookup.

    Args:
        cache_dir: Directory to save cache files
    """
    cache_path = Path(cache_dir)
    cache_path.mkdir(exist_ok=True)

    smiles_lookup, name_lookup = build_smiles_to_chebi_lookup()

    # Save SMILES → ChEBI lookup
    with open(cache_path / "smiles_to_chebi.json", "w") as f:
        json.dump(smiles_lookup, f, indent=2)
    print(f"Saved SMILES→ChEBI lookup to {cache_path / 'smiles_to_chebi.json'}")

    # Save ChEBI → name lookup
    with open(cache_path / "chebi_names.json", "w") as f:
        json.dump(name_lookup, f, indent=2)
    print(f"Saved ChEBI names to {cache_path / 'chebi_names.json'}")


def load_smiles_to_chebi_cache(
    cache_dir: str = "cache",
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Load SMILES→ChEBI lookup from cache.

    Args:
        cache_dir: Directory containing cache files

    Returns:
        Tuple of (smiles_to_chebi, chebi_to_name) dictionaries

    Raises:
        FileNotFoundError: If cache files don't exist
    """
    cache_path = Path(cache_dir)

    smiles_cache = cache_path / "smiles_to_chebi.json"
    if not smiles_cache.exists():
        raise FileNotFoundError(
            f"SMILES→ChEBI cache not found at {smiles_cache}. "
            "Run 'autarch cache-chebi' to build it."
        )

    names_cache = cache_path / "chebi_names.json"

    with open(smiles_cache) as f:
        smiles_to_chebi = json.load(f)

    chebi_to_name: Dict[str, str] = {}
    if names_cache.exists():
        with open(names_cache) as f:
            chebi_to_name = json.load(f)

    return smiles_to_chebi, chebi_to_name


class ChEBILookup:
    """Singleton class for SMILES→ChEBI lookups with caching.

    Examples:
        >>> lookup = ChEBILookup.get_instance()  # doctest: +SKIP
        >>> lookup.get_chebi_id("O")  # doctest: +SKIP
        'CHEBI:15377'
        >>> lookup.get_name("CHEBI:15377")  # doctest: +SKIP
        'water'
    """

    _instance: Optional["ChEBILookup"] = None
    _smiles_to_chebi: Dict[str, str]
    _chebi_to_name: Dict[str, str]

    def __init__(self, cache_dir: str = "cache"):
        self._smiles_to_chebi, self._chebi_to_name = load_smiles_to_chebi_cache(
            cache_dir
        )

    @classmethod
    def get_instance(cls, cache_dir: str = "cache") -> "ChEBILookup":
        """Get singleton instance of ChEBI lookup."""
        if cls._instance is None:
            cls._instance = cls(cache_dir)
        return cls._instance

    def get_chebi_id(self, smiles: str) -> Optional[str]:
        """Look up ChEBI ID for a SMILES string.

        Args:
            smiles: SMILES string (will be canonicalized)

        Returns:
            ChEBI ID or None if not found
        """
        canonical = canonicalize_smiles(smiles)
        if canonical is None:
            return None
        return self._smiles_to_chebi.get(canonical)

    def get_name(self, chebi_id: str) -> Optional[str]:
        """Get name for a ChEBI ID.

        Args:
            chebi_id: ChEBI identifier (e.g., "CHEBI:15377")

        Returns:
            Chemical name or None if not found
        """
        return self._chebi_to_name.get(chebi_id)
