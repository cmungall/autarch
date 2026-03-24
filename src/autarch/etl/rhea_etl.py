"""RHEA ETL - Primary interface for fetching RHEA reaction data.

This module provides the main interface for fetching RHEA reaction data.
It uses TSV files from RHEA's FTP site as the primary data source.
"""

import logging
import json
import csv
import re
from pathlib import Path
from typing import Dict, List, Optional

from autarch.datamodel import RheaTerm
from autarch.etl.rhea_tsv_etl import RheaTSVETL

logger = logging.getLogger(__name__)


def parse_location_from_label(label: str) -> dict[str, str]:
    """Parse location/compartment information from RHEA label.

    Looks for patterns like:
    - ATP(in) + ADP(out) = ATP(out) + ADP(in)
    - H(+)(periplasm) = H(+)(cytoplasm)
    - glucose(extracellular) = glucose(cytoplasm)
    - an [alpha-factor](in) = an [alpha-factor](out)

    Args:
        label: RHEA reaction label

    Returns:
        Dictionary mapping molecule names to their locations

    >>> parse_location_from_label("ATP(in) + H2O = ADP(out) + Pi")
    {'ATP': 'in', 'ADP': 'out'}
    >>> parse_location_from_label("glucose(extracellular) = glucose(cytoplasm)")
    {'glucose': 'cytoplasm'}
    >>> parse_location_from_label("an [alpha-factor](in) + ATP = an [alpha-factor](out) + ADP")
    {'[alpha-factor]': 'out', 'ATP': None, 'ADP': None}
    """
    locations = {}

    # Enhanced pattern to match:
    # - Simple molecules: ATP(in), H(+)(periplasm)
    # - Bracketed molecules: [alpha-factor](in), an [alpha-factor](out)
    # Also handle molecules without locations

    # First pattern for bracketed molecules with optional "a/an" prefix
    bracket_pattern = r"(?:a|an)?\s*(\[[^\]]+\])\s*\(([a-zA-Z]+)\)"
    matches = re.findall(bracket_pattern, label)
    for molecule, location in matches:
        # For bracketed molecules, last occurrence wins
        locations[molecule] = location

    # Second pattern for regular molecules with locations
    regular_pattern = r"([A-Za-z0-9\(\)\+\-]+)\s*\(([a-zA-Z]+)\)"
    matches = re.findall(regular_pattern, label)
    for molecule, location in matches:
        # Skip if this looks like it's part of a bracketed molecule
        if "[" not in label[:label.find(molecule) if molecule in label else 0]:
            molecule = molecule.strip()
            # Always update - last occurrence wins
            locations[molecule] = location

    return locations


def extract_molecules_from_label(label: str) -> dict[str, dict]:
    """Extract molecule names and metadata from RHEA label.

    This parses the equation text to extract molecule names, especially
    for generic molecules without ChEBI IDs like "[alpha-factor]".

    Args:
        label: RHEA reaction equation/label

    Returns:
        Dictionary mapping molecule identifiers to their metadata
        (name, location, etc.)

    >>> extract_molecules_from_label("an [alpha-factor](in) + ATP = an [alpha-factor](out) + ADP")
    {'[alpha-factor]_in': {'name': '[alpha-factor]', 'location': 'in'}, '[alpha-factor]_out': {'name': '[alpha-factor]', 'location': 'out'}, 'ATP': {'name': 'ATP', 'location': None}, 'ADP': {'name': 'ADP', 'location': None}}
    """
    molecules = {}

    # Split equation into left and right sides
    if "=" in label:
        left_side, right_side = label.split("=", 1)
        sides = [("left", left_side), ("right", right_side)]
    else:
        sides = [("both", label)]

    for side_name, side_text in sides:
        # Split by + to get individual molecules
        parts = side_text.split("+")

        for part in parts:
            part = part.strip()
            if not part:
                continue

            # Handle bracketed molecules like "an [alpha-factor](in)"
            bracket_match = re.match(r"(?:a|an)?\s*(\[[^\]]+\])(?:\(([a-zA-Z]+)\))?", part)
            if bracket_match:
                name = bracket_match.group(1)
                location = bracket_match.group(2)
                # Create unique key with location if present
                key = f"{name}_{location}" if location else name
                molecules[key] = {"name": name, "location": location}
                continue

            # Handle regular molecules with optional location
            regular_match = re.match(r"([^()]+?)(?:\(([a-zA-Z]+)\))?$", part)
            if regular_match:
                name = regular_match.group(1).strip()
                location = regular_match.group(2)
                # Skip common small molecules that are well-handled
                if name not in ["H2O", "H(+)", "H+"]:
                    key = f"{name}_{location}" if location else name
                    if key not in molecules:  # Don't override if already found
                        molecules[key] = {"name": name, "location": location}

    return molecules


def get_rhea_ids_with_go_annotations(cache_dir: str = "cache") -> set[str]:
    """Get set of RHEA IDs that have GO annotations.

    This function reads from cached GO term mappings to find which
    RHEA reactions have associated GO terms.

    Args:
        cache_dir: Directory containing cached data files

    Returns:
        Set of RHEA IDs that have GO annotations
    """
    cache_path = Path(cache_dir)
    go_file = cache_path / "go_chebi_rhea.json"

    if not go_file.exists():
        logger.warning(f"GO annotations file not found: {go_file}")
        return set()

    try:
        with open(go_file) as f:
            data = json.load(f)

        rhea_ids = set()
        for entry in data:
            if "rhea_ids" in entry and entry["rhea_ids"]:
                rhea_ids.update(entry["rhea_ids"])

        return rhea_ids
    except Exception as e:
        logger.error(f"Error reading GO annotations: {e}")
        return set()


def fetch_rhea_reactions(
    limit: Optional[int] = None,
    rhea_ids: Optional[List[str]] = None,
    only_with_go: bool = False,
    only_with_ec: bool = False,
    include_participants: bool = True,
    participant_limit: Optional[int] = None,
    cache_dir: str = "cache",
) -> Dict[str, RheaTerm]:
    """Fetch RHEA reactions using TSV files.

    Args:
        limit: Maximum number of reactions to fetch (None = all)
        rhea_ids: Specific RHEA IDs to fetch (overrides limit)
        only_with_go: If True, only fetch reactions with GO annotations
        only_with_ec: If True, only fetch reactions with EC numbers
        include_participants: If True, fetch reaction participants (always True for TSV)
        participant_limit: Max reactions to fetch participants for (ignored for TSV)
        cache_dir: Directory to cache TSV files

    Returns:
        Dictionary mapping RHEA IDs to RheaTerm objects
    """
    etl = RheaTSVETL(cache_dir=Path(cache_dir) / "rhea_tsv")

    # Clean RHEA IDs if provided
    clean_ids = None
    if rhea_ids:
        clean_ids = [rid.replace("RHEA:", "") for rid in rhea_ids]

    # Load reactions
    reactions = etl.load_reactions_batch(rhea_ids=clean_ids, limit=limit)
    rhea_terms = etl.create_rhea_terms(reactions)

    # Convert to dictionary format
    rhea_mappings = {}
    for term in rhea_terms:
        rhea_id = term.rhea_id.replace("RHEA:", "")
        rhea_mappings[f"RHEA:{rhea_id}"] = term

    # Filter by GO if requested
    if only_with_go:
        go_rhea_ids = get_rhea_ids_with_go_annotations(cache_dir)
        rhea_mappings = {
            k: v for k, v in rhea_mappings.items()
            if k in go_rhea_ids or k.replace("RHEA:", "") in go_rhea_ids
        }

    # Filter by EC if requested (TSV ETL already prioritizes EC-annotated reactions)
    if only_with_ec:
        rhea_mappings = {
            k: v for k, v in rhea_mappings.items()
            if v.ec_numbers
        }

    return rhea_mappings


def serialize_rhea_reactions(
    rhea_mappings: Dict[str, RheaTerm], output_file: str, format: str = "json"
) -> None:
    """Serialize RHEA reactions to a file.

    Args:
        rhea_mappings: Dictionary of RHEA ID to RheaTerm objects
        output_file: Path to output file
        format: Output format ('json', 'jsonl', or 'csv')
    """
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if format == "json":
        with open(output_path, "w") as f:
            # Convert to serializable format
            data = {
                rhea_id: term.model_dump() for rhea_id, term in rhea_mappings.items()
            }
            json.dump(data, f, indent=2)

    elif format == "jsonl":
        with open(output_path, "w") as f:
            for rhea_id, term in rhea_mappings.items():
                json.dump(term.model_dump(), f)
                f.write("\n")

    elif format == "csv":
        with open(output_path, "w", newline="") as f:
            writer = csv.writer(f)
            # Write header
            writer.writerow(
                ["rhea_id", "label", "ec_numbers", "go_terms", "has_participants"]
            )

            # Write data
            for rhea_id, term in rhea_mappings.items():
                has_participants = bool(
                    term.reaction
                    and (
                        term.reaction.left_participants
                        or term.reaction.right_participants
                    )
                )
                writer.writerow(
                    [
                        rhea_id,
                        term.label,
                        ";".join(term.ec_numbers) if term.ec_numbers else "",
                        ";".join(term.go_terms) if term.go_terms else "",
                        has_participants,
                    ]
                )

    else:
        raise ValueError(f"Unsupported format: {format}")

    logger.info(f"Serialized {len(rhea_mappings)} reactions to {output_path}")


# For backward compatibility, export the key functions
__all__ = [
    "parse_location_from_label",
    "extract_molecules_from_label",
    "get_rhea_ids_with_go_annotations",
    "fetch_rhea_reactions",
    "serialize_rhea_reactions",
]
