"""ETL functions for mining GO classifications and their mappings to reaction/enzyme databases."""

import csv
import json
import logging
from datetime import datetime
from typing import Optional
from pathlib import Path

from oaklib import get_adapter  # type: ignore[import-untyped]
from oaklib.datamodels.vocabulary import IS_A  # type: ignore[import-untyped]

from autarch.datamodel import GoTerm

logger = logging.getLogger(__name__)


def _dedupe_preserve_order(values: list[str]) -> list[str]:
    """Deduplicate a list while preserving first-seen order."""
    seen: set[str] = set()
    result: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        result.append(value)
    return result


def fetch_go_enzyme_mappings(go_terms: Optional[list[str]] = None) -> dict[str, GoTerm]:
    """Fetch GO terms with their mappings to RHEA reactions and EC numbers.

    Args:
        go_terms: Optional list of GO term IDs to fetch.
                 If None, fetches all molecular function terms.

    Returns:
        Dictionary mapping GO IDs to GoTerm objects.

    Examples:
        >>> terms = fetch_go_enzyme_mappings(go_terms=["GO:0004553"])
        >>> isinstance(terms, dict)
        True
        >>> if terms:
        ...     term = next(iter(terms.values()))
        ...     isinstance(term, GoTerm)
        ... else:
        ...     True
        True
    """
    results = {}

    # Use the Gene Ontology adapter
    adapter = get_adapter("sqlite:obo:go")

    # Get all GO molecular function terms if no filter provided
    if go_terms is None:
        # Get all descendants of molecular_function (GO:0003674) - reflexive by default
        go_terms = list(adapter.descendants("GO:0003674", predicates=[IS_A]))

    # Get all labels and definitions at once - more efficient
    labels = dict(adapter.labels(go_terms))
    definitions = {go_id: defn for go_id, defn, _ in adapter.definitions(go_terms)}

    # Initialize all GO terms with ancestors
    for go_id in go_terms:
        # Get ancestors (closure) for this term
        ancestors = list(adapter.ancestors(go_id, predicates=[IS_A]))

        # filter out non-GO (e.g. BFO) terms
        ancestors = [a for a in ancestors if a.startswith("GO:")]
        ancestors = _dedupe_preserve_order(ancestors)

        results[go_id] = GoTerm(
            go_id=go_id,
            label=labels.get(go_id, ""),
            definition=definitions.get(go_id, ""),
            ancestors=ancestors,
        )

    # Get mappings for all sources - more efficient to pass all curies at once
    sources_to_fetch = ["RHEA", "EC"]  # Can easily add more sources here

    for source in sources_to_fetch:
        for mapping in adapter.sssom_mappings(curies=go_terms, source=source):
            go_id = mapping.subject_id
            if go_id not in results:
                continue

            if (
                source == "RHEA"
                and mapping.object_id
                and mapping.object_id.startswith("RHEA:")
            ):
                rhea_id = mapping.object_id.split("_")[0]  # Remove direction suffix
                if rhea_id not in results[go_id].rhea_ids:
                    results[go_id].rhea_ids.append(rhea_id)
            elif (
                source == "EC"
                and mapping.object_id
                and mapping.object_id.startswith("EC:")
            ):
                ec_number = mapping.object_id.replace("EC:", "")
                if ec_number not in results[go_id].ec_numbers:
                    results[go_id].ec_numbers.append(ec_number)

    # Filter to only terms with mappings if not specifically requested
    if go_terms:
        # Keep all requested terms
        return results
    else:
        # Only keep terms that have some mappings
        return {
            go_id: term
            for go_id, term in results.items()
            if term.rhea_ids or term.ec_numbers
        }


def mine_go_enzyme_classifications(output_file: Optional[str] = None) -> dict:
    """Mine GO enzyme classifications with their RHEA and EC mappings.

    Args:
        output_file: Optional path to save JSON output

    Returns:
        Dictionary containing GO terms and statistics.

    Examples:
        >>> result = mine_go_enzyme_classifications()
        >>> "go_terms" in result
        True
        >>> "total_terms" in result
        True
        >>> "timestamp" in result
        True
    """
    logger.info("Starting GO enzyme classification mining...")

    # Fetch all mappings
    go_terms = fetch_go_enzyme_mappings()
    logger.info(f"Found {len(go_terms)} GO terms with mappings")

    # Convert to serializable format using Pydantic's model_dump
    terms_dict = {
        go_id: {
            **term.model_dump(),
            "rhea_count": len(term.rhea_ids),
            "ec_count": len(term.ec_numbers),
        }
        for go_id, term in go_terms.items()
    }

    # Prepare result
    result = {
        "go_terms": terms_dict,
        "total_terms": len(go_terms),
        "total_rhea_mappings": sum(len(t.rhea_ids) for t in go_terms.values()),
        "total_ec_mappings": sum(len(t.ec_numbers) for t in go_terms.values()),
        "timestamp": datetime.now().isoformat(),
    }

    # Save to file if requested
    if output_file:
        with open(output_file, "w") as f:
            json.dump(result, f, indent=2)
        logger.info(f"Saved results to {output_file}")

    return result


def serialize_go_terms(
    go_terms: dict[str, GoTerm], output_file: str, format: str = "jsonl"
) -> None:
    """Serialize GO terms to JSONL or CSV format.

    Args:
        go_terms: Dictionary of GO terms
        output_file: Path to output file
        format: Output format - "jsonl" or "csv"

    Examples:
        >>> terms = {"GO:0004553": GoTerm(go_id="GO:0004553", label="test")}
        >>> import tempfile
        >>> with tempfile.NamedTemporaryFile(suffix=".jsonl", delete=False) as f:
        ...     serialize_go_terms(terms, f.name, "jsonl")
        >>> Path(f.name).exists()
        True
    """
    output_path = Path(output_file)

    if format == "jsonl":
        with open(output_path, "w") as f:
            for go_id, term in go_terms.items():
                f.write(json.dumps(term.model_dump()) + "\n")
        logger.info(f"Wrote {len(go_terms)} terms to {output_path}")

    elif format == "csv":
        if not go_terms:
            logger.warning("No terms to write")
            return

        # Get all unique fields from the terms
        fieldnames = [
            "go_id",
            "label",
            "definition",
            "rhea_ids",
            "ec_numbers",
            "ancestors",
        ]

        with open(output_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()

            for go_id, term in go_terms.items():
                row = term.model_dump()
                # Convert lists to pipe-separated strings for CSV
                row["rhea_ids"] = "|".join(row.get("rhea_ids", []))
                row["ec_numbers"] = "|".join(row.get("ec_numbers", []))
                row["ancestors"] = "|".join(row.get("ancestors", []))
                writer.writerow(row)

        logger.info(f"Wrote {len(go_terms)} terms to {output_path}")
    else:
        raise ValueError(f"Unsupported format: {format}. Use 'jsonl' or 'csv'")
