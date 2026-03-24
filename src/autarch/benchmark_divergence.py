"""GO versus EC benchmark divergence analysis."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from autarch.classifier import ReactionClassifier
from autarch.datamodel import GoTerm
from autarch.evaluation import ec_matches_prefix


def ec_specificity_depth(
    ec_prefix: str | None,
    ec_numbers: list[str] | None = None,
) -> int:
    """Return the number of fixed EC levels for a class."""
    if ec_numbers:
        return 4
    if not ec_prefix:
        return 0
    return sum(1 for part in ec_prefix.split(".") if part and part != "-")


def classify_divergence_relation(
    go_positive_count: int,
    ec_positive_count: int,
    intersection_count: int,
) -> str:
    """Describe the set relationship between GO- and EC-derived positives."""
    if go_positive_count == 0 and ec_positive_count == 0:
        return "no_support"
    if go_positive_count > 0 and ec_positive_count == 0:
        return "go_only"
    if ec_positive_count > 0 and go_positive_count == 0:
        return "ec_only"
    if intersection_count == go_positive_count == ec_positive_count:
        return "exact_match"
    if intersection_count == go_positive_count:
        return "go_subset_of_ec"
    if intersection_count == ec_positive_count:
        return "ec_subset_of_go"
    return "divergent"


def summarize_set_divergence(
    class_name: str,
    go_id: str | None,
    ec_prefix: str | None,
    ec_numbers: list[str],
    go_positive_ids: set[str],
    ec_positive_ids: set[str],
) -> dict[str, Any]:
    """Summarize GO-versus-EC support divergence for one class."""
    intersection_ids = go_positive_ids & ec_positive_ids
    union_ids = go_positive_ids | ec_positive_ids
    go_only_ids = go_positive_ids - ec_positive_ids
    ec_only_ids = ec_positive_ids - go_positive_ids

    go_positive_count = len(go_positive_ids)
    ec_positive_count = len(ec_positive_ids)
    intersection_count = len(intersection_ids)
    union_count = len(union_ids)

    relation = classify_divergence_relation(
        go_positive_count,
        ec_positive_count,
        intersection_count,
    )

    go_covered_by_ec = (
        intersection_count / go_positive_count if go_positive_count else 0.0
    )
    ec_covered_by_go = (
        intersection_count / ec_positive_count if ec_positive_count else 0.0
    )
    jaccard = intersection_count / union_count if union_count else 0.0

    return {
        "class": class_name,
        "go_id": go_id or "",
        "ec_prefix": ec_prefix or "",
        "ec_numbers": ";".join(ec_numbers),
        "ec_specificity_depth": ec_specificity_depth(ec_prefix, ec_numbers),
        "go_positive_count": go_positive_count,
        "ec_positive_count": ec_positive_count,
        "intersection_count": intersection_count,
        "union_count": union_count,
        "go_only_count": len(go_only_ids),
        "ec_only_count": len(ec_only_ids),
        "go_covered_by_ec": go_covered_by_ec,
        "ec_covered_by_go": ec_covered_by_go,
        "jaccard": jaccard,
        "relation": relation,
    }


def compute_benchmark_divergence(cache_dir: str | Path = "cache") -> pd.DataFrame:
    """Compute GO-versus-EC divergence per classifier from cached RHEA data."""
    cache_path = Path(cache_dir)

    go_terms: dict[str, GoTerm] = {}
    with open(cache_path / "go_terms.jsonl") as handle:
        for line in handle:
            term = GoTerm(**json.loads(line))
            go_terms[term.go_id] = term

    rhea_reactions: dict[str, dict[str, Any]] = {}
    with open(cache_path / "rhea_reactions.jsonl") as handle:
        for line in handle:
            data = json.loads(line)
            rhea_reactions[data["rhea_id"]] = data

    classifier = ReactionClassifier()
    rows: list[dict[str, Any]] = []

    for class_name, reaction_class_obj in classifier.reaction_classes.items():
        go_id = getattr(reaction_class_obj, "GO_ID", None)
        ec_prefix = getattr(reaction_class_obj, "EC_NUMBER_PREFIX", None)
        ec_numbers = list(getattr(reaction_class_obj, "EC_NUMBERS", None) or [])

        if not go_id and not ec_prefix and not ec_numbers:
            continue

        go_positive_ids: set[str] = set()
        ec_positive_ids: set[str] = set()

        if go_id:
            for rhea_id, rhea_data in rhea_reactions.items():
                for reaction_go_term in rhea_data.get("go_terms", []):
                    if reaction_go_term == go_id:
                        go_positive_ids.add(rhea_id)
                        break
                    if (
                        reaction_go_term in go_terms
                        and go_id in go_terms[reaction_go_term].ancestors
                    ):
                        go_positive_ids.add(rhea_id)
                        break

        if ec_numbers:
            for rhea_id, rhea_data in rhea_reactions.items():
                reaction_ec_numbers = rhea_data.get("ec_numbers", [])
                if any(ec_number in ec_numbers for ec_number in reaction_ec_numbers):
                    ec_positive_ids.add(rhea_id)
        elif ec_prefix:
            for rhea_id, rhea_data in rhea_reactions.items():
                reaction_ec_numbers = rhea_data.get("ec_numbers", [])
                if any(
                    ec_matches_prefix(ec_number, ec_prefix)
                    for ec_number in reaction_ec_numbers
                ):
                    ec_positive_ids.add(rhea_id)

        rows.append(
            summarize_set_divergence(
                class_name,
                go_id,
                ec_prefix,
                ec_numbers,
                go_positive_ids,
                ec_positive_ids,
            )
        )

    return pd.DataFrame(rows).sort_values(
        ["relation", "jaccard", "class"],
        ascending=[True, False, True],
    )


def save_benchmark_divergence_report(
    output_dir: str | Path,
    cache_dir: str | Path = "cache",
) -> tuple[Path, Path]:
    """Write GO-versus-EC divergence reports to disk."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    divergence_df = compute_benchmark_divergence(cache_dir)
    csv_path = output_path / "benchmark_divergence.csv"
    divergence_df.to_csv(csv_path, index=False)

    summary_path = output_path / "benchmark_divergence_summary.txt"
    relation_counts = (
        divergence_df["relation"].value_counts().sort_index().to_dict()
        if not divergence_df.empty
        else {}
    )

    with open(summary_path, "w") as handle:
        handle.write("GO versus EC Divergence Summary\n")
        handle.write("=" * 60 + "\n\n")
        handle.write(f"Classes analyzed: {len(divergence_df)}\n")
        handle.write(
            "Classes with both GO and EC support: "
            f"{int(((divergence_df['go_positive_count'] > 0) & (divergence_df['ec_positive_count'] > 0)).sum())}\n\n"
        )

        handle.write("Relation counts\n")
        handle.write("-" * 30 + "\n")
        for relation, count in relation_counts.items():
            handle.write(f"  {relation}: {count}\n")

        def write_top_section(title: str, frame: pd.DataFrame, count_column: str) -> None:
            handle.write(f"\n{title}\n")
            handle.write("-" * 30 + "\n")
            if frame.empty:
                handle.write("  none\n")
                return
            for row in frame.head(15).to_dict(orient="records"):
                handle.write(
                    f"  {row['class']}: "
                    f"{count_column}={row[count_column]}, "
                    f"relation={row['relation']}, "
                    f"GO={row['go_positive_count']}, "
                    f"EC={row['ec_positive_count']}, "
                    f"Jaccard={row['jaccard']:.3f}\n"
                )

        write_top_section(
            "Largest GO-only tails",
            divergence_df.sort_values(
                ["go_only_count", "jaccard"],
                ascending=[False, True],
            ),
            "go_only_count",
        )
        write_top_section(
            "Largest EC-only tails",
            divergence_df.sort_values(
                ["ec_only_count", "jaccard"],
                ascending=[False, True],
            ),
            "ec_only_count",
        )

    return csv_path, summary_path
