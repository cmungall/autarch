"""Curator-facing objective scores for balancing performance and parsimony.

This module combines evaluation metrics with code-complexity reports to produce
an internal curation objective. The goal is to favor classifiers that perform
well, remain simple, and avoid over-interpreting tiny support counts.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import pandas as pd

from autarch.complexity import ComplexityReport

DEFAULT_COMPLEXITY_WEIGHT = 0.15


def normalize_classifier_name(name: str) -> str:
    """Normalize classifier names for joining across reports."""
    return re.sub(r"[^a-z0-9]+", "", name.lower())


def wilson_lower_bound(successes: int, total: int, z: float = 1.96) -> float:
    """Return the Wilson-score lower bound for a binomial proportion.

    This provides a conservative estimate that naturally penalizes small sample
    sizes, which is useful for tiny-support reaction classes.

    Examples:
        >>> round(wilson_lower_bound(1, 1), 3)
        0.207
        >>> round(wilson_lower_bound(10, 10), 3)
        0.722
    """
    if total <= 0:
        return 0.0

    proportion = successes / total
    z_squared = z * z
    denominator = 1.0 + (z_squared / total)
    center = proportion + (z_squared / (2.0 * total))
    margin = z * math.sqrt(
        (proportion * (1.0 - proportion) + (z_squared / (4.0 * total))) / total
    )
    return max(0.0, (center - margin) / denominator)


def harmonic_mean(x: float, y: float) -> float:
    """Return the harmonic mean of two non-negative numbers."""
    if x + y == 0.0:
        return 0.0
    return 2.0 * x * y / (x + y)


def build_curation_objective_df(
    eval_df: pd.DataFrame,
    complexity_report: ComplexityReport,
    complexity_weight: float = DEFAULT_COMPLEXITY_WEIGHT,
) -> pd.DataFrame:
    """Build a curation report that rewards stable performance and parsimony.

    The score intentionally keeps complexity as a mild tiebreaker. Performance
    first passes through a conservative transform using Wilson lower bounds on
    precision and recall, which strongly penalizes tiny-support "perfect"
    classifiers.
    """
    required_columns = {"class", "tp", "fp", "fn", "f1_score", "precision", "recall"}
    missing = required_columns - set(eval_df.columns)
    if missing:
        missing_str = ", ".join(sorted(missing))
        raise ValueError(f"Evaluation results missing required columns: {missing_str}")

    complexity_rows = []
    for classifier in complexity_report.classifiers:
        complexity_rows.append(
            {
                "classifier_key": normalize_classifier_name(classifier.name),
                "complexity_score": classifier.complexity_score,
                "cyclomatic_complexity": classifier.cyclomatic_complexity,
                "cc_rank": classifier.cc_rank,
                "complexity_loc": classifier.loc,
                "complexity_method_loc": classifier.method_loc,
                "complexity_if_statements": classifier.if_statements,
                "complexity_patterns": classifier.patterns,
                "complexity_chebi_ids": classifier.chebi_ids,
                "complexity_exclusions": classifier.exclusions,
            }
        )
    complexity_df = pd.DataFrame(complexity_rows)

    merged = eval_df.copy()
    merged["classifier_key"] = merged["class"].map(normalize_classifier_name)
    merged = merged.merge(complexity_df, on="classifier_key", how="left")

    unresolved = merged[merged["complexity_score"].isna()]["class"].tolist()
    if unresolved:
        unresolved_str = ", ".join(sorted(unresolved))
        raise ValueError(
            f"Could not match complexity metrics for classes: {unresolved_str}"
        )

    merged["positive_support"] = merged["tp"] + merged["fn"]
    merged["predicted_positive_support"] = merged["tp"] + merged["fp"]
    merged["precision_lb95"] = [
        wilson_lower_bound(tp, tp + fp)
        for tp, fp in zip(merged["tp"], merged["fp"])
    ]
    merged["recall_lb95"] = [
        wilson_lower_bound(tp, tp + fn)
        for tp, fn in zip(merged["tp"], merged["fn"])
    ]
    merged["conservative_f1"] = [
        harmonic_mean(precision_lb, recall_lb)
        for precision_lb, recall_lb in zip(
            merged["precision_lb95"], merged["recall_lb95"]
        )
    ]
    merged["complexity_percentile"] = merged["complexity_score"].rank(
        pct=True, method="average"
    )
    merged["f1_optimism_gap"] = merged["f1_score"] - merged["conservative_f1"]
    merged["curation_objective"] = (
        merged["conservative_f1"]
        - complexity_weight * merged["complexity_percentile"]
    )
    merged["small_support_flag"] = merged["positive_support"] < 5
    merged["overfit_risk_flag"] = (
        (merged["positive_support"] <= 10) & (merged["f1_optimism_gap"] >= 0.25)
    )

    objective_median = merged["curation_objective"].median()
    merged["parsimony_pressure_flag"] = (
        (merged["complexity_percentile"] >= 0.75)
        & (merged["curation_objective"] < objective_median)
    )

    return merged.sort_values(
        ["curation_objective", "conservative_f1", "f1_score"],
        ascending=[False, False, False],
    )


def save_curation_objective_report(
    eval_df: pd.DataFrame,
    complexity_report: ComplexityReport,
    output_dir: Path,
    complexity_weight: float = DEFAULT_COMPLEXITY_WEIGHT,
) -> pd.DataFrame:
    """Save curator-facing objective reports as CSV and text summary."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    report_df = build_curation_objective_df(
        eval_df=eval_df,
        complexity_report=complexity_report,
        complexity_weight=complexity_weight,
    )

    csv_path = output_dir / "curation_objective.csv"
    report_df.to_csv(csv_path, index=False)

    summary_lines = [
        "ARCTURUS CURATION OBJECTIVE SUMMARY",
        "=" * 60,
        "",
        f"Complexity weight: {complexity_weight:.2f}",
        "Objective = conservative_f1 - complexity_weight * complexity_percentile",
        "",
        "Top 5 curator objective scores:",
    ]
    for _, row in report_df.head(5).iterrows():
        summary_lines.append(
            f"  - {row['class']}: objective={row['curation_objective']:.3f}, "
            f"conservative_f1={row['conservative_f1']:.3f}, "
            f"raw_f1={row['f1_score']:.3f}, support={int(row['positive_support'])}"
        )

    summary_lines.append("")
    summary_lines.append("Highest overfitting-risk classes:")
    risky_df = report_df.sort_values(
        ["f1_optimism_gap", "positive_support"], ascending=[False, True]
    ).head(5)
    for _, row in risky_df.iterrows():
        summary_lines.append(
            f"  - {row['class']}: gap={row['f1_optimism_gap']:.3f}, "
            f"raw_f1={row['f1_score']:.3f}, conservative_f1={row['conservative_f1']:.3f}, "
            f"support={int(row['positive_support'])}"
        )

    summary_lines.append("")
    summary_lines.append("High-complexity simplification candidates:")
    simplify_df = report_df[report_df["parsimony_pressure_flag"]].sort_values(
        ["complexity_percentile", "curation_objective"],
        ascending=[False, True],
    ).head(5)
    for _, row in simplify_df.iterrows():
        summary_lines.append(
            f"  - {row['class']}: objective={row['curation_objective']:.3f}, "
            f"complexity_score={row['complexity_score']:.1f}, "
            f"cc={int(row['cyclomatic_complexity'])}, raw_f1={row['f1_score']:.3f}"
        )

    summary_path = output_dir / "curation_objective_summary.txt"
    with open(summary_path, "w") as handle:
        handle.write("\n".join(summary_lines) + "\n")

    return report_df
