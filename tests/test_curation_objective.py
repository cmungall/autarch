"""Tests for curator-facing objective scoring."""

import pandas as pd

from autarch.complexity import ClassifierComplexity, ComplexityReport
from autarch.curation_objective import (
    build_curation_objective_df,
    wilson_lower_bound,
)


def test_wilson_lower_bound_penalizes_small_perfect_support() -> None:
    """Tiny-support perfect scores should have much lower confidence."""
    assert wilson_lower_bound(1, 1) < wilson_lower_bound(10, 10)
    assert wilson_lower_bound(10, 10) < 1.0


def test_curation_objective_prefers_simpler_rule_when_f1_matches() -> None:
    """When performance matches, the simpler classifier should rank higher."""
    eval_df = pd.DataFrame(
        [
            {
                "class": "SimpleRule",
                "tp": 10,
                "fp": 2,
                "tn": 90,
                "fn": 2,
                "precision": 10 / 12,
                "recall": 10 / 12,
                "f1_score": 10 / 12,
            },
            {
                "class": "ComplexRule",
                "tp": 10,
                "fp": 2,
                "tn": 90,
                "fn": 2,
                "precision": 10 / 12,
                "recall": 10 / 12,
                "f1_score": 10 / 12,
            },
        ]
    )
    complexity_report = ComplexityReport(
        classifiers=[
            ClassifierComplexity(
                name="simple_rule",
                file_path="simple_rule.py",
                cyclomatic_complexity=8,
                chebi_ids=1,
                exclusions=1,
                if_statements=3,
            ),
            ClassifierComplexity(
                name="complex_rule",
                file_path="complex_rule.py",
                cyclomatic_complexity=30,
                chebi_ids=12,
                exclusions=5,
                if_statements=15,
            ),
        ]
    )

    report_df = build_curation_objective_df(eval_df, complexity_report)
    objective_by_class = dict(
        zip(report_df["class"], report_df["curation_objective"], strict=False)
    )

    assert objective_by_class["SimpleRule"] > objective_by_class["ComplexRule"]


def test_curation_objective_penalizes_tiny_support_perfect_scores() -> None:
    """Conservative F1 should push down one-example perfect classifiers."""
    eval_df = pd.DataFrame(
        [
            {
                "class": "TinyPerfect",
                "tp": 1,
                "fp": 0,
                "tn": 99,
                "fn": 0,
                "precision": 1.0,
                "recall": 1.0,
                "f1_score": 1.0,
            },
            {
                "class": "StablePerfect",
                "tp": 20,
                "fp": 0,
                "tn": 80,
                "fn": 0,
                "precision": 1.0,
                "recall": 1.0,
                "f1_score": 1.0,
            },
        ]
    )
    complexity_report = ComplexityReport(
        classifiers=[
            ClassifierComplexity(
                name="tiny_perfect",
                file_path="tiny_perfect.py",
                cyclomatic_complexity=10,
            ),
            ClassifierComplexity(
                name="stable_perfect",
                file_path="stable_perfect.py",
                cyclomatic_complexity=10,
            ),
        ]
    )

    report_df = build_curation_objective_df(eval_df, complexity_report)
    row_by_class = {
        row["class"]: row for _, row in report_df.iterrows()
    }

    assert (
        row_by_class["TinyPerfect"]["conservative_f1"]
        < row_by_class["StablePerfect"]["conservative_f1"]
    )
    assert (
        row_by_class["TinyPerfect"]["curation_objective"]
        < row_by_class["StablePerfect"]["curation_objective"]
    )
