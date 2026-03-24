"""Tests for GO-versus-EC benchmark divergence analysis."""

from autarch.benchmark_divergence import (
    classify_divergence_relation,
    ec_specificity_depth,
    summarize_set_divergence,
)


def test_ec_specificity_depth_counts_fixed_levels() -> None:
    """EC specificity should count only fixed prefix levels."""
    assert ec_specificity_depth("4.-.-.-") == 1
    assert ec_specificity_depth("4.2.3.-") == 3
    assert ec_specificity_depth(None, ["4.2.3.166"]) == 4
    assert ec_specificity_depth(None, []) == 0


def test_classify_divergence_relation_distinguishes_set_patterns() -> None:
    """Set relationships should map to stable qualitative labels."""
    assert classify_divergence_relation(0, 0, 0) == "no_support"
    assert classify_divergence_relation(3, 0, 0) == "go_only"
    assert classify_divergence_relation(0, 4, 0) == "ec_only"
    assert classify_divergence_relation(5, 5, 5) == "exact_match"
    assert classify_divergence_relation(4, 10, 4) == "go_subset_of_ec"
    assert classify_divergence_relation(10, 4, 4) == "ec_subset_of_go"
    assert classify_divergence_relation(10, 10, 3) == "divergent"


def test_summarize_set_divergence_reports_counts_and_ratios() -> None:
    """Divergence rows should preserve set counts and overlap metrics."""
    row = summarize_set_divergence(
        class_name="TerpeneSynthase",
        go_id="GO:0010333",
        ec_prefix="4.2.3.-",
        ec_numbers=[],
        go_positive_ids={"RHEA:1", "RHEA:2"},
        ec_positive_ids={"RHEA:1", "RHEA:2", "RHEA:3", "RHEA:4"},
    )

    assert row["relation"] == "go_subset_of_ec"
    assert row["ec_specificity_depth"] == 3
    assert row["go_positive_count"] == 2
    assert row["ec_positive_count"] == 4
    assert row["intersection_count"] == 2
    assert row["go_only_count"] == 0
    assert row["ec_only_count"] == 2
    assert row["go_covered_by_ec"] == 1.0
    assert row["ec_covered_by_go"] == 0.5
    assert row["jaccard"] == 0.5
