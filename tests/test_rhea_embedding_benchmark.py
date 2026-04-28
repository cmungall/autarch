"""Tests for the rule-embedding benchmark helpers."""

from collections import Counter

import numpy as np

from autarch.rhea_embedding_benchmark import (
    build_drfp_feature_spec,
    build_rule_benchmark_feature_matrices,
    build_rule_benchmark_feature_specs,
    select_benchmark_classes,
    select_benchmark_values,
)


def test_build_rule_benchmark_feature_matrices_derives_expected_spaces() -> None:
    """Benchmark spaces should preserve reaction, diff, and symmetric features."""
    reaction = np.array([[1.0, 2.0]], dtype=np.float32)
    participant_only = np.array([[0.5, 1.5]], dtype=np.float32)
    lhs = np.array([[2.0, 4.0]], dtype=np.float32)
    rhs = np.array([[5.0, 10.0]], dtype=np.float32)
    base = {
        "reaction": reaction,
        "reaction_participants_only": participant_only,
        "lhs": lhs,
        "rhs": rhs,
        "rhs_minus_lhs": rhs - lhs,
    }

    feature_matrices = build_rule_benchmark_feature_matrices(base)

    assert np.allclose(feature_matrices["reaction"], reaction)
    assert np.allclose(
        feature_matrices["reaction_participants_only"],
        participant_only,
    )
    assert np.allclose(feature_matrices["rhs_minus_lhs"], [[3.0, 6.0]])
    assert np.allclose(
        feature_matrices["symmetric_sum_sqdiff"],
        [[3.5, 7.0, 9.0, 36.0]],
    )


def test_build_rule_benchmark_feature_specs_preserves_row_indices() -> None:
    """Feature specs should carry the shared row subset explicitly."""
    base = {
        "reaction": np.array([[1.0, 2.0]], dtype=np.float32),
        "reaction_participants_only": np.array([[0.5, 1.5]], dtype=np.float32),
        "lhs": np.array([[2.0, 4.0]], dtype=np.float32),
        "rhs": np.array([[5.0, 10.0]], dtype=np.float32),
        "rhs_minus_lhs": np.array([[3.0, 6.0]], dtype=np.float32),
    }

    feature_specs = build_rule_benchmark_feature_specs(
        base,
        row_indices=np.array([7], dtype=np.int32),
    )

    assert np.array_equal(feature_specs["reaction"]["row_indices"], np.array([7]))
    assert np.allclose(feature_specs["reaction"]["matrix"], [[1.0, 2.0]])


def test_build_drfp_feature_spec_uses_valid_reaction_smiles_subset() -> None:
    """DRFP should skip rows without concrete reaction SMILES."""
    import pandas as pd

    df = pd.DataFrame(
        {
            "reaction_smiles": [
                "CCO.O>>CC=O.O",
                None,
                "CC(=O)O.O>>CC(=O)[O-].[H+]",
            ]
        }
    )

    feature_spec = build_drfp_feature_spec(df, n_folded_length=128)

    assert feature_spec["matrix"].shape == (2, 128)
    assert np.array_equal(feature_spec["row_indices"], np.array([0, 2], dtype=np.int32))


def test_select_benchmark_classes_filters_and_sorts_by_support() -> None:
    """Class selection should apply thresholds before max-class truncation."""
    class_counts = Counter(
        {
            "Transferase": 40,
            "Hydrolase": 25,
            "Lyase": 12,
            "Oxidoreductase": 40,
        }
    )

    selected = select_benchmark_classes(
        class_counts,
        min_positives=20,
        max_classes=3,
    )

    assert selected == ["Oxidoreductase", "Transferase", "Hydrolase"]


def test_select_benchmark_values_preserves_declared_order() -> None:
    """Requested model/space subsets should respect the canonical output order."""
    selected = select_benchmark_values(
        requested=["rhs_minus_lhs", "reaction_participants_only"],
        available=[
            "reaction",
            "reaction_participants_only",
            "rhs_minus_lhs",
            "symmetric_sum_sqdiff",
        ],
        value_type="feature space",
    )

    assert selected == ["reaction_participants_only", "rhs_minus_lhs"]
