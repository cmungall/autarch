"""Tests for local embedding explanation helpers."""

import numpy as np

from autarch.rhea_embedding_explanations import (
    binary_classifier_scores,
    build_embedding_explanation_components,
    build_reaction_text_embedding_matrix,
    normalize_rhea_id,
)


def test_normalize_rhea_id_accepts_prefixed_and_unprefixed_values() -> None:
    """RHEA IDs should normalize to the cached prefix form."""
    assert normalize_rhea_id("10000") == "RHEA:10000"
    assert normalize_rhea_id("RHEA:10000") == "RHEA:10000"


def test_build_embedding_explanation_components_covers_label_sides_and_participants() -> None:
    """Occlusion units should include the label, each side, and each participant."""
    components = build_embedding_explanation_components(
        "ATP + H2O = ADP + phosphate",
        [{"name": "ATP"}, {"name": "H2O"}],
        [{"name": "ADP"}, {"name": "phosphate"}],
    )

    keys = [component["key"] for component in components]
    assert keys == ["label", "lhs", "rhs", "lhs_0", "lhs_1", "rhs_0", "rhs_1"]
    assert components[0]["masked_text"].startswith("reactants:")
    assert "products: ADP; phosphate" in components[0]["masked_text"]
    assert "reactants: ATP | products: ADP; phosphate" in components[4]["masked_text"]
    assert "products: ADP" in components[6]["masked_text"]
    assert "products: ADP; phosphate" not in components[6]["masked_text"]


def test_build_embedding_explanation_components_supports_participant_only_space() -> None:
    """Participant-only explanations should not reintroduce the label text."""
    components = build_embedding_explanation_components(
        "ATP + H2O = ADP + phosphate",
        [{"name": "ATP"}, {"name": "H2O"}],
        [{"name": "ADP"}, {"name": "phosphate"}],
        include_label=False,
    )

    assert components[0]["masked_text"] == "reactants: ATP; H2O | products: ADP; phosphate"
    assert not components[1]["masked_text"].startswith("ATP + H2O = ADP + phosphate")


def test_build_reaction_text_embedding_matrix_supports_lexical_backend() -> None:
    """Explanation texts should be embeddable without live API calls."""
    matrix = build_reaction_text_embedding_matrix(
        [
            "ATP + H2O = ADP + phosphate | reactants: ATP; H2O | products: ADP; phosphate",
            "reactants: ATP; H2O | products: ADP; phosphate",
        ],
        use_linkml_store=False,
        n_features=32,
    )

    assert matrix.shape == (2, 32)
    assert np.linalg.norm(matrix[0]) > 0.0


def test_binary_classifier_scores_prefers_decision_function() -> None:
    """Models with a decision function should expose that score directly."""

    class DummyModel:
        def decision_function(self, X: np.ndarray) -> np.ndarray:
            return X.sum(axis=1)

    scores = binary_classifier_scores(DummyModel(), np.array([[1.0, 2.0], [0.5, -1.0]]))

    assert np.allclose(scores, [3.0, -0.5])
