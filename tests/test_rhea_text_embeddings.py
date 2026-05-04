"""Tests for RHEA text embedding utilities."""

import json
from typing import cast

import numpy as np
import pytest
import autarch.rhea_text_embeddings as rte

from autarch.rhea_text_embeddings import (
    build_reaction_smiles,
    build_bidirectional_distance_matrix,
    build_ec_number_entries,
    build_go_annotations,
    build_reaction_embedding_text,
    build_reaction_side_embedding_text,
    build_text_feature_matrix,
    classify_rule_status,
    classify_participant_bucket,
    extract_ec_hierarchy,
    participant_display_text,
    project_embedding_matrix,
    save_rhea_browser_html,
    tokenize_label_text,
)


def test_tokenize_label_text_includes_word_and_character_features() -> None:
    """Tokenization should retain both lexical and local character signal."""
    tokens = tokenize_label_text("ATP + H2O = ADP + phosphate")
    assert "w:atp" in tokens
    assert "w:phosphate" in tokens
    assert any(token.startswith("c:at") for token in tokens)


def test_text_feature_matrix_is_deterministic() -> None:
    """Hashed lexical embeddings must be stable across repeated calls."""
    texts = ["ATP + H2O = ADP + phosphate", "glucose + NAD(+) = gluconolactone + NADH"]
    matrix_a = build_text_feature_matrix(texts, n_features=128)
    matrix_b = build_text_feature_matrix(texts, n_features=128)
    assert matrix_a.shape == (2, 128)
    assert np.allclose(matrix_a, matrix_b)


def test_project_embedding_matrix_returns_2d_coordinates() -> None:
    """Projection should always return the requested number of coordinates."""
    matrix = build_text_feature_matrix(
        ["a = b", "c = d", "e = f"],
        n_features=64,
    )
    coords = project_embedding_matrix(matrix, n_components=2)
    assert coords.shape == (3, 2)


def test_participant_display_text_renders_structured_metadata() -> None:
    """Participant display text should preserve key stoichiometric metadata."""
    text = participant_display_text(
        {
            "name": "RNA",
            "chebi_id": "CHEBI:33697",
            "count": 2,
            "polymer_type": "rna",
            "polymer_index": "n+1",
        }
    )
    assert "2 RNA" in text
    assert "[CHEBI:33697]" in text
    assert "{rna, n+1}" in text


def test_build_reaction_embedding_text_includes_participant_context() -> None:
    """Embedding text should combine the reaction label and participant summaries."""
    text = build_reaction_embedding_text(
        "ATP + H2O = ADP + phosphate",
        [{"name": "ATP"}, {"name": "H2O"}],
        [{"name": "ADP"}, {"name": "phosphate"}],
    )
    assert text.startswith("ATP + H2O = ADP + phosphate")
    assert "reactants: ATP; H2O" in text
    assert "products: ADP; phosphate" in text


def test_build_reaction_embedding_text_can_exclude_label() -> None:
    """Participant-only embedding text should omit the reaction label."""
    text = build_reaction_embedding_text(
        "ATP + H2O = ADP + phosphate",
        [{"name": "ATP"}, {"name": "H2O"}],
        [{"name": "ADP"}, {"name": "phosphate"}],
        include_label=False,
    )
    assert not text.startswith("ATP + H2O = ADP + phosphate")
    assert text == "reactants: ATP; H2O | products: ADP; phosphate"


def test_build_reaction_smiles_expands_integer_stoichiometry() -> None:
    """Reaction SMILES should repeat fixed-count participants and preserve sides."""
    reaction_smiles = build_reaction_smiles(
        [{"smiles": "O", "count": 2}, {"smiles": "CCO"}],
        [{"smiles": "CC=O"}, {"smiles": "O"}],
    )
    assert reaction_smiles == "O.O.CCO>>CC=O.O"


def test_build_reaction_smiles_rejects_symbolic_stoichiometry() -> None:
    """Polymer-style symbolic stoichiometry should not be forced into reaction SMILES."""
    reaction_smiles = build_reaction_smiles(
        [{"smiles": "CCO", "stoichiometry": "n"}],
        [{"smiles": "CC=O"}],
    )
    assert reaction_smiles is None


def test_build_rhea_embedding_space_df_marks_drfp_subset(tmp_path) -> None:
    """DRFP coordinates should only appear for reactions with complete reaction SMILES."""
    pytest.importorskip("drfp")
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    records = [
        {
            "rhea_id": "RHEA:30000",
            "label": "ethanol oxidation",
            "go_terms": [],
            "ec_numbers": [],
            "reaction": {
                "left_participants": [
                    {"name": "ethanol", "smiles": "CCO", "count": 1},
                    {"name": "water", "smiles": "O", "count": 1},
                ],
                "right_participants": [
                    {"name": "acetaldehyde", "smiles": "CC=O", "count": 1},
                    {"name": "water", "smiles": "O", "count": 1},
                ],
            },
        },
        {
            "rhea_id": "RHEA:30004",
            "label": "partial record",
            "go_terms": [],
            "ec_numbers": [],
            "reaction": {
                "left_participants": [
                    {"name": "ethanol", "count": 1},
                ],
                "right_participants": [
                    {"name": "acetaldehyde", "smiles": "CC=O", "count": 1},
                ],
            },
        },
    ]
    (cache_dir / "rhea_reactions.jsonl").write_text(
        "\n".join(json.dumps(record) for record in records) + "\n"
    )
    (cache_dir / "go_terms.jsonl").write_text("")

    df = rte.build_rhea_embedding_space_df(
        cache_dir=cache_dir,
        n_features=64,
        use_linkml_store=False,
    ).set_index("rhea_id")

    assert bool(df.loc["RHEA:30000", "has__reaction_drfp"])
    assert not np.isnan(cast(float, df.loc["RHEA:30000", "x__reaction_drfp"]))
    assert not np.isnan(cast(float, df.loc["RHEA:30000", "y__reaction_drfp"]))
    assert not bool(df.loc["RHEA:30004", "has__reaction_drfp"])
    assert np.isnan(cast(float, df.loc["RHEA:30004", "x__reaction_drfp"]))
    assert np.isnan(cast(float, df.loc["RHEA:30004", "y__reaction_drfp"]))


def test_build_reaction_side_embedding_text_preserves_directional_side() -> None:
    """Side-specific text should preserve the chosen LHS/RHS label."""
    text = build_reaction_side_embedding_text(
        "reactants",
        [{"name": "ATP"}, {"name": "H2O"}],
    )
    assert text == "reactants: ATP; H2O"


def test_build_bidirectional_distance_matrix_is_swap_invariant() -> None:
    """Swap-equivalent reactions should collapse to zero bidirectional distance."""
    lhs = np.array([[1.0, 0.0], [0.0, 1.0]], dtype=np.float32)
    rhs = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.float32)
    distance_matrix = build_bidirectional_distance_matrix(lhs, rhs)
    assert distance_matrix.shape == (2, 2)
    assert np.isclose(distance_matrix[0, 1], 0.0)
    assert np.isclose(distance_matrix[1, 0], 0.0)


def test_extract_ec_hierarchy_returns_major_and_subclass_levels() -> None:
    """EC hierarchy extraction should normalize mixed EC mappings."""
    majors, subclasses = extract_ec_hierarchy(["2.7.1.1", "1.14.13.39", "2.-.-.-"])
    assert majors == ["1", "2"]
    assert subclasses == ["1.14", "2.7"]


def test_classify_participant_bucket_uses_stable_ranges() -> None:
    """Participant buckets should use a small fixed set of ranges."""
    assert classify_participant_bucket(2) == "1-2"
    assert classify_participant_bucket(4) == "3-4"
    assert classify_participant_bucket(6) == "5-6"
    assert classify_participant_bucket(7) == "7+"


def test_build_go_annotations_includes_direct_terms_and_closure_labels() -> None:
    """GO annotation expansion should preserve direct labels and ancestor closure."""
    lookup = {
        "GO:0016740": {"label": "transferase activity", "ancestors": ["GO:0003824"]},
        "GO:0003824": {"label": "catalytic activity", "ancestors": ["GO:0003674"]},
        "GO:0003674": {"label": "molecular_function", "ancestors": []},
    }
    direct_entries, closure_ids, closure_labels = build_go_annotations(
        ["GO:0016740"],
        lookup,
    )
    assert direct_entries == [
        {"id": "GO:0016740", "label": "transferase activity"},
    ]
    assert closure_ids == ["GO:0016740", "GO:0003824", "GO:0003674"]
    assert "transferase activity" in closure_labels


def test_build_ec_number_entries_attaches_labels_when_available() -> None:
    """EC entries should preserve ids and fill labels from the lookup."""
    entries = build_ec_number_entries(
        ["2.5.1.17", "3.6.1.3"],
        {"2.5.1.17": "corrinoid adenosyltransferase"},
    )
    assert entries == [
        {"id": "2.5.1.17", "label": "corrinoid adenosyltransferase"},
        {"id": "3.6.1.3", "label": ""},
    ]


def test_classify_rule_status_separates_unannotated_from_ec_backed_gaps() -> None:
    """GO-missing agent-only matches should separate unannotated from EC-backed."""
    assert (
        classify_rule_status(
            [], ["ATPHydrolysis"], has_go_terms=False, has_ec_numbers=False
        )
        == "unknown_positive"
    )
    assert (
        classify_rule_status(
            [], ["ATPHydrolysis"], has_go_terms=False, has_ec_numbers=True
        )
        == "ec_backed_positive"
    )
    assert (
        classify_rule_status(
            [], ["ATPHydrolysis"], has_go_terms=True, has_ec_numbers=False
        )
        == "inferred_only"
    )


def test_save_rhea_browser_html_writes_browser_page(tmp_path) -> None:
    """Browser export should emit a standalone page with the linked explorer."""
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    results_dir = tmp_path / "results"
    results_dir.mkdir()
    sample_record = {
        "rhea_id": "RHEA:10000",
        "label": "ATP + H2O = ADP + phosphate",
        "go_terms": ["GO:0016887"],
        "ec_numbers": ["3.6.1.3"],
        "reaction": {
            "left_participants": [
                {
                    "name": "ATP",
                    "chebi_id": "CHEBI:30616",
                    "count": 1,
                    "smiles": "OP(=O)(O)OP(=O)(O)O",
                },
                {
                    "name": "H2O",
                    "chebi_id": "CHEBI:15377",
                    "count": 1,
                    "smiles": "O",
                },
            ],
            "right_participants": [
                {
                    "name": "ADP",
                    "chebi_id": "CHEBI:16761",
                    "count": 1,
                    "smiles": "OP(=O)(O)O",
                },
                {
                    "name": "phosphate",
                    "chebi_id": "CHEBI:18367",
                    "count": 1,
                    "smiles": "O=P(O)(O)O",
                },
            ],
        },
    }
    (cache_dir / "rhea_reactions.jsonl").write_text(json.dumps(sample_record) + "\n")
    (cache_dir / "go_terms.jsonl").write_text(
        json.dumps(
            {
                "go_id": "GO:0016887",
                "label": "ATP hydrolysis activity",
                "ancestors": ["GO:0016817", "GO:0003824"],
            }
        )
        + "\n"
    )
    (results_dir / "detailed_predictions.csv").write_text(
        "class,go_term,rhea_id,label,prediction,actual,outcome,explanation\n"
        "ATPHydrolysis,GO:0016887,RHEA:10000,ATP + H2O = ADP + phosphate,positive,positive,TP,matched\n"
    )

    output_file = tmp_path / "rhea_browser.html"
    save_rhea_browser_html(
        output_file,
        cache_dir=cache_dir,
        results_dir=results_dir,
        n_features=64,
        use_linkml_store=False,
    )

    html_text = output_file.read_text()
    assert "RHEA reaction browser with linked text embeddings" in html_text
    assert "RHEA:10000" in html_text
    assert "Structured record" in html_text
    assert "embedding-plot" in html_text
    assert "ATP hydrolysis activity" in html_text
    assert "All rule statuses" in html_text
    assert "Asserted classes:" in html_text
    assert "Agent classification:" in html_text
    assert "Use the Rule status facet to isolate mismatches" in html_text
    assert "ATPHydrolysis" in html_text
    assert "atphydrolysis.html" in html_text
    assert "Bidirectional" in html_text
    assert "Reaction SMILES (DRFP)" in html_text
    assert "embedding-space-coverage" in html_text
    assert "Embedding-space guide" in html_text
    assert "Chemistry-native:" in html_text
    assert "Coverage note:" in html_text
    assert "Reaction SMILES" in html_text
    assert "RHS-LHS diff" in html_text
    assert "Reactant-side embedding text" in html_text


def test_save_rhea_browser_html_surfaces_unannotated_candidate_status(
    tmp_path, monkeypatch
) -> None:
    """GO-missing reactions with no EC support should be labeled as unannotated."""
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    sample_record = {
        "rhea_id": "RHEA:20000",
        "label": "ATP + H2O = ADP + phosphate",
        "go_terms": [],
        "ec_numbers": [],
        "reaction": {
            "left_participants": [
                {"name": "ATP", "chebi_id": "CHEBI:30616", "count": 1},
                {"name": "H2O", "chebi_id": "CHEBI:15377", "count": 1},
            ],
            "right_participants": [
                {"name": "ADP", "chebi_id": "CHEBI:16761", "count": 1},
                {"name": "phosphate", "chebi_id": "CHEBI:18367", "count": 1},
            ],
        },
    }
    (cache_dir / "rhea_reactions.jsonl").write_text(json.dumps(sample_record) + "\n")
    monkeypatch.setattr(
        rte,
        "infer_rule_matches_for_go_missing",
        lambda records: {"RHEA:20000": ["ATPHydrolysis"]},
    )

    output_file = tmp_path / "rhea_browser.html"
    save_rhea_browser_html(
        output_file,
        cache_dir=cache_dir,
        results_dir=None,
        n_features=64,
        use_linkml_store=False,
    )

    html_text = output_file.read_text()
    assert "Unannotated candidate" in html_text
    assert "unannotated candidate (1)" in html_text


def test_save_rhea_browser_html_surfaces_ec_backed_go_missing_status(
    tmp_path, monkeypatch
) -> None:
    """GO-missing reactions with EC support should get their own status."""
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    sample_record = {
        "rhea_id": "RHEA:20001",
        "label": "farnesyl diphosphate + H2O = germacradien-6-ol + diphosphate",
        "go_terms": [],
        "ec_numbers": ["4.2.3.166"],
        "reaction": {
            "left_participants": [
                {
                    "name": "farnesyl diphosphate",
                    "chebi_id": "CHEBI:66914",
                    "count": 1,
                },
                {"name": "H2O", "chebi_id": "CHEBI:15377", "count": 1},
            ],
            "right_participants": [
                {
                    "name": "germacradien-6-ol",
                    "chebi_id": "CHEBI:133670",
                    "count": 1,
                },
                {"name": "diphosphate", "chebi_id": "CHEBI:33019", "count": 1},
            ],
        },
    }
    (cache_dir / "rhea_reactions.jsonl").write_text(json.dumps(sample_record) + "\n")
    monkeypatch.setattr(
        rte,
        "infer_rule_matches_for_go_missing",
        lambda records: {"RHEA:20001": ["TerpeneSynthase"]},
    )

    output_file = tmp_path / "rhea_browser.html"
    save_rhea_browser_html(
        output_file,
        cache_dir=cache_dir,
        results_dir=None,
        n_features=64,
        use_linkml_store=False,
    )

    html_text = output_file.read_text()
    assert "EC-backed, GO-missing" in html_text
    assert (
        "The agent assigned classes to a reaction that already has EC support but is missing from the GO-labeled benchmark."
        in html_text
    )
