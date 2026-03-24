"""Tests for ModelSEED benchmark materialization."""

import json
from pathlib import Path

from autarch.datamodel import ClassificationResult
from autarch.etl.modelseed_etl import ModelSeedETL, cache_modelseed_dataset
from autarch.modelseed_benchmark import build_modelseed_benchmark


COMPOUND_HEADER = [
    "id",
    "abbreviation",
    "name",
    "formula",
    "mass",
    "source",
    "inchikey",
    "charge",
    "is_core",
    "is_obsolete",
    "linked_compound",
    "is_cofactor",
    "deltag",
    "deltagerr",
    "pka",
    "pkb",
    "abstract_compound",
    "comprised_of",
    "aliases",
    "smiles",
    "notes",
]

REACTION_HEADER = [
    "id",
    "abbreviation",
    "name",
    "code",
    "stoichiometry",
    "is_transport",
    "equation",
    "definition",
    "reversibility",
    "direction",
    "abstract_reaction",
    "pathways",
    "aliases",
    "ec_numbers",
    "deltag",
    "deltagerr",
    "compound_ids",
    "status",
    "is_obsolete",
    "linked_reaction",
    "notes",
    "source",
]


def _write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    lines = ["\t".join(header)]
    lines.extend("\t".join(row) for row in rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def _write_chebi_lookup_cache(cache_dir: Path) -> None:
    (cache_dir / "smiles_to_chebi.json").write_text(
        json.dumps(
            {
                "O": "CHEBI:15377",
                "CC(N)C(=O)O": "CHEBI:16449",
                "N": "CHEBI:16134",
            }
        )
    )
    (cache_dir / "chebi_names.json").write_text(
        json.dumps(
            {
                "CHEBI:15377": "water",
                "CHEBI:16449": "alanine",
                "CHEBI:16134": "ammonia",
            }
        )
    )


def _write_modelseed_benchmark_test_files(base_dir: Path) -> None:
    _write_tsv(
        base_dir / "compounds.tsv",
        COMPOUND_HEADER,
        [
            [
                "cpd00001",
                "h2o",
                "water",
                "H2O",
                "18.0",
                "Primary Database",
                "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
                "0",
                "1",
                "0",
                "null",
                "0",
                "0",
                "0",
                "null",
                "null",
                "null",
                "null",
                "Name: water|KEGG: C00001",
                "O",
                "",
            ],
            [
                "cpd00002",
                "ala",
                "alanine",
                "C3H7NO2",
                "89.0",
                "Primary Database",
                "QNAYBMKLOCPYGJ-REOHCLBHSA-N",
                "0",
                "1",
                "0",
                "null",
                "0",
                "0",
                "0",
                "null",
                "null",
                "null",
                "null",
                "Name: alanine",
                "C[C@H](N)C(=O)O",
                "",
            ],
            [
                "cpd00004",
                "nh3",
                "ammonia",
                "NH3",
                "17.0",
                "Primary Database",
                "",
                "0",
                "1",
                "0",
                "null",
                "0",
                "0",
                "0",
                "null",
                "null",
                "null",
                "null",
                "Name: ammonia",
                "N",
                "",
            ],
            [
                "cpd99999",
                "mystery",
                "mystery compound",
                "X",
                "1.0",
                "Primary Database",
                "",
                "0",
                "1",
                "0",
                "null",
                "0",
                "0",
                "0",
                "null",
                "null",
                "null",
                "null",
                "Name: mystery compound",
                "",
                "",
            ],
        ],
    )

    _write_tsv(
        base_dir / "reactions.tsv",
        REACTION_HEADER,
        [
            [
                "rxn00001",
                "watert",
                "water transport",
                "",
                '-1:cpd00001:0:0:"Water";1:cpd00001:1:0:"Water"',
                "1",
                "",
                "water transport",
                "=",
                "=",
                "null",
                "null",
                "",
                "",
                "0",
                "0",
                "cpd00001",
                "OK",
                "0",
                "null",
                "",
                "Primary Database",
            ],
            [
                "rxn00002",
                "alanine-ok",
                "alanine benchmark reaction",
                "",
                '-1:cpd00002:0:0:"alanine";1:cpd00004:0:0:"ammonia"',
                "0",
                "",
                "alanine benchmark reaction",
                "=",
                "=",
                "null",
                "null",
                "",
                "2.1.1.1",
                "0",
                "0",
                "cpd00002;cpd00004",
                "OK",
                "0",
                "null",
                "",
                "Primary Database",
            ],
            [
                "rxn00003",
                "bad-status",
                "bad status reaction",
                "",
                '-1:cpd00002:0:0:"alanine";1:cpd00004:0:0:"ammonia"',
                "0",
                "",
                "bad status reaction",
                "=",
                "=",
                "null",
                "null",
                "",
                "",
                "0",
                "0",
                "cpd00002;cpd00004",
                "CPDFORMERROR",
                "0",
                "null",
                "",
                "Primary Database",
            ],
            [
                "rxn00004",
                "partial-map",
                "partial mapping reaction",
                "",
                '-1:cpd99999:0:0:"mystery compound";1:cpd00004:0:0:"ammonia"',
                "0",
                "",
                "partial mapping reaction",
                "=",
                "=",
                "null",
                "null",
                "",
                "",
                "0",
                "0",
                "cpd99999;cpd00004",
                "OK",
                "0",
                "null",
                "",
                "Primary Database",
            ],
            [
                "rxn00005",
                "fractional",
                "fractional reaction",
                "",
                '-0.5:cpd00001:0:0:"Water";1:cpd00004:0:0:"ammonia"',
                "0",
                "",
                "fractional reaction",
                "=",
                "=",
                "null",
                "null",
                "",
                "1.1.1.1",
                "0",
                "0",
                "cpd00001;cpd00004",
                "OK",
                "0",
                "null",
                "",
                "Primary Database",
            ],
        ],
    )

    _write_tsv(
        base_dir / "Aliases" / "Unique_ModelSEED_Reaction_Aliases.txt",
        ["ModelSEED ID", "External ID", "Source"],
        [["rxn00001", "10000", "rhea"]],
    )

    _write_tsv(
        base_dir / "Aliases" / "Unique_ModelSEED_Reaction_ECs.txt",
        ["ModelSEED ID", "External ID", "Source"],
        [["rxn00002", "2.1.1.1", "Enzyme Class"], ["rxn00005", "1.1.1.1", "Enzyme Class"]],
    )


class _AlwaysPositive:
    def check_membership(self, reaction):
        names = {participant.name for participant in reaction.all_participants()}
        return ClassificationResult(
            is_member="alanine" in names,
            explanation="contains alanine",
        )


class _NeverPositive:
    def check_membership(self, reaction):
        return ClassificationResult(is_member=False, explanation="never")


class _FakeReactionClassifier:
    def __init__(self):
        self.reaction_classes = {
            "AlwaysPositive": _AlwaysPositive,
            "NeverPositive": _NeverPositive,
        }


def test_build_modelseed_benchmark_filters_and_writes_artifacts(
    tmp_path: Path, monkeypatch
):
    """Benchmark builder should materialize the default filtered OOD slice."""
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    _write_chebi_lookup_cache(cache_dir)
    _write_modelseed_benchmark_test_files(cache_dir / "modelseed")

    monkeypatch.setattr(ModelSeedETL, "download_required_files", lambda self: None)
    monkeypatch.setattr(
        "autarch.modelseed_benchmark.ReactionClassifier",
        _FakeReactionClassifier,
    )

    cache_modelseed_dataset(cache_dir=cache_dir)
    summary = build_modelseed_benchmark(
        cache_dir=cache_dir,
        output_dir=cache_dir / "benchmark",
    )

    assert summary["input_reactions"] == 5
    assert summary["selected_reactions"] == 1
    assert summary["selected_with_ec"] == 1
    assert summary["positive_reactions"] == 1
    assert summary["skipped_counts"]["rhea_alias"] == 1
    assert summary["skipped_counts"]["status"] == 1
    assert summary["skipped_counts"]["partial_mapping"] == 1
    assert summary["skipped_counts"]["fractional_stoichiometry"] == 1
    assert summary["top_positive_classes"] == [("AlwaysPositive", 1)]

    candidates_lines = (cache_dir / "benchmark" / "candidates.jsonl").read_text().splitlines()
    predictions_lines = (cache_dir / "benchmark" / "predictions.jsonl").read_text().splitlines()
    benchmark_summary = json.loads((cache_dir / "benchmark" / "summary.json").read_text())

    assert len(candidates_lines) == 1
    assert len(predictions_lines) == 1
    assert json.loads(candidates_lines[0])["modelseed_id"] == "rxn00002"
    assert json.loads(predictions_lines[0])["positive_classes"] == ["AlwaysPositive"]
    assert benchmark_summary["artifacts"]["output_dir"] == str(cache_dir / "benchmark")
