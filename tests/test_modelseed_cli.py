"""CLI tests for the ModelSEED cache command."""

import json
from pathlib import Path
from typing import TypedDict

from typer.testing import CliRunner

from autarch.cli import app


runner = CliRunner()


class CacheModelseedCalls(TypedDict):
    """Captured arguments passed to the cache helper."""

    cache_dir: str
    limit: int | None
    force_download: bool


def test_cache_modelseed_command(monkeypatch, tmp_path: Path):
    """ModelSEED cache command delegates to the ETL helper and reports success."""
    calls: CacheModelseedCalls = {
        "cache_dir": "",
        "limit": None,
        "force_download": False,
    }

    def fake_cache_modelseed_dataset(cache_dir: str, limit: int | None = None, force_download: bool = False):
        calls["cache_dir"] = cache_dir
        calls["limit"] = limit
        calls["force_download"] = force_download
        path = Path(cache_dir)
        path.mkdir(parents=True, exist_ok=True)
        (path / "modelseed_compounds.jsonl").write_text("")
        (path / "modelseed_to_chebi.jsonl").write_text("")
        (path / "modelseed_reactions.jsonl").write_text("")
        return {
            "total_compounds": 4,
            "mapped_compounds": 3,
            "total_reactions": 2,
            "reactions_with_rhea": 1,
            "fully_mapped_reactions": 1,
        }

    monkeypatch.setattr(
        "autarch.etl.modelseed_etl.cache_modelseed_dataset",
        fake_cache_modelseed_dataset,
    )

    result = runner.invoke(
        app,
        ["cache-modelseed", "--cache-dir", str(tmp_path), "--limit", "25", "--force"],
    )

    assert result.exit_code == 0
    assert "Cached 4 ModelSEED compounds" in result.stdout
    assert "Mapped to CHEBI: 3" in result.stdout
    assert "Reactions with RHEA aliases: 1" in result.stdout
    assert calls["cache_dir"] == str(tmp_path)
    assert calls["limit"] == 25
    assert calls["force_download"] is True


def test_summarize_modelseed_command_json(monkeypatch, tmp_path: Path):
    """Summary command returns cached ModelSEED coverage as JSON."""

    def fake_summarize_modelseed_cache(cache_dir: str):
        assert cache_dir == str(tmp_path)
        return {
            "total_compounds": 4,
            "mapped_compounds": 3,
            "total_reactions": 2,
            "reactions_with_rhea": 1,
            "cached_reactions": 2,
            "cached_reactions_with_rhea": 1,
            "cached_fully_mapped_reactions": 1,
            "compound_mapping_methods": {"exact_smiles": 2, "unmapped": 1},
            "reaction_alias_sources": {
                "rhea": {"rows": 2, "distinct_ids": 1},
                "KEGG": {"rows": 1, "distinct_ids": 1},
            },
            "compound_alias_sources": {
                "KEGG": {"rows": 3, "distinct_ids": 2},
            },
            "ec_coverage": {"rows": 5, "distinct_reactions": 2},
        }

    monkeypatch.setattr(
        "autarch.etl.modelseed_etl.summarize_modelseed_cache",
        fake_summarize_modelseed_cache,
    )

    result = runner.invoke(
        app,
        ["summarize-modelseed", "--cache-dir", str(tmp_path), "--format", "json"],
    )

    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert payload["total_compounds"] == 4
    assert payload["reaction_alias_sources"]["rhea"]["rows"] == 2


def test_benchmark_modelseed_command_json(monkeypatch, tmp_path: Path):
    """Benchmark command returns cached slice metadata as JSON."""

    def fake_build_modelseed_benchmark(
        cache_dir: str,
        output_dir: str | None = None,
        limit: int | None = None,
        include_rhea: bool = False,
        include_transport: bool = False,
        all_statuses: bool = False,
        allow_partial_mapping: bool = False,
        require_ec: bool = False,
    ):
        assert cache_dir == str(tmp_path)
        assert output_dir is None
        assert include_rhea is False
        assert include_transport is False
        assert all_statuses is False
        assert allow_partial_mapping is False
        assert require_ec is True
        assert limit == 10
        return {
            "input_reactions": 10,
            "selected_reactions": 4,
            "selected_with_ec": 4,
            "selected_without_ec": 0,
            "positive_reactions": 3,
            "positive_rate": 0.75,
            "classifiers_run": 2,
            "duration_seconds": 0.1,
            "skipped_counts": {"rhea_alias": 3},
            "top_positive_classes": [["Hydrolase", 2]],
            "artifacts": {"output_dir": str(tmp_path / "modelseed_benchmark")},
        }

    monkeypatch.setattr(
        "autarch.modelseed_benchmark.build_modelseed_benchmark",
        fake_build_modelseed_benchmark,
    )

    result = runner.invoke(
        app,
        [
            "benchmark-modelseed",
            "--cache-dir",
            str(tmp_path),
            "--require-ec",
            "--limit",
            "10",
            "--format",
            "json",
        ],
    )

    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert payload["selected_reactions"] == 4
    assert payload["positive_reactions"] == 3
