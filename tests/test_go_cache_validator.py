"""Tests for GO cache validation."""

from __future__ import annotations

from pathlib import Path

from autarch.validation.go_cache_validator import validate_go_cache


def test_validate_go_cache_detects_duplicate_and_missing_ancestors(
    tmp_path: Path, monkeypatch
) -> None:
    cache_path = tmp_path / "go_terms.jsonl"
    cache_path.write_text(
        '{"go_id":"GO:TEST","label":"test term","ancestors":["GO:ROOT","GO:ROOT"]}\n'
    )

    class FakeAdapter:
        def ancestors(self, go_id: str, predicates: list[str]):  # noqa: ARG002
            return ["GO:TEST", "GO:ROOT", "GO:PARENT"]

    monkeypatch.setattr(
        "autarch.validation.go_cache_validator.get_adapter",
        lambda spec: FakeAdapter(),
    )

    report = validate_go_cache(cache_path)

    assert report.total_terms == 1
    assert report.duplicate_issue_count == 1
    assert report.mismatch_issue_count == 1
    assert len(report.issues) == 1
    assert report.issues[0].duplicate_ancestors == ["GO:ROOT"]
    assert report.issues[0].missing_ancestors == ["GO:PARENT", "GO:TEST"]
