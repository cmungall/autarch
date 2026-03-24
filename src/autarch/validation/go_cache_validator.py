"""Validation for cached GO ancestor closure consistency."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Optional

from oaklib import get_adapter  # type: ignore[import-untyped]
from oaklib.datamodels.vocabulary import IS_A  # type: ignore[import-untyped]


def _dedupe_preserve_order(values: list[str]) -> list[str]:
    """Deduplicate a list while preserving first-seen order."""
    seen: set[str] = set()
    result: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        result.append(value)
    return result


@dataclass
class GoCacheIssue:
    """A single cached GO term inconsistency."""

    go_id: str
    label: str
    duplicate_ancestors: list[str]
    missing_ancestors: list[str]
    extra_ancestors: list[str]


@dataclass
class GoCacheValidationReport:
    """Summary of GO cache validation."""

    total_terms: int
    issue_count: int
    duplicate_issue_count: int
    mismatch_issue_count: int
    issues: list[GoCacheIssue]

    @property
    def is_valid(self) -> bool:
        return self.issue_count == 0


def validate_go_cache(
    cache_path: str | Path = "cache/go_terms.jsonl",
    *,
    max_issues: Optional[int] = 25,
) -> GoCacheValidationReport:
    """Validate cached GO ancestor closures against the local GO sqlite."""

    path = Path(cache_path)
    if not path.exists():
        raise FileNotFoundError(f"GO cache not found: {path}")

    adapter = get_adapter("sqlite:obo:go")
    issues: list[GoCacheIssue] = []
    total_terms = 0
    issue_count = 0
    duplicate_issue_count = 0
    mismatch_issue_count = 0

    with path.open() as handle:
        for line in handle:
            row = json.loads(line)
            total_terms += 1
            go_id = row["go_id"]
            label = row.get("label", "")
            raw_cached = [a for a in row.get("ancestors", []) if isinstance(a, str) and a.startswith("GO:")]
            cached = _dedupe_preserve_order(raw_cached)
            duplicates = sorted({a for a in raw_cached if raw_cached.count(a) > 1})
            fresh = _dedupe_preserve_order(
                [a for a in adapter.ancestors(go_id, predicates=[IS_A]) if isinstance(a, str) and a.startswith("GO:")]
            )

            missing = sorted(set(fresh) - set(cached))
            extra = sorted(set(cached) - set(fresh))

            if duplicates:
                duplicate_issue_count += 1
            if missing or extra:
                mismatch_issue_count += 1

            if duplicates or missing or extra:
                issue_count += 1
                if max_issues is None or len(issues) < max_issues:
                    issues.append(
                        GoCacheIssue(
                            go_id=go_id,
                            label=label,
                            duplicate_ancestors=duplicates,
                            missing_ancestors=missing,
                            extra_ancestors=extra,
                        )
                    )

    return GoCacheValidationReport(
        total_terms=total_terms,
        issue_count=issue_count,
        duplicate_issue_count=duplicate_issue_count,
        mismatch_issue_count=mismatch_issue_count,
        issues=issues,
    )


def print_go_cache_validation_report(report: GoCacheValidationReport) -> bool:
    """Print a concise validation report and return success."""

    print("🔍 Validating cached GO ancestor closure against local GO sqlite...\n")
    print(f"Checked terms: {report.total_terms}")
    print(f"Terms with duplicate cached ancestors: {report.duplicate_issue_count}")
    print(f"Terms with ancestor-set mismatches: {report.mismatch_issue_count}")

    if report.duplicate_issue_count == 0 and report.mismatch_issue_count == 0:
        print("\n✅ GO cache is consistent with local GO sqlite.")
        return True

    print("\n❌ GO cache drift detected.")
    print("   Regenerate with: just cache-go-all")

    for issue in report.issues:
        print(f"\n{issue.go_id} {issue.label}")
        if issue.duplicate_ancestors:
            print(f"  duplicate cached ancestors: {', '.join(issue.duplicate_ancestors[:10])}")
        if issue.missing_ancestors:
            print(f"  missing from cache: {', '.join(issue.missing_ancestors[:10])}")
        if issue.extra_ancestors:
            print(f"  extra in cache: {', '.join(issue.extra_ancestors[:10])}")

    return False
