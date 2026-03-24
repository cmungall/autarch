"""Complexity analysis for reaction classifiers.

This module measures various complexity metrics for classifier code to help
identify over-engineered classifiers that may benefit from simplification.

Metrics collected:
- Cyclomatic Complexity (CC): Control flow paths through code (via radon)
- Lines of Code (LOC): Total lines in the classifier file
- ChEBI IDs: Count of hardcoded chemical identifiers
- If Statements: Number of conditional branches
- Exclusions: Number of negative classification returns
"""

import ast
import json
import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


REPO_ROOT = Path(__file__).resolve().parents[2]


@dataclass
class ClassifierComplexity:
    """Complexity metrics for a single classifier."""

    name: str
    file_path: str

    # Basic metrics
    loc: int = 0  # Lines of code in file
    sloc: int = 0  # Source lines (non-blank, non-comment)

    # AST-based metrics
    if_statements: int = 0
    patterns: int = 0  # Number of PATTERNS defined

    # Domain-specific metrics
    chebi_ids: int = 0  # Hardcoded chemical identifiers
    exclusions: int = 0  # is_member=False returns

    # Radon metrics
    cyclomatic_complexity: int = 0
    cc_rank: str = "?"  # A, B, C, D, E, F
    method_loc: int = 0  # LOC of check_membership_impl

    @property
    def complexity_score(self) -> float:
        """Composite complexity score (higher = more complex)."""
        # Weighted combination of metrics
        return (
            self.cyclomatic_complexity * 1.0 +
            self.chebi_ids * 0.5 +
            self.exclusions * 2.0 +
            self.if_statements * 0.3
        )

    @property
    def is_disabled(self) -> bool:
        """Check if classifier is effectively disabled (CC <= 2)."""
        return self.cyclomatic_complexity <= 2

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "name": self.name,
            "file_path": self.file_path,
            "loc": self.loc,
            "sloc": self.sloc,
            "if_statements": self.if_statements,
            "patterns": self.patterns,
            "chebi_ids": self.chebi_ids,
            "exclusions": self.exclusions,
            "cyclomatic_complexity": self.cyclomatic_complexity,
            "cc_rank": self.cc_rank,
            "method_loc": self.method_loc,
            "complexity_score": round(self.complexity_score, 2),
            "is_disabled": self.is_disabled,
        }


@dataclass
class ComplexityReport:
    """Aggregate complexity report for all classifiers."""

    classifiers: list[ClassifierComplexity] = field(default_factory=list)

    @property
    def total_loc(self) -> int:
        return sum(c.loc for c in self.classifiers)

    @property
    def avg_cc(self) -> float:
        if not self.classifiers:
            return 0.0
        return sum(c.cyclomatic_complexity for c in self.classifiers) / len(self.classifiers)

    @property
    def avg_chebi_ids(self) -> float:
        if not self.classifiers:
            return 0.0
        return sum(c.chebi_ids for c in self.classifiers) / len(self.classifiers)

    @property
    def most_complex(self) -> list[ClassifierComplexity]:
        """Top 10 most complex classifiers by CC."""
        return sorted(self.classifiers, key=lambda x: x.cyclomatic_complexity, reverse=True)[:10]

    @property
    def simplest(self) -> list[ClassifierComplexity]:
        """Top 5 simplest classifiers (excluding disabled)."""
        active = [c for c in self.classifiers if not c.is_disabled]
        return sorted(active, key=lambda x: x.cyclomatic_complexity)[:5]

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "summary": {
                "total_classifiers": len(self.classifiers),
                "total_loc": self.total_loc,
                "avg_cyclomatic_complexity": round(self.avg_cc, 2),
                "avg_chebi_ids": round(self.avg_chebi_ids, 2),
                "disabled_count": sum(1 for c in self.classifiers if c.is_disabled),
            },
            "classifiers": [c.to_dict() for c in self.classifiers],
        }


def analyze_file_ast(filepath: Path) -> dict:
    """Analyze a Python file using AST.

    Returns dict with: loc, if_statements, patterns
    """
    with open(filepath) as f:
        content = f.read()
        tree = ast.parse(content)

    metrics = {
        "loc": len(content.splitlines()),
        "sloc": len([line for line in content.splitlines() if line.strip() and not line.strip().startswith("#")]),
        "if_statements": 0,
        "patterns": 0,
    }

    for node in ast.walk(tree):
        # Count if statements
        if isinstance(node, ast.If):
            metrics["if_statements"] += 1

        # Count PATTERNS list items
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == "PATTERNS":
                    if isinstance(node.value, ast.List):
                        metrics["patterns"] = len(node.value.elts)

    return metrics


def count_chebi_ids(filepath: Path) -> int:
    """Count unique ChEBI IDs mentioned in a file."""
    with open(filepath) as f:
        content = f.read()

    # Match CHEBI:12345, "CHEBI:12345", CHEBI_12345, etc.
    chebi_pattern = r'CHEBI[:\d_]+|"CHEBI:\d+"'
    matches = re.findall(chebi_pattern, content)
    return len(set(matches))


def count_exclusions(filepath: Path) -> int:
    """Count is_member=False returns (exclusion rules)."""
    with open(filepath) as f:
        content = f.read()
    return content.count("is_member=False")


def get_radon_complexity(ontology_dir: Path) -> dict[str, dict]:
    """Get cyclomatic complexity from radon for all classifiers.

    Returns dict mapping filename -> {complexity, rank, lineno, endline}
    """
    result = subprocess.run(
        ["uv", "run", "radon", "cc", str(ontology_dir), "-s", "--json"],
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        return {}

    data = json.loads(result.stdout)

    # Extract complexity for check_membership_impl methods
    complexities = {}
    for filepath, items in data.items():
        filename = Path(filepath).stem
        for item in items:
            if item["name"] == "check_membership_impl":
                complexities[filename] = {
                    "complexity": item["complexity"],
                    "rank": item["rank"],
                    "lineno": item["lineno"],
                    "endline": item["endline"],
                    "method_loc": item["endline"] - item["lineno"] + 1,
                }
                break

    return complexities


def analyze_classifier(filepath: Path, radon_data: dict) -> ClassifierComplexity:
    """Analyze a single classifier file."""
    name = filepath.stem
    resolved_path = filepath.resolve()
    try:
        display_path = str(resolved_path.relative_to(REPO_ROOT))
    except ValueError:
        display_path = str(resolved_path)

    # Get AST metrics
    ast_metrics = analyze_file_ast(filepath)

    # Get domain-specific metrics
    chebi_count = count_chebi_ids(filepath)
    exclusion_count = count_exclusions(filepath)

    # Get radon metrics
    radon = radon_data.get(name, {})

    return ClassifierComplexity(
        name=name,
        file_path=display_path,
        loc=ast_metrics["loc"],
        sloc=ast_metrics["sloc"],
        if_statements=ast_metrics["if_statements"],
        patterns=ast_metrics["patterns"],
        chebi_ids=chebi_count,
        exclusions=exclusion_count,
        cyclomatic_complexity=radon.get("complexity", 0),
        cc_rank=radon.get("rank", "?"),
        method_loc=radon.get("method_loc", 0),
    )


def analyze_all_classifiers(ontology_dir: Optional[Path] = None) -> ComplexityReport:
    """Analyze all classifier files and return a complexity report.

    Args:
        ontology_dir: Path to ontology directory. Defaults to src/autarch/ontology.

    Returns:
        ComplexityReport with metrics for all classifiers.

    Example:
        >>> report = analyze_all_classifiers()
        >>> isinstance(report, ComplexityReport)
        True
    """
    if ontology_dir is None:
        # Find the ontology directory relative to this file
        this_file = Path(__file__)
        ontology_dir = this_file.parent / "ontology"

    # Get radon complexity for all files
    radon_data = get_radon_complexity(ontology_dir)

    # Analyze each classifier file
    classifiers = []
    for filepath in sorted(ontology_dir.glob("*.py")):
        # Skip __init__.py and base classes
        if filepath.name.startswith("_") or filepath.name == "reaction.py":
            continue

        complexity = analyze_classifier(filepath, radon_data)
        classifiers.append(complexity)

    return ComplexityReport(classifiers=classifiers)


def print_complexity_report(report: ComplexityReport, verbose: bool = False) -> None:
    """Print a formatted complexity report to stdout."""
    print("=" * 75)
    print("CLASSIFIER COMPLEXITY ANALYSIS")
    print("=" * 75)
    print()

    # Summary
    print("SUMMARY")
    print("-" * 40)
    print(f"Total classifiers: {len(report.classifiers)}")
    print(f"Total LOC: {report.total_loc}")
    print(f"Average Cyclomatic Complexity: {report.avg_cc:.1f}")
    print(f"Average ChEBI IDs per classifier: {report.avg_chebi_ids:.1f}")
    print(f"Disabled classifiers: {sum(1 for c in report.classifiers if c.is_disabled)}")
    print()

    # Most complex
    print("MOST COMPLEX (Top 10 by Cyclomatic Complexity)")
    print("-" * 75)
    print(f"{'Classifier':<35} {'CC':>4} {'Rank':>5} {'ChEBI':>6} {'Excl':>5} {'LOC':>5}")
    print("-" * 75)
    for c in report.most_complex:
        print(f"{c.name:<35} {c.cyclomatic_complexity:>4} {c.cc_rank:>5} {c.chebi_ids:>6} {c.exclusions:>5} {c.method_loc:>5}")
    print()

    # Simplest (good examples)
    print("SIMPLEST ACTIVE CLASSIFIERS (Good Examples)")
    print("-" * 75)
    for c in report.simplest:
        print(f"{c.name:<35} CC={c.cyclomatic_complexity}, {c.chebi_ids} ChEBI IDs")
    print()

    # Complexity ranks explanation
    print("COMPLEXITY RANKS")
    print("-" * 40)
    print("A (1-5): Simple, easy to test")
    print("B (6-10): Moderate complexity")
    print("C (11-20): Complex, consider refactoring")
    print("D (21-30): Very complex")
    print("E (31-40): Highly complex, hard to maintain")
    print("F (41+): Extremely complex, needs refactoring")
    print()

    if verbose:
        print("ALL CLASSIFIERS")
        print("-" * 75)
        print(f"{'Classifier':<35} {'CC':>4} {'Rank':>5} {'ChEBI':>6} {'If':>4} {'Excl':>5} {'LOC':>5}")
        print("-" * 75)
        for c in sorted(report.classifiers, key=lambda x: x.name):
            print(f"{c.name:<35} {c.cyclomatic_complexity:>4} {c.cc_rank:>5} {c.chebi_ids:>6} {c.if_statements:>4} {c.exclusions:>5} {c.loc:>5}")


def save_complexity_report(report: ComplexityReport, output_dir: Path) -> None:
    """Save complexity report as JSON and CSV files."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save as JSON
    json_path = output_dir / "complexity_report.json"
    with open(json_path, "w") as f:
        json.dump(report.to_dict(), f, indent=2)
    print(f"Saved JSON report to {json_path}")

    # Save as CSV
    csv_path = output_dir / "complexity_report.csv"
    with open(csv_path, "w") as f:
        headers = ["name", "cyclomatic_complexity", "cc_rank", "chebi_ids",
                   "exclusions", "if_statements", "method_loc", "loc", "complexity_score"]
        f.write(",".join(headers) + "\n")
        for c in sorted(report.classifiers, key=lambda x: x.cyclomatic_complexity, reverse=True):
            row = [
                c.name,
                str(c.cyclomatic_complexity),
                c.cc_rank,
                str(c.chebi_ids),
                str(c.exclusions),
                str(c.if_statements),
                str(c.method_loc),
                str(c.loc),
                str(round(c.complexity_score, 2)),
            ]
            f.write(",".join(row) + "\n")
    print(f"Saved CSV report to {csv_path}")


if __name__ == "__main__":
    report = analyze_all_classifiers()
    print_complexity_report(report, verbose=True)
