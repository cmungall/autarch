"""Architecture tests for explicit aggregate wrappers."""

from pathlib import Path


def test_public_ontology_classes_do_not_use_ec_prefix_proxy_dispatch() -> None:
    """Only the helper module may define EC-prefix aggregate dispatch helpers.

    Public classifier modules should use explicit child disjunctions instead of
    inferred EC-prefix membership.
    """
    ontology_dir = Path("src/autarch/ontology")
    offenders: list[str] = []
    for path in sorted(ontology_dir.glob("*.py")):
        if path.name == "ec_prefix_aggregate.py":
            continue
        text = path.read_text()
        if "aggregate_ec_prefix_membership(" in text:
            offenders.append(f"{path}:aggregate_ec_prefix_membership")
        if "aggregate_ec_prefix_supports_evaluation(" in text:
            offenders.append(f"{path}:aggregate_ec_prefix_supports_evaluation")

    assert offenders == []
