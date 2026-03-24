"""Tests for GO-missing and unannotated positives in HTML export."""

from types import SimpleNamespace
from typing import Any, cast

import pandas as pd

from autarch.export import build_go_missing_positive_index, generate_class_page
from autarch.ontology.cis_trans_isomerase import CisTransIsomerase


class FakeClassifier:
    """Minimal classifier stub for export tests."""

    def __init__(self) -> None:
        self.reaction_classes: dict[str, type] = {"Hydrolase": FakeHydrolase}

    def classify(self, reaction):  # type: ignore[no-untyped-def]
        label = getattr(reaction, "label", "")
        is_match = "ATP + H2O" in label
        return {
            "Hydrolase": SimpleNamespace(
                is_member=is_match,
                explanation="matched hydrolysis signature" if is_match else "",
            )
        }


class FakeHydrolase:
    """Minimal reaction class metadata for export tests."""

    GO_ID = "GO:0016787"
    EC_NUMBER_PREFIX = "3.-.-.-"
    PATTERNS = ["hydrolysis signature"]

    def check_membership_impl(self, reaction) -> None:  # type: ignore[no-untyped-def]
        return None


def test_build_go_missing_positive_index_collects_positive_matches() -> None:
    """Positive matches on GO-missing reactions should be indexed per class."""
    classifier = FakeClassifier()
    rhea_reactions = {
        "RHEA:10000": {
            "label": "ATP + H2O = ADP + phosphate",
            "go_terms": [],
            "ec_numbers": ["3.6.1.3"],
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
        },
        "RHEA:20000": {
            "label": "annotated reaction",
            "go_terms": ["GO:0000001"],
            "ec_numbers": ["3.1.1.1"],
            "reaction": {
                "left_participants": [],
                "right_participants": [],
            },
        },
    }

    index = build_go_missing_positive_index(rhea_reactions, cast(Any, classifier))
    assert "Hydrolase" in index
    assert index["Hydrolase"] == [
        {
            "rhea_id": "RHEA:10000",
            "label": "ATP + H2O = ADP + phosphate",
            "explanation": "matched hydrolysis signature",
            "ec_numbers": ["3.6.1.3"],
            "annotation_bucket": "ec_aligned",
        }
    ]


def test_generate_class_page_splits_go_missing_positives_by_ec_provenance() -> None:
    """Class pages should distinguish unannotated and EC-backed GO-missing positives."""
    classifier = FakeClassifier()
    html = generate_class_page(
        "Hydrolase",
        {
            "total": 12,
            "f1_score": 0.9,
            "precision": 0.8,
            "recall": 1.0,
            "accuracy": 0.95,
            "mcc": 0.88,
            "specificity": 0.9,
            "tp": 3,
            "fp": 1,
            "fn": 0,
            "tn": 8,
        },
        pd.DataFrame(),
        {},
        cast(Any, classifier),
        {
            "Hydrolase": [
                {
                    "rhea_id": "RHEA:10000",
                    "label": "ATP + H2O = ADP + phosphate",
                    "explanation": "matched hydrolysis signature",
                    "ec_numbers": ["3.6.1.3"],
                    "annotation_bucket": "ec_aligned",
                },
                {
                    "rhea_id": "RHEA:10001",
                    "label": "farnesyl diphosphate + H2O = terpene alcohol + diphosphate",
                    "explanation": "matched hydrolysis signature",
                    "ec_numbers": ["4.2.3.166"],
                    "annotation_bucket": "ec_conflict",
                },
                {
                    "rhea_id": "RHEA:10002",
                    "label": "unannotated hydrolysis candidate",
                    "explanation": "matched hydrolysis signature",
                    "ec_numbers": [],
                    "annotation_bucket": "unannotated",
                }
            ]
        },
    )

    assert "Positives Without Direct GO Support" in html
    assert "Truly unannotated candidate positives" in html
    assert "EC-backed, GO-missing positives supporting this class" in html
    assert "EC-backed, GO-missing positives conflicting with this class" in html
    assert "Class-level declarative context used by the classifier" in html
    assert "PATTERNS" in html
    assert "RHEA:10000" in html
    assert "RHEA:10001" in html
    assert "RHEA:10002" in html
    assert "ATP + H2O = ADP + phosphate" in html
    assert "EC:3.6.1.3" in html
    assert "EC:4.2.3.166" in html
    assert "matched hydrolysis signature" in html


def test_generate_class_page_shows_broad_ec_xrefs_without_exact_prefix() -> None:
    """Broad EC xrefs should be rendered separately from exact EC prefixes."""

    class FakeMonooxygenase:
        GO_ID = "GO:0004497"
        EC_NUMBER_PREFIX = None
        EC_BROAD_XREFS = ["1.-.-.-"]

        def check_membership_impl(self, reaction) -> None:  # type: ignore[no-untyped-def]
            return None

    classifier = FakeClassifier()
    classifier.reaction_classes["Monooxygenase"] = FakeMonooxygenase
    html = generate_class_page(
        "Monooxygenase",
        {
            "total": 12,
            "f1_score": 0.9,
            "precision": 0.8,
            "recall": 1.0,
            "accuracy": 0.95,
            "mcc": 0.88,
            "specificity": 0.9,
            "tp": 3,
            "fp": 1,
            "fn": 0,
            "tn": 8,
        },
        pd.DataFrame(),
        {},
        cast(Any, classifier),
        {},
    )

    assert "GO:0004497" in html
    assert "EC Prefix" not in html
    assert "GO EC Xref Type" in html
    assert "broad match" in html
    assert "GO EC Broad Xrefs" in html
    assert "EC:1.-.-.-" in html


def test_generate_class_page_shows_declared_ec_scope_and_go_broad_xrefs() -> None:
    """Exact class EC scope and GO broad xrefs should both be visible when present."""

    class FakeTerpeneSynthase:
        GO_ID = "GO:0010333"
        EC_NUMBER_PREFIX = "4.2.3.-"
        EC_BROAD_XREFS = ["4.2.3.-"]

        def check_membership_impl(self, reaction) -> None:  # type: ignore[no-untyped-def]
            return None

    classifier = FakeClassifier()
    classifier.reaction_classes["TerpeneSynthase"] = FakeTerpeneSynthase
    html = generate_class_page(
        "TerpeneSynthase",
        {
            "total": 12,
            "f1_score": 0.9,
            "precision": 0.8,
            "recall": 1.0,
            "accuracy": 0.95,
            "mcc": 0.88,
            "specificity": 0.9,
            "tp": 3,
            "fp": 1,
            "fn": 0,
            "tn": 8,
        },
        pd.DataFrame(),
        {},
        cast(Any, classifier),
        {},
    )

    assert "Declared EC Scope" in html
    assert "GO EC Xref Type" in html
    assert "GO EC Broad Xrefs" in html
    assert "EC:4.2.3.-" in html


def test_generate_class_page_explains_aggregate_ec_wrapper_strategy() -> None:
    """Aggregate wrappers should expose child-dispatch semantics."""

    classifier = FakeClassifier()
    classifier.reaction_classes["CisTransIsomerase"] = CisTransIsomerase
    html = generate_class_page(
        "CisTransIsomerase",
        {
            "total": 12,
            "f1_score": 0.9,
            "precision": 0.8,
            "recall": 1.0,
            "accuracy": 0.95,
            "mcc": 0.88,
            "specificity": 0.9,
            "tp": 3,
            "fp": 1,
            "fn": 0,
            "tn": 8,
        },
        pd.DataFrame(),
        {},
        cast(Any, classifier),
        {},
    )

    assert "Aggregate Membership Strategy" in html
    assert "does not inspect reaction GO or EC annotations" in html
    assert "Aggregate Children" in html
    assert "EC:5.2.-.-" in html
    assert "Curated child classifiers" in html
    assert "CisTransIsomerases" in html
    assert "PeptidylProlylCisTransIsomerase" in html
    assert "CONCEPT_PHRASE" not in html
