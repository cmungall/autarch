"""Tests for EC hierarchy report generation."""

from autarch.ec_hierarchy_report import (
    ClassifierReportEntry,
    ECNode,
    annotate_subtree_f1,
    format_ec_with_dashes,
    iter_display_prefixes,
    normalize_ec_for_ontology,
    render_ec_hierarchy_report,
)


def test_normalize_ec_for_ontology_strips_dashes_and_limits_depth() -> None:
    """EC identifiers should normalize to the sqlite:obo:ec ontology shape."""
    assert normalize_ec_for_ontology("EC:2.7.-.-") == "2.7"
    assert normalize_ec_for_ontology("1.1.1.1", max_depth=3) == "1.1.1"
    assert normalize_ec_for_ontology("7.-.-.-") == "7"


def test_iter_display_prefixes_returns_all_displayed_ancestors() -> None:
    """Displayed prefixes should include all EC levels through depth 3."""
    assert iter_display_prefixes("1.1.1.1") == ["1", "1.1", "1.1.1"]
    assert iter_display_prefixes("3.6.-.-") == ["3", "3.6"]


def test_format_ec_with_dashes_round_trips_partial_ec_ids() -> None:
    """Partial EC identifiers should render as four-level CURIEs."""
    assert format_ec_with_dashes("2.7") == "EC:2.7.-.-"
    assert format_ec_with_dashes("4.2.1") == "EC:4.2.1.-"


def test_render_ec_hierarchy_report_includes_go_counts_and_classifier_stats() -> None:
    """Rendered HTML should expose node counts plus classifier metrics."""
    nodes = [
        ECNode(
            ec_id="2.7",
            label="Transferring phosphorus-containing groups",
            depth=2,
            go_preview=[("GO:0016772", "transferase activity, transferring phosphorus-containing groups")],
            hidden_go_count=0,
            subsumed_go_count=9,
            subsumed_rhea_count=41,
            classes=[
                ClassifierReportEntry(
                    class_name="Kinase",
                    page_href="kinase.html",
                    go_id="GO:0016301",
                    go_label="kinase activity",
                    declared_ecs=["2.7.-.-"],
                    declared_ec_labels=["Transferring phosphorus-containing groups"],
                    f1_score=0.912,
                    precision=0.901,
                    recall=0.923,
                    mcc=0.881,
                    total=3444,
                    cyclomatic_complexity=21,
                    loc=164,
                    complexity_score=41.3,
                    declared_depth=2,
                )
            ],
        )
    ]

    annotate_subtree_f1(nodes)
    html = render_ec_hierarchy_report(nodes)

    assert "EC Hierarchy Report" in html
    assert "Subsumed GO 9" in html
    assert "Subsumed RHEAs 41" in html
    assert "GO:0016772" in html
    assert "transferase activity, transferring phosphorus-containing groups" in html
    assert "Kinase" in html
    assert "kinase activity" in html
    assert "tree-row" in html
    assert "detail-panel" in html
    assert ">0.912<" in html
    assert ">21<" in html
    assert ">41.3<" in html
    assert "Avg F1 0.912" in html


def test_annotate_subtree_f1_deduplicates_multi_placed_classifiers() -> None:
    """Subtree F1 averages should deduplicate classifiers by class name."""
    child = ECNode(
        ec_id="2.7.1",
        label="Transferases with an alcohol group as acceptor",
        depth=3,
        go_preview=[],
        hidden_go_count=0,
        subsumed_go_count=0,
        subsumed_rhea_count=0,
        classes=[
            ClassifierReportEntry(
                class_name="SharedKinase",
                page_href="sharedkinase.html",
                go_id=None,
                go_label=None,
                declared_ecs=["2.7.1.-"],
                declared_ec_labels=[""],
                f1_score=0.9,
                declared_depth=3,
            ),
            ClassifierReportEntry(
                class_name="AlcoholKinase",
                page_href="alcoholkinase.html",
                go_id=None,
                go_label=None,
                declared_ecs=["2.7.1.-"],
                declared_ec_labels=[""],
                f1_score=0.7,
                declared_depth=3,
            ),
        ],
    )
    root = ECNode(
        ec_id="2.7",
        label="Transferring phosphorus-containing groups",
        depth=2,
        go_preview=[],
        hidden_go_count=0,
        subsumed_go_count=0,
        subsumed_rhea_count=0,
        classes=[
            ClassifierReportEntry(
                class_name="SharedKinase",
                page_href="sharedkinase.html",
                go_id=None,
                go_label=None,
                declared_ecs=["2.7.-.-"],
                declared_ec_labels=[""],
                f1_score=0.9,
                declared_depth=2,
            )
        ],
        children=[child],
    )

    annotate_subtree_f1([root])

    assert root.subtree_f1_class_count == 2
    assert root.subtree_avg_f1 == 0.8
    assert child.subtree_f1_class_count == 2
    assert child.subtree_avg_f1 == 0.8


def test_render_ec_hierarchy_report_shows_broad_ec_placement_note_in_declared_ec_cell() -> None:
    """Broad GO-to-EC xrefs should be labeled as such in the Declared EC column."""
    nodes = [
        ECNode(
            ec_id="1",
            label="Oxidoreductases",
            depth=1,
            go_preview=[],
            hidden_go_count=0,
            subsumed_go_count=0,
            subsumed_rhea_count=0,
            classes=[
                ClassifierReportEntry(
                    class_name="Monooxygenase",
                    page_href="monooxygenase.html",
                    go_id="GO:0004497",
                    go_label="monooxygenase activity",
                    declared_ecs=["1.-.-.-"],
                    declared_ec_labels=["Oxidoreductases"],
                    placement_note="Placed via GO broad EC xref",
                    f1_score=0.856,
                    declared_depth=1,
                )
            ],
        )
    ]

    annotate_subtree_f1(nodes)
    html = render_ec_hierarchy_report(nodes)

    assert "monooxygenase.html" in html
    assert "EC:1.-.-.-" in html
    assert "Placed via GO broad EC xref" in html
