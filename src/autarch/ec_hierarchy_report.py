"""One-page EC hierarchy report with GO, RHEA, and classifier overlays."""

from __future__ import annotations

import html
import json
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import pandas as pd

from autarch.classifier import ReactionClassifier

TOP_LEVEL_EC_IDS = [str(index) for index in range(1, 8)]
GO_PREVIEW_LIMIT = 8


def go_label_sort_key(go_id: str, go_labels: dict[str, str]) -> tuple[str, str]:
    """Provide a stable sort key for GO identifiers using labels when available."""
    return ((go_labels.get(go_id) or go_id).lower(), go_id)


@dataclass
class ClassifierReportEntry:
    """Classifier metadata rendered within an EC node."""

    class_name: str
    page_href: str
    go_id: Optional[str]
    go_label: Optional[str]
    declared_ecs: list[str]
    declared_ec_labels: list[str]
    placement_note: Optional[str] = None
    f1_score: Optional[float] = None
    precision: Optional[float] = None
    recall: Optional[float] = None
    mcc: Optional[float] = None
    total: Optional[int] = None
    cyclomatic_complexity: Optional[int] = None
    loc: Optional[int] = None
    complexity_score: Optional[float] = None
    declared_depth: int = 0


@dataclass
class ECNode:
    """A rendered EC hierarchy node."""

    ec_id: str
    label: str
    depth: int
    go_preview: list[tuple[str, str]]
    hidden_go_count: int
    subsumed_go_count: int
    subsumed_rhea_count: int
    classes: list[ClassifierReportEntry] = field(default_factory=list)
    children: list["ECNode"] = field(default_factory=list)
    subtree_avg_f1: Optional[float] = None
    subtree_f1_class_count: int = 0


def ec_sort_key(ec_id: str) -> tuple[int, ...]:
    """Sort EC identifiers numerically."""
    return tuple(int(part) for part in ec_id.split("."))


def ec_depth(ec_id: str) -> int:
    """Return EC hierarchy depth."""
    return len(ec_id.split("."))


def bioregistry_link(curie: str, display: Optional[str] = None) -> str:
    """Create a bioregistry link."""
    shown = html.escape(display or curie)
    return f'<a href="https://bioregistry.io/{curie}" target="_blank">{shown}</a>'


def get_stat_class(value: Optional[float], thresholds: tuple[float, float] = (0.7, 0.9)) -> str:
    """Map a metric to a CSS class."""
    if value is None:
        return "neutral"
    if value >= thresholds[1]:
        return "good"
    if value >= thresholds[0]:
        return "medium"
    return "bad"


def format_metric(value: Optional[float], digits: int = 3) -> str:
    """Format an optional numeric metric."""
    if value is None:
        return "n/a"
    return f"{value:.{digits}f}"


def format_integer(value: Optional[int]) -> str:
    """Format an optional integer metric."""
    if value is None:
        return "n/a"
    return str(value)


def format_ec_with_dashes(ec_value: str) -> str:
    """Render an EC identifier as a four-level CURIE."""
    cleaned = ec_value.replace("EC:", "")
    parts = cleaned.split(".")
    while len(parts) < 4:
        parts.append("-")
    return f"EC:{'.'.join(parts)}"


def normalize_ec_for_ontology(ec_value: Optional[str], max_depth: Optional[int] = None) -> Optional[str]:
    """Normalize EC identifiers to the sqlite:obo:ec ontology identifier shape.

    Examples:
        >>> normalize_ec_for_ontology("1.1.1.1")
        '1.1.1.1'
        >>> normalize_ec_for_ontology("1.1.-.-")
        '1.1'
        >>> normalize_ec_for_ontology("EC:2.7.1.-", max_depth=3)
        '2.7.1'
    """
    if not ec_value:
        return None

    parts: list[str] = []
    for part in ec_value.replace("EC:", "").split("."):
        if not part or part == "-":
            break
        parts.append(part)

    if max_depth is not None:
        parts = parts[:max_depth]

    if not parts:
        return None
    return ".".join(parts)


def iter_display_prefixes(ec_value: str, max_depth: int = 3) -> list[str]:
    """Return all displayed hierarchy prefixes covered by an EC identifier.

    Examples:
        >>> iter_display_prefixes("1.1.1.1")
        ['1', '1.1', '1.1.1']
        >>> iter_display_prefixes("3.6.-.-")
        ['3', '3.6']
    """
    normalized = normalize_ec_for_ontology(ec_value)
    if normalized is None:
        return []
    parts = normalized.split(".")
    return [".".join(parts[:depth]) for depth in range(1, min(len(parts), max_depth) + 1)]


def load_go_stats(
    cache_dir: Path,
) -> tuple[dict[str, str], dict[str, list[str]], dict[str, set[str]], dict[str, set[str]]]:
    """Load GO labels plus EC-derived node associations."""
    go_labels: dict[str, str] = {}
    go_ec_numbers: dict[str, list[str]] = {}
    exact_go_by_node: dict[str, set[str]] = defaultdict(set)
    subsumed_go_by_node: dict[str, set[str]] = defaultdict(set)

    cache_path = Path(cache_dir) / "go_terms.jsonl"
    if not cache_path.exists():
        raise FileNotFoundError(f"Missing GO cache: {cache_path}")

    with cache_path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            go_id = row["go_id"]
            go_labels[go_id] = str(row.get("label") or go_id)
            go_ec_numbers[go_id] = [str(ec_number) for ec_number in row.get("ec_numbers", [])]

            display_nodes: set[str] = set()
            for ec_number in row.get("ec_numbers", []):
                for node_id in iter_display_prefixes(ec_number):
                    subsumed_go_by_node[node_id].add(go_id)
                display_node = normalize_ec_for_ontology(ec_number, max_depth=3)
                if display_node:
                    display_nodes.add(display_node)

            for node_id in display_nodes:
                exact_go_by_node[node_id].add(go_id)

    return go_labels, go_ec_numbers, exact_go_by_node, subsumed_go_by_node


def load_rhea_counts(cache_dir: Path) -> dict[str, set[str]]:
    """Load subtree RHEA counts per EC node."""
    subsumed_rheas_by_node: dict[str, set[str]] = defaultdict(set)

    cache_path = Path(cache_dir) / "rhea_reactions.jsonl"
    if not cache_path.exists():
        raise FileNotFoundError(f"Missing RHEA cache: {cache_path}")

    with cache_path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            rhea_id = row["rhea_id"]
            for ec_number in row.get("ec_numbers", []):
                for node_id in iter_display_prefixes(ec_number):
                    subsumed_rheas_by_node[node_id].add(rhea_id)

    return subsumed_rheas_by_node


def build_go_preview(
    node_id: str,
    go_labels: dict[str, str],
    exact_go_by_node: dict[str, set[str]],
    subsumed_go_by_node: dict[str, set[str]],
    limit: int = GO_PREVIEW_LIMIT,
) -> tuple[list[tuple[str, str]], int]:
    """Build a preview list of mapped GO terms for an EC node."""
    exact_ids = sorted(
        exact_go_by_node.get(node_id, set()),
        key=lambda go_id: go_label_sort_key(go_id, go_labels),
    )
    other_ids = sorted(
        subsumed_go_by_node.get(node_id, set()) - set(exact_ids),
        key=lambda go_id: go_label_sort_key(go_id, go_labels),
    )
    preview_ids = (exact_ids + other_ids)[:limit]
    preview = [(go_id, go_labels.get(go_id, go_id)) for go_id in preview_ids]
    hidden_count = max(0, len(exact_ids) + len(other_ids) - len(preview_ids))
    return preview, hidden_count


def load_ec_labels() -> dict[str, str]:
    """Load labels for all EC terms up to level 3 from sqlite:obo:ec."""
    from oaklib import get_adapter  # type: ignore[import-untyped]

    adapter = get_adapter("sqlite:obo:ec")
    labels: dict[str, str] = {}

    for top_level_id in TOP_LEVEL_EC_IDS:
        top_curie = f"EC:{top_level_id}"
        labels[top_level_id] = adapter.label(top_curie) or top_level_id
        for descendant_curie in adapter.descendants(top_curie, predicates=["rdfs:subClassOf"]):
            ec_id = descendant_curie.removeprefix("EC:")
            if ec_depth(ec_id) <= 3:
                labels[ec_id] = adapter.label(descendant_curie) or ec_id

    return labels


def build_classifier_entries(
    metrics_df: pd.DataFrame,
    complexity_df: pd.DataFrame,
    go_labels: dict[str, str],
    go_ec_numbers: dict[str, list[str]],
    ec_labels: dict[str, str],
) -> dict[str, list[ClassifierReportEntry]]:
    """Associate Python classes with EC nodes."""
    classifier = ReactionClassifier()
    eval_lookup: dict[str, dict[str, Any]] = {
        str(row["class"]): row.to_dict() for _, row in metrics_df.iterrows()
    }
    complexity_lookup: dict[str, dict[str, Any]] = {
        str(row["name"]): row.to_dict() for _, row in complexity_df.iterrows()
    }
    entries_by_node: dict[str, list[ClassifierReportEntry]] = defaultdict(list)

    for class_name, cls in classifier.reaction_classes.items():
        go_id = getattr(cls, "GO_ID", None)
        explicit_ec_values: list[str] = []
        broad_ec_values: list[str] = []

        ec_prefix = getattr(cls, "EC_NUMBER_PREFIX", None)
        if ec_prefix:
            explicit_ec_values.append(ec_prefix)

        explicit_numbers = getattr(cls, "EC_NUMBERS", None) or []
        explicit_ec_values.extend(explicit_numbers)
        broad_ec_values.extend(getattr(cls, "EC_BROAD_XREFS", None) or [])

        node_ids: set[str] = set()
        placement_note: Optional[str] = None

        for ec_value in explicit_ec_values:
            node_id = normalize_ec_for_ontology(ec_value, max_depth=3)
            if node_id:
                node_ids.add(node_id)

        if explicit_ec_values and broad_ec_values:
            placement_note = "GO term has broad EC xref"

        if not node_ids and broad_ec_values:
            for ec_value in broad_ec_values:
                node_id = normalize_ec_for_ontology(ec_value, max_depth=3)
                if node_id:
                    node_ids.add(node_id)
            if node_ids:
                placement_note = "Placed via GO broad EC xref"

        if not node_ids and go_id:
            for ec_value in go_ec_numbers.get(go_id, []):
                node_id = normalize_ec_for_ontology(ec_value, max_depth=3)
                if node_id:
                    node_ids.add(node_id)
            if node_ids:
                placement_note = "EC placement inferred from GO xrefs"

        if not node_ids:
            continue

        declared_ecs = explicit_ec_values or broad_ec_values or sorted(node_ids, key=ec_sort_key)
        declared_ec_labels: list[str] = []
        declared_depth = 0
        for ec_value in declared_ecs:
            normalized = normalize_ec_for_ontology(ec_value)
            if normalized:
                declared_depth = max(declared_depth, ec_depth(normalized))
                declared_ec_labels.append(ec_labels.get(normalized, ""))
            else:
                declared_ec_labels.append("")

        eval_row: dict[str, Any] = eval_lookup.get(class_name, {})
        module_name = cls.__module__.split(".")[-1]
        complexity_row: dict[str, Any] = complexity_lookup.get(module_name, {})

        entry = ClassifierReportEntry(
            class_name=class_name,
            page_href=f"{class_name.lower()}.html",
            go_id=go_id,
            go_label=go_labels.get(go_id, go_id) if go_id else None,
            declared_ecs=declared_ecs,
            declared_ec_labels=declared_ec_labels,
            placement_note=placement_note,
            f1_score=eval_row.get("f1_score"),
            precision=eval_row.get("precision"),
            recall=eval_row.get("recall"),
            mcc=eval_row.get("mcc"),
            total=eval_row.get("total"),
            cyclomatic_complexity=complexity_row.get("cyclomatic_complexity"),
            loc=complexity_row.get("loc"),
            complexity_score=complexity_row.get("complexity_score"),
            declared_depth=declared_depth,
        )

        for node_id in sorted(node_ids, key=ec_sort_key):
            entries_by_node[node_id].append(entry)

    for node_entries in entries_by_node.values():
        node_entries.sort(
            key=lambda entry: (
                entry.declared_depth,
                -(entry.f1_score if entry.f1_score is not None else -1.0),
                entry.class_name,
            )
        )

    return entries_by_node


def build_ec_tree(
    ec_labels: dict[str, str],
    go_labels: dict[str, str],
    exact_go_by_node: dict[str, set[str]],
    subsumed_go_by_node: dict[str, set[str]],
    subsumed_rheas_by_node: dict[str, set[str]],
    classes_by_node: dict[str, list[ClassifierReportEntry]],
) -> list[ECNode]:
    """Build the nested EC tree."""
    level2_children: dict[str, list[str]] = defaultdict(list)
    level3_children: dict[str, list[str]] = defaultdict(list)

    for ec_id in ec_labels:
        depth = ec_depth(ec_id)
        if depth == 2:
            level2_children[ec_id.split(".")[0]].append(ec_id)
        elif depth == 3:
            level3_children[".".join(ec_id.split(".")[:2])].append(ec_id)

    for child_list in level2_children.values():
        child_list.sort(key=ec_sort_key)
    for child_list in level3_children.values():
        child_list.sort(key=ec_sort_key)

    def make_node(ec_id: str) -> ECNode:
        preview, hidden_go_count = build_go_preview(
            ec_id,
            go_labels,
            exact_go_by_node,
            subsumed_go_by_node,
        )
        node = ECNode(
            ec_id=ec_id,
            label=ec_labels.get(ec_id, ec_id),
            depth=ec_depth(ec_id),
            go_preview=preview,
            hidden_go_count=hidden_go_count,
            subsumed_go_count=len(subsumed_go_by_node.get(ec_id, set())),
            subsumed_rhea_count=len(subsumed_rheas_by_node.get(ec_id, set())),
            classes=list(classes_by_node.get(ec_id, [])),
        )
        if node.depth == 1:
            node.children = [make_node(child_id) for child_id in level2_children.get(ec_id, [])]
        elif node.depth == 2:
            node.children = [make_node(child_id) for child_id in level3_children.get(ec_id, [])]
        return node

    return [make_node(top_level_id) for top_level_id in TOP_LEVEL_EC_IDS]


def annotate_subtree_f1(nodes: list[ECNode]) -> None:
    """Annotate each node with mean F1 across directly placed and descendant classes.

    Classifiers are deduplicated by class name within each subtree so that
    multi-placed wrappers do not get counted multiple times.
    """

    def annotate_node(node: ECNode) -> dict[str, float]:
        subtree_scores: dict[str, float] = {}
        for entry in node.classes:
            if entry.f1_score is not None:
                subtree_scores.setdefault(entry.class_name, entry.f1_score)

        for child in node.children:
            child_scores = annotate_node(child)
            for class_name, f1_score in child_scores.items():
                subtree_scores.setdefault(class_name, f1_score)

        node.subtree_f1_class_count = len(subtree_scores)
        if subtree_scores:
            node.subtree_avg_f1 = sum(subtree_scores.values()) / len(subtree_scores)
        else:
            node.subtree_avg_f1 = None
        return subtree_scores

    for node in nodes:
        annotate_node(node)


def render_go_pills(node: ECNode) -> str:
    """Render GO preview pills for the detail pane."""
    if not node.go_preview:
        return '<div class="empty-state">No GO preview terms for this node.</div>'

    pills = []
    for go_id, label in node.go_preview:
        pills.append(
            (
                '<span class="pill pill-go">'
                f'{bioregistry_link(go_id, go_id)}'
                f'<span class="pill-label">{html.escape(label)}</span>'
                '</span>'
            )
        )

    if node.hidden_go_count:
        pills.append(
            f'<span class="pill pill-muted">+{node.hidden_go_count} additional mapped GO terms</span>'
        )

    return f'<div class="pill-row">{"".join(pills)}</div>'


def render_classifier_table(node: ECNode) -> str:
    """Render classifier stats in a compact table for the detail pane."""
    if not node.classes:
        return '<div class="empty-state">No Python classifier is placed directly at this EC node.</div>'

    rows: list[str] = []
    for entry in node.classes:
        if entry.go_id:
            go_label = entry.go_label or entry.go_id
            go_html = (
                f'{bioregistry_link(entry.go_id, entry.go_id)} '
                f'<span class="class-go-label">{html.escape(go_label)}</span>'
            )
        else:
            go_html = '<span class="muted">n/a</span>'

        declared_ecs: list[str] = []
        for ec_value, ec_label in zip(entry.declared_ecs, entry.declared_ec_labels, strict=False):
            ec_curie = format_ec_with_dashes(ec_value)
            if ec_label:
                declared_ecs.append(
                    f'{bioregistry_link(ec_curie, ec_curie)} '
                    f'<span class="declared-ec-label">{html.escape(ec_label)}</span>'
                )
            else:
                declared_ecs.append(bioregistry_link(ec_curie, ec_curie))

        note_html = (
            f'<div class="class-note">{html.escape(entry.placement_note)}</div>'
            if entry.placement_note
            else ''
        )
        declared_ec_html = '<br>'.join(declared_ecs) if declared_ecs else '<span class="muted">n/a</span>'
        if note_html:
            declared_ec_html = f"{declared_ec_html}{note_html}"

        rows.append(
            """
            <tr>
                <td><a class="class-link" href="{page_href}">{class_name}</a></td>
                <td>{go_html}</td>
                <td>{declared_ecs}</td>
                <td class="{f1_class}">{f1}</td>
                <td class="{precision_class}">{precision}</td>
                <td class="{recall_class}">{recall}</td>
                <td class="{mcc_class}">{mcc}</td>
                <td>{total}</td>
                <td>{cc}</td>
                <td>{loc}</td>
                <td>{complexity}</td>
            </tr>
            """.format(
                page_href=html.escape(entry.page_href),
                class_name=html.escape(entry.class_name),
                go_html=go_html,
                declared_ecs=declared_ec_html,
                f1=format_metric(entry.f1_score),
                precision=format_metric(entry.precision),
                recall=format_metric(entry.recall),
                mcc=format_metric(entry.mcc),
                total=format_integer(entry.total),
                cc=format_integer(entry.cyclomatic_complexity),
                loc=format_integer(entry.loc),
                complexity=format_metric(entry.complexity_score, digits=1),
                f1_class=get_stat_class(entry.f1_score),
                precision_class=get_stat_class(entry.precision),
                recall_class=get_stat_class(entry.recall),
                mcc_class=get_stat_class(entry.mcc, thresholds=(0.4, 0.7)),
            )
        )

    return """
    <table class="class-table">
        <thead>
            <tr>
                <th>Class</th>
                <th>GO</th>
                <th>Declared EC</th>
                <th>F1</th>
                <th>P</th>
                <th>R</th>
                <th>MCC</th>
                <th>N</th>
                <th>CC</th>
                <th>LOC</th>
                <th>Complexity</th>
            </tr>
        </thead>
        <tbody>{rows}</tbody>
    </table>
    """.format(rows=''.join(rows))


def build_search_text(node: ECNode) -> str:
    """Build a flattened search string for client-side filtering."""
    chunks = [f"EC:{node.ec_id}", node.label]
    for go_id, go_label in node.go_preview:
        chunks.extend([go_id, go_label])
    for entry in node.classes:
        chunks.append(entry.class_name)
        if entry.go_id:
            chunks.append(entry.go_id)
        if entry.go_label:
            chunks.append(entry.go_label)
        chunks.extend(entry.declared_ecs)
        chunks.extend(label for label in entry.declared_ec_labels if label)
    return ' '.join(chunks).lower()


def build_row_go_snippet(node: ECNode) -> str:
    """Build a compact one-line GO snippet for the tree row."""
    if not node.go_preview:
        return '<span class="row-go-snippet muted">no GO preview</span>'
    go_id, label = node.go_preview[0]
    suffix = f" +{node.subsumed_go_count - 1}" if node.subsumed_go_count > 1 else ''
    full_label = f"{go_id} {label}{suffix}"
    return (
        f'<span class="row-go-snippet" title="{html.escape(full_label, quote=True)}">'
        f'{html.escape(full_label)}'
        '</span>'
    )


def render_detail_panel_content(node: ECNode) -> str:
    """Render the selected-node detail pane content."""
    ec_link = bioregistry_link(f"EC:{node.ec_id}", f"EC:{node.ec_id}")
    return (
        f'<div class="detail-header"><div class="detail-ec">{ec_link}</div>'
        f'<div class="detail-label">{html.escape(node.label)}</div></div>'
        '<div class="detail-summary">'
        f'<span class="count-pill">Subsumed GO {node.subsumed_go_count}</span>'
        f'<span class="count-pill">Subsumed RHEAs {node.subsumed_rhea_count}</span>'
        f'<span class="count-pill">Py classes {len(node.classes)}</span>'
        f'<span class="count-pill {get_stat_class(node.subtree_avg_f1)}">Avg F1 {format_metric(node.subtree_avg_f1)} ({node.subtree_f1_class_count})</span>'
        '</div>'
        '<div class="detail-section">'
        '<div class="section-title">GO preview</div>'
        f'{render_go_pills(node)}'
        '</div>'
        '<div class="detail-section">'
        '<div class="section-title">Python classes</div>'
        f'{render_classifier_table(node)}'
        '</div>'
    )


def collect_nodes(nodes: list[ECNode]) -> list[ECNode]:
    """Flatten the EC tree."""
    ordered: list[ECNode] = []
    for node in nodes:
        ordered.append(node)
        ordered.extend(collect_nodes(node.children))
    return ordered


def render_tree_node(node: ECNode) -> str:
    """Render a single tree node row plus nested children."""
    children_html = ''
    if node.children:
        children_html = '<ul class="tree-children">' + ''.join(render_tree_node(child) for child in node.children) + '</ul>'

    toggle = (
        '<button class="toggle-button" type="button" aria-label="Toggle children">▾</button>'
        if node.children
        else '<span class="toggle-spacer"></span>'
    )
    default_class = '' if node.depth == 1 else ' collapsed'
    return (
        f'<li class="tree-node depth-{node.depth}{default_class}" data-depth="{node.depth}" '
        f'data-search="{html.escape(build_search_text(node), quote=True)}" '
        f'data-ec-id="{html.escape(node.ec_id, quote=True)}">'
        f'<div class="tree-row" data-ec-id="{html.escape(node.ec_id, quote=True)}">'
        f'{toggle}'
        f'<span class="row-ec-id">{bioregistry_link(f"EC:{node.ec_id}", f"EC:{node.ec_id}")}</span>'
        f'<span class="row-ec-label" title="{html.escape(node.label, quote=True)}">{html.escape(node.label)}</span>'
        f'{build_row_go_snippet(node)}'
        '<span class="row-metrics">'
        f'<span class="row-pill">GO {node.subsumed_go_count}</span>'
        f'<span class="row-pill">RHEA {node.subsumed_rhea_count}</span>'
        f'<span class="row-pill">Py {len(node.classes)}</span>'
        f'<span class="row-pill {get_stat_class(node.subtree_avg_f1)}">Avg F1 {format_metric(node.subtree_avg_f1)}</span>'
        '</span>'
        '</div>'
        f'{children_html}'
        '</li>'
    )


def render_ec_hierarchy_report(nodes: list[ECNode]) -> str:
    """Render the full one-page EC hierarchy report."""
    total_nodes = 0
    nodes_with_go = 0
    nodes_with_rheas = 0
    nodes_with_classes = 0
    unique_class_names: set[str] = set()

    stack = list(nodes)
    while stack:
        node = stack.pop()
        total_nodes += 1
        if node.subsumed_go_count:
            nodes_with_go += 1
        if node.subsumed_rhea_count:
            nodes_with_rheas += 1
        if node.classes:
            nodes_with_classes += 1
            unique_class_names.update(entry.class_name for entry in node.classes)
        stack.extend(node.children)

    flat_nodes = collect_nodes(nodes)
    tree_html = '<ul class="tree-root">' + ''.join(render_tree_node(node) for node in nodes) + '</ul>'
    detail_map = {node.ec_id: render_detail_panel_content(node) for node in flat_nodes}
    initial_detail = (
        render_detail_panel_content(nodes[0])
        if nodes
        else '<div class="empty-state">No EC nodes loaded.</div>'
    )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Autarch EC Hierarchy Report</title>
    <style>
        :root {{
            --bg-primary: #0d1117;
            --bg-secondary: #161b22;
            --bg-tertiary: #21262d;
            --text-primary: #c9d1d9;
            --text-secondary: #8b949e;
            --text-link: #58a6ff;
            --border-color: #30363d;
            --green: #3fb950;
            --red: #f85149;
            --yellow: #d29922;
        }}
        * {{ box-sizing: border-box; }}
        body {{
            margin: 0;
            padding: 24px;
            background: var(--bg-primary);
            color: var(--text-primary);
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif;
            line-height: 1.5;
        }}
        a {{ color: var(--text-link); text-decoration: none; }}
        a:hover {{ text-decoration: underline; }}
        .container {{ max-width: 1720px; margin: 0 auto; }}
        h1 {{ margin-top: 0; margin-bottom: 10px; }}
        .lede {{
            color: var(--text-secondary);
            max-width: 1100px;
            margin-bottom: 20px;
        }}
        .nav {{ margin-bottom: 16px; }}
        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 12px;
            margin-bottom: 20px;
        }}
        .summary-card {{
            background: var(--bg-secondary);
            border: 1px solid var(--border-color);
            border-radius: 8px;
            padding: 14px;
        }}
        .summary-value {{ font-size: 1.6em; font-weight: 700; }}
        .summary-label {{ color: var(--text-secondary); margin-top: 4px; font-size: 0.9em; }}
        .controls {{
            display: flex;
            gap: 10px;
            flex-wrap: wrap;
            margin-bottom: 18px;
        }}
        .controls input {{
            flex: 1 1 360px;
            min-width: 260px;
            background: var(--bg-secondary);
            color: var(--text-primary);
            border: 1px solid var(--border-color);
            border-radius: 8px;
            padding: 10px 12px;
        }}
        .controls button {{
            background: var(--bg-secondary);
            color: var(--text-primary);
            border: 1px solid var(--border-color);
            border-radius: 8px;
            padding: 10px 14px;
            cursor: pointer;
        }}
        .controls button:hover {{ border-color: var(--text-link); }}
        .note {{ color: var(--text-secondary); margin-bottom: 18px; font-size: 0.95em; }}
        .layout {{
            display: grid;
            grid-template-columns: minmax(820px, 1.55fr) minmax(420px, 1fr);
            gap: 16px;
            align-items: start;
        }}
        .tree-panel, .detail-panel {{
            background: var(--bg-secondary);
            border: 1px solid var(--border-color);
            border-radius: 8px;
        }}
        .tree-panel {{ overflow: auto; }}
        .detail-panel {{ position: sticky; top: 18px; padding: 14px 16px 16px; }}
        .detail-panel-inner {{ max-height: calc(100vh - 80px); overflow: auto; }}
        .tree-root, .tree-children {{ list-style: none; margin: 0; padding: 0; }}
        .tree-children {{ padding-left: 18px; }}
        .tree-node.collapsed > .tree-children {{ display: none; }}
        .tree-row {{
            display: grid;
            grid-template-columns: 22px 78px minmax(260px, 1fr) minmax(260px, 0.95fr) auto;
            gap: 10px;
            align-items: center;
            padding: 6px 10px;
            border-bottom: 1px solid rgba(48, 54, 61, 0.45);
            cursor: pointer;
            min-width: 980px;
            white-space: nowrap;
        }}
        .tree-row:hover {{ background: rgba(88, 166, 255, 0.06); }}
        .tree-row.selected {{ background: rgba(88, 166, 255, 0.14); }}
        .toggle-button, .toggle-spacer {{
            width: 18px;
            height: 18px;
            display: inline-flex;
            align-items: center;
            justify-content: center;
        }}
        .toggle-button {{
            background: transparent;
            color: var(--text-secondary);
            border: none;
            cursor: pointer;
            padding: 0;
            font-size: 13px;
            transform: rotate(0deg);
            transition: transform 0.15s ease;
        }}
        .tree-node.collapsed > .tree-row .toggle-button {{ transform: rotate(-90deg); }}
        .row-ec-id {{ font-weight: 700; }}
        .row-ec-label, .row-go-snippet {{ overflow: hidden; text-overflow: ellipsis; }}
        .row-go-snippet {{ color: var(--text-secondary); font-size: 0.9em; }}
        .row-metrics {{ display: inline-flex; gap: 6px; justify-content: flex-end; }}
        .count-pill, .row-pill {{
            background: var(--bg-tertiary);
            border: 1px solid var(--border-color);
            border-radius: 999px;
            padding: 4px 10px;
            font-size: 0.85em;
            color: var(--text-secondary);
        }}
        .row-pill {{ padding: 3px 9px; }}
        .count-pill.good, .row-pill.good {{ color: var(--green); }}
        .count-pill.medium, .row-pill.medium {{ color: var(--yellow); }}
        .count-pill.bad, .row-pill.bad {{ color: var(--red); }}
        .detail-header {{ margin-bottom: 12px; }}
        .detail-ec {{ font-size: 1.1em; font-weight: 700; margin-bottom: 4px; }}
        .detail-label {{ color: var(--text-primary); }}
        .detail-summary {{ display: flex; gap: 8px; flex-wrap: wrap; margin-bottom: 16px; }}
        .detail-section + .detail-section {{ margin-top: 18px; }}
        .section-title {{
            font-size: 0.9em;
            font-weight: 700;
            color: var(--text-secondary);
            text-transform: uppercase;
            letter-spacing: 0.04em;
            margin-bottom: 8px;
        }}
        .pill-row {{ display: flex; gap: 8px; flex-wrap: wrap; }}
        .pill {{
            display: inline-flex;
            align-items: center;
            gap: 6px;
            background: var(--bg-tertiary);
            border: 1px solid var(--border-color);
            border-radius: 999px;
            padding: 5px 10px;
            font-size: 0.88em;
        }}
        .pill-label, .declared-ec-label, .muted, .empty-state {{ color: var(--text-secondary); }}
        .empty-state {{ font-size: 0.95em; line-height: 1.5; }}
        .class-go-label {{ color: var(--text-primary); }}
        .class-note {{ color: var(--text-secondary); margin-top: 4px; font-size: 0.88em; }}
        .class-table {{ width: 100%; border-collapse: collapse; font-size: 0.9em; }}
        .class-table th, .class-table td {{
            vertical-align: top;
            text-align: left;
            padding: 8px 10px;
            border-bottom: 1px solid rgba(48, 54, 61, 0.45);
        }}
        .class-table th {{
            color: var(--text-secondary);
            font-weight: 700;
            position: sticky;
            top: 0;
            background: var(--bg-secondary);
        }}
        .class-table td.good {{ color: var(--green); }}
        .class-table td.medium {{ color: var(--yellow); }}
        .class-table td.bad {{ color: var(--red); }}
        .class-link {{ font-weight: 700; }}
        @media (max-width: 1180px) {{
            .layout {{ grid-template-columns: 1fr; }}
            .detail-panel {{ position: static; }}
            .detail-panel-inner {{ max-height: none; }}
        }}
        @media (max-width: 720px) {{
            body {{ padding: 16px; }}
            .tree-row {{ min-width: 780px; }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="nav"><a href="index.html">← Back to Index</a></div>
        <h1>EC Hierarchy Report</h1>
        <div class="lede">
            Compact tree view of the EC hierarchy through level 3, with one row per EC term. Select a node to inspect mapped GO terms, subsumed RHEA coverage, and any directly aligned Autarch classifier stats in the side panel.
        </div>
        <div class="summary-grid">
            <div class="summary-card"><div class="summary-value">{total_nodes}</div><div class="summary-label">Displayed EC terms</div></div>
            <div class="summary-card"><div class="summary-value">{nodes_with_go}</div><div class="summary-label">Nodes with GO mappings</div></div>
            <div class="summary-card"><div class="summary-value">{nodes_with_rheas}</div><div class="summary-label">Nodes with RHEA mappings</div></div>
            <div class="summary-card"><div class="summary-value">{nodes_with_classes}</div><div class="summary-label">Nodes with Python classes</div></div>
            <div class="summary-card"><div class="summary-value">{len(unique_class_names)}</div><div class="summary-label">Python classes placed</div></div>
        </div>
        <div class="controls">
            <input id="search-box" type="search" placeholder="Filter by EC, GO, class name, or label">
            <button id="expand-all" type="button">Expand all</button>
            <button id="collapse-all" type="button">Collapse to level 1</button>
        </div>
        <div class="note">Counts are unique within each node and are not additive across siblings because GO terms and RHEA reactions can carry multiple EC mappings. Row-level Avg F1 is the mean F1 across directly placed and descendant Python classifiers under that EC node, deduplicated by class name.</div>
        <div class="layout">
            <div class="tree-panel">{tree_html}</div>
            <aside class="detail-panel"><div class="detail-panel-inner" id="detail-panel">{initial_detail}</div></aside>
        </div>
    </div>
    <script>
        const detailHtmlById = {json.dumps(detail_map)};
        const detailPanel = document.getElementById('detail-panel');
        const searchBox = document.getElementById('search-box');
        const treeNodes = Array.from(document.querySelectorAll('.tree-node'));
        const treeRows = Array.from(document.querySelectorAll('.tree-row'));

        function childNodes(node) {{
            return Array.from(node.querySelectorAll(':scope > .tree-children > .tree-node'));
        }}

        function selectNode(ecId) {{
            treeRows.forEach(row => {{
                row.classList.toggle('selected', row.dataset.ecId === ecId);
            }});
            detailPanel.innerHTML = detailHtmlById[ecId] || '<div class="empty-state">No details for this EC node.</div>';
        }}

        function applyFilter() {{
            const query = searchBox.value.trim().toLowerCase();
            const nodes = [...treeNodes].sort((left, right) => Number(right.dataset.depth) - Number(left.dataset.depth));
            for (const node of nodes) {{
                const ownMatch = query === '' || (node.dataset.search || '').includes(query);
                const visibleChild = childNodes(node).some(child => !child.hidden);
                const visible = ownMatch || visibleChild;
                node.hidden = !visible;
                if (query && visibleChild) {{
                    node.classList.remove('collapsed');
                }}
            }}
            const selected = document.querySelector('.tree-row.selected');
            if (!selected || selected.closest('.tree-node').hidden) {{
                const firstVisible = treeRows.find(row => !row.closest('.tree-node').hidden);
                if (firstVisible) {{
                    selectNode(firstVisible.dataset.ecId);
                }}
            }}
        }}

        document.querySelectorAll('.toggle-button').forEach(button => {{
            button.addEventListener('click', event => {{
                event.stopPropagation();
                button.closest('.tree-node').classList.toggle('collapsed');
            }});
        }});

        treeRows.forEach(row => {{
            row.addEventListener('click', () => {{
                selectNode(row.dataset.ecId);
            }});
        }});

        document.getElementById('expand-all').addEventListener('click', () => {{
            treeNodes.forEach(node => {{
                if (!node.hidden) {{
                    node.classList.remove('collapsed');
                }}
            }});
        }});

        document.getElementById('collapse-all').addEventListener('click', () => {{
            treeNodes.forEach(node => {{
                if (childNodes(node).length) {{
                    node.classList.add('collapsed');
                }}
            }});
        }});

        searchBox.addEventListener('input', applyFilter);

        if (treeRows.length) {{
            selectNode(treeRows[0].dataset.ecId);
        }}
    </script>
</body>
</html>
"""


def save_ec_hierarchy_report_html(
    output_file: Path,
    metrics_df: pd.DataFrame,
    complexity_df: pd.DataFrame,
    cache_dir: Path = Path("cache"),
) -> None:
    """Save the one-page EC hierarchy report."""
    ec_labels = load_ec_labels()
    go_labels, go_ec_numbers, exact_go_by_node, subsumed_go_by_node = load_go_stats(cache_dir)
    subsumed_rheas_by_node = load_rhea_counts(cache_dir)
    classes_by_node = build_classifier_entries(
        metrics_df,
        complexity_df,
        go_labels,
        go_ec_numbers,
        ec_labels,
    )
    nodes = build_ec_tree(
        ec_labels,
        go_labels,
        exact_go_by_node,
        subsumed_go_by_node,
        subsumed_rheas_by_node,
        classes_by_node,
    )
    annotate_subtree_f1(nodes)
    output_file.write_text(render_ec_hierarchy_report(nodes))
