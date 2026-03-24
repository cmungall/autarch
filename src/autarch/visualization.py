#!/usr/bin/env python3
"""
Comprehensive visualization suite for Autarch evaluation results.
Generates beautiful charts and insights from classifier performance data.
"""

import json
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Circle
import seaborn as sns
import plotly.graph_objects as go  # type: ignore[import-untyped]
import plotly.express as px  # type: ignore[import-untyped]
from plotly.subplots import make_subplots  # type: ignore[import-untyped]
import numpy as np
from pathlib import Path
import textwrap
from typing import Any, Optional, TypedDict

from rdkit import Chem
from rdkit.Chem import Draw

from autarch.complexity import analyze_all_classifiers, ComplexityReport
from autarch.curation_objective import (
    build_curation_objective_df,
    save_curation_objective_report,
)

def setup_style():
    plt.style.use('seaborn-v0_8-darkgrid')
    sns.set_palette("husl")
    plt.rcParams['figure.figsize'] = (12, 8)
    plt.rcParams['font.size'] = 11
    plt.rcParams['axes.titlesize'] = 14
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['ytick.labelsize'] = 10
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['figure.dpi'] = 100


def build_complexity_performance_df(
    eval_df: pd.DataFrame,
    complexity_report: ComplexityReport,
) -> pd.DataFrame:
    """Build a merged dataframe for performance/complexity scatter plots."""
    return build_curation_objective_df(eval_df, complexity_report)


def create_complexity_vs_performance_figure(
    eval_df: pd.DataFrame,
    complexity_report: ComplexityReport,
    class_page_links: bool = False,
) -> go.Figure:
    """Create an interactive F1-vs-complexity scatter plot."""
    df = build_complexity_performance_df(eval_df, complexity_report).copy()
    if df.empty:
        return go.Figure()

    df["plot_support"] = df["positive_support"].clip(lower=1)
    df["target_url"] = ""
    if class_page_links:
        df["target_url"] = df["class"].str.lower() + ".html"

    fig = px.scatter(
        df,
        x="cyclomatic_complexity",
        y="f1_score",
        size="plot_support",
        color="conservative_f1",
        custom_data=["target_url"],
        hover_name="class",
        hover_data={
            "positive_support": True,
            "precision": ":.3f",
            "recall": ":.3f",
            "f1_score": ":.3f",
            "conservative_f1": ":.3f",
            "f1_optimism_gap": ":.3f",
            "complexity_score": ":.1f",
            "curation_objective": ":.3f",
            "target_url": False,
        },
        labels={
            "cyclomatic_complexity": "Cyclomatic Complexity",
            "f1_score": "F1 Score",
            "positive_support": "Positive Support",
            "conservative_f1": "Conservative F1",
        },
        title="F1 versus implementation complexity",
        color_continuous_scale="Viridis",
        size_max=28,
    )

    median_cc = df["cyclomatic_complexity"].median()
    median_f1 = df["f1_score"].median()
    fig.add_vline(x=median_cc, line_dash="dash", line_color="gray", opacity=0.5)
    fig.add_hline(y=median_f1, line_dash="dash", line_color="gray", opacity=0.5)

    fig.update_traces(
        marker=dict(
            opacity=0.82,
            line=dict(width=0.6, color="rgba(15, 23, 42, 0.6)"),
        )
    )
    fig.update_layout(
        height=720,
        xaxis_title="Cyclomatic complexity",
        yaxis_title="F1 score",
        coloraxis_colorbar_title="Conservative F1",
    )
    return fig


def create_complexity_vs_performance(
    eval_df: pd.DataFrame,
    complexity_report: ComplexityReport,
    output_dir: str | Path,
    class_page_links: bool = False,
) -> Path:
    """Create an interactive HTML chart showing complexity/performance trade-off."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    html_file = output_path / "complexity_vs_performance.html"
    fig = create_complexity_vs_performance_figure(
        eval_df,
        complexity_report,
        class_page_links=class_page_links,
    )
    post_script = None
    if class_page_links:
        post_script = """
const plot = document.getElementById('{plot_id}');
plot.on('plotly_click', function(event) {
  const url = event?.points?.[0]?.customdata?.[0];
  if (url) {
    window.location.href = url;
  }
});
"""
    fig.write_html(
        str(html_file),
        include_plotlyjs="cdn",
        full_html=True,
        post_script=post_script,
    )
    return html_file


def create_complexity_vs_performance_static(
    eval_df: pd.DataFrame,
    complexity_report: ComplexityReport,
    output_file: str | Path,
) -> Path:
    """Create a static PNG scatter plot for manuscripts and summaries."""
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df = build_complexity_performance_df(eval_df, complexity_report).copy()
    if df.empty:
        raise ValueError("No complexity/performance data available for plotting")

    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    sizes = 25 + 18 * np.sqrt(np.maximum(df["positive_support"], 1))
    scatter = ax.scatter(
        df["cyclomatic_complexity"],
        df["f1_score"],
        c=df["positive_support"],
        s=sizes,
        cmap="viridis",
        alpha=0.82,
        edgecolors="black",
        linewidths=0.35,
    )
    ax.set_xlabel("Cyclomatic complexity")
    ax.set_ylabel("F1 score")
    ax.set_title("Classifier performance versus implementation complexity")
    ax.grid(True, alpha=0.25)

    median_cc = df["cyclomatic_complexity"].median()
    median_f1 = df["f1_score"].median()
    ax.axvline(median_cc, linestyle="--", color="gray", alpha=0.5)
    ax.axhline(median_f1, linestyle="--", color="gray", alpha=0.5)

    label_candidates = pd.concat(
        [
            df.nlargest(3, "cyclomatic_complexity"),
            df[df["positive_support"] >= 10].nlargest(3, "f1_score"),
        ]
    ).drop_duplicates(subset=["class"])
    for _, row in label_candidates.iterrows():
        ax.annotate(
            row["class"],
            (row["cyclomatic_complexity"], row["f1_score"]),
            textcoords="offset points",
            xytext=(4, 4),
            fontsize=8,
        )

    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label("Positive support (TP + FN)")

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path


def _infer_classifier_family(class_name: str) -> str:
    """Infer a broad reaction family label from a classifier name."""
    lowered = class_name.lower()
    if "translocase" in lowered or "transmembrane" in lowered or "transporter" in lowered:
        return "Transport"
    if "oxidoreductase" in lowered or "oxidase" in lowered or "dehydrogenase" in lowered or "peroxidase" in lowered:
        return "Oxidoreductase"
    if "transferase" in lowered or "kinase" in lowered:
        return "Transferase"
    if "hydrolase" in lowered or "phosphatase" in lowered or "peptidase" in lowered or "nuclease" in lowered:
        return "Hydrolase"
    if "lyase" in lowered or "synthase" in lowered:
        return "Lyase"
    if "ligase" in lowered:
        return "Ligase"
    if "isomerase" in lowered or "mutase" in lowered or "epimerase" in lowered or "racemase" in lowered:
        return "Isomerase"
    return "Other"


class OverviewBoxSpec(TypedDict):
    """Typed figure specification for a system-overview box."""

    xy: tuple[float, float]
    w: float
    h: float
    fc: str
    ec: str
    title: str
    lines: list[str]


def create_f1_vs_support_static(
    eval_df: pd.DataFrame,
    output_file: str | Path,
) -> Path:
    """Create a manuscript-ready scatter plot of per-class F1 versus positive support."""
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df = eval_df.copy()
    df["positive_support"] = df["tp"] + df["fn"]
    df["family"] = df["class"].map(_infer_classifier_family)
    df["plot_support"] = df["positive_support"].clip(lower=0) + 1

    palette = {
        "Oxidoreductase": "#a61e4d",
        "Transferase": "#0b7285",
        "Hydrolase": "#2b8a3e",
        "Lyase": "#e67700",
        "Ligase": "#5f3dc4",
        "Isomerase": "#c2255c",
        "Transport": "#1c7ed6",
        "Other": "#495057",
    }

    fig, ax = plt.subplots(figsize=(8.7, 6.4))
    for family, family_df in df.groupby("family"):
        family_label = str(family)
        ax.scatter(
            family_df["plot_support"],
            family_df["f1_score"],
            s=36,
            alpha=0.82,
            c=palette.get(family_label, "#495057"),
            edgecolors="white",
            linewidths=0.4,
            label=family_label,
        )

    ax.set_xscale("log")
    ax.set_xlabel("Positive support (TP + FN, log scale)")
    ax.set_ylabel("F1 score")
    ax.set_title("Per-class F1 versus positive support")
    ax.grid(True, alpha=0.25)
    ax.set_ylim(-0.02, 1.04)
    ax.axhline(df["f1_score"].median(), linestyle="--", color="gray", alpha=0.5)

    label_candidates = set(df.nlargest(3, "positive_support")["class"])
    label_candidates.update(df[(df["positive_support"] >= 10)].nlargest(4, "f1_score")["class"])
    label_candidates.update(df[(df["positive_support"] >= 20)].nsmallest(4, "f1_score")["class"])
    label_candidates.update(
        df[df["class"].isin(["AminoacylTRNALigase", "RNAPolymerase", "PolysialicAcidOAcetyltransferase"])]["class"]
    )

    offsets = [(5, 5), (5, -10), (-40, 6), (8, 10), (-18, -12), (6, 14), (-48, 10), (8, -14)]
    for idx, class_name in enumerate(sorted(label_candidates)):
        row = df[df["class"] == class_name].iloc[0]
        dx, dy = offsets[idx % len(offsets)]
        ax.annotate(
            class_name,
            (row["plot_support"], row["f1_score"]),
            textcoords="offset points",
            xytext=(dx, dy),
            fontsize=8,
        )

    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax.legend(
        unique.values(),
        unique.keys(),
        title="Family",
        loc="lower right",
        frameon=True,
        ncol=2,
    )

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path


def load_benchmark_snapshot_counts(
    cache_dir: str | Path = "cache",
    results_dir: str | Path = "eval-results",
) -> dict[str, int]:
    """Load the current benchmark snapshot counts used in the manuscript figures."""
    cache_path = Path(cache_dir) / "rhea_reactions.jsonl"
    results_path = Path(results_dir) / "evaluation_results.csv"
    if not cache_path.exists():
        raise FileNotFoundError(f"Missing RHEA cache: {cache_path}")
    if not results_path.exists():
        raise FileNotFoundError(f"Missing evaluation results: {results_path}")

    master_reactions = 0
    go_linked_reactions = 0
    with cache_path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            master_reactions += 1
            if row.get("go_terms"):
                go_linked_reactions += 1

    eval_df = pd.read_csv(results_path)
    return {
        "master_reactions": master_reactions,
        "go_linked_reactions": go_linked_reactions,
        "evaluated_classifiers": len(eval_df),
    }


def load_resource_overlap_counts(
    cache_dir: str | Path = "cache",
) -> dict[str, int]:
    """Load overlap counts for the assembled RHEA/GO/EC reaction corpus."""
    cache_path = Path(cache_dir) / "rhea_reactions.jsonl"
    if not cache_path.exists():
        raise FileNotFoundError(f"Missing RHEA cache: {cache_path}")

    total = 0
    go_count = 0
    ec_count = 0
    both_count = 0
    go_only = 0
    ec_only = 0

    with cache_path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            has_go = bool(row.get("go_terms"))
            has_ec = bool(row.get("ec_numbers"))
            if not (has_go or has_ec):
                continue
            total += 1
            if has_go:
                go_count += 1
            if has_ec:
                ec_count += 1
            if has_go and has_ec:
                both_count += 1
            elif has_go:
                go_only += 1
            elif has_ec:
                ec_only += 1

    return {
        "total": total,
        "go_count": go_count,
        "ec_count": ec_count,
        "both_count": both_count,
        "go_only": go_only,
        "ec_only": ec_only,
    }


def create_resource_overlap_figure(
    output_file: str | Path,
    cache_dir: str | Path = "cache",
) -> Path:
    """Create a manuscript figure summarizing RHEA/GO/EC resource overlap."""
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    counts = load_resource_overlap_counts(cache_dir=cache_dir)

    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    outer = FancyBboxPatch(
        (0.05, 0.08),
        0.90,
        0.80,
        boxstyle="round,pad=0.02,rounding_size=0.03",
        linewidth=1.6,
        edgecolor="#495057",
        facecolor="#f8f9fa",
    )
    ax.add_patch(outer)
    ax.text(
        0.5,
        0.84,
        f"RHEA master reactions retained in the assembled corpus ({counts['total']:,})",
        ha="center",
        va="center",
        fontsize=13,
        fontweight="bold",
        color="#212529",
    )

    go_circle = Circle((0.42, 0.45), 0.22, facecolor="#d0ebff", edgecolor="#1c7ed6", alpha=0.80, linewidth=1.6)
    ec_circle = Circle((0.58, 0.45), 0.22, facecolor="#d3f9d8", edgecolor="#2b8a3e", alpha=0.80, linewidth=1.6)
    ax.add_patch(go_circle)
    ax.add_patch(ec_circle)

    ax.text(0.31, 0.66, "GO-linked\ncatalytic activity", ha="center", va="center", fontsize=12, fontweight="bold", color="#1864ab")
    ax.text(0.69, 0.66, "EC-linked\nreactions", ha="center", va="center", fontsize=12, fontweight="bold", color="#2b8a3e")

    ax.text(0.31, 0.44, f"{counts['go_only']:,}\nGO only", ha="center", va="center", fontsize=14, fontweight="bold", color="#1864ab")
    ax.text(0.50, 0.45, f"{counts['both_count']:,}\nGO + EC", ha="center", va="center", fontsize=15, fontweight="bold", color="#212529")
    ax.text(0.69, 0.44, f"{counts['ec_only']:,}\nEC only", ha="center", va="center", fontsize=14, fontweight="bold", color="#2b8a3e")

    ax.text(
        0.5,
        0.15,
        (
            "The benchmark corpus is assembled from the union of GO-linked and EC-linked "
            "master RHEA reactions; GO-linked reactions drive the primary evaluation."
        ),
        ha="center",
        va="center",
        fontsize=10.5,
        color="#495057",
        wrap=True,
    )

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path


def _draw_class_box(
    ax: plt.Axes,
    x: float,
    y: float,
    w: float,
    h: float,
    title: str,
    lines: list[str],
    edge_color: str,
    face_color: str,
) -> None:
    """Draw a UML-like class box."""
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.01,rounding_size=0.02",
        linewidth=1.4,
        edgecolor=edge_color,
        facecolor=face_color,
    )
    ax.add_patch(patch)
    ax.text(
        x + w / 2,
        y + h - 0.05,
        title,
        ha="center",
        va="center",
        fontsize=11.5,
        fontweight="bold",
        color=edge_color,
    )
    ax.plot([x + 0.02, x + w - 0.02], [y + h - 0.085, y + h - 0.085], color=edge_color, linewidth=1.0, alpha=0.6)
    ax.text(
        x + 0.02,
        y + h - 0.11,
        "\n".join(lines),
        ha="left",
        va="top",
        fontsize=9.2,
        color="#212529",
        linespacing=1.35,
    )


def create_classifier_architecture_figure(
    output_file: str | Path,
) -> Path:
    """Create a UML-like overview of the classifier architecture."""
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10.5, 6.4))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    _draw_class_box(
        ax,
        0.03,
        0.56,
        0.22,
        0.26,
        "Reaction / Participant",
        [
            "+ left_participants",
            "+ right_participants",
            "+ polymer_index / location",
            "+ is_transport_reaction()",
        ],
        "#1c7ed6",
        "#e7f5ff",
    )
    _draw_class_box(
        ax,
        0.35,
        0.60,
        0.28,
        0.28,
        "ReactionClass",
        [
            "+ GO_ID / EC_NUMBER_PREFIX",
            "+ supports_evaluation(reaction)",
            "+ check_membership(reaction)",
            "# check_membership_impl(reaction)",
        ],
        "#2b8a3e",
        "#ebfbee",
    )
    _draw_class_box(
        ax,
        0.33,
        0.28,
        0.22,
        0.18,
        "Kinase",
        [
            "+ PATTERNS = [...]",
            "+ direct chemistry rules",
        ],
        "#0b7285",
        "#e3fafc",
    )
    _draw_class_box(
        ax,
        0.60,
        0.28,
        0.24,
        0.18,
        "CisTransIsomerase",
        [
            "+ CHILD_CLASSES = (...)",
            "+ explicit union wrapper",
        ],
        "#9c36b5",
        "#f8f0fc",
    )
    _draw_class_box(
        ax,
        0.58,
        0.05,
        0.18,
        0.14,
        "CisTransIsomerases",
        ["+ EC 5.2.1.- wrapper"],
        "#9c36b5",
        "#f8f0fc",
    )
    _draw_class_box(
        ax,
        0.79,
        0.05,
        0.18,
        0.14,
        "PeptidylProlyl\nCisTransIsomerase",
        ["+ specific child class"],
        "#9c36b5",
        "#f8f0fc",
    )
    _draw_class_box(
        ax,
        0.72,
        0.60,
        0.24,
        0.22,
        "Evaluation / Reports",
        [
            "+ classifier-specific evidence",
            "+ TP / FP / FN / TN",
            "+ HTML reports",
            "+ manuscript figures",
        ],
        "#e67700",
        "#fff4e6",
    )

    arrows = [
        ((0.25, 0.69), (0.35, 0.74), "#1c7ed6"),
        ((0.49, 0.60), (0.44, 0.46), "#2b8a3e"),
        ((0.53, 0.60), (0.69, 0.46), "#2b8a3e"),
        ((0.63, 0.74), (0.72, 0.71), "#e67700"),
        ((0.72, 0.28), (0.67, 0.19), "#9c36b5"),
        ((0.72, 0.28), (0.88, 0.19), "#9c36b5"),
    ]
    for start, end, color in arrows:
        ax.add_patch(
            FancyArrowPatch(
                start,
                end,
                arrowstyle="-|>",
                mutation_scale=14,
                linewidth=1.6,
                color=color,
                shrinkA=8,
                shrinkB=8,
            )
        )

    ax.text(0.50, 0.95, "Classifier architecture in Autarch", ha="center", va="center", fontsize=14, fontweight="bold")
    ax.text(0.50, 0.91, "Concrete chemistry classes and explicit union wrappers share the same evaluation interface.", ha="center", va="center", fontsize=10.5, color="#495057")

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path


def _load_rhea_reaction_record(
    rhea_id: int,
    cache_dir: str | Path = "cache",
) -> dict[str, Any]:
    """Load a single cached RHEA reaction record by master ID."""
    cache_path = Path(cache_dir) / "rhea_reactions.jsonl"
    if not cache_path.exists():
        raise FileNotFoundError(f"Missing RHEA cache: {cache_path}")
    with cache_path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            current_id = int(str(row["rhea_id"]).split(":")[-1])
            if current_id == rhea_id:
                return row
    raise ValueError(f"RHEA:{rhea_id} not found in {cache_path}")


def _participant_display_name(name: str) -> str:
    replacements = {
        "ATP(4-)": "ATP",
        "ADP(3-)": "ADP",
        "NADPH(4-)": "NADPH",
        "NADP(3-)": "NADP+",
        "dioxygen": "O2",
        "hydron": "H+",
        "2'-deoxyadenosine 5'-monophosphate(2-)": "dAMP",
    }
    return replacements.get(name, name)


def _draw_participant_image(ax: plt.Axes, participant: dict[str, Any], x: float, y: float, w: float, h: float) -> None:
    """Draw one participant as structure image plus text label."""
    smiles = participant.get("smiles")
    name = _participant_display_name(participant.get("name") or participant.get("chebi_id") or "participant")
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            image = Draw.MolToImage(mol, size=(220, 160))
            ax.imshow(np.asarray(image), extent=(x, x + w, y, y + h), aspect="auto", zorder=2)
    ax.text(
        x + w / 2,
        y - 0.015,
        "\n".join(textwrap.wrap(name, width=16)),
        ha="center",
        va="top",
        fontsize=8.5,
        color="#212529",
    )


def _draw_reaction_panel(
    ax: plt.Axes,
    row: dict[str, Any],
    title: str,
    go_text: str,
    ec_text: str,
) -> None:
    """Draw a single reaction example panel."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    left = row["reaction"]["left_participants"]
    right = row["reaction"]["right_participants"]

    ax.text(0.02, 0.96, title, ha="left", va="top", fontsize=12, fontweight="bold", color="#212529")

    left_positions = np.linspace(0.03, 0.34, len(left))
    right_positions = np.linspace(0.58, 0.89, len(right))
    y = 0.40
    w = 0.12
    h = 0.32

    for idx, participant in enumerate(left):
        x = float(left_positions[idx])
        _draw_participant_image(ax, participant, x, y, w, h)
        if idx < len(left) - 1:
            ax.text(x + w + 0.015, y + h / 2, "+", fontsize=16, fontweight="bold", ha="center", va="center", color="#495057")

    ax.add_patch(
        FancyArrowPatch(
            (0.44, y + h / 2),
            (0.56, y + h / 2),
            arrowstyle="simple",
            mutation_scale=18,
            linewidth=0,
            facecolor="#495057",
            edgecolor="#495057",
        )
    )

    for idx, participant in enumerate(right):
        x = float(right_positions[idx])
        _draw_participant_image(ax, participant, x, y, w, h)
        if idx < len(right) - 1:
            ax.text(x + w + 0.015, y + h / 2, "+", fontsize=16, fontweight="bold", ha="center", va="center", color="#495057")

    go_box = FancyBboxPatch((0.04, 0.03), 0.42, 0.20, boxstyle="round,pad=0.01,rounding_size=0.02", linewidth=1.2, edgecolor="#1c7ed6", facecolor="#e7f5ff")
    ec_box = FancyBboxPatch((0.54, 0.03), 0.42, 0.20, boxstyle="round,pad=0.01,rounding_size=0.02", linewidth=1.2, edgecolor="#2b8a3e", facecolor="#ebfbee")
    ax.add_patch(go_box)
    ax.add_patch(ec_box)
    ax.text(0.05, 0.205, "GO route", ha="left", va="center", fontsize=9.5, fontweight="bold", color="#1864ab")
    ax.text(0.55, 0.205, "EC route", ha="left", va="center", fontsize=9.5, fontweight="bold", color="#2b8a3e")
    ax.text(0.05, 0.168, "\n".join(textwrap.wrap(go_text, width=34)), ha="left", va="top", fontsize=8.4, color="#212529")
    ax.text(0.55, 0.168, "\n".join(textwrap.wrap(ec_text, width=34)), ha="left", va="top", fontsize=8.4, color="#212529")


def create_scope_conflict_examples_figure(
    output_file: str | Path,
    cache_dir: str | Path = "cache",
) -> Path:
    """Create manuscript figure showing concrete GO/EC scope-conflict reactions."""
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    kinase_row = _load_rhea_reaction_record(23452, cache_dir=cache_dir)
    monooxygenase_row = _load_rhea_reaction_record(26085, cache_dir=cache_dir)

    fig, axes = plt.subplots(2, 1, figsize=(12.5, 9.0))
    _draw_reaction_panel(
        axes[0],
        kinase_row,
        "A. RHEA:23452  2'-deoxyadenosine + ATP = dAMP + ADP + H+",
        "GO groups this reaction under deoxyadenosine kinase activity and related nucleobase-containing-compound kinase branches.",
        "EC routes the same chemistry through EC 2.7.1.-, the alcohol-group phosphotransferase branch.",
    )
    _draw_reaction_panel(
        axes[1],
        monooxygenase_row,
        "B. RHEA:26085  (4S)-limonene + NADPH + O2 + H+ = limonene 1,2-epoxide + NADP+ + H2O",
        "GO treats this as limonene 1,2-monooxygenase [NAD(P)H] activity, emphasizing oxygen incorporation into substrate.",
        "EC places it in EC 1.14.13.-, emphasizing paired-donor oxygenase chemistry with NAD(P)H as electron donor.",
    )
    fig.suptitle("Concrete reactions illustrating GO/EC scope differences", fontsize=14, fontweight="bold", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path


def create_system_overview_figure(
    output_file: str | Path,
    cache_dir: str | Path = "cache",
    results_dir: str | Path = "eval-results",
) -> Path:
    """Create a manuscript schematic summarizing the Autarch workflow."""
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    counts = load_benchmark_snapshot_counts(cache_dir=cache_dir, results_dir=results_dir)

    fig, ax = plt.subplots(figsize=(11.2, 4.2))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    boxes: list[OverviewBoxSpec] = [
        {
            "xy": (0.03, 0.18),
            "w": 0.2,
            "h": 0.64,
            "fc": "#e7f5ff",
            "ec": "#1c7ed6",
            "title": "Inputs",
            "lines": [
                "RHEA reactions",
                "GO / EC mappings",
                "ChEBI structures",
            ],
        },
        {
            "xy": (0.28, 0.18),
            "w": 0.2,
            "h": 0.64,
            "fc": "#ebfbee",
            "ec": "#2b8a3e",
            "title": "Normalization",
            "lines": [
                "Master RHEA IDs",
                "Participant normalization",
                "Polymer / location parsing",
                "Evidence-aware evaluability",
            ],
        },
        {
            "xy": (0.53, 0.18),
            "w": 0.2,
            "h": 0.64,
            "fc": "#fff4e6",
            "ec": "#e67700",
            "title": "Classifier Layer",
            "lines": [
                "Declarative patterns",
                "Curated Python rules",
                "Rule explanations",
                f"{counts['evaluated_classifiers']} evaluated classifiers",
            ],
        },
        {
            "xy": (0.78, 0.18),
            "w": 0.19,
            "h": 0.64,
            "fc": "#f8f0fc",
            "ec": "#9c36b5",
            "title": "Outputs",
            "lines": [
                "GO / EC class predictions",
                "Per-class metrics",
                "Complexity analysis",
                "HTML browser + reports",
            ],
        },
    ]

    for box in boxes:
        x, y = box["xy"]
        patch = FancyBboxPatch(
            (x, y),
            box["w"],
            box["h"],
            boxstyle="round,pad=0.012,rounding_size=0.02",
            linewidth=1.6,
            edgecolor=box["ec"],
            facecolor=box["fc"],
        )
        ax.add_patch(patch)
        ax.text(
            x + box["w"] / 2,
            y + box["h"] - 0.08,
            box["title"],
            ha="center",
            va="center",
            fontsize=13,
            fontweight="bold",
            color=box["ec"],
        )
        body = "\n".join(f"- {textwrap.fill(line, width=24)}" for line in box["lines"])
        ax.text(
            x + 0.025,
            y + box["h"] - 0.16,
            body,
            ha="left",
            va="top",
            fontsize=10.5,
            color="#212529",
            linespacing=1.5,
        )

    arrow_y = 0.5
    for x0, x1, color in [(0.235, 0.28, "#1c7ed6"), (0.485, 0.53, "#2b8a3e"), (0.735, 0.78, "#e67700")]:
        arrow = FancyArrowPatch(
            (x0, arrow_y),
            (x1, arrow_y),
            arrowstyle="simple",
            mutation_scale=18,
            linewidth=0,
            facecolor=color,
            edgecolor=color,
            alpha=0.9,
        )
        ax.add_patch(arrow)

    ax.text(
        0.5,
        0.92,
        "Autarch workflow for declarative biochemical reaction classification",
        ha="center",
        va="center",
        fontsize=14,
        fontweight="bold",
        color="#212529",
    )
    ax.text(
        0.5,
        0.08,
        (
            "Evaluation snapshot: "
            f"{counts['master_reactions']:,} master RHEA reactions, "
            f"{counts['go_linked_reactions']:,} GO-linked reactions, "
            f"{counts['evaluated_classifiers']:,} evaluated classifiers"
        ),
        ha="center",
        va="center",
        fontsize=10,
        color="#495057",
    )

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path

def load_data(results_dir="results"):
    results_path = Path(results_dir)
    
    # Load main evaluation results
    eval_df = pd.read_csv(results_path / "evaluation_results.csv")
    
    # Load detailed predictions if available
    detailed_df = None
    if (results_path / "detailed_predictions.csv").exists():
        detailed_df = pd.read_csv(results_path / "detailed_predictions.csv")
    
    return eval_df, detailed_df

def create_performance_barchart(df, output_dir):
    """Create comprehensive performance bar charts."""
    
    # Sort by F1 score for better visualization
    df_sorted = df.sort_values('f1_score', ascending=True)
    
    # Create subplots for different metrics
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Classifier Performance Metrics', fontsize=16, fontweight='bold')
    
    # 1. F1 Score comparison
    ax1 = axes[0, 0]
    bars1 = ax1.barh(df_sorted['class'], df_sorted['f1_score'])
    ax1.set_xlabel('F1 Score')
    ax1.set_title('F1 Score by Classifier')
    ax1.set_xlim(0, 1)
    # Color bars based on performance
    colors = ['red' if x < 0.3 else 'orange' if x < 0.6 else 'green' for x in df_sorted['f1_score']]
    for bar, color in zip(bars1, colors):
        bar.set_color(color)
    ax1.grid(True, alpha=0.3)
    
    # 2. Precision vs Recall
    ax2 = axes[0, 1]
    x = np.arange(len(df_sorted))
    width = 0.35
    ax2.barh(x - width/2, df_sorted['precision'], width, label='Precision', color='skyblue')
    ax2.barh(x + width/2, df_sorted['recall'], width, label='Recall', color='lightcoral')
    ax2.set_yticks(x)
    ax2.set_yticklabels(df_sorted['class'])
    ax2.set_xlabel('Score')
    ax2.set_title('Precision vs Recall')
    ax2.legend()
    ax2.set_xlim(0, 1)
    ax2.grid(True, alpha=0.3)
    
    # 3. Accuracy and MCC
    ax3 = axes[1, 0]
    ax3.barh(x - width/2, df_sorted['accuracy'], width, label='Accuracy', color='lightgreen')
    ax3.barh(x + width/2, df_sorted['mcc'], width, label='MCC', color='gold')
    ax3.set_yticks(x)
    ax3.set_yticklabels(df_sorted['class'])
    ax3.set_xlabel('Score')
    ax3.set_title('Accuracy vs Matthews Correlation Coefficient')
    ax3.legend()
    ax3.set_xlim(0, 1)
    ax3.grid(True, alpha=0.3)
    
    # 4. Sample distribution (TP, FP, TN, FN)
    ax4 = axes[1, 1]
    categories = ['True Positive', 'False Positive', 'True Negative', 'False Negative']
    # Calculate averages
    avg_values = [
        df['tp'].mean(),
        df['fp'].mean(),
        df['tn'].mean(),
        df['fn'].mean()
    ]
    bars4 = ax4.bar(categories, avg_values, color=['green', 'orange', 'blue', 'red'])
    ax4.set_ylabel('Average Count')
    ax4.set_title('Average Prediction Distribution')
    ax4.tick_params(axis='x', rotation=45)
    
    # Add value labels on bars
    for bar in bars4:
        height = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.1f}',
                ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/performance_metrics.png", dpi=300, bbox_inches='tight')
    plt.close()
    
def create_interactive_heatmap(df, output_dir):
    """Create an interactive heatmap using Plotly."""
    
    # Prepare data for heatmap
    metrics = ['precision', 'recall', 'f1_score', 'accuracy', 'mcc']
    heatmap_data = df[['class'] + metrics].set_index('class')[metrics].T
    
    # Create interactive heatmap
    fig = go.Figure(data=go.Heatmap(
        z=heatmap_data.values,
        x=heatmap_data.columns,
        y=heatmap_data.index,
        colorscale='RdYlGn',
        text=np.round(heatmap_data.values, 3),
        texttemplate="%{text}",
        textfont={"size": 10},
        colorbar=dict(title="Score"),
        hovertemplate='Classifier: %{x}<br>Metric: %{y}<br>Value: %{z:.3f}<extra></extra>'
    ))
    
    fig.update_layout(
        title='Performance Metrics Heatmap',
        xaxis_title='Classifier',
        yaxis_title='Metric',
        height=600,
        xaxis={'tickangle': 45}
    )
    
    fig.write_html(f"{output_dir}/interactive_heatmap.html")

def create_confusion_matrix_grid(df, output_dir):
    """Create a grid of confusion matrices for each classifier."""
    
    n_classifiers = len(df)
    cols = 4
    rows = (n_classifiers + cols - 1) // cols
    
    fig, axes = plt.subplots(rows, cols, figsize=(16, rows * 4))
    fig.suptitle('Confusion Matrices for All Classifiers', fontsize=16, fontweight='bold')
    
    # Flatten axes array for easier iteration
    axes = axes.flatten() if n_classifiers > 1 else [axes]
    
    for idx, (_, row) in enumerate(df.iterrows()):
        if idx < len(axes):
            ax = axes[idx]
            
            # Create confusion matrix
            cm = np.array([[row['tp'], row['fp']], 
                          [row['fn'], row['tn']]])
            
            # Plot
            sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', 
                       ax=ax, cbar=False,
                       xticklabels=['Predicted +', 'Predicted -'],
                       yticklabels=['Actual +', 'Actual -'])
            ax.set_title(f"{row['class']}\n(F1: {row['f1_score']:.3f})")
    
    # Hide unused subplots
    for idx in range(n_classifiers, len(axes)):
        axes[idx].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/confusion_matrices.png", dpi=300, bbox_inches='tight')
    plt.close()

def create_radar_chart(df, output_dir):
    """Create radar charts for top performers."""
    
    # Select top 5 classifiers by F1 score
    top_classifiers = df.nlargest(5, 'f1_score')
    
    # Metrics for radar chart
    metrics = ['precision', 'recall', 'specificity', 'accuracy', 'f1_score']
    
    # Create subplot with multiple radar charts
    fig = make_subplots(
        rows=2, cols=3,
        specs=[[{'type': 'polar'}] * 3] * 2,
        subplot_titles=top_classifiers['class'].tolist() + ['']
    )
    
    colors = px.colors.qualitative.Set1
    
    for idx, (_, row) in enumerate(top_classifiers.iterrows()):
        r_idx = idx // 3
        c_idx = idx % 3
        
        values = [row[metric] for metric in metrics]
        
        fig.add_trace(
            go.Scatterpolar(
                r=values + [values[0]],  # Close the polygon
                theta=metrics + [metrics[0]],
                fill='toself',
                name=row['class'],
                line_color=colors[idx],
                fillcolor=colors[idx],
                opacity=0.6
            ),
            row=r_idx + 1,
            col=c_idx + 1
        )
        
        fig.update_polars(
            radialaxis=dict(
                visible=True,
                range=[0, 1]
            ),
            row=r_idx + 1,
            col=c_idx + 1
        )
    
    fig.update_layout(
        title_text="Top 5 Classifiers - Performance Radar Charts",
        showlegend=False,
        height=800
    )
    
    fig.write_html(f"{output_dir}/radar_charts.html")
    # fig.show() - disabled for batch processing

def create_ec_class_analysis(df, output_dir):
    """Analyze performance by EC class patterns."""
    
    # Map classifiers to EC classes (approximate)
    ec_mapping = {
        'Oxidoreductase': 'EC 1',
        'OxidoreductaseActingOnTheCHOHGroupOfDonorsNADOrNADPAsAcceptor': 'EC 1',
        'OxidoreductaseActingOnAHemeGroupOfDonors': 'EC 1',
        'OxidoreductaseActingOnNADHOrNADPHQuinoneOrSimilarCompoundAsAcceptor': 'EC 1',
        'Transferase': 'EC 2',
        'Kinase': 'EC 2',
        'Transaminase': 'EC 2',
        'Hydrolase': 'EC 3',
        'Phosphatase': 'EC 3',
        'ATPHydrolysis': 'EC 3',
        'HydrolaseActingOnAcidSulfurNitrogenBonds': 'EC 3',
        'HydrolaseActingOnSulfurNitrogenBonds': 'EC 3',
        'Lyase': 'EC 4',
        'CarboxyLyase': 'EC 4',
        'Isomerase': 'EC 5',
        'Ligase': 'EC 6',
        'LigaseFormingCarbonNitrogenBonds': 'EC 6',
        'PrimaryActiveTransmembraneTransporter': 'EC 7'
    }
    
    df['ec_class'] = df['class'].map(ec_mapping).fillna('Other')
    
    # Group by EC class and calculate average metrics
    ec_grouped = df.groupby('ec_class')[['precision', 'recall', 'f1_score', 'accuracy']].mean()
    
    # Create grouped bar chart
    fig, ax = plt.subplots(figsize=(12, 6))
    ec_grouped.plot(kind='bar', ax=ax)
    ax.set_title('Average Performance by Enzyme Commission (EC) Class', fontsize=14, fontweight='bold')
    ax.set_xlabel('EC Class')
    ax.set_ylabel('Score')
    ax.legend(title='Metrics')
    ax.grid(True, alpha=0.3)
    plt.xticks(rotation=45)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/ec_class_analysis.png", dpi=300, bbox_inches='tight')
    plt.close()

def create_performance_distribution(df, output_dir):
    """Create distribution plots for performance metrics."""
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle('Performance Metric Distributions', fontsize=16, fontweight='bold')
    
    metrics = ['precision', 'recall', 'f1_score', 'accuracy', 'mcc', 'specificity']
    colors = ['blue', 'green', 'red', 'purple', 'orange', 'brown']
    
    for idx, (metric, color) in enumerate(zip(metrics, colors)):
        ax = axes[idx // 3, idx % 3]
        
        # Create violin plot with box plot overlay
        parts = ax.violinplot([df[metric].dropna()], positions=[1], 
                              showmeans=True, showmedians=True)
        
        # Color the violin
        for pc in parts['bodies']:
            pc.set_facecolor(color)
            pc.set_alpha(0.7)
        
        # Add scatter points
        y = df[metric].dropna()
        x = np.random.normal(1, 0.04, size=len(y))
        ax.scatter(x, y, alpha=0.5, s=30, color=color)
        
        # Add horizontal lines for reference
        ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)
        ax.axhline(y=0.8, color='green', linestyle='--', alpha=0.3)
        
        ax.set_title(f'{metric.capitalize()} Distribution')
        ax.set_ylabel('Score')
        ax.set_xlim(0.5, 1.5)
        ax.set_ylim(-0.05, 1.05)
        ax.set_xticks([])
        ax.grid(True, alpha=0.3)
        
        # Add statistics text
        mean_val = df[metric].mean()
        median_val = df[metric].median()
        std_val = df[metric].std()
        ax.text(0.55, 0.95, f'μ={mean_val:.3f}\nM={median_val:.3f}\nσ={std_val:.3f}',
                transform=ax.transAxes, fontsize=9,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/performance_distributions.png", dpi=300, bbox_inches='tight')
    plt.close()

def create_interactive_3d_scatter(df, output_dir):
    """Create an interactive 3D scatter plot of classifier performance."""

    # Ensure MCC is non-negative for size parameter (plotly requires positive values)
    df_plot = df.copy()
    df_plot['mcc_size'] = df_plot['mcc'].clip(lower=0.01)  # Minimum size of 0.01

    fig = px.scatter_3d(df_plot,
                        x='precision',
                        y='recall',
                        z='f1_score',
                        color='accuracy',
                        size='mcc_size',
                        hover_data=['class', 'tp', 'fp', 'tn', 'fn', 'mcc'],
                        text='class',
                        color_continuous_scale='Viridis',
                        title='3D Classifier Performance Space')
    
    fig.update_traces(textposition='top center', textfont_size=8)
    
    fig.update_layout(
        scene=dict(
            xaxis_title='Precision',
            yaxis_title='Recall',
            zaxis_title='F1 Score',
            camera=dict(
                eye=dict(x=1.5, y=1.5, z=1.5)
            )
        ),
        height=700
    )
    
    fig.write_html(f"{output_dir}/3d_performance_scatter.html")
    # fig.show() - disabled for batch processing

def create_sunburst_chart(df, output_dir):
    """Create a hierarchical sunburst chart of classifier performance."""
    
    # Categorize performance
    def categorize_performance(f1):
        if f1 >= 0.8:
            return 'Excellent'
        elif f1 >= 0.6:
            return 'Good'
        elif f1 >= 0.4:
            return 'Fair'
        elif f1 >= 0.2:
            return 'Poor'
        else:
            return 'Very Poor'
    
    df['performance_category'] = df['f1_score'].apply(categorize_performance)
    df['ec_class'] = df['class'].apply(lambda x: 'EC ' + str(hash(x) % 6 + 1))  # Simplified EC mapping
    
    # Prepare data for sunburst
    sunburst_data = []
    for _, row in df.iterrows():
        sunburst_data.append({
            'labels': row['class'],
            'parents': row['performance_category'],
            'values': row['f1_score'],
            'text': f"F1: {row['f1_score']:.3f}<br>Acc: {row['accuracy']:.3f}"
        })
    
    # Add category level
    for category in df['performance_category'].unique():
        sunburst_data.append({
            'labels': category,
            'parents': '',
            'values': 0,
            'text': category
        })
    
    sunburst_df = pd.DataFrame(sunburst_data)
    
    fig = go.Figure(go.Sunburst(
        labels=sunburst_df['labels'],
        parents=sunburst_df['parents'],
        values=sunburst_df['values'],
        text=sunburst_df['text'],
        branchvalues="total",
        marker=dict(
            colorscale='RdYlGn',
            cmid=0.5
        ),
        hovertemplate='<b>%{label}</b><br>%{text}<br>F1 Score: %{value:.3f}<extra></extra>'
    ))
    
    fig.update_layout(
        title='Classifier Performance Hierarchy',
        height=600
    )
    
    fig.write_html(f"{output_dir}/sunburst_performance.html")
    # fig.show() - disabled for batch processing

def create_complexity_chart(complexity_report: ComplexityReport, output_dir: str) -> None:
    """Create complexity visualization charts."""

    # Convert to DataFrame for easier plotting
    data = []
    for c in complexity_report.classifiers:
        data.append({
            'classifier': c.name,
            'cyclomatic_complexity': c.cyclomatic_complexity,
            'chebi_ids': c.chebi_ids,
            'exclusions': c.exclusions,
            'method_loc': c.method_loc,
            'cc_rank': c.cc_rank,
            'complexity_score': c.complexity_score,
        })

    df = pd.DataFrame(data)
    df = df.sort_values('cyclomatic_complexity', ascending=True)

    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Classifier Complexity Analysis', fontsize=16, fontweight='bold')

    # 1. Cyclomatic Complexity bar chart
    ax1 = axes[0, 0]
    colors = []
    for cc in df['cyclomatic_complexity']:
        if cc <= 10:
            colors.append('green')
        elif cc <= 20:
            colors.append('yellowgreen')
        elif cc <= 30:
            colors.append('orange')
        else:
            colors.append('red')

    bars = ax1.barh(df['classifier'], df['cyclomatic_complexity'], color=colors)
    ax1.set_xlabel('Cyclomatic Complexity')
    ax1.set_title('Cyclomatic Complexity by Classifier')
    ax1.axvline(x=10, color='green', linestyle='--', alpha=0.5, label='Target (<10)')
    ax1.axvline(x=20, color='orange', linestyle='--', alpha=0.5, label='Warning (>20)')
    ax1.legend()

    # 2. ChEBI IDs vs Exclusions scatter
    ax2 = axes[0, 1]
    scatter = ax2.scatter(df['chebi_ids'], df['exclusions'],
                         c=df['cyclomatic_complexity'], cmap='RdYlGn_r',
                         s=100, alpha=0.7)
    ax2.set_xlabel('ChEBI IDs (Hardcoded Chemical Knowledge)')
    ax2.set_ylabel('Exclusion Rules')
    ax2.set_title('Knowledge Complexity vs Rule Complexity')
    plt.colorbar(scatter, ax=ax2, label='Cyclomatic Complexity')

    # Add classifier labels for outliers
    for _, row in df.iterrows():
        if row['cyclomatic_complexity'] > 30 or row['chebi_ids'] > 20:
            ax2.annotate(row['classifier'], (row['chebi_ids'], row['exclusions']),
                        fontsize=8, alpha=0.7)

    # 3. Complexity distribution by rank
    ax3 = axes[1, 0]
    rank_order = ['A', 'B', 'C', 'D', 'E', 'F', '?']
    rank_colors = {'A': 'green', 'B': 'yellowgreen', 'C': 'yellow',
                   'D': 'orange', 'E': 'orangered', 'F': 'red', '?': 'gray'}
    rank_counts = df['cc_rank'].value_counts().reindex(rank_order, fill_value=0)
    bars = ax3.bar(rank_counts.index, rank_counts.values,
                  color=[rank_colors.get(r, 'gray') for r in rank_counts.index])
    ax3.set_xlabel('Complexity Rank')
    ax3.set_ylabel('Number of Classifiers')
    ax3.set_title('Distribution by Complexity Rank')

    # Add value labels
    for bar in bars:
        height = bar.get_height()
        if height > 0:
            ax3.text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(height)}', ha='center', va='bottom')

    # 4. Top 15 most complex
    ax4 = axes[1, 1]
    top_15 = df.nlargest(15, 'cyclomatic_complexity')
    ax4.barh(top_15['classifier'], top_15['cyclomatic_complexity'], color='coral')
    ax4.set_xlabel('Cyclomatic Complexity')
    ax4.set_title('Top 15 Most Complex Classifiers')

    plt.tight_layout()
    plt.savefig(f"{output_dir}/complexity_analysis.png", dpi=300, bbox_inches='tight')
    plt.close()


def generate_summary_report(df, output_dir, complexity_report: Optional[ComplexityReport] = None):
    """Generate a comprehensive summary report."""

    output_path = Path(output_dir)
    report_lines = []
    report_lines.append("=" * 80)
    report_lines.append("ARCTURUS CLASSIFIER EVALUATION - COMPREHENSIVE REPORT")
    report_lines.append("=" * 80)
    report_lines.append("")
    
    # Overall statistics
    report_lines.append("OVERALL PERFORMANCE STATISTICS")
    report_lines.append("-" * 40)
    report_lines.append(f"Total Classifiers Evaluated: {len(df)}")
    report_lines.append(f"Average F1 Score: {df['f1_score'].mean():.4f} (±{df['f1_score'].std():.4f})")
    report_lines.append(f"Average Accuracy: {df['accuracy'].mean():.4f} (±{df['accuracy'].std():.4f})")
    report_lines.append(f"Average Precision: {df['precision'].mean():.4f} (±{df['precision'].std():.4f})")
    report_lines.append(f"Average Recall: {df['recall'].mean():.4f} (±{df['recall'].std():.4f})")
    report_lines.append(f"Average MCC: {df['mcc'].mean():.4f} (±{df['mcc'].std():.4f})")
    report_lines.append("")
    
    # Top performers
    report_lines.append("TOP 5 PERFORMERS (by F1 Score)")
    report_lines.append("-" * 40)
    top_5 = df.nlargest(5, 'f1_score')[['class', 'f1_score', 'precision', 'recall', 'accuracy']]
    for idx, row in top_5.iterrows():
        report_lines.append(f"{row['class']:30} | F1: {row['f1_score']:.3f} | P: {row['precision']:.3f} | R: {row['recall']:.3f} | A: {row['accuracy']:.3f}")
    report_lines.append("")
    
    # Bottom performers
    report_lines.append("BOTTOM 5 PERFORMERS (by F1 Score)")
    report_lines.append("-" * 40)
    bottom_5 = df.nsmallest(5, 'f1_score')[['class', 'f1_score', 'precision', 'recall', 'accuracy']]
    for idx, row in bottom_5.iterrows():
        report_lines.append(f"{row['class']:30} | F1: {row['f1_score']:.3f} | P: {row['precision']:.3f} | R: {row['recall']:.3f} | A: {row['accuracy']:.3f}")
    report_lines.append("")
    
    # Classification challenges
    report_lines.append("CLASSIFICATION CHALLENGES")
    report_lines.append("-" * 40)
    high_fp = df.nlargest(3, 'fp')[['class', 'fp', 'precision']]
    report_lines.append("Highest False Positives:")
    for _, row in high_fp.iterrows():
        report_lines.append(f"  - {row['class']}: {row['fp']} FP (Precision: {row['precision']:.3f})")
    report_lines.append("")
    
    high_fn = df.nlargest(3, 'fn')[['class', 'fn', 'recall']]
    report_lines.append("Highest False Negatives:")
    for _, row in high_fn.iterrows():
        report_lines.append(f"  - {row['class']}: {row['fn']} FN (Recall: {row['recall']:.3f})")
    report_lines.append("")
    
    # Perfect classifiers
    perfect = df[(df['precision'] == 1.0) | (df['recall'] == 1.0)]
    if not perfect.empty:
        report_lines.append("PERFECT SCORES")
        report_lines.append("-" * 40)
        for _, row in perfect.iterrows():
            if row['precision'] == 1.0:
                report_lines.append(f"  - {row['class']}: Perfect Precision (no false positives)")
            if row['recall'] == 1.0:
                report_lines.append(f"  - {row['class']}: Perfect Recall (no false negatives)")
        report_lines.append("")
    
    # Recommendations
    report_lines.append("RECOMMENDATIONS FOR IMPROVEMENT")
    report_lines.append("-" * 40)
    
    low_recall = df[df['recall'] < 0.3]
    if not low_recall.empty:
        report_lines.append("Classifiers needing recall improvement (< 0.3):")
        for _, row in low_recall.iterrows():
            report_lines.append(f"  - {row['class']}: Current recall {row['recall']:.3f}")
    
    low_precision = df[df['precision'] < 0.3]
    if not low_precision.empty:
        report_lines.append("Classifiers needing precision improvement (< 0.3):")
        for _, row in low_precision.iterrows():
            report_lines.append(f"  - {row['class']}: Current precision {row['precision']:.3f}")

    # Complexity analysis section
    if complexity_report:
        curation_df = build_curation_objective_df(df, complexity_report)
        save_curation_objective_report(df, complexity_report, output_path)

        report_lines.append("")
        report_lines.append("CODE COMPLEXITY ANALYSIS")
        report_lines.append("-" * 40)
        report_lines.append(f"Total LOC: {complexity_report.total_loc}")
        report_lines.append(f"Average Cyclomatic Complexity: {complexity_report.avg_cc:.1f}")
        report_lines.append(f"Average ChEBI IDs: {complexity_report.avg_chebi_ids:.1f}")
        report_lines.append(f"Disabled classifiers: {sum(1 for c in complexity_report.classifiers if c.is_disabled)}")
        report_lines.append("")

        report_lines.append("Most Complex Classifiers (may need refactoring):")
        for c in complexity_report.most_complex[:5]:
            report_lines.append(f"  - {c.name}: CC={c.cyclomatic_complexity} (Rank {c.cc_rank}), {c.chebi_ids} ChEBI IDs, {c.exclusions} exclusions")
        report_lines.append("")

        report_lines.append("Simplest Active Classifiers (good examples):")
        for c in complexity_report.simplest[:5]:
            report_lines.append(f"  - {c.name}: CC={c.cyclomatic_complexity}, {c.chebi_ids} ChEBI IDs")
        report_lines.append("")

        report_lines.append("Complexity Ranks: A(1-5)=Simple, B(6-10)=Moderate, C(11-20)=Complex,")
        report_lines.append("                  D(21-30)=Very Complex, E(31-40)=Highly Complex, F(41+)=Extreme")

        report_lines.append("")
        report_lines.append("CURATION OBJECTIVE")
        report_lines.append("-" * 40)
        report_lines.append(
            "Objective = conservative_F1 - 0.15 * complexity_percentile"
        )
        report_lines.append(
            "Conservative F1 uses Wilson lower bounds on precision and recall to penalize small-support optimism."
        )
        report_lines.append("")

        report_lines.append("Top 5 curator objective scores:")
        for _, row in curation_df.head(5).iterrows():
            report_lines.append(
                f"  - {row['class']}: objective={row['curation_objective']:.3f}, "
                f"conservative_F1={row['conservative_f1']:.3f}, "
                f"raw_F1={row['f1_score']:.3f}, support={int(row['positive_support'])}"
            )

        report_lines.append("")
        report_lines.append("Highest overfitting-risk classes:")
        risk_df = curation_df.sort_values(
            ["f1_optimism_gap", "positive_support"], ascending=[False, True]
        ).head(5)
        for _, row in risk_df.iterrows():
            report_lines.append(
                f"  - {row['class']}: gap={row['f1_optimism_gap']:.3f}, "
                f"raw_F1={row['f1_score']:.3f}, conservative_F1={row['conservative_f1']:.3f}, "
                f"support={int(row['positive_support'])}"
            )

    report_lines.append("")
    report_lines.append("=" * 80)
    report_lines.append(f"Report generated: {pd.Timestamp.now()}")
    
    # Write report
    report_text = "\n".join(report_lines)
    with open(output_path / "comprehensive_report.txt", 'w') as f:
        f.write(report_text)
    
    print(report_text)
    
    return report_text

def main():
    setup_style()

    # Create output directory for visualizations
    output_dir = Path("results/visualizations")
    output_dir.mkdir(exist_ok=True, parents=True)

    # Load data
    print("Loading evaluation results...")
    eval_df, detailed_df = load_data()

    print(f"Loaded {len(eval_df)} classifier results")
    print("\nGenerating visualizations...\n")

    # Generate complexity analysis
    print("0. Analyzing code complexity...")
    complexity_report = analyze_all_classifiers()
    print(f"   Analyzed {len(complexity_report.classifiers)} classifiers")
    print(f"   Average Cyclomatic Complexity: {complexity_report.avg_cc:.1f}")

    # Generate all visualizations
    print("1. Creating performance bar charts...")
    create_performance_barchart(eval_df, output_dir)

    print("2. Creating interactive heatmap...")
    create_interactive_heatmap(eval_df, output_dir)

    print("3. Creating confusion matrix grid...")
    create_confusion_matrix_grid(eval_df, output_dir)

    print("4. Creating radar charts for top performers...")
    create_radar_chart(eval_df, output_dir)

    print("5. Analyzing by EC class...")
    create_ec_class_analysis(eval_df, output_dir)

    print("6. Creating performance distribution plots...")
    create_performance_distribution(eval_df, output_dir)

    print("7. Creating 3D performance scatter plot...")
    create_interactive_3d_scatter(eval_df, output_dir)

    print("8. Creating sunburst hierarchy chart...")
    create_sunburst_chart(eval_df, output_dir)

    print("9. Creating complexity analysis charts...")
    create_complexity_chart(complexity_report, output_dir)

    print("10. Creating complexity vs performance chart...")
    create_complexity_vs_performance(eval_df, complexity_report, output_dir)

    print("11. Generating comprehensive report...")
    generate_summary_report(eval_df, output_dir, complexity_report)

    # Save complexity report
    from autarch.complexity import save_complexity_report
    save_complexity_report(complexity_report, output_dir)

    print(f"\n All visualizations saved to {output_dir}")
    print("\nInteractive visualizations:")
    print(f"  - {output_dir}/interactive_heatmap.html")
    print(f"  - {output_dir}/radar_charts.html")
    print(f"  - {output_dir}/3d_performance_scatter.html")
    print(f"  - {output_dir}/sunburst_performance.html")
    print(f"  - {output_dir}/complexity_vs_performance.html")

    print("\nStatic visualizations:")
    print(f"  - {output_dir}/performance_metrics.png")
    print(f"  - {output_dir}/confusion_matrices.png")
    print(f"  - {output_dir}/ec_class_analysis.png")
    print(f"  - {output_dir}/performance_distributions.png")
    print(f"  - {output_dir}/complexity_analysis.png")

    print("\nReports:")
    print(f"  - {output_dir}/comprehensive_report.txt")
    print(f"  - {output_dir}/complexity_report.json")
    print(f"  - {output_dir}/complexity_report.csv")

if __name__ == "__main__":
    main()
