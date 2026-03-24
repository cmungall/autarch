"""HTML export utility for reaction classifiers.

Generates standalone HTML pages for each reaction class with:
- Source code of the classification routine
- Summary statistics (F1, precision, recall, etc.)
- Drilldown sections for TP/FP/TN/FN results
"""

import ast
import inspect
import json
import re
from pathlib import Path
from typing import Any, Optional

import pandas as pd

from autarch.classifier import ReactionClassifier
from autarch.complexity import analyze_all_classifiers
from autarch.datamodel import Participant, Reaction
from autarch.ec_hierarchy_report import save_ec_hierarchy_report_html
from autarch.rhea_text_embeddings import save_rhea_browser_html
from autarch.visualization import create_complexity_vs_performance

# Maximum number of items to show in each outcome section (TP, FP, FN, TN)
MAX_ITEMS_PER_SECTION = 50


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
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
            --blue: #58a6ff;
            --purple: #a371f7;
        }}
        * {{ box-sizing: border-box; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif;
            background: var(--bg-primary);
            color: var(--text-primary);
            margin: 0;
            padding: 20px;
            line-height: 1.6;
        }}
        .container {{ max-width: 1200px; margin: 0 auto; }}
        h1 {{ color: var(--text-primary); border-bottom: 1px solid var(--border-color); padding-bottom: 10px; }}
        h2 {{ color: var(--text-secondary); margin-top: 30px; }}
        a {{ color: var(--text-link); text-decoration: none; }}
        a:hover {{ text-decoration: underline; }}
        .nav {{ margin-bottom: 20px; padding: 10px 0; border-bottom: 1px solid var(--border-color); }}
        .nav a {{ margin-right: 15px; }}

        /* Metadata section */
        .metadata {{
            background: var(--bg-secondary);
            border: 1px solid var(--border-color);
            border-radius: 6px;
            padding: 15px;
            margin-bottom: 20px;
        }}
        .metadata-item {{ display: inline-block; margin-right: 30px; }}
        .metadata-label {{ color: var(--text-secondary); font-size: 0.9em; }}
        .metadata-value {{ color: var(--text-primary); font-weight: 600; }}

        /* Stats grid */
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
            gap: 15px;
            margin-bottom: 30px;
        }}
        .stat-card {{
            background: var(--bg-secondary);
            border: 1px solid var(--border-color);
            border-radius: 6px;
            padding: 15px;
            text-align: center;
        }}
        .stat-value {{ font-size: 1.8em; font-weight: 600; }}
        .stat-label {{ color: var(--text-secondary); font-size: 0.85em; margin-top: 5px; }}
        .stat-value.good {{ color: var(--green); }}
        .stat-value.medium {{ color: var(--yellow); }}
        .stat-value.bad {{ color: var(--red); }}

        /* Confusion matrix */
        .confusion-matrix {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 10px;
            max-width: 400px;
            margin: 20px 0;
        }}
        .cm-cell {{
            background: var(--bg-secondary);
            border: 1px solid var(--border-color);
            border-radius: 6px;
            padding: 15px;
            text-align: center;
        }}
        .cm-cell.tp {{ border-color: var(--green); }}
        .cm-cell.tn {{ border-color: var(--green); }}
        .cm-cell.fp {{ border-color: var(--red); }}
        .cm-cell.fn {{ border-color: var(--yellow); }}
        .cm-count {{ font-size: 1.5em; font-weight: 600; }}
        .cm-label {{ color: var(--text-secondary); font-size: 0.85em; }}

        /* Code block */
        .code-section {{ margin: 20px 0; }}
        pre {{
            background: var(--bg-secondary);
            border: 1px solid var(--border-color);
            border-radius: 6px;
            padding: 15px;
            overflow-x: auto;
            font-family: 'SFMono-Regular', Consolas, 'Liberation Mono', Menlo, monospace;
            font-size: 0.9em;
            line-height: 1.5;
        }}
        code {{ color: var(--text-primary); }}
        .code-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 10px;
        }}
        .file-path {{ color: var(--text-secondary); font-size: 0.85em; }}

        /* Syntax highlighting */
        .kw {{ color: var(--red); }}
        .str {{ color: #a5d6ff; }}
        .num {{ color: #79c0ff; }}
        .comment {{ color: var(--text-secondary); font-style: italic; }}
        .func {{ color: var(--purple); }}
        .cls {{ color: var(--yellow); }}

        /* Collapsible sections */
        .collapsible {{
            background: var(--bg-secondary);
            border: 1px solid var(--border-color);
            border-radius: 6px;
            margin: 10px 0;
        }}
        .collapsible-header {{
            padding: 15px;
            cursor: pointer;
            display: flex;
            justify-content: space-between;
            align-items: center;
            user-select: none;
        }}
        .collapsible-header:hover {{ background: var(--bg-tertiary); }}
        .collapsible-title {{ font-weight: 600; }}
        .collapsible-count {{
            background: var(--bg-tertiary);
            padding: 2px 10px;
            border-radius: 20px;
            font-size: 0.85em;
        }}
        .collapsible-content {{
            display: none;
            padding: 0 15px 15px;
            border-top: 1px solid var(--border-color);
        }}
        .collapsible.open .collapsible-content {{ display: block; }}
        .collapsible-icon {{ transition: transform 0.2s; }}
        .collapsible.open .collapsible-icon {{ transform: rotate(90deg); }}

        /* Reaction items */
        .reaction-item {{
            padding: 10px 0;
            border-bottom: 1px solid var(--border-color);
        }}
        .reaction-item:last-child {{ border-bottom: none; }}
        .reaction-id {{ color: var(--blue); font-weight: 600; }}
        .reaction-label {{ color: var(--text-primary); }}
        .go-tag {{
            color: var(--purple);
            font-size: 0.85em;
            margin-left: 5px;
        }}
        .go-tag a {{ color: var(--purple); }}
        .reaction-explanation {{
            color: var(--text-secondary);
            font-size: 0.9em;
            margin-top: 5px;
            padding-left: 15px;
            border-left: 2px solid var(--border-color);
        }}
        .reaction-meta {{
            color: var(--text-secondary);
            font-size: 0.85em;
            margin-top: 4px;
        }}
        .reaction-meta a {{
            color: var(--text-link);
        }}

        /* Category colors */
        .tp-header {{ border-left: 3px solid var(--green); }}
        .tn-header {{ border-left: 3px solid var(--green); }}
        .fp-header {{ border-left: 3px solid var(--red); }}
        .fn-header {{ border-left: 3px solid var(--yellow); }}

        /* Index page */
        .class-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(300px, 1fr));
            gap: 15px;
        }}
        .class-card {{
            background: var(--bg-secondary);
            border: 1px solid var(--border-color);
            border-radius: 6px;
            padding: 15px;
            transition: border-color 0.2s;
        }}
        .class-card:hover {{ border-color: var(--text-link); }}
        .class-name {{ font-size: 1.1em; font-weight: 600; margin-bottom: 10px; }}
        .class-stats {{ display: flex; gap: 15px; }}
        .class-stat {{ font-size: 0.85em; }}
        .class-stat-label {{ color: var(--text-secondary); }}

        /* Definition section */
        .definition {{
            margin: 20px 0;
        }}
        .definition-content {{
            background: var(--bg-secondary);
            border: 1px solid var(--border-color);
            border-radius: 6px;
            padding: 15px;
            color: var(--text-primary);
            line-height: 1.7;
        }}
        .definition-content p {{
            margin: 0 0 10px 0;
        }}
        .definition-content p:last-child {{
            margin-bottom: 0;
        }}

        /* Truncation notice */
        .truncation-notice {{
            padding: 10px;
            margin-top: 10px;
            background: var(--bg-tertiary);
            border-radius: 4px;
            color: var(--text-secondary);
            font-style: italic;
            font-size: 0.9em;
        }}

        /* Footer */
        .footer {{
            margin-top: 40px;
            padding-top: 20px;
            border-top: 1px solid var(--border-color);
            color: var(--text-secondary);
            font-size: 0.85em;
        }}
        .plot-section {{
            margin: 30px 0;
        }}
        .section-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            gap: 10px;
            flex-wrap: wrap;
        }}
        .section-note {{
            color: var(--text-secondary);
            margin-top: 8px;
            margin-bottom: 14px;
        }}
        .plot-frame {{
            width: 100%;
            min-height: 760px;
            border: 1px solid var(--border-color);
            border-radius: 8px;
            background: var(--bg-secondary);
        }}
        .viz-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
            gap: 15px;
            margin: 20px 0 30px;
        }}
        .viz-card {{
            display: block;
            background: var(--bg-secondary);
            border: 1px solid var(--border-color);
            border-radius: 6px;
            padding: 16px;
        }}
        .viz-card:hover {{
            border-color: var(--text-link);
            text-decoration: none;
        }}
        .viz-card-title {{
            font-size: 1.05em;
            font-weight: 600;
            margin-bottom: 8px;
            color: var(--text-primary);
        }}
        .viz-card-text {{
            color: var(--text-secondary);
        }}
    </style>
</head>
<body>
    <div class="container">
        {content}
    </div>
    <script>
        document.querySelectorAll('.collapsible-header').forEach(header => {{
            header.addEventListener('click', () => {{
                header.parentElement.classList.toggle('open');
            }});
        }});
    </script>
</body>
</html>
"""


def get_stat_class(value: float, thresholds: tuple[float, float] = (0.7, 0.9)) -> str:
    """Return CSS class based on metric value."""
    if value >= thresholds[1]:
        return "good"
    elif value >= thresholds[0]:
        return "medium"
    return "bad"


def format_ec_with_dashes(ec: str) -> str:
    """Format EC number with dashes for missing parts.

    >>> format_ec_with_dashes("3")
    'EC:3.-.-.-'
    >>> format_ec_with_dashes("2.7")
    'EC:2.7.-.-'
    >>> format_ec_with_dashes("1.1.1")
    'EC:1.1.1.-'
    >>> format_ec_with_dashes("3.2.1.1")
    'EC:3.2.1.1'
    >>> format_ec_with_dashes("2.7.-.-")
    'EC:2.7.-.-'
    """
    if not ec or ec == "N/A":
        return ec

    # Remove EC: prefix if present
    ec = ec.replace("EC:", "")

    # Split into parts
    parts = ec.split(".")

    # Pad with dashes to make 4 parts
    while len(parts) < 4:
        parts.append("-")

    return f"EC:{'.'.join(parts)}"


def bioregistry_link(curie: str, display: Optional[str] = None) -> str:
    """Create a hyperlink to bioregistry.io for a CURIE.

    >>> bioregistry_link("GO:0016787")
    '<a href="https://bioregistry.io/GO:0016787" target="_blank">GO:0016787</a>'
    >>> bioregistry_link("RHEA:10000", "RHEA:10000 - some reaction")
    '<a href="https://bioregistry.io/RHEA:10000" target="_blank">RHEA:10000 - some reaction</a>'
    """
    if not curie or curie == "N/A":
        return curie or "N/A"

    display = display or curie
    return f'<a href="https://bioregistry.io/{curie}" target="_blank">{display}</a>'


def format_rhea_id(rhea_id: str) -> str:
    """Format RHEA ID as a bioregistry link.

    Handles both 'RHEA:10000' and '10000' formats.
    """
    if not rhea_id:
        return ""

    # Normalize to RHEA:XXXXX format
    if not rhea_id.startswith("RHEA:"):
        rhea_id = f"RHEA:{rhea_id}"

    return bioregistry_link(rhea_id)


def get_class_docstring(cls: type) -> str:
    """Extract the docstring from a class, cleaned up for display."""
    if not cls or not cls.__doc__:
        return ""

    doc = cls.__doc__.strip()

    # Remove excessive whitespace but preserve paragraph breaks
    lines = doc.split("\n")
    cleaned_lines = []
    for line in lines:
        cleaned_lines.append(line.strip())

    return "\n".join(cleaned_lines)


def simple_syntax_highlight(code: str) -> str:
    """Apply basic syntax highlighting to Python code.

    Uses a token-based approach to avoid nested span issues.
    """
    import html

    # First escape HTML
    code = html.escape(code)

    keywords = {
        'def', 'class', 'return', 'if', 'elif', 'else', 'for', 'while',
        'in', 'not', 'and', 'or', 'True', 'False', 'None', 'import',
        'from', 'as', 'try', 'except', 'finally', 'with', 'raise',
        'pass', 'break', 'continue', 'lambda', 'yield', 'assert', 'self',
    }

    # Use placeholders to avoid re-processing
    placeholder_id = 0
    placeholders: dict[str, str] = {}

    def make_placeholder(content: str) -> str:
        nonlocal placeholder_id
        key = f"\x00PH{placeholder_id}\x00"
        placeholder_id += 1
        placeholders[key] = content
        return key

    # Highlight triple-quoted strings first
    def replace_triple_string(m: re.Match) -> str:
        return make_placeholder(f'<span class="str">{m.group(0)}</span>')

    code = re.sub(r'&quot;&quot;&quot;.*?&quot;&quot;&quot;', replace_triple_string, code, flags=re.DOTALL)
    code = re.sub(r"&#x27;&#x27;&#x27;.*?&#x27;&#x27;&#x27;", replace_triple_string, code, flags=re.DOTALL)

    # Highlight comments
    def replace_comment(m: re.Match) -> str:
        return make_placeholder(f'<span class="comment">{m.group(0)}</span>')

    code = re.sub(r'#.*$', replace_comment, code, flags=re.MULTILINE)

    # Highlight single-quoted strings
    def replace_string(m: re.Match) -> str:
        return make_placeholder(f'<span class="str">{m.group(0)}</span>')

    code = re.sub(r'&quot;[^&]*?&quot;', replace_string, code)
    code = re.sub(r"&#x27;[^&]*?&#x27;", replace_string, code)

    # Highlight keywords (only whole words not already in placeholders)
    def replace_keyword(m: re.Match) -> str:
        word = m.group(1)
        if word in keywords:
            return make_placeholder(f'<span class="kw">{word}</span>')
        return word

    code = re.sub(r'\b(\w+)\b', replace_keyword, code)

    # Highlight numbers
    def replace_number(m: re.Match) -> str:
        return make_placeholder(f'<span class="num">{m.group(0)}</span>')

    code = re.sub(r'\b(\d+\.?\d*)\b', replace_number, code)

    # Restore placeholders
    for key, value in placeholders.items():
        code = code.replace(key, value)

    return code


def extract_method_source(cls: type, method_name: str) -> tuple[str, str, int]:
    """Extract source code for a method from a class.

    Returns:
        Tuple of (source_code, file_path, start_line)
    """
    try:
        source_file = inspect.getfile(cls)
        source_lines, class_start = inspect.getsourcelines(cls)

        # Parse the class to find the method
        with open(source_file) as f:
            tree = ast.parse(f.read())

        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef) and node.name == cls.__name__:
                for item in node.body:
                    if isinstance(item, ast.FunctionDef) and item.name == method_name:
                        start_line = item.lineno
                        end_line = item.end_lineno or start_line

                        with open(source_file) as f:
                            all_lines = f.readlines()

                        method_source = ''.join(all_lines[start_line - 1:end_line])
                        return method_source, source_file, start_line

        # Fallback: use inspect
        method = getattr(cls, method_name, None)
        if method:
            source = inspect.getsource(method)
            return source, source_file, class_start

    except Exception as e:
        return f"# Could not extract source: {e}", "", 0

    return "# Method not found", "", 0


def extract_class_context_source(
    cls: type, method_name: str
) -> tuple[str, str, int]:
    """Extract declarative class-level assignments used as classifier context.

    This surfaces the class attributes defined before the membership method,
    such as GO/EC identifiers, evaluation evidence declarations, and pattern
    lists. The goal is to make the declarative part of the classifier visible
    without dumping the entire module boilerplate into the HTML report.

    Returns:
        Tuple of (source_code, file_path, start_line)
    """
    try:
        source_file = inspect.getfile(cls)
        with open(source_file) as f:
            all_lines = f.readlines()
        with open(source_file) as f:
            tree = ast.parse(f.read())

        def assignment_target_names(node: ast.Assign | ast.AnnAssign) -> list[str]:
            if isinstance(node, ast.Assign):
                names: list[str] = []
                for target in node.targets:
                    if isinstance(target, ast.Name):
                        names.append(target.id)
                return names
            if isinstance(node.target, ast.Name):
                return [node.target.id]
            return []

        for node in ast.walk(tree):
            if not isinstance(node, ast.ClassDef) or node.name != cls.__name__:
                continue

            target_method_line = None
            for item in node.body:
                if isinstance(item, ast.FunctionDef) and item.name == method_name:
                    target_method_line = item.lineno
                    break
            if target_method_line is None:
                return "", source_file, 0

            snippets: list[str] = []
            start_line = 0
            for item in node.body:
                if item.lineno >= target_method_line:
                    break
                if not isinstance(item, (ast.Assign, ast.AnnAssign)):
                    continue

                target_names = assignment_target_names(item)
                if not target_names:
                    continue
                if not any(re.fullmatch(r"[A-Z0-9_]+", name) for name in target_names):
                    continue

                snippet_start = item.lineno
                snippet_end = item.end_lineno or snippet_start
                if start_line == 0:
                    start_line = snippet_start
                snippets.append("".join(all_lines[snippet_start - 1 : snippet_end]))

            return "\n".join(snippets).strip(), source_file, start_line

    except Exception as e:
        return f"# Could not extract class context: {e}", "", 0

    return "", "", 0


def get_aggregate_membership_context(cls: type | None) -> dict[str, Any] | None:
    """Return report metadata for aggregate wrapper classifiers.

    Aggregate wrappers are legitimate only if they dispatch over curated child
    classifiers. This helper makes that explicit on the HTML page so readers do
    not mistake declared EC/GO scope for benchmark leakage.
    """
    if cls is None:
        return None

    method_source, _, _ = extract_method_source(cls, "check_membership_impl")
    if "aggregate_ec_prefix_membership(" in method_source:
        from autarch.ontology.ec_prefix_aggregate import aggregate_child_classes

        ec_prefix = getattr(cls, "EC_NUMBER_PREFIX", None)
        if not isinstance(ec_prefix, str):
            return None
        child_classes = aggregate_child_classes(ec_prefix, cls.__name__)
        return {
            "title": "Aggregate Membership Strategy",
            "summary": (
                "This wrapper does not inspect reaction EC annotations. "
                "It delegates membership to curated child classifiers chosen "
                "statically from the ontology by declared EC scope."
            ),
            "scope_label": "Aggregate EC Scope",
            "scope_value": format_ec_with_dashes(ec_prefix),
            "children": [
                {
                    "name": child_cls.__name__,
                    "ec": format_ec_with_dashes(
                        getattr(child_cls, "EC_NUMBER_PREFIX", None)
                        or getattr(child_cls, "EC_NUMBER", None)
                        or ""
                    ),
                    "go_id": getattr(child_cls, "GO_ID", None),
                }
                for child_cls in child_classes
            ],
        }

    explicit_child_classes = getattr(cls, "CHILD_CLASSES", None)
    if isinstance(explicit_child_classes, tuple) and explicit_child_classes:
        return {
            "title": "Aggregate Membership Strategy",
            "summary": (
                "This wrapper does not inspect reaction GO or EC annotations. "
                "It delegates membership to the explicit curated child "
                "classifiers listed below."
            ),
            "scope_label": "Aggregate Children",
            "scope_value": str(len(explicit_child_classes)),
            "children": [
                {
                    "name": child_cls.__name__,
                    "ec": format_ec_with_dashes(
                        getattr(child_cls, "EC_NUMBER_PREFIX", None)
                        or getattr(child_cls, "EC_NUMBER", None)
                        or ""
                    ),
                    "go_id": getattr(child_cls, "GO_ID", None),
                }
                for child_cls in explicit_child_classes
            ],
        }

    return None


def load_evaluation_data(
    results_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, dict]]:
    """Load evaluation results from CSV files.

    Returns:
        Tuple of (metrics_df, detailed_df, rhea_reactions)
    """
    metrics_df = pd.read_csv(results_dir / "evaluation_results.csv")

    detailed_df = pd.DataFrame()
    detailed_file = results_dir / "detailed_predictions.csv"
    if detailed_file.exists():
        detailed_df = pd.read_csv(detailed_file)

    # Load RHEA reactions for labels
    rhea_reactions = {}
    cache_dir = Path("cache")
    rhea_cache = cache_dir / "rhea_reactions.jsonl"
    if rhea_cache.exists():
        with open(rhea_cache) as f:
            for line in f:
                data = json.loads(line)
                rhea_reactions[data["rhea_id"]] = data

    return metrics_df, detailed_df, rhea_reactions


def build_go_missing_positive_index(
    rhea_reactions: dict[str, dict[str, Any]],
    classifier: ReactionClassifier,
) -> dict[str, list[dict[str, Any]]]:
    """Index positive classifier matches for reactions lacking GO annotations."""
    from autarch.evaluation import ec_matches_prefix

    unknown_positive_index: dict[str, list[dict[str, Any]]] = {}

    def classify_annotation_bucket(
        cls: type | None, reaction_record: dict[str, Any]
    ) -> str:
        ec_numbers = reaction_record.get("ec_numbers") or []
        if not ec_numbers:
            return "unannotated"
        if cls is None:
            return "ec_conflict"

        class_ec_numbers = getattr(cls, "EC_NUMBERS", None) or []
        if any(ec_number in class_ec_numbers for ec_number in ec_numbers):
            return "ec_aligned"

        ec_prefix = getattr(cls, "EC_NUMBER_PREFIX", None)
        if ec_prefix and any(
            ec_matches_prefix(ec_number, ec_prefix) for ec_number in ec_numbers
        ):
            return "ec_aligned"

        return "ec_conflict"

    for rhea_id, reaction_record in sorted(rhea_reactions.items()):
        if reaction_record.get("go_terms"):
            continue
        reaction_data = reaction_record.get("reaction")
        if not reaction_data:
            continue

        reaction = Reaction(
            left_participants=[
                Participant(**participant)
                for participant in reaction_data.get("left_participants") or []
            ],
            right_participants=[
                Participant(**participant)
                for participant in reaction_data.get("right_participants") or []
            ],
            label=reaction_record.get("label", ""),
        )
        results = classifier.classify(reaction)
        for class_name, result in results.items():
            if not result.is_member:
                continue
            cls = classifier.reaction_classes.get(class_name)
            unknown_positive_index.setdefault(class_name, []).append(
                {
                    "rhea_id": rhea_id,
                    "label": reaction_record.get("label", "") or "No label",
                    "explanation": result.explanation or "",
                    "ec_numbers": reaction_record.get("ec_numbers") or [],
                    "annotation_bucket": classify_annotation_bucket(cls, reaction_record),
                }
            )

    return unknown_positive_index


def generate_class_page(
    class_name: str,
    metrics: dict,
    detailed_df: pd.DataFrame,
    rhea_reactions: dict[str, dict],
    classifier: ReactionClassifier,
    unknown_positive_index: dict[str, list[dict[str, Any]]],
) -> str:
    """Generate HTML content for a single reaction class page."""
    import html as html_module

    # Get class object and metadata
    cls = classifier.reaction_classes.get(class_name)
    go_id = getattr(cls, "GO_ID", "N/A") if cls else "N/A"
    ec_prefix = getattr(cls, "EC_NUMBER_PREFIX", None) if cls else None
    ec_broad_xrefs = getattr(cls, "EC_BROAD_XREFS", []) if cls else []

    # Format GO and EC with links
    go_html = bioregistry_link(go_id) if go_id and go_id != "N/A" else "N/A"
    ec_formatted = format_ec_with_dashes(ec_prefix) if ec_prefix else "N/A"
    ec_html = bioregistry_link(ec_formatted) if ec_prefix else "N/A"
    ec_broad_html = (
        ", ".join(
            bioregistry_link(format_ec_with_dashes(ec_value))
            for ec_value in ec_broad_xrefs
        )
        if ec_broad_xrefs
        else "N/A"
    )
    ec_prefix_metadata_html = ""
    if ec_prefix:
        ec_prefix_metadata_html = f'''
            <div class="metadata-item">
                <div class="metadata-label">Declared EC Scope</div>
                <div class="metadata-value">{ec_html}</div>
            </div>
        '''

    broad_xref_metadata_html = ""
    if ec_broad_xrefs:
        broad_xref_metadata_html = f'''
            <div class="metadata-item">
                <div class="metadata-label">GO EC Xref Type</div>
                <div class="metadata-value">broad match</div>
            </div>
            <div class="metadata-item">
                <div class="metadata-label">GO EC Broad Xrefs</div>
                <div class="metadata-value">{ec_broad_html}</div>
            </div>
        '''

    # Get class docstring
    class_docstring = get_class_docstring(cls) if cls else ""
    docstring_html = ""
    if class_docstring:
        # Escape HTML and convert newlines to <br>
        escaped_doc = html_module.escape(class_docstring)
        # Convert double newlines to paragraph breaks
        escaped_doc = escaped_doc.replace("\n\n", "</p><p>")
        escaped_doc = escaped_doc.replace("\n", "<br>")
        docstring_html = f'''
        <div class="definition">
            <h2>Definition</h2>
            <div class="definition-content">
                <p>{escaped_doc}</p>
            </div>
        </div>
        '''

    # Extract source code
    source_code = ""
    class_context_code = ""
    file_path = ""
    context_start_line = 0
    start_line = 0
    aggregate_context = get_aggregate_membership_context(cls) if cls else None
    if cls:
        class_context_code, file_path, context_start_line = extract_class_context_source(
            cls, "check_membership_impl"
        )
        source_code, file_path, start_line = extract_method_source(cls, "check_membership_impl")
        # Make file path relative
        if file_path:
            file_path = Path(file_path).name

    context_html = ""
    if class_context_code:
        context_html = f'''
        <div class="code-section">
            <div class="code-header">
                <span class="file-path">{file_path}:{context_start_line}</span>
            </div>
            <div class="section-note">
                Class-level declarative context used by the classifier, including identifiers, evidence requirements, and pattern definitions.
            </div>
            <pre><code>{simple_syntax_highlight(class_context_code)}</code></pre>
        </div>
        '''

    aggregate_context_html = ""
    if aggregate_context:
        child_items = ""
        for child in aggregate_context["children"]:
            child_ec_html = (
                bioregistry_link(child["ec"]) if child["ec"] else "No declared EC scope"
            )
            child_go_html = (
                bioregistry_link(child["go_id"])
                if child["go_id"]
                else "No declared GO term"
            )
            child_items += f"""
            <li>
                <strong>{html_module.escape(child['name'])}</strong>
                <span>({child_ec_html}; {child_go_html})</span>
            </li>
            """

        aggregate_context_html = f'''
        <div class="code-section">
            <div class="code-header">
                <span class="file-path">{html_module.escape(aggregate_context["title"])}</span>
            </div>
            <div class="section-note">
                {html_module.escape(aggregate_context["summary"])}
            </div>
            <div class="metadata-grid" style="margin-top: 1rem;">
                <div class="metadata-item">
                    <div class="metadata-label">{html_module.escape(aggregate_context["scope_label"])}</div>
                    <div class="metadata-value">{html_module.escape(aggregate_context["scope_value"])}</div>
                </div>
                <div class="metadata-item">
                    <div class="metadata-label">Curated child classifiers</div>
                    <div class="metadata-value">{len(aggregate_context["children"])}</div>
                </div>
            </div>
            <div class="section-note" style="margin-top: 1rem;">
                The wrapper returns positive only if one of these child classifiers matches the reaction.
            </div>
            <ul>
                {child_items}
            </ul>
        </div>
        '''

    # Get detailed predictions for this class
    class_details = detailed_df[detailed_df["class"] == class_name] if len(detailed_df) > 0 else pd.DataFrame()

    # Build sections for TP/FP/TN/FN
    def build_reaction_section(outcome: str, title: str, css_class: str) -> str:
        items = class_details[class_details["outcome"] == outcome] if len(class_details) > 0 else pd.DataFrame()
        if len(items) == 0:
            return ""

        total_count = len(items)
        truncated = total_count > MAX_ITEMS_PER_SECTION
        display_items = items.head(MAX_ITEMS_PER_SECTION)

        reactions_html = ""
        for _, row in display_items.iterrows():
            rhea_id = row.get("rhea_id", "")
            label = row.get("label", "") or rhea_reactions.get(rhea_id, {}).get("label", "No label")
            explanation = row.get("explanation", "")
            go_term = row.get("go_term", "")

            # Create linked RHEA ID
            rhea_link = format_rhea_id(rhea_id)

            # Create linked GO term
            go_link = bioregistry_link(go_term) if go_term else ""
            go_html = f' <span class="go-tag">[{go_link}]</span>' if go_link else ""

            reactions_html += f'''
            <div class="reaction-item">
                <span class="reaction-id">{rhea_link}</span>{go_html}:
                <span class="reaction-label">{html_module.escape(str(label))}</span>
                {f'<div class="reaction-explanation">{html_module.escape(str(explanation))}</div>' if explanation else ''}
            </div>
            '''

        # Add truncation notice if needed
        if truncated:
            reactions_html += f'''
            <div class="truncation-notice">
                ... and {total_count - MAX_ITEMS_PER_SECTION} more (showing {MAX_ITEMS_PER_SECTION} of {total_count})
            </div>
            '''

        return f'''
        <div class="collapsible">
            <div class="collapsible-header {css_class}-header">
                <span class="collapsible-title">
                    <span class="collapsible-icon">▶</span> {title}
                </span>
                <span class="collapsible-count">{total_count}</span>
            </div>
            <div class="collapsible-content">
                {reactions_html}
            </div>
        </div>
        '''

    def build_unknown_positive_section() -> str:
        items = unknown_positive_index.get(class_name, [])
        if not items:
            return ""

        bucket_order = [
            (
                "unannotated",
                "Truly unannotated candidate positives",
                "No GO mapping and no EC annotation are present for these reactions.",
                "fn-header",
            ),
            (
                "ec_aligned",
                "EC-backed, GO-missing positives supporting this class",
                "These reactions lack direct GO support in the benchmark, but their EC annotations already support this class.",
                "tp-header",
            ),
            (
                "ec_conflict",
                "EC-backed, GO-missing positives conflicting with this class",
                "These reactions lack direct GO support for this class and already carry a different EC assignment. They are useful false-positive review candidates.",
                "fp-header",
            ),
        ]

        def render_ec_numbers(ec_numbers: list[str]) -> str:
            if not ec_numbers:
                return ""
            links = ", ".join(
                bioregistry_link(format_ec_with_dashes(ec_number))
                for ec_number in ec_numbers
            )
            return f'<div class="reaction-meta">EC annotation: {links}</div>'

        sections_html = ""
        for bucket, title, note, css_class in bucket_order:
            bucket_items = [
                item for item in items if item.get("annotation_bucket") == bucket
            ]
            if not bucket_items:
                continue

            total_count = len(bucket_items)
            truncated = total_count > MAX_ITEMS_PER_SECTION
            display_items = bucket_items[:MAX_ITEMS_PER_SECTION]

            reactions_html = ""
            for item in display_items:
                rhea_link = format_rhea_id(item["rhea_id"])
                reactions_html += f'''
                <div class="reaction-item">
                    <span class="reaction-id">{rhea_link}</span>:
                    <span class="reaction-label">{html_module.escape(item["label"])}</span>
                    {render_ec_numbers(item.get("ec_numbers") or [])}
                    {f'<div class="reaction-explanation">{html_module.escape(item["explanation"])}</div>' if item["explanation"] else ''}
                </div>
                '''
            if truncated:
                reactions_html += f'''
                <div class="truncation-notice">
                    ... and {total_count - MAX_ITEMS_PER_SECTION} more (showing {MAX_ITEMS_PER_SECTION} of {total_count})
                </div>
                '''

            sections_html += f'''
            <div class="collapsible">
                <div class="collapsible-header {css_class}">
                    <span class="collapsible-title">
                        <span class="collapsible-icon">▶</span> {title}
                    </span>
                    <span class="collapsible-count">{total_count}</span>
                </div>
                <div class="collapsible-content">
                    <div class="section-note">{note}</div>
                    {reactions_html}
                </div>
            </div>
            '''

        return f'''
        <h2>Positives Without Direct GO Support</h2>
        <div class="section-note">
            Positive rule matches on reactions without direct GO support for this class. These matches are partitioned into truly unannotated candidates, EC-backed reactions already consistent with this class, and EC-backed reactions that conflict with this class.
        </div>
        {sections_html}
        '''

    # Build the page content
    content = f'''
        <nav class="nav">
            <a href="index.html">← Back to Index</a>
        </nav>

        <h1>{class_name}</h1>

        <div class="metadata">
            <div class="metadata-item">
                <div class="metadata-label">GO Term</div>
                <div class="metadata-value">{go_html}</div>
            </div>
            {ec_prefix_metadata_html}
            {broad_xref_metadata_html}
            <div class="metadata-item">
                <div class="metadata-label">Total Evaluated</div>
                <div class="metadata-value">{metrics.get("total", 0)}</div>
            </div>
        </div>

        {docstring_html}

        <h2>Performance Metrics</h2>
        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-value {get_stat_class(metrics.get('f1_score', 0))}">{metrics.get('f1_score', 0):.3f}</div>
                <div class="stat-label">F1 Score</div>
            </div>
            <div class="stat-card">
                <div class="stat-value {get_stat_class(metrics.get('precision', 0))}">{metrics.get('precision', 0):.3f}</div>
                <div class="stat-label">Precision</div>
            </div>
            <div class="stat-card">
                <div class="stat-value {get_stat_class(metrics.get('recall', 0))}">{metrics.get('recall', 0):.3f}</div>
                <div class="stat-label">Recall</div>
            </div>
            <div class="stat-card">
                <div class="stat-value {get_stat_class(metrics.get('accuracy', 0))}">{metrics.get('accuracy', 0):.3f}</div>
                <div class="stat-label">Accuracy</div>
            </div>
            <div class="stat-card">
                <div class="stat-value {get_stat_class(metrics.get('mcc', 0), (0.4, 0.7))}">{metrics.get('mcc', 0):.3f}</div>
                <div class="stat-label">MCC</div>
            </div>
            <div class="stat-card">
                <div class="stat-value {get_stat_class(metrics.get('specificity', 0))}">{metrics.get('specificity', 0):.3f}</div>
                <div class="stat-label">Specificity</div>
            </div>
        </div>

        <h2>Confusion Matrix</h2>
        <div class="confusion-matrix">
            <div class="cm-cell tp">
                <div class="cm-count">{metrics.get('tp', 0)}</div>
                <div class="cm-label">True Positives</div>
            </div>
            <div class="cm-cell fp">
                <div class="cm-count">{metrics.get('fp', 0)}</div>
                <div class="cm-label">False Positives</div>
            </div>
            <div class="cm-cell fn">
                <div class="cm-count">{metrics.get('fn', 0)}</div>
                <div class="cm-label">False Negatives</div>
            </div>
            <div class="cm-cell tn">
                <div class="cm-count">{metrics.get('tn', 0)}</div>
                <div class="cm-label">True Negatives</div>
            </div>
        </div>

        <h2>Classification Code</h2>
        {context_html}
        {aggregate_context_html}
        <div class="code-section">
            <div class="code-header">
                <span class="file-path">{file_path}:{start_line}</span>
            </div>
            <div class="section-note">
                Imperative membership logic. This often calls the declarative class attributes shown above.
            </div>
            <pre><code>{simple_syntax_highlight(source_code)}</code></pre>
        </div>

        <h2>Detailed Results</h2>
        {build_reaction_section("TP", "True Positives", "tp")}
        {build_reaction_section("FP", "False Positives", "fp")}
        {build_reaction_section("FN", "False Negatives", "fn")}
        {build_reaction_section("TN", "True Negatives", "tn")}
        {build_unknown_positive_section()}

        <div class="footer">
            Generated by autarch export utility
        </div>
    '''

    return HTML_TEMPLATE.format(title=f"{class_name} - Autarch", content=content)


def generate_index_page(
    metrics_df: pd.DataFrame,
    has_complexity_plot: bool = False,
    has_embedding_page: bool = False,
    has_ec_hierarchy_page: bool = False,
) -> str:
    """Generate HTML content for the index page."""

    # Sort by F1 score descending
    sorted_df = metrics_df.sort_values("f1_score", ascending=False)

    # Calculate summary stats
    avg_f1 = metrics_df["f1_score"].mean()
    avg_precision = metrics_df["precision"].mean()
    avg_recall = metrics_df["recall"].mean()
    total_classes = len(metrics_df)

    # Build class cards
    class_cards = ""
    for _, row in sorted_df.iterrows():
        class_name = row["class"]
        f1 = row["f1_score"]
        precision = row["precision"]
        recall = row["recall"]

        class_cards += f'''
        <a href="{class_name.lower()}.html" class="class-card">
            <div class="class-name">{class_name}</div>
            <div class="class-stats">
                <div class="class-stat">
                    <span class="class-stat-label">F1:</span>
                    <span class="{get_stat_class(f1)}">{f1:.3f}</span>
                </div>
                <div class="class-stat">
                    <span class="class-stat-label">P:</span>
                    <span class="{get_stat_class(precision)}">{precision:.3f}</span>
                </div>
                <div class="class-stat">
                    <span class="class-stat-label">R:</span>
                    <span class="{get_stat_class(recall)}">{recall:.3f}</span>
                </div>
            </div>
        </a>
        '''

    viz_cards = ""
    if has_complexity_plot:
        viz_cards += '''
        <a href="complexity_vs_performance.html" class="viz-card">
            <div class="viz-card-title">Complexity versus F1</div>
            <div class="viz-card-text">Interactive scatter plot of classifier F1, complexity, support, and conservative confidence.</div>
        </a>
        '''
    if has_embedding_page:
        viz_cards += '''
        <a href="rhea_browser.html" class="viz-card">
            <div class="viz-card-title">RHEA Browser</div>
            <div class="viz-card-text">Browse structured RHEA reactions, inspect participant slots, and navigate a linked text embedding view built from reaction definitions and participant descriptors.</div>
        </a>
        '''
    if has_ec_hierarchy_page:
        viz_cards += '''
        <a href="ec_hierarchy.html" class="viz-card">
            <div class="viz-card-title">EC Hierarchy</div>
            <div class="viz-card-text">One-page EC level 1-3 hierarchy with GO mappings, subsumed RHEA counts, and Python classifier evaluation plus complexity overlays.</div>
        </a>
        '''

    interactive_views_html = ""
    if viz_cards:
        interactive_views_html += f'''
        <h2>Interactive Views</h2>
        <div class="viz-grid">
            {viz_cards}
        </div>
        '''
    if has_complexity_plot:
        interactive_views_html += '''
        <div class="plot-section">
            <div class="section-header">
                <h2>Complexity versus F1</h2>
                <a href="complexity_vs_performance.html">Open full view</a>
            </div>
            <div class="section-note">
                Points are sized by positive support and colored by conservative F1, which penalizes low-support optimism.
            </div>
            <iframe src="complexity_vs_performance.html" class="plot-frame" loading="lazy"></iframe>
        </div>
        '''

    content = f'''
        <h1>Autarch Reaction Classification Report</h1>

        <div class="metadata">
            <div class="metadata-item">
                <div class="metadata-label">Total Classes</div>
                <div class="metadata-value">{total_classes}</div>
            </div>
            <div class="metadata-item">
                <div class="metadata-label">Avg F1</div>
                <div class="metadata-value {get_stat_class(avg_f1)}">{avg_f1:.3f}</div>
            </div>
            <div class="metadata-item">
                <div class="metadata-label">Avg Precision</div>
                <div class="metadata-value {get_stat_class(avg_precision)}">{avg_precision:.3f}</div>
            </div>
            <div class="metadata-item">
                <div class="metadata-label">Avg Recall</div>
                <div class="metadata-value {get_stat_class(avg_recall)}">{avg_recall:.3f}</div>
            </div>
        </div>

        {interactive_views_html}

        <h2>Reaction Classes (sorted by F1 score)</h2>
        <div class="class-grid">
            {class_cards}
        </div>

        <div class="footer">
            Generated by autarch export utility
        </div>
    '''

    return HTML_TEMPLATE.format(title="Autarch Classification Report", content=content)


def export_html(
    results_dir: Path,
    output_dir: Optional[Path] = None,
    verbose: bool = False,
) -> Path:
    """Export evaluation results as standalone HTML pages.

    Args:
        results_dir: Directory containing evaluation_results.csv and detailed_predictions.csv
        output_dir: Directory to write HTML files (default: results_dir/html)
        verbose: Print progress messages

    Returns:
        Path to the output directory
    """
    results_dir = Path(results_dir)
    if output_dir is None:
        output_dir = results_dir / "html"
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for html_file in output_dir.glob("*.html"):
        html_file.unlink()

    if verbose:
        print(f"Loading evaluation data from {results_dir}...")

    # Load data
    metrics_df, detailed_df, rhea_reactions = load_evaluation_data(results_dir)
    complexity_report = analyze_all_classifiers()

    # Generate interactive report-level visualizations
    if verbose:
        print("Generating complexity/performance scatter...")
    create_complexity_vs_performance(
        metrics_df,
        complexity_report,
        output_dir,
        class_page_links=True,
    )

    has_embedding_page = False
    try:
        if verbose:
            print("Generating RHEA reaction browser...")
        save_rhea_browser_html(
            output_dir / "rhea_browser.html",
            results_dir=results_dir,
        )
        has_embedding_page = True
    except FileNotFoundError:
        has_embedding_page = False

    has_ec_hierarchy_page = False
    try:
        if verbose:
            print("Generating EC hierarchy report...")
        complexity_df = pd.DataFrame([classifier.to_dict() for classifier in complexity_report.classifiers])
        save_ec_hierarchy_report_html(
            output_dir / "ec_hierarchy.html",
            metrics_df,
            complexity_df,
        )
        has_ec_hierarchy_page = True
    except FileNotFoundError:
        has_ec_hierarchy_page = False

    # Initialize classifier for source code extraction
    classifier = ReactionClassifier()
    unknown_positive_index = build_go_missing_positive_index(
        rhea_reactions,
        classifier,
    )

    # Generate index page
    if verbose:
        print("Generating index page...")
    index_html = generate_index_page(
        metrics_df,
        has_complexity_plot=True,
        has_embedding_page=has_embedding_page,
        has_ec_hierarchy_page=has_ec_hierarchy_page,
    )
    (output_dir / "index.html").write_text(index_html)

    # Generate page for each class
    for _, row in metrics_df.iterrows():
        class_name = row["class"]
        if verbose:
            print(f"Generating page for {class_name}...")

        metrics = row.to_dict()
        class_html = generate_class_page(
            class_name,
            metrics,
            detailed_df,
            rhea_reactions,
            classifier,
            unknown_positive_index,
        )
        (output_dir / f"{class_name.lower()}.html").write_text(class_html)

    if verbose:
        print(f"Exported {len(metrics_df)} class pages to {output_dir}")

    return output_dir
