"""Text embeddings and browser views for RHEA reactions.

This is a deterministic lexical baseline intended for exploration. It uses
hashed TF-IDF features over reaction descriptions assembled from the curated
reaction equation plus participant metadata, then projects them to 2D with PCA.
"""

from __future__ import annotations

import ast
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Iterable, cast

import numpy as np
import pandas as pd
import plotly.express as px  # type: ignore[import-untyped]


DEFAULT_N_FEATURES = 512
EC_MAJOR_NAMES = {
    "1": "Oxidoreductases",
    "2": "Transferases",
    "3": "Hydrolases",
    "4": "Lyases",
    "5": "Isomerases",
    "6": "Ligases",
    "7": "Translocases",
}
PARTICIPANT_BUCKET_ORDER = ["1-2", "3-4", "5-6", "7+"]
GO_FACET_EXCLUDE_IDS = {"GO:0003674", "GO:0003824"}
GO_FACET_LIMIT = 18
RULE_STATUS_LABELS = {
    "agreement": "Consistent",
    "asserted_only": "Asserted only",
    "inferred_only": "Agent only",
    "unknown_positive": "Unannotated candidate",
    "ec_backed_positive": "EC-backed, GO-missing",
    "mixed": "Mismatch",
    "none": "No rule coverage",
}


def ec_major_label(ec_major: str) -> str:
    """Render an EC major class key into a compact label."""
    name = EC_MAJOR_NAMES.get(ec_major)
    if name is None:
        return ec_major
    return f"{ec_major} {name}"


def ec_sort_key(ec_fragment: str) -> tuple[int, ...]:
    """Sort EC fragments numerically where possible."""
    return tuple(int(part) for part in ec_fragment.split("."))


def extract_ec_hierarchy(ec_numbers: list[str]) -> tuple[list[str], list[str]]:
    """Extract EC major classes and subclasses from mapped EC numbers.

    Examples:
        >>> extract_ec_hierarchy(["2.7.1.1", "1.14.13.39", "2.-.-.-"])
        (['1', '2'], ['1.14', '2.7'])
    """
    majors: set[str] = set()
    subclasses: set[str] = set()
    for ec_number in ec_numbers:
        parts = ec_number.split(".")
        if len(parts) >= 1 and parts[0] and parts[0] != "-":
            majors.add(parts[0])
        if (
            len(parts) >= 2
            and parts[0]
            and parts[1]
            and parts[0] != "-"
            and parts[1] != "-"
        ):
            subclasses.add(f"{parts[0]}.{parts[1]}")
    return sorted(majors, key=ec_sort_key), sorted(subclasses, key=ec_sort_key)


def has_polymer_context(participants: list[dict[str, Any]]) -> bool:
    """Return True if any participant carries polymer annotations."""
    return any(
        participant.get("polymer_type")
        or participant.get("polymer_index")
        or participant.get("monomer")
        for participant in participants
    )


def has_location_context(participants: list[dict[str, Any]]) -> bool:
    """Return True if any participant carries a location annotation."""
    return any(participant.get("location") for participant in participants)


def classify_participant_bucket(participant_total: int) -> str:
    """Bucket reactions by total participant count."""
    if participant_total <= 2:
        return "1-2"
    if participant_total <= 4:
        return "3-4"
    if participant_total <= 6:
        return "5-6"
    return "7+"


def dedupe_preserve_order(values: Iterable[str]) -> list[str]:
    """Deduplicate a sequence while preserving first-seen order."""
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return ordered


def load_go_term_lookup(cache_dir: Path = Path("cache")) -> dict[str, dict[str, Any]]:
    """Load GO labels and ancestor links from the local cache."""
    cache_path = Path(cache_dir) / "go_terms.jsonl"
    if not cache_path.exists():
        return {}

    lookup: dict[str, dict[str, Any]] = {}
    with open(cache_path) as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            go_id = row["go_id"]
            lookup[go_id] = {
                "label": row.get("label", go_id),
                "ancestors": dedupe_preserve_order(row.get("ancestors", [])),
            }
    return lookup


def build_go_annotations(
    go_terms: list[str],
    go_lookup: dict[str, dict[str, Any]],
) -> tuple[list[dict[str, str]], list[str], list[str]]:
    """Build direct GO entries plus closure ids/labels for browser use."""
    direct_ids = dedupe_preserve_order(go_terms)
    direct_entries = [
        {"id": go_id, "label": go_lookup.get(go_id, {}).get("label", go_id)}
        for go_id in direct_ids
    ]
    closure_ids: list[str] = []
    seen: set[str] = set()

    def visit(go_id: str) -> None:
        if go_id in seen:
            return
        seen.add(go_id)
        closure_ids.append(go_id)
        for ancestor_id in go_lookup.get(go_id, {}).get("ancestors", []):
            visit(ancestor_id)

    for go_id in direct_ids:
        visit(go_id)
    closure_labels = [
        go_lookup.get(go_id, {}).get("label", go_id) for go_id in closure_ids
    ]
    return direct_entries, closure_ids, closure_labels


def load_ec_label_lookup(ec_numbers: Iterable[str]) -> dict[str, str]:
    """Load EC labels from the local OBO SQLite adapter."""
    unique_ec_numbers = sorted(set(ec_numbers))
    if not unique_ec_numbers:
        return {}

    try:
        from oaklib import get_adapter  # type: ignore[import-untyped]

        adapter = get_adapter("sqlite:obo:ec")
        return {
            ec_curie.removeprefix("EC:"): label
            for ec_curie, label in adapter.labels(
                [f"EC:{ec_number}" for ec_number in unique_ec_numbers]
            )
            if label
        }
    except Exception:
        return {}


def build_ec_number_entries(
    ec_numbers: list[str],
    ec_label_lookup: dict[str, str],
) -> list[dict[str, str]]:
    """Build EC entries with labels when available."""
    return [
        {"id": ec_number, "label": ec_label_lookup.get(ec_number, "")}
        for ec_number in dedupe_preserve_order(ec_numbers)
    ]


def tokenize_label_text(text: str) -> list[str]:
    """Tokenize reaction text into lexical features.

    Includes word-level tokens and character trigrams to retain signal from
    biochemical names and reaction equations.

    Examples:
        >>> toks = tokenize_label_text("ATP + H2O = ADP + phosphate")
        >>> "w:atp" in toks
        True
        >>> any(tok.startswith("c:at") for tok in toks)
        True
    """
    normalized = text.lower().strip()
    word_tokens = re.findall(r"[a-z0-9][a-z0-9+()/-]*", normalized)
    features = [f"w:{token}" for token in word_tokens]

    compact = re.sub(r"\s+", " ", normalized)
    for i in range(max(0, len(compact) - 2)):
        trigram = compact[i : i + 3]
        if trigram.strip():
            features.append(f"c:{trigram}")
    return features


def stable_feature_index(feature: str, n_features: int) -> int:
    """Map a feature string to a stable hashed index."""
    digest = hashlib.blake2b(feature.encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(digest, "little") % n_features


def build_text_feature_matrix(
    texts: Iterable[str],
    n_features: int = DEFAULT_N_FEATURES,
) -> np.ndarray:
    """Build a deterministic hashed TF-IDF matrix from reaction texts."""
    text_list = list(texts)
    matrix = np.zeros((len(text_list), n_features), dtype=np.float32)
    document_frequency = np.zeros(n_features, dtype=np.int32)

    for row_index, text in enumerate(text_list):
        token_counts: dict[int, float] = {}
        seen_indices: set[int] = set()
        for feature in tokenize_label_text(text):
            index = stable_feature_index(feature, n_features)
            weight = 0.35 if feature.startswith("c:") else 1.0
            token_counts[index] = token_counts.get(index, 0.0) + weight
            seen_indices.add(index)

        for index, count in token_counts.items():
            matrix[row_index, index] = count
        for index in seen_indices:
            document_frequency[index] += 1

    if len(text_list) == 0:
        return matrix

    idf = np.log((1.0 + len(text_list)) / (1.0 + document_frequency)) + 1.0
    matrix *= idf

    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    norms[norms == 0.0] = 1.0
    matrix /= norms
    return matrix


def project_embedding_matrix(matrix: np.ndarray, n_components: int = 2) -> np.ndarray:
    """Project a feature matrix to low-dimensional coordinates with PCA."""
    if matrix.size == 0:
        return np.zeros((0, n_components), dtype=np.float32)

    centered = matrix - matrix.mean(axis=0, keepdims=True)
    covariance = centered.T @ centered
    eigvals, eigvecs = np.linalg.eigh(covariance)
    order = np.argsort(eigvals)[::-1][:n_components]
    components = eigvecs[:, order]
    coords = centered @ components

    if coords.shape[1] < n_components:
        pad = np.zeros(
            (coords.shape[0], n_components - coords.shape[1]),
            dtype=coords.dtype,
        )
        coords = np.hstack([coords, pad])
    return coords.astype(np.float32, copy=False)


def classify_annotation_group(go_terms: list[str], ec_numbers: list[str]) -> str:
    """Classify reactions by available annotation source."""
    if go_terms and ec_numbers:
        return "GO+EC"
    if go_terms:
        return "GO only"
    if ec_numbers:
        return "EC only"
    return "Unannotated"


def participant_browser_record(participant: dict[str, Any]) -> dict[str, Any]:
    """Trim a participant record to the fields needed in the browser."""
    return {
        "name": participant.get("name"),
        "chebi_id": participant.get("chebi_id"),
        "count": participant.get("count"),
        "stoichiometry": participant.get("stoichiometry"),
        "polymer_index": participant.get("polymer_index"),
        "polymer_type": participant.get("polymer_type"),
        "location": participant.get("location"),
        "monomer": participant.get("monomer"),
    }


def participant_display_text(participant: dict[str, Any]) -> str:
    """Render a participant into a compact textual descriptor.

    Examples:
        >>> participant_display_text({"name": "ATP", "chebi_id": "CHEBI:30616", "count": 2})
        '2 ATP [CHEBI:30616]'
        >>> participant_display_text({"name": "RNA", "polymer_index": "n+1", "polymer_type": "rna"})
        'RNA {rna, n+1}'
    """
    name = (
        participant.get("name")
        or participant.get("chebi_id")
        or "unnamed participant"
    )
    prefix = ""
    if participant.get("stoichiometry"):
        prefix = f"{participant['stoichiometry']} "
    elif participant.get("count") and participant.get("count") != 1:
        prefix = f"{participant['count']} "

    pieces = [f"{prefix}{name}"]
    if participant.get("chebi_id"):
        pieces.append(f"[{participant['chebi_id']}]")

    polymer_bits = []
    if participant.get("polymer_type"):
        polymer_bits.append(str(participant["polymer_type"]))
    if participant.get("polymer_index"):
        polymer_bits.append(str(participant["polymer_index"]))
    if polymer_bits:
        pieces.append("{" + ", ".join(polymer_bits) + "}")

    if participant.get("location"):
        pieces.append(f"@{participant['location']}")
    return " ".join(pieces)


def build_reaction_embedding_text(
    label: str,
    left_participants: list[dict[str, Any]],
    right_participants: list[dict[str, Any]],
) -> str:
    """Build embedding text from reaction definition and participants.

    Examples:
        >>> text = build_reaction_embedding_text(
        ...     "ATP + H2O = ADP + phosphate",
        ...     [{"name": "ATP"}, {"name": "H2O"}],
        ...     [{"name": "ADP"}, {"name": "phosphate"}],
        ... )
        >>> "reactants" in text and "products" in text
        True
    """
    left_text = "; ".join(participant_display_text(p) for p in left_participants)
    right_text = "; ".join(participant_display_text(p) for p in right_participants)
    parts = [label.strip()]
    if left_text:
        parts.append(f"reactants: {left_text}")
    if right_text:
        parts.append(f"products: {right_text}")
    return " | ".join(part for part in parts if part)


def build_reaction_search_text(
    rhea_id: str,
    label: str,
    left_summary: str,
    right_summary: str,
    go_terms: list[str],
    ec_numbers: list[str],
    extra_terms: list[str] | None = None,
) -> str:
    """Build search text for browser filtering."""
    return " ".join(
        [
            rhea_id,
            label,
            left_summary,
            right_summary,
            " ".join(go_terms),
            " ".join(ec_numbers),
            " ".join(extra_terms or []),
        ]
    ).lower()


def load_rhea_text_df(cache_dir: Path = Path("cache")) -> pd.DataFrame:
    """Load RHEA reactions plus rich browser/embedding metadata from cache."""
    cache_path = Path(cache_dir) / "rhea_reactions.jsonl"
    if not cache_path.exists():
        raise FileNotFoundError(f"RHEA cache not found at {cache_path}")

    go_lookup = load_go_term_lookup(cache_dir)
    raw_rows = []
    with open(cache_path) as handle:
        for line in handle:
            if not line.strip():
                continue
            raw_rows.append(json.loads(line))
    ec_label_lookup = load_ec_label_lookup(
        ec_number
        for row in raw_rows
        for ec_number in row.get("ec_numbers", [])
    )
    rows = []
    for row in raw_rows:
        go_terms = row.get("go_terms", [])
        ec_numbers = row.get("ec_numbers", [])
        go_term_entries, go_closure_ids, go_closure_labels = build_go_annotations(
            go_terms,
            go_lookup,
        )
        ec_number_entries = build_ec_number_entries(ec_numbers, ec_label_lookup)
        reaction = row.get("reaction") or {}
        left_raw = reaction.get("left_participants") or []
        right_raw = reaction.get("right_participants") or []
        left_participants = [participant_browser_record(p) for p in left_raw]
        right_participants = [participant_browser_record(p) for p in right_raw]
        all_participants = left_participants + right_participants
        left_summary = "; ".join(
            participant_display_text(p) for p in left_participants
        )
        right_summary = "; ".join(
            participant_display_text(p) for p in right_participants
        )
        label = row.get("label", "") or row["rhea_id"]
        ec_major_classes, ec_subclasses = extract_ec_hierarchy(ec_numbers)
        primary_ec_major = ec_major_classes[0] if ec_major_classes else "Unclassified"
        participant_total = len(all_participants)
        rows.append(
            {
                "rhea_id": row["rhea_id"],
                "label": label,
                "left_participants": left_participants,
                "right_participants": right_participants,
                "left_summary": left_summary,
                "right_summary": right_summary,
                "embedding_text": build_reaction_embedding_text(
                    label,
                    left_participants,
                    right_participants,
                ),
                "search_text": build_reaction_search_text(
                    row["rhea_id"],
                    label,
                    left_summary,
                    right_summary,
                    go_terms,
                    ec_numbers,
                    extra_terms=(
                        go_closure_ids
                        + go_closure_labels
                        + [entry["label"] for entry in go_term_entries]
                        + [entry["label"] for entry in ec_number_entries if entry["label"]]
                        + [ec_major_label(ec_major) for ec_major in ec_major_classes]
                    ),
                ),
                "go_terms": go_terms,
                "go_term_entries": go_term_entries,
                "go_closure_ids": go_closure_ids,
                "go_closure_labels": go_closure_labels,
                "ec_numbers": ec_numbers,
                "ec_number_entries": ec_number_entries,
                "ec_major_classes": ec_major_classes,
                "ec_major_labels": [
                    ec_major_label(ec_major) for ec_major in ec_major_classes
                ],
                "ec_subclasses": ec_subclasses,
                "primary_ec_major": primary_ec_major,
                "primary_ec_major_label": (
                    ec_major_label(primary_ec_major)
                    if primary_ec_major != "Unclassified"
                    else "Unclassified"
                ),
                "go_count": len(go_terms),
                "ec_count": len(ec_numbers),
                "participant_total": participant_total,
                "participant_bucket": classify_participant_bucket(participant_total),
                "has_polymer_context": has_polymer_context(all_participants),
                "has_location_context": has_location_context(all_participants),
                "has_multiple_ec": len(ec_numbers) > 1,
                "annotation_group": classify_annotation_group(go_terms, ec_numbers),
                "url": f"https://bioregistry.io/{row['rhea_id']}",
            }
        )
    return pd.DataFrame(rows)


def build_rhea_text_embedding_df(
    cache_dir: Path = Path("cache"),
    n_features: int = DEFAULT_N_FEATURES,
) -> pd.DataFrame:
    """Build a dataframe with 2D coordinates for RHEA text embeddings."""
    df = load_rhea_text_df(cache_dir)
    matrix = build_text_feature_matrix(
        df["embedding_text"].tolist(),
        n_features=n_features,
    )
    coords = project_embedding_matrix(matrix, n_components=2)
    df = df.copy()
    df["x"] = coords[:, 0]
    df["y"] = coords[:, 1]
    return df


def literal_string_value(node: ast.AST) -> str | None:
    """Extract a string literal value from an AST node."""
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    return None


def literal_string_list(node: ast.AST) -> list[str]:
    """Extract a list of string literals from an AST node."""
    if not isinstance(node, (ast.List, ast.Tuple, ast.Set)):
        return []
    values: list[str] = []
    for element in node.elts:
        value = literal_string_value(element)
        if value is not None:
            values.append(value)
    return values


def load_rule_class_metadata() -> dict[str, dict[str, Any]]:
    """Load GO/EC metadata for Python reaction classes from ontology source files."""
    ontology_dir = Path(__file__).resolve().parent / "ontology"
    metadata_by_class: dict[str, dict[str, Any]] = {}
    for source_path in sorted(ontology_dir.glob("*.py")):
        if source_path.name in {"__init__.py", "reaction.py"}:
            continue
        module = ast.parse(source_path.read_text(), filename=str(source_path))
        for node in module.body:
            if not isinstance(node, ast.ClassDef):
                continue
            go_id: str | None = None
            ec_prefix: str | None = None
            ec_numbers: list[str] = []
            for statement in node.body:
                if not isinstance(statement, (ast.Assign, ast.AnnAssign)):
                    continue
                target_names: list[str] = []
                value_node: ast.AST | None = None
                if isinstance(statement, ast.Assign):
                    target_names = [
                        target.id
                        for target in statement.targets
                        if isinstance(target, ast.Name)
                    ]
                    value_node = statement.value
                else:
                    if isinstance(statement.target, ast.Name):
                        target_names = [statement.target.id]
                    value_node = statement.value

                if not value_node:
                    continue
                if "GO_ID" in target_names:
                    go_id = literal_string_value(value_node)
                if "EC_NUMBER_PREFIX" in target_names:
                    ec_prefix = literal_string_value(value_node)
                if "EC_NUMBERS" in target_names:
                    ec_numbers = literal_string_list(value_node)

            if go_id or ec_prefix or ec_numbers:
                metadata_by_class[node.name] = {
                    "go_id": go_id,
                    "ec_prefix": ec_prefix,
                    "ec_numbers": ec_numbers,
                }
    return metadata_by_class


def build_asserted_rule_support(
    go_closure_ids: list[str],
    ec_numbers: list[str],
    rule_class_metadata: dict[str, dict[str, Any]],
) -> dict[str, list[str]]:
    """Map Python rule classes to the asserted GO/EC evidence supporting them."""
    from autarch.evaluation import ec_matches_prefix

    support_by_class: dict[str, list[str]] = {}
    for class_name, metadata in rule_class_metadata.items():
        support: list[str] = []
        go_id = metadata["go_id"]
        if go_id and go_id in go_closure_ids:
            support.append("GO")

        matched_ec = any(ec in metadata["ec_numbers"] for ec in ec_numbers)
        if not matched_ec and metadata["ec_prefix"]:
            matched_ec = any(
                ec_matches_prefix(ec_number, metadata["ec_prefix"])
                for ec_number in ec_numbers
            )
        if matched_ec:
            support.append("EC")

        if support:
            support_by_class[class_name] = support
    return support_by_class


def format_asserted_rule_label(class_name: str, support: list[str]) -> str:
    """Render a class name with compact GO/EC evidence tags."""
    return f"{class_name} ({'+'.join(support)})"


def class_report_href(class_name: str) -> str:
    """Return the relative report page href for a classifier."""
    return f"{class_name.lower()}.html"


def classify_rule_status(
    asserted_rule_classes: list[str],
    inferred_rule_classes: list[str],
    has_go_terms: bool,
    has_ec_numbers: bool,
) -> str:
    """Classify the relationship between asserted and inferred classes."""
    asserted_set = set(asserted_rule_classes)
    inferred_set = set(inferred_rule_classes)

    if not asserted_set and not inferred_set:
        return "none"
    if inferred_set and not asserted_set and not has_go_terms:
        if has_ec_numbers:
            return "ec_backed_positive"
        return "unknown_positive"
    if asserted_set == inferred_set:
        return "agreement"
    if asserted_set and not inferred_set:
        return "asserted_only"
    if inferred_set and not asserted_set:
        return "inferred_only"
    return "mixed"


def rule_badge_text(
    rule_status: str,
    missing_rule_classes: list[str],
    extra_rule_classes: list[str],
) -> str:
    """Render a compact status badge for the search results list."""
    if rule_status == "agreement":
        return "consistent"
    if rule_status == "asserted_only":
        return f"asserted only ({len(missing_rule_classes)})"
    if rule_status == "inferred_only":
        return f"agent only ({len(extra_rule_classes)})"
    if rule_status == "unknown_positive":
        return f"unannotated candidate ({len(extra_rule_classes)})"
    if rule_status == "ec_backed_positive":
        return f"ec-backed, GO-missing ({len(extra_rule_classes)})"
    if rule_status == "mixed":
        return "mismatch"
    return "no coverage"


def load_inferred_rule_matches(results_dir: Path) -> dict[str, list[str]]:
    """Load positive rule matches from detailed evaluation predictions."""
    predictions_path = Path(results_dir) / "detailed_predictions.csv"
    if not predictions_path.exists():
        return {}

    prediction_df = pd.read_csv(
        predictions_path,
        usecols=["rhea_id", "class", "prediction"],
    )
    positive_df = prediction_df[
        prediction_df["prediction"] == "positive"
    ].drop_duplicates(subset=["rhea_id", "class"])

    inferred_by_rhea: dict[str, list[str]] = {}
    for row in positive_df.to_dict(orient="records"):
        inferred_by_rhea.setdefault(row["rhea_id"], []).append(row["class"])

    return {
        rhea_id: sorted(class_names)
        for rhea_id, class_names in inferred_by_rhea.items()
    }


def infer_rule_matches_for_go_missing(
    records: list[dict[str, Any]],
) -> dict[str, list[str]]:
    """Directly classify GO-missing reactions for browser-only comparison.

    The benchmark detailed predictions only cover GO-linked reactions. To make
    GO-missing reactions visible in the browser, we classify those reactions
    directly and surface positive matches as either unannotated candidates or
    EC-backed, GO-missing reactions.
    """
    from autarch.classifier import ReactionClassifier
    from autarch.datamodel import Participant, Reaction

    classifier = ReactionClassifier()
    inferred_by_rhea: dict[str, list[str]] = {}
    for record in records:
        if record.get("go_terms"):
            continue
        reaction = Reaction(
            left_participants=[
                Participant(**participant)
                for participant in record["left_participants"]
            ],
            right_participants=[
                Participant(**participant)
                for participant in record["right_participants"]
            ],
            label=record["label"],
        )
        matches = sorted(classifier.get_matching_classes(reaction))
        if matches:
            inferred_by_rhea[record["rhea_id"]] = matches
    return inferred_by_rhea


def build_browser_facet_options(records: list[dict[str, Any]]) -> dict[str, Any]:
    """Build facet option metadata for the browser UI."""
    ec_major_counts: dict[str, int] = {}
    ec_subclass_counts: dict[str, int] = {}
    participant_bucket_counts = {bucket: 0 for bucket in PARTICIPANT_BUCKET_ORDER}
    go_facet_counts: dict[str, int] = {}
    go_facet_labels: dict[str, str] = {}
    rule_status_counts = {
        "covered": 0,
        "agreement": 0,
        "asserted_only": 0,
        "inferred_only": 0,
        "unknown_positive": 0,
        "ec_backed_positive": 0,
        "mixed": 0,
    }

    for record in records:
        for ec_major in record["ec_major_classes"]:
            ec_major_counts[ec_major] = ec_major_counts.get(ec_major, 0) + 1
        for ec_subclass in record["ec_subclasses"]:
            ec_subclass_counts[ec_subclass] = ec_subclass_counts.get(ec_subclass, 0) + 1
        for go_id, go_label in zip(
            record["go_closure_ids"],
            record["go_closure_labels"],
            strict=False,
        ):
            if go_id in GO_FACET_EXCLUDE_IDS:
                continue
            go_facet_counts[go_id] = go_facet_counts.get(go_id, 0) + 1
            go_facet_labels.setdefault(go_id, go_label)
        participant_bucket = record["participant_bucket"]
        participant_bucket_counts[participant_bucket] = (
            participant_bucket_counts.get(participant_bucket, 0) + 1
        )
        if record.get("has_rule_coverage"):
            rule_status_counts["covered"] += 1
        rule_status = record.get("rule_status", "none")
        if rule_status in rule_status_counts:
            rule_status_counts[rule_status] += 1

    ec_majors = [
        {
            "value": ec_major,
            "label": ec_major_label(ec_major),
            "count": ec_major_counts[ec_major],
        }
        for ec_major in sorted(ec_major_counts, key=ec_sort_key)
    ]
    ec_subclasses = [
        {
            "value": ec_subclass,
            "label": ec_subclass,
            "major": ec_subclass.split(".")[0],
            "count": ec_subclass_counts[ec_subclass],
        }
        for ec_subclass in sorted(ec_subclass_counts, key=ec_sort_key)
    ]
    participant_buckets = [
        {
            "value": bucket,
            "label": f"{bucket} participants",
            "count": participant_bucket_counts[bucket],
        }
        for bucket in PARTICIPANT_BUCKET_ORDER
        if participant_bucket_counts.get(bucket, 0) > 0
    ]
    go_facets = [
        {
            "id": go_id,
            "label": go_facet_labels[go_id],
            "count": go_facet_counts[go_id],
        }
        for go_id in sorted(
            go_facet_counts,
            key=lambda go_id: (-go_facet_counts[go_id], go_facet_labels[go_id]),
        )[:GO_FACET_LIMIT]
    ]
    rule_statuses = [
        {
            "value": "covered",
            "label": "Has rule coverage",
            "count": rule_status_counts["covered"],
        },
        {
            "value": "agreement",
            "label": RULE_STATUS_LABELS["agreement"],
            "count": rule_status_counts["agreement"],
        },
        {
            "value": "asserted_only",
            "label": RULE_STATUS_LABELS["asserted_only"],
            "count": rule_status_counts["asserted_only"],
        },
        {
            "value": "inferred_only",
            "label": RULE_STATUS_LABELS["inferred_only"],
            "count": rule_status_counts["inferred_only"],
        },
        {
            "value": "unknown_positive",
            "label": RULE_STATUS_LABELS["unknown_positive"],
            "count": rule_status_counts["unknown_positive"],
        },
        {
            "value": "ec_backed_positive",
            "label": RULE_STATUS_LABELS["ec_backed_positive"],
            "count": rule_status_counts["ec_backed_positive"],
        },
        {
            "value": "mixed",
            "label": RULE_STATUS_LABELS["mixed"],
            "count": rule_status_counts["mixed"],
        },
    ]
    return {
        "ec_majors": ec_majors,
        "ec_subclasses": ec_subclasses,
        "participant_buckets": participant_buckets,
        "go_facets": go_facets,
        "rule_statuses": [
            option
            for option in rule_statuses
            if option["value"] in {"unknown_positive", "ec_backed_positive"}
            or cast(int, option["count"]) > 0
        ],
        "ec_major_labels": {
            ec_major: ec_major_label(ec_major) for ec_major in ec_major_counts
        },
    }


def save_rhea_text_embedding_explorer(
    output_file: Path,
    cache_dir: Path = Path("cache"),
    n_features: int = DEFAULT_N_FEATURES,
) -> Path:
    """Save a simple interactive HTML explorer for RHEA text embeddings."""
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df = build_rhea_text_embedding_df(cache_dir=cache_dir, n_features=n_features)
    fig = px.scatter(
        df,
        x="x",
        y="y",
        color="annotation_group",
        hover_name="rhea_id",
        hover_data={
            "label": True,
            "go_count": True,
            "ec_count": True,
            "url": False,
        },
        custom_data=["url"],
        title="RHEA text embedding explorer",
        labels={"x": "PC1", "y": "PC2"},
        opacity=0.78,
    )
    fig.update_traces(
        marker=dict(size=8, line=dict(width=0.25, color="rgba(0,0,0,0.25)"))
    )
    fig.update_layout(height=820)
    fig.write_html(
        str(output_path),
        include_plotlyjs="cdn",
        full_html=True,
        post_script="""
const plot = document.getElementById('{plot_id}');
plot.on('plotly_click', function(event) {
  const url = event?.points?.[0]?.customdata?.[0];
  if (url) {
    window.open(url, '_blank');
  }
});
""",
    )
    return output_path


def save_rhea_browser_html(
    output_file: Path,
    cache_dir: Path = Path("cache"),
    results_dir: Path | None = None,
    n_features: int = DEFAULT_N_FEATURES,
) -> Path:
    """Save a browser-style RHEA explorer with linked text embeddings."""
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df = build_rhea_text_embedding_df(cache_dir=cache_dir, n_features=n_features)
    row_records = cast(list[dict[str, Any]], df.to_dict(orient="records"))
    rule_class_metadata = load_rule_class_metadata()
    inferred_rule_matches: dict[str, list[str]] = {}
    if results_dir is not None:
        inferred_rule_matches = load_inferred_rule_matches(results_dir)
    go_missing_rule_matches = infer_rule_matches_for_go_missing(row_records)
    for rhea_id, class_names in go_missing_rule_matches.items():
        merged = set(inferred_rule_matches.get(rhea_id, []))
        merged.update(class_names)
        inferred_rule_matches[rhea_id] = sorted(merged)

    records = []
    for row in row_records:
        asserted_rule_support = build_asserted_rule_support(
            row["go_closure_ids"],
            row["ec_numbers"],
            rule_class_metadata,
        )
        asserted_rule_classes = sorted(asserted_rule_support)
        inferred_rule_classes = inferred_rule_matches.get(row["rhea_id"], [])
        missing_rule_classes = sorted(
            set(asserted_rule_classes) - set(inferred_rule_classes)
        )
        extra_rule_classes = sorted(
            set(inferred_rule_classes) - set(asserted_rule_classes)
        )
        rule_status = classify_rule_status(
            asserted_rule_classes,
            inferred_rule_classes,
            has_go_terms=bool(row["go_terms"]),
            has_ec_numbers=bool(row["ec_numbers"]),
        )
        asserted_rule_display = [
            format_asserted_rule_label(class_name, asserted_rule_support[class_name])
            for class_name in asserted_rule_classes
        ]
        rule_search_terms = (
            asserted_rule_display
            + inferred_rule_classes
            + missing_rule_classes
            + extra_rule_classes
            + [RULE_STATUS_LABELS[rule_status]]
        )
        records.append(
            {
                "rhea_id": row["rhea_id"],
                "label": row["label"],
                "left_summary": row["left_summary"],
                "right_summary": row["right_summary"],
                "left_participants": row["left_participants"],
                "right_participants": row["right_participants"],
                "go_terms": row["go_terms"],
                "go_term_entries": row["go_term_entries"],
                "go_closure_ids": row["go_closure_ids"],
                "go_closure_labels": row["go_closure_labels"],
                "ec_numbers": row["ec_numbers"],
                "ec_number_entries": row["ec_number_entries"],
                "ec_major_classes": row["ec_major_classes"],
                "ec_major_labels": row["ec_major_labels"],
                "ec_subclasses": row["ec_subclasses"],
                "primary_ec_major": row["primary_ec_major"],
                "primary_ec_major_label": row["primary_ec_major_label"],
                "go_count": row["go_count"],
                "ec_count": row["ec_count"],
                "participant_total": row["participant_total"],
                "participant_bucket": row["participant_bucket"],
                "has_polymer_context": row["has_polymer_context"],
                "has_location_context": row["has_location_context"],
                "has_multiple_ec": row["has_multiple_ec"],
                "annotation_group": row["annotation_group"],
                "embedding_text": row["embedding_text"],
                "search_text": f"{row['search_text']} {' '.join(rule_search_terms).lower()}",
                "url": row["url"],
                "asserted_rule_classes": asserted_rule_classes,
                "asserted_rule_display": asserted_rule_display,
                "asserted_rule_entries": [
                    {
                        "name": class_name,
                        "display": display,
                        "href": class_report_href(class_name),
                    }
                    for class_name, display in zip(
                        asserted_rule_classes,
                        asserted_rule_display,
                        strict=False,
                    )
                ],
                "inferred_rule_classes": inferred_rule_classes,
                "inferred_rule_entries": [
                    {
                        "name": class_name,
                        "display": class_name,
                        "href": class_report_href(class_name),
                    }
                    for class_name in inferred_rule_classes
                ],
                "missing_rule_classes": missing_rule_classes,
                "extra_rule_classes": extra_rule_classes,
                "rule_status": rule_status,
                "rule_status_label": RULE_STATUS_LABELS[rule_status],
                "rule_badge": rule_badge_text(
                    rule_status,
                    missing_rule_classes,
                    extra_rule_classes,
                ),
                "has_rule_coverage": bool(
                    asserted_rule_classes or inferred_rule_classes
                ),
                "x": round(float(row["x"]), 5),
                "y": round(float(row["y"]), 5),
            }
        )

    stats = {
        "count": len(records),
        "go_linked": sum(1 for record in records if record["go_count"] > 0),
        "ec_linked": sum(1 for record in records if record["ec_count"] > 0),
        "polymer_context": sum(
            1 for record in records if record["has_polymer_context"]
        ),
        "location_context": sum(
            1 for record in records if record["has_location_context"]
        ),
        "multi_ec": sum(1 for record in records if record["has_multiple_ec"]),
    }
    facet_options = build_browser_facet_options(records)
    top_go_facet_ids = {facet["id"] for facet in facet_options["go_facets"]}
    for record in records:
        record["go_facet_ids"] = [
            go_id for go_id in record["go_closure_ids"] if go_id in top_go_facet_ids
        ]
        del record["go_closure_labels"]
    records_json = json.dumps(records, separators=(",", ":"))
    facet_options_json = json.dumps(facet_options, separators=(",", ":"))
    annotation_colors_json = json.dumps(
        {
            "GO+EC": "#1a5f5b",
            "GO only": "#2d7dd2",
            "EC only": "#c2553d",
            "Unannotated": "#8b7e74",
        },
        separators=(",", ":"),
    )
    ec_major_colors_json = json.dumps(
        {
            "1": "#b23a48",
            "2": "#355070",
            "3": "#588157",
            "4": "#bc6c25",
            "5": "#7b2cbf",
            "6": "#2a9d8f",
            "7": "#6c757d",
            "Unclassified": "#8b7e74",
        },
        separators=(",", ":"),
    )

    html_text = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>RHEA Browser</title>
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <style>
    :root {
      --bg: #f5f1e8;
      --panel: rgba(255, 251, 241, 0.96);
      --panel-strong: #fffaf0;
      --ink: #1f2328;
      --muted: #6a5f52;
      --line: #d9ccb8;
      --accent: #1a5f5b;
      --accent-2: #c2553d;
      --shadow: 0 12px 32px rgba(61, 43, 14, 0.08);
      --radius: 18px;
      --mono: "SFMono-Regular", Consolas, "Liberation Mono", Menlo, monospace;
      --sans: "Iowan Old Style", "Palatino Linotype", "Book Antiqua", Palatino, Georgia, serif;
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      background:
        radial-gradient(circle at top left, rgba(194, 85, 61, 0.08), transparent 28%),
        radial-gradient(circle at top right, rgba(26, 95, 91, 0.12), transparent 32%),
        linear-gradient(180deg, #f7f3eb 0%, #f2ede1 100%);
      color: var(--ink);
      font-family: var(--sans);
    }
    .page {
      max-width: 1720px;
      margin: 0 auto;
      padding: 28px 24px 40px;
    }
    .hero {
      display: grid;
      grid-template-columns: 1.7fr 1fr;
      gap: 18px;
      margin-bottom: 20px;
    }
    .hero-card, .stat-card, .panel {
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: var(--radius);
      box-shadow: var(--shadow);
    }
    .hero-card {
      padding: 22px 24px;
    }
    .kicker {
      color: var(--accent);
      text-transform: uppercase;
      letter-spacing: 0.08em;
      font-size: 0.82rem;
      margin-bottom: 10px;
      font-weight: 700;
    }
    h1 {
      margin: 0 0 10px;
      font-size: clamp(2rem, 4vw, 3.4rem);
      line-height: 1.02;
      letter-spacing: -0.03em;
    }
    .lead {
      margin: 0;
      color: var(--muted);
      font-size: 1.02rem;
      line-height: 1.5;
    }
    .hero-note {
      margin: 14px 0 0;
      color: var(--accent-2);
      font-size: 0.95rem;
      line-height: 1.45;
      font-weight: 600;
    }
    .hero-links {
      display: flex;
      flex-wrap: wrap;
      gap: 14px;
      margin-top: 16px;
    }
    .hero-link {
      color: var(--accent);
      text-decoration: none;
      font-weight: 700;
    }
    .hero-link:hover {
      text-decoration: underline;
    }
    .stats {
      display: grid;
      grid-template-columns: repeat(3, minmax(0, 1fr));
      gap: 14px;
    }
    .stat-card {
      padding: 18px;
    }
    .stat-label {
      color: var(--muted);
      font-size: 0.85rem;
      text-transform: uppercase;
      letter-spacing: 0.08em;
      margin-bottom: 6px;
    }
    .stat-value {
      font-size: 2rem;
      font-weight: 700;
    }
    .controls {
      display: grid;
      grid-template-columns: minmax(260px, 2fr) minmax(190px, 1fr) minmax(220px, 1fr) auto;
      gap: 12px;
      margin-bottom: 18px;
    }
    input, select, button {
      width: 100%;
      padding: 12px 14px;
      border-radius: 12px;
      border: 1px solid var(--line);
      background: var(--panel-strong);
      color: var(--ink);
      font: inherit;
    }
    button {
      background: var(--accent);
      color: white;
      font-weight: 700;
      cursor: pointer;
      border-color: transparent;
    }
    .layout {
      display: grid;
      grid-template-columns: 410px 1fr;
      gap: 18px;
      align-items: start;
    }
    .sidebar {
      display: grid;
      gap: 18px;
    }
    .panel {
      overflow: hidden;
    }
    .panel-header {
      padding: 16px 18px 10px;
      border-bottom: 1px solid var(--line);
      display: flex;
      justify-content: space-between;
      align-items: baseline;
      gap: 12px;
    }
    .panel-title {
      font-size: 1.05rem;
      font-weight: 700;
    }
    .panel-meta {
      color: var(--muted);
      font-size: 0.88rem;
    }
    .facet-body {
      padding: 16px 18px 18px;
      display: grid;
      gap: 18px;
    }
    .facet-group {
      display: grid;
      gap: 10px;
    }
    .facet-label {
      color: var(--muted);
      font-size: 0.82rem;
      text-transform: uppercase;
      letter-spacing: 0.08em;
      font-weight: 700;
    }
    .facet-controls {
      display: grid;
      gap: 10px;
    }
    .toggle-row {
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
    }
    .facet-chip-row {
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
    }
    .toggle-button {
      width: auto;
      background: var(--panel-strong);
      color: var(--ink);
      border: 1px solid var(--line);
      font-weight: 600;
      padding: 8px 12px;
      text-align: left;
      white-space: normal;
    }
    .toggle-button.active {
      background: var(--accent);
      color: #fff;
      border-color: transparent;
    }
    .results-list {
      max-height: 980px;
      overflow: auto;
    }
    .result-item {
      padding: 14px 18px;
      border-top: 1px solid rgba(217, 204, 184, 0.7);
      cursor: pointer;
    }
    .result-item.discrepancy {
      background: rgba(194, 85, 61, 0.04);
    }
    .result-item:hover {
      background: rgba(26, 95, 91, 0.05);
    }
    .result-item.selected {
      background: rgba(26, 95, 91, 0.1);
      border-left: 4px solid var(--accent);
      padding-left: 14px;
    }
    .result-id {
      color: var(--accent);
      font-family: var(--mono);
      font-size: 0.88rem;
      margin-bottom: 4px;
    }
    .result-label {
      font-weight: 700;
      line-height: 1.35;
      margin-bottom: 8px;
    }
    .result-top {
      display: flex;
      justify-content: space-between;
      align-items: start;
      gap: 10px;
    }
    .result-status-badge {
      display: inline-flex;
      align-items: center;
      border-radius: 999px;
      padding: 4px 10px;
      background: rgba(26, 95, 91, 0.12);
      color: var(--accent);
      border: 1px solid rgba(26, 95, 91, 0.16);
      font-size: 0.76rem;
      font-weight: 700;
      text-transform: uppercase;
      letter-spacing: 0.06em;
      white-space: nowrap;
    }
    .result-status-badge.mixed,
    .result-status-badge.asserted_only,
    .result-status-badge.inferred_only {
      background: rgba(194, 85, 61, 0.12);
      color: var(--accent-2);
      border-color: rgba(194, 85, 61, 0.16);
    }
    .result-status-badge.unknown_positive {
      background: rgba(230, 119, 0, 0.14);
      color: #9a5b00;
      border-color: rgba(230, 119, 0, 0.22);
    }
    .result-status-badge.ec_backed_positive {
      background: rgba(32, 80, 160, 0.12);
      color: #204f9e;
      border-color: rgba(32, 80, 160, 0.2);
    }
    .result-meta {
      color: var(--muted);
      font-size: 0.86rem;
      display: flex;
      gap: 10px;
      flex-wrap: wrap;
    }
    .result-flag {
      color: var(--accent-2);
      font-weight: 700;
    }
    .result-annotations {
      display: flex;
      flex-wrap: wrap;
      gap: 6px;
      margin-top: 10px;
    }
    .result-rule-lines {
      display: grid;
      gap: 4px;
      margin-top: 10px;
      font-size: 0.82rem;
      color: var(--muted);
      line-height: 1.35;
    }
    .result-rule-line strong {
      color: var(--ink);
    }
    .result-rule-heading {
      color: var(--ink);
      font-size: 0.76rem;
      font-weight: 700;
      letter-spacing: 0.08em;
      text-transform: uppercase;
      margin-bottom: 2px;
    }
    .result-rule-alert {
      margin-top: 4px;
      border-left: 3px solid rgba(194, 85, 61, 0.4);
      padding-left: 8px;
      color: var(--accent-2);
      font-weight: 600;
    }
    .result-pill {
      display: inline-flex;
      border: 1px solid var(--line);
      border-radius: 999px;
      padding: 3px 9px;
      background: var(--panel-strong);
      font-size: 0.78rem;
      font-family: var(--mono);
    }
    .class-inline-link {
      color: var(--accent);
      font-weight: 700;
      text-decoration: none;
    }
    .class-inline-link:hover {
      text-decoration: underline;
    }
    .main-stack {
      display: grid;
      grid-template-rows: minmax(500px, 58vh) auto;
      gap: 18px;
    }
    .plot-wrap {
      padding: 12px 12px 0;
    }
    #embedding-plot {
      width: 100%;
      height: 100%;
      min-height: 560px;
    }
    .detail-grid {
      display: grid;
      grid-template-columns: 1.15fr 1fr;
      gap: 18px;
      padding: 18px;
    }
    .slot-grid {
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 12px;
      margin-top: 14px;
      margin-bottom: 18px;
    }
    .slot {
      background: rgba(255, 251, 241, 0.7);
      border: 1px solid var(--line);
      border-radius: 14px;
      padding: 12px 14px;
    }
    .slot-name {
      color: var(--muted);
      text-transform: uppercase;
      letter-spacing: 0.08em;
      font-size: 0.76rem;
      margin-bottom: 6px;
    }
    .slot-value {
      font-size: 0.98rem;
      line-height: 1.45;
      word-break: break-word;
    }
    .classification-callout {
      border-radius: 12px;
      border: 1px solid var(--line);
      background: rgba(26, 95, 91, 0.06);
      padding: 10px 12px;
    }
    .classification-callout.discrepancy {
      background: rgba(194, 85, 61, 0.08);
      border-color: rgba(194, 85, 61, 0.18);
    }
    .classification-callout-status {
      color: var(--ink);
      font-weight: 700;
      margin-bottom: 4px;
    }
    .classification-callout-detail {
      margin-top: 6px;
      color: var(--accent-2);
      font-weight: 600;
    }
    .pill-row {
      display: flex;
      gap: 8px;
      flex-wrap: wrap;
    }
    .pill {
      display: inline-flex;
      align-items: center;
      border: 1px solid var(--line);
      border-radius: 999px;
      padding: 4px 10px;
      background: var(--panel-strong);
      font-size: 0.82rem;
      font-family: var(--mono);
    }
    .pill-link {
      text-decoration: none;
      color: inherit;
    }
    .pill-link:hover {
      color: var(--accent);
      border-color: var(--accent);
    }
    .participant-section {
      margin-top: 10px;
    }
    .participant-section h3 {
      margin: 0 0 10px;
      font-size: 0.98rem;
    }
    .participant-list {
      display: grid;
      gap: 8px;
    }
    .participant-item {
      border: 1px solid var(--line);
      border-radius: 12px;
      padding: 10px 12px;
      background: rgba(255, 251, 241, 0.72);
    }
    .participant-name {
      font-weight: 700;
      margin-bottom: 4px;
    }
    .participant-meta {
      color: var(--muted);
      font-size: 0.84rem;
      display: flex;
      flex-wrap: wrap;
      gap: 10px;
    }
    .detail-title {
      display: flex;
      justify-content: space-between;
      align-items: start;
      gap: 12px;
      margin-bottom: 10px;
    }
    .detail-title h2 {
      margin: 0;
      font-size: 1.5rem;
      line-height: 1.15;
    }
    .detail-link {
      color: var(--accent);
      text-decoration: none;
      font-weight: 700;
      white-space: nowrap;
    }
    .detail-link:hover {
      text-decoration: underline;
    }
    .detail-subtitle {
      color: var(--muted);
      margin: 0 0 10px;
    }
    .embedding-box {
      border: 1px solid var(--line);
      border-radius: 14px;
      background: rgba(255, 251, 241, 0.72);
      padding: 12px 14px;
      font-size: 0.9rem;
      line-height: 1.5;
      max-height: 240px;
      overflow: auto;
    }
    .json-box {
      max-height: 320px;
      white-space: pre-wrap;
      font-family: var(--mono);
      font-size: 0.84rem;
    }
    .empty-state {
      padding: 20px 18px;
      color: var(--muted);
    }
    @media (max-width: 1100px) {
      .hero, .layout, .detail-grid, .controls {
        grid-template-columns: 1fr;
      }
      .stats {
        grid-template-columns: repeat(2, minmax(0, 1fr));
      }
      .main-stack {
        grid-template-rows: minmax(420px, 55vh) auto;
      }
    }
    @media (max-width: 680px) {
      .stats {
        grid-template-columns: 1fr;
      }
    }
  </style>
</head>
<body>
  <div class="page">
    <section class="hero">
      <div class="hero-card">
        <div class="kicker">Autarch Browser</div>
        <h1>RHEA reaction browser with linked text embeddings</h1>
        <p class="lead">
          Browse curated RHEA reactions as structured objects, search across
          reaction definitions and participant metadata, and inspect a linked
          2D embedding view built from reaction equation text plus participant
          descriptors.
        </p>
        <p class="hero-note">
          Use the Rule status facet to isolate mismatches between asserted classes and agent classification.
        </p>
        <div class="hero-links">
          <a class="hero-link" href="index.html">Report index</a>
          <a class="hero-link" href="../index.html">Docs home</a>
        </div>
      </div>
      <div class="stats">
        <div class="stat-card">
          <div class="stat-label">Reactions</div>
          <div class="stat-value">__STAT_COUNT__</div>
        </div>
        <div class="stat-card">
          <div class="stat-label">GO-linked</div>
          <div class="stat-value">__STAT_GO_LINKED__</div>
        </div>
        <div class="stat-card">
          <div class="stat-label">EC-linked</div>
          <div class="stat-value">__STAT_EC_LINKED__</div>
        </div>
        <div class="stat-card">
          <div class="stat-label">Polymer Context</div>
          <div class="stat-value">__STAT_POLYMER_CONTEXT__</div>
        </div>
        <div class="stat-card">
          <div class="stat-label">Location Context</div>
          <div class="stat-value">__STAT_LOCATION_CONTEXT__</div>
        </div>
        <div class="stat-card">
          <div class="stat-label">Multi-EC</div>
          <div class="stat-value">__STAT_MULTI_EC__</div>
        </div>
      </div>
    </section>

    <div class="controls">
      <input id="search-input" type="search" placeholder="Search RHEA ID, definition, participants, GO or EC">
      <select id="annotation-filter">
        <option value="all">All annotation groups</option>
        <option value="GO+EC">GO+EC</option>
        <option value="GO only">GO only</option>
        <option value="EC only">EC only</option>
        <option value="Unannotated">Unannotated</option>
      </select>
      <select id="color-mode">
        <option value="annotation">Color by annotation source</option>
        <option value="ec_major">Color by EC major class</option>
      </select>
      <button id="reset-button" type="button">Reset</button>
    </div>

    <section class="layout">
      <aside class="sidebar">
        <section class="panel">
          <div class="panel-header">
            <div class="panel-title">Facets</div>
            <div class="panel-meta" id="facet-summary"></div>
          </div>
          <div class="facet-body">
            <div class="facet-group">
              <div class="facet-label">EC hierarchy</div>
              <div class="facet-controls">
                <select id="ec-major-filter"></select>
                <select id="ec-subclass-filter"></select>
              </div>
            </div>
            <div class="facet-group">
              <div class="facet-label">Operational</div>
              <div class="facet-controls">
                <select id="rule-status-filter"></select>
                <select id="participant-filter"></select>
              </div>
              <div class="toggle-row">
                <button class="toggle-button" type="button" data-flag="polymerOnly">Polymer context</button>
                <button class="toggle-button" type="button" data-flag="locationOnly">Location annotated</button>
                <button class="toggle-button" type="button" data-flag="multiEcOnly">Multiple EC mappings</button>
              </div>
            </div>
            <div class="facet-group">
              <div class="facet-label">GO closure</div>
              <div id="go-facet-chips" class="facet-chip-row"></div>
            </div>
          </div>
        </section>

        <section class="panel">
          <div class="panel-header">
            <div class="panel-title">Search results</div>
            <div class="panel-meta" id="result-count"></div>
          </div>
          <div id="results-list" class="results-list"></div>
        </section>
      </aside>

      <div class="main-stack">
        <section class="panel plot-wrap">
          <div class="panel-header">
            <div class="panel-title">Text embedding map</div>
            <div class="panel-meta" id="plot-meta">definition/equation + participant descriptors</div>
          </div>
          <div id="embedding-plot"></div>
        </section>

        <section class="panel">
          <div class="panel-header">
            <div class="panel-title">Reaction detail</div>
            <div class="panel-meta" id="detail-meta"></div>
          </div>
          <div id="detail-pane"></div>
        </section>
      </div>
    </section>
  </div>

  <script>
    const records = __RECORDS_JSON__;
    const facetOptions = __FACET_OPTIONS_JSON__;
    const annotationColors = __ANNOTATION_COLORS_JSON__;
    const ecMajorColors = __EC_MAJOR_COLORS_JSON__;
    const state = {
      query: "",
      annotationGroup: "all",
      goFacet: "all",
      ecMajor: "all",
      ecSubclass: "all",
      participantBucket: "all",
      ruleStatus: "all",
      polymerOnly: false,
      locationOnly: false,
      multiEcOnly: false,
      colorMode: "annotation",
      selectedId: records.length ? records[0].rhea_id : null,
    };

    function escapeHtml(value) {
      return String(value)
        .replaceAll("&", "&amp;")
        .replaceAll("<", "&lt;")
        .replaceAll(">", "&gt;")
        .replaceAll('"', "&quot;");
    }

    function optionHtml(value, label, selected) {
      const selectedAttr = selected ? " selected" : "";
      return `<option value="${escapeHtml(value)}"${selectedAttr}>${escapeHtml(label)}</option>`;
    }

    function populateFacetControls() {
      const majorSelect = document.getElementById("ec-major-filter");
      majorSelect.innerHTML = [
        optionHtml("all", "All EC major classes", state.ecMajor === "all"),
        ...facetOptions.ec_majors.map((option) =>
          optionHtml(option.value, `${option.label} (${option.count})`, state.ecMajor === option.value)
        ),
      ].join("");

      const subclassOptions = facetOptions.ec_subclasses.filter(
        (option) => state.ecMajor === "all" || option.major === state.ecMajor
      );
      if (
        state.ecSubclass !== "all" &&
        !subclassOptions.some((option) => option.value === state.ecSubclass)
      ) {
        state.ecSubclass = "all";
      }
      const subclassSelect = document.getElementById("ec-subclass-filter");
      subclassSelect.innerHTML = [
        optionHtml("all", "All EC subclasses", state.ecSubclass === "all"),
        ...subclassOptions.map((option) =>
          optionHtml(option.value, `${option.label} (${option.count})`, state.ecSubclass === option.value)
        ),
      ].join("");

      const participantSelect = document.getElementById("participant-filter");
      participantSelect.innerHTML = [
        optionHtml("all", "All participant-count buckets", state.participantBucket === "all"),
        ...facetOptions.participant_buckets.map((option) =>
          optionHtml(option.value, `${option.label} (${option.count})`, state.participantBucket === option.value)
        ),
      ].join("");

      const ruleStatusSelect = document.getElementById("rule-status-filter");
      ruleStatusSelect.innerHTML = [
        optionHtml("all", "All rule statuses", state.ruleStatus === "all"),
        ...facetOptions.rule_statuses.map((option) =>
          optionHtml(option.value, `${option.label} (${option.count})`, state.ruleStatus === option.value)
        ),
      ].join("");
    }

    function syncToggleButtons() {
      document.querySelectorAll(".toggle-button").forEach((button) => {
        const flagName = button.dataset.flag;
        button.classList.toggle("active", Boolean(state[flagName]));
      });
    }

    function renderGoFacetChips() {
      const container = document.getElementById("go-facet-chips");
      container.innerHTML = facetOptions.go_facets
        .map((facet) => `
          <button
            class="toggle-button ${state.goFacet === facet.id ? "active" : ""}"
            type="button"
            data-go-facet="${facet.id}"
          >
            ${escapeHtml(facet.label)} (${facet.count})
          </button>
        `)
        .join("");
      container.querySelectorAll("[data-go-facet]").forEach((button) => {
        button.addEventListener("click", () => {
          state.goFacet =
            state.goFacet === button.dataset.goFacet ? "all" : button.dataset.goFacet;
          renderAll();
        });
      });
    }

    function activeFilterCount() {
      return [
        state.annotationGroup !== "all",
        state.goFacet !== "all",
        state.ecMajor !== "all",
        state.ecSubclass !== "all",
        state.participantBucket !== "all",
        state.ruleStatus !== "all",
        state.polymerOnly,
        state.locationOnly,
        state.multiEcOnly,
        state.query.trim() !== "",
      ].filter(Boolean).length;
    }

    function updateFacetSummary(filteredCount) {
      const activeCount = activeFilterCount();
      const summary = activeCount
        ? `${activeCount} active · ${filteredCount} match`
        : `Browsing ${records.length} reactions`;
      document.getElementById("facet-summary").textContent = summary;
    }

    function filteredRecords() {
      const query = state.query.trim().toLowerCase();
      return records.filter((record) => {
        const annotationOk =
          state.annotationGroup === "all" ||
          record.annotation_group === state.annotationGroup;
        const goFacetOk =
          state.goFacet === "all" ||
          record.go_facet_ids.includes(state.goFacet);
        const ecMajorOk =
          state.ecMajor === "all" ||
          record.ec_major_classes.includes(state.ecMajor);
        const ecSubclassOk =
          state.ecSubclass === "all" ||
          record.ec_subclasses.includes(state.ecSubclass);
        const participantOk =
          state.participantBucket === "all" ||
          record.participant_bucket === state.participantBucket;
        const ruleStatusOk =
          state.ruleStatus === "all" ||
          (state.ruleStatus === "covered"
            ? record.has_rule_coverage
            : record.rule_status === state.ruleStatus);
        const polymerOk = !state.polymerOnly || record.has_polymer_context;
        const locationOk = !state.locationOnly || record.has_location_context;
        const multiEcOk = !state.multiEcOnly || record.has_multiple_ec;
        const searchOk = !query || record.search_text.includes(query);
        return (
          annotationOk &&
          goFacetOk &&
          ecMajorOk &&
          ecSubclassOk &&
          participantOk &&
          ruleStatusOk &&
          polymerOk &&
          locationOk &&
          multiEcOk &&
          searchOk
        );
      });
    }

    function getSelectedRecord(filtered) {
      if (!filtered.length) {
        return null;
      }
      const selected =
        filtered.find((record) => record.rhea_id === state.selectedId) || filtered[0];
      state.selectedId = selected.rhea_id;
      return selected;
    }

    function participantHtml(participant) {
      const meta = [];
      if (participant.chebi_id) meta.push(participant.chebi_id);
      if (participant.count && participant.count !== 1) meta.push(`count=${participant.count}`);
      if (participant.stoichiometry) meta.push(`stoich=${participant.stoichiometry}`);
      if (participant.polymer_type) meta.push(`polymer=${participant.polymer_type}`);
      if (participant.polymer_index) meta.push(`index=${participant.polymer_index}`);
      if (participant.location) meta.push(`location=${participant.location}`);
      if (participant.monomer) meta.push(`monomer=${participant.monomer}`);
      return `
        <div class="participant-item">
          <div class="participant-name">${escapeHtml(participant.name || participant.chebi_id || "unnamed participant")}</div>
          <div class="participant-meta">${meta.map((item) => `<span>${escapeHtml(item)}</span>`).join("")}</div>
        </div>
      `;
    }

    function compactRuleList(items) {
      if (!items.length) {
        return "none";
      }
      const limit = 3;
      if (items.length <= limit) {
        return items.join(", ");
      }
      return `${items.slice(0, limit).join(", ")} +${items.length - limit} more`;
    }

    function classLinkHtml(entry, extraClass = "class-inline-link") {
      return `<a class="${extraClass}" href="${encodeURI(entry.href)}">${escapeHtml(entry.display)}</a>`;
    }

    function compactRuleLinksHtml(entries) {
      if (!entries.length) {
        return "none";
      }
      const limit = 3;
      const visible = entries.slice(0, limit).map((entry) => classLinkHtml(entry)).join(", ");
      if (entries.length <= limit) {
        return visible;
      }
      return `${visible} +${entries.length - limit} more`;
    }

    function discrepancySummary(record) {
      const parts = [];
      if (record.missing_rule_classes.length) {
        parts.push(`Asserted only: ${compactRuleList(record.missing_rule_classes)}`);
      }
      if (record.extra_rule_classes.length) {
        parts.push(`Agent only: ${compactRuleList(record.extra_rule_classes)}`);
      }
      return parts.join(" | ");
    }

    function classificationSummaryText(record) {
      if (record.rule_status === "agreement") {
        return "Asserted classes and agent classification are consistent.";
      }
      if (record.rule_status === "asserted_only") {
        return "Some asserted classes are not present in the agent classification.";
      }
      if (record.rule_status === "inferred_only") {
        return "The agent assigned classes with no asserted support.";
      }
      if (record.rule_status === "unknown_positive") {
        return "The agent assigned classes to a reaction with neither GO nor EC support.";
      }
      if (record.rule_status === "ec_backed_positive") {
        return "The agent assigned classes to a reaction that already has EC support but is missing from the GO-labeled benchmark.";
      }
      if (record.rule_status === "mixed") {
        return "The asserted classes and agent classification only partially overlap.";
      }
      return "No Python-rule classification is available for this reaction.";
    }

    function pillRowHtml(items) {
      return items.length
        ? items.map((item) => `<span class="pill">${escapeHtml(item)}</span>`).join("")
        : '<span class="pill">none</span>';
    }

    function classPillRowHtml(entries) {
      return entries.length
        ? entries
            .map((entry) => `<a class="pill pill-link" href="${encodeURI(entry.href)}">${escapeHtml(entry.display)}</a>`)
            .join("")
        : '<span class="pill">none</span>';
    }

    function renderResults(filtered) {
      const resultsList = document.getElementById("results-list");
      const resultCount = document.getElementById("result-count");
      const topRecords = filtered.slice(0, 250);
      resultCount.textContent =
        filtered.length > topRecords.length
          ? `${topRecords.length} of ${filtered.length} shown`
          : `${filtered.length} shown`;

      if (!filtered.length) {
        resultsList.innerHTML = '<div class="empty-state">No reactions match the current search.</div>';
        return;
      }

      resultsList.innerHTML = topRecords
        .map((record) => `
          <div class="result-item ${record.rhea_id === state.selectedId ? "selected" : ""} ${record.rule_status !== "agreement" && record.has_rule_coverage ? "discrepancy" : ""}" data-rhea-id="${record.rhea_id}">
            <div class="result-id">${escapeHtml(record.rhea_id)}</div>
            <div class="result-top">
              <div class="result-label">${escapeHtml(record.label)}</div>
              <div class="result-status-badge ${record.rule_status}">${escapeHtml(record.rule_badge)}</div>
            </div>
            <div class="result-meta">
              <span>${escapeHtml(record.annotation_group)}</span>
              <span>${escapeHtml(record.primary_ec_major_label)}</span>
              <span>${record.participant_total} participants</span>
              ${record.has_polymer_context ? '<span class="result-flag">polymer</span>' : ""}
              ${record.has_location_context ? '<span class="result-flag">location</span>' : ""}
            </div>
            <div class="result-rule-lines">
              <div class="result-rule-heading">Classification</div>
              <div class="result-rule-line"><strong>Asserted classes:</strong> ${compactRuleLinksHtml(record.asserted_rule_entries)}</div>
              <div class="result-rule-line"><strong>Agent classification:</strong> ${compactRuleLinksHtml(record.inferred_rule_entries)}</div>
              ${record.rule_status === "agreement" || !discrepancySummary(record)
                ? ""
                : `<div class="result-rule-alert">${escapeHtml(discrepancySummary(record))}</div>`}
            </div>
            <div class="result-annotations">
              ${record.go_term_entries.slice(0, 2).map((entry) => `<span class="result-pill">${escapeHtml(entry.id)} ${escapeHtml(entry.label)}</span>`).join("")}
              ${record.ec_number_entries.slice(0, 2).map((entry) => `<span class="result-pill">${escapeHtml(entry.id)}${entry.label ? ` ${escapeHtml(entry.label)}` : ""}</span>`).join("")}
            </div>
          </div>
        `)
        .join("");

      resultsList.querySelectorAll(".result-item").forEach((element) => {
        element.addEventListener("click", () => {
          state.selectedId = element.dataset.rheaId;
          renderAll();
        });
      });
    }

    function renderDetail(record) {
      const detailPane = document.getElementById("detail-pane");
      const detailMeta = document.getElementById("detail-meta");
      if (!record) {
        detailMeta.textContent = "";
        detailPane.innerHTML = '<div class="empty-state">Select a reaction to inspect its slots and embedding text.</div>';
        return;
      }

      detailMeta.textContent = `${record.annotation_group} · ${record.primary_ec_major_label} · ${record.rule_status_label}`;
      const goTerms = record.go_term_entries.length
        ? record.go_term_entries
            .map((entry) => `<span class="pill">${escapeHtml(entry.id)} ${escapeHtml(entry.label)}</span>`)
            .join("")
        : '<span class="pill">none</span>';
      const ecNumbers = record.ec_number_entries.length
        ? record.ec_number_entries
            .map((entry) => `<span class="pill">${escapeHtml(entry.id)}${entry.label ? ` ${escapeHtml(entry.label)}` : ""}</span>`)
            .join("")
        : '<span class="pill">none</span>';
      const ecHierarchy = [
        ...record.ec_major_labels.map((term) => `<span class="pill">${escapeHtml(term)}</span>`),
        ...record.ec_subclasses.map((term) => `<span class="pill">${escapeHtml(term)}</span>`),
      ].join("") || '<span class="pill">unclassified</span>';
      const assertedClasses = classPillRowHtml(record.asserted_rule_entries);
      const agentClasses = classPillRowHtml(record.inferred_rule_entries);
      const discrepancyText = discrepancySummary(record);
      const classificationCallout = `
        <div class="classification-callout ${record.rule_status !== "agreement" && record.has_rule_coverage ? "discrepancy" : ""}">
          <div class="classification-callout-status">${escapeHtml(record.rule_status_label)}</div>
          <div>${escapeHtml(classificationSummaryText(record))}</div>
          ${discrepancyText ? `<div class="classification-callout-detail">${escapeHtml(discrepancyText)}</div>` : ""}
        </div>
      `;
      const contextPills = [
        `<span class="pill">${escapeHtml(record.participant_bucket)} participants</span>`,
        record.has_polymer_context ? '<span class="pill">polymer context</span>' : "",
        record.has_location_context ? '<span class="pill">location annotated</span>' : "",
        record.has_multiple_ec ? '<span class="pill">multiple EC mappings</span>' : "",
      ].join("");
      const rawRecord = {
        rhea_id: record.rhea_id,
        label: record.label,
        annotation_group: record.annotation_group,
        go_terms: record.go_term_entries,
        ec_numbers: record.ec_number_entries,
        rule_status: record.rule_status_label,
        asserted_classes: record.asserted_rule_entries,
        agent_classification: record.inferred_rule_entries,
        discrepancy: discrepancyText || null,
        missing_rule_classes: record.missing_rule_classes,
        extra_rule_classes: record.extra_rule_classes,
        ec_major_classes: record.ec_major_classes,
        ec_subclasses: record.ec_subclasses,
        participant_bucket: record.participant_bucket,
        has_polymer_context: record.has_polymer_context,
        has_location_context: record.has_location_context,
        left_participants: record.left_participants,
        right_participants: record.right_participants,
      };

      detailPane.innerHTML = `
        <div class="detail-grid">
          <div>
            <div class="detail-title">
              <div>
                <h2>${escapeHtml(record.rhea_id)}</h2>
                <p class="detail-subtitle">${escapeHtml(record.label)}</p>
              </div>
              <a class="detail-link" href="${record.url}" target="_blank" rel="noreferrer">Open RHEA</a>
            </div>

            <div class="slot-grid">
              <div class="slot">
                <div class="slot-name">Reactants</div>
                <div class="slot-value">${escapeHtml(record.left_summary || "none")}</div>
              </div>
              <div class="slot">
                <div class="slot-name">Products</div>
                <div class="slot-value">${escapeHtml(record.right_summary || "none")}</div>
              </div>
              <div class="slot">
                <div class="slot-name">GO Terms</div>
                <div class="slot-value"><div class="pill-row">${goTerms}</div></div>
              </div>
              <div class="slot">
                <div class="slot-name">EC Numbers</div>
                <div class="slot-value"><div class="pill-row">${ecNumbers}</div></div>
              </div>
              <div class="slot">
                <div class="slot-name">EC Hierarchy</div>
                <div class="slot-value"><div class="pill-row">${ecHierarchy}</div></div>
              </div>
              <div class="slot">
                <div class="slot-name">Asserted Classes</div>
                <div class="slot-value"><div class="pill-row">${assertedClasses}</div></div>
              </div>
              <div class="slot">
                <div class="slot-name">Agent Classification</div>
                <div class="slot-value"><div class="pill-row">${agentClasses}</div></div>
              </div>
              <div class="slot">
                <div class="slot-name">Classification Comparison</div>
                <div class="slot-value">${classificationCallout}</div>
              </div>
              <div class="slot">
                <div class="slot-name">Context</div>
                <div class="slot-value"><div class="pill-row">${contextPills}</div></div>
              </div>
            </div>

            <div class="participant-section">
              <h3>Left participants</h3>
              <div class="participant-list">
                ${record.left_participants.length ? record.left_participants.map(participantHtml).join("") : '<div class="empty-state">No left participants recorded.</div>'}
              </div>
            </div>
          </div>

          <div>
            <div class="participant-section">
              <h3>Right participants</h3>
              <div class="participant-list">
                ${record.right_participants.length ? record.right_participants.map(participantHtml).join("") : '<div class="empty-state">No right participants recorded.</div>'}
              </div>
            </div>

            <div class="participant-section">
              <h3>Embedding text</h3>
              <div class="embedding-box">${escapeHtml(record.embedding_text)}</div>
            </div>

            <div class="participant-section">
              <h3>Structured record</h3>
              <div class="embedding-box json-box">${escapeHtml(JSON.stringify(rawRecord, null, 2))}</div>
            </div>
          </div>
        </div>
      `;
    }

    function colorDescriptor(record) {
      if (state.colorMode === "ec_major") {
        const key = record.primary_ec_major;
        return {
          key,
          label: record.primary_ec_major_label,
          color: ecMajorColors[key] || ecMajorColors.Unclassified,
        };
      }
      return {
        key: record.annotation_group,
        label: record.annotation_group,
        color: annotationColors[record.annotation_group] || annotationColors.Unannotated,
      };
    }

    function renderPlot(filtered, selected) {
      const plotMeta = document.getElementById("plot-meta");
      plotMeta.textContent =
        state.colorMode === "ec_major"
          ? "colored by EC major class; definition/equation + participant descriptors"
          : "colored by annotation source; definition/equation + participant descriptors";

      if (!filtered.length) {
        Plotly.react(
          "embedding-plot",
          [],
          {
            margin: { t: 10, r: 10, b: 50, l: 50 },
            paper_bgcolor: "rgba(0,0,0,0)",
            plot_bgcolor: "rgba(255,255,255,0)",
            annotations: [
              {
                text: "No reactions match the current filters.",
                showarrow: false,
                x: 0.5,
                y: 0.5,
                xref: "paper",
                yref: "paper",
                font: { size: 16, color: "#6a5f52" },
              },
            ],
          },
          { responsive: true, displaylogo: false }
        );
        return;
      }

      const grouped = new Map();
      filtered.forEach((record) => {
        const descriptor = colorDescriptor(record);
        if (!grouped.has(descriptor.key)) {
          grouped.set(descriptor.key, { ...descriptor, items: [] });
        }
        grouped.get(descriptor.key).items.push(record);
      });

      const traces = Array.from(grouped.values())
        .sort((left, right) => left.label.localeCompare(right.label, undefined, { numeric: true }))
        .map((group) => ({
          type: "scattergl",
          mode: "markers",
          name: group.label,
          x: group.items.map((record) => record.x),
          y: group.items.map((record) => record.y),
          customdata: group.items.map((record) => [
            record.rhea_id,
            record.annotation_group,
            record.label,
            record.go_count,
            record.ec_count,
            record.primary_ec_major_label,
          ]),
          hovertemplate:
            "<b>%{customdata[0]}</b><br>%{customdata[2]}<br>Annotation=%{customdata[1]}<br>EC major=%{customdata[5]}<br>GO=%{customdata[3]} EC=%{customdata[4]}<extra>%{fullData.name}</extra>",
          marker: {
            size: group.items.map(
              (record) => 7 + Math.sqrt(Math.max(record.go_count + record.ec_count, 1)) * 2.5
            ),
            color: group.color,
            opacity: 0.78,
            line: { width: 0.4, color: "rgba(33, 28, 20, 0.5)" },
          },
        }));

      if (selected) {
        traces.push({
          type: "scattergl",
          mode: "markers+text",
          x: [selected.x],
          y: [selected.y],
          text: [selected.rhea_id],
          textposition: "top center",
          textfont: { size: 11, color: "#1f2328" },
          hoverinfo: "skip",
          showlegend: false,
          marker: {
            size: 22,
            color: "#ffd166",
            line: { width: 2.0, color: "#1f2328" },
            symbol: "diamond",
          },
        });
      }

      Plotly.react(
        "embedding-plot",
        traces,
        {
          margin: { t: 10, r: 10, b: 50, l: 50 },
          paper_bgcolor: "rgba(0,0,0,0)",
          plot_bgcolor: "rgba(255,255,255,0)",
          xaxis: { title: "PC1", zeroline: false, gridcolor: "rgba(217, 204, 184, 0.6)" },
          yaxis: { title: "PC2", zeroline: false, gridcolor: "rgba(217, 204, 184, 0.6)" },
          legend: { orientation: "h", yanchor: "bottom", y: 1.01, x: 0.0 },
        },
        { responsive: true, displaylogo: false }
      );

      const plot = document.getElementById("embedding-plot");
      if (plot.dataset.boundClick !== "true") {
        plot.on("plotly_click", (event) => {
          const rheaId = event?.points?.[0]?.customdata?.[0];
          if (rheaId) {
            state.selectedId = rheaId;
            renderAll();
          }
        });
        plot.dataset.boundClick = "true";
      }
    }

    function renderAll() {
      populateFacetControls();
      syncToggleButtons();
      renderGoFacetChips();
      const filtered = filteredRecords();
      const selected = getSelectedRecord(filtered);
      updateFacetSummary(filtered.length);
      renderResults(filtered);
      renderDetail(selected);
      renderPlot(filtered, selected);
    }

    document.getElementById("search-input").addEventListener("input", (event) => {
      state.query = event.target.value;
      renderAll();
    });

    document.getElementById("annotation-filter").addEventListener("change", (event) => {
      state.annotationGroup = event.target.value;
      renderAll();
    });

    document.getElementById("color-mode").addEventListener("change", (event) => {
      state.colorMode = event.target.value;
      renderAll();
    });

    document.getElementById("ec-major-filter").addEventListener("change", (event) => {
      state.ecMajor = event.target.value;
      if (state.ecMajor === "all") {
        state.ecSubclass = "all";
      }
      renderAll();
    });

    document.getElementById("ec-subclass-filter").addEventListener("change", (event) => {
      state.ecSubclass = event.target.value;
      renderAll();
    });

    document.getElementById("participant-filter").addEventListener("change", (event) => {
      state.participantBucket = event.target.value;
      renderAll();
    });

    document.getElementById("rule-status-filter").addEventListener("change", (event) => {
      state.ruleStatus = event.target.value;
      renderAll();
    });

    document.querySelectorAll(".toggle-button").forEach((button) => {
      button.addEventListener("click", () => {
        const flagName = button.dataset.flag;
        state[flagName] = !state[flagName];
        renderAll();
      });
    });

    document.getElementById("reset-button").addEventListener("click", () => {
      state.query = "";
      state.annotationGroup = "all";
      state.goFacet = "all";
      state.ecMajor = "all";
      state.ecSubclass = "all";
      state.participantBucket = "all";
      state.ruleStatus = "all";
      state.polymerOnly = false;
      state.locationOnly = false;
      state.multiEcOnly = false;
      state.colorMode = "annotation";
      document.getElementById("search-input").value = "";
      document.getElementById("annotation-filter").value = "all";
      document.getElementById("color-mode").value = "annotation";
      renderAll();
    });

    renderAll();
  </script>
</body>
</html>
"""
    replacements = {
        "__RECORDS_JSON__": records_json,
        "__FACET_OPTIONS_JSON__": facet_options_json,
        "__ANNOTATION_COLORS_JSON__": annotation_colors_json,
        "__EC_MAJOR_COLORS_JSON__": ec_major_colors_json,
        "__STAT_COUNT__": str(stats["count"]),
        "__STAT_GO_LINKED__": str(stats["go_linked"]),
        "__STAT_EC_LINKED__": str(stats["ec_linked"]),
        "__STAT_POLYMER_CONTEXT__": str(stats["polymer_context"]),
        "__STAT_LOCATION_CONTEXT__": str(stats["location_context"]),
        "__STAT_MULTI_EC__": str(stats["multi_ec"]),
    }
    for token, value in replacements.items():
        html_text = html_text.replace(token, value)

    output_path.write_text(html_text)
    return output_path
