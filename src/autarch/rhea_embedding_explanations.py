"""Local explanations for embedding-based rule classifiers."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, TypedDict

import numpy as np
import pandas as pd

from autarch.classifier import ReactionClassifier
from autarch.datamodel import Participant, Reaction
from autarch.rhea_embedding_benchmark import (
    build_rule_benchmark_feature_matrices,
    build_rule_target_data,
    fit_binary_rule_classifier,
)
from autarch.rhea_text_embeddings import (
    DEFAULT_EMBEDDING_MODEL_NAME,
    DEFAULT_N_FEATURES,
    build_asserted_rule_support,
    build_lexical_embedding_matrices,
    build_linkml_store_embedding_matrix,
    build_linkml_store_embedding_matrices,
    build_reaction_embedding_text,
    load_rhea_text_df,
    load_rule_class_metadata,
    participant_display_text,
    build_text_feature_matrix,
)


class EmbeddingExplanationComponent(TypedDict):
    """One interpretable occlusion unit for a reaction embedding explanation."""

    key: str
    kind: str
    side: str
    label: str
    removed_text: str
    masked_text: str


class EmbeddingExplanationComponentScore(EmbeddingExplanationComponent):
    """Occlusion unit plus score deltas for one class."""

    masked_score: float
    score_delta: float
    supports_class: bool


class EmbeddingRuleComparison(TypedDict):
    """Structured local explanation comparing SVM and rule outputs."""

    rhea_id: str
    reaction_label: str
    reaction_class: str
    model_name: str
    feature_space: str
    embedding_backend: str
    embedding_model_name: str
    asserted_positive: bool
    positive_support: int
    negative_support: int
    svm_score: float
    svm_predicted_positive: bool
    rule_predicted_positive: bool
    rule_explanation: str
    embedding_text: str
    component_scores: list[EmbeddingExplanationComponentScore]


def normalize_rhea_id(rhea_id: str) -> str:
    """Normalize user input into the cached master-RHEA identifier format."""
    clean_id = rhea_id.replace("RHEA:", "").strip()
    return f"RHEA:{clean_id}"


def load_rhea_cache_record(
    rhea_id: str,
    cache_dir: str | Path = "cache",
) -> dict[str, Any]:
    """Load one raw cached RHEA record from the JSONL cache."""
    normalized_id = normalize_rhea_id(rhea_id)
    cache_path = Path(cache_dir) / "rhea_reactions.jsonl"
    if not cache_path.exists():
        raise FileNotFoundError(f"RHEA cache not found at {cache_path}")

    with open(cache_path) as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            if row.get("rhea_id") == normalized_id:
                return row
    raise KeyError(f"{normalized_id} not found in {cache_path}")


def reaction_from_rhea_cache_record(record: dict[str, Any]) -> Reaction:
    """Convert one cached RHEA record into a runtime Reaction object."""
    reaction_dict = record.get("reaction") or {}
    if not reaction_dict:
        raise ValueError(f"{record.get('rhea_id', 'reaction')} has no reaction payload.")

    return Reaction(
        left_participants=[
            Participant(**participant)
            for participant in reaction_dict.get("left_participants", [])
        ],
        right_participants=[
            Participant(**participant)
            for participant in reaction_dict.get("right_participants", [])
        ],
        label=record.get("label", ""),
    )


def build_embedding_explanation_components(
    reaction_label: str,
    left_participants: list[dict[str, Any]],
    right_participants: list[dict[str, Any]],
    include_label: bool = True,
) -> list[EmbeddingExplanationComponent]:
    """Build interpretable masking units for one reaction embedding."""
    components: list[EmbeddingExplanationComponent] = [
        {
            "key": "label",
            "kind": "label",
            "side": "meta",
            "label": "Label text",
            "removed_text": reaction_label.strip(),
            "masked_text": build_reaction_embedding_text(
                "",
                left_participants,
                right_participants,
                include_label=include_label,
            ),
        },
        {
            "key": "lhs",
            "kind": "side",
            "side": "lhs",
            "label": "All reactants",
            "removed_text": "; ".join(
                participant_display_text(participant) for participant in left_participants
            ),
            "masked_text": build_reaction_embedding_text(
                reaction_label,
                [],
                right_participants,
                include_label=include_label,
            ),
        },
        {
            "key": "rhs",
            "kind": "side",
            "side": "rhs",
            "label": "All products",
            "removed_text": "; ".join(
                participant_display_text(participant)
                for participant in right_participants
            ),
            "masked_text": build_reaction_embedding_text(
                reaction_label,
                left_participants,
                [],
                include_label=include_label,
            ),
        },
    ]

    for index, participant in enumerate(left_participants):
        components.append(
            {
                "key": f"lhs_{index}",
                "kind": "participant",
                "side": "lhs",
                "label": f"Reactant {index + 1}",
                "removed_text": participant_display_text(participant),
                "masked_text": build_reaction_embedding_text(
                    reaction_label,
                    left_participants[:index] + left_participants[index + 1 :],
                    right_participants,
                    include_label=include_label,
                ),
            }
        )
    for index, participant in enumerate(right_participants):
        components.append(
            {
                "key": f"rhs_{index}",
                "kind": "participant",
                "side": "rhs",
                "label": f"Product {index + 1}",
                "removed_text": participant_display_text(participant),
                "masked_text": build_reaction_embedding_text(
                    reaction_label,
                    left_participants,
                    right_participants[:index] + right_participants[index + 1 :],
                    include_label=include_label,
                ),
            }
        )
    return components


def build_reaction_text_embedding_matrix(
    texts: list[str],
    cache_dir: str | Path = "cache",
    use_linkml_store: bool = True,
    embedding_model_name: str = DEFAULT_EMBEDDING_MODEL_NAME,
    n_features: int = DEFAULT_N_FEATURES,
    space_key: str = "reaction_explain",
) -> np.ndarray:
    """Embed arbitrary reaction texts using the configured backend."""
    if use_linkml_store:
        df = pd.DataFrame(
            {
                "rhea_id": [f"EXPLAIN:{index}" for index in range(len(texts))],
                "embedding_text": texts,
            }
        )
        return build_linkml_store_embedding_matrix(
            df,
            cache_dir=Path(cache_dir),
            text_field="embedding_text",
            space_key=space_key,
            embedding_model_name=embedding_model_name,
        )
    return build_text_feature_matrix(texts, n_features=n_features)


def binary_classifier_scores(model: Any, X: np.ndarray) -> np.ndarray:
    """Score a fitted binary model using its natural continuous output."""
    if hasattr(model, "decision_function"):
        return np.asarray(model.decision_function(X), dtype=np.float32)
    if hasattr(model, "predict_proba"):
        return np.asarray(model.predict_proba(X)[:, 1], dtype=np.float32)
    raise TypeError("Model must provide decision_function() or predict_proba().")


def binary_classifier_threshold(model: Any) -> float:
    """Return the positive/negative threshold for one fitted benchmark model."""
    if hasattr(model, "decision_function"):
        return 0.0
    if hasattr(model, "predict_proba"):
        return 0.5
    raise TypeError("Model must provide decision_function() or predict_proba().")


def explain_reaction_class_embedding(
    rhea_id: str,
    reaction_class: str,
    cache_dir: str | Path = "cache",
    model_name: str = "linear_svm",
    feature_space: str = "reaction",
    use_linkml_store: bool = True,
    embedding_model_name: str = DEFAULT_EMBEDDING_MODEL_NAME,
    n_features: int = DEFAULT_N_FEATURES,
    random_state: int = 42,
) -> EmbeddingRuleComparison:
    """Explain one learned class prediction and compare it with the rule output."""
    if feature_space not in {"reaction", "reaction_participants_only"}:
        raise ValueError(
            "Participant-level occlusion explanations currently support only "
            "'reaction' and 'reaction_participants_only'."
        )

    cache_path = Path(cache_dir)
    normalized_id = normalize_rhea_id(rhea_id)
    df = load_rhea_text_df(cache_path)
    row_df = df[df["rhea_id"] == normalized_id]
    if row_df.empty:
        raise KeyError(f"{normalized_id} not found in {cache_path / 'rhea_reactions.jsonl'}")
    row = row_df.iloc[0]

    _, labels_by_class, _ = build_rule_target_data(df)
    if reaction_class not in labels_by_class:
        raise KeyError(f"Unknown benchmark class '{reaction_class}'.")
    y = labels_by_class[reaction_class]
    positive_support = int(y.sum())
    negative_support = int(len(y) - positive_support)
    if positive_support == 0 or negative_support == 0:
        raise ValueError(
            f"Class '{reaction_class}' does not have both positive and negative examples."
        )

    if use_linkml_store:
        base_matrices = build_linkml_store_embedding_matrices(
            df,
            cache_dir=cache_path,
            embedding_model_name=embedding_model_name,
        )
        embedding_backend = "linkml_store"
    else:
        base_matrices = build_lexical_embedding_matrices(
            df,
            n_features=n_features,
        )
        embedding_backend = "lexical"

    feature_matrices = build_rule_benchmark_feature_matrices(base_matrices)
    model = fit_binary_rule_classifier(
        feature_matrices[feature_space],
        y,
        model_name=model_name,
        random_state=random_state,
    )

    include_label = feature_space == "reaction"
    components = build_embedding_explanation_components(
        row["label"],
        row["left_participants"],
        row["right_participants"],
        include_label=include_label,
    )
    full_text = (
        row["embedding_text"]
        if feature_space == "reaction"
        else row["participant_embedding_text"]
    )
    variant_texts = [full_text] + [component["masked_text"] for component in components]
    variant_matrix = build_reaction_text_embedding_matrix(
        variant_texts,
        cache_dir=cache_path,
        use_linkml_store=use_linkml_store,
        embedding_model_name=embedding_model_name,
        n_features=n_features,
        space_key=f"{feature_space}_explain",
    )
    scores = binary_classifier_scores(model, variant_matrix)
    full_score = float(scores[0])
    threshold = binary_classifier_threshold(model)

    component_scores: list[EmbeddingExplanationComponentScore] = []
    for offset, component in enumerate(components, start=1):
        masked_score = float(scores[offset])
        score_delta = full_score - masked_score
        component_scores.append(
            {
                **component,
                "masked_score": masked_score,
                "score_delta": score_delta,
                "supports_class": score_delta > 0.0,
            }
        )

    raw_record = load_rhea_cache_record(normalized_id, cache_dir=cache_path)
    reaction = reaction_from_rhea_cache_record(raw_record)
    classifier = ReactionClassifier()
    if reaction_class not in classifier.reaction_classes:
        raise KeyError(f"Rule classifier '{reaction_class}' is not available.")
    rule_result = classifier.classify_with_class(
        reaction,
        classifier.reaction_classes[reaction_class],
    )

    asserted_support = build_asserted_rule_support(
        row["go_closure_ids"],
        row["ec_numbers"],
        load_rule_class_metadata(),
    )

    return {
        "rhea_id": normalized_id,
        "reaction_label": row["label"],
        "reaction_class": reaction_class,
        "model_name": model_name,
        "feature_space": feature_space,
        "embedding_backend": embedding_backend,
        "embedding_model_name": embedding_model_name,
        "asserted_positive": reaction_class in asserted_support,
        "positive_support": positive_support,
        "negative_support": negative_support,
        "svm_score": full_score,
        "svm_predicted_positive": full_score >= threshold,
        "rule_predicted_positive": rule_result.is_member,
        "rule_explanation": rule_result.explanation,
        "embedding_text": full_text,
        "component_scores": component_scores,
    }
