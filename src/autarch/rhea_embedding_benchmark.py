"""Benchmark autarch rule-class prediction from RHEA embedding spaces."""

from __future__ import annotations

import json
import time
from collections import Counter
from pathlib import Path
from typing import Any, TypedDict, cast

import numpy as np
import pandas as pd

from autarch.rhea_text_embeddings import (
    DEFAULT_EMBEDDING_MODEL_NAME,
    DEFAULT_N_FEATURES,
    build_asserted_rule_support,
    build_lexical_embedding_matrices,
    build_linkml_store_embedding_matrices,
    load_rhea_text_df,
    load_rule_class_metadata,
)


FEATURE_SPACE_ORDER = [
    "reaction",
    "reaction_participants_only",
    "reaction_drfp",
    "rhs_minus_lhs",
    "symmetric_sum_sqdiff",
]
FEATURE_SPACE_LABELS = {
    "reaction": "Reaction",
    "reaction_participants_only": "Reaction participants only",
    "reaction_drfp": "Reaction SMILES (DRFP)",
    "rhs_minus_lhs": "RHS-LHS diff",
    "symmetric_sum_sqdiff": "Symmetric side summary",
}
MODEL_ORDER = [
    "logistic_regression",
    "linear_svm",
    "random_forest",
    "xgboost",
]
MODEL_LABELS = {
    "logistic_regression": "Logistic regression",
    "linear_svm": "Linear SVM",
    "random_forest": "Random forest",
    "xgboost": "XGBoost",
}


class RuleEmbeddingBenchmarkArtifacts(TypedDict):
    """Paths written by the rule-embedding benchmark."""

    output_dir: str
    summary_json: str
    per_class_csv: str
    summary_csv: str


class RuleEmbeddingBenchmarkSummary(TypedDict):
    """Summary payload returned after benchmarking rule prediction."""

    reaction_count: int
    candidate_classes: int
    benchmarked_classes: int
    min_positives: int
    max_classes: int | None
    cv_splits: int
    models: list[str]
    feature_spaces: list[str]
    feature_dimensions: dict[str, int]
    feature_row_counts: dict[str, int]
    embedding_backend: str
    embedding_model_name: str
    duration_seconds: float
    overall_metrics: list[dict[str, Any]]
    versus_reaction: list[dict[str, Any]]
    artifacts: RuleEmbeddingBenchmarkArtifacts


class FeatureMatrixSpec(TypedDict):
    """Feature matrix plus the dataframe rows it covers."""

    matrix: np.ndarray
    row_indices: np.ndarray


def build_rule_benchmark_feature_matrices(
    base_matrices: dict[str, np.ndarray],
) -> dict[str, np.ndarray]:
    """Build benchmark feature spaces from reaction/LHS/RHS embeddings."""
    reaction_matrix = base_matrices["reaction"]
    lhs_matrix = base_matrices["lhs"]
    rhs_matrix = base_matrices["rhs"]
    rhs_minus_lhs = base_matrices["rhs_minus_lhs"]
    symmetric_sum_sqdiff = np.hstack(
        [
            0.5 * (lhs_matrix + rhs_matrix),
            np.square(rhs_matrix - lhs_matrix),
        ]
    )
    return {
        "reaction": reaction_matrix.astype(np.float32, copy=False),
        "reaction_participants_only": base_matrices[
            "reaction_participants_only"
        ].astype(np.float32, copy=False),
        "rhs_minus_lhs": rhs_minus_lhs.astype(np.float32, copy=False),
        "symmetric_sum_sqdiff": symmetric_sum_sqdiff.astype(np.float32, copy=False),
    }


def build_rule_benchmark_feature_specs(
    base_matrices: dict[str, np.ndarray],
    row_indices: np.ndarray | None = None,
) -> dict[str, FeatureMatrixSpec]:
    """Build benchmark feature specs with explicit row coverage."""
    feature_matrices = build_rule_benchmark_feature_matrices(base_matrices)
    if row_indices is None:
        row_indices = np.arange(
            feature_matrices["reaction"].shape[0],
            dtype=np.int32,
        )
    return {
        space_name: {
            "matrix": matrix,
            "row_indices": row_indices,
        }
        for space_name, matrix in feature_matrices.items()
    }


def build_drfp_feature_spec(
    df: pd.DataFrame,
    n_folded_length: int = 2048,
) -> FeatureMatrixSpec:
    """Build a DRFP feature spec from the valid reaction-SMILES subset."""
    try:
        from drfp import DrfpEncoder
    except ImportError as e:
        raise ValueError(
            "drfp is not installed. Add it to the embeddings dependency group "
            "before selecting the reaction_drfp feature space."
        ) from e

    reaction_smiles = df["reaction_smiles"].tolist()
    row_indices = np.array(
        [
            row_index
            for row_index, smiles in enumerate(reaction_smiles)
            if isinstance(smiles, str) and smiles
        ],
        dtype=np.int32,
    )
    if len(row_indices) == 0:
        raise ValueError("No valid reaction SMILES are available for DRFP encoding.")

    valid_reaction_smiles = [cast(str, reaction_smiles[row_index]) for row_index in row_indices]
    fingerprints = DrfpEncoder.encode(
        valid_reaction_smiles,
        n_folded_length=n_folded_length,
    )
    matrix = np.asarray(fingerprints, dtype=np.float32)
    return {
        "matrix": matrix,
        "row_indices": row_indices,
    }


def select_benchmark_classes(
    class_counts: Counter[str],
    min_positives: int,
    max_classes: int | None = None,
) -> list[str]:
    """Select benchmark classes by support count."""
    selected = [
        class_name
        for class_name, count in sorted(
            class_counts.items(),
            key=lambda item: (-item[1], item[0]),
        )
        if count >= min_positives
    ]
    if max_classes is not None:
        return selected[:max_classes]
    return selected


def select_benchmark_values(
    requested: list[str] | None,
    available: list[str],
    value_type: str,
) -> list[str]:
    """Validate and normalize requested model/space selections."""
    if not requested:
        return available.copy()
    invalid = sorted(set(requested) - set(available))
    if invalid:
        raise ValueError(
            f"Unknown {value_type}(s): {', '.join(invalid)}. "
            f"Choose from: {', '.join(available)}."
        )
    return [value for value in available if value in requested]


def build_rule_target_data(
    df: pd.DataFrame,
) -> tuple[list[str], dict[str, np.ndarray], Counter[str]]:
    """Build one-vs-rest rule labels from cached GO/EC support."""
    rule_class_metadata = load_rule_class_metadata()
    rows = df.to_dict(orient="records")
    class_counts: Counter[str] = Counter()
    labels_by_class: dict[str, np.ndarray] = {}

    for row_index, row in enumerate(rows):
        asserted_support = build_asserted_rule_support(
            row["go_closure_ids"],
            row["ec_numbers"],
            rule_class_metadata,
        )
        for class_name in asserted_support:
            class_counts[class_name] += 1
            if class_name not in labels_by_class:
                labels_by_class[class_name] = np.zeros(len(rows), dtype=np.int8)
            labels_by_class[class_name][row_index] = 1

    return [row["rhea_id"] for row in rows], labels_by_class, class_counts


def _build_binary_classifier(
    model_name: str,
    positive_count: int,
    negative_count: int,
    random_state: int = 42,
) -> tuple[Any, str, float]:
    """Build a binary classifier and its score/threshold conventions."""
    if model_name == "logistic_regression":
        from sklearn.linear_model import LogisticRegression

        return (
            LogisticRegression(
                max_iter=3000,
                class_weight="balanced",
                solver="lbfgs",
            ),
            "predict_proba",
            0.5,
        )
    if model_name == "linear_svm":
        from sklearn.svm import LinearSVC

        return (
            LinearSVC(
                class_weight="balanced",
                dual="auto",
                max_iter=10000,
                random_state=random_state,
            ),
            "decision_function",
            0.0,
        )
    if model_name == "random_forest":
        from sklearn.ensemble import RandomForestClassifier

        return (
            RandomForestClassifier(
                class_weight="balanced_subsample",
                max_depth=None,
                n_estimators=300,
                n_jobs=1,
                random_state=random_state,
            ),
            "predict_proba",
            0.5,
        )
    if model_name == "xgboost":
        try:
            from xgboost import XGBClassifier
        except ImportError as e:
            raise ValueError(
                "xgboost is not installed. Add it to the environment before "
                "selecting the xgboost benchmark model."
            ) from e

        scale_pos_weight = negative_count / positive_count
        return (
            XGBClassifier(
                colsample_bytree=0.8,
                eval_metric="logloss",
                learning_rate=0.05,
                max_depth=4,
                n_estimators=300,
                n_jobs=1,
                random_state=random_state,
                reg_lambda=1.0,
                scale_pos_weight=scale_pos_weight,
                subsample=0.8,
                tree_method="hist",
                verbosity=0,
            ),
            "predict_proba",
            0.5,
        )
    raise ValueError(
        f"Unknown benchmark model '{model_name}'. "
        f"Choose from: {', '.join(MODEL_ORDER)}."
    )


def _evaluate_binary_rule_classifier(
    X: np.ndarray,
    y: np.ndarray,
    cv_splits: int,
    model_name: str = "logistic_regression",
    n_jobs: int = -1,
    random_state: int = 42,
) -> dict[str, float]:
    """Evaluate one binary one-vs-rest task with cross-validated predictions."""
    from sklearn.metrics import average_precision_score, f1_score, roc_auc_score
    from sklearn.model_selection import StratifiedKFold, cross_val_predict

    positive_count = int(y.sum())
    negative_count = int(len(y) - positive_count)
    model, score_method, decision_threshold = _build_binary_classifier(
        model_name=model_name,
        positive_count=positive_count,
        negative_count=negative_count,
        random_state=random_state,
    )
    cv = StratifiedKFold(
        n_splits=cv_splits,
        shuffle=True,
        random_state=random_state,
    )
    scores = cross_val_predict(
        model,
        X,
        y,
        cv=cv,
        method=score_method,
        n_jobs=n_jobs,
    )
    if score_method == "predict_proba":
        continuous_scores = scores[:, 1]
    else:
        continuous_scores = scores
    predictions = (continuous_scores >= decision_threshold).astype(np.int8)
    return {
        "f1": float(f1_score(y, predictions, zero_division=0)),
        "average_precision": float(average_precision_score(y, continuous_scores)),
        "roc_auc": float(roc_auc_score(y, continuous_scores)),
    }


def fit_binary_rule_classifier(
    X: np.ndarray,
    y: np.ndarray,
    model_name: str = "logistic_regression",
    random_state: int = 42,
) -> Any:
    """Fit one benchmark binary classifier on the full dataset."""
    positive_count = int(y.sum())
    negative_count = int(len(y) - positive_count)
    model, _, _ = _build_binary_classifier(
        model_name=model_name,
        positive_count=positive_count,
        negative_count=negative_count,
        random_state=random_state,
    )
    model.fit(X, y)
    return model


def benchmark_rule_embedding_spaces(
    cache_dir: str | Path = "cache",
    output_dir: str | Path | None = None,
    min_positives: int = 100,
    max_classes: int | None = 25,
    cv_splits: int = 5,
    model_names: list[str] | None = None,
    feature_spaces: list[str] | None = None,
    n_jobs: int = -1,
    use_linkml_store: bool = True,
    embedding_model_name: str = DEFAULT_EMBEDDING_MODEL_NAME,
    n_features: int = DEFAULT_N_FEATURES,
) -> RuleEmbeddingBenchmarkSummary:
    """Benchmark one-vs-rest rule prediction from cached embedding spaces."""
    cache_path = Path(cache_dir)
    output_path = (
        Path(output_dir)
        if output_dir is not None
        else cache_path / "rhea_rule_embedding_benchmark"
    )
    output_path.mkdir(parents=True, exist_ok=True)

    started_at = time.perf_counter()
    df = load_rhea_text_df(cache_path)
    rhea_ids, labels_by_class, class_counts = build_rule_target_data(df)
    selected_models = select_benchmark_values(
        requested=model_names,
        available=MODEL_ORDER,
        value_type="model",
    )
    selected_classes = select_benchmark_classes(
        class_counts,
        min_positives=min_positives,
        max_classes=max_classes,
    )
    if not selected_classes:
        raise ValueError(
            "No classes meet the benchmark threshold. "
            f"Try lowering min_positives below {min_positives}."
        )

    selected_spaces = select_benchmark_values(
        requested=feature_spaces,
        available=FEATURE_SPACE_ORDER,
        value_type="feature space",
    )

    feature_specs: dict[str, FeatureMatrixSpec] = {}
    embedding_backend_parts: list[str] = []

    requested_text_spaces = [
        space_name for space_name in selected_spaces if space_name != "reaction_drfp"
    ]
    if requested_text_spaces:
        if use_linkml_store:
            base_matrices = build_linkml_store_embedding_matrices(
                df,
                cache_dir=cache_path,
                embedding_model_name=embedding_model_name,
            )
            embedding_backend_parts.append("linkml_store")
        else:
            base_matrices = build_lexical_embedding_matrices(
                df,
                n_features=n_features,
            )
            embedding_backend_parts.append("lexical")
        feature_specs.update(
            build_rule_benchmark_feature_specs(
                base_matrices,
                row_indices=np.arange(len(df), dtype=np.int32),
            )
        )
    if "reaction_drfp" in selected_spaces:
        feature_specs["reaction_drfp"] = build_drfp_feature_spec(df)
        embedding_backend_parts.append("drfp")
    embedding_backend = "+".join(embedding_backend_parts)
    detail_rows: list[dict[str, Any]] = []

    for class_name in selected_classes:
        y = labels_by_class[class_name]
        for model_name in selected_models:
            for space_name in selected_spaces:
                feature_spec = feature_specs[space_name]
                X = feature_spec["matrix"]
                row_indices = feature_spec["row_indices"]
                y_space = y[row_indices]
                positive_count = int(y_space.sum())
                negative_count = int(len(y_space) - positive_count)
                if positive_count < min_positives:
                    continue
                if positive_count < cv_splits or negative_count < cv_splits:
                    continue
                prevalence = positive_count / len(y_space) if len(y_space) else 0.0
                metrics = _evaluate_binary_rule_classifier(
                    X,
                    y_space,
                    cv_splits=cv_splits,
                    model_name=model_name,
                    n_jobs=n_jobs,
                )
                detail_rows.append(
                    {
                        "class": class_name,
                        "model": model_name,
                        "model_label": MODEL_LABELS[model_name],
                        "space": space_name,
                        "space_label": FEATURE_SPACE_LABELS[space_name],
                        "positives": positive_count,
                        "negatives": negative_count,
                        "samples": len(y_space),
                        "prevalence": prevalence,
                        **metrics,
                    }
                )

    if not detail_rows:
        raise ValueError("No class/space evaluations were completed.")

    detail_df = pd.DataFrame(detail_rows).sort_values(
        ["model", "space", "f1", "average_precision", "class"],
        ascending=[True, True, False, False, True],
    )
    summary_df = (
        detail_df.groupby(["model", "model_label", "space", "space_label"], as_index=False)
        .agg(
            classes=("class", "count"),
            mean_f1=("f1", "mean"),
            median_f1=("f1", "median"),
            mean_average_precision=("average_precision", "mean"),
            mean_roc_auc=("roc_auc", "mean"),
            mean_prevalence=("prevalence", "mean"),
        )
        .sort_values("mean_f1", ascending=False)
    )

    comparison_rows: list[dict[str, Any]] = []
    if "reaction" in selected_spaces:
        for model_name in selected_models:
            reaction_scores = (
                detail_df[
                    (detail_df["model"] == model_name) & (detail_df["space"] == "reaction")
                ][["class", "f1"]]
                .rename(columns={"f1": "reaction_f1"})
                .set_index("class")
            )
            for space_name in selected_spaces:
                if space_name == "reaction":
                    continue
                space_scores = (
                    detail_df[
                        (detail_df["model"] == model_name) & (detail_df["space"] == space_name)
                    ][["class", "f1"]]
                    .rename(columns={"f1": "space_f1"})
                    .set_index("class")
                )
                merged = reaction_scores.join(space_scores, how="inner")
                deltas = merged["space_f1"] - merged["reaction_f1"]
                comparison_rows.append(
                    {
                        "model": model_name,
                        "model_label": MODEL_LABELS[model_name],
                        "space": space_name,
                        "space_label": FEATURE_SPACE_LABELS[space_name],
                        "classes_compared": int(len(merged)),
                        "better_than_reaction": int((deltas > 0).sum()),
                        "worse_than_reaction": int((deltas < 0).sum()),
                        "mean_f1_delta": float(deltas.mean()),
                    }
                )

    detail_csv_path = output_path / "per_class_metrics.csv"
    summary_csv_path = output_path / "summary_metrics.csv"
    summary_json_path = output_path / "summary.json"
    detail_df.to_csv(detail_csv_path, index=False)
    summary_df.to_csv(summary_csv_path, index=False)

    duration_seconds = round(time.perf_counter() - started_at, 3)
    overall_metrics = cast(list[dict[str, Any]], summary_df.to_dict(orient="records"))
    summary_payload: RuleEmbeddingBenchmarkSummary = {
        "reaction_count": len(rhea_ids),
        "candidate_classes": len(class_counts),
        "benchmarked_classes": len(detail_df["class"].unique()),
        "min_positives": min_positives,
        "max_classes": max_classes,
        "cv_splits": cv_splits,
        "models": selected_models,
        "feature_spaces": selected_spaces,
        "feature_dimensions": {
            space_name: int(feature_specs[space_name]["matrix"].shape[1])
            for space_name in selected_spaces
        },
        "feature_row_counts": {
            space_name: int(feature_specs[space_name]["matrix"].shape[0])
            for space_name in selected_spaces
        },
        "embedding_backend": embedding_backend,
        "embedding_model_name": embedding_model_name,
        "duration_seconds": duration_seconds,
        "overall_metrics": overall_metrics,
        "versus_reaction": comparison_rows,
        "artifacts": {
            "output_dir": str(output_path),
            "summary_json": str(summary_json_path),
            "per_class_csv": str(detail_csv_path),
            "summary_csv": str(summary_csv_path),
        },
    }
    summary_json_path.write_text(json.dumps(summary_payload, indent=2))
    return summary_payload
