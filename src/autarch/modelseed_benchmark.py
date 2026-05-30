"""ModelSEED benchmark materialization and cached prediction utilities."""

from __future__ import annotations

import json
import time
from collections import Counter
from pathlib import Path
from typing import Optional, TypedDict

from pydantic import BaseModel, Field

from autarch.classifier import ReactionClassifier
from autarch.datamodel import ClassificationResult, Reaction
from autarch.etl.modelseed_etl import (
    ModelSeedCacheSummary,
    ModelSeedReactionRecord,
    cache_modelseed_dataset,
    load_modelseed_reaction_records,
    serialize_modelseed_reactions,
    summarize_modelseed_cache,
)
from autarch.ontology.reaction import ReactionClass


class ModelSeedBenchmarkFilters(TypedDict):
    """Filter settings used to materialize one benchmark slice."""

    include_rhea: bool
    include_transport: bool
    all_statuses: bool
    allow_partial_mapping: bool
    require_ec: bool
    limit: Optional[int]


class ModelSeedBenchmarkArtifacts(TypedDict):
    """Paths to benchmark artifacts written to disk."""

    output_dir: str
    candidates: str
    predictions: str
    summary: str


class ModelSeedBenchmarkSample(TypedDict):
    """Compact sample of one positive benchmark prediction."""

    modelseed_id: str
    name: str
    status: str
    ec_numbers: list[str]
    positive_classes: list[str]


class ModelSeedBenchmarkSummary(TypedDict):
    """Typed summary returned after benchmark materialization."""

    input_reactions: int
    upstream_total_reactions: int
    classifiers_run: int
    selected_reactions: int
    selected_with_ec: int
    selected_without_ec: int
    positive_reactions: int
    positive_rate: float
    duration_seconds: float
    filters: ModelSeedBenchmarkFilters
    skipped_counts: dict[str, int]
    selected_status_counts: list[tuple[str, int]]
    top_positive_classes: list[tuple[str, int]]
    sample_positive_predictions: list[ModelSeedBenchmarkSample]
    artifacts: ModelSeedBenchmarkArtifacts


class ModelSeedBenchmarkPrediction(BaseModel):
    """Cached classifier output for one ModelSEED reaction."""

    modelseed_id: str
    positive_classes: list[str] = Field(default_factory=list)
    positive_explanations: dict[str, str] = Field(default_factory=dict)


def serialize_modelseed_benchmark_predictions(
    predictions: list[ModelSeedBenchmarkPrediction],
    output_file: str | Path,
) -> None:
    """Serialize ModelSEED benchmark predictions to JSONL."""
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as handle:
        for prediction in predictions:
            json.dump(prediction.model_dump(), handle)
            handle.write("\n")


def reaction_has_fractional_stoichiometry(reaction: Optional[Reaction]) -> bool:
    """Return True when a cached reaction contains numeric fractional stoichiometry."""
    if reaction is None:
        return False

    for participant in reaction.all_participants():
        if not participant.stoichiometry:
            continue
        try:
            value = float(participant.stoichiometry)
        except ValueError:
            continue
        if not value.is_integer():
            return True

    return False


def _ensure_modelseed_cache(cache_dir: Path) -> ModelSeedCacheSummary:
    """Ensure the transformed ModelSEED cache exists and is not limited."""
    reactions_path = cache_dir / "modelseed_reactions.jsonl"

    if not reactions_path.exists():
        cache_modelseed_dataset(cache_dir=cache_dir)
        return summarize_modelseed_cache(cache_dir)

    summary = summarize_modelseed_cache(cache_dir)
    if summary["cached_reactions"] < summary["total_reactions"]:
        cache_modelseed_dataset(cache_dir=cache_dir)
        summary = summarize_modelseed_cache(cache_dir)

    return summary


def _classify_reaction(
    reaction: Reaction,
    classifier_instances: list[tuple[str, ReactionClass]],
) -> dict[str, ClassificationResult]:
    """Classify one reaction using pre-instantiated classifiers."""
    results: dict[str, ClassificationResult] = {}
    for class_name, classifier_instance in classifier_instances:
        results[class_name] = classifier_instance.check_membership(reaction)
    return results


def build_modelseed_benchmark(
    cache_dir: str | Path = "cache",
    output_dir: str | Path | None = None,
    limit: Optional[int] = None,
    include_rhea: bool = False,
    include_transport: bool = False,
    all_statuses: bool = False,
    allow_partial_mapping: bool = False,
    require_ec: bool = False,
) -> ModelSeedBenchmarkSummary:
    """Materialize a filtered ModelSEED benchmark slice and cached predictions."""
    cache_path = Path(cache_dir)
    summary = _ensure_modelseed_cache(cache_path)

    reactions = load_modelseed_reaction_records(cache_path / "modelseed_reactions.jsonl")
    if output_dir is not None:
        output_path = Path(output_dir)
    elif require_ec:
        output_path = cache_path / "modelseed_benchmark_ec"
    else:
        output_path = cache_path / "modelseed_benchmark"
    output_path.mkdir(parents=True, exist_ok=True)

    classifier = ReactionClassifier()
    classifier_instances = [
        (class_name, classifier.reaction_classes[class_name]())
        for class_name in sorted(classifier.reaction_classes)
    ]

    selected_reactions: dict[str, ModelSeedReactionRecord] = {}
    predictions: list[ModelSeedBenchmarkPrediction] = []
    skipped_counts: Counter[str] = Counter()
    selected_status_counts: Counter[str] = Counter()
    positive_class_counts: Counter[str] = Counter()
    positive_reactions = 0
    selected_with_ec = 0
    sample_positive_predictions: list[ModelSeedBenchmarkSample] = []

    started_at = time.perf_counter()

    for reaction_record in reactions.values():
        if not include_rhea and reaction_record.rhea_ids:
            skipped_counts["rhea_alias"] += 1
            continue
        if not all_statuses and reaction_record.status != "OK":
            skipped_counts["status"] += 1
            continue
        if not include_transport and reaction_record.is_transport:
            skipped_counts["transport"] += 1
            continue
        if not allow_partial_mapping and not reaction_record.fully_mapped_to_chebi:
            skipped_counts["partial_mapping"] += 1
            continue
        if reaction_has_fractional_stoichiometry(reaction_record.reaction):
            skipped_counts["fractional_stoichiometry"] += 1
            continue
        if require_ec and not reaction_record.ec_numbers:
            skipped_counts["missing_ec"] += 1
            continue
        if reaction_record.reaction is None:
            skipped_counts["missing_reaction"] += 1
            continue

        selected_reactions[reaction_record.modelseed_id] = reaction_record
        selected_status_counts[reaction_record.status or "<blank>"] += 1
        if reaction_record.ec_numbers:
            selected_with_ec += 1

        classification_results = _classify_reaction(
            reaction_record.reaction,
            classifier_instances,
        )
        positive_explanations = {
            class_name: result.explanation
            for class_name, result in classification_results.items()
            if result.is_member
        }
        positive_classes = sorted(positive_explanations)

        predictions.append(
            ModelSeedBenchmarkPrediction(
                modelseed_id=reaction_record.modelseed_id,
                positive_classes=positive_classes,
                positive_explanations=positive_explanations,
            )
        )

        if positive_classes:
            positive_reactions += 1
            positive_class_counts.update(positive_classes)
            if len(sample_positive_predictions) < 20:
                sample_positive_predictions.append(
                    {
                        "modelseed_id": reaction_record.modelseed_id,
                        "name": reaction_record.name,
                        "status": reaction_record.status,
                        "ec_numbers": reaction_record.ec_numbers,
                        "positive_classes": positive_classes,
                    }
                )

        if limit is not None and len(selected_reactions) >= limit:
            break

    duration_seconds = round(time.perf_counter() - started_at, 3)

    candidates_path = output_path / "candidates.jsonl"
    predictions_path = output_path / "predictions.jsonl"
    summary_path = output_path / "summary.json"

    serialize_modelseed_reactions(selected_reactions, candidates_path)
    serialize_modelseed_benchmark_predictions(predictions, predictions_path)

    selected_count = len(selected_reactions)
    summary_payload: ModelSeedBenchmarkSummary = {
        "input_reactions": len(reactions),
        "upstream_total_reactions": summary["total_reactions"],
        "classifiers_run": len(classifier_instances),
        "selected_reactions": selected_count,
        "selected_with_ec": selected_with_ec,
        "selected_without_ec": selected_count - selected_with_ec,
        "positive_reactions": positive_reactions,
        "positive_rate": round(positive_reactions / selected_count, 4)
        if selected_count
        else 0.0,
        "duration_seconds": duration_seconds,
        "filters": {
            "include_rhea": include_rhea,
            "include_transport": include_transport,
            "all_statuses": all_statuses,
            "allow_partial_mapping": allow_partial_mapping,
            "require_ec": require_ec,
            "limit": limit,
        },
        "skipped_counts": dict(sorted(skipped_counts.items())),
        "selected_status_counts": selected_status_counts.most_common(),
        "top_positive_classes": positive_class_counts.most_common(25),
        "sample_positive_predictions": sample_positive_predictions,
        "artifacts": {
            "output_dir": str(output_path),
            "candidates": str(candidates_path),
            "predictions": str(predictions_path),
            "summary": str(summary_path),
        },
    }

    with open(summary_path, "w") as handle:
        json.dump(summary_payload, handle, indent=2)

    return summary_payload
