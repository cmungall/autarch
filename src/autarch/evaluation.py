"""Evaluation metrics and utilities for reaction classification."""

from dataclasses import dataclass
from typing import Optional
import json
from pathlib import Path

from autarch.datamodel import Reaction, Participant, GoTerm
from autarch.classifier import ReactionClassifier
from autarch.moiety import save_moiety_cache


def ec_matches_prefix(ec_number: str, ec_prefix: str) -> bool:
    """Check if an EC number matches a prefix pattern.

    Handles wildcard patterns like "3.-.-.-" or "3.5.-.-".

    Args:
        ec_number: Full EC number like "3.5.1.50"
        ec_prefix: Prefix pattern like "3.-.-.-" or "3.5.-.-"

    Returns:
        True if the EC number matches the prefix pattern

    Examples:
        >>> ec_matches_prefix("3.5.1.50", "3.-.-.-")
        True
        >>> ec_matches_prefix("3.5.1.50", "3.5.-.-")
        True
        >>> ec_matches_prefix("3.5.1.50", "3.4.-.-")
        False
        >>> ec_matches_prefix("1.1.1.1", "3.-.-.-")
        False
    """
    number_parts = ec_number.split(".")
    prefix_parts = ec_prefix.split(".")

    for number_part, prefix_part in zip(number_parts, prefix_parts):
        if prefix_part == "-":
            return True
        if number_part != prefix_part:
            return False

    return len(number_parts) == len(prefix_parts)

# Suppress RDKit warnings about hydrogen atoms and other parsing issues
try:
    from rdkit import RDLogger  # type: ignore

    RDLogger.DisableLog("rdApp.*")  # type: ignore
except Exception:
    pass  # RDLogger might not have DisableLog in all versions


@dataclass
class EvaluationMetrics:
    """Metrics for evaluating classification performance."""

    true_positives: int = 0
    false_positives: int = 0
    true_negatives: int = 0
    false_negatives: int = 0

    # Lists to track individual predictions
    tp_reactions: Optional[list[str]] = None
    fp_reactions: Optional[list[str]] = None
    tn_reactions: Optional[list[str]] = None
    fn_reactions: Optional[list[str]] = None

    # Store explanations for debugging
    tp_explanations: Optional[dict[str, str]] = None
    fp_explanations: Optional[dict[str, str]] = None
    tn_explanations: Optional[dict[str, str]] = None
    fn_explanations: Optional[dict[str, str]] = None
    
    # Track unlabeled reactions (no GO/EC annotations)
    unlabeled_predicted_positive: Optional[list[tuple[str, str]]] = None
    unlabeled_predicted_negative: Optional[list[tuple[str, str]]] = None

    def __post_init__(self):
        """Initialize reaction lists and explanation dicts if not provided."""
        if self.tp_reactions is None:
            self.tp_reactions = []
        if self.fp_reactions is None:
            self.fp_reactions = []
        if self.tn_reactions is None:
            self.tn_reactions = []
        if self.fn_reactions is None:
            self.fn_reactions = []

        if self.tp_explanations is None:
            self.tp_explanations = {}
        if self.fp_explanations is None:
            self.fp_explanations = {}
        if self.tn_explanations is None:
            self.tn_explanations = {}
        if self.fn_explanations is None:
            self.fn_explanations = {}

    @property
    def total(self) -> int:
        """Total number of predictions."""
        return (
            self.true_positives
            + self.false_positives
            + self.true_negatives
            + self.false_negatives
        )

    @property
    def precision(self) -> float:
        """Precision: TP / (TP + FP)."""
        denominator = self.true_positives + self.false_positives
        if denominator == 0:
            return 0.0
        return self.true_positives / denominator

    @property
    def recall(self) -> float:
        """Recall (sensitivity): TP / (TP + FN)."""
        denominator = self.true_positives + self.false_negatives
        if denominator == 0:
            return 0.0
        return self.true_positives / denominator

    @property
    def specificity(self) -> float:
        """Specificity: TN / (TN + FP)."""
        denominator = self.true_negatives + self.false_positives
        if denominator == 0:
            return 0.0
        return self.true_negatives / denominator

    @property
    def accuracy(self) -> float:
        """Accuracy: (TP + TN) / total."""
        if self.total == 0:
            return 0.0
        return (self.true_positives + self.true_negatives) / self.total

    @property
    def f1_score(self) -> float:
        """F1 score: harmonic mean of precision and recall."""
        if self.precision + self.recall == 0:
            return 0.0
        return 2 * (self.precision * self.recall) / (self.precision + self.recall)

    @property
    def matthews_correlation_coefficient(self) -> float:
        """Matthews correlation coefficient for binary classification."""
        tp, tn, fp, fn = (
            self.true_positives,
            self.true_negatives,
            self.false_positives,
            self.false_negatives,
        )

        denominator = ((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) ** 0.5
        if denominator == 0:
            return 0.0

        return (tp * tn - fp * fn) / denominator


def evaluate_reaction_class(
    reaction_class: str,
    go_term: Optional[str] = None,
    ec_prefix: Optional[str] = None,
    cache_dir: str = "cache",
    verbose: bool = False,
    go_only: bool = True,
) -> EvaluationMetrics:
    """Evaluate a reaction class against cached GO and RHEA data.

    Args:
        reaction_class: Name of the reaction class to evaluate
        go_term: GO term associated with this reaction class
        ec_prefix: EC number prefix (e.g., "3.1.1" for esterases)
        cache_dir: Directory containing cached data
        verbose: Whether to print detailed information
        go_only: If True, only evaluate reactions that have GO term mappings.
            This ensures we have proper ground truth for evaluation.
            If False, also use EC numbers as ground truth (may inflate metrics).

    Returns:
        EvaluationMetrics with classification results

    Examples:
        >>> # This would evaluate Hydrolase if cache files existed
        >>> # metrics = evaluate_reaction_class("Hydrolase", "GO:0016787", "cache")
        >>> # isinstance(metrics, EvaluationMetrics)
        >>> # True
        >>>
        >>> # For testing, just check the function exists
        >>> callable(evaluate_reaction_class)
        True
    """
    cache_path = Path(cache_dir)
    go_cache = cache_path / "go_terms.jsonl"
    rhea_cache = cache_path / "rhea_reactions.jsonl"
    chebi_smiles_cache = cache_path / "chebi_smiles.json"

    if not go_cache.exists() or not rhea_cache.exists():
        raise FileNotFoundError(f"Cache files not found in {cache_dir}")

    if not chebi_smiles_cache.exists():
        raise FileNotFoundError(
            f"CHEBI SMILES cache not found at {chebi_smiles_cache}. "
            f"Run 'autarch cache-chebi --cache-dir {cache_dir}' first."
        )

    # Load CHEBI to SMILES mapping
    with open(chebi_smiles_cache) as f:
        chebi_to_smiles = json.load(f)

    # Load GO terms (for ancestor information)
    go_terms = {}
    with open(go_cache) as f:
        for line in f:
            term_data = json.loads(line)
            term = GoTerm(**term_data)
            go_terms[term.go_id] = term

    if verbose:
        print(f"Loaded {len(go_terms)} GO terms")

    # Load RHEA reactions (now with go_terms directly in the data)
    rhea_reactions = {}
    with open(rhea_cache) as f:
        for line in f:
            reaction_data = json.loads(line)
            rhea_reactions[reaction_data["rhea_id"]] = reaction_data

    # If go_only mode, filter to only reactions with GO mappings
    if go_only:
        original_count = len(rhea_reactions)
        rhea_reactions = {
            rhea_id: data
            for rhea_id, data in rhea_reactions.items()
            if data.get("go_terms")  # Has direct GO mappings
        }
        if verbose:
            print(f"GO-only mode: filtered {original_count} -> {len(rhea_reactions)} reactions with GO mappings")

    # Initialize classifier registry and target class once.
    classifier = ReactionClassifier()
    if reaction_class not in classifier.reaction_classes:
        raise ValueError(f"Unknown reaction class: {reaction_class}")
    reaction_class_obj = classifier.reaction_classes.get(reaction_class)
    if reaction_class_obj is None:
        raise ValueError(f"Unknown reaction class: {reaction_class}")
    reaction_class_instance = reaction_class_obj()

    # If no GO term or EC prefix provided, try to get them from the class
    if not go_term and reaction_class_obj and hasattr(reaction_class_obj, "GO_ID"):
        go_term = reaction_class_obj.GO_ID
        if verbose:
            print(f"Using GO_ID from class: {go_term}")
    
    if not ec_prefix and reaction_class_obj and hasattr(reaction_class_obj, "EC_NUMBER_PREFIX"):
        ec_prefix = reaction_class_obj.EC_NUMBER_PREFIX
        if verbose:
            print(f"Using EC_NUMBER_PREFIX from class: {ec_prefix}")

    # Check for specific EC_NUMBERS list (for classifiers that target specific EC numbers)
    ec_numbers_list = None
    if reaction_class_obj and hasattr(reaction_class_obj, "EC_NUMBERS"):
        ec_numbers_list = reaction_class_obj.EC_NUMBERS
        if verbose:
            print(f"Using EC_NUMBERS from class: {ec_numbers_list}")

    # Helper function to check if a GO term is or is a descendant of the target
    def is_positive_go_term(reaction_go_term: str, target_go_term: str) -> bool:
        """Check if reaction_go_term matches target or is a descendant of target."""
        if reaction_go_term == target_go_term:
            return True
        if reaction_go_term in go_terms:
            # Check if target is an ancestor of this reaction's GO term
            return target_go_term in go_terms[reaction_go_term].ancestors
        return False

    # Find RHEA reactions that should be positive (belong to this GO term or EC prefix)
    positive_rhea_ids = set()
    negative_rhea_ids = set()
    unlabeled_rhea_ids = set()  # Reactions without any GO/EC annotations

    # First, check specific EC numbers list if provided (skip in go_only mode)
    if ec_numbers_list and not go_only:
        for rhea_id, rhea_data in rhea_reactions.items():
            reaction_ec_numbers = rhea_data.get("ec_numbers", [])
            for ec in reaction_ec_numbers:
                if ec in ec_numbers_list:
                    positive_rhea_ids.add(rhea_id)
                    break

        if verbose and positive_rhea_ids:
            print(f"Found {len(positive_rhea_ids)} positive RHEA IDs for EC numbers {ec_numbers_list}")

    # Then check EC prefix if provided (and no specific list, skip in go_only mode)
    if ec_prefix and not ec_numbers_list and not go_only:
        for rhea_id, rhea_data in rhea_reactions.items():
            ec_numbers = rhea_data.get("ec_numbers", [])
            for ec in ec_numbers:
                if ec_matches_prefix(ec, ec_prefix):
                    positive_rhea_ids.add(rhea_id)
                    break

        if verbose and positive_rhea_ids:
            print(f"Found {len(positive_rhea_ids)} positive RHEA IDs for EC prefix {ec_prefix}")

    # Use GO mappings from RHEA cache directly (much simpler!)
    if go_term:
        for rhea_id, rhea_data in rhea_reactions.items():
            # Check if any of the reaction's GO terms match the target
            for reaction_go_term in rhea_data.get("go_terms", []):
                if is_positive_go_term(reaction_go_term, go_term):
                    positive_rhea_ids.add(rhea_id)
                    break

        if verbose and positive_rhea_ids:
            print(
                f"Found {len(positive_rhea_ids)} positive RHEA IDs for {go_term} and its descendants"
            )

    # Categorize remaining reactions as negative (has GO/EC but not matching) or unlabeled
    for rhea_id, rhea_data in rhea_reactions.items():
        if rhea_id not in positive_rhea_ids:
            # Check if this reaction has any GO or EC annotations
            has_go = bool(rhea_data.get("go_terms", []))
            has_ec = bool(rhea_data.get("ec_numbers", []))

            if has_go or has_ec:
                # Has annotations but not for our target class - true negative
                negative_rhea_ids.add(rhea_id)
            else:
                # No annotations - we don't know the ground truth
                unlabeled_rhea_ids.add(rhea_id)

    # Evaluate each RHEA reaction
    metrics = EvaluationMetrics()
    skipped_insufficient_evidence = 0
    
    # Track predictions on unlabeled data separately
    unlabeled_predicted_positive = []
    unlabeled_predicted_negative = []

    for rhea_id, rhea_data in rhea_reactions.items():
        # Handle both old and new format
        # New format: reaction object with left_participants/right_participants
        # Old format: direct inputs/outputs arrays

        left_participants = []
        right_participants = []

        if "reaction" in rhea_data and rhea_data["reaction"]:
            # New format from SPARQL ETL
            reaction_obj = rhea_data["reaction"]

            # Extract participants from the reaction object
            for p in reaction_obj.get("left_participants", []):
                if p:
                    chebi_id = p.get("chebi_id")
                    # Use SMILES from reaction object, fall back to ChEBI cache
                    smiles = p.get("smiles") or (chebi_to_smiles.get(chebi_id) if chebi_id else None)
                    inchi = p.get("inchi")  # Include InChI for stereochemistry detection
                    location = p.get("location")  # Include location for transport detection
                    name = p.get("name")  # Include name for better debugging
                    polymer_index = p.get("polymer_index")  # For polymer reactions
                    polymer_type = p.get("polymer_type")  # For polymer type classification
                    # Include ALL participants, even without ChEBI IDs
                    left_participants.append(
                        Participant(
                            chebi_id=chebi_id, smiles=smiles, inchi=inchi, location=location,
                            name=name, polymer_index=polymer_index, polymer_type=polymer_type
                        )
                    )

            for p in reaction_obj.get("right_participants", []):
                if p:
                    chebi_id = p.get("chebi_id")
                    # Use SMILES from reaction object, fall back to ChEBI cache
                    smiles = p.get("smiles") or (chebi_to_smiles.get(chebi_id) if chebi_id else None)
                    inchi = p.get("inchi")  # Include InChI for stereochemistry detection
                    location = p.get("location")  # Include location for transport detection
                    name = p.get("name")  # Include name for better debugging
                    polymer_index = p.get("polymer_index")  # For polymer reactions
                    polymer_type = p.get("polymer_type")  # For polymer type classification
                    # Include ALL participants, even without ChEBI IDs
                    right_participants.append(
                        Participant(
                            chebi_id=chebi_id, smiles=smiles, inchi=inchi, location=location,
                            name=name, polymer_index=polymer_index, polymer_type=polymer_type
                        )
                    )
        else:
            # Old format - backward compatibility
            inputs = rhea_data.get("inputs", [])
            outputs = rhea_data.get("outputs", [])

            if not inputs or not outputs:
                continue

            for chebi_id in inputs:
                smiles = chebi_to_smiles.get(chebi_id)
                if smiles:
                    left_participants.append(
                        Participant(chebi_id=chebi_id, smiles=smiles)
                    )
                else:
                    left_participants.append(Participant(chebi_id=chebi_id))

            for chebi_id in outputs:
                smiles = chebi_to_smiles.get(chebi_id)
                if smiles:
                    right_participants.append(
                        Participant(chebi_id=chebi_id, smiles=smiles)
                    )
                else:
                    right_participants.append(Participant(chebi_id=chebi_id))

        # Skip if no participants found
        if not left_participants or not right_participants:
            continue

        # Get the human-readable label from RHEA for pattern matching
        rhea_label = rhea_data.get("label", "")

        reaction = Reaction(
            left_participants=left_participants, right_participants=right_participants,
            label=rhea_label
        )

        if not reaction_class_obj.supports_evaluation(reaction):
            skipped_insufficient_evidence += 1
            continue

        # Evaluate only the target classifier. This is faster and respects
        # classifier-specific evidence requirements.
        classification_result = reaction_class_instance.check_membership(reaction)
        is_predicted_positive = classification_result.is_member
        explanation = classification_result.explanation

        # Check if this reaction has a known label
        if rhea_id in unlabeled_rhea_ids:
            # This reaction has no GO/EC annotations - track separately
            if is_predicted_positive:
                unlabeled_predicted_positive.append((rhea_id, explanation))
            else:
                unlabeled_predicted_negative.append((rhea_id, explanation))
        else:
            # This reaction has annotations - include in metrics
            is_actually_positive = rhea_id in positive_rhea_ids
            
            # Update metrics and store explanations
            if is_actually_positive and is_predicted_positive:
                metrics.true_positives += 1
                assert metrics.tp_reactions is not None
                assert metrics.tp_explanations is not None
                metrics.tp_reactions.append(rhea_id)
                metrics.tp_explanations[rhea_id] = explanation
            elif is_actually_positive and not is_predicted_positive:
                metrics.false_negatives += 1
                assert metrics.fn_reactions is not None
                assert metrics.fn_explanations is not None
                metrics.fn_reactions.append(rhea_id)
                metrics.fn_explanations[rhea_id] = explanation
            elif not is_actually_positive and is_predicted_positive:
                metrics.false_positives += 1
                assert metrics.fp_reactions is not None
                assert metrics.fp_explanations is not None
                metrics.fp_reactions.append(rhea_id)
                metrics.fp_explanations[rhea_id] = explanation
            else:  # not is_actually_positive and not is_predicted_positive
                metrics.true_negatives += 1
                assert metrics.tn_reactions is not None
                assert metrics.tn_explanations is not None
                metrics.tn_reactions.append(rhea_id)
                metrics.tn_explanations[rhea_id] = explanation

    if verbose:
        if skipped_insufficient_evidence > 0:
            print(
                "Warning: Skipped "
                f"{skipped_insufficient_evidence} reactions without enough evidence "
                f"for {reaction_class}"
            )
        
        # Report on unlabeled reactions
        if unlabeled_predicted_positive or unlabeled_predicted_negative:
            print("\n=== Unlabeled Reactions (no GO/EC annotations) ===")
            print(f"Predicted as {reaction_class}: {len(unlabeled_predicted_positive)}")
            print(f"Predicted as not {reaction_class}: {len(unlabeled_predicted_negative)}")
            
            if unlabeled_predicted_positive:
                print(f"\nTop 10 unlabeled reactions predicted as {reaction_class}:")
                for rhea_id, explanation in unlabeled_predicted_positive[:10]:
                    rhea_data = rhea_reactions.get(rhea_id, {})
                    label = rhea_data.get("label", "Unknown")
                    print(f"  • {rhea_id}: {label}")
                    print(f"    → {explanation}")

    if metrics.total == 0:
        raise RuntimeError(
            f"No labeled reactions could be evaluated. "
            f"{skipped_insufficient_evidence} reactions were skipped due to insufficient "
            f"evaluation evidence for {reaction_class}. "
            f"{len(unlabeled_rhea_ids)} reactions had no GO/EC annotations to determine ground truth. "
            f"Ensure CHEBI SMILES are available in {cache_dir}/chebi_smiles.json"
        )

    # Store unlabeled predictions in metrics for external access if needed
    metrics.unlabeled_predicted_positive = unlabeled_predicted_positive
    metrics.unlabeled_predicted_negative = unlabeled_predicted_negative

    # Save moiety cache (incremental - only saves if modified)
    save_moiety_cache()

    return metrics
