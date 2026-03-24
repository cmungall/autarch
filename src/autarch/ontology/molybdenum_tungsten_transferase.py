"""Molybdenum/tungsten transferase reaction classification.

Transferases (EC 2.10) that transfer molybdenum- or tungsten-containing
groups. The primary reaction is molybdopterin molybdotransferase, which
inserts molybdenum into the molybdopterin cofactor.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.transferase import Transferase


class MolybdenumTungstenTransferase(Transferase):
    """Catalysis of the reaction adenylyl-molybdopterin + molybdate = molybdenum cofactor + AMP."""

    GO_ID = "GO:0061599"  # molybdopterin molybdotransferase activity
    EC_NUMBER_PREFIX = "2.10.-.-"  # Transferring molybdenum or tungsten

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a molybdenum/tungsten transferase.

        Strategy:
        1. Look for molybdenum/tungsten indicators in labels
        2. Look for molybdopterin-related ChEBI participants
        """
        label_lower = reaction.label.lower() if reaction.label else ""

        # Look for molybdenum/tungsten indicators in label
        metal_indicators = [
            "molybdopterin",
            "molybdenum",
            "molybdo",
            "molybdate",
            "tungsten",
            "tungsto",
        ]

        for indicator in metal_indicators:
            if indicator in label_lower:
                return ClassificationResult(
                    is_member=True,
                    explanation=f"Molybdenum/tungsten transferase: label contains '{indicator}'",
                )

        # Check for known molybdopterin-related ChEBI IDs
        molybdopterin_ids = {
            "CHEBI:62727",  # molybdopterin adenine dinucleotide
            "CHEBI:36264",  # molybdate
            "CHEBI:44733",  # molybdopterin
            "CHEBI:71302",  # molybdenum cofactor
        }

        has_molybdopterin = any(
            p.chebi_id in molybdopterin_ids
            for p in reaction.left_participants + reaction.right_participants
        )

        if has_molybdopterin:
            return ClassificationResult(
                is_member=True,
                explanation="Molybdenum/tungsten transferase: molybdopterin participant detected",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No molybdenum/tungsten transfer indicators detected",
        )
