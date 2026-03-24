"""RNA polymerase reaction classification.

RNA polymerases catalyze template-directed RNA synthesis:
    RNA(n) + NTP -> RNA(n+1) + PPi
    RNA(n) + NDP <=> RNA(n+1) + Pi

This classifier handles polymer reactions using (n) notation because these
reactions are typically not represented as fully resolved SMILES.

EC numbers covered:
- 2.7.7.6: DNA-directed RNA polymerase
- 2.7.7.8: polynucleotide phosphorylase
- 2.7.7.48: RNA-directed RNA polymerase
"""

from typing import Optional

from autarch.datamodel import ClassificationResult, Reaction, PolymerType
from autarch.ontology.reaction import ReactionClass
from autarch.molecules import (
    CHEBI_ATP, CHEBI_GTP, CHEBI_CTP, CHEBI_UTP,
    CHEBI_ADP, CHEBI_GDP, CHEBI_CDP, CHEBI_UDP,
    CHEBI_PHOSPHATE, CHEBI_DIPHOSPHATE,
)

# RNA nucleotide phosphates used by RNA polymerases
RNA_NUCLEOTIDE_CHEBI = {
    CHEBI_ATP,   # ATP
    CHEBI_GTP,   # GTP
    CHEBI_CTP,   # CTP
    CHEBI_UTP,   # UTP
    CHEBI_ADP,   # ADP
    CHEBI_GDP,   # GDP
    CHEBI_CDP,   # CDP
    CHEBI_UDP,   # UDP
}


class RNAPolymerase(ReactionClass):
    """Classifier for RNA polymerase reactions.

    RNA polymerases synthesize RNA by adding ribonucleotides:
        RNA(n) + NTP -> RNA(n+1) + PPi (diphosphate)
        RNA(n) + NDP <=> RNA(n+1) + Pi (phosphate)

    Key features detected:
    1. RNA polymer on both sides with growth (n -> n+1)
    2. RNA nucleotide substrate (NTP or NDP)
    3. Pyrophosphate or phosphate product
    """

    GO_ID = "GO:0097747"  # RNA polymerase activity
    EC_NUMBERS = ["2.7.7.6", "2.7.7.8", "2.7.7.48"]
    EC_NUMBER_PREFIX = None  # Don't use prefix matching

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is catalyzed by an RNA polymerase.

        Detects:
        - RNA polymer growth (n -> n+1 or n+1 -> n for reverse)
        - RNA nucleotide triphosphate/diphosphate substrate
        - Pyrophosphate/phosphate product
        """
        # Use RHEA label for pattern matching
        label_lower = reaction.label.lower() if reaction.label else ""

        # Check for RNA polymer participants
        left_polymer = self._find_rna_polymer(reaction.left_participants, label_lower)
        right_polymer = self._find_rna_polymer(reaction.right_participants, label_lower)

        if not left_polymer or not right_polymer:
            return ClassificationResult(
                is_member=False,
                explanation="No RNA polymer found on both sides"
            )

        # Check for polymer growth pattern
        left_idx = left_polymer.polymer_index
        right_idx = right_polymer.polymer_index
        # Determine polymer name from label or polymer_type
        if "rna" in label_lower:
            polymer_name = "RNA"
        elif left_polymer.polymer_type:
            polymer_name = str(left_polymer.polymer_type.value)
        else:
            polymer_name = "RNA"

        if not left_idx or not right_idx:
            return ClassificationResult(
                is_member=False,
                explanation="Polymer missing (n) index notation"
            )

        # Determine growth direction
        growth = self._check_polymer_growth(left_idx, right_idx)

        if growth == "forward":
            # n → n+1: synthesis
            explanation = f"RNA polymerase: {polymer_name}({left_idx}) + nucleotide -> {polymer_name}({right_idx})"
        elif growth == "reverse":
            # n+1 → n: degradation (reverse reaction)
            explanation = f"RNA polymerase (reverse): {polymer_name}({left_idx}) -> {polymer_name}({right_idx}) + nucleotide"
        else:
            return ClassificationResult(
                is_member=False,
                explanation=f"No polymer growth pattern: {left_idx} vs {right_idx}"
            )

        # Check for nucleotide substrate/product
        has_nucleotide = self._has_nucleotide(reaction)
        if not has_nucleotide:
            return ClassificationResult(
                is_member=False,
                explanation="No RNA nucleotide (NTP/NDP) found"
            )

        # Check for phosphate/pyrophosphate
        has_phosphate = self._has_phosphate_product(reaction)
        if not has_phosphate:
            return ClassificationResult(
                is_member=False,
                explanation="No phosphate/diphosphate product found"
            )

        return ClassificationResult(
            is_member=True,
            explanation=explanation
        )

    def _find_rna_polymer(self, participants, label_lower: str = ""):
        """Find RNA polymer participant."""
        for p in participants:
            # Check by polymer_type (primary)
            if p.polymer_type in (
                PolymerType.RNA,
                PolymerType.TRNA,
                PolymerType.MRNA,
                PolymerType.RRNA,
            ):
                return p
            # Name-based fallback when polymer_type is missing
            if p.name and "rna" in p.name.lower():
                return p
            # Check for polymer with index (n) notation
            if p.polymer_index and p.name and "rna" in p.name.lower():
                return p

        # Label-based fallback for RNA reactions
        if "rna" in label_lower:
            for p in participants:
                if p.polymer_index:
                    return p
        return None

    def _check_polymer_growth(
        self, left_idx: str, right_idx: str
    ) -> Optional[str]:
        """Check if polymer grows from left to right.

        Returns:
            'forward' if n → n+1 (synthesis)
            'reverse' if n+1 → n (degradation)
            None if no growth pattern
        """
        # Normalize indices
        left_norm = left_idx.replace(" ", "")
        right_norm = right_idx.replace(" ", "")

        # Forward: n -> n+1
        if left_norm == "n" and right_norm == "n+1":
            return "forward"
        # Reverse: n+1 -> n
        if left_norm == "n+1" and right_norm == "n":
            return "reverse"
        # Also handle variations like n-1 → n
        if left_norm == "n-1" and right_norm == "n":
            return "forward"
        if left_norm == "n" and right_norm == "n-1":
            return "reverse"

        return None

    def _has_nucleotide(self, reaction: Reaction) -> bool:
        """Check if reaction has RNA nucleotide substrate (using ChEBI IDs)."""
        return any(
            p.chebi_id in RNA_NUCLEOTIDE_CHEBI
            for p in reaction.left_participants + reaction.right_participants
        )

    def _has_phosphate_product(self, reaction: Reaction) -> bool:
        """Check if reaction has phosphate or pyrophosphate product (using ChEBI IDs)."""
        return any(
            p.chebi_id in {CHEBI_PHOSPHATE, CHEBI_DIPHOSPHATE}
            for p in reaction.right_participants + reaction.left_participants
        )
