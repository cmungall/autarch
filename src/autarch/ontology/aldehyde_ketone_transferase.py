"""Aldehyde/ketone transferase reaction classification.

Transferases that transfer aldehyde or ketonic groups (EC 2.2).
This includes transketolases and transaldolases that transfer
2-carbon or 3-carbon ketol/aldol units between sugar phosphates.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.transferase import Transferase


class AldehydeKetoneTransferase(Transferase):
    """Catalysis of the transfer of an aldehyde or ketonic group from one compound (donor) to another (acceptor)."""

    GO_ID = "GO:0016744"  # transketolase or transaldolase activity
    EC_NUMBER_PREFIX = "2.2.-.-"  # Transferring aldehyde or ketonic groups

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is an aldehyde/ketone transferase.

        Strategy:
        1. Must be a transferase (parent class)
        2. Look for transketolase or transaldolase indicators in labels
        3. Look for dihydroxyacetone involvement (common intermediate)
        """
        # First check if it's a transferase
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not a transferase: {parent_result.explanation}",
            )

        label_lower = reaction.label.lower() if reaction.label else ""

        # Look for transketolase/transaldolase indicators
        ketol_indicators = [
            "transketolase",
            "transaldolase",
            "dihydroxyacetone",
            "glycolaldehyde",
            "sedoheptulose",
        ]

        for indicator in ketol_indicators:
            if indicator in label_lower:
                return ClassificationResult(
                    is_member=True,
                    explanation=f"Aldehyde/ketone transferase: label contains '{indicator}'",
                )

        return ClassificationResult(
            is_member=False,
            explanation="No aldehyde/ketone transfer indicators detected",
        )
