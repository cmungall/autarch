"""Carbon-phosphorus lyase reaction classification.

EC 4.7: Carbon-phosphorus lyases catalyze the cleavage of C-P bonds.
These enzymes degrade phosphonate compounds by breaking the direct
carbon-phosphorus bond.

Pattern: R-PO3 → R-H + phosphorus species
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.lyase import Lyase


class CarbonPhosphorusLyase(Lyase):
    """carbon-phosphorus lyase

    EC 4.7.x.x includes:
    - Phosphonoacetaldehyde hydrolase (4.7.1.1) - not a true hydrolase
    - Phosphonoacetate hydrolase (4.7.1.2)
    - Methylphosphonate synthase

    Examples:
    - Phosphonoacetaldehyde = acetaldehyde + phosphate
    - Methylphosphonate degradation
    """

    GO_ID = "GO:0018835"  # carbon-phosphorus lyase activity
    EC_NUMBER_PREFIX = "4.7.-.-"

    CP_LABEL_PATTERNS = [
        "phosphonate", "phosphonic",
        "methylphosphonate", "phosphonoacet",
        "c-p lyase", "carbon-phosphorus",
        "phosphonatase",
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a carbon-phosphorus lyase.

        Strategy:
        1. Must pass parent Lyase checks
        2. Look for phosphonate substrates or C-P lyase label patterns
        """
        parent_result = super().check_membership_impl(reaction)

        label_lower = reaction.label.lower() if reaction.label else ""

        has_cp_label = any(
            pattern in label_lower
            for pattern in self.CP_LABEL_PATTERNS
        )

        if has_cp_label:
            if parent_result.is_member:
                return ClassificationResult(
                    is_member=True,
                    explanation="Carbon-phosphorus lyase: C-P bond cleavage"
                )
            # Phosphonate patterns are strong enough signals even without parent match
            return ClassificationResult(
                is_member=True,
                explanation="Carbon-phosphorus lyase: phosphonate substrate"
            )

        if parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation="Lyase but no C-P bond cleavage pattern"
            )

        return ClassificationResult(
            is_member=False,
            explanation=f"Not a lyase: {parent_result.explanation}"
        )
