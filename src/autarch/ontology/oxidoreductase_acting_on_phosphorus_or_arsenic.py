"""Oxidoreductase acting on phosphorus or arsenic reaction classification.

EC 1.20: Oxidoreductases acting on phosphorus or arsenic in donors.
These enzymes catalyze redox reactions where the donor contains phosphorus
or arsenic. Examples include arsenate reductase (glutaredoxin),
phosphonate dehydrogenase, and arsenite oxidase.

Key features:
- Arsenate/arsenite interconversion (arsenic redox)
- Phosphonate oxidation (C-P bond cleavage with electron transfer)
- Hypophosphite oxidation
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.oxidoreductase import Oxidoreductase

# Arsenic species
CHEBI_ARSENATE = "CHEBI:22633"  # arsenate
CHEBI_ARSENITE = "CHEBI:29242"  # arsenite


class OxidoreductaseActingOnPhosphorusOrArsenic(Oxidoreductase):
    """Catalysis of an oxidation-reduction (redox) reaction in which phosphorus or arsenic in the donor group acts as hydrogen or electron donor and reduces an acceptor."""

    GO_ID = "GO:0030613"  # oxidoreductase activity, acting on phosphorus or arsenic in donors
    EC_NUMBER_PREFIX = "1.20.-.-"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is an oxidoreductase acting on phosphorus or arsenic.

        Strategy:
        1. Must be an oxidoreductase (parent class check)
        2. Look for arsenate/arsenite or phosphonate as participants
        3. Label-based detection for key substrate names
        """
        # First check parent
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not an oxidoreductase: {parent_result.explanation}",
            )

        all_participants = reaction.left_participants + reaction.right_participants
        label_lower = reaction.label.lower() if reaction.label else ""

        # Check for arsenate/arsenite by ChEBI ID
        has_arsenic_chebi = any(
            p.chebi_id in {CHEBI_ARSENATE, CHEBI_ARSENITE}
            for p in all_participants
        )

        # Label-based detection for arsenic and phosphorus donors
        arsenic_phosphorus_labels = [
            "arsenate",
            "arsenite",
            "phosphonate",
            "hypophosphite",
            "phosphite",
        ]
        has_label_match = any(pat in label_lower for pat in arsenic_phosphorus_labels)

        if has_arsenic_chebi:
            return ClassificationResult(
                is_member=True,
                explanation="Phosphorus/arsenic oxidoreductase: arsenate/arsenite redox",
            )

        if has_label_match:
            return ClassificationResult(
                is_member=True,
                explanation=f"Phosphorus/arsenic oxidoreductase: label match in '{label_lower}'",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No phosphorus or arsenic donor involvement detected",
        )
