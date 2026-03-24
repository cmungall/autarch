"""Oxidoreductase acting on reduced flavodoxin reaction classification.

EC 1.19: Oxidoreductases acting on reduced flavodoxin as donor.
These enzymes use reduced flavodoxin as the electron donor.
Flavodoxin is a small electron transfer protein containing FMN as
a prosthetic group. Examples include flavodoxin-based nitrogenase
reactions.

Key features:
- Flavodoxin as electron donor
- Often overlaps with nitrogenase reactions
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.oxidoreductase import Oxidoreductase


class OxidoreductaseActingOnReducedFlavodoxin(Oxidoreductase):
    """Catalysis of an oxidation-reduction (redox) reaction in which reduced flavodoxin acts as hydrogen or electron donor and reduces an acceptor."""

    GO_ID = "GO:0016737"  # oxidoreductase activity, acting on reduced flavodoxin as donor
    EC_NUMBER_PREFIX = "1.19.-.-"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is an oxidoreductase acting on reduced flavodoxin.

        Strategy:
        1. Must be an oxidoreductase (parent class check)
        2. Look for flavodoxin as participant (by label or participant name)
        """
        # First check parent
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not an oxidoreductase: {parent_result.explanation}",
            )

        label_lower = reaction.label.lower() if reaction.label else ""

        # Check for flavodoxin in label
        has_flavodoxin_label = "flavodoxin" in label_lower

        # Check participant names for flavodoxin
        has_flavodoxin_participant = any(
            p.name and "flavodoxin" in p.name.lower()
            for p in reaction.left_participants + reaction.right_participants
        )

        if has_flavodoxin_label or has_flavodoxin_participant:
            return ClassificationResult(
                is_member=True,
                explanation="Flavodoxin-dependent oxidoreductase",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No flavodoxin involvement detected",
        )
