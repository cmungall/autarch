"""Oxidoreductase acting on iron-sulfur proteins reaction classification.

EC 1.18: Oxidoreductases acting on iron-sulfur proteins as donors.
These enzymes use iron-sulfur proteins (ferredoxin, rubredoxin) as
electron donors. Examples include ferredoxin-NADP+ reductase,
nitrogenase, and rubredoxin-NAD+ reductase.

Key features:
- Ferredoxin as electron donor/acceptor
- Rubredoxin as electron donor/acceptor
- Iron-sulfur cluster proteins in electron transfer chains
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.oxidoreductase import Oxidoreductase

# Ferredoxin ChEBI IDs
CHEBI_FERREDOXIN_OXIDIZED = "CHEBI:33737"  # oxidized [2Fe-2S] ferredoxin
CHEBI_FERREDOXIN_REDUCED = "CHEBI:33738"  # reduced [2Fe-2S] ferredoxin


class OxidoreductaseActingOnIronSulfurProteins(Oxidoreductase):
    """Catalysis of an oxidation-reduction (redox) reaction in which an iron-sulfur protein acts as hydrogen or electron donor and reduces an acceptor."""

    GO_ID = "GO:0016730"  # oxidoreductase activity, acting on iron-sulfur proteins as donors
    EC_NUMBER_PREFIX = "1.18.-.-"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is an oxidoreductase acting on iron-sulfur proteins.

        Strategy:
        1. Must be an oxidoreductase (parent class check)
        2. Look for ferredoxin or rubredoxin as participants (by ChEBI or label)
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

        # Check for ferredoxin by ChEBI ID
        has_ferredoxin_chebi = any(
            p.chebi_id in {CHEBI_FERREDOXIN_OXIDIZED, CHEBI_FERREDOXIN_REDUCED}
            for p in all_participants
        )

        # Check for ferredoxin or rubredoxin in label
        has_ferredoxin_label = "ferredoxin" in label_lower
        has_rubredoxin_label = "rubredoxin" in label_lower

        if has_ferredoxin_chebi or has_ferredoxin_label:
            return ClassificationResult(
                is_member=True,
                explanation="Iron-sulfur protein oxidoreductase: ferredoxin-dependent",
            )

        if has_rubredoxin_label:
            return ClassificationResult(
                is_member=True,
                explanation="Iron-sulfur protein oxidoreductase: rubredoxin-dependent",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No iron-sulfur protein (ferredoxin/rubredoxin) involvement detected",
        )
