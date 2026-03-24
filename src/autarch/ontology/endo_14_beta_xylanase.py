"""Endo-1,4-beta-xylanase reaction classification using pattern DSL.

Xylanases are hydrolases that cleave xylan (β-1,4-linked xylose chains).
EC 3.2.1.8 classification.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.hydrolase_acting_on_glycosyl_bonds import HydrolaseActingOnGlycosylBonds
from autarch.ontology.reaction_diff import ReactionDiff
from autarch.pattern_dsl import var, optional, match_patterns
from autarch.molecules import CHEBI_H2O, h_plus, p, water


class Endo14BetaXylanase(HydrolaseActingOnGlycosylBonds):
    """endo-1,4-beta-xylanase

    Examples:
    - Endoxylanase: xylan + H2O → xylobiose + shorter xylan chains
    - β-xylosidase: xylobiose + H2O → 2 xylose
    - Arabinoxylanase: arabinoxylan + H2O → xylose + arabinose oligomers
    - Glucuronoxylanase: glucuronoxylan + H2O → xylose + glucuronic acid
    """

    GO_ID = "GO:0031176"  # endo-1,4-beta-xylanase activity
    EC_NUMBER_PREFIX = "3.2.1.8"  # Endo-1,4-beta-xylanase

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # Xylan hydrolysis: xylan + H2O → oligomers + xylose
        var("xylan") + p(water) >> var("oligomers") + var("xylose") + optional(h_plus),
        # Xylobiose hydrolysis: xylobiose + H2O → 2 xylose
        var("xylobiose") + p(water) >> var("xylose1") + var("xylose2") + optional(h_plus),
        # Arabinoxylan hydrolysis: arabinoxylan + H2O → xylose + arabinose
        var("arabinoxylan") + p(water) >> var("xylose") + var("arabinose") + optional(h_plus),
        # General xylan cleavage: xylan polymer + H2O → smaller fragments
        var("xylan_polymer") + p(water) >> var("smaller_fragments") + optional(h_plus),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a xylanase using pattern matching.

        Strategy:
        1. Must be a glycosidase (parent class)
        2. Must involve xylan or xylose polymers
        3. Must use water for hydrolysis
        4. Must produce xylose or xylose oligomers
        """
        # First check if it's a glycosidase
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not a glycosidase: {parent_result.explanation}",
            )

        ReactionDiff(reaction)

        # Must use water as reactant
        has_water = any(
            p.chebi_id == CHEBI_H2O
            for p in reaction.left_participants
        )

        if not has_water:
            return ClassificationResult(
                is_member=False,
                explanation="No water reactant - not endo-1,4-beta-xylanase hydrolysis"
            )

        # Look for xylan and related substrates

        xylan_chebi = {
            "CHEBI:27560",  # xylose
            "CHEBI:18222",  # xylan
        }

        has_xylan_substrate = any(
            p.chebi_id in xylan_chebi
            for p in reaction.left_participants
        )

        if not has_xylan_substrate:
            return ClassificationResult(
                is_member=False,
                explanation="No xylan or xylose substrate detected"
            )

        # Look for xylose products

        xylose_chebi = {
            "CHEBI:27560",  # xylose
            "CHEBI:15940",  # L-arabinose
        }

        any(
            p.chebi_id in xylose_chebi
            for p in reaction.right_participants
        )

        # Try to match patterns
        match = match_patterns(reaction, self.PATTERNS, strict=False)

        explanation = "Endo-1,4-beta-xylanase: beta-1,4-xylan hydrolysis"

        if match.matched:
            explanation += " [pattern matched]"

        return ClassificationResult(is_member=True, explanation=explanation)
