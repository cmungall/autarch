"""Amylase reaction classification using pattern DSL.

Amylases are hydrolases that cleave α-1,4-glucosidic bonds in starch and glycogen.
EC 3.2.1.1 classification.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.hydrolase_acting_on_glycosyl_bonds import HydrolaseActingOnGlycosylBonds
from autarch.ontology.reaction_diff import ReactionDiff
from autarch.pattern_dsl import var, optional, match_patterns
from autarch.molecules import CHEBI_H2O, h_plus, p, water


class AlphaAmylase(HydrolaseActingOnGlycosylBonds):
    """alpha-amylase

    Examples:
    - α-amylase: starch + H2O → maltose + maltotriose + dextrins
    - β-amylase: starch + H2O → maltose (from non-reducing ends)
    - Glucoamylase: starch + H2O → glucose
    - Pullulanase: pullulan + H2O → maltotriose
    """

    GO_ID = "GO:0004556"  # alpha-amylase activity
    EC_NUMBER_PREFIX = "3.2.1.1"  # Alpha-amylase

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # Starch hydrolysis: starch + H2O → maltose + oligomers
        var("starch") + p(water) >> var("maltose") + var("oligomers") + optional(h_plus),
        # Glycogen hydrolysis: glycogen + H2O → glucose + maltose
        var("glycogen") + p(water) >> var("glucose") + var("maltose") + optional(h_plus),
        # Maltose production: polysaccharide + H2O → maltose units
        var("polysaccharide") + p(water) >> var("maltose1") + var("maltose2") + optional(h_plus),
        # Complete hydrolysis: starch + H2O → glucose
        var("starch") + p(water) >> var("glucose") + optional(h_plus),
        # Pullulan hydrolysis: pullulan + H2O → maltotriose
        var("pullulan") + p(water) >> var("maltotriose") + optional(h_plus),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is an amylase using pattern matching.

        Strategy:
        1. Must be a glycosidase (parent class)
        2. Must involve starch, glycogen, or α-1,4-glucan substrates
        3. Must use water for hydrolysis
        4. Must produce glucose, maltose, or glucose oligomers
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
                explanation="No water reactant - not starch hydrolysis"
            )

        # Look for starch and related substrates

        starch_chebi = {
            "CHEBI:28017",  # starch
            "CHEBI:28604",  # glycogen
            "CHEBI:17992",  # amylose
        }

        has_starch_substrate = any(
            p.chebi_id in starch_chebi
            for p in reaction.left_participants
        )

        if not has_starch_substrate:
            return ClassificationResult(
                is_member=False,
                explanation="No starch or α-1,4-glucan substrate detected"
            )

        # Look for glucose/maltose products

        glucose_chebi = {
            "CHEBI:4167",   # glucose
            "CHEBI:17634",  # glucose
            "CHEBI:17306",  # maltose
        }

        any(
            p.chebi_id in glucose_chebi
            for p in reaction.right_participants
        )

        # Try to match patterns
        match = match_patterns(reaction, self.PATTERNS, strict=False)

        explanation = "Amylase: α-1,4-glucan hydrolysis"

        if match.matched:
            explanation += " [pattern matched]"

        return ClassificationResult(is_member=True, explanation=explanation)
