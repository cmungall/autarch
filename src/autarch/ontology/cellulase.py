"""Cellulase reaction classification using pattern DSL.

Cellulases are hydrolases that cleave cellulose and related β-1,4-glucan chains.
EC 3.2.1.4 classification.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.hydrolase_acting_on_glycosyl_bonds import HydrolaseActingOnGlycosylBonds
from autarch.ontology.reaction_diff import ReactionDiff
from autarch.pattern_dsl import var, optional, match_patterns
from autarch.molecules import CHEBI_H2O, h_plus, p, water


class Cellulase(HydrolaseActingOnGlycosylBonds):
    """cellulase

    Examples:
    - Endocellulase: cellulose + H2O → cellobiose + shorter cellulose chains
    - Exocellulase: cellulose + H2O → cellobiose (from chain ends)
    - β-glucosidase: cellobiose + H2O → 2 glucose
    - Mixed cellulase: cellulose + H2O → glucose + oligosaccharides
    """

    GO_ID = "GO:0008810"  # cellulase activity
    EC_NUMBER_PREFIX = "3.2.1.4"  # Cellulase

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # Cellulose hydrolysis: cellulose + H2O → oligomers + glucose
        var("cellulose") + p(water) >> var("oligomers") + var("glucose") + optional(h_plus),
        # Cellobiose hydrolysis: cellobiose + H2O → 2 glucose
        var("cellobiose") + p(water) >> var("glucose1") + var("glucose2") + optional(h_plus),
        # General β-glucan hydrolysis: β-glucan + H2O → smaller glucans
        var("beta_glucan") + p(water) >> var("smaller_glucan") + var("glucose_unit") + optional(h_plus),
        # Polysaccharide cleavage: polysaccharide + H2O → oligosaccharides
        var("polysaccharide") + p(water) >> var("oligosaccharides") + optional(h_plus),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a cellulase using pattern matching.

        Strategy:
        1. Must be a glycosidase (parent class)
        2. Must involve cellulose or related β-1,4-glucans
        3. Must use water for hydrolysis
        4. Must produce glucose or glucose oligomers
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
                explanation="No water reactant - not cellulase hydrolysis"
            )

        # Look for cellulose and related substrates

        cellulose_chebi = {
            "CHEBI:16991",  # cellobiose
            "CHEBI:17062",  # cellulose
        }

        has_cellulose_substrate = any(
            p.chebi_id in cellulose_chebi
            for p in reaction.left_participants
        )

        if not has_cellulose_substrate:
            return ClassificationResult(
                is_member=False,
                explanation="No cellulose or β-glucan substrate detected"
            )

        # Look for glucose products (characteristic of cellulase activity)

        glucose_chebi = {
            "CHEBI:4167",   # glucose
            "CHEBI:17634",  # glucose
            "CHEBI:16991",  # cellobiose
        }

        any(
            p.chebi_id in glucose_chebi
            for p in reaction.right_participants
        )

        # Try to match patterns
        match = match_patterns(reaction, self.PATTERNS, strict=False)

        explanation = "Cellulase: β-1,4-glucan hydrolysis"

        if match.matched:
            explanation += " [pattern matched]"

        return ClassificationResult(is_member=True, explanation=explanation)
