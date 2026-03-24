"""Chitinase reaction classification using pattern DSL.

Chitinases are hydrolases that cleave chitin (β-1,4-linked N-acetylglucosamine chains).
EC 3.2.1.14 classification.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.hydrolase_acting_on_glycosyl_bonds import HydrolaseActingOnGlycosylBonds
from autarch.ontology.reaction_diff import ReactionDiff
from autarch.pattern_dsl import var, optional, match_patterns
from autarch.molecules import CHEBI_H2O, h_plus, p, water


class Chitinase(HydrolaseActingOnGlycosylBonds):
    """chitinase

    Examples:
    - Endochitinase: chitin + H2O → chitobiose + shorter chitin chains
    - Exochitinase: chitin + H2O → chitobiose (from chain ends)
    - Chitobiosidase: chitobiose + H2O → 2 N-acetylglucosamine
    - Bacterial chitinase: chitin + H2O → GlcNAc oligomers
    """

    GO_ID = "GO:0004568"  # chitinase activity
    EC_NUMBER_PREFIX = "3.2.1.14"  # Chitinase

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # Chitin hydrolysis: chitin + H2O → oligomers + GlcNAc
        var("chitin") + p(water) >> var("oligomers") + var("glcnac") + optional(h_plus),
        # Chitobiose hydrolysis: chitobiose + H2O → 2 GlcNAc
        var("chitobiose") + p(water) >> var("glcnac1") + var("glcnac2") + optional(h_plus),
        # General chitin cleavage: chitin polymer + H2O → smaller fragments
        var("chitin_polymer") + p(water) >> var("smaller_fragments") + optional(h_plus),
        # N-acetylglucosamine polymer: NAG polymer + H2O → NAG units
        var("nag_polymer") + p(water) >> var("nag_units") + optional(h_plus),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a chitinase using pattern matching.

        Strategy:
        1. Must be a glycosidase (parent class)
        2. Must involve chitin or N-acetylglucosamine polymers
        3. Must use water for hydrolysis
        4. Must produce N-acetylglucosamine or its oligomers
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
                explanation="No water reactant - not chitin hydrolysis"
            )

        # Look for chitin and related substrates

        chitin_chebi = {
            "CHEBI:16336",  # chitin
            "CHEBI:506227", # chitobiose
            "CHEBI:28007",  # N-acetyl-D-glucosamine
        }

        has_chitin_substrate = any(
            p.chebi_id in chitin_chebi
            for p in reaction.left_participants
        )

        if not has_chitin_substrate:
            return ClassificationResult(
                is_member=False,
                explanation="No chitin or N-acetylglucosamine substrate detected"
            )

        # Look for N-acetylglucosamine products

        glcnac_chebi = {
            "CHEBI:28007",  # N-acetyl-D-glucosamine
            "CHEBI:506227", # chitobiose
        }

        any(
            p.chebi_id in glcnac_chebi
            for p in reaction.right_participants
        )

        # Try to match patterns
        match = match_patterns(reaction, self.PATTERNS, strict=False)

        explanation = "Chitinase: β-1,4-N-acetylglucosamine hydrolysis"

        if match.matched:
            explanation += " [pattern matched]"

        return ClassificationResult(is_member=True, explanation=explanation)
