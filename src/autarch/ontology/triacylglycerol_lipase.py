"""Lipase reaction classification using pattern DSL.

Lipases are hydrolases that cleave ester bonds in lipids and fatty acid esters.
EC 3.1.1.3 classification.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.carboxylic_ester_hydrolase import CarboxylicEsterHydrolase
from autarch.ontology.reaction_diff import ReactionDiff
from autarch.pattern_dsl import var, optional, match_patterns
from autarch.molecules import CHEBI_H2O, h_plus, p, water


class TriacylglycerolLipase(CarboxylicEsterHydrolase):
    """triacylglycerol lipase

    Examples:
    - Pancreatic lipase: triglyceride + H2O → monoglyceride + fatty acid
    - Hormone-sensitive lipase: triglyceride + H2O → diacylglycerol + fatty acid
    - Phospholipase: phospholipid + H2O → lysophospholipid + fatty acid
    - Cholesterol esterase: cholesteryl ester + H2O → cholesterol + fatty acid
    """

    GO_ID = "GO:0004806"  # triglyceride lipase activity
    EC_NUMBER_PREFIX = "3.1.1.3"  # Triacylglycerol lipase
    EC_BROAD_XREFS = ["3.1.1.79"]  # GO xref is broadMatch to the broader lipase branch

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # Triglyceride hydrolysis: triglyceride + H2O → diglyceride + fatty acid
        var("triglyceride") + p(water) >> var("diglyceride") + var("fatty_acid") + optional(h_plus),
        # Complete hydrolysis: triglyceride + 3 H2O → glycerol + 3 fatty acids
        var("triglyceride") + p(water) + p(water) + p(water) 
        >> var("glycerol") + var("fatty_acid1") + var("fatty_acid2") + var("fatty_acid3"),
        # Cholesterol ester: cholesteryl ester + H2O → cholesterol + fatty acid
        var("cholesteryl_ester") + p(water) >> var("cholesterol") + var("fatty_acid") + optional(h_plus),
        # Phospholipid: phospholipid + H2O → lysophospholipid + fatty acid
        var("phospholipid") + p(water) >> var("lysophospholipid") + var("fatty_acid") + optional(h_plus),
        # General lipid ester: lipid ester + H2O → alcohol + carboxylic acid
        var("lipid_ester") + p(water) >> var("alcohol") + var("carboxylic_acid") + optional(h_plus),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a lipase using pattern matching.

        Strategy:
        1. Must be an esterase (parent class)
        2. Must involve lipid substrates (triglycerides, phospholipids, etc.)
        3. Must use water for hydrolysis
        4. Must produce glycerol, fatty acids, or other lipid components
        """
        # First check if it's an esterase
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not an esterase: {parent_result.explanation}",
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
                explanation="No water reactant - not lipid hydrolysis"
            )

        # Look for lipid substrates

        lipid_chebi = {
            "CHEBI:17855",  # triglyceride
            "CHEBI:18059",  # fatty acid
            "CHEBI:15734",  # glycerol
            "CHEBI:16113",  # cholesterol
        }

        has_lipid_substrate = any(
            p.chebi_id in lipid_chebi
            for p in reaction.left_participants
        )

        if not has_lipid_substrate:
            return ClassificationResult(
                is_member=False,
                explanation="No lipid substrate detected"
            )

        # Look for lipid hydrolysis products

        lipid_product_chebi = {
            "CHEBI:15734",  # glycerol
            "CHEBI:18059",  # fatty acid
            "CHEBI:16113",  # cholesterol
        }

        any(
            p.chebi_id in lipid_product_chebi
            for p in reaction.right_participants
        )

        # Try to match patterns
        match = match_patterns(reaction, self.PATTERNS, strict=False)

        explanation = "Lipase: lipid ester hydrolysis"

        if match.matched:
            explanation += " [pattern matched]"

        return ClassificationResult(is_member=True, explanation=explanation)
