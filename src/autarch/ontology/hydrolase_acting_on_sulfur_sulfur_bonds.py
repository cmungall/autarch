"""Hydrolase acting on sulfur-sulfur bonds (GO:0016828, EC 3.12.-.-).

Catalysis of the hydrolysis of any sulfur-sulfur bond.
These reactions cleave S-S bonds using water.
"""

import re

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.hydrolase import Hydrolase
from autarch.ontology.reaction_diff import ReactionDiff


# Name patterns indicating S-S bond hydrolysis
SS_NAME_PATTERNS = [
    "trithionat",
    "disulfide",
    "sulfur-sulfur",
    "s-s bond",
    "thiosulfat.*hydro",
    "polysulfide",
]


class HydrolaseActingOnSulfurSulfurBonds(Hydrolase):
    """Catalysis of the hydrolysis of any sulfur-sulfur bond.

    These enzymes cleave S-S bonds using water, releasing sulfur-containing
    species as products.

    Examples include:
    - Trithionate hydrolase (EC 3.12.1.1)
    """

    GO_ID = "GO:0016828"  # hydrolase activity, acting on sulfur-sulfur bonds
    EC_NUMBER_PREFIX = "3.12.-.-"  # EC prefix for S-S bond hydrolases

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is hydrolysis of an S-S bond.

        Strategy:
        1. Must be a hydrolase (water-mediated bond cleavage)
        2. Check reaction/participant names for S-S bond patterns
        3. Check for sulfur in both reactants and products
        4. Look for multiple sulfur species in products (from S-S cleavage)
        """
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not a hydrolase reaction: {parent_result.explanation}",
            )

        # Check names for S-S bond patterns
        has_ss_name = False
        matched_name = ""
        all_names = []
        for p in reaction.left_participants + reaction.right_participants:
            if p.name:
                all_names.append(p.name.lower())

        for name in all_names:
            for pattern in SS_NAME_PATTERNS:
                if re.search(pattern, name):
                    has_ss_name = True
                    matched_name = name
                    break
            if has_ss_name:
                break

        # Strong name-based evidence
        if has_ss_name:
            return ClassificationResult(
                is_member=True,
                explanation=f"S-S bond hydrolysis: pattern in '{matched_name}'",
            )

        # Structural: sulfur in reactants, multiple sulfur species in products
        diff = ReactionDiff(reaction)
        has_sulfur_reactant = "S" in diff.reactant_elements
        has_sulfur_product = "S" in diff.product_elements

        if has_sulfur_reactant and has_sulfur_product:
            # Count sulfur-containing products (excluding water/H+)
            sulfur_products = [
                p for p in reaction.right_participants
                if p.name and "sulf" in p.name.lower()
            ]
            if len(sulfur_products) >= 2:
                return ClassificationResult(
                    is_member=True,
                    explanation="S-S bond hydrolysis: multiple sulfur species in products",
                )

        return ClassificationResult(
            is_member=False,
            explanation="No evidence of S-S bond hydrolysis",
        )
