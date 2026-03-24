"""Hydrolase acting on carbon-sulfur bonds (GO:0046508, EC 3.13.-.-).

Catalysis of the hydrolysis of any carbon-sulfur bond, C-S.
These reactions cleave C-S bonds using water.
"""

import re

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.hydrolase import Hydrolase
from autarch.ontology.reaction_diff import ReactionDiff


# Name patterns indicating C-S bond hydrolysis
CS_NAME_PATTERNS = [
    "thioether",
    "thioester",
    "carbon-sulfur",
    "carbon-sulphur",
    "c-s bond",
    "sulfoquinovos",
    "sulfonat.*hydro",
    "desulfinat",
    "desulfonat",
    r"s-alkyl",
    "thioglycosid",
]


class HydrolaseActingOnCarbonSulfurBonds(Hydrolase):
    """Catalysis of the hydrolysis of any carbon-sulfur bond, C-S.

    These enzymes cleave C-S bonds using water, releasing sulfur-containing
    products (sulfite, sulfate, thiol, etc.).

    Examples include:
    - UDP-sulfoquinovose synthase
    - S-formylglutathione hydrolase (EC 3.1.2.12, though classified under ester hydrolases)
    - Carbon-sulfur lyases that use water
    """

    GO_ID = "GO:0046508"  # hydrolase activity, acting on carbon-sulfur bonds
    EC_NUMBER_PREFIX = "3.13.-.-"  # EC prefix for C-S bond hydrolases

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is hydrolysis of a C-S bond.

        Strategy:
        1. Must be a hydrolase (water-mediated bond cleavage)
        2. Check reaction/participant names for C-S bond patterns
        3. Check for sulfur in reactants and sulfur species in products
        """
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not a hydrolase reaction: {parent_result.explanation}",
            )

        # Check names for C-S bond patterns
        has_cs_name = False
        matched_name = ""
        all_names = []
        for p in reaction.left_participants + reaction.right_participants:
            if p.name:
                all_names.append(p.name.lower())

        for name in all_names:
            for pattern in CS_NAME_PATTERNS:
                if re.search(pattern, name):
                    has_cs_name = True
                    matched_name = name
                    break
            if has_cs_name:
                break

        # Strong name-based evidence
        if has_cs_name:
            return ClassificationResult(
                is_member=True,
                explanation=f"C-S bond hydrolysis: pattern in '{matched_name}'",
            )

        # Structural: sulfur in reactants, carbon in reactants, sulfur species in products
        diff = ReactionDiff(reaction)
        has_sulfur = "S" in diff.reactant_elements
        has_carbon = "C" in diff.reactant_elements

        if has_sulfur and has_carbon:
            # Look for sulfur-containing products (sulfite, sulfate, thiol, H2S, etc.)
            has_sulfur_product = any(
                p.name and (
                    "sulfit" in p.name.lower()
                    or "sulfid" in p.name.lower()
                    or "thiol" in p.name.lower()
                    or "mercapto" in p.name.lower()
                    or "hydrogen sulfide" in p.name.lower()
                )
                for p in reaction.right_participants
            )

            # Also check for sulfur-containing reactant that loses sulfur
            has_sulfur_substrate = any(
                p.name and (
                    "thio" in p.name.lower()
                    or "sulfo" in p.name.lower()
                    or "mercapto" in p.name.lower()
                )
                for p in reaction.left_participants
            )

            if has_sulfur_substrate and has_sulfur_product:
                return ClassificationResult(
                    is_member=True,
                    explanation="C-S bond hydrolysis: sulfur-containing substrate yields sulfur product",
                )

        return ClassificationResult(
            is_member=False,
            explanation="No evidence of C-S bond hydrolysis",
        )
