"""Hydrolase acting on carbon-phosphorus bonds (GO:0016827, EC 3.11.-.-).

Catalysis of the hydrolysis of any carbon-phosphorus bond.
These reactions cleave C-P bonds (phosphonate bonds) using water.
"""

import re

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.hydrolase import Hydrolase
from autarch.ontology.reaction_diff import ReactionDiff


# Name patterns indicating C-P bond hydrolysis
CP_NAME_PATTERNS = [
    "phosphonat",
    "phosphonoacetaldehyde",
    "phosphonoacetate",
    "phosphonopyruvat",
    "carbon-phosphorus",
    "c-p bond",
    "c-p lyase",  # sometimes misnamed but related
    "phosphonoformate",
    "phosphonomethyl",
    "aminoethylphosphonate",
]


class HydrolaseActingOnCarbonPhosphorusBonds(Hydrolase):
    """Catalysis of the hydrolysis of any carbon-phosphorus bond.

    These enzymes cleave C-P bonds (phosphonate bonds) using water.
    Phosphonates contain a direct C-P bond, unlike phosphate esters
    which have C-O-P linkages.

    Examples include:
    - Phosphonoacetaldehyde hydrolase (EC 3.11.1.1)
    - Phosphonoacetate hydrolase (EC 3.11.1.2)
    - Phosphonopyruvate hydrolase (EC 3.11.1.3)
    """

    GO_ID = "GO:0016827"  # hydrolase activity, acting on carbon-phosphorus bonds
    EC_NUMBER_PREFIX = "3.11.-.-"  # EC prefix for C-P bond hydrolases

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is hydrolysis of a C-P bond.

        Strategy:
        1. Must be a hydrolase (water-mediated bond cleavage)
        2. Check reaction/participant names for phosphonate patterns
        3. Check for phosphorus in reactants (C-P bond substrates)
        4. Look for phosphate product (P released from C-P cleavage)
        """
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not a hydrolase reaction: {parent_result.explanation}",
            )

        # Check names for C-P bond patterns
        has_cp_name = False
        matched_name = ""
        all_names = []
        for p in reaction.left_participants + reaction.right_participants:
            if p.name:
                all_names.append(p.name.lower())

        for name in all_names:
            for pattern in CP_NAME_PATTERNS:
                if re.search(pattern, name):
                    has_cp_name = True
                    matched_name = name
                    break
            if has_cp_name:
                break

        # Strong name-based evidence
        if has_cp_name:
            return ClassificationResult(
                is_member=True,
                explanation=f"C-P bond hydrolysis: phosphonate pattern in '{matched_name}'",
            )

        # Structural check: phosphorus in reactants with phosphate product
        diff = ReactionDiff(reaction)
        has_phosphorus = "P" in diff.reactant_elements
        has_carbon = "C" in diff.reactant_elements

        if has_phosphorus and has_carbon:
            has_phosphate_product = any(
                p.chebi_id == "CHEBI:43474"
                or (p.name and "phosphate" in p.name.lower() and "phosphonat" not in p.name.lower())
                for p in reaction.right_participants
            )
            # Look for phosphonate in reactants specifically
            has_phosphonate_reactant = any(
                p.name and "phosphonat" in p.name.lower()
                for p in reaction.left_participants
            )

            if has_phosphonate_reactant and has_phosphate_product:
                return ClassificationResult(
                    is_member=True,
                    explanation="C-P bond hydrolysis: phosphonate substrate yields phosphate",
                )

        return ClassificationResult(
            is_member=False,
            explanation="No evidence of C-P bond hydrolysis",
        )
