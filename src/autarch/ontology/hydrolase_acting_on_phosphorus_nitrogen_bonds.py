"""Hydrolase acting on phosphorus-nitrogen bonds (GO:0016825, EC 3.9.-.-).

Catalysis of the hydrolysis of any phosphorus-nitrogen bond.
These reactions cleave P-N bonds using water, such as phosphoamidase activity.
"""

import re

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.hydrolase import Hydrolase
from autarch.ontology.reaction_diff import ReactionDiff


# Name patterns indicating P-N bond hydrolysis
PN_NAME_PATTERNS = [
    "phosphoamid",
    "phosphoramid",
    "phosphoramidat",
    "phosphonamidat",
    r"phospho.*amine",
    "phosphorus-nitrogen",
    "p-n bond",
]


class HydrolaseActingOnPhosphorusNitrogenBonds(Hydrolase):
    """Catalysis of the hydrolysis of any phosphorus-nitrogen bond.

    Phosphoamidases cleave P-N bonds using water, releasing
    phosphate/phosphonate and an amine.

    Examples include:
    - Phosphoamidase (EC 3.9.1.1)
    - Phosphonamidate hydrolase
    """

    GO_ID = "GO:0016825"  # hydrolase activity, acting on phosphorus-nitrogen bonds
    EC_NUMBER_PREFIX = "3.9.-.-"  # EC prefix for P-N bond hydrolases

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is hydrolysis of a P-N bond.

        Strategy:
        1. Must be a hydrolase (water-mediated bond cleavage)
        2. Check for both phosphorus and nitrogen in reactants
        3. Check reaction/participant names for phosphoamidase patterns
        4. Look for phosphate product + amine/ammonia product
        """
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not a hydrolase reaction: {parent_result.explanation}",
            )

        diff = ReactionDiff(reaction)
        has_phosphorus = "P" in diff.reactant_elements
        has_nitrogen = "N" in diff.reactant_elements

        # Check names for P-N bond patterns
        has_pn_name = False
        all_names = []
        for p in reaction.left_participants + reaction.right_participants:
            if p.name:
                all_names.append(p.name.lower())

        for name in all_names:
            for pattern in PN_NAME_PATTERNS:
                if re.search(pattern, name):
                    has_pn_name = True
                    break
            if has_pn_name:
                break

        # Strong name-based evidence
        if has_pn_name:
            return ClassificationResult(
                is_member=True,
                explanation="P-N bond hydrolysis: phosphoamidase-related name pattern",
            )

        # Structural: both P and N in reactant, with phosphate and amine products
        if has_phosphorus and has_nitrogen:
            has_phosphate_product = any(
                p.chebi_id == "CHEBI:43474"
                or (p.name and "phosphat" in p.name.lower())
                for p in reaction.right_participants
            )
            has_amine_product = any(
                p.chebi_id in ("CHEBI:16134", "CHEBI:28938")  # ammonia, ammonium
                or (p.name and ("amine" in p.name.lower() or "ammoni" in p.name.lower()))
                for p in reaction.right_participants
            )

            if has_phosphate_product and has_amine_product:
                return ClassificationResult(
                    is_member=True,
                    explanation="P-N bond hydrolysis: P and N in substrate, phosphate + amine products",
                )

        return ClassificationResult(
            is_member=False,
            explanation="No evidence of P-N bond hydrolysis",
        )
