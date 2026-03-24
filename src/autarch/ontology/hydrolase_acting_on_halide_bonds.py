"""Hydrolase acting on halide bonds (GO:0016824, EC 3.8.-.-).

Catalysis of the hydrolysis of any halide bond, i.e. C-halogen bonds.
These are dehalogenase reactions where a halogen substituent (F, Cl, Br, I)
is removed from a substrate via hydrolysis, releasing a halide ion.
"""

import re

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.hydrolase import Hydrolase
from autarch.ontology.reaction_diff import ReactionDiff


# ChEBI IDs for halide ions
HALIDE_CHEBI_IDS = {
    "CHEBI:17996",  # chloride
    "CHEBI:24061",  # fluoride
    "CHEBI:29227",  # bromide
    "CHEBI:16788",  # iodide
}

# Halogen element symbols
HALOGEN_ELEMENTS = {"F", "Cl", "Br", "I"}

# Name patterns indicating halide bond hydrolysis
HALIDE_NAME_PATTERNS = [
    "dehalogenase",
    "chloro",
    "bromo",
    "fluoro",
    "iodo",
    "halide",
    "halogenase",  # some reverse reactions named this way
    "dechlorin",
    "debrom",
    "defluorin",
    "deiodin",
]


class HydrolaseActingOnHalideBonds(Hydrolase):
    """Catalysis of the hydrolysis of any halide bond.

    Dehalogenases cleave C-X bonds (where X is F, Cl, Br, or I)
    using water, releasing a halide ion as product.

    Examples include:
    - Haloalkane dehalogenases
    - Haloacetate dehalogenases
    - Atrazine chlorohydrolase
    """

    GO_ID = "GO:0016824"  # hydrolase activity, acting on halide bonds
    EC_NUMBER_PREFIX = "3.8.-.-"  # EC prefix for halide bond hydrolases

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is hydrolysis of a halide bond.

        Strategy:
        1. Must be a hydrolase (water-mediated bond cleavage)
        2. Check for halide ions in products (Cl-, Br-, F-, I-)
        3. Check for halogens in reactant elements
        4. Check reaction/participant names for dehalogenase patterns
        """
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not a hydrolase reaction: {parent_result.explanation}",
            )

        # Check for halide ion products by ChEBI ID
        has_halide_product = any(
            p.chebi_id in HALIDE_CHEBI_IDS
            for p in reaction.right_participants
        )

        # Check for halide ion products by name
        if not has_halide_product:
            has_halide_product = any(
                p.name and p.name.lower() in ("chloride", "bromide", "fluoride", "iodide")
                for p in reaction.right_participants
            )

        # Check for halogen elements in reactants
        diff = ReactionDiff(reaction)
        has_halogen_in_reactants = bool(HALOGEN_ELEMENTS & diff.reactant_elements)

        # Check names for dehalogenase patterns
        has_halide_name = False
        all_names = []
        for p in reaction.left_participants + reaction.right_participants:
            if p.name:
                all_names.append(p.name.lower())

        for name in all_names:
            for pattern in HALIDE_NAME_PATTERNS:
                if re.search(pattern, name):
                    has_halide_name = True
                    break
            if has_halide_name:
                break

        # Positive if we have structural evidence (halide product + halogen reactant)
        # or strong name-based evidence
        if has_halide_product and has_halogen_in_reactants:
            return ClassificationResult(
                is_member=True,
                explanation="Halide bond hydrolysis: halogen in reactant, halide ion in product",
            )

        if has_halide_name and (has_halide_product or has_halogen_in_reactants):
            return ClassificationResult(
                is_member=True,
                explanation="Halide bond hydrolysis: dehalogenase-related name pattern",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No evidence of halide bond hydrolysis",
        )
