"""Peroxidase reaction classification using pattern DSL.

Peroxidases are oxidoreductases that use hydrogen peroxide as electron acceptor.
EC 1.11.1.x classification.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.oxidoreductase import Oxidoreductase
from autarch.ontology.reaction_diff import ReactionDiff
from autarch.pattern_dsl import var, match_patterns
from autarch.molecules import (
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_O2,
    h_plus,
    p,
    water,
)


class Peroxidase(Oxidoreductase):
    """peroxidase

    Examples:
    - Horseradish peroxidase: phenol + H2O2 → phenoxyl radical + 2 H2O
    - Myeloperoxidase: Cl- + H2O2 + H+ → HOCl + H2O
    - Catalase: 2 H2O2 → 2 H2O + O2
    - Glutathione peroxidase: GSH + H2O2 → GSSG + 2 H2O
    """

    GO_ID = "GO:0004601"  # peroxidase activity
    EC_NUMBER_PREFIX = "1.11.1.-"  # Peroxidases

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # General peroxidase: substrate + H2O2 → oxidized substrate + 2 H2O
        var("substrate") + var("h2o2") >> var("oxidized_substrate") + p(water) + p(water),
        # With proton: substrate + H2O2 + H+ → product + 2 H2O
        var("substrate") + var("h2o2") + p(h_plus) >> var("product") + p(water) + p(water),
        # Catalase: 2 H2O2 → 2 H2O + O2
        var("h2o2_1") + var("h2o2_2") >> p(water) + p(water) + var("oxygen"),
        # Haloperoxidase: halide + H2O2 + H+ → hypohalous acid + H2O
        var("halide") + var("h2o2") + p(h_plus) >> var("hypohalous_acid") + p(water),
        # Glutathione peroxidase: 2 GSH + H2O2 → GSSG + 2 H2O
        var("gsh_1") + var("gsh_2") + var("h2o2") >> var("gssg") + p(water) + p(water),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a peroxidase using pattern matching.

        Strategy:
        1. Must be an oxidoreductase (parent class)
        2. Must use hydrogen peroxide (H2O2) as oxidant
        3. Must produce water as product
        4. Look for characteristic peroxidase substrates
        """
        # First check if it's an oxidoreductase
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not an oxidoreductase: {parent_result.explanation}",
            )

        ReactionDiff(reaction)

        # Must use hydrogen peroxide as substrate
        h2o2_indicators = {
            CHEBI_H2O2,  # hydrogen peroxide
        }

        has_h2o2 = any(
            p.chebi_id in h2o2_indicators
            for p in reaction.left_participants
        )

        if not has_h2o2:
            return ClassificationResult(
                is_member=False,
                explanation="No hydrogen peroxide - not peroxidase"
            )

        # Must produce water (most peroxidases reduce H2O2 to H2O)
        has_water_product = any(
            p.chebi_id == CHEBI_H2O
            for p in reaction.right_participants
        )

        if not has_water_product:
            return ClassificationResult(
                is_member=False,
                explanation="No water product - H2O2 not reduced"
            )

        # Look for peroxidase substrates

        peroxidase_chebi = {
            "CHEBI:15882",  # phenol
            "CHEBI:16856",  # glutathione
            "CHEBI:29985",  # L-glutamate (in some peroxidases)
        }

        has_peroxidase_substrate = any(
            p.chebi_id in peroxidase_chebi
            for p in reaction.left_participants
        )

        # Apply exclusions before classifying as peroxidase

        # Exclude oxygenases (use O2, not H2O2)
        has_oxygen_substrate = any(
            p.chebi_id == CHEBI_O2
            for p in reaction.left_participants
        )

        if has_oxygen_substrate:
            return ClassificationResult(
                is_member=False,
                explanation="Uses O2 - oxygenase, not peroxidase"
            )

        # Exclude oxidases (produce H2O2, don't consume it)
        # Check if H2O2 is produced rather than consumed
        h2o2_product_count = sum(
            1 for p in reaction.right_participants
            if p.chebi_id in h2o2_indicators
        )

        if h2o2_product_count > 0:
            return ClassificationResult(
                is_member=False,
                explanation="Produces H2O2 - oxidase, not peroxidase"
            )

        # Exclude hydrolases (primarily use water as substrate, not H2O2)
        has_water_substrate = any(
            p.chebi_id == CHEBI_H2O
            for p in reaction.left_participants
        )

        if has_water_substrate and not has_peroxidase_substrate:
            return ClassificationResult(
                is_member=False,
                explanation="Hydrolase - not peroxidase"
            )

        # Try to match patterns
        match = match_patterns(reaction, self.PATTERNS, strict=False)

        explanation = "Peroxidase: H2O2-dependent oxidation"

        # Check for special catalase reaction (H2O2 → H2O + O2)
        has_oxygen_product = any(
            p.chebi_id == CHEBI_O2
            for p in reaction.right_participants
        )

        if has_oxygen_product and len(reaction.left_participants) <= 2:
            explanation += " [catalase dismutation]"

        if match.matched:
            explanation += " [pattern matched]"

        return ClassificationResult(is_member=True, explanation=explanation)
