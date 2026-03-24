"""Superoxide Dismutase reaction classification using pattern DSL.

Superoxide dismutases are ultra-specific oxidoreductases that catalyze the
disproportionation of superoxide radicals to oxygen and hydrogen peroxide.
This is a critical antioxidant enzyme.

Key biochemical signature: 2 O2•- + 2 H+ → O2 + H2O2

This is one of the most specific and important enzymatic reactions - only
one substrate (superoxide), unique disproportionation mechanism.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.oxidoreductase import Oxidoreductase
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE
from autarch.pattern_dsl import var
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_ATP,
    CHEBI_H2O2,
    CHEBI_H_PLUS,
    CHEBI_NADPH,
    CHEBI_O2,
    h_plus,
    oxygen,
    p,
)


class SuperoxideDismutase(Oxidoreductase):
    """superoxide dismutase

    Examples:
    - SOD1: Cu/Zn superoxide dismutase (cytosolic)
    - SOD2: Mn superoxide dismutase (mitochondrial)
    - SOD3: extracellular Cu/Zn superoxide dismutase

    This is ULTIMATE SPECIFICITY:
    - Only 1 reaction type in biochemistry
    - Requires superoxide as substrate
    - Produces both O2 and H2O2 (disproportionation)
    - Essential for oxidative stress protection
    """

    GO_ID = "GO:0004784"  # superoxide dismutase activity
    EC_NUMBER_PREFIX = "1.15.1.1"  # Superoxide dismutase
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # Classic SOD: 2 O2•- + 2 H+ → O2 + H2O2
        var("superoxide") + p(h_plus) >> p(oxygen) + var("hydrogen_peroxide"),
        # Simplified without explicit stoichiometry
        var("superoxide") >> p(oxygen) + var("hydrogen_peroxide"),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is superoxide dismutase using SUPEROXIDE-ULTIMATE-SPECIFICITY.

        SUPEROXIDE DISMUTASE SIGNATURE:
        1. Must involve superoxide radical (O2•-)
        2. Must produce both O2 and H2O2 (disproportionation)
        3. Should consume H+ (protonation)
        4. Should be 2→2 or 3→3 reaction
        5. Should NOT require cofactors in the reaction
        
        This is the most specific antioxidant enzyme!
        """
        # Must involve superoxide radical
        has_superoxide = any(
            p.chebi_id == "CHEBI:18421"  # superoxide
            for p in reaction.left_participants
        )
        
        if not has_superoxide:
            return ClassificationResult(
                is_member=False,
                explanation="No superoxide substrate - SOD requires superoxide radical"
            )
        
        # Must produce O2 (molecular oxygen)
        has_oxygen_product = any(
            p.chebi_id == CHEBI_O2 for p in reaction.right_participants
        )
        
        if not has_oxygen_product:
            return ClassificationResult(
                is_member=False,
                explanation="No O2 product - SOD produces molecular oxygen"
            )
        
        # Must produce H2O2 (hydrogen peroxide)
        has_hydrogen_peroxide = any(
            p.chebi_id == CHEBI_H2O2
            for p in reaction.right_participants
        )
        
        if not has_hydrogen_peroxide:
            return ClassificationResult(
                is_member=False,
                explanation="No H2O2 product - SOD produces hydrogen peroxide"
            )
        
        # Should involve proton consumption
        has_proton = any(
            p.chebi_id == CHEBI_H_PLUS for p in reaction.left_participants
        )
        
        if not has_proton:
            return ClassificationResult(
                is_member=False,
                explanation="No proton consumption - SOD requires H+ for disproportionation"
            )
        
        # Should NOT involve major cofactors (metal cofactors are part of enzyme, not reaction)
        has_external_cofactors = any(
            p.chebi_id in [CHEBI_ATP, CHEBI_NADPH, CHEBI_ADP] for p in (reaction.left_participants + reaction.right_participants)
        )
        
        if has_external_cofactors:
            return ClassificationResult(
                is_member=False,
                explanation="Uses external cofactors - SOD is cofactor-independent for reaction"
            )
        
        return ClassificationResult(
            is_member=True,
            explanation="Superoxide Dismutase: ultimate specificity for superoxide disproportionation"
        )
