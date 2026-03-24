"""Carbonate dehydratase reaction classification using pattern DSL.

Carbonate dehydratases are ultra-specific lyases that catalyze the reversible
hydration of carbon dioxide to bicarbonate and protons. This is one of the
most specific and fastest enzymatic reactions known.

Key biochemical signature: CO2 + H2O ⇌ HCO3- + H+

This is THE most specific enzymatic reaction - only one substrate (CO2),
only one product pair (bicarbonate + proton), and extremely high catalytic rates.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.lyase import Lyase
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE
from autarch.pattern_dsl import var
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_ATP,
    CHEBI_CO2,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_NADPH,
    co2,
    h_plus,
    p,
    water,
)


class CarbonateDehydratase(Lyase):
    """carbonate dehydratase

    Examples:
    - Carbonic anhydrase I: cytosolic, low activity
    - Carbonic anhydrase II: cytosolic, very high activity
    - Carbonic anhydrase IV: membrane-bound, extracellular
    - Carbonic anhydrase IX: tumor-associated
    
    This is ULTIMATE SPECIFICITY:
    - Only 1 reaction type in all of biochemistry
    - Requires CO2 as the only true substrate
    - Forms bicarbonate as the only major product
    - Fastest known enzyme (kcat > 10^6 s^-1)
    """

    GO_ID = "GO:0004089"  # carbonate dehydratase activity
    EC_NUMBER_PREFIX = "4.2.1.1"  # Carbonic anhydrase
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # Forward: CO2 + H2O → HCO3- + H+
        p(co2) + p(water) >> var("bicarbonate") + p(h_plus),
        # Reverse: HCO3- + H+ → CO2 + H2O
        var("bicarbonate") + p(h_plus) >> p(co2) + p(water),
        # Simplified without explicit water
        p(co2) >> var("bicarbonate") + p(h_plus),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is carbonate dehydratase using CO2-ULTIMATE-SPECIFICITY.

        CARBONIC ANHYDRASE SIGNATURE:
        1. Must involve CO2 as substrate or product
        2. Must involve bicarbonate/hydrogen carbonate
        3. Should be 2→2 reaction (CO2 + H2O ⇌ HCO3- + H+)
        4. Should NOT involve other cofactors
        5. Should be reversible (can go both directions)
        
        This is the most specific enzyme in biochemistry!
        """
        # Must involve CO2
        has_co2 = any(
            p.chebi_id == CHEBI_CO2 for p in (reaction.left_participants + reaction.right_participants)
        )
        
        if not has_co2:
            return ClassificationResult(
                is_member=False,
                explanation="No CO2 - carbonate dehydratase requires CO2 as substrate"
            )
        
        # Must involve bicarbonate/hydrogen carbonate
        has_bicarbonate = any(
            p.chebi_id == "CHEBI:17544"  # bicarbonate
            for p in (reaction.left_participants + reaction.right_participants)
        )
        
        if not has_bicarbonate:
            return ClassificationResult(
                is_member=False,
                explanation="No bicarbonate - carbonate dehydratase produces HCO3-"
            )
        
        # Should involve proton (H+)
        has_proton = any(
            p.chebi_id == CHEBI_H_PLUS for p in (reaction.left_participants + reaction.right_participants)
        )
        
        if not has_proton:
            return ClassificationResult(
                is_member=False,
                explanation="No proton - carbonate dehydratase produces H+"
            )
        
        # Should be simple 2→2 reaction
        reactant_count = len(reaction.left_participants)
        product_count = len(reaction.right_participants)
        
        if reactant_count not in [1, 2] or product_count not in [2, 3]:
            return ClassificationResult(
                is_member=False,
                explanation=f"Wrong stoichiometry ({reactant_count}→{product_count}) - expected CO2 hydration"
            )
        
        # Should NOT involve major cofactors (this is cofactor-free)
        has_cofactors = any(
            p.chebi_id in [CHEBI_ATP, CHEBI_NADPH, CHEBI_ADP] for p in (reaction.left_participants + reaction.right_participants)
        )
        
        if has_cofactors:
            return ClassificationResult(
                is_member=False,
                explanation="Uses cofactors - carbonate dehydratase is cofactor-free"
            )
        
        # Should involve water (may be implicit)
        any(
            p.chebi_id == CHEBI_H2O for p in (reaction.left_participants + reaction.right_participants)
        )
        
        return ClassificationResult(
            is_member=True,
            explanation="Carbonate dehydratase: ultimate specificity for CO2 ⇌ HCO3- + H+"
        )
