"""Simplified Oxidoreductase reaction classification.

Oxidoreductases catalyze electron transfer reactions (oxidation-reduction).
"""

from typing import ClassVar, Optional

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import ReactionClass

from autarch.molecules import (
    CHEBI_FAD,
    CHEBI_FADH2,
    CHEBI_H2O2,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
)


class Oxidoreductase(ReactionClass):
    """oxidoreductase"""
    
    GO_ID = "GO:0016491"  # oxidoreductase activity
    EC_NUMBER_PREFIX: ClassVar[Optional[str]] = "1.-.-.-"  # All oxidoreductases
    
    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is an oxidoreductase using cofactor presence.
        
        Simple criteria: Presence of redox cofactors or O2/H2O2.
        """
        # Common redox cofactors and their reduced/oxidized forms
        REDOX_COFACTORS = {
            # NAD+/NADH
            CHEBI_NAD_PLUS,  # NAD+
            CHEBI_NADH,  # NADH
            # NADP+/NADPH  
            CHEBI_NADP_PLUS,  # NADP+
            CHEBI_NADPH,  # NADPH
            # FAD/FADH2
            CHEBI_FAD,  # FAD
            CHEBI_FADH2,  # FADH2
            # Oxygen species
            CHEBI_O2,  # O2
            CHEBI_H2O2,  # H2O2
            # FMN
            "CHEBI:58210",  # FMN
            "CHEBI:58307",  # FMNH2
            # Quinones
            "CHEBI:16389",  # ubiquinone
            "CHEBI:17976",  # ubiquinol
            # Cytochromes (generic)
            "CHEBI:18070",  # iron(2+)
            "CHEBI:29034",  # iron(3+)
        }
        
        # Check for redox cofactors
        has_redox_cofactor = any(
            p.chebi_id in REDOX_COFACTORS 
            for p in reaction.left_participants + reaction.right_participants
        )
        
        if has_redox_cofactor:
            # Determine the type of oxidoreductase
            explanation = "Oxidoreductase: "
            
            # Check specific cofactors
            has_nad = any(
                p.chebi_id in {CHEBI_NAD_PLUS, CHEBI_NADH}
                for p in reaction.left_participants + reaction.right_participants
            )
            has_nadp = any(
                p.chebi_id in {CHEBI_NADP_PLUS, CHEBI_NADPH}
                for p in reaction.left_participants + reaction.right_participants
            )
            has_oxygen = any(
                p.chebi_id == CHEBI_O2
                for p in reaction.left_participants + reaction.right_participants
            )
            has_peroxide = any(
                p.chebi_id == CHEBI_H2O2
                for p in reaction.left_participants + reaction.right_participants
            )
            
            if has_nad:
                explanation += "NAD+/NADH-dependent"
            elif has_nadp:
                explanation += "NADP+/NADPH-dependent"
            elif has_oxygen:
                explanation += "O2-dependent oxidation"
            elif has_peroxide:
                explanation += "H2O2-involved"
            else:
                explanation += "electron transfer reaction"

            return ClassificationResult(
                is_member=True,
                explanation=explanation
            )

        return ClassificationResult(
            is_member=False,
            explanation="No redox cofactors or patterns detected"
        )
