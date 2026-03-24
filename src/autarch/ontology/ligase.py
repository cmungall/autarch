"""Simplified Ligase reaction classification.

Ligases form bonds using ATP energy.
"""

from typing import ClassVar, Optional

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import ReactionClass

from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_AMP,
    CHEBI_ATP,
    CHEBI_CO2,
    CHEBI_DIPHOSPHATE,
    CHEBI_GDP,
    CHEBI_GTP,
    CHEBI_H2O,
    CHEBI_PHOSPHATE,
)


class Ligase(ReactionClass):
    """ligase"""
    
    GO_ID = "GO:0016874"  # ligase activity
    EC_NUMBER_PREFIX: ClassVar[Optional[str]] = "6.-.-.-"  # All ligases
    
    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a ligase.
        
        Simple criteria:
        - ATP consumption to ADP or AMP
        - Multiple substrates combining (bond formation)
        - Not simple hydrolysis
        """
        # Must have ATP or other nucleotide triphosphate as reactant
        has_ntp = any(
            p.chebi_id in {CHEBI_ATP, CHEBI_GTP} for p in reaction.left_participants
        )
        
        if not has_ntp:
            return ClassificationResult(
                is_member=False,
                explanation="No nucleotide triphosphate - ligases require ATP/GTP"
            )
        
        # Must produce nucleotide diphosphate (NTP consumed)
        has_adp = any(
            p.chebi_id == CHEBI_ADP
            for p in reaction.right_participants
        )
        has_amp = any(
            p.chebi_id == CHEBI_AMP
            for p in reaction.right_participants
        )
        has_gdp = any(
            p.chebi_id == CHEBI_GDP
            for p in reaction.right_participants
        )
        
        if not (has_adp or has_amp or has_gdp):
            return ClassificationResult(
                is_member=False,
                explanation="Nucleotide triphosphate not consumed to diphosphate/monophosphate"
            )
        
        # Exclude simple ATP hydrolysis (ATP + H2O → ADP + Pi)
        has_water = any(
            p.chebi_id == CHEBI_H2O
            for p in reaction.left_participants
        )
        
        # Count non-ATP/water substrates
        non_atp_water_reactants = [
            p for p in reaction.left_participants
            if p.chebi_id not in {CHEBI_ATP, CHEBI_H2O}
        ]
        
        if has_water and len(non_atp_water_reactants) == 0:
            # Just ATP + H2O → products (simple hydrolysis)
            return ClassificationResult(
                is_member=False,
                explanation="Simple ATP hydrolysis - not ligase"
            )
        
        # Look for bond formation signatures
        if len(non_atp_water_reactants) >= 2:
            # Multiple substrates with ATP → bond formation
            return ClassificationResult(
                is_member=True,
                explanation=f"Ligase: ATP-dependent joining of {len(non_atp_water_reactants)} substrates"
            )
        
        # Check for CO2 fixation (carboxylases)
        has_co2 = any(
            p.chebi_id == CHEBI_CO2
            for p in reaction.left_participants
        )
        
        if has_co2 and len(non_atp_water_reactants) >= 1:
            return ClassificationResult(
                is_member=True,
                explanation="Ligase: ATP-dependent carboxylation"
            )

        
        # Exclude kinase reactions (simple phosphorylation with single substrate)
        # Pattern: substrate + ATP → substrate-phosphate + ADP
        if len(non_atp_water_reactants) == 1 and has_adp:
            # Simple phosphorylation is kinase, not ligase
            return ClassificationResult(
                is_member=False,
                explanation="Single substrate + ATP → ADP pattern - likely kinase not ligase"
            )
        
        
        # Check for phosphate/diphosphate release (indicates bond formation) - more stringent
        has_phosphate = any(
            p.chebi_id in {CHEBI_PHOSPHATE, "CHEBI:18367"}
            for p in reaction.right_participants
        )
        has_diphosphate = any(
            p.chebi_id == CHEBI_DIPHOSPHATE
            for p in reaction.right_participants
        )
        
        # Only classify as ligase if multiple substrates or clear biosynthetic pattern
        if (has_phosphate or has_diphosphate) and len(non_atp_water_reactants) >= 2:
            nucleotide = "AMP + PPi" if has_amp else ("GDP + Pi" if has_gdp else "ADP + Pi")
            return ClassificationResult(
                is_member=True,
                explanation=f"Ligase: ATP → {nucleotide} with bond formation"
            )
        
        return ClassificationResult(
            is_member=False,
            explanation="No clear ligase pattern"
        )
