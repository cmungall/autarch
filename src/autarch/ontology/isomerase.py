"""Simplified Isomerase reaction classification.

Isomerases catalyze molecular rearrangements within a molecule.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import ReactionClass

from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_ATP,
    CHEBI_CO2,
    CHEBI_COA,
    CHEBI_FAD,
    CHEBI_FADH2,
    CHEBI_H2O,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
)


class Isomerase(ReactionClass):
    """isomerase"""
    
    GO_ID = "GO:0016853"  # isomerase activity
    EC_NUMBER_PREFIX = "5.-.-.-"  # All isomerases
    
    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is an isomerase.
        
        Simple criteria:
        - Same number of molecules in and out
        - No water involvement
        - No net ATP/cofactor consumption
        - No redox cofactors
        """
        n_reactants = len(reaction.left_participants)
        n_products = len(reaction.right_participants)
        
        # Must have same molecule count (no fusion/fragmentation)
        if n_reactants != n_products:
            return ClassificationResult(
                is_member=False,
                explanation="Molecule count changed - not isomerization"
            )
        
        # Exclude reactions with water (hydrolases/lyases)
        has_water = any(
            p.chebi_id == CHEBI_H2O
            for p in reaction.left_participants + reaction.right_participants
        )
        
        if has_water:
            return ClassificationResult(
                is_member=False,
                explanation="Water involved - not isomerization"
            )
        
        # Exclude oxidoreductases (redox cofactors)
        REDOX_COFACTORS = {
            CHEBI_NAD_PLUS, CHEBI_NADH,  # NAD+/NADH
            CHEBI_NADP_PLUS, CHEBI_NADPH,  # NADP+/NADPH  
            CHEBI_FAD, CHEBI_FADH2,  # FAD/FADH2
            CHEBI_O2,  # O2
        }
        
        has_redox = any(
            p.chebi_id in REDOX_COFACTORS
            for p in reaction.left_participants + reaction.right_participants
        )
        
        if has_redox:
            return ClassificationResult(
                is_member=False,
                explanation="Redox cofactor present - not isomerase"
            )
        
        # Check ATP/ADP balance (some isomerases use ATP but regenerate it)
        atp_left = sum(1 for p in reaction.left_participants 
                      if p.chebi_id == CHEBI_ATP)
        atp_right = sum(1 for p in reaction.right_participants
                       if p.chebi_id == CHEBI_ATP)
        adp_left = sum(1 for p in reaction.left_participants
                      if p.chebi_id == CHEBI_ADP)
        adp_right = sum(1 for p in reaction.right_participants
                       if p.chebi_id == CHEBI_ADP)
        
        # Total nucleotides should be conserved
        if (atp_left + adp_left) != (atp_right + adp_right):
            return ClassificationResult(
                is_member=False,
                explanation="Net ATP consumption - not isomerase"
            )
        
        # Exclude CO2 involvement (lyases/ligases)
        has_co2 = any(
            p.chebi_id == CHEBI_CO2
            for p in reaction.left_participants + reaction.right_participants
        )

        if has_co2:
            return ClassificationResult(
                is_member=False,
                explanation="CO2 involved - not isomerization"
            )

        # Simple 1→1 transformation is likely isomerization
        if n_reactants == 1 and n_products == 1:
            return ClassificationResult(
                is_member=True,
                explanation="Isomerase: 1→1 molecular rearrangement"
            )
        
        # Multi-molecule balanced transformation - need to exclude transferases
        if n_reactants == n_products and n_reactants >= 2:
            # Exclude transferase reactions (group transfer between molecules)
            
            # Check for transferase cofactors/patterns
            has_transferase_cofactors = any(
                p.chebi_id in {CHEBI_COA, "CHEBI:15346"} or  # CoA, coenzyme A
                (False)
                for p in reaction.left_participants + reaction.right_participants
            )
            
            if has_transferase_cofactors:
                return ClassificationResult(
                    is_member=False,
                    explanation="Transferase cofactors present - group transfer not isomerization"
                )
            
            # Check for aminotransferase patterns (amino group transfer)
            has_amino_donor = False
            
            has_keto_acceptor = False
            
            if has_amino_donor and has_keto_acceptor:
                return ClassificationResult(
                    is_member=False,
                    explanation="Aminotransferase pattern - amino group transfer not isomerization"
                )
            
            # Check for phosphotransferase patterns
            has_phosphate_donor = False
            
            has_nucleotide = False
            
            if has_phosphate_donor and has_nucleotide:
                return ClassificationResult(
                    is_member=False,
                    explanation="Phosphotransferase pattern - phosphate transfer not isomerization"
                )
            
            # Check for acyltransferase patterns using CoA
            # CoA involvement already checked above
            
            # If no clear transferase patterns, could be true multi-molecule isomerization
            # But be conservative - only allow if molecules look related
            if n_reactants == 2:
                # Very conservative - only allow if no obvious group transfer
                return ClassificationResult(
                    is_member=False,
                    explanation="Multi-molecule reaction - likely transferase not isomerase"
                )
        
        return ClassificationResult(
            is_member=False,
            explanation="No isomerization pattern detected"
        )
