"""Biotin Ligase reaction classification using pattern DSL.

Biotin Ligases are ultra-specific ligases that attach biotin to target proteins,
typically to form holoenzymes from apoenzymes. This is an extremely specific
post-translational modification essential for carboxylase function.

Key biochemical signature: apoprotein + biotin + ATP → biotinylated_protein + AMP + PPi

Only a few dozen reactions in all of biochemistry involve biotin ligation.
This tests our ability to classify ultra-specific cofactor dependencies.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.ligase import Ligase
from autarch.pattern_dsl import var, optional
from autarch.molecules import (
    CHEBI_AMP,
    CHEBI_ATP,
    CHEBI_DIPHOSPHATE,
    amp,
    atp,
    diphosphate,
    h_plus,
    p,
)


class BiotinBiotinCarboxylCarrierProteinLigase(Ligase):
    """biotin biotin carboxyl-carrier protein ligase

    Examples:
    - Biotin protein ligase: attaches biotin to lysine residues
    - Holoenzyme formation: apo-acetyl-CoA carboxylase → holo-enzyme

    This is ULTRA-SPECIFIC:
    - Only ~2 reactions in our dataset  
    - Requires biotin as substrate
    - Forms covalent biotin-lysine bond
    - Essential for carboxylase function
    """

    GO_ID = "GO:0004077"  # biotin-protein ligase activity
    EC_NUMBER_PREFIX = "6.3.4.-"  # Ligases forming C-N bonds (biotin attachment)

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # Classic biotin ligase: apoprotein + biotin + ATP → holoprotein + AMP + PPi
        var("apoprotein") + var("biotin") + p(atp) 
        >> var("holoprotein") + p(amp) + p(diphosphate) + optional(h_plus),
        # Alternative order
        var("biotin") + var("apoprotein") + p(atp)
        >> var("holoprotein") + p(amp) + p(diphosphate) + optional(h_plus),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a biotin ligase using ULTRA-SPECIFIC criteria.

        BIOTIN SIGNATURE:
        1. Must have biotin as substrate (extremely specific)
        2. Must use ATP → AMP + PPi energetics  
        3. Should involve protein biotinylation
        4. Should be exactly 3→3 or 3→4 reaction
        
        This is the most specific enzyme class possible!
        """
        # Must have biotin as substrate - THE KEY SPECIFICITY
        label_lower = reaction.label.lower() if reaction.label else ""
        label_parts = label_lower.split("=")
        substrate_str = label_parts[0] if label_parts else ""
        product_str = label_parts[1] if len(label_parts) > 1 else ""

        has_biotin = any(
            p.chebi_id == "CHEBI:57586"  # biotin
            for p in reaction.left_participants
        )
        
        if not has_biotin:
            return ClassificationResult(
                is_member=False,
                explanation="No biotin substrate - only biotin ligases use biotin"
            )
            
        # Must use ATP → AMP + PPi energetics (like tRNA synthetases)
        has_atp = any(
            p.chebi_id == CHEBI_ATP
            for p in reaction.left_participants
        )
        
        if not has_atp:
            return ClassificationResult(
                is_member=False,
                explanation="No ATP - biotin ligases require ATP for activation"
            )
            
        # Must produce AMP 
        has_amp = any(
            p.chebi_id == CHEBI_AMP
            for p in reaction.right_participants
        )
        
        if not has_amp:
            return ClassificationResult(
                is_member=False,
                explanation="No AMP product - biotin ligases use ATP → AMP + PPi energetics"
            )
        
        # Must produce diphosphate (PPi)
        has_diphosphate = any(
            p.chebi_id == CHEBI_DIPHOSPHATE or 
            (False)
            for p in reaction.right_participants
        )
        
        if not has_diphosphate:
            return ClassificationResult(
                is_member=False,
                explanation="No diphosphate - biotin ligases produce AMP + PPi"
            )

        # GO:0004077 is specific to protein biotinylation.
        # Exclude biotin-CoA ligase and other non-protein biotin ligations.
        has_protein_context = (
            "protein" in substrate_str
            or "protein" in product_str
            or "carboxyl-carrier" in label_lower
            or "lysyl" in label_lower
        )
        if not has_protein_context:
            return ClassificationResult(
                is_member=False,
                explanation="No protein/carboxyl-carrier context - likely non-protein biotin ligase",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Biotin Ligase: ultra-specific biotin attachment to proteins"
        )
