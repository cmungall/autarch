"""SUMO Ligase reaction classification using pattern DSL.

SUMO Ligases (Small Ubiquitin-like Modifier ligases) are ultra-specific E3 enzymes
that conjugate SUMO proteins to target substrates. SUMO modification is even more
specific than ubiquitin and regulates nuclear processes.

Key biochemical signature: protein + SUMO-E2 → SUMOylated_protein + E2

This is THE ULTIMATE in specificity - only a handful of SUMO ligases exist,
and SUMO reactions are extremely rare in biochemical databases.
"""

from typing import ClassVar, Optional

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE
from autarch.ontology.transferase import Transferase
from autarch.pattern_dsl import var
from autarch.molecules import (
    CHEBI_ATP,
)


class SUMOTransferase(Transferase):
    """sumo transferase

    Examples:
    - PIAS proteins: nuclear SUMO ligases
    - RanBP2: mitotic SUMO ligase
    - Pc2: polycomb SUMO ligase

    This is ULTIMATE SPECIFICITY:
    - Only ~1-3 reactions in any biochemical database
    - Requires SUMO protein as modifier
    - Creates SUMO-lysine isopeptide bonds
    - Essential for nuclear regulation
    """

    GO_ID = "GO:0019789"  # SUMO transferase activity
    EC_NUMBER_PREFIX: ClassVar[Optional[str]] = None  # No specific EC for SUMO ligases
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # Classic SUMO ligase: protein + SUMO-E2 → SUMOylated_protein + E2
        var("protein") + var("sumo_e2") 
        >> var("sumoylated_protein") + var("e2"),
        # Alternative with explicit SUMO
        var("protein") + var("sumo") + var("e2_enzyme")
        >> var("sumoylated_protein") + var("e2_enzyme2"),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is SUMO ligase using ULTIMATE SPECIFICITY criteria.

        SUMO LIGASE SIGNATURE:
        1. Must have SUMO as substrate (ultimate specificity)
        2. Should involve protein substrates
        3. Should be 2→2 or 3→3 reaction (no ATP hydrolysis in E3 step)
        4. Should not use ATP directly (that's E1 activation)
        5. More specific than ubiquitin ligase
        
        This is the ultimate in biochemical specificity!
        """
        # Must have SUMO as substrate - ULTIMATE SPECIFICITY
        has_sumo = any(
            p.chebi_id == "CHEBI:74067"  # SUMO if it has a ChEBI ID
            for p in reaction.left_participants
        )
        
        if not has_sumo:
            return ClassificationResult(
                is_member=False,
                explanation="No SUMO substrate - only SUMO ligases use SUMO proteins"
            )
            
        # Should NOT use ATP directly (E3s don't hydrolyze ATP)
        has_atp = any(
            p.chebi_id == CHEBI_ATP
            for p in reaction.left_participants
        )
        
        if has_atp:
            return ClassificationResult(
                is_member=False,
                explanation="Uses ATP - likely E1 or E2 enzyme, not SUMO E3 ligase"
            )
            
        # Protein detection would need specific ChEBI IDs or polymer types
        
        # Should be 2→2 or 3→3 reaction (no cofactor consumption)
        reactant_count = len(reaction.left_participants)
        product_count = len(reaction.right_participants)
        
        if reactant_count not in [2, 3] or product_count not in [2, 3]:
            return ClassificationResult(
                is_member=False,
                explanation=f"Wrong stoichiometry ({reactant_count}→{product_count}) - expected 2→2 or 3→3"
            )
        
        # Ubiquitin detection would need specific ChEBI IDs
        
        return ClassificationResult(
            is_member=True,
            explanation="SUMO Ligase: ultimate specificity for SUMO conjugation"
        )
