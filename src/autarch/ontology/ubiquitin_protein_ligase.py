"""Ubiquitin protein ligase reaction classification using pattern DSL.

Ubiquitin Ligases (E3 enzymes) are the final step in ubiquitin conjugation,
transferring ubiquitin from E2 enzymes to target proteins. This is a critical
post-translational modification for protein degradation and signaling.

Key biochemical signature: protein + ubiquitin-E2 → ubiquitinated_protein + E2

This is ULTRA-SPECIFIC - only a handful of reactions involve ubiquitin attachment.
Tests our ability to classify highly specific protein modifications.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE
from autarch.ontology.transferase import Transferase
from autarch.pattern_dsl import var
from autarch.molecules import (
    CHEBI_ATP,
)


class UbiquitinProteinLigase(Transferase):
    """ubiquitin protein ligase

    Examples:
    - MDM2: ubiquitinates p53 for degradation
    - RING E3s: transfer ubiquitin from E2 to substrate
    - HECT E3s: form intermediate with ubiquitin before transfer

    This is ULTRA-SPECIFIC:
    - Only ~5-10 reactions in biochemical databases
    - Requires ubiquitin as modification group
    - Creates isopeptide bonds (Lys-Gly linkage)
    - Central to proteostasis and signaling
    """

    GO_ID = "GO:0061630"  # ubiquitin protein ligase activity
    EC_NUMBER_PREFIX = "2.3.2.-"  # Aminoacyltransferases (includes ubiquitin transfer)
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # Classic ubiquitin ligase: protein + ubiquitin-E2 → ubiquitinated_protein + E2
        var("protein") + var("ubiquitin_e2") 
        >> var("ubiquitinated_protein") + var("e2"),
        # Alternative with explicit ubiquitin
        var("protein") + var("ubiquitin") + var("e2_enzyme")
        >> var("ubiquitinated_protein") + var("e2_enzyme2"),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a ubiquitin protein ligase using ULTRA-SPECIFIC criteria.

        UBIQUITIN SIGNATURE:
        1. Must have ubiquitin as substrate (extremely specific)
        2. Should involve protein substrates
        3. Should be 2→2 or 3→3 reaction (no ATP hydrolysis in E3 step)
        4. Should not use ATP directly (that's E1 activation)
        
        This is one of the most specific protein modifications!
        """
        # Must have ubiquitin as substrate - THE KEY SPECIFICITY
        has_ubiquitin = any(
            p.chebi_id == "CHEBI:29801"  # ubiquitin
            for p in reaction.left_participants
        )
        
        if not has_ubiquitin:
            return ClassificationResult(
                is_member=False,
                explanation="No ubiquitin substrate - only ubiquitin protein ligases use ubiquitin"
            )
            
        # Should NOT use ATP directly (E3s don't hydrolyze ATP)
        has_atp = any(
            p.chebi_id == CHEBI_ATP
            for p in reaction.left_participants
        )
        
        if has_atp:
            return ClassificationResult(
                is_member=False,
                explanation="Uses ATP - likely E1 or E2 enzyme, not E3 ligase"
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
        
        return ClassificationResult(
            is_member=True,
            explanation="Ubiquitin protein ligase: ultra-specific ubiquitin conjugation to proteins"
        )
