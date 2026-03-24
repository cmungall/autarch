"""Ribonuclease reaction classification using pattern DSL.

Ribonucleases are ultra-specific nucleases that cleave RNA chains, either
endonucleolytically (internally) or exonucleolytically (from ends). These
enzymes are essential for RNA processing and degradation.

Key biochemical signature: RNA + H2O → RNA fragments + H+

This tests our ability to classify RNA-specific hydrolases and distinguish
them from DNases and general phosphodiesterases.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import ReactionClass
from autarch.pattern_dsl import var, optional
from autarch.molecules import (
    CHEBI_H2O,
    h_plus,
    p,
    phosphate,
    water,
)


class RNANuclease(ReactionClass):
    """RNA nuclease

    Examples:
    - RNase A: endonuclease cleaving at pyrimidines
    - RNase H: RNA-DNA hybrid specific
    - RNase III: double-stranded RNA endonuclease
    - RNase T1: guanosine-specific endonuclease
    
    This is ULTRA-SPECIFIC:
    - Only ~50-100 reactions in databases
    - Requires RNA substrates (not DNA)
    - Produces ribonucleotides (not deoxyribonucleotides)
    - Essential for RNA metabolism
    """

    GO_ID = "GO:0004540"  # ribonuclease activity
    EC_NUMBER_PREFIX = "3.1.27.-"  # Endoribonucleases producing 5'-phosphomonoesters

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # RNA hydrolysis: RNA + H2O → fragments
        var("RNA") + p(water) >> var("RNA_fragments"),
        # With phosphate production
        var("RNA") + p(water) >> var("RNA_fragments") + p(phosphate) + optional(h_plus),
        # With nucleotide production
        var("RNA") + p(water) >> var("RNA_fragment") + var("nucleotide"),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is ribonuclease using RNA-SPECIFIC criteria.

        RIBONUCLEASE SIGNATURE:
        1. Must consume H2O (hydrolytic cleavage)
        2. Should involve RNA substrates (not DNA)
        3. Should produce ribonucleotides or RNA fragments
        4. Should be 2→2 or 2→3 reaction typically
        5. May produce phosphate from backbone cleavage
        
        This tests RNA vs DNA specificity!
        """
        # Must consume H2O for hydrolysis
        has_water = any(
            p.chebi_id == CHEBI_H2O for p in reaction.left_participants
        )
        
        if not has_water:
            return ClassificationResult(
                is_member=False,
                explanation="No H2O reactant - ribonucleases require water for hydrolysis"
            )
            
        # RNA-related ChEBI IDs (ribonucleosides and ribonucleotides)
        RNA_CHEBI = {
            "CHEBI:17596",  # CMP (cytidine 5'-monophosphate)
            "CHEBI:17659",  # AMP (adenosine 5'-monophosphate)
            "CHEBI:17345",  # GMP (guanosine 5'-monophosphate)
            "CHEBI:16695",  # UMP (uridine 5'-monophosphate)
            "CHEBI:15047",  # cytidine
            "CHEBI:16335",  # adenosine
            "CHEBI:16750",  # guanosine
            "CHEBI:16704",  # uridine
        }

        has_rna_substrate = any(
            p.chebi_id in RNA_CHEBI
            for p in reaction.left_participants
        )

        has_rna_product = any(
            p.chebi_id in RNA_CHEBI
            for p in reaction.right_participants
        )

        if not (has_rna_substrate or has_rna_product):
            return ClassificationResult(
                is_member=False,
                explanation="No RNA ChEBI IDs - ribonucleases are RNA-specific"
            )
        
        # DNA-specific ChEBI IDs (deoxyribonucleosides)
        DNA_CHEBI = {
            "CHEBI:17239",  # dCMP
            "CHEBI:16890",  # dAMP
            "CHEBI:16192",  # dGMP
            "CHEBI:17013",  # dTMP
        }

        has_dna_components = any(
            p.chebi_id in DNA_CHEBI
            for p in (reaction.left_participants + reaction.right_participants)
        )

        if has_dna_components:
            return ClassificationResult(
                is_member=False,
                explanation="Contains DNA ChEBI IDs - this is likely a DNase, not RNase"
            )
        
        # Should be reasonable stoichiometry for nuclease
        reactant_count = len(reaction.left_participants)
        product_count = len(reaction.right_participants)
        
        if reactant_count not in [2, 3] or product_count not in [2, 3, 4]:
            return ClassificationResult(
                is_member=False,
                explanation=f"Wrong stoichiometry ({reactant_count}→{product_count}) - expected nuclease pattern"
            )
        
        # May produce phosphate (from RNA backbone cleavage)
        any(
            p.chebi_id == "CHEBI:18367" or  # phosphate
            (False)
            for p in reaction.right_participants
        )
        
        return ClassificationResult(
            is_member=True,
            explanation="Ribonuclease: ultra-specific RNA hydrolysis"
        )
