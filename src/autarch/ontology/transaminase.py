"""Aminotransferase reaction classification using pattern DSL.

Aminotransferases (transaminases) are transferases that move amino groups.

NOTE: Previous complex implementation had 0% precision (30 FPs, 0 TPs, 11 FNs).
Simplified to focus on core aminotransferase pattern: amino acid + keto acid exchange.
Pattern was too broad and missed actual aminotransferases. Staying SIMPLE and DECLARATIVE.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.transferase import Transferase
from autarch.pattern_dsl import var
from autarch.molecules import (
    CHEBI_H2O,
    alanine,
    alpha_ketoglutarate,
    aspartate,
    glutamate,
    oxaloacetate,
    p,
    pyruvate,
)


class Transaminase(Transferase):
    """transaminase

    Examples:
    - ALT (alanine aminotransferase): alanine + α-ketoglutarate → pyruvate + glutamate
    - AST (aspartate aminotransferase): aspartate + α-ketoglutarate → oxaloacetate + glutamate
    - Requires PLP (pyridoxal phosphate) as cofactor (usually regenerated)
    """

    GO_ID = "GO:0008483"  # transaminase activity
    EC_NUMBER_PREFIX = "2.6.1.-"

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # Classic aminotransferase: amino_acid1 + keto_acid1 → keto_acid2 + amino_acid2
        var("amino_acid1") + var("keto_acid1")
        >> var("keto_acid2") + var("amino_acid2"),
        # ALT: alanine + α-ketoglutarate → pyruvate + glutamate
        p(alanine) + p(alpha_ketoglutarate) >> p(pyruvate) + p(glutamate),
        # AST: aspartate + α-ketoglutarate → oxaloacetate + glutamate
        p(aspartate) + p(alpha_ketoglutarate) >> p(oxaloacetate) + p(glutamate),
        # With cofactor (PLP usually regenerated so not shown): aa1 + ka1 + cofactor → ka2 + aa2 + cofactor
        var("amino_acid1") + var("keto_acid1") + var("cofactor")
        >> var("keto_acid2") + var("amino_acid2") + var("cofactor"),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is an aminotransferase using SIMPLE pattern matching.

        SIMPLE APPROACH: Look for classic aminotransferase signature:
        - 2-oxoglutarate as co-substrate (universal amino donor/acceptor) 
        - Glutamate as co-product (from 2-oxoglutarate)
        - No water involvement (not hydrolytic)
        - 2 reactants → 2 products (simple exchange)
        
        This catches the core biochemical signature without complex heuristics.
        """
        # Must have 2-oxoglutarate (α-ketoglutarate) - the universal amino acceptor/donor
        has_oxoglutarate = any(
            p.chebi_id == "CHEBI:16810"
            for p in reaction.left_participants
        )
        
        if not has_oxoglutarate:
            return ClassificationResult(
                is_member=False,
                explanation="No 2-oxoglutarate - not classic aminotransferase"
            )
            
        # Must have glutamate as product (from 2-oxoglutarate + amino group)
        has_glutamate_product = any(
            p.chebi_id == "CHEBI:29985"
            for p in reaction.right_participants
        )
        
        if not has_glutamate_product:
            return ClassificationResult(
                is_member=False,
                explanation="No glutamate product - not classic aminotransferase"
            )
        
        # Must be 2→2 reaction (amino acid exchange, not synthesis/degradation)
        if len(reaction.left_participants) != 2 or len(reaction.right_participants) != 2:
            return ClassificationResult(
                is_member=False,
                explanation="Not 2→2 reaction - aminotransferases exchange amino groups"
            )
            
        # Must not involve water (not hydrolytic)
        has_water = any(
            p.chebi_id == CHEBI_H2O
            for p in reaction.left_participants + reaction.right_participants
        )
        
        if has_water:
            return ClassificationResult(
                is_member=False,
                explanation="Water involved - amino transfer is not hydrolytic"
            )
        
        return ClassificationResult(
            is_member=True,
            explanation="Aminotransferase: 2-oxoglutarate + amino acid → glutamate + keto acid"
        )
