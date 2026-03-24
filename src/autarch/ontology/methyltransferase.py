"""Methyltransferase reaction classification using pattern DSL.

Methyltransferases are a specific class of transferases that transfer methyl groups,
typically from S-adenosyl-L-methionine (SAM) to various acceptor molecules.

Key biochemical signature: SAM + acceptor → SAH + methylated_acceptor
This is one of the most common and well-characterized enzyme reaction types.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE
from autarch.ontology.transferase import Transferase
from autarch.pattern_dsl import var, optional
from autarch.molecules import (
    CHEBI_H2O,
    CHEBI_SAH,
    CHEBI_SAM,
    h_plus,
    p,
    sah,
    sam,
)


class Methyltransferase(Transferase):
    """methyltransferase

    Examples:
    - DNA methyltransferases: cytosine methylation
    - Histone methyltransferases: chromatin modification  
    - Small molecule methyltransferases: creatine synthesis
    - O-methyltransferases: catechol methylation
    - N-methyltransferases: neurotransmitter metabolism
    """

    GO_ID = "GO:0008168"  # methyltransferase activity
    EC_NUMBER_PREFIX = "2.1.1.-"  # Transferases transferring one-carbon groups (methyltransferases)
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    # Define patterns using DSL - SAM is the universal methyl donor
    PATTERNS: list[Reaction] = [
        # Classic methyltransferase: SAM + acceptor → SAH + product + H+
        p(sam) + var("acceptor") >> p(sah) + var("product") + optional(h_plus),
        # Without explicit H+
        p(sam) + var("acceptor") >> p(sah) + var("product"),
        # With multiple acceptors (some methyltransferases have complex substrates)
        p(sam) + var("acceptor1") + var("acceptor2") 
        >> p(sah) + var("product1") + var("product2") + optional(h_plus),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a methyltransferase using SIMPLE, DECLARATIVE approach.

        BIOCHEMICAL SIGNATURE:
        1. Must have SAM (S-adenosyl-L-methionine) as methyl donor
        2. Must produce SAH (S-adenosyl-L-homocysteine) as co-product
        3. Should not involve water (not hydrolytic)
        4. Should not involve CoA (different transferase type)
        
        This captures the core biochemical mechanism without complex heuristics.
        """
        # Must have SAM as methyl donor
        has_sam = any(
            p.chebi_id == CHEBI_SAM
            for p in reaction.left_participants
        )
        
        if not has_sam:
            return ClassificationResult(
                is_member=False,
                explanation="No SAM (S-adenosyl-L-methionine) - not SAM-dependent methyltransferase"
            )
            
        # Must produce SAH as co-product
        has_sah = any(
            p.chebi_id == CHEBI_SAH
            for p in reaction.right_participants
        )
        
        if not has_sah:
            return ClassificationResult(
                is_member=False,
                explanation="No SAH (S-adenosyl-L-homocysteine) product - not methyl transfer"
            )
        
        # Should not involve water (methylation is not hydrolytic)
        has_water = any(
            p.chebi_id == CHEBI_H2O
            for p in reaction.left_participants + reaction.right_participants
        )
        
        if has_water:
            return ClassificationResult(
                is_member=False,
                explanation="Water involved - methyl transfer is not hydrolytic"
            )
        
        # Should not involve CoA (that's acyltransferase territory)
        has_coa = any(
            False
            for p in reaction.left_participants + reaction.right_participants
        )
        
        if has_coa:
            return ClassificationResult(
                is_member=False,
                explanation="CoA involved - likely acyltransferase not methyltransferase"
            )
        
        return ClassificationResult(
            is_member=True,
            explanation="Methyltransferase: SAM-dependent methyl group transfer"
        )
