"""Sulfotransferase reaction classification using pattern DSL.

Sulfotransferases transfer sulfate groups from PAPS to acceptor molecules.
EC 2.8.2.x classification.

History
-------

## 2025-12-22 (v2)

Fixed false positives by:
1. Excluding PAPS hydrolysis (PAPS + H2O → PAP + sulfate) - this is NOT sulfotransfer
2. Excluding CoA transfers that produce PAP (phosphopantetheinyl transfers)
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.transferase import Transferase
from autarch.pattern_dsl import var, optional, match_patterns
from autarch.molecules import (
    CHEBI_ATP,
    CHEBI_H2O,
    CHEBI_PHOSPHATE,
    h_plus,
    p,
)


class Sulfotransferase(Transferase):
    """sulfotransferase

    Examples:
    - Estrogen sulfotransferase: estradiol + PAPS → estradiol sulfate + PAP
    - Heparan sulfate sulfotransferase: heparan + PAPS → heparan sulfate + PAP
    - Tyrosine sulfotransferase: tyrosine + PAPS → tyrosine sulfate + PAP
    - Phenol sulfotransferase: phenol + PAPS → phenyl sulfate + PAP
    """

    GO_ID = "GO:0008146"  # sulfotransferase activity
    EC_NUMBER_PREFIX = "2.8.2.-"  # Sulfotransferases

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # PAPS-dependent sulfation: substrate + PAPS → sulfated-product + PAP
        var("substrate") + var("paps") >> var("sulfated_product") + var("pap") + optional(h_plus),
        # General sulfotransferase: donor + acceptor → modified-donor + sulfated-acceptor
        var("sulfate_donor") + var("acceptor") >> var("modified_donor") + var("sulfated_acceptor"),
        # With proton: substrate + PAPS + H+ → sulfated-product + PAP
        var("substrate") + var("paps") + p(h_plus) >> var("sulfated_product") + var("pap"),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a sulfotransferase using pattern matching.

        Strategy:
        1. Must be a transferase (parent class)
        2. Must involve PAPS (3'-phosphoadenosine-5'-phosphosulfate) as sulfate donor
        3. Must produce PAP (3'-phosphoadenosine-5'-phosphate) as product
        4. Look for sulfation patterns
        """
        # Use RHEA label for pattern matching
        label_lower = reaction.label.lower() if reaction.label else ""
        label_parts = label_lower.split("=")
        product_str = label_parts[1] if len(label_parts) > 1 else ""

        # First check if it's a transferase
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not a transferase: {parent_result.explanation}",
            )

        # Look for PAPS (sulfate donor)
        paps_indicators = {
            "CHEBI:58339",  # PAPS (3'-phosphoadenosine-5'-phosphosulfate)
        }


        has_paps = any(
            p.chebi_id in paps_indicators
            for p in reaction.left_participants
        )

        # Look for PAP (sulfate transfer product)
        pap_indicators = {
            "CHEBI:58343",  # PAP (3'-phosphoadenosine-5'-phosphate)
        }


        has_pap = any(
            p.chebi_id in pap_indicators
            for p in reaction.right_participants
        )

        # Must have PAPS → PAP pattern
        if not (has_paps or has_pap):
            return ClassificationResult(
                is_member=False,
                explanation="No PAPS/PAP pattern - not sulfotransferase"
            )

        # Apply exclusions before classifying as sulfotransferase

        # Exclude PAPS hydrolysis: PAPS + H2O → PAP + sulfate
        # This is NOT sulfotransfer - the sulfate goes to water, not an acceptor
        has_water = any(
            p.chebi_id == CHEBI_H2O
            for p in reaction.left_participants
        )

        CHEBI_SULFATE = "CHEBI:16189"
        has_sulfate_product = any(
            p.chebi_id == CHEBI_SULFATE
            for p in reaction.right_participants
        ) or "sulfate" in product_str

        if has_paps and has_water and has_sulfate_product:
            return ClassificationResult(
                is_member=False,
                explanation="PAPS hydrolysis - not sulfotransfer to acceptor"
            )

        # Exclude CoA/ACP transfers that happen to produce PAP
        # (phosphopantetheinyl transfers)
        CHEBI_COA = "CHEBI:57287"
        has_coa = any(
            p.chebi_id == CHEBI_COA
            for p in reaction.left_participants + reaction.right_participants
        ) or "coa" in label_lower or "acp" in label_lower

        if has_coa:
            return ClassificationResult(
                is_member=False,
                explanation="CoA/ACP involved - not sulfotransferase"
            )

        # Exclude ATP-dependent sulfate activation (sulfate + ATP → APS/PAPS)
        has_atp = any(p.chebi_id == CHEBI_ATP for p in reaction.left_participants)
        has_sulfate_substrate = any(
            p.chebi_id == "CHEBI:16189"
            for p in reaction.left_participants
        )
        
        if has_atp and has_sulfate_substrate:
            return ClassificationResult(
                is_member=False,
                explanation="Sulfate activation - not sulfotransferase"
            )

        # Exclude phosphotransferase reactions
        has_phosphate = any(
            p.chebi_id == CHEBI_PHOSPHATE
            for p in reaction.left_participants + reaction.right_participants
        )
        
        if has_phosphate and not (has_paps or has_pap):
            return ClassificationResult(
                is_member=False,
                explanation="Phosphotransferase - not sulfotransferase"
            )

        # Try to match patterns
        match = match_patterns(reaction, self.PATTERNS, strict=False)

        explanation = "Sulfotransferase: sulfate group transfer"

        # Check for PAPS dependency
        if has_paps:
            explanation += " [PAPS-dependent]"

        if match.matched:
            explanation += " [pattern matched]"

        return ClassificationResult(is_member=True, explanation=explanation)
