"""Protein deglycase (EC 3.5.1.124, GO:0036524) classifier.

This enzyme removes methylglyoxal and glyoxal modifications from proteins,
particularly from arginine, lysine, and cysteine residues.

History
-------

## 2025-12-21

Reviewed for SMARTS conversion. Name-based patterns retained as appropriate
since classification depends on protein modification nomenclature (e.g.,
"[protein]", "1-hydroxy-2-oxopropyl") which cannot be detected structurally.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.hydrolase import Hydrolase
from autarch.pattern_dsl import var, match_patterns
from autarch.molecules import (
    CHEBI_H2O,
    CHEBI_H_PLUS,
    h_plus,
    p,
    water,
)


class ProteinDeglycase(Hydrolase):
    """protein deglycase

    Examples:
    - N(6)-(1-hydroxy-2-oxopropyl)-L-lysyl-[protein] + H2O → L-lysyl-[protein] + lactate + H+
    - S-(1-hydroxy-2-oxoethyl)-L-cysteinyl-[protein] + H2O → L-cysteinyl-[protein] + glycolate + H+
    """
    
    GO_ID = "GO:0036524"
    EC_NUMBER_PREFIX = "3.5.1.124"
    
    # Known products of deglycase activity
    DEGLYCASE_PRODUCTS = {
        "CHEBI:24996",  # lactate
        "CHEBI:16891",  # L-lactate
        "CHEBI:29805",  # glycolate
        "CHEBI:17497",  # glycolic acid
    }
    
    # Modification patterns in compound names
    MODIFICATION_PATTERNS = [
        "1-hydroxy-2-oxopropyl",  # methylglyoxal adduct
        "1-hydroxy-2-oxoethyl",    # glyoxal adduct
        "methylglyoxal",
        "glyoxal",
        "hemithioacetal",
        "hemiaminal",
    ]
    
    # Modified amino acid patterns
    MODIFIED_AA_PATTERNS = [
        r"N\(omega\).*arginyl",   # Modified arginine
        r"N\(6\).*lysyl",         # Modified lysine  
        r"S-.*cysteinyl",         # Modified cysteine
        r"N\(\w+\)-.*\[protein\]",
        r"S-.*\[protein\]",
    ]
    
    # Pattern DSL for deglycase reactions
    PATTERNS: list[Reaction] = [
        # Modified protein + H2O → protein + lactate + H+
        var("modified_protein") + p(water)
        >> var("protein") + var("lactate") + h_plus,
        
        # Modified protein + H2O → lactate + protein + H+
        var("modified_protein") + p(water)
        >> var("lactate") + var("protein") + h_plus,
        
        # With glycolate instead of lactate
        var("modified_protein") + p(water)
        >> var("protein") + var("glycolate") + h_plus,
        
        # Sometimes without explicit H+
        var("modified_protein") + p(water)
        >> var("protein") + var("product"),
    ]
    
    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a protein deglycase reaction.
        
        Strategy:
        1. Check for deglycase products (lactate/glycolate) first
        2. Look for modified protein substrates (methylglyoxal/glyoxal adducts)
        3. Check reaction label for deglycase patterns
        4. Verify water consumption
        5. Pattern matching
        """
        # Use RHEA label for pattern matching (not p.name which has ChEBI generic names)
        label_lower = reaction.label.lower() if reaction.label else ""
        label_parts = label_lower.split("=")
        substrate_str = label_parts[0] if label_parts else ""
        product_str = label_parts[1] if len(label_parts) > 1 else ""

        # First check for deglycase products (lactate or glycolate)
        has_deglycase_product = False
        product_name = None

        # Check by ChEBI ID first
        for participant in reaction.right_participants:
            if participant.chebi_id in self.DEGLYCASE_PRODUCTS:
                has_deglycase_product = True
                product_name = "lactate/glycolate"
                break

        # Fallback to label pattern
        if not has_deglycase_product:
            if "lactate" in product_str or "lactic" in product_str:
                has_deglycase_product = True
                product_name = "lactate"
            elif "glycolate" in product_str or "glycolic" in product_str:
                has_deglycase_product = True
                product_name = "glycolate"

        # If we have deglycase products, check for water and bypass parent validation
        if has_deglycase_product:
            # Check basic requirements: water consumption (ChEBI ID or SMILES)
            has_water = any(
                participant.chebi_id == CHEBI_H2O or participant.smiles == "O"
                for participant in reaction.left_participants
            ) or "h2o" in substrate_str or "water" in substrate_str

            if not has_water:
                return ClassificationResult(
                    is_member=False,
                    explanation="No water consumption for protein deglycase"
                )

            # If we have water and lactate/glycolate, this is deglycase
            explanation = "Protein deglycase: removes glycation modification"
            if product_name:
                explanation += f" → {product_name}"

            # Additional validation: check for H+ production (common in deglycase)
            has_proton = any(
                participant.chebi_id == CHEBI_H_PLUS
                for participant in reaction.right_participants
            ) or "h+" in product_str or "h(+)" in product_str

            if has_proton:
                explanation += " + H+"

            return ClassificationResult(is_member=True, explanation=explanation)
        
        # Otherwise check parent hydrolase requirements
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not a hydrolase: {parent_result.explanation}"
            )
        
        # Check for modified protein substrate (use label, not p.name)
        has_modified_protein = self._has_modified_protein_label(substrate_str)

        # Check for deglycase products via label (already checked ChEBI IDs above)
        if not has_deglycase_product:
            for participant in reaction.right_participants:
                if participant.chebi_id in self.DEGLYCASE_PRODUCTS:
                    has_deglycase_product = True
                    product_name = "lactate/glycolate"
                    break
            if not has_deglycase_product:
                if "lactate" in product_str or "lactic" in product_str:
                    has_deglycase_product = True
                    product_name = "lactate"
                elif "glycolate" in product_str or "glycolic" in product_str:
                    has_deglycase_product = True
                    product_name = "glycolate"

        # Check for protein in products (should be unmodified)
        has_protein_product = "[protein]" in product_str and not self._has_modified_protein_label(product_str)

        # Try pattern matching
        match = match_patterns(reaction, self.PATTERNS, strict=False)

        # Determine if this is deglycase activity
        if has_modified_protein and (has_deglycase_product or has_protein_product):
            explanation = "Protein deglycase: removes "

            # Identify modification type from label
            if "1-hydroxy-2-oxopropyl" in substrate_str:
                explanation += "methylglyoxal"
            elif "1-hydroxy-2-oxoethyl" in substrate_str:
                explanation += "glyoxal"
            else:
                explanation += "glycation"
            explanation += " modification from protein"

            if product_name:
                explanation += f" → {product_name}"

            if match.matched:
                explanation += " [pattern matched]"

            return ClassificationResult(is_member=True, explanation=explanation)

        # Check for general glycation repair pattern (in label)
        glycation_terms = ["glycat", "amadori", "schiff", "maillard"]
        if any(term in substrate_str for term in glycation_terms):
            return ClassificationResult(
                is_member=True,
                explanation="Protein deglycase: glycation repair"
            )

        return ClassificationResult(
            is_member=False,
            explanation="Not a protein deglycase reaction"
        )

    def _has_modified_protein_label(self, label_str: str) -> bool:
        """Check if label indicates modified protein substrate."""
        import re

        # Must contain [protein] or similar
        if "[protein]" not in label_str:
            return False

        # Check for modification patterns
        for pattern in self.MODIFICATION_PATTERNS:
            if pattern.lower() in label_str:
                return True

        # Check for modified amino acid patterns
        for pattern in self.MODIFIED_AA_PATTERNS:
            if re.search(pattern, label_str, re.IGNORECASE):
                return True

        return False
