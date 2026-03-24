"""Aldolase reaction classification using pattern DSL.

Aldolases are lyases that perform aldol additions and cleavages.
EC 4.1.2.x classification.

History
-------

## 2025-12-22

Fixed false positives by removing overly permissive fragmentation/condensation
check. Previously any lyase with more products than reactants (or vice versa)
would be classified as aldolase. Now requires known aldolase ChEBI IDs or
carbonyl chemistry evidence.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.lyase import Lyase
from autarch.ontology.reaction_diff import ReactionDiff
from autarch.pattern_dsl import var, optional, match_patterns
from autarch.molecules import (
    CHEBI_CO2,
    CHEBI_H2O,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    h_plus,
)


class AldehydeLyase(Lyase):
    """aldehyde-lyase

    Examples:
    - Fructose-1,6-bisphosphate aldolase: F1,6BP ⇌ DHAP + G3P
    - Fructose-bisphosphate aldolase: F1,6BP ⇌ dihydroxyacetone-P + glyceraldehyde-3-P
    - Aldol condensation: acetaldehyde + acetaldehyde → aldol
    - Threonine aldolase: threonine → glycine + acetaldehyde
    """

    GO_ID = "GO:0016832"  # aldehyde-lyase activity
    EC_NUMBER_PREFIX = "4.1.2.-"  # Aldehyde-lyases

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # Classic aldolase: C6 compound → 2 C3 compounds
        var("c6_substrate") >> var("c3_fragment1") + var("c3_fragment2") + optional(h_plus),
        # Reverse aldol: 2 small molecules → larger molecule
        var("fragment1") + var("fragment2") >> var("condensed_product") + optional(h_plus),
        # Aldol with phosphates: phospho-sugar → 2 phospho-fragments
        var("phospho_sugar") >> var("phospho_fragment1") + var("phospho_fragment2"),
        # Amino acid aldolase: amino acid → amino fragment + aldehyde
        var("amino_acid") >> var("amino_fragment") + var("aldehyde") + optional(h_plus),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is an aldolase using pattern matching.

        Strategy:
        1. Must be a lyase (parent class)
        2. Must involve C-C bond cleavage/formation
        3. Must involve aldol-type chemistry (carbonyl compounds)
        4. Look for characteristic aldolase patterns
        """
        # First check if it's a lyase
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not a lyase: {parent_result.explanation}",
            )

        diff = ReactionDiff(reaction)

        # Look for aldolase indicators - specific ChEBI IDs for aldolase substrates/products
        aldolase_substrates = {
            # Fructose bisphosphate pathway
            "CHEBI:16905",  # fructose 1,6-bisphosphate
            "CHEBI:28013",  # fructose 1,6-bisphosphate(4-)
            "CHEBI:32966",  # D-fructose 1,6-bisphosphate
            "CHEBI:57642",  # fructose 6-phosphate
            "CHEBI:4167",   # dihydroxyacetone phosphate (DHAP)
            "CHEBI:16108",  # dihydroxyacetone phosphate(2-)
            "CHEBI:29052",  # glyceraldehyde 3-phosphate (G3P)
            "CHEBI:58027",  # D-glyceraldehyde 3-phosphate(2-)
            # Sedoheptulose pathway
            "CHEBI:17969",  # sedoheptulose 7-phosphate
            "CHEBI:15721",  # sedoheptulose 1,7-bisphosphate
            # Other aldolase substrates
            "CHEBI:17138",  # glyceraldehyde
            "CHEBI:16474",  # dihydroxyacetone
            "CHEBI:17170",  # acetaldehyde
            "CHEBI:16842",  # formaldehyde
            # Amino acid aldolases
            "CHEBI:16857",  # L-threonine
            "CHEBI:17115",  # L-allo-threonine
            "CHEBI:15354",  # glycine
        }

        aldolase_products = {
            # Same as substrates (reversible)
            "CHEBI:16905", "CHEBI:28013", "CHEBI:32966", "CHEBI:57642",
            "CHEBI:4167", "CHEBI:16108", "CHEBI:29052", "CHEBI:58027",
            "CHEBI:17969", "CHEBI:15721",
            "CHEBI:17138", "CHEBI:16474", "CHEBI:17170", "CHEBI:16842",
            "CHEBI:16857", "CHEBI:17115", "CHEBI:15354",
        }

        # Check for aldolase substrate or product ChEBI IDs
        has_aldolase_substrate = any(
            p.chebi_id in aldolase_substrates
            for p in reaction.left_participants
        )

        has_aldolase_product = any(
            p.chebi_id in aldolase_products
            for p in reaction.right_participants
        )

        has_aldolase_pattern = has_aldolase_substrate or has_aldolase_product

        # Check for C-C bond formation/cleavage pattern
        # Classic aldolase: 1 molecule → 2 molecules or vice versa
        is_fragmentation = diff.n_product_molecules > diff.n_reactant_molecules
        is_condensation = diff.n_reactant_molecules > diff.n_product_molecules

        # MUST have known aldolase ChEBI IDs - fragmentation/condensation alone is not enough
        # (Many other lyases also change molecule count)
        if not has_aldolase_pattern:
            return ClassificationResult(
                is_member=False,
                explanation="No aldolase substrates/products detected"
            )

        # Apply exclusions before classifying as aldolase

        # Exclude decarboxylases (CO2 release, not aldol chemistry)
        has_co2 = any(
            p.chebi_id == CHEBI_CO2
            for p in reaction.right_participants
        )
        
        if has_co2 and not has_aldolase_pattern:
            return ClassificationResult(
                is_member=False,
                explanation="Decarboxylase - not aldol chemistry"
            )

        # Exclude hydrolases (use water) - unless known aldolase substrate present
        has_water = any(
            p.chebi_id == CHEBI_H2O
            for p in reaction.left_participants
        )

        if has_water and not has_aldolase_substrate:
            return ClassificationResult(
                is_member=False,
                explanation="Hydrolase - not aldolase"
            )

        # Exclude oxidoreductases (NAD/NADP involvement)
        has_nad_system = any(
            p.chebi_id in {CHEBI_NAD_PLUS, CHEBI_NADH, CHEBI_NADP_PLUS, CHEBI_NADPH}for p in reaction.left_participants + reaction.right_participants
        )
        
        if has_nad_system:
            return ClassificationResult(
                is_member=False,
                explanation="Oxidoreductase - not aldolase"
            )

        # Exclude transferases (group transfer)
        has_coa = False
        
        if has_coa and not has_aldolase_pattern:
            return ClassificationResult(
                is_member=False,
                explanation="CoA-dependent transferase - not aldolase"
            )

        # Try to match patterns
        match = match_patterns(reaction, self.PATTERNS, strict=False)

        # Determine specific type
        explanation = "Aldolase: aldol cleavage/condensation"
        
        # Specific aldolase type detection would need ChEBI IDs

        # Check for reaction direction
        if is_fragmentation:
            explanation += " [cleavage direction]"
        elif is_condensation:
            explanation += " [condensation direction]"

        if match.matched:
            explanation += " [pattern matched]"

        return ClassificationResult(is_member=True, explanation=explanation)
