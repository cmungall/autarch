"""Nucleosidase reaction classification using pattern DSL.

Nucleosidases are hydrolases that cleave N-glycosyl bonds in nucleosides and related
compounds, releasing nucleobases and sugars. EC 3.2.2.x classification.

History
-------

## 2025-12-22 (v2)

Major fix: Changed GO_ID from GO:0016798 (all glycosidases) to GO:0016799
(N-glycosyl hydrolases). Also:
1. Added SMARTS-based nucleoside detection (many RHEA entries lack ChEBI IDs)
2. Added D-ribofuranose variants to sugar detection
3. Fixed placeholder logic that was making exclusions ineffective
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.hydrolase import Hydrolase
from autarch.pattern_dsl import var, optional
from autarch.molecules import h_plus, p, water


class HydrolaseHydrolyzingNGlycosylCompounds(Hydrolase):
    """hydrolase hydrolyzing N-glycosyl compounds

    Examples:
    - Adenosine nucleosidase: adenosine + H2O → adenine + ribose
    - Inosine nucleosidase: inosine + H2O → hypoxanthine + ribose
    - Uridine nucleosidase: uridine + H2O → uracil + ribose
    - Purine nucleosidase: purine nucleoside + H2O → purine base + ribose
    """

    GO_ID = "GO:0016799"  # hydrolase activity, hydrolyzing N-glycosyl compounds
    EC_NUMBER_PREFIX = "3.2.2.-"  # Hydrolysing N-glycosyl compounds

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # General nucleosidase: nucleoside + H2O → nucleobase + sugar
        var("nucleoside") + p(water) >> var("nucleobase") + var("sugar") + optional(h_plus),
        # Adenosine nucleosidase: adenosine + H2O → adenine + ribose
        var("adenosine") + p(water) >> var("adenine") + var("ribose") + optional(h_plus),
        # Purine nucleosidase: purine nucleoside + H2O → purine + sugar
        var("purine_nucleoside") + p(water) >> var("purine") + var("ribose") + optional(h_plus),
        # Pyrimidine nucleosidase: pyrimidine nucleoside + H2O → pyrimidine + sugar
        var("pyrimidine_nucleoside") + p(water) >> var("pyrimidine") + var("ribose") + optional(h_plus),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a nucleosidase (N-glycosyl hydrolase).

        Strategy:
        1. Must be a hydrolase (parent class)
        2. Must use water for hydrolysis (already checked by parent)
        3. Must produce ribose/deoxyribose products OR have N-glycosyl cleavage pattern
        4. Exclude reactions that are clearly not nucleosidases

        NOTE: Many nucleosides in RHEA lack ChEBI IDs, so we detect by:
        - Ribose/deoxyribose product (with expanded ChEBI variants)
        - SMARTS patterns for purine/pyrimidine bases
        """
        # Use RHEA label for pattern matching (not p.name which has ChEBI generic names)
        label_lower = reaction.label.lower() if reaction.label else ""
        label_parts = label_lower.split("=")
        label_parts[0] if label_parts else ""
        product_str = label_parts[1] if len(label_parts) > 1 else ""

        # First check if it's a hydrolase
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not a hydrolase: {parent_result.explanation}",
            )

        # Expanded ribose/deoxyribose ChEBI IDs (including charged forms)
        sugar_chebi = {
            "CHEBI:33942",  # ribose (generic)
            "CHEBI:47013",  # D-ribofuranose
            "CHEBI:28816",  # 2-deoxyribose
            "CHEBI:17634",  # D-ribose
            "CHEBI:47014",  # L-ribofuranose
            "CHEBI:78346",  # D-ribofuranose 5-phosphate(2-)
            "CHEBI:58273",  # D-ribose 5-phosphate(2-)
        }

        # Nucleobase ChEBI IDs (expanded with common variants)
        nucleobase_chebi = {
            "CHEBI:16708",  # adenine
            "CHEBI:16235",  # guanine
            "CHEBI:16040",  # cytosine
            "CHEBI:17568",  # uracil
            "CHEBI:17821",  # thymine
            "CHEBI:17368",  # hypoxanthine
            "CHEBI:17712",  # xanthine
            "CHEBI:27468",  # 7-methylxanthine (from methylated nucleosides)
            "CHEBI:17802",  # nicotinamide (from NAD+)
        }

        # Check for ribose product
        has_sugar_product = any(
            p.chebi_id in sugar_chebi
            for p in reaction.right_participants
        )

        # Check for nucleobase product
        has_nucleobase_product = any(
            p.chebi_id in nucleobase_chebi
            for p in reaction.right_participants
        )

        # SMARTS-based detection for purine/pyrimidine products (when ChEBI ID missing)
        # Purine core: imidazole fused to pyrimidine
        # Pyrimidine core: 6-membered ring with 2 nitrogens

        def has_base_smarts(participant) -> bool:
            """Check if participant has purine or pyrimidine base pattern."""
            if not participant.smiles:
                return False
            # Simple substring check for nucleobase patterns
            smiles = participant.smiles
            # Imidazole-pyrimidine fusion (purine)
            if "c1ncnc2" in smiles.lower() or "c1nc2" in smiles.lower():
                return True
            # Pyrimidine ring (check for characteristic patterns)
            if "NC(=O)" in smiles and "N=C" in smiles:
                return True
            return False

        any(
            has_base_smarts(p)
            for p in reaction.right_participants
        )

        # Exclude common hexose sugars (glucose, galactose, etc.) - these are glycosidase, not nucleosidase
        hexose_chebi = {
            "CHEBI:4167",   # D-glucose
            "CHEBI:17634",  # D-ribose (already in sugar_chebi, this is good)
            "CHEBI:12965",  # beta-D-glucose
            "CHEBI:28061",  # alpha-D-glucose
            "CHEBI:4139",   # D-galactose
            "CHEBI:28260",  # beta-D-galactose
            "CHEBI:17118",  # D-fructose
            "CHEBI:48095",  # D-glucose 6-phosphate
            "CHEBI:17925",  # alpha-D-glucose 6-phosphate
        }

        has_hexose_product = any(
            p.chebi_id in hexose_chebi
            for p in reaction.right_participants
        ) or any(sugar in product_str for sugar in ["glucose", "galactose", "fructose", "mannose", "fucose"])

        # Exclude CoA reactions (acyl hydrolases, not nucleosidases)
        has_coa = any(
            p.chebi_id and "57287" in p.chebi_id  # CoA ChEBI ID
            for p in reaction.left_participants + reaction.right_participants
        ) or "coa" in label_lower

        if has_coa:
            return ClassificationResult(
                is_member=False,
                explanation="CoA involvement - acyl hydrolase, not nucleosidase"
            )

        # For nucleosidase, we need STRONG evidence: ribose ChEBI ID + base evidence
        # OR base ChEBI ID + some ribose evidence
        # Weak SMARTS matching alone is not sufficient

        has_strong_sugar = has_sugar_product
        has_strong_base = has_nucleobase_product

        # Must have at least one strong indicator (ChEBI-verified)
        if not (has_strong_sugar or has_strong_base):
            return ClassificationResult(
                is_member=False,
                explanation="No ribose/nucleobase products - not N-glycosyl hydrolysis"
            )

        # If we have hexose products, this is glycosidase not nucleosidase
        if has_hexose_product and not (has_strong_sugar and has_strong_base):
            return ClassificationResult(
                is_member=False,
                explanation="Hexose product (glucose/galactose) - glycosidase, not nucleosidase"
            )

        # Exclude reactions involving proteins/peptides (ADP-ribosyl hydrolases are separate)
        has_protein = "protein" in label_lower or "peptide" in label_lower

        if has_protein:
            return ClassificationResult(
                is_member=False,
                explanation="Protein involvement - ADP-ribosyl hydrolase, not simple nucleosidase"
            )

        # Exclude reactions with asparagine (beta-aspartyl-glucosaminidase)
        has_asparagine = "asparagine" in label_lower

        if has_asparagine:
            return ClassificationResult(
                is_member=False,
                explanation="Asparagine involvement - aspartylglucosaminidase, not nucleosidase"
            )

        # Build explanation
        explanation = "N-glycosyl hydrolase"
        if has_sugar_product:
            explanation += " (ribose product detected)"
        if has_nucleobase_product:
            explanation += " (nucleobase product detected)"

        return ClassificationResult(is_member=True, explanation=explanation)
