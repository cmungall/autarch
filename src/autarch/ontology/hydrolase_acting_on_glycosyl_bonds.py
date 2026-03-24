"""Glycosidase reaction classification using pattern DSL.

Glycosidases are hydrolases that cleave glycosidic bonds in carbohydrates.
EC 3.2.x.x classification.

History
-------

## 2025-12-23

Added exclusions for non-glycosidase reactions that have sugar-like substrates:
- Deacetylase (acetate product)
- Lyases (pyruvate product)
- Asparaginase (aspartate product)

## 2025-12-22

Fixed false positives by excluding reactions that release:
- Sulfate (sulfatase, not glycosidase)
- Phosphate (phosphatase, not glycosidase)
- Fatty acids (esterase, not glycosidase)
These were matching because substrates contained sugar moieties but the
bond being hydrolyzed wasn't glycosidic.

Previously expanded sugar_indicators set with more ChEBI IDs found in RHEA database:
- N-acetyl-D-glucosamine (CHEBI:506227), beta-D-galactose (CHEBI:27667)
- L-fucose (CHEBI:18079), D-ribose (CHEBI:47013)
- N-acetylneuraminic acid (CHEBI:35418), glucuronic acid (CHEBI:4178)
- D-glucosamine (CHEBI:5417), D-galactosamine (CHEBI:27903)

## 2025-12-21

Reviewed for SMARTS conversion. Added Participant.is_glycoside() check for
structural glycosidic bond detection. ChEBI ID detection retained as primary
method since many glycosides share common sugar product ChEBIs.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.hydrolase import Hydrolase
from autarch.pattern_dsl import var, optional, match_patterns
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_GDP,
    CHEBI_H2O,
    CHEBI_PHOSPHATE,
    CHEBI_UDP,
    h_plus,
    p,
    water,
)


class HydrolaseActingOnGlycosylBonds(Hydrolase):
    """hydrolase acting on glycosyl bonds

    Examples:
    - Alpha-amylase: starch + H2O → oligosaccharides
    - Beta-glucosidase: cellobiose + H2O → 2 glucose
    - Sucrase: sucrose + H2O → glucose + fructose
    - Lysozyme: peptidoglycan + H2O → fragments
    """

    GO_ID = "GO:0016798"  # hydrolase activity, acting on glycosyl bonds
    EC_NUMBER_PREFIX = "3.2.-.-"  # Glycosidases

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # Simple glycoside hydrolysis: glycoside + H2O → sugar + aglycon
        var("glycoside") + p(water) >> var("sugar") + var("aglycon") + optional(h_plus),
        # Polysaccharide hydrolysis: polysaccharide + H2O → oligosaccharides
        var("polysaccharide") + p(water) >> var("product1") + var("product2") + optional(h_plus),
        # Multiple water molecules: substrate + n H2O → multiple products
        var("substrate") + var("waters") >> var("product1") + var("product2") + var("product3") + optional(h_plus),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a glycosidase using pattern matching.

        Strategy:
        1. Must be a hydrolase (parent class)
        2. Must involve carbohydrate/sugar substrates or products
        3. Must use water for hydrolysis
        4. Check for glycosidic bond cleavage patterns
        """
        # Use RHEA label for pattern matching (not p.name which has ChEBI generic names)
        label_lower = reaction.label.lower() if reaction.label else ""
        label_parts = label_lower.split("=")
        substrate_str = label_parts[0] if label_parts else ""
        product_str = label_parts[1] if len(label_parts) > 1 else ""

        # First check if it's a hydrolase
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not a hydrolase: {parent_result.explanation}",
            )

        # Must use water as reactant
        has_water = any(
            p.chebi_id == CHEBI_H2O
            for p in reaction.left_participants
        )

        if not has_water:
            return ClassificationResult(
                is_member=False,
                explanation="No water reactant - not glycosidic hydrolysis"
            )

        # Look for carbohydrate/sugar indicators
        sugar_indicators = {
            # Common monosaccharides
            "CHEBI:17234",  # D-glucose
            "CHEBI:4167",   # D-glucose (alternate)
            "CHEBI:15903",  # D-fructose
            "CHEBI:16551",  # D-galactose
            "CHEBI:27667",  # beta-D-galactose
            "CHEBI:17218",  # D-mannose
            "CHEBI:26667",  # D-xylose
            "CHEBI:28645",  # L-arabinose
            "CHEBI:47013",  # D-ribose
            "CHEBI:18079",  # L-fucose
            # Amino sugars
            "CHEBI:506227",  # N-acetyl-D-glucosamine (very common in RHEA)
            "CHEBI:28009",   # N-acetyl-D-glucosamine (alternate)
            "CHEBI:17411",   # N-acetyl-D-galactosamine
            "CHEBI:5417",    # D-glucosamine
            "CHEBI:27903",   # D-galactosamine
            # Acidic sugars
            "CHEBI:35418",  # N-acetylneuraminic acid (sialic acid)
            "CHEBI:4178",   # D-glucuronic acid
            "CHEBI:15948",  # D-glucuronate
            # Disaccharides and oligosaccharides
            "CHEBI:17716",  # lactose
            "CHEBI:17992",  # sucrose
            "CHEBI:16646",  # cellobiose
            "CHEBI:28757",  # maltose
            "CHEBI:17306",  # trehalose
        }

        # Check for sugar/carbohydrate substrates or products
        has_sugar_substrate = any(
            p.chebi_id in sugar_indicators
            for p in reaction.left_participants
        )

        has_sugar_product = any(
            p.chebi_id in sugar_indicators
            for p in reaction.right_participants
        )

        # SMARTS-based glycoside detection for substrates
        has_glycoside_substrate = any(
            p.is_glycoside() for p in reaction.left_participants
        )

        if not (has_sugar_substrate or has_sugar_product or has_glycoside_substrate):
            return ClassificationResult(
                is_member=False,
                explanation="No carbohydrate/sugar substrates or products detected"
            )

        # Apply exclusions before classifying as glycosidase
        
        # Exclude phosphorylase reactions (use phosphate instead of water)
        has_phosphate = any(
            p.chebi_id == CHEBI_PHOSPHATE
            for p in reaction.left_participants
        )
        
        if has_phosphate:
            return ClassificationResult(
                is_member=False,
                explanation="Uses phosphate - phosphorylase not glycosidase"
            )

        # Exclude transferase reactions - detect by nucleotide products
        has_nucleotide_product = any(
            p.chebi_id in {CHEBI_UDP, CHEBI_GDP, CHEBI_ADP}
            for p in reaction.right_participants
        )

        if has_nucleotide_product:
            return ClassificationResult(
                is_member=False,
                explanation="Nucleotide product - likely glycosyltransferase not glycosidase"
            )

        # Exclude sulfatase reactions (sulfate product from sulfate ester, not glycosidic bond)
        CHEBI_SULFATE = "CHEBI:16189"
        has_sulfate_product = any(
            p.chebi_id == CHEBI_SULFATE
            for p in reaction.right_participants
        )
        if has_sulfate_product:
            return ClassificationResult(
                is_member=False,
                explanation="Sulfate product - sulfatase not glycosidase"
            )

        # Exclude phosphatase reactions (phosphate product from phosphate ester)
        # RHEA uses multiple ChEBI IDs for phosphate forms
        PHOSPHATE_PRODUCTS = {
            CHEBI_PHOSPHATE,  # CHEBI:43474
            "CHEBI:16838",    # polyphosphate (used in RHEA)
            "CHEBI:18367",    # phosphate (alternate)
        }
        has_phosphate_product = any(
            p.chebi_id in PHOSPHATE_PRODUCTS
            for p in reaction.right_participants
        ) or "phosphate" in product_str
        if has_phosphate_product:
            return ClassificationResult(
                is_member=False,
                explanation="Phosphate product - phosphatase not glycosidase"
            )

        # Exclude esterase reactions (fatty acid product)
        FATTY_ACID_INDICATORS = {
            "CHEBI:35366",  # fatty acid
            "CHEBI:28868",  # fatty acid anion
        }
        has_fatty_acid_product = any(
            p.chebi_id in FATTY_ACID_INDICATORS
            for p in reaction.right_participants
        ) or "fatty acid" in product_str
        if has_fatty_acid_product:
            return ClassificationResult(
                is_member=False,
                explanation="Fatty acid product - esterase not glycosidase"
            )

        # Exclude deacetylase reactions (acetate product from N-acetyl group, not glycosidic bond)
        CHEBI_ACETATE = "CHEBI:30089"  # acetate
        has_acetate_product = any(
            p.chebi_id == CHEBI_ACETATE
            for p in reaction.right_participants
        ) or "acetate" in product_str
        if has_acetate_product:
            return ClassificationResult(
                is_member=False,
                explanation="Acetate product - deacetylase not glycosidase"
            )

        # Exclude lyases that release pyruvate (isochorismate lyase, etc.)
        CHEBI_PYRUVATE = "CHEBI:15361"
        has_pyruvate_product = any(
            p.chebi_id == CHEBI_PYRUVATE
            for p in reaction.right_participants
        ) or "pyruvate" in product_str
        if has_pyruvate_product:
            return ClassificationResult(
                is_member=False,
                explanation="Pyruvate product - lyase not glycosidase"
            )

        # Exclude asparaginase/amidase reactions (aspartate product from amide bond)
        CHEBI_ASPARTATE = "CHEBI:29991"  # L-aspartate
        has_aspartate_product = any(
            p.chebi_id == CHEBI_ASPARTATE
            for p in reaction.right_participants
        ) or "aspartate" in product_str
        if has_aspartate_product:
            return ClassificationResult(
                is_member=False,
                explanation="Aspartate product - asparaginase not glycosidase"
            )

        # Exclude amidase reactions (ammonia/ammonium product from amide bond)
        AMMONIUM_CHEBIS = {"CHEBI:28938", "CHEBI:16134"}  # NH4+, NH3
        has_ammonia_product = any(
            p.chebi_id in AMMONIUM_CHEBIS
            for p in reaction.right_participants
        )
        if has_ammonia_product:
            return ClassificationResult(
                is_member=False,
                explanation="Ammonia product - amidase not glycosidase"
            )

        # Exclude aminoacyl-tRNA hydrolases (act on ester bond, not glycosidic)
        has_trna_substrate = "trna" in substrate_str or "t-rna" in substrate_str
        if has_trna_substrate:
            return ClassificationResult(
                is_member=False,
                explanation="tRNA substrate - aminoacyl-tRNA hydrolase not glycosidase"
            )

        # Exclude cyclic phosphodiesterases (act on P-O-C bond, not glycosidic)
        cyclic_nucleotide_patterns = [
            "cyclic gmp", "cyclic amp", "cgmp", "camp",
            "3',5'-cyclic", "2',3'-cyclic",
        ]
        has_cyclic_nucleotide = any(pat in substrate_str for pat in cyclic_nucleotide_patterns)
        if has_cyclic_nucleotide:
            return ClassificationResult(
                is_member=False,
                explanation="Cyclic nucleotide substrate - cyclic phosphodiesterase not glycosidase"
            )

        # Exclude urea product (amidase/urease type reactions)
        UREA_CHEBIS = {"CHEBI:16199", "CHEBI:57845"}  # urea
        has_urea_product = any(
            p.chebi_id in UREA_CHEBIS
            for p in reaction.right_participants
        ) or "urea" in product_str
        if has_urea_product:
            return ClassificationResult(
                is_member=False,
                explanation="Urea product - amidase/urease not glycosidase"
            )

        # Note: peptidoglycan synthesis exclusion removed - would need ChEBI-based detection

        # Try to match patterns
        match = match_patterns(reaction, self.PATTERNS, strict=False)

        explanation = "Glycosidase: glycosidic bond hydrolysis"

        if match.matched:
            explanation += " [pattern matched]"

        return ClassificationResult(is_member=True, explanation=explanation)
