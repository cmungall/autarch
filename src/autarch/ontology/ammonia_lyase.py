"""Ammonia lyase reaction classification.

EC 4.3.1: Ammonia-lyases catalyze the release of ammonia from substrates,
often forming a double bond (C=C or C=N).

Pattern: R-CH(-NH2)-R' = R-CH=R' + NH3
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.lyase import Lyase
from autarch.molecules import (
    CHEBI_H2O,
    CHEBI_NH3,
    CHEBI_NH4,
    CHEBI_NAD_PLUS,
    CHEBI_NADH,
    CHEBI_NADP_PLUS,
    CHEBI_NADPH,
    CHEBI_O2,
    CHEBI_H2O2,
)


class AmmoniaLyase(Lyase):
    """ammonia lyase

    EC 4.3.1.x includes:
    - Aspartate ammonia-lyase: L-aspartate → fumarate + NH3
    - Histidine ammonia-lyase: L-histidine → urocanate + NH3
    - Phenylalanine ammonia-lyase: L-phenylalanine → trans-cinnamate + NH3
    - Serine ammonia-lyase: L-serine → pyruvate + NH3
    - Threonine ammonia-lyase: L-threonine → 2-oxobutanoate + NH3

    Examples:
    - L-aspartate → fumarate + NH3
    - L-phenylalanine → trans-cinnamate + NH3
    """

    GO_ID = "GO:0016841"  # ammonia-lyase activity
    EC_NUMBER_PREFIX = "4.3.1.-"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is an ammonia lyase.

        Strategy:
        1. Must release ammonia (NH3 or NH4+)
        2. Must NOT use water as reactant (would be hydrolase)
        3. Must NOT use O2 or redox cofactors (would be oxidoreductase)
        4. Look for characteristic amino acid substrates/products
        """
        # Use RHEA label for pattern matching (not p.name which has ChEBI generic names)
        label_lower = reaction.label.lower() if reaction.label else ""
        label_parts = label_lower.split("=")
        substrate_str = label_parts[0] if label_parts else ""
        product_str = label_parts[1] if len(label_parts) > 1 else ""

        # Must release ammonia or ammonium
        has_ammonia_product = any(
            p.chebi_id in {CHEBI_NH3, CHEBI_NH4}
            for p in reaction.right_participants
        ) or any(x in product_str for x in ["ammonia", "ammonium", "nh3", "nh4"])

        if not has_ammonia_product:
            return ClassificationResult(
                is_member=False,
                explanation="No ammonia release"
            )

        # Must NOT use water as reactant (that would be deaminase/hydrolase)
        has_water_reactant = any(
            p.chebi_id == CHEBI_H2O
            for p in reaction.left_participants
        )

        if has_water_reactant:
            return ClassificationResult(
                is_member=False,
                explanation="Uses water - hydrolase/deaminase, not ammonia lyase"
            )

        # Ammonia lyases typically: 1 substrate → multiple products
        # NOT 2 → 2 which is often a transferase releasing ammonia
        left_count = len(reaction.left_participants)
        right_count = len(reaction.right_participants)

        # Exclude 2 → 2 reactions (transferase pattern)
        if left_count == 2 and right_count == 2:
            return ClassificationResult(
                is_member=False,
                explanation="2→2 stoichiometry - likely transferase, not lyase"
            )

        # Must NOT use redox cofactors (that would be amino acid oxidoreductase)
        has_redox = any(
            p.chebi_id in {CHEBI_NAD_PLUS, CHEBI_NADH, CHEBI_NADP_PLUS, CHEBI_NADPH,
                          CHEBI_O2, CHEBI_H2O2}
            for p in reaction.left_participants + reaction.right_participants
        )

        if has_redox:
            return ClassificationResult(
                is_member=False,
                explanation="Uses redox cofactors - oxidoreductase, not lyase"
            )

        # Characteristic amino acid substrates
        amino_acid_patterns = [
            "aspartate", "aspartic",
            "histidine", "phenylalanine", "tyrosine",
            "serine", "threonine", "tryptophan",
            "alanine", "glutamate", "glutamic",
            "arginine", "ornithine", "lysine",
            "cysteine", "methionine",
            "aminobutyr", "aminopentan",
        ]

        # Characteristic unsaturated products
        unsaturated_products = [
            "fumarate", "cinnamate", "coumarate", "urocanate",
            "crotonoyl", "acryloyl", "mesaconate",
            "butenoyl", "pentenoyl",
            "pyruvate", "oxobutano", "oxopentan",  # keto-products from Schiff base rearrangement
        ]

        has_amino_substrate = any(aa in substrate_str for aa in amino_acid_patterns)

        has_unsaturated_product = any(up in product_str for up in unsaturated_products)

        # CoA-linked substrates (3-aminobutyryl-CoA, beta-alanyl-CoA)
        coa_patterns = ["aminobutyryl-coa", "aminobutanoyl-coa", "alanyl-coa", "aminoacyl-coa"]
        has_coa_substrate = any(cp in substrate_str for cp in coa_patterns)

        if has_amino_substrate and has_unsaturated_product:
            # Derive substrate name from label
            substrate_name = next(
                (aa for aa in amino_acid_patterns if aa in substrate_str),
                "amino acid"
            )

            return ClassificationResult(
                is_member=True,
                explanation=f"Ammonia-lyase: {substrate_name} → unsaturated + NH3"
            )

        if has_coa_substrate:
            return ClassificationResult(
                is_member=True,
                explanation="Ammonia-lyase: aminoacyl-CoA → enoyl-CoA + NH3"
            )

        # Cyclodeaminase pattern (e.g., ornithine cyclodeaminase: ornithine → proline + NH3)
        cyclic_products = ["proline", "pyrrolidine", "piperidine", "lactam"]
        has_cyclic_product = any(cp in product_str for cp in cyclic_products)

        if has_amino_substrate and has_cyclic_product:
            return ClassificationResult(
                is_member=True,
                explanation="Ammonia-lyase: cyclodeaminase (ring formation + NH3)"
            )

        # Generic pattern: amino-containing substrate releases ammonia
        if has_amino_substrate:
            return ClassificationResult(
                is_member=True,
                explanation="Ammonia-lyase: amino acid deamination"
            )

        return ClassificationResult(
            is_member=False,
            explanation="No ammonia-lyase pattern"
        )
