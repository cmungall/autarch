"""Deaminase (cyclic amidine hydrolase) reaction classification.

Deaminases (EC 3.5.4) are hydrolases acting on cyclic amidines.
They catalyze the hydrolytic removal of amino groups, releasing ammonia.

Pattern: amino-compound + H2O → hydroxyl-compound + NH4+
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.hydrolase import Hydrolase
from autarch.molecules import (
    CHEBI_H2O,
)

# Ammonia/ammonium (product of deamination)
CHEBI_AMMONIA = "CHEBI:16134"
CHEBI_AMMONIUM = "CHEBI:28938"

# Common deaminase substrates (nucleosides/nucleotides)
CHEBI_ADENOSINE = "CHEBI:16335"
CHEBI_ADENINE = "CHEBI:16708"
CHEBI_CYTIDINE = "CHEBI:17562"
CHEBI_CYTOSINE = "CHEBI:16040"
CHEBI_GUANINE = "CHEBI:16235"
CHEBI_AMP = "CHEBI:16027"
CHEBI_CMP = "CHEBI:17361"

# Deamination products
CHEBI_INOSINE = "CHEBI:17596"
CHEBI_HYPOXANTHINE = "CHEBI:17368"
CHEBI_URIDINE = "CHEBI:16704"
CHEBI_URACIL = "CHEBI:17568"
CHEBI_XANTHINE = "CHEBI:17712"


class Deaminase(Hydrolase):
    """deaminase

    EC 3.5.4.x includes:
    - Adenosine deaminase: adenosine + H2O → inosine + NH4+
    - Cytidine deaminase: cytidine + H2O → uridine + NH4+
    - Guanine deaminase: guanine + H2O → xanthine + NH4+
    - AMP deaminase: AMP + H2O → IMP + NH4+

    Examples:
    - Adenosine + H2O → inosine + NH4+
    - Cytosine + H2O → uracil + NH4+
    """

    GO_ID = "GO:0019239"  # deaminase activity
    EC_NUMBER_PREFIX = "3.5.4.-"  # Acting on cyclic amidines

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a deaminase.

        Strategy:
        1. Must use water as substrate
        2. Must produce ammonia/ammonium
        3. Look for characteristic deamination substrates/products
        """
        # Must use water
        has_water = any(
            p.chebi_id == CHEBI_H2O
            for p in reaction.left_participants
        )

        if not has_water:
            return ClassificationResult(
                is_member=False,
                explanation="No water substrate - deaminases are hydrolases"
            )

        # Use RHEA label for pattern matching (not p.name which has ChEBI generic names)
        label_lower = reaction.label.lower() if reaction.label else ""
        label_parts = label_lower.split("=")
        substrate_str = label_parts[0] if label_parts else ""
        product_str = label_parts[1] if len(label_parts) > 1 else ""

        # Must produce ammonia or ammonium - ChEBI ID primary, label fallback
        has_ammonia = any(
            p.chebi_id in {CHEBI_AMMONIA, CHEBI_AMMONIUM}
            for p in reaction.right_participants
        ) or any(x in product_str for x in ["ammonia", "ammonium", "nh4", "nh3"])

        if not has_ammonia:
            return ClassificationResult(
                is_member=False,
                explanation="No ammonia/ammonium product - not deamination"
            )

        # Check for characteristic deaminase substrates (ChEBI IDs)
        deaminase_substrates = {
            CHEBI_ADENOSINE: "adenosine",
            CHEBI_ADENINE: "adenine",
            CHEBI_CYTIDINE: "cytidine",
            CHEBI_CYTOSINE: "cytosine",
            CHEBI_GUANINE: "guanine",
            CHEBI_AMP: "AMP",
            CHEBI_CMP: "CMP",
        }

        # Check for characteristic deaminase products (ChEBI IDs)
        deaminase_products = {
            CHEBI_INOSINE: "inosine",
            CHEBI_HYPOXANTHINE: "hypoxanthine",
            CHEBI_URIDINE: "uridine",
            CHEBI_URACIL: "uracil",
            CHEBI_XANTHINE: "xanthine",
        }

        # Check substrates by ChEBI ID first, then label pattern
        substrate_name = None
        for p in reaction.left_participants:
            if p.chebi_id in deaminase_substrates:
                substrate_name = deaminase_substrates[p.chebi_id]
                break

        if not substrate_name:
            substrate_patterns = ["adenosine", "adenine", "cytidine", "cytosine", "guanine"]
            for pat in substrate_patterns:
                if pat in substrate_str:
                    substrate_name = pat
                    break

        # Check products by ChEBI ID first, then label pattern
        product_name = None
        for p in reaction.right_participants:
            if p.chebi_id in deaminase_products:
                product_name = deaminase_products[p.chebi_id]
                break

        if not product_name:
            product_patterns = ["inosine", "hypoxanthine", "uridine", "uracil", "xanthine"]
            for pat in product_patterns:
                if pat in product_str:
                    product_name = pat
                    break

        # Exclude amino acid oxidases/dehydrogenases (EC 1.4 - use O2 or NAD/NADP)
        has_o2 = any(x in substrate_str for x in ["dioxygen", " o2", "oxygen"])
        has_nad = "nad" in label_lower

        if has_o2 or has_nad:
            return ClassificationResult(
                is_member=False,
                explanation="Amino acid oxidase/dehydrogenase - not cyclic amidine deaminase"
            )

        # Exclude simple amidases (EC 3.5.1) - linear amide hydrolysis
        simple_amide_patterns = ["amide", "amidomethylene", "carbamoyl", "carbamide"]
        has_simple_amide = any(pat in substrate_str for pat in simple_amide_patterns)

        if has_simple_amide and not substrate_name:
            return ClassificationResult(
                is_member=False,
                explanation="Simple amidase (EC 3.5.1) - not cyclic amidine deaminase"
            )

        # Exclude reactions producing CO2 (carbamoyl hydrolases)
        has_co2 = "carbon dioxide" in product_str or " co2" in product_str

        if has_co2:
            return ClassificationResult(
                is_member=False,
                explanation="Carbamoyl hydrolase - produces CO2, not cyclic amidine deaminase"
            )

        # Must have nucleoside/nucleotide pattern to be a deaminase (EC 3.5.4)
        if substrate_name and product_name:
            return ClassificationResult(
                is_member=True,
                explanation=f"Deaminase: {substrate_name} → {product_name} + NH4+"
            )

        if substrate_name:
            return ClassificationResult(
                is_member=True,
                explanation=f"Deaminase: {substrate_name} deamination"
            )

        # Without specific nucleoside patterns, reject
        return ClassificationResult(
            is_member=False,
            explanation="No nucleoside/nucleotide deamination pattern"
        )
