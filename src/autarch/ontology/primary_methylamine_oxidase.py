"""Primary amine oxidase reaction classification using pattern DSL.

GO:0008131 "primary methylamine oxidase activity" includes oxidases that
catalyze oxidative deamination of primary amines (including amino acids)
producing aldehydes/ketoacids, ammonia, and hydrogen peroxide.

Key biochemical signature: R-NH2 + O2 + H2O → R=O + NH4+ + H2O2

History
-------

## 2025-12-22 (v2)

Changed parent from Oxidase (heme-dependent) to Oxidoreductase because:
1. Primary amine oxidases use FAD, not heme
2. GO:0008131 includes amino acid oxidases as children
3. Fixed substrate detection to include amino acids
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.oxidoreductase import Oxidoreductase
from autarch.pattern_dsl import var, optional
from autarch.molecules import (
    CHEBI_H2O2,
    CHEBI_NH3,
    CHEBI_O2,
    ammonia,
    oxygen,
    p,
    water,
)


class PrimaryMethylamineOxidase(Oxidoreductase):
    """primary methylamine oxidase

    GO:0008131 "primary methylamine oxidase activity" covers:
    - Monoamine oxidases (MAO-A, MAO-B): serotonin, dopamine, etc.
    - Amino acid oxidases: D/L-glutamate, D/L-aspartate oxidases
    - Amine oxidases: cyclohexylamine oxidase, etc.

    All produce H2O2 as characteristic product (flavin-dependent).
    """

    GO_ID = "GO:0008131"  # primary methylamine oxidase activity
    EC_NUMBER_PREFIX = "1.4.3.21"  # primary methylamine oxidase activity

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # Classic MAO: amine + FAD + O2 + H2O → aldehyde + NH3 + FADH2 + H2O2
        var("monoamine") + var("FAD") + p(oxygen) + p(water)
        >> var("aldehyde") + p(ammonia) + var("FADH2") + var("hydrogen_peroxide"),
        # Simplified without explicit cofactors
        var("monoamine") + p(oxygen) + p(water)
        >> var("aldehyde") + p(ammonia) + var("hydrogen_peroxide"),
        # With different cofactor states
        var("monoamine") + var("cofactor") + p(oxygen) + optional(water)
        >> var("aldehyde") + p(ammonia) + var("cofactor_reduced") + var("hydrogen_peroxide"),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a primary amine oxidase.

        PRIMARY AMINE OXIDASE SIGNATURE (GO:0008131):
        1. Must consume O2 (oxidative reaction)
        2. Must produce H2O2 (characteristic of flavin-dependent oxidases)
        3. Must produce NH4+/NH3 (from deamination)

        This covers monoamine oxidases AND amino acid oxidases.
        """
        # Use RHEA label for pattern matching
        label_lower = reaction.label.lower() if reaction.label else ""
        label_parts = label_lower.split("=")
        product_str = label_parts[1] if len(label_parts) > 1 else ""

        # Must consume O2
        has_oxygen = any(
            p.chebi_id == CHEBI_O2 for p in reaction.left_participants
        )

        if not has_oxygen:
            return ClassificationResult(
                is_member=False,
                explanation="No O2 reactant - requires molecular oxygen"
            )

        # Must produce H2O2 (THE characteristic product of flavin-dependent amine oxidases)
        has_hydrogen_peroxide = any(
            p.chebi_id == CHEBI_H2O2
            for p in reaction.right_participants
        )

        if not has_hydrogen_peroxide:
            return ClassificationResult(
                is_member=False,
                explanation="No H2O2 product - amine oxidases produce hydrogen peroxide"
            )

        # Must produce NH4+ or NH3 (from deamination)
        # Use ChEBI IDs primarily, label as fallback
        CHEBI_NH4 = "CHEBI:28938"  # ammonium
        has_ammonia_product = any(
            p.chebi_id in {CHEBI_NH3, CHEBI_NH4}
            for p in reaction.right_participants
        ) or any(x in product_str for x in ["ammonia", "ammonium", "nh4(+)", "nh3"])

        if not has_ammonia_product:
            return ClassificationResult(
                is_member=False,
                explanation="No NH4+/NH3 product - no deamination"
            )

        # The combination of O2 + H2O2 + NH4+ is highly specific for amine oxidases
        return ClassificationResult(
            is_member=True,
            explanation="Primary amine oxidase: O2 + amine → H2O2 + NH4+ + keto product"
        )
