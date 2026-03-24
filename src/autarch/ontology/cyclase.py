"""Cyclase reaction classification using pattern DSL.

Cyclases are lyases that form rings by eliminating small molecules.
EC 4.6.x.x classification.

History
-------

## 2025-12-22

Fixed false positives by:
1. Requiring specific cyclic nucleotide products (cAMP, cGMP, cCMP, cUMP)
   OR "cycl" in product name
2. Removed diphosphate-only matching - many lyases produce PPi without being cyclases
3. Terpene synthases are lyases (GO:0016829) but NOT cyclases (GO:0009975)
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import ReactionClass
from autarch.pattern_dsl import var, optional
from autarch.molecules import (
    CHEBI_ATP,
    CHEBI_CTP,
    CHEBI_DIPHOSPHATE,
    CHEBI_GTP,
    diphosphate,
    h_plus,
    p,
    water,
)


class Cyclase(ReactionClass):
    """cyclase

    Examples:
    - Adenylyl cyclase: ATP → cAMP + diphosphate
    - Guanylyl cyclase: GTP → cGMP + diphosphate
    - Farnesyl diphosphate synthase: precursors → farnesyl-PP + PP
    - Isopentenyl diphosphate isomerase: linear → cyclic + PP
    - Cyclase reactions forming lactones, cyclic peptides, etc.
    """

    GO_ID = "GO:0009975"  # cyclase activity
    EC_NUMBER_PREFIX = "4.6.-.-"  # Cyclases

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # Nucleotide cyclase: NTP → cyclic NMP + diphosphate
        var("nucleotide_triphosphate") >> var("cyclic_nucleotide") + p(diphosphate),
        # General cyclization: linear substrate → cyclic product + small molecule
        var("linear_substrate") >> var("cyclic_product") + var("small_molecule"),
        # With water elimination: substrate → cyclic product + H2O
        var("substrate") >> var("cyclic_product") + p(water) + optional(h_plus),
        # Terpene cyclase: linear terpene → cyclic terpene + PP
        var("linear_terpene") >> var("cyclic_terpene") + p(diphosphate),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a cyclase using pattern matching.

        Strategy:
        1. Must produce cyclic nucleotide products (cAMP, cGMP, etc.)
           OR have "cyclic" in product name
        2. NOT just any reaction producing diphosphate (terpene synthases
           produce PPi but aren't cyclases in GO)
        """
        # Specific cyclic nucleotide products (the definitive cyclase signature)
        CYCLIC_NUCLEOTIDE_PRODUCTS = {
            "CHEBI:17489",  # cAMP
            "CHEBI:16356",  # cGMP
            "CHEBI:17929",  # cCMP
            "CHEBI:16695",  # cUMP
            "CHEBI:58165",  # 3',5'-cyclic AMP(1-)
            "CHEBI:58115",  # 3',5'-cyclic GMP(1-)
        }

        # Check for cyclic nucleotide products by ChEBI ID
        has_cyclic_product = any(
            p.chebi_id in CYCLIC_NUCLEOTIDE_PRODUCTS
            for p in reaction.right_participants
        )

        if has_cyclic_product:
            return ClassificationResult(
                is_member=True,
                explanation="Cyclase: produces cyclic nucleotide (cAMP/cGMP/cCMP/cUMP)"
            )

        # Check for specific cyclic patterns in label (product side)
        # Must be "cyclic AMP", "cyclic GMP", "cyclic phosphate" etc.
        # NOT just any product with "cyclic" anywhere in name
        label_lower = reaction.label.lower() if reaction.label else ""
        product_str = label_lower.split("=")[1] if "=" in label_lower else ""

        def has_cyclic_nucleotide_pattern(text: str) -> bool:
            """Check if text contains cyclic nucleotide pattern."""
            # Match "cyclic amp", "cyclic gmp", "cyclic cmp", "cyclic ump"
            if any(f"cyclic {x}" in text for x in ["amp", "gmp", "cmp", "ump"]):
                return True
            # Match "3',5'-cyclic" pattern (nucleotide cyclases)
            if "3',5'-cyclic" in text:
                return True
            return False

        has_cyclic_product_label = has_cyclic_nucleotide_pattern(product_str)

        if has_cyclic_product_label:
            return ClassificationResult(
                is_member=True,
                explanation="Cyclase: produces cyclic nucleotide product"
            )

        # Check for "cyclase" pattern: NTP → cyclic product + PPi
        # Must have NTP substrate AND diphosphate product
        cyclase_substrates = {
            CHEBI_ATP,  # ATP (→ cAMP)
            CHEBI_GTP,  # GTP (→ cGMP)
            CHEBI_CTP,  # CTP (→ cCMP)
            "CHEBI:37568",  # UTP (→ cUMP)
        }

        has_ntp_substrate = any(
            p.chebi_id in cyclase_substrates
            for p in reaction.left_participants
        )

        has_diphosphate = any(
            p.chebi_id == CHEBI_DIPHOSPHATE
            for p in reaction.right_participants
        )

        # For NTP → product + PPi, check if it's likely cyclase vs other reactions
        if has_ntp_substrate and has_diphosphate:
            # Cyclases have 1→2 stoichiometry (NTP → cyclic + PPi)
            # NOT 3→3 like ligases (substrate + substrate + NTP → product + AMP + PPi)
            if len(reaction.left_participants) == 1 and len(reaction.right_participants) in [2, 3]:
                return ClassificationResult(
                    is_member=True,
                    explanation="Cyclase: NTP → product + diphosphate (1→2/3 stoichiometry)"
                )

        # NOTE: We do NOT match just on diphosphate output - terpene synthases
        # produce PPi but are lyases (GO:0016829), not cyclases (GO:0009975)

        return ClassificationResult(
            is_member=False,
            explanation="No cyclic nucleotide product or cyclase pattern detected"
        )
