"""Heme-dependent oxidase reaction classification.

GO:0016675 is specifically "oxidoreductase activity, acting on a heme group
of donors" - this covers cytochrome oxidases and similar heme-dependent reactions.

History
-------

## 2025-12-22 (v2)

Major fix: GO:0016675 is a SPECIFIC term for heme-dependent oxidoreductases,
not general oxidases. The classifier was matching any O2→H2O/H2O2 reaction
(165 FP, 0 TP). Fixed by:
1. Requiring heme/cytochrome indicators in participant names or ChEBI IDs
2. This is a specialized classifier, not a general oxidase detector

## 2025-12-22 (v1)

Fixed false positives by excluding oxygenases. Key distinction:
- Oxidases: O2 is terminal electron acceptor → H2O or H2O2 (no NAD(P)H needed)
- Oxygenases: O2 atoms are incorporated into substrate (use NAD(P)H as reductant)
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.oxidoreductase import Oxidoreductase
from autarch.pattern_dsl import var, optional
from autarch.molecules import (
    oxygen, water, hydrogen_peroxide, h_plus, p,
    CHEBI_O2, CHEBI_H2O,
)


class OxidoreductaseActingOnAHemeGroupOfDonors(Oxidoreductase):
    """oxidoreductase acting on a heme group of donors

    Examples:
    - Cytochrome c oxidase: cytochrome c + O2 → oxidized cytochrome c + H2O
    - Glucose oxidase: glucose + O2 → gluconate + H2O2
    - Monoamine oxidase: amine + O2 → aldehyde + NH3 + H2O2
    """

    GO_ID = "GO:0016675"  # oxidoreductase activity, acting on a heme group of donors
    EC_NUMBER_PREFIX = "1.9.-.-"  # oxygen as acceptor

    # Define patterns using DSL with operator overloading
    # Type as list[Reaction] for mypy compatibility
    PATTERNS: list[Reaction] = [
        # Terminal oxidation: substrate + O2 → product + H2O + H+?
        var("substrate") + p(oxygen) >> var("product") + p(water) + optional(h_plus),
        # Peroxide formation: substrate + O2 → product + H2O2
        var("substrate") + p(oxygen) >> var("product") + p(hydrogen_peroxide),
        # Complex oxidation: substrate1 + substrate2 + O2 → product1 + product2 + H2O + H+?
        var("substrate1") + var("substrate2") + p(oxygen)
        >> var("product1") + var("product2") + p(water) + optional(h_plus),
        # Multiple substrate oxidation with H2O2: substrate1 + substrate2 + O2 → product1 + product2 + H2O2
        var("substrate1") + var("substrate2") + p(oxygen)
        >> var("product1") + var("product2") + p(hydrogen_peroxide),
        # Simple oxidation without H+ production: substrate + O2 → product + H2O
        var("substrate") + p(oxygen) >> var("product") + p(water),
        # General O2-dependent reaction: substrate + O2 → product + byproduct
        var("substrate") + p(oxygen) >> var("product") + var("byproduct"),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a heme-dependent oxidase.

        GO:0016675 is "oxidoreductase activity, acting on a heme group of donors"
        This is a SPECIFIC class for cytochrome oxidases and similar enzymes.

        Strategy:
        1. MUST involve heme/cytochrome (by ChEBI ID or label pattern)
        2. Must be an oxidoreductase (parent class)
        3. O2 involvement is optional for this GO term (some are non-O2)
        """
        # Use RHEA label for pattern matching
        label_lower = reaction.label.lower() if reaction.label else ""

        # CRITICAL: This GO term specifically requires heme/cytochrome involvement
        # Without heme indicators, we cannot classify as GO:0016675

        # Heme/cytochrome ChEBI IDs
        # NOTE: CHEBI:29033/29034 are iron ions, NOT cytochromes (data quality issue)
        # The cytochrome context is in RHEA labels but lost in participant data
        HEME_CHEBI_IDS = {
            "CHEBI:17627",   # heme (generic)
            "CHEBI:60344",   # heme b
            "CHEBI:30413",   # heme c
            "CHEBI:18070",   # cytochrome c (protein)
            "CHEBI:18367",   # cytochrome (generic)
            # NOT CHEBI:29033/29034 - those are just Fe2+/Fe3+ ions
        }

        # Check for heme/cytochrome by ChEBI ID
        has_heme_chebi = any(
            p.chebi_id in HEME_CHEBI_IDS
            for p in reaction.left_participants + reaction.right_participants
        )

        # Check for heme/cytochrome in label (not p.name)
        heme_patterns = ["cytochrome", "heme", "haem", "ferricytochrome", "ferrocytochrome"]
        has_heme_label = any(pat in label_lower for pat in heme_patterns)

        if not (has_heme_chebi or has_heme_label):
            return ClassificationResult(
                is_member=False,
                explanation="No heme/cytochrome involvement - GO:0016675 requires heme group donors"
            )

        # Check if it's an oxidoreductase (parent class)
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not an oxidoreductase: {parent_result.explanation}",
            )

        # At this point we have heme involvement and oxidoreductase chemistry
        # Build explanation based on specific indicators
        explanation = "Heme-dependent oxidoreductase"

        # Check for O2 involvement (common in cytochrome oxidases)
        has_o2 = any(
            p.chebi_id == CHEBI_O2
            for p in reaction.left_participants
        )

        if has_o2:
            produces_water = any(
                p.chebi_id == CHEBI_H2O for p in reaction.right_participants
            )
            if produces_water:
                explanation = "Cytochrome oxidase: heme + O2 → H2O"
            else:
                explanation = "Heme-dependent oxidase with O2"
        else:
            explanation = "Heme-dependent oxidoreductase (non-O2)"

        return ClassificationResult(is_member=True, explanation=explanation)
