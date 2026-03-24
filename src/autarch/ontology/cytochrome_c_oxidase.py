"""Cytochrome C Oxidase reaction classification using pattern DSL.

Cytochrome C Oxidase (Complex IV) is the terminal enzyme of the electron transport
chain, catalyzing the four-electron reduction of O2 to water using cytochrome c
as the electron donor. This is THE most specific respiratory enzyme.

Key biochemical signature: 4 cytochrome c(Fe2+) + O2 + 4H+ → 4 cytochrome c(Fe3+) + 2H2O

Only one major reaction in all of biochemistry - the terminal step of respiration.
Tests our ability to classify the most specific biochemical transformations.

History
-------

## 2025-12-22

Fixed false positives by:
1. Requiring cytochrome c-specific ChEBI IDs (was matching any O2→H2O reaction)
2. Excluding 2-oxoglutarate-dependent dioxygenases
3. Excluding CO2-producing reactions (oxygenases, not oxidases)
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.primary_active_transmembrane_transporter import (
    PrimaryActiveTransmembraneTransporter,
)
from autarch.pattern_dsl import var
from autarch.molecules import (
    CHEBI_CO2,
    CHEBI_H2O,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_O2,
    h_plus,
    oxygen,
    p,
    water,
)


class CytochromeCOxidase(PrimaryActiveTransmembraneTransporter):
    """cytochrome c oxidase

    Examples:
    - Complex IV: terminal electron transport
    - Mitochondrial respiration: final oxygen reduction
    - Bacterial cytochrome oxidases: respiratory variants

    This is THE ULTIMATE SPECIFICITY:
    - Only 1-2 reactions in entire biochemistry
    - Requires cytochrome c as electron donor
    - Four-electron reduction of O2
    - Essential for aerobic life
    """

    GO_ID = "GO:0004129"  # cytochrome-c oxidase activity
    EC_NUMBER_PREFIX = "7.1.1.9"  # cytochrome-c oxidase activity

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # Classic cytochrome oxidase: 4 cyt c + O2 + 4H+ → 4 cyt c+ + 2H2O
        var("cytochrome_c_reduced") + p(oxygen) + p(h_plus)
        >> var("cytochrome_c_oxidized") + p(water),
        # Without explicit H+
        var("cytochrome_c_reduced") + p(oxygen)
        >> var("cytochrome_c_oxidized") + p(water),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is cytochrome c oxidase using ULTRA-SPECIFIC criteria.

        CYTOCHROME C OXIDASE SIGNATURE:
        1. Must consume O2 (molecular oxygen)
        2. Must produce H2O (four-electron reduction)
        3. Must involve cytochrome c (heme protein electron donor)
        4. Should NOT produce CO2 (that indicates oxygenase)
        5. Should NOT involve NAD(P)H or 2-oxoglutarate
        """
        # Cytochrome c ChEBI IDs (specific electron donors)
        CYTOCHROME_C_CHEBIS = {
            "CHEBI:29034",  # ferrocytochrome c (Fe2+)
            "CHEBI:29033",  # ferricytochrome c (Fe3+)
            "CHEBI:4067",   # cytochrome c
            "CHEBI:18063",  # cytochrome c2
        }

        # Check for cytochrome c involvement
        has_cytochrome_c = any(
            p.chebi_id in CYTOCHROME_C_CHEBIS
            for p in reaction.left_participants + reaction.right_participants
        )

        if not has_cytochrome_c:
            return ClassificationResult(
                is_member=False,
                explanation="No cytochrome c detected - CcO requires cytochrome c as electron donor"
            )

        # Must consume O2 as terminal electron acceptor
        has_oxygen = any(
            p.chebi_id == CHEBI_O2 for p in reaction.left_participants
        )

        if not has_oxygen:
            return ClassificationResult(
                is_member=False,
                explanation="No O2 reactant - cytochrome oxidase requires molecular oxygen"
            )

        # Must produce H2O (four-electron reduction of oxygen)
        has_water_product = any(
            p.chebi_id == CHEBI_H2O for p in reaction.right_participants
        )

        if not has_water_product:
            return ClassificationResult(
                is_member=False,
                explanation="No H2O product - terminal oxidases produce water"
            )

        # Should NOT produce CO2 (that indicates oxygenase/decarboxylase)
        has_co2_product = any(
            p.chebi_id == CHEBI_CO2 for p in reaction.right_participants
        )

        if has_co2_product:
            return ClassificationResult(
                is_member=False,
                explanation="Produces CO2 - oxygenase not cytochrome oxidase"
            )

        # Should NOT involve NAD(P)H (that's upstream in electron transport)
        has_nadh_cofactor = any(
            p.chebi_id in [CHEBI_NADPH, CHEBI_NADH]
            for p in reaction.left_participants
        )

        if has_nadh_cofactor:
            return ClassificationResult(
                is_member=False,
                explanation="Uses NAD(P)H - likely Complex I/III, not terminal oxidase"
            )

        return ClassificationResult(
            is_member=True,
            explanation="Cytochrome C Oxidase: terminal respiratory enzyme (cyt c + O2 → H2O)"
        )
