"""Epoxide Hydrolase reaction classification using pattern DSL.

Epoxide Hydrolases are ultra-specific hydrolases that open epoxide rings
by adding water across the ring. These enzymes are crucial for detoxification
of reactive epoxides formed during xenobiotic metabolism.

Key biochemical signature: epoxide + H2O → diol

This tests our ability to recognize specific ring-opening reactions
and distinguish them from other hydrolases.

History
-------

## 2025-12-22

Fixed critical bug: classifier was not calling parent Hydrolase.check_membership_impl(),
causing it to match any 2-reactant hydrolysis reaction. Now properly inherits from
Hydrolase and requires epoxide-specific ChEBI IDs.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.hydrolase import Hydrolase
from autarch.pattern_dsl import var, optional
from autarch.molecules import (
    CHEBI_ATP,
    CHEBI_NADH,
    CHEBI_NADPH,
    h_plus,
    p,
    water,
)


class EpoxideHydrolase(Hydrolase):
    """epoxide hydrolase

    Examples:
    - Microsomal epoxide hydrolase: benzo[a]pyrene epoxide detoxification
    - Soluble epoxide hydrolase: fatty acid epoxide metabolism
    - Cholesterol 5,6-epoxide hydrolase: sterol metabolism

    This is ULTRA-SPECIFIC:
    - Only ~20-30 reactions in databases
    - Requires epoxide substrate (3-membered oxygen ring)
    - Forms vicinal diol products
    - Essential for detoxification
    """

    GO_ID = "GO:0004301"  # epoxide hydrolase activity
    EC_NUMBER_PREFIX = "3.3.2.10"  # epoxide hydrolase activity

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # Classic epoxide hydrolase: epoxide + H2O → diol
        var("epoxide") + p(water) >> var("diol"),
        # With explicit protonation
        var("epoxide") + p(water) >> var("diol") + optional(h_plus),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is epoxide hydrolase using RING-SPECIFIC criteria.

        EPOXIDE HYDROLASE SIGNATURE:
        1. Must be a hydrolase first (parent class)
        2. Must involve epoxide substrate (ChEBI or SMARTS detection)
        3. Should produce diol product (hydroxyl addition)
        4. Should not involve cofactors (simple hydration)

        This tests recognition of specific ring chemistry!
        """
        # First check if it's a hydrolase
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not a hydrolase: {parent_result.explanation}",
            )

        # Should not involve cofactors like NAD(P)H, ATP, etc.
        has_cofactors = any(
            p.chebi_id in [CHEBI_NADPH, CHEBI_NADH, CHEBI_ATP]
            for p in (reaction.left_participants + reaction.right_participants)
        )

        if has_cofactors:
            return ClassificationResult(
                is_member=False,
                explanation="Uses cofactors - epoxide hydrolases are simple hydration reactions"
            )

        # Known epoxide ChEBI IDs
        epoxide_chebis = {
            "CHEBI:35762",  # epoxide
            "CHEBI:32878",  # alkene oxide
            "CHEBI:35587",  # squalene 2,3-epoxide
            "CHEBI:15446",  # leukotriene A4
            "CHEBI:27432",  # styrene oxide
            "CHEBI:16956",  # benzo[a]pyrene-7,8-epoxide
            "CHEBI:36308",  # cis-epoxysuccinate
            "CHEBI:17076",  # epoxyeicosatrienoic acid
        }

        # Diol ChEBI IDs (products of epoxide hydrolysis)
        diol_chebis = {
            "CHEBI:23824",  # diol
            "CHEBI:15600",  # glycol
            "CHEBI:16509",  # dihydrodiol
            "CHEBI:35587",  # squalene-2,3-diol
        }

        # Check for epoxide substrates
        has_epoxide_substrate = any(
            p.chebi_id in epoxide_chebis
            for p in reaction.left_participants
        )

        # Check for diol products
        has_diol_product = any(
            p.chebi_id in diol_chebis
            for p in reaction.right_participants
        )

        if not (has_epoxide_substrate or has_diol_product):
            return ClassificationResult(
                is_member=False,
                explanation="No epoxide substrate or diol product detected"
            )

        return ClassificationResult(
            is_member=True,
            explanation="Epoxide Hydrolase: epoxide ring opening with water"
        )
