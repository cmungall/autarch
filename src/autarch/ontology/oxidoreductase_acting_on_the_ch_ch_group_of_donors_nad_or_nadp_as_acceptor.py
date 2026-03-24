"""Enoyl reductase / CH-CH oxidoreductase reaction classification.

EC 1.3.1.x are oxidoreductases acting on the CH-CH group with NAD+/NADP+ as acceptor.
These catalyze the interconversion between saturated and unsaturated C-C bonds.

This classifier inherits from OxidoreductaseActingOnTheCHCHGroupOfDonors which handles:
- Detection of CH-CH patterns (enoyl, saturated/unsaturated)
- Exclusion of CH-OH patterns (→ Dehydrogenase EC 1.1)

Pattern: saturated-substrate + NAD(P)+ ⇌ unsaturated-substrate + NAD(P)H + H+
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.oxidoreductase_acting_on_the_ch_ch_group_of_donors import (
    OxidoreductaseActingOnTheCHCHGroupOfDonors,
)
from autarch.molecules import (
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_FAD,
    CHEBI_FADH2,
)

# ACP (acyl carrier protein) indicators
CHEBI_ACP = "CHEBI:78785"  # saturated acyl-ACP
CHEBI_ENOYL_ACP = "CHEBI:78784"  # enoyl-ACP


class OxidoreductaseActingOnTheCHCHGroupOfDonorsNADOrNADPAsAcceptor(
    OxidoreductaseActingOnTheCHCHGroupOfDonors
):
    """oxidoreductase acting on the CH-CH group of donors, NAD or NADP as acceptor

    EC 1.3.1.x includes:
    - Enoyl-ACP reductase: enoyl-ACP + NADH → acyl-ACP + NAD+
    - Dihydroorotate dehydrogenase: dihydroorotate + NAD+ → orotate + NADH
    - Butyryl-CoA dehydrogenase: butanoyl-CoA + NAD+ → crotonyl-CoA + NADH

    Key distinctions (now handled by parent OxidoreductaseActingOnTheCHCHGroupOfDonors):
    - NOT EC 1.1 (CH-OH group) → Dehydrogenase
    - NOT EC 1.3.8 (FAD-dependent) → handled here by excluding FAD
    """

    GO_ID = "GO:0016628"  # oxidoreductase activity, acting on CH-CH with NAD/NADP
    EC_NUMBER_PREFIX = "1.3.1.-"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a CH-CH oxidoreductase with NAD/NADP (EC 1.3.1).

        Strategy:
        1. Must pass OxidoreductaseActingOnTheCHCHGroupOfDonors checks (CH-CH pattern detected)
        2. Must use NAD+/NADP+ system (not FAD which is EC 1.3.8)
        3. Identify specific substrate patterns
        """
        # Check parent class (handles CH-CH pattern detection)
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        # Must use NAD+/NADP+ system, NOT FAD
        nad_system = {CHEBI_NAD_PLUS, CHEBI_NADH, CHEBI_NADP_PLUS, CHEBI_NADPH}
        fad_system = {CHEBI_FAD, CHEBI_FADH2}

        has_nad = any(
            p.chebi_id in nad_system
            for p in reaction.left_participants + reaction.right_participants
        )

        has_fad = any(
            p.chebi_id in fad_system
            for p in reaction.left_participants + reaction.right_participants
        )

        if has_fad:
            return ClassificationResult(
                is_member=False,
                explanation="FAD-dependent - EC 1.3.8, not EC 1.3.1"
            )

        if not has_nad:
            return ClassificationResult(
                is_member=False,
                explanation="No NAD/NADP system - EC 1.3.1 requires NAD/NADP"
            )

        # Check for ACP involvement (fatty acid biosynthesis) by ChEBI only.
        has_acp = any(
            p.chebi_id in {CHEBI_ACP, CHEBI_ENOYL_ACP}
            for p in reaction.left_participants + reaction.right_participants
        )

        if has_acp:
            return ClassificationResult(
                is_member=True,
                explanation="Enoyl-ACP reductase: CH-CH reduction in fatty acid biosynthesis"
            )

        # Determine direction based on cofactor position
        has_oxidized_left = any(
            p.chebi_id in [CHEBI_NAD_PLUS, CHEBI_NADP_PLUS]
            for p in reaction.left_participants
        )
        has_reduced_right = any(
            p.chebi_id in [CHEBI_NADH, CHEBI_NADPH]
            for p in reaction.right_participants
        )

        if has_oxidized_left and has_reduced_right:
            # Desaturation direction
            cofactor = (
                "NAD+"
                if any(p.chebi_id == CHEBI_NAD_PLUS for p in reaction.left_participants)
                else "NADP+"
            )
            return ClassificationResult(
                is_member=True,
                explanation=f"Enoyl reductase: {cofactor}-dependent desaturation (C-C → C=C)"
            )

        # Reduction direction
        has_reduced_left = any(
            p.chebi_id in [CHEBI_NADH, CHEBI_NADPH]
            for p in reaction.left_participants
        )
        has_oxidized_right = any(
            p.chebi_id in [CHEBI_NAD_PLUS, CHEBI_NADP_PLUS]
            for p in reaction.right_participants
        )

        if has_reduced_left and has_oxidized_right:
            cofactor = (
                "NADH"
                if any(p.chebi_id == CHEBI_NADH for p in reaction.left_participants)
                else "NADPH"
            )
            return ClassificationResult(
                is_member=True,
                explanation=f"Enoyl reductase: {cofactor}-dependent saturation (C=C → C-C)"
            )

        # Default case - parent already confirmed CH-CH pattern
        return ClassificationResult(
            is_member=True,
            explanation="Enoyl reductase: NAD(P)-dependent CH-CH oxidoreduction"
        )
