"""oxidoreductase acting on the CH-CH group of donors with a flavin as acceptor.

Catalysis of an oxidation-reduction (redox) reaction in which a CH-CH group
acts as a hydrogen or electron donor and reduces a flavin.
"""

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_FAD,
    CHEBI_FADH2,
    CHEBI_FMN,
    CHEBI_FMNH2,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
)
from autarch.ontology.oxidoreductase_acting_on_the_ch_ch_group_of_donors import (
    OxidoreductaseActingOnTheCHCHGroupOfDonors,
)

CHEBI_ETF_OXIDIZED = "CHEBI:83723"
CHEBI_ETF_REDUCED = "CHEBI:83726"
CHEBI_REDUCED_FAD_LIKE = "CHEBI:58307"


class OxidoreductaseActingOnTheCHCHGroupOfDonorsWithAFlavinAsAcceptor(
    OxidoreductaseActingOnTheCHCHGroupOfDonors
):
    """oxidoreductase acting on the CH-CH group of donors with a flavin as acceptor.

    Catalysis of an oxidation-reduction (redox) reaction in which a CH-CH group
    acts as a hydrogen or electron donor and reduces a flavin.
    """

    GO_ID = "GO:0052890"
    EC_NUMBER_PREFIX = "1.3.8.-"

    OXIDIZED_FLAVIN_CHEBIS = {CHEBI_FAD, CHEBI_FMN, CHEBI_ETF_OXIDIZED}
    REDUCED_FLAVIN_CHEBIS = {
        CHEBI_FADH2,
        CHEBI_FMNH2,
        CHEBI_ETF_REDUCED,
        CHEBI_REDUCED_FAD_LIKE,
    }
    NAD_CHEBIS = {CHEBI_NAD_PLUS, CHEBI_NADH, CHEBI_NADP_PLUS, CHEBI_NADPH}

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for CH-CH redox with flavin/ETF as the acceptor."""
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        forward = self._check_direction(
            reaction.left_participants,
            reaction.right_participants,
        )
        if forward.is_member:
            return forward

        reverse = self._check_direction(
            reaction.right_participants,
            reaction.left_participants,
        )
        if reverse.is_member:
            return ClassificationResult(
                is_member=True,
                explanation=f"{reverse.explanation} (reverse reaction orientation)",
            )

        return forward

    @classmethod
    def _check_direction(
        cls,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> ClassificationResult:
        """Evaluate one reaction orientation."""
        has_oxidized_flavin = any(
            participant.chebi_id in cls.OXIDIZED_FLAVIN_CHEBIS for participant in left_participants
        )
        has_reduced_flavin = any(
            participant.chebi_id in cls.REDUCED_FLAVIN_CHEBIS for participant in right_participants
        )
        if not (has_oxidized_flavin and has_reduced_flavin):
            return ClassificationResult(
                is_member=False,
                explanation="Missing oxidized-to-reduced flavin conversion",
            )

        if any(
            participant.chebi_id in cls.NAD_CHEBIS
            for participant in left_participants + right_participants
        ):
            return ClassificationResult(
                is_member=False,
                explanation="Nicotinamide cofactors indicate a different CH-CH oxidoreductase branch",
            )

        return ClassificationResult(
            is_member=True,
            explanation="CH-CH oxidoreductase with a flavin as acceptor",
        )
