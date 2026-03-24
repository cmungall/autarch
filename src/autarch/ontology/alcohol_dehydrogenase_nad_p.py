"""alcohol dehydrogenase (NADP+) activity."""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_NADH, CHEBI_NADPH, CHEBI_NADP_PLUS, CHEBI_NAD_PLUS
from autarch.ontology.oxidoreductase_acting_on_the_ch_oh_group_of_donors_nad_or_nadp_as_acceptor import (
    OxidoreductaseActingOnTheCHOHGroupOfDonorsNADOrNADPAsAcceptor,
)


class AlcoholDehydrogenaseNADP(
    OxidoreductaseActingOnTheCHOHGroupOfDonorsNADOrNADPAsAcceptor
):
    """alcohol dehydrogenase (NADP+) activity."""

    GO_ID = "GO:0008106"
    EC_NUMBER_PREFIX = "1.1.1.2"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Restrict alcohol dehydrogenase chemistry to the NADP branch."""
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        has_nad = any(
            participant.chebi_id in {CHEBI_NAD_PLUS, CHEBI_NADH}
            for participant in reaction.all_participants()
        )
        if has_nad:
            return ClassificationResult(
                is_member=False,
                explanation="Uses NAD(H) cofactor - not NADP+-specific alcohol dehydrogenase",
            )

        has_nadp_left = any(
            participant.chebi_id == CHEBI_NADP_PLUS
            for participant in reaction.left_participants
        )
        has_nadph_right = any(
            participant.chebi_id == CHEBI_NADPH
            for participant in reaction.right_participants
        )
        has_nadph_left = any(
            participant.chebi_id == CHEBI_NADPH
            for participant in reaction.left_participants
        )
        has_nadp_right = any(
            participant.chebi_id == CHEBI_NADP_PLUS
            for participant in reaction.right_participants
        )

        if not ((has_nadp_left and has_nadph_right) or (has_nadph_left and has_nadp_right)):
            return ClassificationResult(
                is_member=False,
                explanation="No NADP+/NADPH conversion pattern for alcohol dehydrogenase",
            )

        ignored_ids = {
            CHEBI_NAD_PLUS,
            CHEBI_NADH,
            CHEBI_NADP_PLUS,
            CHEBI_NADPH,
            "CHEBI:15378",
        }
        left_non_cofactors = [
            participant for participant in reaction.left_participants if participant.chebi_id not in ignored_ids
        ]
        right_non_cofactors = [
            participant for participant in reaction.right_participants if participant.chebi_id not in ignored_ids
        ]
        if len(left_non_cofactors) != 1 or len(right_non_cofactors) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Not a simple 1->1 alcohol dehydrogenase substrate/product pair",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Alcohol dehydrogenase (NADP+): NADP-dependent alcohol ⇌ carbonyl chemistry",
        )
