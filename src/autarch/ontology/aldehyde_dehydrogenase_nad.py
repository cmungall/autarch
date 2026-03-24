"""aldehyde dehydrogenase (NAD+) activity."""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_NADH, CHEBI_NADPH, CHEBI_NADP_PLUS, CHEBI_NAD_PLUS
from autarch.ontology.oxidoreductase_acting_on_the_aldehyde_or_oxo_group_of_donors_nad_or_nadp_as_acceptor import (
    OxidoreductaseActingOnTheAldehydeOrOxoGroupOfDonorsNADOrNADPAsAcceptor,
)


class AldehydeDehydrogenaseNAD(
    OxidoreductaseActingOnTheAldehydeOrOxoGroupOfDonorsNADOrNADPAsAcceptor
):
    """aldehyde dehydrogenase (NAD+) activity."""

    GO_ID = "GO:0004029"
    EC_NUMBER_PREFIX = "1.2.1.3"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Restrict aldehyde/oxo dehydrogenase chemistry to the NAD branch."""
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        has_nadp = any(
            participant.chebi_id in {CHEBI_NADP_PLUS, CHEBI_NADPH}
            for participant in reaction.all_participants()
        )
        if has_nadp:
            return ClassificationResult(
                is_member=False,
                explanation="Uses NADP(H) cofactor - not NAD+-specific aldehyde dehydrogenase",
            )

        has_nad_left = any(
            participant.chebi_id == CHEBI_NAD_PLUS
            for participant in reaction.left_participants
        )
        has_nadh_right = any(
            participant.chebi_id == CHEBI_NADH
            for participant in reaction.right_participants
        )
        has_nadh_left = any(
            participant.chebi_id == CHEBI_NADH
            for participant in reaction.left_participants
        )
        has_nad_right = any(
            participant.chebi_id == CHEBI_NAD_PLUS
            for participant in reaction.right_participants
        )

        if not ((has_nad_left and has_nadh_right) or (has_nadh_left and has_nad_right)):
            return ClassificationResult(
                is_member=False,
                explanation="No NAD+/NADH conversion pattern for aldehyde dehydrogenase",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Aldehyde dehydrogenase (NAD+): NAD-dependent aldehyde/oxo oxidation-reduction chemistry",
        )
