"""NAD(P)H dehydrogenase (quinone) activity.

Catalysis of quinone reduction by NADH or NADPH, excluding the broader
"similar acceptor" branch covered by GO:0016655.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.oxidoreductase_acting_on_nadh_or_nadph_quinone_or_similar_compound_as_acceptor import (
    OxidoreductaseActingOnNADHOrNADPHQuinoneOrSimilarCompoundAsAcceptor,
)


class NADHOrNADPHDehydrogenaseQuinone(
    OxidoreductaseActingOnNADHOrNADPHQuinoneOrSimilarCompoundAsAcceptor
):
    """NAD(P)H dehydrogenase (quinone) activity.

    Catalysis of quinone reduction by NADH or NADPH, excluding the broader
    "similar acceptor" branch covered by GO:0016655.
    """

    GO_ID = "GO:0003955"
    EC_NUMBER_PREFIX = "1.6.5.2"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Restrict the parent class to bona fide quinone acceptors."""
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        forward = self._is_specific_quinone_branch(
            reaction.left_participants,
            reaction.right_participants,
        )
        if forward:
            return ClassificationResult(
                is_member=True,
                explanation="NAD(P)H dehydrogenase (quinone): NAD(P)H-dependent reduction of a quinone acceptor",
            )

        reverse = self._is_specific_quinone_branch(
            reaction.right_participants,
            reaction.left_participants,
        )
        if reverse:
            return ClassificationResult(
                is_member=True,
                explanation="NAD(P)H dehydrogenase (quinone): NAD(P)H-dependent reduction of a quinone acceptor (reverse reaction orientation)",
            )

        return ClassificationResult(
            is_member=False,
            explanation="Parent oxidoreductase matched only the broader similar-acceptor branch, not a quinone acceptor",
        )

    def _is_specific_quinone_branch(self, left, right) -> bool:
        donor_pair = self._find_cofactor_pair(left, right)
        if donor_pair is None:
            return False
        donor_id, oxidized_id = donor_pair
        left_reactive = [
            participant
            for participant in self._drop_one_chebi(left, donor_id)
            if participant.chebi_id not in self.LEFT_SPECTATOR_CHEBIS
        ]
        right_reactive = [
            participant
            for participant in self._drop_one_chebi(right, oxidized_id)
            if participant.chebi_id not in self.RIGHT_SPECTATOR_CHEBIS
        ]
        substrate = self._collapse_uniform_participants(left_reactive)
        product = self._collapse_uniform_participants(right_reactive)
        if substrate is None or product is None:
            return False
        if self._is_special_similar_acceptor_pair(substrate, product):
            return False
        return self._is_quinone_like_acceptor(substrate)
