"""acyl-CoA dehydrogenase activity."""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.lipid_utils import double_bond_count, is_fatty_acyl_coa_like
from autarch.ontology.oxidoreductase_acting_on_the_ch_ch_group_of_donors_with_a_flavin_as_acceptor import (
    OxidoreductaseActingOnTheCHCHGroupOfDonorsWithAFlavinAsAcceptor,
)


class AcylCoADehydrogenaseActivity(
    OxidoreductaseActingOnTheCHCHGroupOfDonorsWithAFlavinAsAcceptor
):
    """acyl-CoA dehydrogenase activity."""

    GO_ID = "GO:0003995"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        parent = super().check_membership_impl(reaction)
        if not parent.is_member:
            return parent

        if self._has_matching_pair(reaction.left_participants, reaction.right_participants) or self._has_matching_pair(reaction.right_participants, reaction.left_participants):
            return ClassificationResult(is_member=True, explanation="Acyl-CoA dehydrogenase activity: flavin/ETF-dependent dehydrogenation of an acyl-CoA substrate")
        return ClassificationResult(is_member=False, explanation="No acyl-CoA substrate/product pair with the expected dehydrogenation pattern detected")

    @staticmethod
    def _has_matching_pair(left_participants, right_participants) -> bool:
        left_coa = [participant for participant in left_participants if is_fatty_acyl_coa_like(participant)]
        right_coa = [participant for participant in right_participants if is_fatty_acyl_coa_like(participant)]
        for substrate in left_coa:
            for product in right_coa:
                if double_bond_count(product) == double_bond_count(substrate) + 1:
                    return True
        return False
