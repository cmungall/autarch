"""steroid dehydrogenase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor.

Catalysis of the oxidation of a steroid alcohol group with NAD(P)+ as the
acceptor.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.biochemical_context_utils import is_steroid_like
from autarch.ontology.oxidoreductase_acting_on_the_ch_oh_group_of_donors_nad_or_nadp_as_acceptor import (
    OxidoreductaseActingOnTheCHOHGroupOfDonorsNADOrNADPAsAcceptor,
)


class SteroidDehydrogenaseActivityActingOnTheCHOHGroupOfDonorsNADOrNADPAsAcceptor(
    OxidoreductaseActingOnTheCHOHGroupOfDonorsNADOrNADPAsAcceptor
):
    """steroid dehydrogenase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor.

    Catalysis of the oxidation of a steroid alcohol group with NAD(P)+ as the
    acceptor.
    """

    GO_ID = "GO:0033764"
    EC_BROAD_XREFS = ["1.1.1.-"]  # GO xref is broadMatch, not exact EC equivalence

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        parent = super().check_membership_impl(reaction)
        if not parent.is_member:
            return parent
        if not self._has_steroid_pair(reaction.left_participants, reaction.right_participants) and not self._has_steroid_pair(reaction.right_participants, reaction.left_participants):
            return ClassificationResult(is_member=False, explanation="No conserved steroid-like substrate/product pair detected")
        return ClassificationResult(
            is_member=True,
            explanation="Steroid dehydrogenase activity with NAD(P)+ acceptor: steroid-like alcohol/carbonyl interconversion",
        )

    @staticmethod
    def _has_steroid_pair(left_participants, right_participants) -> bool:
        return any(is_steroid_like(left) for left in left_participants) and any(is_steroid_like(right) for right in right_participants)
