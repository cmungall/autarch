"""UDP-xylosyltransferase activity."""

from autarch.molecules import CHEBI_UDP
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE
from autarch.ontology.udp_glycosyltransferase_activity import (
    UDPGlycosyltransferaseActivity,
)


class UDPXylosyltransferaseActivity(UDPGlycosyltransferaseActivity):
    """UDP-xylosyltransferase activity."""

    GO_ID = "GO:0035252"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE
    DONOR_IDS = frozenset({"CHEBI:57632"})
    PRODUCT_IDS = frozenset({CHEBI_UDP})

    def check_membership_impl(self, reaction):
        result = super().check_membership_impl(reaction)
        if not result.is_member:
            return result
        return type(result)(
            is_member=True,
            explanation="UDP-xylosyltransferase activity: UDP-xylose donor transfers a xylosyl group",
        )
