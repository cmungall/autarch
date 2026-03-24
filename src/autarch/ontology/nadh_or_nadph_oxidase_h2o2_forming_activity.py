"""NAD(P)H oxidase H2O2-forming activity.

Catalysis of hydrogen-peroxide-forming oxidation of NADH or NADPH by dioxygen.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H2O2, CHEBI_H_PLUS, CHEBI_NADH, CHEBI_NADPH, CHEBI_NADP_PLUS, CHEBI_NAD_PLUS, CHEBI_O2
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass


class NADHOrNADPHOxidaseH2O2FormingActivity(ReactionClass):
    """NAD(P)H oxidase H2O2-forming activity."""

    GO_ID = "GO:0016174"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE
    COFACTOR_PAIRS = {
        CHEBI_NADPH: CHEBI_NADP_PLUS,
        CHEBI_NADH: CHEBI_NAD_PLUS,
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = [p.chebi_id for p in reaction.left_participants if p.chebi_id]
        right_ids = [p.chebi_id for p in reaction.right_participants if p.chebi_id]
        if CHEBI_O2 not in left_ids or CHEBI_H_PLUS not in left_ids:
            return ClassificationResult(is_member=False, explanation="Requires dioxygen and hydron on the substrate side")
        if CHEBI_H2O2 not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires hydrogen peroxide as product")
        for reduced, oxidized in self.COFACTOR_PAIRS.items():
            if reduced in left_ids and oxidized in right_ids:
                return ClassificationResult(is_member=True, explanation="NAD(P)H oxidase H2O2-forming activity: oxidation of NADH/NADPH by dioxygen to hydrogen peroxide")
        return ClassificationResult(is_member=False, explanation="Requires NADH or NADPH oxidation to the matching oxidized cofactor")
