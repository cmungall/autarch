"""3-hydroxyphenylacetate 6-hydroxylase activity.

Catalysis of NADH/NADPH-dependent hydroxylation of 3-hydroxyphenylacetate to
homogentisate.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS, CHEBI_NADH, CHEBI_NADPH, CHEBI_NADP_PLUS, CHEBI_NAD_PLUS, CHEBI_O2
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_3_HYDROXYPHENYLACETATE = "CHEBI:58149"
CHEBI_HOMOGENTISATE = "CHEBI:16169"


class ThreeHydroxyphenylacetate6HydroxylaseActivity(ReactionClass):
    """3-hydroxyphenylacetate 6-hydroxylase activity."""

    GO_ID = "GO:0047094"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE
    COFACTOR_PAIRS = {
        CHEBI_NADPH: CHEBI_NADP_PLUS,
        CHEBI_NADH: CHEBI_NAD_PLUS,
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = [p.chebi_id for p in reaction.left_participants if p.chebi_id]
        right_ids = [p.chebi_id for p in reaction.right_participants if p.chebi_id]
        if CHEBI_3_HYDROXYPHENYLACETATE not in left_ids or CHEBI_HOMOGENTISATE not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires 3-hydroxyphenylacetate substrate and homogentisate product")
        if CHEBI_O2 not in left_ids or CHEBI_H_PLUS not in left_ids or CHEBI_H2O not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires dioxygen, hydron, and water in the hydroxylase stoichiometry")
        if any(reduced in left_ids and oxidized in right_ids for reduced, oxidized in self.COFACTOR_PAIRS.items()):
            return ClassificationResult(is_member=True, explanation="3-hydroxyphenylacetate 6-hydroxylase activity: NAD(P)H-dependent hydroxylation of 3-hydroxyphenylacetate")
        return ClassificationResult(is_member=False, explanation="Requires NADH/NADPH oxidation to the matching oxidized cofactor")
