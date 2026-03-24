"""malonate-semialdehyde dehydrogenase (acetylating) activity.

Catalysis of NAD(P)-dependent oxidative decarboxylative acetyl-CoA formation
from malonate semialdehyde.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_CO2, CHEBI_COA, CHEBI_NADH, CHEBI_NADPH, CHEBI_NADP_PLUS, CHEBI_NAD_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_MALONATE_SEMIALDEHYDE = "CHEBI:33190"
CHEBI_ACETYL_COA = "CHEBI:57288"


class MalonateSemialdehydeDehydrogenaseAcetylatingActivity(ReactionClass):
    """malonate-semialdehyde dehydrogenase (acetylating) activity."""

    GO_ID = "GO:0018478"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE
    COFACTOR_PAIRS = {
        CHEBI_NADP_PLUS: CHEBI_NADPH,
        CHEBI_NAD_PLUS: CHEBI_NADH,
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = [p.chebi_id for p in reaction.left_participants if p.chebi_id]
        right_ids = [p.chebi_id for p in reaction.right_participants if p.chebi_id]
        if CHEBI_MALONATE_SEMIALDEHYDE not in left_ids or CHEBI_COA not in left_ids:
            return ClassificationResult(is_member=False, explanation="Requires malonate semialdehyde and coenzyme A substrates")
        if CHEBI_ACETYL_COA not in right_ids or CHEBI_CO2 not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires acetyl-CoA and carbon dioxide products")
        if any(oxidized in left_ids and reduced in right_ids for oxidized, reduced in self.COFACTOR_PAIRS.items()):
            return ClassificationResult(is_member=True, explanation="malonate-semialdehyde dehydrogenase (acetylating) activity: NAD(P)-dependent acetyl-CoA formation from malonate semialdehyde")
        return ClassificationResult(is_member=False, explanation="Requires matched NAD(P) cofactor conversion")
