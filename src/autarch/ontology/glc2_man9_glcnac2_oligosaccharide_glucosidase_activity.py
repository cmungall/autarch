"""Glc2Man9GlcNAc2 oligosaccharide glucosidase activity.

Catalysis of glucosidic trimming steps in Glc2Man9GlcNAc2-type N-glycan
substrates.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H2O
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_BETA_D_GLUCOSE = "CHEBI:15903"
SUBSTRATE_TO_PRODUCT = {
    "CHEBI:59082": "CHEBI:59080",
    "CHEBI:59080": "CHEBI:139493",
}


class Glc2Man9GlcNAc2OligosaccharideGlucosidaseActivity(ReactionClass):
    """Glc2Man9GlcNAc2 oligosaccharide glucosidase activity."""

    GO_ID = "GO:0106407"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = [p.chebi_id for p in reaction.left_participants if p.chebi_id]
        right_ids = [p.chebi_id for p in reaction.right_participants if p.chebi_id]
        if CHEBI_H2O not in left_ids or CHEBI_BETA_D_GLUCOSE not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires hydrolysis with beta-D-glucose release")
        left_core = [chebi_id for chebi_id in left_ids if chebi_id != CHEBI_H2O]
        right_core = [chebi_id for chebi_id in right_ids if chebi_id != CHEBI_BETA_D_GLUCOSE]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(is_member=False, explanation="Expected one glycoprotein substrate and one trimmed glycoprotein product")
        if SUBSTRATE_TO_PRODUCT.get(left_core[0]) != right_core[0]:
            return ClassificationResult(is_member=False, explanation="Glycan branch does not match a supported Glc2Man9GlcNAc2 glucosidase step")
        return ClassificationResult(is_member=True, explanation="Glc2Man9GlcNAc2 oligosaccharide glucosidase activity: glucosidic trimming of an N-glycan substrate")
