"""acetylajmaline esterase activity.

Catalysis of hydrolytic deacetylation of acetylajmaline alkaloids.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_ACETATE = "CHEBI:30089"
SUBSTRATE_TO_PRODUCT = {
    "CHEBI:58679": "CHEBI:58567",
    "CHEBI:77725": "CHEBI:77618",
}


class AcetylajmalineEsteraseActivity(ReactionClass):
    """acetylajmaline esterase activity."""

    GO_ID = "GO:0033879"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = [p.chebi_id for p in reaction.left_participants if p.chebi_id]
        right_ids = [p.chebi_id for p in reaction.right_participants if p.chebi_id]
        if CHEBI_H2O not in left_ids or CHEBI_ACETATE not in right_ids or CHEBI_H_PLUS not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires hydrolysis with acetate and hydron release")
        left_core = [chebi_id for chebi_id in left_ids if chebi_id != CHEBI_H2O]
        right_core = [chebi_id for chebi_id in right_ids if chebi_id not in {CHEBI_ACETATE, CHEBI_H_PLUS}]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(is_member=False, explanation="Expected one acetylajmaline substrate and one deacetylated alkaloid product")
        if SUBSTRATE_TO_PRODUCT.get(left_core[0]) != right_core[0]:
            return ClassificationResult(is_member=False, explanation="Alkaloid branch does not match a supported acetylajmaline esterase pair")
        return ClassificationResult(is_member=True, explanation="acetylajmaline esterase activity: hydrolytic deacetylation of acetylajmaline alkaloids")
