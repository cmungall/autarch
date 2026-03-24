"""serine-tRNA ligase activity.

Catalysis of ATP-dependent loading of serine onto tRNA(Ser) or tRNA(Sec).
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_AMP, CHEBI_ATP, CHEBI_H_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_DIPHOSPHATE = "CHEBI:33019"
CHEBI_SERINE = "CHEBI:33384"
TRNA_TO_AMINOACYL = {
    "CHEBI:78442": "CHEBI:78533",
}


class SerineTRNALigaseActivity(ReactionClass):
    """serine-tRNA ligase activity."""

    GO_ID = "GO:0004828"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = [p.chebi_id for p in reaction.left_participants if p.chebi_id]
        right_ids = [p.chebi_id for p in reaction.right_participants if p.chebi_id]
        if CHEBI_SERINE not in left_ids or CHEBI_ATP not in left_ids:
            return ClassificationResult(is_member=False, explanation="Requires serine and ATP substrates")
        if CHEBI_AMP not in right_ids or CHEBI_DIPHOSPHATE not in right_ids or CHEBI_H_PLUS not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires AMP, diphosphate, and hydron products")
        left_core = [chebi_id for chebi_id in left_ids if chebi_id not in {CHEBI_SERINE, CHEBI_ATP}]
        right_core = [chebi_id for chebi_id in right_ids if chebi_id not in {CHEBI_AMP, CHEBI_DIPHOSPHATE, CHEBI_H_PLUS}]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(is_member=False, explanation="Expected one tRNA substrate and one seryl-tRNA product")
        if TRNA_TO_AMINOACYL.get(left_core[0]) != right_core[0]:
            return ClassificationResult(is_member=False, explanation="tRNA branch does not match a supported serine ligation pair")
        return ClassificationResult(is_member=True, explanation="serine-tRNA ligase activity: ATP-dependent loading of serine onto tRNA")
