"""protein histidine phosphatase activity.

Catalysis of hydrolysis of phosphohistidine residues in proteins.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_PHOSPHATE
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_HISTIDYL_PROTEIN = "CHEBI:29979"
PHOSPHOHISTIDINE_IDS = {"CHEBI:83586", "CHEBI:64837"}


class ProteinHistidinePhosphataseActivity(ReactionClass):
    """protein histidine phosphatase activity."""

    GO_ID = "GO:0101006"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = {p.chebi_id for p in reaction.left_participants if p.chebi_id}
        right_ids = {p.chebi_id for p in reaction.right_participants if p.chebi_id}
        if not (PHOSPHOHISTIDINE_IDS & left_ids):
            return ClassificationResult(is_member=False, explanation="Requires phosphohistidyl-protein substrate")
        if CHEBI_H2O not in left_ids or CHEBI_HISTIDYL_PROTEIN not in right_ids or CHEBI_PHOSPHATE not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires water consumption with histidyl-protein and phosphate release")
        return ClassificationResult(is_member=True, explanation="protein histidine phosphatase activity: hydrolysis of a phosphohistidine residue in protein")
