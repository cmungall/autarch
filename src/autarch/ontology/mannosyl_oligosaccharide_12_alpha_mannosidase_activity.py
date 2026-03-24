"""mannosyl-oligosaccharide 1,2-alpha-mannosidase activity.

Catalysis of hydrolytic trimming of terminal alpha-1,2-mannose residues from
oligomannose N-glycan substrates.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H2O
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

SUBSTRATE_IDS = {"CHEBI:139493", "CHEBI:60628"}
PRODUCT_ID = "CHEBI:59087"
CHEBI_BETA_D_MANNOSE = "CHEBI:28563"


class MannosylOligosaccharide12AlphaMannosidaseActivity(ReactionClass):
    """mannosyl-oligosaccharide 1,2-alpha-mannosidase activity."""

    GO_ID = "GO:0004571"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left = reaction.left_participants
        right = reaction.right_participants
        left_ids = {p.chebi_id for p in left if p.chebi_id}
        right_ids = {p.chebi_id for p in right if p.chebi_id}
        if not (SUBSTRATE_IDS & left_ids):
            return ClassificationResult(is_member=False, explanation="Requires a supported oligomannose glycoprotein substrate")
        if PRODUCT_ID not in right_ids or CHEBI_BETA_D_MANNOSE not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires trimmed oligomannose product and released beta-D-mannose")
        if CHEBI_H2O not in left_ids:
            return ClassificationResult(is_member=False, explanation="Requires hydrolytic water consumption")
        return ClassificationResult(is_member=True, explanation="mannosyl-oligosaccharide 1,2-alpha-mannosidase activity: hydrolytic trimming of terminal alpha-1,2-mannose residues")
