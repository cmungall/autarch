"""leucyl-tRNA--protein transferase activity.

Catalysis of transfer of a leucyl residue from leucyl-tRNA onto the N terminus
of a lysyl- or arginyl-protein substrate.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_LEUCYL_TRNA = "CHEBI:78494"
CHEBI_TRNA_LEU = "CHEBI:78442"
SUBSTRATE_TO_PRODUCT = {
    "CHEBI:65249": "CHEBI:133043",
    "CHEBI:64719": "CHEBI:133044",
}


class LeucylTRNAProteinTransferaseActivity(ReactionClass):
    """leucyl-tRNA--protein transferase activity."""

    GO_ID = "GO:0008914"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = [p.chebi_id for p in reaction.left_participants if p.chebi_id]
        right_ids = [p.chebi_id for p in reaction.right_participants if p.chebi_id]
        if CHEBI_LEUCYL_TRNA not in left_ids or CHEBI_TRNA_LEU not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires leucyl-tRNA substrate and released tRNA(Leu)")
        if CHEBI_H_PLUS not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires hydron release")
        left_core = [chebi_id for chebi_id in left_ids if chebi_id != CHEBI_LEUCYL_TRNA]
        right_core = [chebi_id for chebi_id in right_ids if chebi_id not in {CHEBI_TRNA_LEU, CHEBI_H_PLUS}]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(is_member=False, explanation="Expected one protein acceptor and one leucylated protein product")
        if SUBSTRATE_TO_PRODUCT.get(left_core[0]) != right_core[0]:
            return ClassificationResult(is_member=False, explanation="Protein branch does not match a supported leucyl transfer pair")
        return ClassificationResult(is_member=True, explanation="leucyl-tRNA--protein transferase activity: transfer of leucine from leucyl-tRNA to an N-terminal protein residue")
