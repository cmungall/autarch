"""protein kinase activity.

Catalysis of the phosphorylation of an amino acid residue in a protein, usually
according to the reaction: a protein + ATP = a phosphoprotein + ADP.
"""

from autarch.datamodel import ClassificationResult, PolymerType, Reaction
from autarch.molecules import CHEBI_ADP, CHEBI_ATP
from autarch.ontology.reaction import (
    IDENTIFIER_AND_POLYMER_FRIENDLY_EVIDENCE,
    ReactionClass,
)

PROTEIN_POLYMERS = {
    PolymerType.PROTEIN,
    PolymerType.PEPTIDE,
    PolymerType.POLYPEPTIDE,
}


class ProteinKinaseActivity(ReactionClass):
    """protein kinase activity.

    Catalysis of the phosphorylation of an amino acid residue in a protein,
    usually according to the reaction: a protein + ATP = a phosphoprotein +
    ADP.
    """

    GO_ID = "GO:0004672"
    EVALUATION_EVIDENCE = IDENTIFIER_AND_POLYMER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = {participant.chebi_id for participant in reaction.left_participants if participant.chebi_id}
        right_ids = {participant.chebi_id for participant in reaction.right_participants if participant.chebi_id}
        if CHEBI_ATP not in left_ids or CHEBI_ADP not in right_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Requires ATP consumption with ADP production",
            )

        left_proteins = [participant for participant in reaction.left_participants if participant.polymer_type in PROTEIN_POLYMERS]
        right_proteins = [participant for participant in reaction.right_participants if participant.polymer_type in PROTEIN_POLYMERS]
        if not left_proteins or not right_proteins:
            return ClassificationResult(
                is_member=False,
                explanation="Requires protein or peptide participants on both sides of the reaction",
            )

        if not any(
            left.polymer_type == right.polymer_type and left.location == right.location
            for left in left_proteins
            for right in right_proteins
        ):
            return ClassificationResult(
                is_member=False,
                explanation="Protein substrate and product do not form a conserved polymer branch",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Protein kinase activity: ATP-dependent phosphorylation of a protein or peptide substrate",
        )
