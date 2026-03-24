"""guanylyltransferase activity."""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_CDP, CHEBI_CTP, CHEBI_GDP, CHEBI_GTP, CHEBI_UDP, CHEBI_UTP
from autarch.ontology.nucleotidyltransferase import Nucleotidyltransferase
from autarch.ontology.reaction import IDENTIFIER_AND_POLYMER_FRIENDLY_EVIDENCE


class Guanylyltransferase(Nucleotidyltransferase):
    """guanylyltransferase activity."""

    GO_ID = "GO:0070568"
    EVALUATION_EVIDENCE = IDENTIFIER_AND_POLYMER_FRIENDLY_EVIDENCE
    SPECIFIC_DONORS = {CHEBI_GTP, CHEBI_GDP}
    OTHER_DONORS = {CHEBI_CTP, CHEBI_CDP, CHEBI_UTP, CHEBI_UDP}

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Restrict nucleotidyl transfer to guanylyl donors."""
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        all_ids = {participant.chebi_id for participant in reaction.all_participants() if participant.chebi_id}
        if not (all_ids & self.SPECIFIC_DONORS):
            return ClassificationResult(is_member=False, explanation="No guanylyl donor detected")
        if all_ids & self.OTHER_DONORS:
            return ClassificationResult(is_member=False, explanation="Contains a non-guanylyl nucleotide donor")

        return ClassificationResult(
            is_member=True,
            explanation="Guanylyltransferase: nucleotidyl transfer using a guanine nucleotide donor",
        )
