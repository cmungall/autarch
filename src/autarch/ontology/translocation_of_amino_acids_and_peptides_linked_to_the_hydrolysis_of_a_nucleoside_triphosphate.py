"""translocation of amino acids and peptides linked to the hydrolysis of a nucleoside triphosphate.

Precision-biased EC-level translocase classifier for ATP/GTP-coupled transport
of amino acids, peptides, and directly represented peptide-like cargo.
"""

from typing import ClassVar, Optional

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.translocase_linked_to_hydrolysis import (
    TranslocaseLinkedToHydrolysis,
)
from autarch.ontology.transport_utils import (
    is_amino_acid_or_peptide,
    reactive_transported_pairs,
)


class TranslocationOfAminoAcidsAndPeptidesLinkedToTheHydrolysisOfANucleosideTriphosphate(
    TranslocaseLinkedToHydrolysis
):
    """translocation of amino acids and peptides linked to the hydrolysis of a nucleoside triphosphate.

    Precision-biased EC-level translocase classifier for ATP/GTP-coupled
    transport of amino acids, peptides, and directly represented peptide-like
    cargo.
    """

    GO_ID: ClassVar[Optional[str]] = None
    EC_NUMBER_PREFIX: ClassVar[Optional[str]] = "7.4.2.-"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for ATP/GTP-coupled translocation of amino acids or peptides."""
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        transported_pairs = reactive_transported_pairs(reaction)
        if not transported_pairs:
            return ClassificationResult(
                is_member=False,
                explanation="No transported substrate beyond ATPase-coupling participants",
            )

        if not all(is_amino_acid_or_peptide(left_participant) for left_participant, _ in transported_pairs):
            return ClassificationResult(
                is_member=False,
                explanation="Transported substrate is not consistently amino-acid or peptide-like",
            )

        transport_info = [
            f"{left_participant.chebi_id or 'unknown'}: "
            f"{left_participant.location or 'unknown'}→{right_participant.location or 'unknown'}"
            for left_participant, right_participant in transported_pairs
        ]
        return ClassificationResult(
            is_member=True,
            explanation=(
                "Translocation of amino acids and peptides linked to nucleoside-triphosphate "
                f"hydrolysis: {', '.join(transport_info)}"
            ),
        )
