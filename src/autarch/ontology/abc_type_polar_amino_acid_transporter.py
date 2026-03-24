"""ABC-type polar-amino-acid transporter.

ATP-dependent transmembrane transport of supported polar amino acids.
"""

from typing import ClassVar, Optional

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.translocase_linked_to_hydrolysis import TranslocaseLinkedToHydrolysis
from autarch.ontology.transport_utils import is_supported_polar_amino_acid, reactive_transported_pairs


class ABCTypePolarAminoAcidTransporter(TranslocaseLinkedToHydrolysis):
    """ABC-type polar-amino-acid transporter.

    ATP-dependent transmembrane transport of supported polar amino acids.
    """

    GO_ID: ClassVar[Optional[str]] = None
    EC_NUMBER_PREFIX: ClassVar[Optional[str]] = "7.4.2.1"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        transported_pairs = reactive_transported_pairs(reaction)
        if not transported_pairs:
            return ClassificationResult(is_member=False, explanation="No transported substrate beyond ATPase-coupling participants")

        if not all(is_supported_polar_amino_acid(left_participant) for left_participant, _ in transported_pairs):
            return ClassificationResult(is_member=False, explanation="Transported substrate is not a supported polar amino acid")

        return ClassificationResult(is_member=True, explanation="ABC-type polar-amino-acid transporter: ATP-dependent transport of a supported polar amino acid")
