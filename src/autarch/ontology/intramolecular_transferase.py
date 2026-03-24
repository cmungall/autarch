"""intramolecular transferase activity."""

from autarch.datamodel import Reaction
from autarch.ontology.ec_prefix_aggregate import (
    aggregate_explicit_child_membership,
    aggregate_explicit_child_supports_evaluation,
)
from autarch.ontology.intramolecular_phosphotransferase import IntramolecularPhosphotransferase
from autarch.ontology.phosphotransferases_phosphomutases import PhosphotransferasesPhosphomutases
from autarch.ontology.reaction import ReactionClass
from autarch.ontology.transferring_amino_groups import TransferringAminoGroups
from autarch.ontology.transferring_other_groups import TransferringOtherGroups


class IntramolecularTransferase(ReactionClass):
    """intramolecular transferase activity."""

    GO_ID = "GO:0016866"
    EC_NUMBER_PREFIX = "5.4.-.-"
    CHILD_CLASSES = (
        IntramolecularPhosphotransferase,
        PhosphotransferasesPhosphomutases,
        TransferringAminoGroups,
        TransferringOtherGroups,
    )

    @classmethod
    def supports_evaluation(cls, reaction: Reaction) -> bool:
        return aggregate_explicit_child_supports_evaluation(reaction, cls.CHILD_CLASSES)

    def check_membership_impl(self, reaction: Reaction):
        return aggregate_explicit_child_membership(
            reaction,
            self.CHILD_CLASSES,
            "Intramolecular transferase",
        )
