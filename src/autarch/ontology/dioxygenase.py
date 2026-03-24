"""dioxygenase activity."""

from autarch.datamodel import Reaction
from autarch.ontology.ec_prefix_aggregate import (
    aggregate_explicit_child_membership,
    aggregate_explicit_child_supports_evaluation,
)
from autarch.ontology.nitroarene_dioxygenase import NitroareneDioxygenase
from autarch.ontology.oxidoreductase_acting_on_paired_donors_with_incorporation_or_reduction_of_molecular_oxygen_nadh_or_nadph_as_one_donor_and_incorporation_of_two_atoms_of_oxygen_into_one_donor import (
    OxidoreductaseActingOnPairedDonorsWithIncorporationOrReductionOfMolecularOxygenNADHOrNADPHAsOneDonorAndIncorporationOfTwoAtomsOfOxygenIntoOneDonor,
)
from autarch.ontology.oxidoreductase_acting_on_single_donors_with_incorporation_of_molecular_oxygen_incorporation_of_two_atoms_of_oxygen import (
    OxidoreductaseActingOnSingleDonorsWithIncorporationOfMolecularOxygenIncorporationOfTwoAtomsOfOxygen,
)
from autarch.ontology.reaction import ReactionClass
from autarch.ontology.two_oxoglutarate_dependent_dioxygenase import (
    TwoOxoglutarateDependentDioxygenase,
)


class Dioxygenase(ReactionClass):
    """dioxygenase activity."""

    GO_ID = "GO:0051213"
    CHILD_CLASSES = (
        TwoOxoglutarateDependentDioxygenase,
        OxidoreductaseActingOnSingleDonorsWithIncorporationOfMolecularOxygenIncorporationOfTwoAtomsOfOxygen,
        OxidoreductaseActingOnPairedDonorsWithIncorporationOrReductionOfMolecularOxygenNADHOrNADPHAsOneDonorAndIncorporationOfTwoAtomsOfOxygenIntoOneDonor,
        NitroareneDioxygenase,
    )

    @classmethod
    def supports_evaluation(cls, reaction: Reaction) -> bool:
        return aggregate_explicit_child_supports_evaluation(reaction, cls.CHILD_CLASSES)

    def check_membership_impl(self, reaction: Reaction):
        return aggregate_explicit_child_membership(
            reaction,
            self.CHILD_CLASSES,
            "Dioxygenase",
        )
