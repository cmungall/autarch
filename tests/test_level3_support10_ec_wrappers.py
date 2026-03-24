"""Focused tests for level-3 EC wrappers with transitive support >= 10."""

from __future__ import annotations

from autarch.classifier import ReactionClassifier
from autarch.ontology.acting_on_the_aldehyde_or_oxo_group_of_donors_with_nad_or_nadp_as_acceptor import (
    ActingOnTheAldehydeOrOxoGroupOfDonorsWithNADOrNADPAsAcceptor,
)
from autarch.ontology.acting_on_the_ch_ch_group_of_donors_with_other_acceptors import (
    ActingOnTheCHCHGroupOfDonorsWithOtherAcceptors,
)
from autarch.ontology.carbon_sulfur_lyases import CarbonSulfurLyases
from autarch.ontology.cis_trans_isomerases import CisTransIsomerases
from autarch.ontology.intramolecular_lyases import IntramolecularLyases
from autarch.ontology.phosphotransferases_phosphomutases import (
    PhosphotransferasesPhosphomutases,
)
from autarch.ontology.transketolases_and_transaldolases import (
    TransketolasesAndTransaldolases,
)
from autarch.ontology.transposing_c_c_bonds import TransposingCCBonds


ALL_LEVEL3_SUPPORT10_WRAPPERS = {
    "ActingOnTheAldehydeOrOxoGroupOfDonorsWithNADOrNADPAsAcceptor",
    "TransferringOtherGroups",
    "ProteinSerineThreonineKinases",
    "ActingOnTheCHCHGroupOfDonorsWithOtherAcceptors",
    "CarbonSulfurLyases",
    "ActingOnTheCHOHGroupOfDonorsWithOtherAcceptors",
    "ActingOnPairedDonorsWithIncorporationOrReductionOfMolecularOxygenMiscellaneous",
    "IntramolecularLyases",
    "TransposingCCBonds",
    "OtherCarbonCarbonLyases",
    "TransketolasesAndTransaldolases",
    "TranslocationOfOtherCompoundsLinkedToTheHydrolysisOfANucleosideTriphosphate",
    "PhosphotransferasesPhosphomutases",
    "ActingOnCarbonNitrogenBondsOtherThanPeptideBondsInOtherCompounds",
    "OtherIntramolecularOxidoreductases",
    "CisTransIsomerases",
    "ReactionXHYHXYWithOxygenAsAcceptor",
    "TransferringAminoGroups",
    "ActingOnCHOrCH2GroupsWithNADOrNADPAsAcceptor",
    "ActingOnPairedDonorsWithIncorporationOrReductionOfMolecularOxygenWith2OxoglutarateAsOneDonorAndTheOtherDehydrogenated",
}


def test_level3_support10_wrappers_are_discoverable() -> None:
    classifier = ReactionClassifier()
    assert ALL_LEVEL3_SUPPORT10_WRAPPERS <= set(classifier.reaction_classes)


def test_ec_121_wrapper_has_expected_child() -> None:
    assert (
        ActingOnTheAldehydeOrOxoGroupOfDonorsWithNADOrNADPAsAcceptor.CHILD_CLASSES
        == (
            ActingOnTheAldehydeOrOxoGroupOfDonorsWithNADOrNADPAsAcceptor.CHILD_CLASSES[0],
        )
    )
    assert (
        ActingOnTheAldehydeOrOxoGroupOfDonorsWithNADOrNADPAsAcceptor.CHILD_CLASSES[0].__name__
        == "OxidoreductaseActingOnTheAldehydeOrOxoGroupOfDonorsNADOrNADPAsAcceptor"
    )


def test_ec_1399_wrapper_has_expected_child() -> None:
    assert (
        ActingOnTheCHCHGroupOfDonorsWithOtherAcceptors.CHILD_CLASSES[0].__name__
        == "Quinoline2Oxidoreductase"
    )


def test_ec_441_wrapper_has_expected_child() -> None:
    assert CarbonSulfurLyases.CHILD_CLASSES[0].__name__ == "CarbonSulfurLyase"


def test_ec_551_wrapper_has_expected_child() -> None:
    assert IntramolecularLyases.CHILD_CLASSES[0].__name__ == "IntramolecularLyase"


def test_ec_533_wrapper_has_expected_child() -> None:
    assert TransposingCCBonds.CHILD_CLASSES[0].__name__ == "Delta3Delta2EnoylCoAIsomerase"


def test_ec_221_wrapper_has_expected_child() -> None:
    assert (
        TransketolasesAndTransaldolases.CHILD_CLASSES[0].__name__
        == "AldehydeKetoneTransferase"
    )


def test_ec_542_wrapper_has_expected_child() -> None:
    assert (
        PhosphotransferasesPhosphomutases.CHILD_CLASSES[0].__name__
        == "IntramolecularPhosphotransferase"
    )


def test_ec_521_wrapper_has_expected_child() -> None:
    assert CisTransIsomerases.CHILD_CLASSES[0].__name__ == "PeptidylProlylCisTransIsomerase"
