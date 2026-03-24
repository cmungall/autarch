"""Focused tests for missing EC level-3 shells with support >= 6."""

from __future__ import annotations

from autarch.classifier import ReactionClassifier
from autarch.ontology.ec_prefix_aggregate import aggregate_child_classes


ALL_LEVEL3_SUPPORT6_EC_SHELLS = {
    "ActingOnGTPInvolvedInCellularAndSubcellularMovement",
    "ActingOnIronSulfurProteinsAsDonorsWithNADOrNADPAsAcceptor",
    "ActingOnOtherCompounds",
    "ActingOnOtherNitrogenousCompoundsAsDonorsWithACytochromeAsAcceptor",
    "ActingOnTheCHNHGroupOfDonorsWithOtherAcceptors",
    "CycloLigases",
    "FormingCarbonCarbonBonds",
    "InterconvertingKetoAndEnolGroups",
    "LyasesActingOnAmidesAmidinesEtc",
    "OxidizingMetalIonsWithNADOrNADPAsAcceptor",
    "PhosphorusOxygenLyases",
    "TransferringHydroxyGroups",
}


def test_level3_support6_ec_shells_are_discoverable() -> None:
    classifier = ReactionClassifier()
    assert ALL_LEVEL3_SUPPORT6_EC_SHELLS <= set(classifier.reaction_classes)


def test_level3_support6_ec_shells_currently_have_no_curated_children() -> None:
    expected_empty = {
        "1.16.1.-",
        "1.18.1.-",
        "1.5.99.-",
        "1.7.2.-",
        "3.6.5.-",
        "4.3.2.-",
        "4.6.1.-",
        "5.1.99.-",
        "5.3.2.-",
        "5.4.4.-",
        "6.3.3.-",
        "6.4.1.-",
    }

    for ec_prefix in expected_empty:
        assert aggregate_child_classes(ec_prefix, "Dummy") == ()
