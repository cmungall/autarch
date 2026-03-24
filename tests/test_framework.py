"""Tests for the classifier framework using fake classifiers.

These tests ensure the framework works correctly without depending on
specific chemistry implementations that may change.
"""

import pytest

from autarch.classifier import ReactionClassifier
from autarch.datamodel import Reaction, Participant
from tests.fake_classifiers import (
    FakeHydrolase,
    FakeOxidoreductase,
    FakeTransferase,
    FakeLyase,
    AlwaysTrue,
    AlwaysFalse,
)


class TestClassifierFramework:
    """Test the classifier framework with fake classifiers."""

    @pytest.fixture
    def fake_classifier(self):
        """Create a classifier with only fake classes."""
        classifier = ReactionClassifier()
        # Clear real classes and add fake ones
        classifier.reaction_classes = {
            "FakeHydrolase": FakeHydrolase,
            "FakeOxidoreductase": FakeOxidoreductase,
            "FakeTransferase": FakeTransferase,
            "FakeLyase": FakeLyase,
            "AlwaysTrue": AlwaysTrue,
            "AlwaysFalse": AlwaysFalse,
        }
        return classifier

    @pytest.fixture
    def fake_water_reaction(self):
        """Reaction with FAKE_WATER."""
        return Reaction(
            left_participants=[
                Participant(smiles="FAKE_WATER"),
                Participant(smiles="ABC"),
            ],
            right_participants=[Participant(smiles="AB"), Participant(smiles="C")],
        )

    @pytest.fixture
    def fake_nad_reaction(self):
        """Reaction with FAKE_NAD."""
        return Reaction(
            left_participants=[
                Participant(smiles="SUBSTRATE"),
                Participant(smiles="FAKE_NAD+"),
            ],
            right_participants=[
                Participant(smiles="PRODUCT"),
                Participant(smiles="FAKE_NADH"),
            ],
        )

    @pytest.fixture
    def balanced_reaction(self):
        """Reaction with equal participants."""
        return Reaction(
            left_participants=[Participant(smiles="A"), Participant(smiles="B")],
            right_participants=[Participant(smiles="C"), Participant(smiles="D")],
        )

    @pytest.fixture
    def fragmentation_reaction(self):
        """Reaction with fragmentation."""
        return Reaction(
            left_participants=[Participant(smiles="ABC")],
            right_participants=[
                Participant(smiles="A"),
                Participant(smiles="B"),
                Participant(smiles="C"),
            ],
        )

    def test_classifier_discovers_fake_classes(self, fake_classifier):
        """Test that fake classes are discovered."""
        assert "FakeHydrolase" in fake_classifier.reaction_classes
        assert "FakeOxidoreductase" in fake_classifier.reaction_classes
        assert len(fake_classifier.reaction_classes) == 6

    def test_fake_hydrolase_classification(self, fake_classifier, fake_water_reaction):
        """Test FakeHydrolase classification."""
        results = fake_classifier.classify(fake_water_reaction)
        assert "FakeHydrolase" in results
        assert results["FakeHydrolase"].is_member is True
        assert results["FakeHydrolase"].explanation == "Has FAKE_WATER reactant"

    def test_fake_oxidoreductase_classification(
        self, fake_classifier, fake_nad_reaction
    ):
        """Test FakeOxidoreductase classification."""
        results = fake_classifier.classify(fake_nad_reaction)
        assert "FakeOxidoreductase" in results
        assert results["FakeOxidoreductase"].is_member is True
        assert results["FakeOxidoreductase"].explanation == "Has FAKE_NAD cofactor"

    def test_fake_transferase_classification(self, fake_classifier, balanced_reaction):
        """Test FakeTransferase classification."""
        results = fake_classifier.classify(balanced_reaction)
        assert "FakeTransferase" in results
        assert results["FakeTransferase"].is_member is True
        assert results["FakeTransferase"].explanation == "Equal participant count"

    def test_fake_lyase_classification(self, fake_classifier, fragmentation_reaction):
        """Test FakeLyase classification."""
        results = fake_classifier.classify(fragmentation_reaction)
        assert "FakeLyase" in results
        assert results["FakeLyase"].is_member is True
        assert results["FakeLyase"].explanation == "Fragmentation pattern"

    def test_always_true_classifier(self, fake_classifier):
        """Test AlwaysTrue classifier."""
        reaction = Reaction(left_participants=[], right_participants=[])
        results = fake_classifier.classify(reaction)
        assert "AlwaysTrue" in results
        assert results["AlwaysTrue"].is_member is True

    def test_always_false_classifier(self, fake_classifier):
        """Test AlwaysFalse classifier."""
        reaction = Reaction(left_participants=[], right_participants=[])
        results = fake_classifier.classify(reaction)
        assert "AlwaysFalse" in results
        assert results["AlwaysFalse"].is_member is False

    def test_classify_with_specific_class(self, fake_classifier, fake_water_reaction):
        """Test classification with a specific class."""
        result = fake_classifier.classify_with_class(fake_water_reaction, FakeHydrolase)
        assert result.is_member is True

        result = fake_classifier.classify_with_class(
            fake_water_reaction, FakeOxidoreductase
        )
        assert result.is_member is False

    def test_get_matching_classes(self, fake_classifier, fake_water_reaction):
        """Test getting all matching classes."""
        matching = fake_classifier.get_matching_classes(fake_water_reaction)
        assert "FakeHydrolase" in matching
        assert "AlwaysTrue" in matching
        assert "AlwaysFalse" not in matching
        assert "FakeOxidoreductase" not in matching

    def test_multiple_classifications(self, fake_classifier):
        """Test a reaction matching multiple classes."""
        # This reaction has FAKE_WATER and equal participants
        reaction = Reaction(
            left_participants=[
                Participant(smiles="FAKE_WATER"),
                Participant(smiles="X"),
            ],
            right_participants=[Participant(smiles="Y"), Participant(smiles="Z")],
        )

        fake_classifier.classify(reaction)
        matching = fake_classifier.get_matching_classes(reaction)

        # Should match FakeHydrolase, FakeTransferase, and AlwaysTrue
        assert len(matching) == 3
        assert "FakeHydrolase" in matching
        assert "FakeTransferase" in matching
        assert "AlwaysTrue" in matching

    def test_no_matching_classes(self, fake_classifier):
        """Test a reaction matching no classes (except AlwaysTrue)."""
        reaction = Reaction(
            left_participants=[
                Participant(smiles="A"),
                Participant(smiles="B"),
                Participant(smiles="C"),
            ],
            right_participants=[Participant(smiles="D")],
        )

        matching = fake_classifier.get_matching_classes(reaction)
        # Only AlwaysTrue should match
        assert matching == ["AlwaysTrue"]

    def test_empty_reaction(self, fake_classifier):
        """Test classification of empty reaction."""
        reaction = Reaction(left_participants=[], right_participants=[])
        results = fake_classifier.classify(reaction)

        # FakeTransferase should match (0 == 0)
        assert results["FakeTransferase"].is_member is True
        # FakeLyase should not match (0 < 0 is false)
        assert results["FakeLyase"].is_member is False
