"""Tests for ProteinDeglycase classifier (EC 3.5.1.124, GO:0036524)."""

from pathlib import Path

import pytest

from autarch.io import load_reaction_from_yaml
from autarch.ontology import ProteinDeglycase

TEST_DATA_DIR = Path(__file__).parent / "input"


class TestProteinDeglycase:
    """Test protein deglycase classification."""
    
    @pytest.fixture
    def classifier(self):
        """Create a protein deglycase classifier instance."""
        return ProteinDeglycase()
    
    def test_methylglyoxal_removal_from_lysine(self, classifier):
        """Test RHEA:49552 - methylglyoxal removal from lysine."""
        reaction = load_reaction_from_yaml(
            TEST_DATA_DIR / "rhea_49552_deglycase_lysine.yaml"
        )
        result = classifier.check_membership_impl(reaction)

        assert result.is_member, (
            f"Should classify methylglyoxal removal as deglycase: {result.explanation}"
        )
        # Check explanation mentions deglycase or glycation
        assert "deglycase" in result.explanation.lower() or "glycation" in result.explanation.lower()
    
    def test_glyoxal_removal_from_lysine(self, classifier):
        """Test RHEA:57192 - glyoxal removal from lysine."""
        reaction = load_reaction_from_yaml(
            TEST_DATA_DIR / "rhea_57192_deglycase_glycolate.yaml"
        )
        result = classifier.check_membership_impl(reaction)
        
        assert result.is_member, (
            f"Should classify glyoxal removal as deglycase: {result.explanation}"
        )
        assert "glyoxal" in result.explanation.lower() or "deglycase" in result.explanation.lower()
        assert "glycolate" in result.explanation.lower()
    
    def test_regular_hydrolysis_not_deglycase(self, classifier):
        """Test that regular protein hydrolysis is NOT deglycase."""
        reaction = load_reaction_from_yaml(TEST_DATA_DIR / "protease.yaml")
        result = classifier.check_membership_impl(reaction)
        
        assert not result.is_member, (
            f"Regular proteolysis should NOT be deglycase: {result.explanation}"
        )
    
    def test_ester_hydrolysis_not_deglycase(self, classifier):
        """Test that ester hydrolysis is NOT deglycase."""
        reaction = load_reaction_from_yaml(TEST_DATA_DIR / "ester_hydrolysis.yaml")
        result = classifier.check_membership_impl(reaction)
        
        assert not result.is_member, (
            f"Ester hydrolysis should NOT be deglycase: {result.explanation}"
        )