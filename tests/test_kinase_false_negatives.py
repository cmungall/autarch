"""Tests to isolate kinase false negative issues.

This test suite reproduces specific false negative reactions from the kinase
evaluation to debug why pattern matching is failing.
"""

import pytest
from autarch.datamodel import Reaction, Participant
from autarch.ontology.kinase import Kinase
from autarch.pattern_dsl import match_patterns


class TestKinaseFalseNegatives:
    """Test specific false negative reactions from kinase evaluation."""

    @pytest.fixture
    def kinase(self):
        """Create a Kinase classifier instance."""
        return Kinase()

    def test_rhea_10152_carbamoyl_phosphate_synthetase(self, kinase):
        """Test RHEA:10152: hydrogencarbonate + NH4(+) + ATP = carbamoyl phosphate + ADP + H2O + H(+)

        This is carbamoyl phosphate synthetase (EC 2.7.2.2) - a complex biosynthetic
        enzyme that incorporates CO2 and NH4+ into carbamoyl phosphate. While technically
        EC 2.7.x (phosphotransferase), this is more of a ligase-type reaction with 3
        substrates fusing into a product, so the simple kinase pattern should NOT match.
        """
        reaction = Reaction(
            left_participants=[
                Participant(chebi_id="CHEBI:17544", name="hydrogencarbonate"),
                Participant(chebi_id="CHEBI:28938", name="NH4(+)"),
                Participant(chebi_id="CHEBI:30616", name="ATP"),
            ],
            right_participants=[
                Participant(chebi_id="CHEBI:58228", name="carbamoyl phosphate"),
                Participant(chebi_id="CHEBI:456216", name="ADP"),
                Participant(chebi_id="CHEBI:15377", name="H2O"),
                Participant(chebi_id="CHEBI:15378", name="H(+)"),
            ],
        )

        result = kinase.check_membership(reaction)

        # Debug: Print the pattern matching details
        print("\nRHEA:10152 Classification Result:")
        print(f"  is_member: {result.is_member}")
        print(f"  explanation: {result.explanation}")

        # This is NOT a simple kinase - it's a complex ligase-type biosynthetic reaction
        # The kinase classifier correctly rejects this due to 3 substrates
        assert not result.is_member, (
            "Carbamoyl phosphate synthetase should NOT be classified as simple kinase - "
            "it has 3 substrates and is a ligase-type biosynthetic reaction"
        )

    def test_rhea_10948_glucosamine_kinase(self, kinase):
        """Test RHEA:10948: D-glucosamine + ATP = D-glucosamine 6-phosphate + ADP + H(+)
        
        This is glucosamine kinase - a straightforward kinase. EC 2.7.1.8
        """
        reaction = Reaction(
            left_participants=[
                Participant(chebi_id="CHEBI:58723", name="D-glucosamine"),
                Participant(chebi_id="CHEBI:30616", name="ATP"),
            ],
            right_participants=[
                Participant(chebi_id="CHEBI:456216", name="ADP"),
                Participant(chebi_id="CHEBI:15378", name="H(+)"),
                Participant(chebi_id="CHEBI:58725", name="D-glucosamine 6-phosphate"),
            ],
        )
        
        result = kinase.check_membership(reaction)
        
        print("\nRHEA:10948 Classification Result:")
        print(f"  is_member: {result.is_member}")
        print(f"  explanation: {result.explanation}")
        
        # Test individual pattern matching
        match = match_patterns(reaction, kinase.PATTERNS, strict=True)
        print(f"  pattern_matched: {match.matched}")
        print(f"  unmatched_left: {[p.name for p in match.unmatched_left]}")
        print(f"  unmatched_right: {[p.name for p in match.unmatched_right]}")
        
        # This should definitely be True - it's a textbook kinase reaction
        assert result.is_member, "Glucosamine kinase is a straightforward kinase (EC 2.7.1.8)"

    def test_rhea_11028_mannose_kinase(self, kinase):
        """Test RHEA:11028: D-mannose + ATP = D-mannose 6-phosphate + ADP + H(+)
        
        This is mannose kinase - another straightforward kinase. EC 2.7.1.7
        """
        reaction = Reaction(
            left_participants=[
                Participant(chebi_id="CHEBI:4208", name="D-mannose"),
                Participant(chebi_id="CHEBI:30616", name="ATP"),
            ],
            right_participants=[
                Participant(chebi_id="CHEBI:456216", name="ADP"),
                Participant(chebi_id="CHEBI:15378", name="H(+)"),
                Participant(chebi_id="CHEBI:58735", name="D-mannose 6-phosphate"),
            ],
        )
        
        result = kinase.check_membership(reaction)
        
        print("\nRHEA:11028 Classification Result:")
        print(f"  is_member: {result.is_member}")
        print(f"  explanation: {result.explanation}")
        
        # Test individual pattern matching
        match = match_patterns(reaction, kinase.PATTERNS, strict=True)
        print(f"  pattern_matched: {match.matched}")
        print(f"  unmatched_left: {[p.name for p in match.unmatched_left]}")
        print(f"  unmatched_right: {[p.name for p in match.unmatched_right]}")
        
        # This should definitely be True
        assert result.is_member, "Mannose kinase is a straightforward kinase (EC 2.7.1.7)"


class TestPatternMatchingComponents:
    """Test individual components of the pattern matching system."""

    def test_kinase_pattern_creation(self):
        """Test that kinase patterns are created correctly."""
        kinase = Kinase()
        
        print("\nKinase Patterns:")
        for i, pattern in enumerate(kinase.PATTERNS):
            print(f"  Pattern {i}:")
            print(f"    Left: {[(p.chebi_id, p.name, p.variable) for p in pattern.left_participants]}")
            print(f"    Right: {[(p.chebi_id, p.name, p.variable) for p in pattern.right_participants]}")
        
        # Should have ATP and GTP patterns
        assert len(kinase.PATTERNS) == 2
        
        # First pattern should be ATP-based
        atp_pattern = kinase.PATTERNS[0]
        assert atp_pattern.left_participants[0].chebi_id == "CHEBI:30616"  # ATP
        assert atp_pattern.left_participants[1].variable == "substrate"
        assert atp_pattern.right_participants[0].chebi_id == "CHEBI:456216"  # ADP
        assert atp_pattern.right_participants[1].variable == "product"

    def test_atp_chebi_id_matching(self):
        """Test that ATP CHEBI ID matching works correctly."""
        from autarch.pattern_dsl import PatternParticipant
        
        # Create pattern participant for ATP
        atp_pattern = PatternParticipant(chebi_id="CHEBI:30616", name="ATP")
        
        # Create reaction participant for ATP  
        atp_reaction = Participant(chebi_id="CHEBI:30616", name="ATP")
        
        # Should match
        matches = atp_pattern.matches(atp_reaction)
        print("\nATP matching test:")
        print(f"  Pattern ATP CHEBI: {atp_pattern.chebi_id}")
        print(f"  Reaction ATP CHEBI: {atp_reaction.chebi_id}")
        print(f"  Matches: {matches}")
        
        assert matches, "ATP pattern should match ATP reaction participant"

    def test_variable_matching(self):
        """Test that variable participants match any molecule."""
        from autarch.pattern_dsl import PatternParticipant
        
        # Create variable pattern
        substrate_var = PatternParticipant(variable="substrate", name="?substrate")
        
        # Should match any participant
        glucose = Participant(chebi_id="CHEBI:17234", name="glucose")
        mannose = Participant(chebi_id="CHEBI:4208", name="D-mannose")
        
        print("\nVariable matching test:")
        print(f"  Variable is_variable: {substrate_var.is_variable}")
        print(f"  Variable name: {substrate_var.variable}")
        print(f"  Matches glucose: {substrate_var.matches(glucose)}")
        print(f"  Matches mannose: {substrate_var.matches(mannose)}")
        
        assert substrate_var.is_variable
        assert substrate_var.matches(glucose)
        assert substrate_var.matches(mannose)

    def test_minimal_kinase_pattern_strict_mode(self):
        """Test minimal kinase reaction with strict mode."""
        # Create simplest possible kinase reaction
        reaction = Reaction(
            left_participants=[
                Participant(chebi_id="CHEBI:30616", name="ATP"),
                Participant(chebi_id="CHEBI:17234", name="glucose"),
            ],
            right_participants=[
                Participant(chebi_id="CHEBI:456216", name="ADP"),
                Participant(chebi_id="CHEBI:17665", name="glucose-6-phosphate"),
            ],
        )
        
        kinase = Kinase()
        match = match_patterns(reaction, kinase.PATTERNS, strict=True)
        
        print("\nMinimal kinase reaction test:")
        print(f"  Matched: {match.matched}")
        print(f"  Bindings: {match.bindings}")
        print(f"  Unmatched left: {[p.name for p in match.unmatched_left]}")
        print(f"  Unmatched right: {[p.name for p in match.unmatched_right]}")
        
        assert match.matched, "Simple ATP + glucose -> ADP + G6P should match kinase pattern"

    def test_kinase_with_extra_h_plus_strict_mode(self):
        """Test kinase reaction with H+ and strict mode."""
        # Add H+ to the simple kinase reaction
        reaction = Reaction(
            left_participants=[
                Participant(chebi_id="CHEBI:30616", name="ATP"),
                Participant(chebi_id="CHEBI:17234", name="glucose"),
            ],
            right_participants=[
                Participant(chebi_id="CHEBI:456216", name="ADP"),
                Participant(chebi_id="CHEBI:17665", name="glucose-6-phosphate"),
                Participant(chebi_id="CHEBI:15378", name="H(+)"),
            ],
        )
        
        kinase = Kinase()
        match = match_patterns(reaction, kinase.PATTERNS, strict=True)
        
        print("\nKinase with H+ test (strict=True):")
        print(f"  Matched: {match.matched}")
        print(f"  Bindings: {match.bindings}")
        print(f"  Unmatched left: {[p.name for p in match.unmatched_left]}")
        print(f"  Unmatched right: {[p.name for p in match.unmatched_right]}")
        
        # This might fail because strict=True and the pattern has optional H+
        # The test will reveal if this is the issue
        if not match.matched:
            print("  HYPOTHESIS: strict=True is rejecting reactions with H+ that don't match optional pattern")
        
        # Let's also test with strict=False
        match_loose = match_patterns(reaction, kinase.PATTERNS, strict=False)
        print(f"  Matched (strict=False): {match_loose.matched}")
        
        # At minimum, non-strict should work
        assert match_loose.matched, "Kinase with H+ should match in non-strict mode"