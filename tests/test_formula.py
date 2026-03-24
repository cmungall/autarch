"""Test the chemical reaction formula DSL."""

from autarch.formula import (
    ParticipantList,
    ReactionBuilder,
    molecule,
    ATP,
    ADP,
    Pi,
    H2O,
    H,
    NAD,
    NADH,
    glucose,
    O2,
    CO2,
)
from autarch.datamodel import Reaction


class TestParticipantBuilder:
    """Test the ParticipantBuilder class and its operators."""

    def test_molecule_creation(self):
        """Test creating molecules."""
        mol = molecule("CHEBI:12345", name="test")
        assert mol.chebi_id == "CHEBI:12345"
        assert mol.name == "test"
        assert mol.count == 1
        assert mol.location is None

    def test_stoichiometry_multiplication(self):
        """Test stoichiometry with * operator."""
        # Right multiplication
        mol2 = ATP * 2
        assert mol2.count == 2
        assert mol2.chebi_id == ATP.chebi_id

        # Left multiplication
        mol3 = 3 * H2O
        assert mol3.count == 3
        assert mol3.chebi_id == H2O.chebi_id

        # Chained multiplication
        mol4 = 2 * (ATP * 2)
        assert mol4.count == 4

    def test_location_subscript(self):
        """Test setting location with [] operator."""
        atp_in = ATP["in"]
        assert atp_in.location == "in"
        assert atp_in.chebi_id == ATP.chebi_id

        atp_cytoplasm = ATP["cytoplasm"]
        assert atp_cytoplasm.location == "cytoplasm"

        # With stoichiometry and location
        atp2_mito = (ATP * 2)["mitochondria"]
        assert atp2_mito.count == 2
        assert atp2_mito.location == "mitochondria"

    def test_molecule_addition(self):
        """Test adding molecules with + operator."""
        # Two molecules
        participants = ATP + H2O
        assert isinstance(participants, ParticipantList)
        assert len(participants) == 2
        assert participants.molecules[0].name == "ATP"
        assert participants.molecules[1].name == "H2O"

        # Three molecules
        participants = ATP + H2O + Pi
        assert len(participants) == 3

        # With stoichiometry
        participants = (ATP * 2) + H2O
        assert participants.molecules[0].count == 2


class TestParticipantList:
    """Test the ParticipantList class."""

    def test_participant_list_creation(self):
        """Test creating and using participant lists."""
        plist = ParticipantList([ATP, H2O])
        assert len(plist) == 2

        # Should be iterable
        molecules = list(plist)
        assert molecules[0] == ATP
        assert molecules[1] == H2O

    def test_participant_list_addition(self):
        """Test adding to participant lists."""
        plist1 = ATP + H2O
        plist2 = plist1 + Pi
        assert len(plist2) == 3

        # Adding two lists
        left = ATP + H2O
        right = ADP + Pi
        combined = left + right
        assert len(combined) == 4

    def test_to_participants(self):
        """Test converting to Participant objects."""
        plist = (ATP * 2) + H2O["cytoplasm"]
        participants = plist.to_participants()

        assert len(participants) == 2
        assert participants[0].chebi_id == ATP.chebi_id
        assert participants[0].count == 2
        assert participants[1].chebi_id == H2O.chebi_id
        assert participants[1].location == "cytoplasm"


class TestReactionBuilder:
    """Test building reactions with operators."""

    def test_simple_reaction(self):
        """Test creating a simple reaction."""
        # Using >> for left-to-right
        rxn_builder = ATP + H2O >> ADP + Pi
        assert isinstance(rxn_builder, ReactionBuilder)
        assert rxn_builder.direction == "left-to-right"
        assert len(rxn_builder.left) == 2
        assert len(rxn_builder.right) == 2

        # Build to Reaction
        reaction = rxn_builder.build()
        assert isinstance(reaction, Reaction)
        assert len(reaction.left_participants) == 2
        assert len(reaction.right_participants) == 2

    def test_bidirectional_reaction(self):
        """Test bidirectional reaction with | operator."""
        rxn_builder = NAD + H | NADH
        assert rxn_builder.direction == "bidirectional"

        reaction = rxn_builder.build()
        assert len(reaction.left_participants) == 2
        assert len(reaction.right_participants) == 1

    def test_reverse_reaction(self):
        """Test reverse reaction with << operator."""
        rxn_builder = ADP + Pi << ATP + H2O
        assert rxn_builder.direction == "right-to-left"
        assert len(rxn_builder.left) == 2  # ATP + H2O
        assert len(rxn_builder.right) == 2  # ADP + Pi

    def test_complex_reaction(self):
        """Test a complex reaction with stoichiometry."""
        # Cellular respiration
        rxn_builder = glucose + (O2 * 6) >> (CO2 * 6) + (H2O * 6)

        reaction = rxn_builder.build()
        assert len(reaction.left_participants) == 2
        assert len(reaction.right_participants) == 2

        # Check stoichiometry
        assert reaction.left_participants[1].count == 6  # O2
        assert reaction.right_participants[0].count == 6  # CO2
        assert reaction.right_participants[1].count == 6  # H2O

    def test_transport_reaction(self):
        """Test a transport reaction with locations."""
        # ATP/ADP translocase
        rxn_builder = (ATP["out"] + ADP["in"]) >> (ATP["in"] + ADP["out"])

        reaction = rxn_builder.build()
        assert reaction.left_participants[0].location == "out"
        assert reaction.left_participants[1].location == "in"
        assert reaction.right_participants[0].location == "in"
        assert reaction.right_participants[1].location == "out"

        # Should be detected as transport
        assert reaction.is_transport_reaction()

    def test_reaction_string(self):
        """Test string representation of reactions."""
        rxn_builder = (ATP * 2) + H2O >> (ADP * 2) + (Pi * 2)
        rxn_str = str(rxn_builder)

        assert "2 ATP" in rxn_str
        assert "2 ADP" in rxn_str
        assert "→" in rxn_str

        # Bidirectional
        rxn_builder2 = NAD | NADH
        rxn_str2 = str(rxn_builder2)
        assert "⇌" in rxn_str2


class TestIntegration:
    """Test integration with the classifier system."""

    def test_hydrolase_reaction(self):
        """Test creating a hydrolase reaction."""
        # ATP hydrolysis
        rxn_builder = ATP + H2O >> ADP + Pi + H
        reaction = rxn_builder.build()

        # Check it was built correctly
        assert len(reaction.left_participants) == 2
        assert len(reaction.right_participants) == 3
        assert reaction.left_participants[0].chebi_id == "CHEBI:15422"  # ATP
        assert reaction.left_participants[1].chebi_id == "CHEBI:15377"  # H2O

    def test_custom_molecules(self):
        """Test creating custom molecules."""
        # Create custom molecules
        substrate = molecule("CHEBI:99999", name="substrate")
        product = molecule("CHEBI:88888", name="product")

        # Use in reaction
        rxn_builder = substrate + ATP >> product + ADP + Pi
        reaction = rxn_builder.build()

        assert reaction.left_participants[0].chebi_id == "CHEBI:99999"
        assert reaction.left_participants[0].name == "substrate"

    def test_parentheses_grouping(self):
        """Test that parentheses work for grouping."""
        # Should be able to use parentheses for clarity
        rxn_builder = (glucose + O2 * 6) >> (CO2 * 6 + H2O * 6)

        assert len(rxn_builder.left) == 2
        assert len(rxn_builder.right) == 2
        assert rxn_builder.left.molecules[1].count == 6

    def test_complex_stoichiometry(self):
        """Test various stoichiometry patterns."""
        # Different ways to express stoichiometry
        rxn1 = (2 * H + O2) >> (2 * H2O)
        rxn2 = (H * 2 + O2) >> (H2O * 2)

        # Both should give same result
        r1 = rxn1.build()
        r2 = rxn2.build()

        assert r1.left_participants[0].count == 2
        assert r2.left_participants[0].count == 2
        assert r1.right_participants[0].count == 2
        assert r2.right_participants[0].count == 2
