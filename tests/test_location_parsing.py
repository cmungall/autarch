"""Test location/compartment parsing from RHEA labels."""

from autarch.etl.rhea_etl import parse_location_from_label
from autarch.datamodel import Reaction, Participant


class TestLocationParsing:
    """Test parsing of cellular locations from RHEA labels."""

    def test_parse_simple_in_out_locations(self):
        """Test parsing simple in/out locations."""
        label = "ATP(in) + ADP(out) = ATP(out) + ADP(in)"
        locations = parse_location_from_label(label)

        assert "ATP" in locations
        assert locations["ATP"] == "out"  # Last occurrence wins
        assert "ADP" in locations
        assert locations["ADP"] == "in"  # Last occurrence wins

    def test_parse_periplasm_cytoplasm(self):
        """Test parsing specific compartments."""
        label = "H(+)(periplasm) + glucose(extracellular) = H(+)(cytoplasm) + glucose(cytoplasm)"
        locations = parse_location_from_label(label)

        assert "H(+)" in locations
        assert locations["H(+)"] == "cytoplasm"  # Last occurrence
        assert "glucose" in locations
        assert locations["glucose"] == "cytoplasm"  # Last occurrence

    def test_parse_no_locations(self):
        """Test label without location information."""
        label = "ATP + H2O = ADP + phosphate"
        locations = parse_location_from_label(label)

        assert len(locations) == 0

    def test_parse_mixed_format(self):
        """Test mixed format with some molecules having locations."""
        label = "ATP(mitochondria) + glucose = ADP(mitochondria) + glucose-6-phosphate"
        locations = parse_location_from_label(label)

        assert "ATP" in locations
        assert locations["ATP"] == "mitochondria"
        assert "ADP" in locations
        assert locations["ADP"] == "mitochondria"
        assert "glucose" not in locations
        assert "glucose-6-phosphate" not in locations

    def test_transport_reaction_detection(self):
        """Test that transport reactions are correctly identified."""
        # ATP/ADP translocase
        reaction = Reaction(
            left_participants=[
                Participant(chebi_id="CHEBI:15422", location="in"),  # ATP
                Participant(chebi_id="CHEBI:16761", location="out"),  # ADP
            ],
            right_participants=[
                Participant(chebi_id="CHEBI:15422", location="out"),  # ATP
                Participant(chebi_id="CHEBI:16761", location="in"),  # ADP
            ],
        )

        assert reaction.is_transport_reaction()

        transported = reaction.get_transported_molecules()
        assert len(transported) == 2

        # Check ATP transport
        atp_transport = [t for t in transported if "15422" in t[0]]
        assert len(atp_transport) == 1
        assert atp_transport[0][1] == "in"
        assert atp_transport[0][2] == "out"

        # Check ADP transport
        adp_transport = [t for t in transported if "16761" in t[0]]
        assert len(adp_transport) == 1
        assert adp_transport[0][1] == "out"
        assert adp_transport[0][2] == "in"

    def test_non_transport_reaction(self):
        """Test that regular reactions are not identified as transport."""
        # ATP hydrolysis - not a transport reaction
        reaction = Reaction(
            left_participants=[
                Participant(chebi_id="CHEBI:15422"),  # ATP
                Participant(chebi_id="CHEBI:15377"),  # water
            ],
            right_participants=[
                Participant(chebi_id="CHEBI:16761"),  # ADP
                Participant(chebi_id="CHEBI:43474"),  # phosphate
            ],
        )

        assert not reaction.is_transport_reaction()

    def test_same_location_not_transport(self):
        """Test that molecules in same location are not transport."""
        reaction = Reaction(
            left_participants=[
                Participant(chebi_id="CHEBI:15422", location="cytoplasm"),  # ATP
                Participant(chebi_id="CHEBI:15377", location="cytoplasm"),  # water
            ],
            right_participants=[
                Participant(chebi_id="CHEBI:16761", location="cytoplasm"),  # ADP
                Participant(chebi_id="CHEBI:43474", location="cytoplasm"),  # phosphate
            ],
        )

        assert not reaction.is_transport_reaction()
