"""Tests for RHEA ETL functions."""

import json
import csv

from autarch.datamodel import RheaTerm, Reaction, Participant
from autarch.etl.rhea_etl import (
    fetch_rhea_reactions,
    serialize_rhea_reactions,
)


def test_rhea_reaction_dataclass():
    """Test the RheaTerm dataclass with Reaction composition."""
    reaction = Reaction(
        left_participants=[
            Participant(chebi_id="CHEBI:83725"),
            Participant(chebi_id="CHEBI:15378"),
        ],
        right_participants=[
            Participant(chebi_id="CHEBI:18310"),
            Participant(chebi_id="CHEBI:16793"),
        ],
    )

    rhea_term = RheaTerm(
        rhea_id="RHEA:18777",
        label="an alkylmercury + H(+) = an alkane + Hg(2+)",
        reaction=reaction,
        ec_numbers=["4.99.1.2"],
        go_terms=["GO:0018836"],
        direction="bidirectional",
    )
    assert rhea_term.rhea_id == "RHEA:18777"
    assert len(rhea_term.reaction.left_participants) == 2
    assert len(rhea_term.reaction.right_participants) == 2
    assert rhea_term.ec_numbers[0] == "4.99.1.2"


def test_rhea_reaction_defaults():
    """Test RheaTerm with default values."""
    rhea_term = RheaTerm(rhea_id="RHEA:00001")
    assert rhea_term.rhea_id == "RHEA:00001"
    assert rhea_term.label == ""
    assert rhea_term.reaction is None
    assert rhea_term.ec_numbers == []
    assert rhea_term.direction == "bidirectional"


def test_fetch_rhea_reactions():
    """Test fetching RHEA reactions."""
    # Test with a specific reaction ID from your example
    reactions = fetch_rhea_reactions(rhea_ids=["RHEA:18777"])
    assert isinstance(reactions, dict)

    if reactions and "RHEA:18777" in reactions:
        rhea_term = reactions["RHEA:18777"]
        assert isinstance(rhea_term, RheaTerm)
        assert rhea_term.rhea_id == "RHEA:18777"
        # The reaction should have some properties based on your example
        assert isinstance(rhea_term.label, str)


def test_serialize_rhea_reactions_jsonl(tmp_path):
    """Test serializing reactions to JSONL format."""
    reaction = Reaction(
        left_participants=[
            Participant(chebi_id="CHEBI:1"),
            Participant(chebi_id="CHEBI:2"),
        ],
        right_participants=[Participant(chebi_id="CHEBI:3")],
    )

    reactions = {
        "RHEA:18777": RheaTerm(
            rhea_id="RHEA:18777", label="test reaction", reaction=reaction
        )
    }

    output_file = tmp_path / "reactions.jsonl"
    serialize_rhea_reactions(reactions, str(output_file), format="jsonl")

    assert output_file.exists()

    # Read and verify the content
    with open(output_file) as f:
        lines = f.readlines()
    assert len(lines) == 1

    data = json.loads(lines[0])
    assert data["rhea_id"] == "RHEA:18777"
    assert data["label"] == "test reaction"
    assert len(data["reaction"]["left_participants"]) == 2
    assert data["reaction"]["left_participants"][0]["chebi_id"] == "CHEBI:1"
    assert data["reaction"]["left_participants"][1]["chebi_id"] == "CHEBI:2"
    assert len(data["reaction"]["right_participants"]) == 1
    assert data["reaction"]["right_participants"][0]["chebi_id"] == "CHEBI:3"


def test_serialize_rhea_reactions_csv(tmp_path):
    """Test serializing reactions to CSV format."""
    reaction = Reaction(
        left_participants=[
            Participant(chebi_id="CHEBI:1"),
            Participant(chebi_id="CHEBI:2"),
        ],
        right_participants=[Participant(chebi_id="CHEBI:3")],
    )

    reactions = {
        "RHEA:18777": RheaTerm(
            rhea_id="RHEA:18777",
            label="test reaction",
            reaction=reaction,
            ec_numbers=["1.1.1.1", "2.2.2.2"],
        )
    }

    output_file = tmp_path / "reactions.csv"
    serialize_rhea_reactions(reactions, str(output_file), format="csv")

    assert output_file.exists()

    # Read and verify the CSV
    with open(output_file) as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    assert len(rows) == 1
    assert rows[0]["rhea_id"] == "RHEA:18777"
    assert rows[0]["label"] == "test reaction"
    assert rows[0]["ec_numbers"] == "1.1.1.1;2.2.2.2"
    assert rows[0]["has_participants"] == "True"
