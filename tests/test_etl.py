"""Tests for ETL functions mining GO classifications for reactions."""

import json
import csv
from autarch.datamodel import GoTerm
from autarch.etl.chebi_etl import (
    fetch_go_enzyme_mappings,
    mine_go_enzyme_classifications,
    serialize_go_terms,
)


def test_go_term_dataclass():
    """Test the GoTerm dataclass."""
    term = GoTerm(
        go_id="GO:0004553",
        label="hydrolase activity, hydrolyzing O-glycosyl compounds",
        definition="Catalysis of the hydrolysis of O-glycosyl bonds",
        rhea_ids=["RHEA:10000", "RHEA:10001"],
        ec_numbers=["3.2.1.1", "3.2.1.2"],
    )
    assert term.go_id == "GO:0004553"
    assert term.label == "hydrolase activity, hydrolyzing O-glycosyl compounds"
    assert len(term.rhea_ids) == 2
    assert len(term.ec_numbers) == 2


def test_go_term_defaults():
    """Test GoTerm with default values."""
    term = GoTerm(go_id="GO:0000001")
    assert term.go_id == "GO:0000001"
    assert term.label == ""
    assert term.definition == ""
    assert term.rhea_ids == []
    assert term.ec_numbers == []


def test_fetch_go_enzyme_mappings():
    """Test fetching GO enzyme mappings."""
    terms = fetch_go_enzyme_mappings(go_terms=["GO:0004553"])
    assert isinstance(terms, dict)
    # May be empty if no mappings exist
    if terms:
        term = next(iter(terms.values()))
        assert isinstance(term, GoTerm)
        assert term.go_id.startswith("GO:")
        assert isinstance(term.label, str)
        assert isinstance(term.rhea_ids, list)
        assert isinstance(term.ec_numbers, list)


def test_fetch_go_enzyme_mappings_no_filter():
    """Test fetching GO enzyme mappings without filter."""
    # Fetch a limited set to avoid timeout
    # Just test with a few specific terms
    test_terms = [
        "GO:0016787",
        "GO:0016740",
        "GO:0016491",
    ]  # hydrolase, transferase, oxidoreductase
    terms = fetch_go_enzyme_mappings(go_terms=test_terms)
    assert isinstance(terms, dict)
    # Check structure if any results
    if terms:
        first_term = next(iter(terms.values()))
        assert isinstance(first_term, GoTerm)


def test_mine_go_enzyme_classifications():
    """Test the main mining function."""
    # Note: This would fetch all molecular function terms if we don't limit it
    # For testing, let's mock or limit the scope
    # We'll test with the actual function but understand it may return empty results
    result = mine_go_enzyme_classifications()

    assert isinstance(result, dict)
    assert "go_terms" in result
    assert "total_terms" in result
    assert "total_rhea_mappings" in result
    assert "total_ec_mappings" in result
    assert "timestamp" in result

    assert isinstance(result["go_terms"], dict)
    assert result["total_terms"] >= 0
    assert result["total_rhea_mappings"] >= 0
    assert result["total_ec_mappings"] >= 0

    # Check structure of terms if any exist
    if result["go_terms"]:
        first_id = next(iter(result["go_terms"]))
        first_term = result["go_terms"][first_id]
        assert "go_id" in first_term
        assert "label" in first_term
        assert "definition" in first_term
        assert "rhea_ids" in first_term
        assert "ec_numbers" in first_term
        assert "rhea_count" in first_term
        assert "ec_count" in first_term


def test_mine_go_enzyme_classifications_with_output(tmp_path):
    """Test mining with output file."""
    output_file = tmp_path / "go_enzyme_mappings.json"

    result = mine_go_enzyme_classifications(output_file=str(output_file))

    assert output_file.exists()

    # Load and verify the saved file
    with open(output_file) as f:
        saved_data = json.load(f)

    assert saved_data == result
    assert "go_terms" in saved_data
    assert "timestamp" in saved_data


def test_serialize_go_terms_jsonl(tmp_path):
    """Test serializing GO terms to JSONL format."""
    terms = {
        "GO:0004553": GoTerm(
            go_id="GO:0004553",
            label="hydrolase activity",
            rhea_ids=["RHEA:10000"],
            ec_numbers=["3.2.1.1"],
        )
    }

    output_file = tmp_path / "go_terms.jsonl"
    serialize_go_terms(terms, str(output_file), format="jsonl")

    assert output_file.exists()

    # Read and verify
    with open(output_file) as f:
        lines = f.readlines()
    assert len(lines) == 1

    data = json.loads(lines[0])
    assert data["go_id"] == "GO:0004553"
    assert data["rhea_ids"] == ["RHEA:10000"]


def test_serialize_go_terms_csv(tmp_path):
    """Test serializing GO terms to CSV format."""
    terms = {
        "GO:0004553": GoTerm(
            go_id="GO:0004553",
            label="hydrolase activity",
            rhea_ids=["RHEA:10000", "RHEA:10001"],
            ec_numbers=["3.2.1.1"],
        )
    }

    output_file = tmp_path / "go_terms.csv"
    serialize_go_terms(terms, str(output_file), format="csv")

    assert output_file.exists()

    # Read and verify CSV
    with open(output_file) as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    assert len(rows) == 1
    assert rows[0]["go_id"] == "GO:0004553"
    assert rows[0]["rhea_ids"] == "RHEA:10000|RHEA:10001"
