"""Tests for the CLI interface."""

from pathlib import Path
import json

import pytest
from typer.testing import CliRunner

from autarch.cli import app


runner = CliRunner()
INPUT_DIR = Path(__file__).parent / "input"


class TestCLICommands:
    """Test CLI commands."""

    def test_list_classes_command(self):
        """Test the list-classes command."""
        result = runner.invoke(app, ["list-classes"])
        assert result.exit_code == 0
        assert "Hydrolase" in result.stdout
        assert "Oxidoreductase" in result.stdout
        assert "Available Reaction Classes" in result.stdout

    def test_classify_ester_hydrolysis(self):
        """Test classifying ester hydrolysis reaction."""
        filepath = INPUT_DIR / "ester_hydrolysis.yaml"
        result = runner.invoke(app, ["classify", str(filepath)])

        assert result.exit_code == 0
        assert "Hydrolase" in result.stdout
        assert "✓ Yes" in result.stdout
        # Check that reaction matches something (could be Hydrolase and/or Transferase)
        assert "✓ Reaction matches:" in result.stdout

    def test_classify_with_verbose(self):
        """Test classify command with verbose flag."""
        filepath = INPUT_DIR / "ester_hydrolysis.yaml"
        result = runner.invoke(app, ["classify", str(filepath), "-v"])

        assert result.exit_code == 0
        assert "Hydrolase" in result.stdout
        assert "Explanation" in result.stdout  # Verbose mode shows explanation column
        assert "hydrolysis" in result.stdout.lower() or "water" in result.stdout.lower()

    def test_classify_with_json_format(self):
        """Test classify command with JSON output."""
        filepath = INPUT_DIR / "alcohol_oxidation.yaml"
        result = runner.invoke(app, ["classify", str(filepath), "--format", "json"])

        assert result.exit_code == 0

        # Parse JSON output
        # Extract JSON from output (skip any warnings)
        json_start = result.stdout.find("{")
        json_str = result.stdout[json_start:]
        data = json.loads(json_str)

        # Just check that JSON parsing worked and has expected structure
        assert isinstance(data, dict)
        assert "Oxidoreductase" in data
        assert "is_member" in data["Oxidoreductase"]
        assert "explanation" in data["Oxidoreductase"]

    def test_classify_with_yaml_format(self):
        """Test classify command with YAML output."""
        filepath = INPUT_DIR / "glucose_oxidation.yaml"
        result = runner.invoke(app, ["classify", str(filepath), "--format", "yaml"])

        assert result.exit_code == 0
        assert "Oxidoreductase:" in result.stdout
        # Just check YAML structure is present
        assert "is_member:" in result.stdout
        assert "explanation:" in result.stdout

    def test_classify_specific_class(self):
        """Test classifying with a specific reaction class."""
        filepath = INPUT_DIR / "ester_hydrolysis.yaml"
        result = runner.invoke(app, ["classify", str(filepath), "-c", "Hydrolase"])

        assert result.exit_code == 0
        assert "Hydrolase" in result.stdout
        # Should only show results for Hydrolase, not Oxidoreductase
        assert result.stdout.count("Hydrolase") >= 1
        # Oxidoreductase should not appear in the table (might appear in summary)
        lines = result.stdout.split("\n")
        table_lines = [line for line in lines if "│" in line]
        assert not any("Oxidoreductase" in line for line in table_lines)

    def test_classify_invalid_class(self):
        """Test classifying with an invalid reaction class."""
        filepath = INPUT_DIR / "ester_hydrolysis.yaml"
        result = runner.invoke(app, ["classify", str(filepath), "-c", "InvalidClass"])

        assert result.exit_code == 1
        assert "Unknown reaction class 'InvalidClass'" in result.stdout
        assert "Available classes:" in result.stdout

    def test_classify_nonexistent_file(self):
        """Test classifying with a non-existent file."""
        result = runner.invoke(app, ["classify", "nonexistent.yaml"])

        # Typer validates file existence
        assert result.exit_code != 0

    def test_classify_all_input_files(self):
        """Test that all input files can be classified via CLI."""
        yaml_files = list(INPUT_DIR.glob("*.yaml"))
        assert len(yaml_files) > 0, "No YAML files found"

        for filepath in yaml_files:
            result = runner.invoke(app, ["classify", str(filepath)])
            assert result.exit_code == 0, f"Failed to classify {filepath.name}"
            assert "Reaction Class" in result.stdout
            assert "Is Member" in result.stdout


@pytest.mark.parametrize(
    "filename,expected_contains",
    [
        ("ester_hydrolysis.yaml", "Hydrolase"),
        ("alcohol_oxidation.yaml", "Oxidoreductase"),
        ("glucose_oxidation.yaml", "Oxidoreductase"),
        ("combustion.yaml", "Classification Results"),  # Just check it runs
    ],
)
def test_cli_classification_results(filename, expected_contains):
    """Parametrized test for CLI classification results."""
    filepath = INPUT_DIR / filename
    result = runner.invoke(app, ["classify", str(filepath)])

    assert result.exit_code == 0
    assert expected_contains in result.stdout
