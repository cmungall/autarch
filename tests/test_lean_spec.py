"""Tests for the Lean scaffold bridge."""

import json

from autarch.datamodel import Participant, Reaction
from autarch.lean_spec import (
    kinase_patterns_to_lean_ir,
    kinase_reference_bundle,
    match_problem_to_lean_ir,
    pattern_to_lean_ir,
    reaction_to_lean_ir,
    write_lean_ir_json,
)
from autarch.molecules import adp, atp, h_plus
from autarch.pattern_dsl import optional, pattern, var


def test_reaction_to_lean_ir_serializes_symbolic_fields() -> None:
    """Concrete reactions should serialize to stable keys and metadata."""
    reaction = Reaction(
        left_participants=[
            Participant(
                chebi_id="CHEBI:30616",
                name="ATP",
                count=2,
                location="cytoplasm",
            ),
            Participant(name="glucose"),
        ],
        right_participants=[
            Participant(smiles="O", name="water"),
        ],
        label="example",
    )

    ir = reaction_to_lean_ir(reaction)

    assert ir["lhs"][0]["key"] == {"kind": "chebi", "value": "CHEBI:30616"}
    assert ir["lhs"][0]["count"] == 2
    assert ir["lhs"][0]["location"] == "cytoplasm"
    assert ir["lhs"][1]["key"] == {"kind": "name", "value": "glucose"}
    assert ir["rhs"][0]["key"] == {"kind": "smiles", "value": "O"}
    assert ir["label"] == "example"
    assert ir["transport"] is False


def test_pattern_to_lean_ir_preserves_variables_and_optional_bounds() -> None:
    """Pattern serialization should preserve variables and count bounds."""
    kinase_like = pattern(
        [atp, var("substrate")],
        [adp, var("product"), optional(h_plus)],
    )

    ir = pattern_to_lean_ir(kinase_like)

    assert ir["lhs"][0]["atom"] == {
        "kind": "exact",
        "key": {"kind": "chebi", "value": "CHEBI:30616"},
    }
    assert ir["lhs"][1]["atom"] == {"kind": "var", "name": "substrate"}
    assert ir["rhs"][1]["atom"] == {"kind": "var", "name": "product"}
    assert ir["rhs"][2]["count"] == {"lo": 0, "hi": 1}
    assert ir["rhs"][2]["atom"] == {
        "kind": "exact",
        "key": {"kind": "chebi", "value": "CHEBI:15378"},
    }


def test_match_problem_to_lean_ir_wraps_reaction_and_pattern() -> None:
    """A match problem should contain strictness plus serialized payloads."""
    reaction = Reaction(
        left_participants=[Participant(chebi_id="CHEBI:30616", name="ATP")],
        right_participants=[Participant(chebi_id="CHEBI:456216", name="ADP")],
    )
    pat = pattern([atp], [adp])

    payload = match_problem_to_lean_ir(reaction, pat, strict=True)

    assert payload["strict"] is True
    assert payload["reaction"]["lhs"][0]["key"]["value"] == "CHEBI:30616"
    assert payload["pattern"]["rhs"][0]["atom"]["key"]["value"] == "CHEBI:456216"


def test_kinase_patterns_to_lean_ir_exports_both_reference_patterns() -> None:
    """Kinase export should stay aligned with the current Python classifier."""
    patterns = kinase_patterns_to_lean_ir()

    assert len(patterns) == 2
    assert patterns[0]["lhs"][0]["atom"]["key"]["value"] == "CHEBI:30616"
    assert patterns[0]["rhs"][0]["atom"]["key"]["value"] == "CHEBI:456216"
    assert patterns[1]["lhs"][0]["atom"]["key"]["value"] == "CHEBI:37565"
    assert patterns[1]["rhs"][0]["atom"]["key"]["value"] == "CHEBI:58189"


def test_kinase_reference_bundle_contains_classifier_metadata() -> None:
    """The export bundle should include classifier metadata plus patterns."""
    bundle = kinase_reference_bundle()

    assert bundle["classifier"] == "Kinase"
    assert bundle["go_id"] == "GO:0016301"
    assert bundle["ec_number_prefix"] == "2.7.-.-"
    assert len(bundle["patterns"]) == 2


def test_write_lean_ir_json_round_trips(tmp_path) -> None:
    """Lean IR JSON output should be stable and machine-readable."""
    destination = tmp_path / "kinase.json"

    write_lean_ir_json(kinase_reference_bundle(), destination)

    payload = json.loads(destination.read_text())
    assert payload["classifier"] == "Kinase"
    assert payload["patterns"][0]["lhs"][0]["atom"]["key"]["value"] == "CHEBI:30616"
