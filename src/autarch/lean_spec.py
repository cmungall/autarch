"""Lean-oriented symbolic IR bridge.

This module exports the symbolic subset of autarch reactions and pattern DSL
terms into a small JSON-friendly representation that can be consumed by a Lean
reference implementation. It intentionally ignores RDKit state and focuses on
stable identifiers plus simple multiplicity constraints.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from autarch.datamodel import Participant, Reaction
from autarch.ontology.kinase import Kinase
from autarch.pattern_dsl import PatternParticipant

JsonDict = dict[str, Any]


def participant_key_to_lean_ir(participant: Participant) -> JsonDict:
    """Serialize the best available participant identity key.

    The precedence matches the current Python matcher semantics: CHEBI ID,
    then SMILES, then name.
    """
    if participant.chebi_id:
        return {"kind": "chebi", "value": participant.chebi_id}
    if participant.smiles:
        return {"kind": "smiles", "value": participant.smiles}
    if participant.name:
        return {"kind": "name", "value": participant.name}
    return {"kind": "anonymous", "value": str(participant)}


def participant_to_lean_ir(participant: Participant) -> JsonDict:
    """Serialize a concrete participant for Lean."""
    return {
        "key": participant_key_to_lean_ir(participant),
        "count": participant.count,
        "location": participant.location,
        "stoichiometry": participant.stoichiometry,
        "polymer_index": participant.polymer_index,
    }


def pattern_term_to_lean_ir(participant: Participant) -> JsonDict:
    """Serialize a pattern participant or exact participant for Lean."""
    if isinstance(participant, PatternParticipant) and participant.is_variable:
        atom: JsonDict = {
            "kind": "var",
            "name": participant.variable,
        }
    else:
        atom = {
            "kind": "exact",
            "key": participant_key_to_lean_ir(participant),
        }

    min_count = 1
    max_count: int | None = 1
    if isinstance(participant, PatternParticipant):
        min_count = participant.min_count
        max_count = participant.max_count

    return {
        "atom": atom,
        "count": {
            "lo": min_count,
            "hi": max_count,
        },
        "location": participant.location,
    }


def reaction_to_lean_ir(reaction: Reaction) -> JsonDict:
    """Serialize a concrete reaction for Lean."""
    return {
        "lhs": [participant_to_lean_ir(p) for p in reaction.left_participants],
        "rhs": [participant_to_lean_ir(p) for p in reaction.right_participants],
        "label": reaction.label,
        "transport": reaction.is_transport_reaction(),
    }


def pattern_to_lean_ir(pattern: Reaction) -> JsonDict:
    """Serialize a pattern reaction for Lean."""
    return {
        "lhs": [pattern_term_to_lean_ir(p) for p in pattern.left_participants],
        "rhs": [pattern_term_to_lean_ir(p) for p in pattern.right_participants],
        "label": pattern.label,
    }


def match_problem_to_lean_ir(
    reaction: Reaction,
    pattern: Reaction,
    *,
    strict: bool,
) -> JsonDict:
    """Bundle a reaction and pattern into one JSON-friendly match problem."""
    return {
        "strict": strict,
        "reaction": reaction_to_lean_ir(reaction),
        "pattern": pattern_to_lean_ir(pattern),
    }


def kinase_patterns_to_lean_ir() -> list[JsonDict]:
    """Serialize the current kinase patterns for Lean."""
    return [pattern_to_lean_ir(pattern) for pattern in Kinase.PATTERNS]


def kinase_reference_bundle() -> JsonDict:
    """Build a small export bundle for the current kinase classifier."""
    return {
        "classifier": "Kinase",
        "go_id": Kinase.GO_ID,
        "ec_number_prefix": Kinase.EC_NUMBER_PREFIX,
        "patterns": kinase_patterns_to_lean_ir(),
    }


def write_lean_ir_json(data: Any, path: str | Path) -> None:
    """Write Lean IR data to disk as formatted JSON."""
    target = Path(path)
    target.write_text(
        json.dumps(data, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


__all__ = [
    "participant_key_to_lean_ir",
    "participant_to_lean_ir",
    "pattern_term_to_lean_ir",
    "reaction_to_lean_ir",
    "pattern_to_lean_ir",
    "match_problem_to_lean_ir",
    "kinase_patterns_to_lean_ir",
    "kinase_reference_bundle",
    "write_lean_ir_json",
]
