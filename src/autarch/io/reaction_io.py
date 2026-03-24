"""IO operations for chemical reactions."""

import json
from pathlib import Path
from typing import Any, Dict

import yaml

from autarch.datamodel import Reaction, Participant


def load_reaction_from_yaml(file_path: Path | str) -> Reaction:
    """Load a reaction from a YAML file.

    Args:
        file_path: Path to the YAML file

    Returns:
        Reaction object

    Expected YAML format:
    ```yaml
    left_participants:
      - smiles: "CC(=O)OCC"
        count: 1
      - smiles: "O"
        count: 1
    right_participants:
      - smiles: "CC(=O)O"
        count: 1
      - smiles: "CCO"
        count: 1
    ```
    """
    with open(file_path, "r") as f:
        data = yaml.safe_load(f)

    return _dict_to_reaction(data)


def save_reaction_to_yaml(reaction: Reaction, file_path: Path | str) -> None:
    """Save a reaction to a YAML file.

    Args:
        reaction: Reaction object to save
        file_path: Path to save the YAML file
    """
    data = _reaction_to_dict(reaction)

    with open(file_path, "w") as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)


def load_reaction_from_json(file_path: Path | str) -> Reaction:
    """Load a reaction from a JSON file.

    Args:
        file_path: Path to the JSON file

    Returns:
        Reaction object
    """
    with open(file_path, "r") as f:
        data = json.load(f)

    return _dict_to_reaction(data)


def save_reaction_to_json(reaction: Reaction, file_path: Path | str) -> None:
    """Save a reaction to a JSON file.

    Args:
        reaction: Reaction object to save
        file_path: Path to save the JSON file
    """
    data = _reaction_to_dict(reaction)

    with open(file_path, "w") as f:
        json.dump(data, f, indent=2)


def _dict_to_reaction(data: Dict[str, Any]) -> Reaction:
    """Convert a dictionary to a Reaction object.

    Args:
        data: Dictionary containing reaction data

    Returns:
        Reaction object
    """
    left = [
        Participant(
            smiles=p.get("smiles"),
            chebi_id=p.get("chebi_id"),
            name=p.get("name"),
            count=p.get("count", 1)
        )
        for p in data.get("left_participants", [])
    ]

    right = [
        Participant(
            smiles=p.get("smiles"),
            chebi_id=p.get("chebi_id"), 
            name=p.get("name"),
            count=p.get("count", 1)
        )
        for p in data.get("right_participants", [])
    ]

    return Reaction(left_participants=left, right_participants=right)


def _reaction_to_dict(reaction: Reaction) -> Dict[str, Any]:
    """Convert a Reaction object to a dictionary.

    Args:
        reaction: Reaction object

    Returns:
        Dictionary representation of the reaction
    """
    def participant_to_dict(p):
        d = {}
        if p.smiles:
            d["smiles"] = p.smiles
        if p.chebi_id:
            d["chebi_id"] = p.chebi_id
        if p.name:
            d["name"] = p.name
        d["count"] = p.count
        return d
    
    return {
        "left_participants": [
            participant_to_dict(stoi)
            for stoi in reaction.left_participants
        ],
        "right_participants": [
            participant_to_dict(stoi)
            for stoi in reaction.right_participants
        ],
    }
