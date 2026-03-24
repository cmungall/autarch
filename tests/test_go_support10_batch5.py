"""Focused tests for additional transport GO classes with support above 10."""

from __future__ import annotations

from typing_extensions import TypedDict

from autarch.datamodel import Participant, Reaction
from autarch.ontology.atpase_coupled_monoatomic_cation_transmembrane_transporter_activity import (
    ATPaseCoupledMonoatomicCationTransmembraneTransporterActivity,
)
from autarch.ontology.metal_cation_monoatomic_cation_antiporter_activity import (
    MetalCationMonoatomicCationAntiporterActivity,
)
from autarch.ontology.solute_proton_symporter_activity import (
    SoluteProtonSymporterActivity,
)
from autarch.ontology.sugar_transmembrane_transporter_activity import (
    SugarTransmembraneTransporterActivity,
)
from autarch.ontology.sulfur_compound_transmembrane_transporter_activity import (
    SulfurCompoundTransmembraneTransporterActivity,
)

CHEBI_ADP = "CHEBI:456216"
CHEBI_ATP = "CHEBI:30616"
CHEBI_GLUCOSE = "CHEBI:4167"
CHEBI_H2O = "CHEBI:15377"
CHEBI_H_PLUS = "CHEBI:15378"
CHEBI_METHIONINE = "CHEBI:57844"
CHEBI_SODIUM = "CHEBI:29101"
CHEBI_FERROUS = "CHEBI:29033"


class ParticipantInput(TypedDict, total=False):
    chebi_id: str
    smiles: str
    name: str
    location: str


def _reaction(left: list[ParticipantInput], right: list[ParticipantInput]) -> Reaction:
    return Reaction(
        left_participants=[Participant(**participant) for participant in left],
        right_participants=[Participant(**participant) for participant in right],
    )


def _atpase_cation_transport() -> Reaction:
    return _reaction(
        left=[
            {"chebi_id": CHEBI_SODIUM, "name": "sodium ion", "location": "out"},
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
            {"chebi_id": CHEBI_H2O, "name": "water"},
        ],
        right=[
            {"chebi_id": CHEBI_SODIUM, "name": "sodium ion", "location": "in"},
            {"chebi_id": CHEBI_ADP, "name": "ADP"},
        ],
    )


def _proton_symport() -> Reaction:
    return _reaction(
        left=[
            {
                "chebi_id": CHEBI_METHIONINE,
                "smiles": "CSCC[C@H]([NH3+])C(=O)[O-]",
                "name": "L-methionine",
                "location": "out",
            },
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron", "location": "out"},
        ],
        right=[
            {
                "chebi_id": CHEBI_METHIONINE,
                "smiles": "CSCC[C@H]([NH3+])C(=O)[O-]",
                "name": "L-methionine",
                "location": "in",
            },
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron", "location": "in"},
        ],
    )


def _glucose_transport() -> Reaction:
    return _reaction(
        left=[
            {
                "chebi_id": CHEBI_GLUCOSE,
                "smiles": "OC[C@H]1O[C@H](O)[C@@H](O)[C@H](O)C1",
                "name": "D-glucose",
                "location": "out",
            }
        ],
        right=[
            {
                "chebi_id": CHEBI_GLUCOSE,
                "smiles": "OC[C@H]1O[C@H](O)[C@@H](O)[C@H](O)C1",
                "name": "D-glucose",
                "location": "in",
            }
        ],
    )


def _sulfur_transport() -> Reaction:
    return _reaction(
        left=[
            {
                "chebi_id": CHEBI_METHIONINE,
                "smiles": "CSCC[C@H]([NH3+])C(=O)[O-]",
                "name": "L-methionine",
                "location": "out",
            }
        ],
        right=[
            {
                "chebi_id": CHEBI_METHIONINE,
                "smiles": "CSCC[C@H]([NH3+])C(=O)[O-]",
                "name": "L-methionine",
                "location": "in",
            }
        ],
    )


def _metal_antiport() -> Reaction:
    return _reaction(
        left=[
            {"chebi_id": CHEBI_FERROUS, "name": "Fe2+", "location": "out"},
            {"chebi_id": CHEBI_SODIUM, "name": "sodium ion", "location": "in"},
        ],
        right=[
            {"chebi_id": CHEBI_FERROUS, "name": "Fe2+", "location": "in"},
            {"chebi_id": CHEBI_SODIUM, "name": "sodium ion", "location": "out"},
        ],
    )


def test_atpase_coupled_monoatomic_cation_transmembrane_transporter_activity_positive() -> None:
    assert ATPaseCoupledMonoatomicCationTransmembraneTransporterActivity().check_membership_impl(
        _atpase_cation_transport()
    ).is_member is True


def test_solute_proton_symporter_activity_positive() -> None:
    assert SoluteProtonSymporterActivity().check_membership_impl(_proton_symport()).is_member is True


def test_sugar_transmembrane_transporter_activity_positive() -> None:
    assert SugarTransmembraneTransporterActivity().check_membership_impl(_glucose_transport()).is_member is True


def test_sulfur_compound_transmembrane_transporter_activity_positive() -> None:
    assert SulfurCompoundTransmembraneTransporterActivity().check_membership_impl(_sulfur_transport()).is_member is True


def test_metal_cation_monoatomic_cation_antiporter_activity_positive() -> None:
    assert MetalCationMonoatomicCationAntiporterActivity().check_membership_impl(_metal_antiport()).is_member is True
