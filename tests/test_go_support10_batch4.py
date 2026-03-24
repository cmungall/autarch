"""Focused tests for transport GO classes with support above 10."""

from __future__ import annotations

from typing_extensions import TypedDict

from autarch.datamodel import Participant, Reaction
from autarch.ontology.carboxylic_acid_transmembrane_transporter_activity import (
    CarboxylicAcidTransmembraneTransporterActivity,
)
from autarch.ontology.metal_ion_transmembrane_transporter_activity import (
    MetalIonTransmembraneTransporterActivity,
)
from autarch.ontology.monosaccharide_transmembrane_transporter_activity import (
    MonosaccharideTransmembraneTransporterActivity,
)
from autarch.ontology.sodium_ion_transmembrane_transporter_activity import (
    SodiumIonTransmembraneTransporterActivity,
)
from autarch.ontology.solute_monoatomic_cation_symporter_activity import (
    SoluteMonoatomicCationSymporterActivity,
)
from autarch.ontology.solute_sodium_symporter_activity import (
    SoluteSodiumSymporterActivity,
)

CHEBI_GLUCOSE = "CHEBI:4167"
CHEBI_METHIONINE = "CHEBI:57844"
CHEBI_SODIUM = "CHEBI:29101"


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


def _sodium_transport() -> Reaction:
    return _reaction(
        left=[{"chebi_id": CHEBI_SODIUM, "name": "sodium ion", "location": "out"}],
        right=[{"chebi_id": CHEBI_SODIUM, "name": "sodium ion", "location": "in"}],
    )


def _sodium_symport() -> Reaction:
    return _reaction(
        left=[
            {
                "chebi_id": CHEBI_METHIONINE,
                "smiles": "CSCC[C@H]([NH3+])C(=O)[O-]",
                "name": "L-methionine",
                "location": "out",
            },
            {"chebi_id": CHEBI_SODIUM, "name": "sodium ion", "location": "out"},
        ],
        right=[
            {
                "chebi_id": CHEBI_METHIONINE,
                "smiles": "CSCC[C@H]([NH3+])C(=O)[O-]",
                "name": "L-methionine",
                "location": "in",
            },
            {"chebi_id": CHEBI_SODIUM, "name": "sodium ion", "location": "in"},
        ],
    )


def _succinate_transport() -> Reaction:
    return _reaction(
        left=[{"smiles": "O=C([O-])CCC(=O)[O-]", "name": "succinate", "location": "out"}],
        right=[{"smiles": "O=C([O-])CCC(=O)[O-]", "name": "succinate", "location": "in"}],
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


def test_metal_ion_transmembrane_transporter_activity_positive() -> None:
    assert MetalIonTransmembraneTransporterActivity().check_membership_impl(_sodium_transport()).is_member is True


def test_sodium_ion_transmembrane_transporter_activity_positive() -> None:
    assert SodiumIonTransmembraneTransporterActivity().check_membership_impl(_sodium_transport()).is_member is True


def test_carboxylic_acid_transmembrane_transporter_activity_positive() -> None:
    assert CarboxylicAcidTransmembraneTransporterActivity().check_membership_impl(_succinate_transport()).is_member is True


def test_solute_monoatomic_cation_symporter_activity_positive() -> None:
    assert SoluteMonoatomicCationSymporterActivity().check_membership_impl(_sodium_symport()).is_member is True


def test_solute_sodium_symporter_activity_positive() -> None:
    assert SoluteSodiumSymporterActivity().check_membership_impl(_sodium_symport()).is_member is True


def test_monosaccharide_transmembrane_transporter_activity_positive() -> None:
    assert MonosaccharideTransmembraneTransporterActivity().check_membership_impl(_glucose_transport()).is_member is True
