"""Focused tests for donor-specific GO classes with support above 10."""

from __future__ import annotations

from typing_extensions import TypedDict

from autarch.datamodel import Participant, Reaction
from autarch.molecules import CHEBI_GDP, CHEBI_UDP
from autarch.ontology.acetylglucosaminyltransferase_activity import (
    AcetylglucosaminyltransferaseActivity,
)
from autarch.ontology.antiporter_activity import AntiporterActivity
from autarch.ontology.galactosyltransferase_activity import (
    GalactosyltransferaseActivity,
)
from autarch.ontology.glucosyltransferase_activity import (
    GlucosyltransferaseActivity,
)
from autarch.ontology.mannosyltransferase_activity import (
    MannosyltransferaseActivity,
)
from autarch.ontology.secondary_active_transmembrane_transporter_activity import (
    SecondaryActiveTransmembraneTransporterActivity,
)
from autarch.ontology.symporter_activity import SymporterActivity
from autarch.ontology.udp_galactosyltransferase_activity import (
    UDPGalactosyltransferaseActivity,
)
from autarch.ontology.udp_glucosyltransferase_activity import (
    UDPGlucosyltransferaseActivity,
)
from autarch.ontology.udp_glycosyltransferase_activity import (
    UDPGlycosyltransferaseActivity,
)

CHEBI_UDP_GLUCOSE = "CHEBI:18066"
CHEBI_UDP_GALACTOSE = "CHEBI:18307"
CHEBI_UDP_GLCNAC = "CHEBI:57705"
CHEBI_GDP_MANNOSE = "CHEBI:57527"
CHEBI_ADP = "CHEBI:456216"
CHEBI_ATP = "CHEBI:30616"
CHEBI_H2O = "CHEBI:15377"
CHEBI_H_PLUS = "CHEBI:15378"
CHEBI_DIPHOSPHATE = "CHEBI:33019"
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


def _glycosyl_transfer(donor_chebi: str, product_chebi: str) -> Reaction:
    return _reaction(
        left=[
            {"chebi_id": donor_chebi, "name": "sugar donor"},
            {"smiles": "Oc1ccccc1", "name": "phenol"},
        ],
        right=[
            {"chebi_id": product_chebi, "name": "released nucleotide"},
            {
                "smiles": "OC[C@H]1O[C@H](Oc2ccccc2)[C@@H](O)[C@H](O)[C@H]1O",
                "name": "aryl glycoside",
            },
        ],
    )


def _symport_reaction() -> Reaction:
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


def _antiport_reaction() -> Reaction:
    return _reaction(
        left=[
            {
                "chebi_id": CHEBI_METHIONINE,
                "smiles": "CSCC[C@H]([NH3+])C(=O)[O-]",
                "name": "L-methionine",
                "location": "out",
            },
            {"chebi_id": CHEBI_SODIUM, "name": "sodium ion", "location": "in"},
        ],
        right=[
            {
                "chebi_id": CHEBI_METHIONINE,
                "smiles": "CSCC[C@H]([NH3+])C(=O)[O-]",
                "name": "L-methionine",
                "location": "in",
            },
            {"chebi_id": CHEBI_SODIUM, "name": "sodium ion", "location": "out"},
        ],
    )


def _primary_active_transport_reaction() -> Reaction:
    return _reaction(
        left=[
            {
                "chebi_id": CHEBI_METHIONINE,
                "smiles": "CSCC[C@H]([NH3+])C(=O)[O-]",
                "name": "L-methionine",
                "location": "out",
            },
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
            {"chebi_id": CHEBI_H2O, "name": "water"},
        ],
        right=[
            {
                "chebi_id": CHEBI_METHIONINE,
                "smiles": "CSCC[C@H]([NH3+])C(=O)[O-]",
                "name": "L-methionine",
                "location": "in",
            },
            {"chebi_id": CHEBI_ADP, "name": "ADP"},
            {"chebi_id": CHEBI_DIPHOSPHATE, "name": "diphosphate"},
        ],
    )


def test_udp_glycosyltransferase_activity_positive() -> None:
    assert UDPGlycosyltransferaseActivity().check_membership_impl(
        _glycosyl_transfer(CHEBI_UDP_GLUCOSE, CHEBI_UDP)
    ).is_member is True


def test_glucosyltransferase_activity_positive() -> None:
    assert GlucosyltransferaseActivity().check_membership_impl(
        _glycosyl_transfer(CHEBI_UDP_GLUCOSE, CHEBI_UDP)
    ).is_member is True


def test_udp_glucosyltransferase_activity_rejects_galactose_donor() -> None:
    cls = UDPGlucosyltransferaseActivity()
    assert cls.check_membership_impl(_glycosyl_transfer(CHEBI_UDP_GLUCOSE, CHEBI_UDP)).is_member is True
    assert cls.check_membership_impl(_glycosyl_transfer(CHEBI_UDP_GALACTOSE, CHEBI_UDP)).is_member is False


def test_galactosyltransferase_activity_positive() -> None:
    assert GalactosyltransferaseActivity().check_membership_impl(
        _glycosyl_transfer(CHEBI_UDP_GALACTOSE, CHEBI_UDP)
    ).is_member is True


def test_udp_galactosyltransferase_activity_rejects_glucose_donor() -> None:
    cls = UDPGalactosyltransferaseActivity()
    assert cls.check_membership_impl(_glycosyl_transfer(CHEBI_UDP_GALACTOSE, CHEBI_UDP)).is_member is True
    assert cls.check_membership_impl(_glycosyl_transfer(CHEBI_UDP_GLUCOSE, CHEBI_UDP)).is_member is False


def test_acetylglucosaminyltransferase_activity_positive() -> None:
    assert AcetylglucosaminyltransferaseActivity().check_membership_impl(
        _glycosyl_transfer(CHEBI_UDP_GLCNAC, CHEBI_UDP)
    ).is_member is True


def test_mannosyltransferase_activity_positive() -> None:
    cls = MannosyltransferaseActivity()
    assert cls.check_membership_impl(_glycosyl_transfer(CHEBI_GDP_MANNOSE, CHEBI_GDP)).is_member is True
    assert cls.check_membership_impl(_glycosyl_transfer(CHEBI_UDP_GLUCOSE, CHEBI_UDP)).is_member is False


def test_symporter_activity_positive() -> None:
    assert SymporterActivity().check_membership_impl(_symport_reaction()).is_member is True


def test_antiporter_activity_positive() -> None:
    assert AntiporterActivity().check_membership_impl(_antiport_reaction()).is_member is True


def test_secondary_active_transmembrane_transporter_activity() -> None:
    cls = SecondaryActiveTransmembraneTransporterActivity()
    assert cls.check_membership_impl(_symport_reaction()).is_member is True
    assert cls.check_membership_impl(_antiport_reaction()).is_member is True
    assert cls.check_membership_impl(_primary_active_transport_reaction()).is_member is False
