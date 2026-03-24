"""Focused tests for batch 14 specific classifiers."""

from __future__ import annotations

from typing_extensions import TypedDict

from autarch.datamodel import Participant, Reaction
from autarch.molecules import (
    CHEBI_ATP,
    CHEBI_DIPHOSPHATE,
    CHEBI_FAD,
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_NH4,
    CHEBI_O2,
    CHEBI_SAH,
    CHEBI_SAM,
)
from autarch.ontology.ampylase import AMPylase
from autarch.ontology.corrinoid_adenosyltransferase import CorrinoidAdenosyltransferase
from autarch.ontology.nicotinamide_n_methyltransferase import (
    NicotinamideNMethyltransferase,
)
from autarch.ontology.pyridoxamine_phosphate_oxidase import (
    PyridoxaminePhosphateOxidase,
)
from autarch.ontology.spermine_oxidase import SpermineOxidase

CHEBI_REDUCED_ETF = "CHEBI:58307"
CHEBI_TRIPHOSPHATE = "CHEBI:18036"
CHEBI_COB_II_ALAMIN = "CHEBI:16304"
CHEBI_ADENOSYLCOBALAMIN = "CHEBI:18408"
CHEBI_COB_II_INAMIDE = "CHEBI:141013"
CHEBI_ADENOSYLCOBINAMIDE = "CHEBI:2480"

CHEBI_TYROSYL_PROTEIN = "CHEBI:46858"
CHEBI_AMPYL_TYROSYL_PROTEIN = "CHEBI:83624"
CHEBI_FREE_TYROSINE = "CHEBI:17895"

CHEBI_SPERMINE = "CHEBI:45725"
CHEBI_SPERMIDINE = "CHEBI:57834"
CHEBI_NORSPERMINE = "CHEBI:58704"
CHEBI_NORSPERMIDINE = "CHEBI:57920"
CHEBI_AMINOPROPANAL = "CHEBI:133427"
CHEBI_N1_ACETYLSPERMINE = "CHEBI:58101"

CHEBI_NICOTINAMIDE = "CHEBI:17154"
CHEBI_1_METHYLNICOTINAMIDE = "CHEBI:16797"
CHEBI_NICOTINIC_ACID = "CHEBI:15940"

CHEBI_PYRIDOXAMINE_PHOSPHATE = "CHEBI:58451"
CHEBI_PYRIDOXAL_PHOSPHATE = "CHEBI:597326"
CHEBI_PYRIDOXINE_PHOSPHATE = "CHEBI:43097"


class ParticipantInput(TypedDict, total=False):
    """Subset of participant fields used in classifier tests."""

    chebi_id: str
    smiles: str
    inchi: str
    name: str
    count: int


def _reaction(left: list[ParticipantInput], right: list[ParticipantInput]) -> Reaction:
    return Reaction(
        left_participants=[Participant(**participant) for participant in left],
        right_participants=[Participant(**participant) for participant in right],
    )


def test_go_0008817_positive_corrinoid_adenosyltransferase() -> None:
    cls = CorrinoidAdenosyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_COB_II_ALAMIN, "name": "cob(II)alamin"},
            {"chebi_id": CHEBI_COB_II_ALAMIN, "name": "cob(II)alamin"},
            {"chebi_id": CHEBI_REDUCED_ETF, "name": "reduced ETF"},
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
        ],
        right=[
            {"chebi_id": CHEBI_ADENOSYLCOBALAMIN, "name": "adenosylcob(III)alamin"},
            {"chebi_id": CHEBI_ADENOSYLCOBALAMIN, "name": "adenosylcob(III)alamin"},
            {"chebi_id": CHEBI_TRIPHOSPHATE, "name": "triphosphate"},
            {"chebi_id": CHEBI_TRIPHOSPHATE, "name": "triphosphate"},
            {"chebi_id": CHEBI_FAD, "name": "oxidized ETF"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "adenosyltransferase" in result.explanation.lower()


def test_go_0008817_rejects_non_triphosphate_corrinoid_reaction() -> None:
    cls = CorrinoidAdenosyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_COB_II_INAMIDE, "name": "cob(II)inamide"},
            {"chebi_id": CHEBI_REDUCED_ETF, "name": "reduced ETF"},
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
        ],
        right=[
            {"chebi_id": CHEBI_ADENOSYLCOBINAMIDE, "name": "adenosylcobinamide"},
            {"chebi_id": CHEBI_DIPHOSPHATE, "name": "diphosphate"},
            {"chebi_id": CHEBI_FAD, "name": "oxidized ETF"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_go_0070733_positive_ampylase() -> None:
    cls = AMPylase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_TYROSYL_PROTEIN, "name": "L-tyrosyl-[protein]"},
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
        ],
        right=[
            {"chebi_id": CHEBI_AMPYL_TYROSYL_PROTEIN, "name": "O-(5'-adenylyl)-L-tyrosyl-[protein]"},
            {"chebi_id": CHEBI_DIPHOSPHATE, "name": "diphosphate"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "ampylase" in result.explanation.lower() or "adenylyl" in result.explanation.lower()


def test_go_0070733_rejects_ampylation_of_free_tyrosine() -> None:
    cls = AMPylase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_FREE_TYROSINE, "name": "L-tyrosine"},
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
        ],
        right=[
            {"chebi_id": CHEBI_AMPYL_TYROSYL_PROTEIN, "name": "O-(5'-adenylyl)-L-tyrosyl-[protein]"},
            {"chebi_id": CHEBI_DIPHOSPHATE, "name": "diphosphate"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_go_0052901_positive_spermine_oxidase() -> None:
    cls = SpermineOxidase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_SPERMINE, "name": "spermine"},
            {"chebi_id": CHEBI_O2, "name": "dioxygen"},
            {"chebi_id": CHEBI_H2O, "name": "water"},
        ],
        right=[
            {"chebi_id": CHEBI_AMINOPROPANAL, "name": "3-aminopropanal"},
            {"chebi_id": CHEBI_SPERMIDINE, "name": "spermidine"},
            {"chebi_id": CHEBI_H2O2, "name": "hydrogen peroxide"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "spermine oxidase" in result.explanation.lower()


def test_go_0052901_rejects_n1_acetylpolyamine_branch() -> None:
    cls = SpermineOxidase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_N1_ACETYLSPERMINE, "name": "N(1)-acetylspermine"},
            {"chebi_id": CHEBI_O2, "name": "dioxygen"},
            {"chebi_id": CHEBI_H2O, "name": "water"},
        ],
        right=[
            {"chebi_id": CHEBI_AMINOPROPANAL, "name": "3-aminopropanal"},
            {"chebi_id": CHEBI_SPERMIDINE, "name": "spermidine"},
            {"chebi_id": CHEBI_H2O2, "name": "hydrogen peroxide"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_go_0008112_positive_nicotinamide_n_methyltransferase() -> None:
    cls = NicotinamideNMethyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_NICOTINAMIDE, "name": "nicotinamide"},
            {"chebi_id": CHEBI_SAM, "name": "SAM"},
        ],
        right=[
            {"chebi_id": CHEBI_1_METHYLNICOTINAMIDE, "name": "1-methylnicotinamide"},
            {"chebi_id": CHEBI_SAH, "name": "SAH"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "nicotinamide" in result.explanation.lower()


def test_go_0008112_rejects_nicotinic_acid_methylation() -> None:
    cls = NicotinamideNMethyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_NICOTINIC_ACID, "name": "nicotinic acid"},
            {"chebi_id": CHEBI_SAM, "name": "SAM"},
        ],
        right=[
            {"chebi_id": CHEBI_1_METHYLNICOTINAMIDE, "name": "1-methylnicotinamide"},
            {"chebi_id": CHEBI_SAH, "name": "SAH"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_go_0004733_positive_pyridoxamine_phosphate_oxidase() -> None:
    cls = PyridoxaminePhosphateOxidase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_PYRIDOXAMINE_PHOSPHATE, "name": "pyridoxamine 5'-phosphate"},
            {"chebi_id": CHEBI_O2, "name": "dioxygen"},
            {"chebi_id": CHEBI_H2O, "name": "water"},
        ],
        right=[
            {"chebi_id": CHEBI_PYRIDOXAL_PHOSPHATE, "name": "pyridoxal 5'-phosphate"},
            {"chebi_id": CHEBI_H2O2, "name": "hydrogen peroxide"},
            {"chebi_id": CHEBI_NH4, "name": "ammonium"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "pyridoxamine phosphate oxidase" in result.explanation.lower()


def test_go_0004733_rejects_pyridoxine_phosphate_oxidation() -> None:
    cls = PyridoxaminePhosphateOxidase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_PYRIDOXINE_PHOSPHATE, "name": "pyridoxine 5'-phosphate"},
            {"chebi_id": CHEBI_O2, "name": "dioxygen"},
            {"chebi_id": CHEBI_H2O, "name": "water"},
        ],
        right=[
            {"chebi_id": CHEBI_PYRIDOXAL_PHOSPHATE, "name": "pyridoxal 5'-phosphate"},
            {"chebi_id": CHEBI_H2O2, "name": "hydrogen peroxide"},
            {"chebi_id": CHEBI_NH4, "name": "ammonium"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False
