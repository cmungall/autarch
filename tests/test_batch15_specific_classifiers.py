"""Focused tests for batch 15 specific classifiers."""

from __future__ import annotations

from typing_extensions import TypedDict

from autarch.datamodel import Participant, Reaction
from autarch.molecules import CHEBI_H2O2, CHEBI_H_PLUS, CHEBI_O2, CHEBI_SAH, CHEBI_SAM
from autarch.ontology.galactose_oxidase import GalactoseOxidase
from autarch.ontology.methionine_s_methyltransferase import MethionineSMethyltransferase
from autarch.ontology.precorrin_3b_c17_methyltransferase import (
    Precorrin3BC17Methyltransferase,
)
from autarch.ontology.ten_hydroxydihydrosanguinarine_10_o_methyltransferase import (
    TenHydroxydihydrosanguinarine10OMethyltransferase,
)
from autarch.ontology.tyramine_n_methyltransferase import TyramineNMethyltransferase

CHEBI_10_HYDROXYDIHYDROSANGUINARINE = "CHEBI:15878"
CHEBI_DIHYDROCHELIRUBINE = "CHEBI:17789"
CHEBI_PRECORRIN_3B = "CHEBI:77870"
CHEBI_PRECORRIN_4 = "CHEBI:57769"
CHEBI_METHIONINE = "CHEBI:57844"
CHEBI_S_METHYL_METHIONINE = "CHEBI:58252"
CHEBI_HOMOCYSTEINE = "CHEBI:58199"
CHEBI_TYRAMINE = "CHEBI:327995"
CHEBI_N_METHYLTYRAMINE = "CHEBI:58155"
CHEBI_GALACTOSE = "CHEBI:4139"
CHEBI_GALACTO_HEXODIALDOSE = "CHEBI:16222"
CHEBI_GLUCOSE = "CHEBI:4167"
CHEBI_2_DEHYDRO_GLUCOSE = "CHEBI:16609"


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


def test_go_0030779_positive_ten_hydroxydihydrosanguinarine_10_o_methyltransferase() -> None:
    cls = TenHydroxydihydrosanguinarine10OMethyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_10_HYDROXYDIHYDROSANGUINARINE, "name": "10-hydroxydihydrosanguinarine"},
            {"chebi_id": CHEBI_SAM, "name": "SAM"},
        ],
        right=[
            {"chebi_id": CHEBI_DIHYDROCHELIRUBINE, "name": "dihydrochelirubine"},
            {"chebi_id": CHEBI_SAH, "name": "SAH"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "methyltransferase" in result.explanation.lower()


def test_go_0030779_rejects_unmethylated_product() -> None:
    cls = TenHydroxydihydrosanguinarine10OMethyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_10_HYDROXYDIHYDROSANGUINARINE, "name": "10-hydroxydihydrosanguinarine"},
            {"chebi_id": CHEBI_SAM, "name": "SAM"},
        ],
        right=[
            {"chebi_id": CHEBI_10_HYDROXYDIHYDROSANGUINARINE, "name": "10-hydroxydihydrosanguinarine"},
            {"chebi_id": CHEBI_SAH, "name": "SAH"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_go_0030789_positive_precorrin_3b_c17_methyltransferase() -> None:
    cls = Precorrin3BC17Methyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_PRECORRIN_3B, "name": "precorrin-3B"},
            {"chebi_id": CHEBI_SAM, "name": "SAM"},
        ],
        right=[
            {"chebi_id": CHEBI_PRECORRIN_4, "name": "precorrin-4"},
            {"chebi_id": CHEBI_SAH, "name": "SAH"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "precorrin" in result.explanation.lower()


def test_go_0030789_rejects_missing_precorrin_4_product() -> None:
    cls = Precorrin3BC17Methyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_PRECORRIN_3B, "name": "precorrin-3B"},
            {"chebi_id": CHEBI_SAM, "name": "SAM"},
        ],
        right=[
            {"chebi_id": CHEBI_PRECORRIN_3B, "name": "precorrin-3B"},
            {"chebi_id": CHEBI_SAH, "name": "SAH"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_go_0030732_positive_methionine_s_methyltransferase() -> None:
    cls = MethionineSMethyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_METHIONINE, "name": "L-methionine"},
            {"chebi_id": CHEBI_SAM, "name": "SAM"},
        ],
        right=[
            {"chebi_id": CHEBI_S_METHYL_METHIONINE, "name": "S-methyl-L-methionine"},
            {"chebi_id": CHEBI_SAH, "name": "SAH"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "methionine" in result.explanation.lower()


def test_go_0030732_rejects_homocysteine_methylation() -> None:
    cls = MethionineSMethyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_HOMOCYSTEINE, "name": "L-homocysteine"},
            {"chebi_id": CHEBI_SAM, "name": "SAM"},
        ],
        right=[
            {"chebi_id": CHEBI_METHIONINE, "name": "L-methionine"},
            {"chebi_id": CHEBI_SAH, "name": "SAH"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_go_0030738_positive_tyramine_n_methyltransferase() -> None:
    cls = TyramineNMethyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_TYRAMINE, "name": "tyramine"},
            {"chebi_id": CHEBI_SAM, "name": "SAM"},
        ],
        right=[
            {"chebi_id": CHEBI_N_METHYLTYRAMINE, "name": "N-methyltyramine"},
            {"chebi_id": CHEBI_SAH, "name": "SAH"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "tyramine" in result.explanation.lower()


def test_go_0030738_rejects_unmethylated_tyramine_product() -> None:
    cls = TyramineNMethyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_TYRAMINE, "name": "tyramine"},
            {"chebi_id": CHEBI_SAM, "name": "SAM"},
        ],
        right=[
            {"chebi_id": CHEBI_TYRAMINE, "name": "tyramine"},
            {"chebi_id": CHEBI_SAH, "name": "SAH"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_go_0045480_positive_galactose_oxidase() -> None:
    cls = GalactoseOxidase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_GALACTOSE, "name": "D-galactose"},
            {"chebi_id": CHEBI_O2, "name": "dioxygen"},
        ],
        right=[
            {"chebi_id": CHEBI_GALACTO_HEXODIALDOSE, "name": "D-galacto-hexodialdose"},
            {"chebi_id": CHEBI_H2O2, "name": "hydrogen peroxide"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "galactose oxidase" in result.explanation.lower()


def test_go_0045480_rejects_glucose_oxidation() -> None:
    cls = GalactoseOxidase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_GLUCOSE, "name": "D-glucose"},
            {"chebi_id": CHEBI_O2, "name": "dioxygen"},
        ],
        right=[
            {"chebi_id": CHEBI_2_DEHYDRO_GLUCOSE, "name": "2-dehydro-D-glucose"},
            {"chebi_id": CHEBI_H2O2, "name": "hydrogen peroxide"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False
