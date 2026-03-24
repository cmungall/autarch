"""Focused tests for high-support batch 18 classifiers."""

from __future__ import annotations

from typing_extensions import TypedDict

from autarch.datamodel import Participant, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_ATP,
    CHEBI_CO2,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
    CHEBI_SAH,
    CHEBI_SAM,
)
from autarch.ontology.demethylmenaquinone_methyltransferase import (
    DemethylmenaquinoneMethyltransferase,
)
from autarch.ontology.gibberellin_20_oxidase import Gibberellin20Oxidase
from autarch.ontology.hexokinase import Hexokinase
from autarch.ontology.limonene_12_monooxygenase import (
    Limonene12Monooxygenase,
)
from autarch.ontology.nitroarene_dioxygenase import NitroareneDioxygenase

CHEBI_2_OXOGLUTARATE = "CHEBI:16810"
CHEBI_SUCCINATE = "CHEBI:30031"
CHEBI_GA12 = "CHEBI:58627"
CHEBI_GA9 = "CHEBI:73255"
CHEBI_GA24 = "CHEBI:143957"

CHEBI_D_GLUCOSE = "CHEBI:4167"
CHEBI_D_GLUCOSE_6_PHOSPHATE = "CHEBI:61548"
CHEBI_D_FRUCTOSE = "CHEBI:37721"

CHEBI_NITROBENZENE = "CHEBI:27798"
CHEBI_CATECHOL = "CHEBI:18135"
CHEBI_NITRITE = "CHEBI:16301"

CHEBI_DEMETHYLMENAQUINOL = "CHEBI:55437"
CHEBI_MENAQUINOL = "CHEBI:18151"
CHEBI_MENAQUINONE = "CHEBI:16374"

CHEBI_R_LIMONENE = "CHEBI:15382"
CHEBI_LIMONENE_12_EPOXIDE = "CHEBI:16431"


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


def test_go_0045544_positive_gibberellin_20_oxidase() -> None:
    cls = Gibberellin20Oxidase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_GA12, "name": "GA12"},
            {"chebi_id": CHEBI_2_OXOGLUTARATE, "name": "2-oxoglutarate", "count": 2},
            {"chebi_id": CHEBI_O2, "name": "dioxygen", "count": 3},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
        right=[
            {"chebi_id": CHEBI_GA9, "name": "GA9"},
            {"chebi_id": CHEBI_SUCCINATE, "name": "succinate", "count": 2},
            {"chebi_id": CHEBI_CO2, "name": "carbon dioxide", "count": 3},
            {"chebi_id": CHEBI_H2O, "name": "water", "count": 2},
        ],
    )

    result = cls.check_membership_impl(reaction)
    assert result.is_member is True
    assert "gibberellin 20-oxidase" in result.explanation.lower()


def test_go_0045544_rejects_wrong_gibberellin_product_branch() -> None:
    cls = Gibberellin20Oxidase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_GA12, "name": "GA12"},
            {"chebi_id": CHEBI_2_OXOGLUTARATE, "name": "2-oxoglutarate"},
            {"chebi_id": CHEBI_O2, "name": "dioxygen"},
        ],
        right=[
            {"chebi_id": CHEBI_GA24, "name": "GA24"},
            {"chebi_id": CHEBI_SUCCINATE, "name": "succinate"},
            {"chebi_id": CHEBI_CO2, "name": "carbon dioxide"},
            {"chebi_id": CHEBI_H2O, "name": "water"},
        ],
    )

    assert cls.check_membership_impl(reaction).is_member is False


def test_ec_2711_positive_hexokinase() -> None:
    cls = Hexokinase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_D_GLUCOSE, "name": "D-glucose"},
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
        ],
        right=[
            {"chebi_id": CHEBI_D_GLUCOSE_6_PHOSPHATE, "name": "D-glucose 6-phosphate"},
            {"chebi_id": CHEBI_ADP, "name": "ADP"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)
    assert result.is_member is True
    assert "hexokinase" in result.explanation.lower()


def test_ec_2711_rejects_wrong_hexose_product_branch() -> None:
    cls = Hexokinase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_D_GLUCOSE, "name": "D-glucose"},
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
        ],
        right=[
            {"chebi_id": CHEBI_D_FRUCTOSE, "name": "D-fructose"},
            {"chebi_id": CHEBI_ADP, "name": "ADP"},
        ],
    )

    assert cls.check_membership_impl(reaction).is_member is False


def test_ec_11413107_positive_limonene_12_monooxygenase() -> None:
    cls = Limonene12Monooxygenase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_R_LIMONENE, "name": "limonene"},
            {"chebi_id": CHEBI_NADPH, "name": "NADPH"},
            {"chebi_id": CHEBI_O2, "name": "dioxygen"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
        right=[
            {"chebi_id": CHEBI_LIMONENE_12_EPOXIDE, "name": "limonene 1,2-epoxide"},
            {"chebi_id": CHEBI_NADP_PLUS, "name": "NADP+"},
            {"chebi_id": CHEBI_H2O, "name": "water"},
        ],
    )

    result = cls.check_membership_impl(reaction)
    assert result.is_member is True
    assert "limonene 1,2-monooxygenase" in result.explanation.lower()


def test_ec_11413107_rejects_missing_epoxide_product() -> None:
    cls = Limonene12Monooxygenase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_R_LIMONENE, "name": "limonene"},
            {"chebi_id": CHEBI_NADPH, "name": "NADPH"},
            {"chebi_id": CHEBI_O2, "name": "dioxygen"},
        ],
        right=[
            {"chebi_id": CHEBI_CATECHOL, "name": "catechol"},
            {"chebi_id": CHEBI_NADP_PLUS, "name": "NADP+"},
            {"chebi_id": CHEBI_H2O, "name": "water"},
        ],
    )

    assert cls.check_membership_impl(reaction).is_member is False


def test_ec_1141223_positive_nitroarene_dioxygenase() -> None:
    cls = NitroareneDioxygenase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_NITROBENZENE, "name": "nitrobenzene"},
            {"chebi_id": CHEBI_NADH, "name": "NADH"},
            {"chebi_id": CHEBI_O2, "name": "dioxygen"},
        ],
        right=[
            {"chebi_id": CHEBI_CATECHOL, "name": "catechol"},
            {"chebi_id": CHEBI_NITRITE, "name": "nitrite"},
            {"chebi_id": CHEBI_NAD_PLUS, "name": "NAD+"},
        ],
    )

    result = cls.check_membership_impl(reaction)
    assert result.is_member is True
    assert "nitroarene dioxygenase" in result.explanation.lower()


def test_ec_1141223_rejects_nadph_branch() -> None:
    cls = NitroareneDioxygenase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_NITROBENZENE, "name": "nitrobenzene"},
            {"chebi_id": CHEBI_NADPH, "name": "NADPH"},
            {"chebi_id": CHEBI_O2, "name": "dioxygen"},
        ],
        right=[
            {"chebi_id": CHEBI_CATECHOL, "name": "catechol"},
            {"chebi_id": CHEBI_NITRITE, "name": "nitrite"},
            {"chebi_id": CHEBI_NADP_PLUS, "name": "NADP+"},
        ],
    )

    assert cls.check_membership_impl(reaction).is_member is False


def test_go_0043770_positive_demethylmenaquinone_methyltransferase() -> None:
    cls = DemethylmenaquinoneMethyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_DEMETHYLMENAQUINOL, "name": "demethylmenaquinol"},
            {"chebi_id": CHEBI_SAM, "name": "SAM"},
        ],
        right=[
            {"chebi_id": CHEBI_MENAQUINOL, "name": "menaquinol"},
            {"chebi_id": CHEBI_SAH, "name": "SAH"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)
    assert result.is_member is True
    assert "demethylmenaquinone" in result.explanation.lower()


def test_go_0043770_rejects_already_methylated_menaquinone() -> None:
    cls = DemethylmenaquinoneMethyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_MENAQUINONE, "name": "menaquinone"},
            {"chebi_id": CHEBI_SAM, "name": "SAM"},
        ],
        right=[
            {"chebi_id": CHEBI_MENAQUINOL, "name": "menaquinol"},
            {"chebi_id": CHEBI_SAH, "name": "SAH"},
        ],
    )

    assert cls.check_membership_impl(reaction).is_member is False
