"""Focused tests for benchmark-supported exact classifiers."""

from __future__ import annotations

from typing_extensions import TypedDict

from autarch.datamodel import Participant, Reaction
from autarch.molecules import CHEBI_ADP, CHEBI_ATP, CHEBI_H2O, CHEBI_H_PLUS, CHEBI_NADPH, CHEBI_NADP_PLUS, CHEBI_O2
from autarch.ontology.monocyclic_monoterpene_ketone_monooxygenase import MonocyclicMonoterpeneKetoneMonooxygenase
from autarch.ontology.one_piperideine_2_carboxylate_one_pyrroline_2_carboxylate_reductase_nadh_or_nadph import (
    OnePiperideine2CarboxylateOnePyrroline2CarboxylateReductaseNADHOrNADPH,
)
from autarch.ontology.quinate_shikimate_dehydrogenase_nad_p import QuinateShikimateDehydrogenaseNADP
from autarch.ontology.two_deamino_two_hydroxyneamine_1_alpha_d_kanosaminyltransferase import (
    TwoDeaminoTwoHydroxyneamine1AlphaDKanosaminyltransferase,
)
from autarch.ontology.uridine_cytidine_kinase import UridineCytidineKinase

CHEBI_SHIKIMATE = "CHEBI:36208"
CHEBI_3_DEHYDROSHIKIMATE = "CHEBI:16630"
CHEBI_PIPECOLATE = "CHEBI:61185"
CHEBI_PIPERIDEINE_2_CARBOXYLATE = "CHEBI:77631"
CHEBI_URIDINE = "CHEBI:16704"
CHEBI_UMP = "CHEBI:57865"
CHEBI_DEOXYTHYMIDINE = "CHEBI:17748"
CHEBI_DTMP = "CHEBI:63528"
CHEBI_NEAMINE = "CHEBI:65076"
CHEBI_UDP_ALPHA_D_KANOSAMINE = "CHEBI:71964"
CHEBI_KANAMYCIN_B = "CHEBI:58549"
CHEBI_UDP = "CHEBI:58223"
CHEBI_GDP_BETA_L_FUCOSE = "CHEBI:57273"
CHEBI_MENTHONE = "CHEBI:15410"
CHEBI_MENTHONE_LACTONE = "CHEBI:50250"


class ParticipantInput(TypedDict, total=False):
    chebi_id: str
    name: str
    count: int


def _reaction(left: list[ParticipantInput], right: list[ParticipantInput]) -> Reaction:
    return Reaction(
        left_participants=[Participant(**participant) for participant in left],
        right_participants=[Participant(**participant) for participant in right],
    )


def test_ec_111282_positive_quinate_shikimate_dehydrogenase() -> None:
    cls = QuinateShikimateDehydrogenaseNADP()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_SHIKIMATE},
            {"chebi_id": CHEBI_NADP_PLUS},
        ],
        right=[
            {"chebi_id": CHEBI_3_DEHYDROSHIKIMATE},
            {"chebi_id": CHEBI_NADPH},
            {"chebi_id": CHEBI_H_PLUS},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_1511_positive_piperideine_pyrroline_reductase() -> None:
    cls = OnePiperideine2CarboxylateOnePyrroline2CarboxylateReductaseNADHOrNADPH()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_PIPECOLATE},
            {"chebi_id": CHEBI_NADP_PLUS},
        ],
        right=[
            {"chebi_id": CHEBI_PIPERIDEINE_2_CARBOXYLATE},
            {"chebi_id": CHEBI_NADPH},
            {"chebi_id": CHEBI_H_PLUS},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_27148_positive_uridine_cytidine_kinase() -> None:
    cls = UridineCytidineKinase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_URIDINE},
            {"chebi_id": CHEBI_ATP},
        ],
        right=[
            {"chebi_id": CHEBI_UMP},
            {"chebi_id": CHEBI_ADP},
            {"chebi_id": CHEBI_H_PLUS},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_27148_rejects_deoxynucleoside_kinase() -> None:
    cls = UridineCytidineKinase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_DEOXYTHYMIDINE},
            {"chebi_id": CHEBI_ATP},
        ],
        right=[
            {"chebi_id": CHEBI_DTMP},
            {"chebi_id": CHEBI_ADP},
            {"chebi_id": CHEBI_H_PLUS},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_ec_241301_positive_kanosaminyltransferase() -> None:
    cls = TwoDeaminoTwoHydroxyneamine1AlphaDKanosaminyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_NEAMINE},
            {"chebi_id": CHEBI_UDP_ALPHA_D_KANOSAMINE},
        ],
        right=[
            {"chebi_id": CHEBI_KANAMYCIN_B},
            {"chebi_id": CHEBI_UDP},
            {"chebi_id": CHEBI_H_PLUS},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_241301_rejects_wrong_sugar_donor() -> None:
    cls = TwoDeaminoTwoHydroxyneamine1AlphaDKanosaminyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_NEAMINE},
            {"chebi_id": CHEBI_GDP_BETA_L_FUCOSE},
        ],
        right=[
            {"chebi_id": CHEBI_KANAMYCIN_B},
            {"chebi_id": CHEBI_UDP},
            {"chebi_id": CHEBI_H_PLUS},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_ec_11413105_positive_monoterpene_ketone_monooxygenase() -> None:
    cls = MonocyclicMonoterpeneKetoneMonooxygenase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_MENTHONE},
            {"chebi_id": CHEBI_NADPH},
            {"chebi_id": CHEBI_O2},
            {"chebi_id": CHEBI_H_PLUS},
        ],
        right=[
            {"chebi_id": CHEBI_MENTHONE_LACTONE},
            {"chebi_id": CHEBI_NADP_PLUS},
            {"chebi_id": CHEBI_H2O},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_11413105_requires_oxygen() -> None:
    cls = MonocyclicMonoterpeneKetoneMonooxygenase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_MENTHONE},
            {"chebi_id": CHEBI_NADPH},
            {"chebi_id": CHEBI_H_PLUS},
        ],
        right=[
            {"chebi_id": CHEBI_MENTHONE_LACTONE},
            {"chebi_id": CHEBI_NADP_PLUS},
            {"chebi_id": CHEBI_H2O},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is False
