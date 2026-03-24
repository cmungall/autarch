"""Focused tests for a larger support-driven classifier batch."""

from __future__ import annotations

from typing_extensions import TypedDict

from autarch.datamodel import Participant, Reaction
from autarch.molecules import CHEBI_ADP, CHEBI_ATP, CHEBI_H2O, CHEBI_H_PLUS, CHEBI_PHOSPHATE
from autarch.ontology.abc_type_polar_amino_acid_transporter import ABCTypePolarAminoAcidTransporter
from autarch.ontology.aromatic_amino_acid_transaminase import AromaticAminoAcidTransaminase
from autarch.ontology.gellan_tetrasaccharide_unsaturated_glucuronosyl_hydrolase import (
    GellanTetrasaccharideUnsaturatedGlucuronosylHydrolase,
)
from autarch.ontology.gibberellin_a4_carboxyl_methyltransferase import GibberellinA4CarboxylMethyltransferase
from autarch.ontology.gibberellin_a9_o_methyltransferase import GibberellinA9OMethyltransferase
from autarch.ontology.quinoline_2_oxidoreductase import Quinoline2Oxidoreductase

CHEBI_SAM = "CHEBI:59789"
CHEBI_SAH = "CHEBI:57856"
CHEBI_GIBBERELLIN_A4 = "CHEBI:73251"
CHEBI_GIBBERELLIN_A4_METHYL_ESTER = "CHEBI:73252"
CHEBI_GIBBERELLIN_A34 = "CHEBI:73258"
CHEBI_GIBBERELLIN_A34_METHYL_ESTER = "CHEBI:73260"
CHEBI_UNSATURATED_GELLAN = "CHEBI:134390"
CHEBI_GELLAN_TRISACCHARIDE = "CHEBI:134389"
CHEBI_5_DEHYDRO_4_DEOXY_GLUCURONATE = "CHEBI:17117"
CHEBI_ARGININE = "CHEBI:32682"
CHEBI_LYSINE = "CHEBI:32551"
CHEBI_QUINOLINE = "CHEBI:17362"
CHEBI_QUINOLINONE = "CHEBI:18289"
CHEBI_GENERIC_OXIDANT = "CHEBI:15377"
CHEBI_GENERIC_REDUCED_OXIDANT = "CHEBI:17499"
CHEBI_TYROSINE = "CHEBI:58315"
CHEBI_2_OXOGLUTARATE = "CHEBI:16810"
CHEBI_TYROSINE_KETO_ACID = "CHEBI:36242"
CHEBI_GLUTAMATE = "CHEBI:29985"


class ParticipantInput(TypedDict, total=False):
    chebi_id: str
    name: str
    count: int
    location: str


def _reaction(left: list[ParticipantInput], right: list[ParticipantInput]) -> Reaction:
    return Reaction(
        left_participants=[Participant(**participant) for participant in left],
        right_participants=[Participant(**participant) for participant in right],
    )


def test_ec_211275_positive_gibberellin_a9_o_methyltransferase() -> None:
    cls = GibberellinA9OMethyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_GIBBERELLIN_A4},
            {"chebi_id": CHEBI_SAM},
        ],
        right=[
            {"chebi_id": CHEBI_GIBBERELLIN_A4_METHYL_ESTER},
            {"chebi_id": CHEBI_SAH},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_211276_positive_gibberellin_a4_carboxyl_methyltransferase() -> None:
    cls = GibberellinA4CarboxylMethyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_GIBBERELLIN_A34},
            {"chebi_id": CHEBI_SAM},
        ],
        right=[
            {"chebi_id": CHEBI_GIBBERELLIN_A34_METHYL_ESTER},
            {"chebi_id": CHEBI_SAH},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_321179_positive_gellan_unsaturated_glucuronosyl_hydrolase() -> None:
    cls = GellanTetrasaccharideUnsaturatedGlucuronosylHydrolase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_UNSATURATED_GELLAN},
            {"chebi_id": CHEBI_H2O},
        ],
        right=[
            {"chebi_id": CHEBI_GELLAN_TRISACCHARIDE},
            {"chebi_id": CHEBI_5_DEHYDRO_4_DEOXY_GLUCURONATE},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_7421_positive_abc_polar_amino_acid_transporter() -> None:
    cls = ABCTypePolarAminoAcidTransporter()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_ARGININE, "location": "out"},
            {"chebi_id": CHEBI_ATP},
            {"chebi_id": CHEBI_H2O},
        ],
        right=[
            {"chebi_id": CHEBI_ARGININE, "location": "in"},
            {"chebi_id": CHEBI_ADP},
            {"chebi_id": CHEBI_PHOSPHATE},
            {"chebi_id": CHEBI_H_PLUS},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_7421_positive_lysine_transport() -> None:
    cls = ABCTypePolarAminoAcidTransporter()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_LYSINE, "location": "out"},
            {"chebi_id": CHEBI_ATP},
            {"chebi_id": CHEBI_H2O},
        ],
        right=[
            {"chebi_id": CHEBI_LYSINE, "location": "in"},
            {"chebi_id": CHEBI_ADP},
            {"chebi_id": CHEBI_PHOSPHATE},
            {"chebi_id": CHEBI_H_PLUS},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_139917_positive_quinoline_2_oxidoreductase() -> None:
    cls = Quinoline2Oxidoreductase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_QUINOLINE},
            {"chebi_id": "CHEBI:13193"},
            {"chebi_id": CHEBI_GENERIC_OXIDANT, "count": 2},
        ],
        right=[
            {"chebi_id": CHEBI_QUINOLINONE},
            {"chebi_id": CHEBI_GENERIC_REDUCED_OXIDANT},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_26157_positive_aromatic_amino_acid_transaminase() -> None:
    cls = AromaticAminoAcidTransaminase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_TYROSINE},
            {"chebi_id": CHEBI_2_OXOGLUTARATE},
        ],
        right=[
            {"chebi_id": CHEBI_TYROSINE_KETO_ACID},
            {"chebi_id": CHEBI_GLUTAMATE},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True
