"""Focused tests for the next GO support batch."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from typing_extensions import TypedDict

from autarch.datamodel import Participant, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS, CHEBI_O2
from autarch.ontology.arachidonate_5_lipoxygenase_activity import (
    Arachidonate5LipoxygenaseActivity,
)
from autarch.ontology.cinnamyl_alcohol_dehydrogenase_activity import (
    CinnamylAlcoholDehydrogenaseActivity,
)
from autarch.ontology.galactosylceramide_sulfotransferase_activity import (
    GalactosylceramideSulfotransferaseActivity,
)
from autarch.ontology.sterol_12_alpha_hydroxylase_activity import (
    Sterol12AlphaHydroxylaseActivity,
)
from autarch.ontology.sterol_14_demethylase_activity import (
    Sterol14DemethylaseActivity,
)
from autarch.ontology.three_beta_hydroxysteroid_3_dehydrogenase_nad_p_activity import (
    ThreeBetaHydroxysteroid3DehydrogenaseNADPActivity,
)
from autarch.ontology.two_hydroxyacyl_coa_lyase_activity import (
    TwoHydroxyacylCoALyaseActivity,
)
from autarch.ontology.zeaxanthin_epoxidase_activity import ZeaxanthinEpoxidaseActivity

CHEBI_FORMYL_COA = "CHEBI:57376"
CHEBI_ACETYL_COA = "CHEBI:57288"
CHEBI_ZEAXANTHIN = "CHEBI:27547"
CHEBI_ANTHERAXANTHIN = "CHEBI:27867"
CHEBI_VIOLAXANTHIN = "CHEBI:35288"
CHEBI_REDUCED_FERREDOXIN = "CHEBI:33738"
CHEBI_OXIDIZED_FERREDOXIN = "CHEBI:33737"
CHEBI_CINNAMYL_ALCOHOL = "CHEBI:33227"
CHEBI_CINNAMALDEHYDE = "CHEBI:16731"
CHEBI_NADP_PLUS = "CHEBI:58349"
CHEBI_NADPH = "CHEBI:57783"
CHEBI_STEROL_SUBSTRATE = "CHEBI:17791"
CHEBI_STEROL_PRODUCT = "CHEBI:30109"
CHEBI_REDUCED_HEMOPROTEIN_REDUCTASE = "CHEBI:57618"
CHEBI_OXIDIZED_HEMOPROTEIN_REDUCTASE = "CHEBI:58210"
CHEBI_FORMATE = "CHEBI:15740"
CHEBI_ARACHIDONATE = "CHEBI:32395"
CHEBI_5_HPETE = "CHEBI:57450"
CHEBI_LTA4 = "CHEBI:57463"
CHEBI_HYDROXYSTEROID = "CHEBI:18378"
CHEBI_OXOSTEROID = "CHEBI:16495"
CHEBI_GALACTOSYLCERAMIDE = "CHEBI:18390"
CHEBI_SULFO_GALACTOSYLCERAMIDE = "CHEBI:75956"
CHEBI_PAPS = "CHEBI:58339"
CHEBI_PAP = "CHEBI:58343"


class ParticipantInput(TypedDict, total=False):
    chebi_id: str
    smiles: str
    inchi: str
    name: str
    count: int



def _load_rhea_reaction(rhea_id: str) -> Reaction:
    chebi_to_smiles = json.loads(Path("cache/chebi_smiles.json").read_text())
    for line in Path("cache/rhea_reactions.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        term = json.loads(line)
        if term["rhea_id"] == rhea_id:
            reaction = term["reaction"]
            left_participants = []
            right_participants = []
            for participant in reaction["left_participants"]:
                chebi_id = participant.get("chebi_id")
                smiles = participant.get("smiles") or (
                    chebi_to_smiles.get(chebi_id) if chebi_id else None
                )
                left_participants.append(
                    Participant(
                        chebi_id=chebi_id,
                        smiles=smiles,
                        inchi=participant.get("inchi"),
                        location=participant.get("location"),
                        name=participant.get("name"),
                        polymer_index=participant.get("polymer_index"),
                        polymer_type=participant.get("polymer_type"),
                    )
                )
            for participant in reaction["right_participants"]:
                chebi_id = participant.get("chebi_id")
                smiles = participant.get("smiles") or (
                    chebi_to_smiles.get(chebi_id) if chebi_id else None
                )
                right_participants.append(
                    Participant(
                        chebi_id=chebi_id,
                        smiles=smiles,
                        inchi=participant.get("inchi"),
                        location=participant.get("location"),
                        name=participant.get("name"),
                        polymer_index=participant.get("polymer_index"),
                        polymer_type=participant.get("polymer_type"),
                    )
                )
            return Reaction(
                left_participants=left_participants,
                right_participants=right_participants,
                label=term.get("label", ""),
            )
    raise AssertionError(f"Missing cached RHEA reaction: {rhea_id}")



def _reaction(left: list[ParticipantInput], right: list[ParticipantInput]) -> Reaction:
    return Reaction(
        left_participants=[Participant(**participant) for participant in left],
        right_participants=[Participant(**participant) for participant in right],
    )


@pytest.mark.parametrize(
    ("classifier", "rhea_id"),
    [
        (TwoHydroxyacylCoALyaseActivity, "RHEA:25379"),
        (ZeaxanthinEpoxidaseActivity, "RHEA:24084"),
        (CinnamylAlcoholDehydrogenaseActivity, "RHEA:10392"),
        (Sterol14DemethylaseActivity, "RHEA:14917"),
        (Sterol12AlphaHydroxylaseActivity, "RHEA:15261"),
        (Arachidonate5LipoxygenaseActivity, "RHEA:17485"),
        (ThreeBetaHydroxysteroid3DehydrogenaseNADPActivity, "RHEA:18409"),
        (GalactosylceramideSulfotransferaseActivity, "RHEA:20613"),
    ],
)
def test_batch27_positive_cached_examples(classifier, rhea_id) -> None:
    reaction = _load_rhea_reaction(rhea_id)
    assert classifier().check_membership(reaction).is_member is True


def test_two_hydroxyacyl_coa_lyase_rejects_simple_thioester_cleavage() -> None:
    cls = TwoHydroxyacylCoALyaseActivity()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_ACETYL_COA, "name": "acetyl-CoA"}],
        right=[{"chebi_id": CHEBI_FORMYL_COA, "name": "formyl-CoA"}],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_zeaxanthin_epoxidase_rejects_missing_ferredoxin() -> None:
    cls = ZeaxanthinEpoxidaseActivity()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_ZEAXANTHIN, "name": "zeaxanthin"},
            {"chebi_id": CHEBI_O2, "name": "dioxygen"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
        right=[
            {"chebi_id": CHEBI_ANTHERAXANTHIN, "name": "antheraxanthin"},
            {"chebi_id": CHEBI_H2O, "name": "water"},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_cinnamyl_alcohol_dehydrogenase_rejects_nad_branch() -> None:
    cls = CinnamylAlcoholDehydrogenaseActivity()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_CINNAMYL_ALCOHOL}, {"chebi_id": CHEBI_NADPH}],
        right=[{"chebi_id": CHEBI_CINNAMALDEHYDE}, {"chebi_id": CHEBI_NADP_PLUS}],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_sterol_14_demethylase_rejects_missing_formate() -> None:
    cls = Sterol14DemethylaseActivity()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_STEROL_SUBSTRATE},
            {"chebi_id": CHEBI_REDUCED_HEMOPROTEIN_REDUCTASE},
            {"chebi_id": CHEBI_REDUCED_HEMOPROTEIN_REDUCTASE},
            {"chebi_id": CHEBI_REDUCED_HEMOPROTEIN_REDUCTASE},
            {"chebi_id": CHEBI_O2},
            {"chebi_id": CHEBI_O2},
            {"chebi_id": CHEBI_O2},
        ],
        right=[
            {"chebi_id": CHEBI_STEROL_PRODUCT},
            {"chebi_id": CHEBI_OXIDIZED_HEMOPROTEIN_REDUCTASE},
            {"chebi_id": CHEBI_OXIDIZED_HEMOPROTEIN_REDUCTASE},
            {"chebi_id": CHEBI_OXIDIZED_HEMOPROTEIN_REDUCTASE},
            {"chebi_id": CHEBI_H2O},
            {"chebi_id": CHEBI_H2O},
            {"chebi_id": CHEBI_H2O},
            {"chebi_id": CHEBI_H_PLUS},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_arachidonate_5_lipoxygenase_rejects_unrelated_water_release() -> None:
    cls = Arachidonate5LipoxygenaseActivity()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_5_HPETE}],
        right=[{"chebi_id": CHEBI_ZEAXANTHIN}, {"chebi_id": CHEBI_H2O}],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_steroid_dehydrogenase_rejects_without_nadp_coupling() -> None:
    cls = ThreeBetaHydroxysteroid3DehydrogenaseNADPActivity()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_HYDROXYSTEROID}],
        right=[{"chebi_id": CHEBI_OXOSTEROID}],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_galactosylceramide_sulfotransferase_rejects_missing_pap() -> None:
    cls = GalactosylceramideSulfotransferaseActivity()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_GALACTOSYLCERAMIDE}, {"chebi_id": CHEBI_PAPS}],
        right=[{"chebi_id": CHEBI_SULFO_GALACTOSYLCERAMIDE}, {"chebi_id": CHEBI_H_PLUS}],
    )
    assert cls.check_membership_impl(reaction).is_member is False
