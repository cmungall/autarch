"""Focused tests for the next uncovered GO support batch."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from typing_extensions import TypedDict

from autarch.datamodel import Participant, Reaction
from autarch.molecules import CHEBI_ACETYL_COA, CHEBI_COA, CHEBI_H2O, CHEBI_H_PLUS
from autarch.ontology.c_acyltransferase_activity import CAcyltransferaseActivity
from autarch.ontology.carbohydrate_phosphatase_activity import (
    CarbohydratePhosphataseActivity,
)
from autarch.ontology.hydroxycinnamoyltransferase_activity import (
    HydroxycinnamoyltransferaseActivity,
)
from autarch.ontology.limonene_monooxygenase_activity import (
    LimoneneMonooxygenaseActivity,
)
from autarch.ontology.phosphoprotein_phosphatase_activity import (
    PhosphoproteinPhosphataseActivity,
)
from autarch.ontology.protein_n_acetyltransferase_activity import (
    ProteinNAcetyltransferaseActivity,
)
from autarch.ontology.protein_n_acyltransferase_activity import (
    ProteinNAcyltransferaseActivity,
)
from autarch.ontology.ribonucleoside_triphosphate_phosphatase_activity import (
    RibonucleosideTriphosphatePhosphataseActivity,
)
from autarch.ontology.succinyltransferase_activity import SuccinyltransferaseActivity
from autarch.ontology.sugar_phosphatase_activity import SugarPhosphataseActivity
from autarch.ontology.udp_xylosyltransferase_activity import (
    UDPXylosyltransferaseActivity,
)
from autarch.ontology.xylosyltransferase_activity import XylosyltransferaseActivity

CHEBI_WATER = CHEBI_H2O
CHEBI_PHOSPHATE = "CHEBI:16838"
CHEBI_FRUCTOSE_BISPHOSPHATE = "CHEBI:32966"
CHEBI_FRUCTOSE_PHOSPHATE = "CHEBI:57634"
CHEBI_UDP_XYLOSE = "CHEBI:57632"
CHEBI_UDP = "CHEBI:58223"
CHEBI_UDP_GLUCOSE = "CHEBI:18066"
CHEBI_GLUCOSE = "CHEBI:4167"
CHEBI_LIMONENE_S = "CHEBI:15383"
CHEBI_O2 = "CHEBI:15379"
CHEBI_OXIDIZED_REDUCTASE = "CHEBI:58210"
CHEBI_ISOPIPERITENOL = "CHEBI:15406"
CHEBI_PHOSPHOPROTEIN = "CHEBI:61978"
CHEBI_PROTEIN = "CHEBI:46858"
CHEBI_SUCCINYL_COA = "CHEBI:57292"
CHEBI_ARGININE = "CHEBI:32682"
CHEBI_SUCCINYL_ARGININE = "CHEBI:58241"
CHEBI_ATP = "CHEBI:30616"
CHEBI_ADP = "CHEBI:456216"
CHEBI_DTTP = "CHEBI:37568"
CHEBI_DTDP = "CHEBI:58369"
CHEBI_PROTEIN_AMINO = "CHEBI:29969"
CHEBI_N_ACETYL_PROTEIN = "CHEBI:61930"
CHEBI_HYDROSULFIDE = "CHEBI:29919"
CHEBI_THIOACETATE = "CHEBI:30320"


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
        (CAcyltransferaseActivity, "RHEA:21036"),
        (SugarPhosphataseActivity, "RHEA:16689"),
        (CarbohydratePhosphataseActivity, "RHEA:16689"),
        (HydroxycinnamoyltransferaseActivity, "RHEA:12124"),
        (UDPXylosyltransferaseActivity, "RHEA:14721"),
        (XylosyltransferaseActivity, "RHEA:14721"),
        (LimoneneMonooxygenaseActivity, "RHEA:17945"),
        (PhosphoproteinPhosphataseActivity, "RHEA:10684"),
        (ProteinNAcyltransferaseActivity, "RHEA:21992"),
        (ProteinNAcetyltransferaseActivity, "RHEA:21992"),
        (SuccinyltransferaseActivity, "RHEA:15185"),
        (RibonucleosideTriphosphatePhosphataseActivity, "RHEA:19669"),
    ],
)
def test_batch26_positive_cached_examples(classifier, rhea_id) -> None:
    reaction = _load_rhea_reaction(rhea_id)
    assert classifier().check_membership(reaction).is_member is True


def test_sugar_phosphatase_rejects_phosphoprotein_hydrolysis() -> None:
    cls = SugarPhosphataseActivity()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_PHOSPHOPROTEIN, "name": "phosphoprotein"},
            {"chebi_id": CHEBI_WATER, "name": "water"},
        ],
        right=[
            {"chebi_id": CHEBI_PROTEIN, "name": "protein"},
            {"chebi_id": CHEBI_PHOSPHATE, "name": "phosphate"},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_udp_xylosyltransferase_rejects_udp_glucose_donor() -> None:
    cls = UDPXylosyltransferaseActivity()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_UDP_GLUCOSE, "name": "UDP-glucose"},
            {"chebi_id": CHEBI_GLUCOSE, "name": "glucose"},
        ],
        right=[
            {"chebi_id": CHEBI_UDP, "name": "UDP"},
            {"chebi_id": CHEBI_GLUCOSE, "name": "glycosylated product"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_limonene_monooxygenase_rejects_missing_oxygen() -> None:
    cls = LimoneneMonooxygenaseActivity()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_LIMONENE_S, "name": "(4S)-limonene"},
        ],
        right=[
            {"chebi_id": CHEBI_ISOPIPERITENOL, "name": "isopiperitenol"},
            {"chebi_id": CHEBI_OXIDIZED_REDUCTASE, "name": "oxidized reductase"},
            {"chebi_id": CHEBI_WATER, "name": "water"},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_protein_n_acetyltransferase_rejects_thioacetate_transfer() -> None:
    cls = ProteinNAcetyltransferaseActivity()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_ACETYL_COA, "name": "acetyl-CoA"},
            {"chebi_id": CHEBI_HYDROSULFIDE, "name": "hydrosulfide"},
        ],
        right=[
            {"chebi_id": CHEBI_COA, "name": "CoA"},
            {"chebi_id": CHEBI_THIOACETATE, "name": "thioacetate"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_ribonucleoside_triphosphate_phosphatase_rejects_deoxy_triphosphate() -> None:
    cls = RibonucleosideTriphosphatePhosphataseActivity()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_DTTP, "name": "dTTP"},
            {"chebi_id": CHEBI_WATER, "name": "water"},
        ],
        right=[
            {"chebi_id": CHEBI_DTDP, "name": "dTDP"},
            {"chebi_id": CHEBI_PHOSPHATE, "name": "phosphate"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is False
