"""Focused tests for uncovered EC/GO classes added after the support-threshold scan."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from typing_extensions import TypedDict

from autarch.datamodel import Participant, Reaction
from autarch.evaluation import evaluate_reaction_class
from autarch.molecules import CHEBI_H_PLUS, CHEBI_SAH, CHEBI_SAM
from autarch.ontology.histone_methyltransferase_activity import (
    HistoneMethyltransferaseActivity,
)
from autarch.ontology.histone_modifying_activity import HistoneModifyingActivity
from autarch.ontology.lysine_n_methyltransferase_activity import (
    LysineNMethyltransferaseActivity,
)
from autarch.ontology.oxidoreductase_acting_on_the_ch_oh_group_of_donors_quinone_or_similar_compound_as_acceptor import (
    OxidoreductaseActingOnTheCHOHGroupOfDonorsQuinoneOrSimilarCompoundAsAcceptor,
)
from autarch.ontology.oxidosqualene_cyclase_activity import (
    OxidosqualeneCyclaseActivity,
)
from autarch.ontology.protein_lysine_n_methyltransferase_activity import (
    ProteinLysineNMethyltransferaseActivity,
)
from autarch.ontology.pseudouridine_synthase_activity import (
    PseudouridineSynthaseActivity,
)
from autarch.ontology.s_acyltransferase_activity import SAcyltransferaseActivity
from autarch.ontology.trna_guanine_methyltransferase_activity import (
    TRNAGuanineMethyltransferaseActivity,
)

CHEBI_PROTEIN_LYSINE = "CHEBI:29969"
CHEBI_METHYL_PROTEIN_LYSINE = "CHEBI:61929"
CHEBI_URIDINE_IN_RNA = "CHEBI:65315"
CHEBI_PSEUDOURIDINE_IN_RNA = "CHEBI:65314"
CHEBI_GUANOSINE_IN_TRNA = "CHEBI:74269"
CHEBI_N1_METHYLGUANOSINE_IN_TRNA = "CHEBI:73542"
CHEBI_QUINONE_GENERIC = "CHEBI:132124"
CHEBI_QUINOL_GENERIC = "CHEBI:24646"
CHEBI_QUINATE = "CHEBI:29751"
CHEBI_3_DEHYDROQUINATE = "CHEBI:32364"
CHEBI_HYDROSULFIDE = "CHEBI:29919"
CHEBI_ACETYL_COA = "CHEBI:57288"
CHEBI_THIOACETATE = "CHEBI:30320"
CHEBI_COA = "CHEBI:57287"
CHEBI_OXIDOSQUALENE = "CHEBI:15441"
CHEBI_LANOSTEROL = "CHEBI:16521"
CHEBI_GLYCEROL_3_PHOSPHATE = "CHEBI:57597"
CHEBI_DHAP = "CHEBI:57642"
CHEBI_NAD_PLUS = "CHEBI:57540"
CHEBI_NADH = "CHEBI:57945"


class ParticipantInput(TypedDict, total=False):
    """Subset of participant fields used in classifier tests."""

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
        (
            OxidoreductaseActingOnTheCHOHGroupOfDonorsQuinoneOrSimilarCompoundAsAcceptor,
            "RHEA:23672",
        ),
        (ProteinLysineNMethyltransferaseActivity, "RHEA:21556"),
        (LysineNMethyltransferaseActivity, "RHEA:21556"),
        (HistoneMethyltransferaseActivity, "RHEA:10024"),
        (HistoneModifyingActivity, "RHEA:10024"),
        (PseudouridineSynthaseActivity, "RHEA:42532"),
        (SAcyltransferaseActivity, "RHEA:16625"),
        (OxidosqualeneCyclaseActivity, "RHEA:14621"),
        (TRNAGuanineMethyltransferaseActivity, "RHEA:36899"),
    ],
)
def test_batch25_positive_cached_examples(classifier, rhea_id) -> None:
    reaction = _load_rhea_reaction(rhea_id)
    assert classifier().check_membership(reaction).is_member is True


def test_go_0016901_rejects_nad_branch() -> None:
    cls = OxidoreductaseActingOnTheCHOHGroupOfDonorsQuinoneOrSimilarCompoundAsAcceptor()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_QUINATE, "name": "quinate"},
            {"chebi_id": CHEBI_NAD_PLUS, "name": "NAD+"},
        ],
        right=[
            {"chebi_id": CHEBI_3_DEHYDROQUINATE, "name": "3-dehydroquinate"},
            {"chebi_id": CHEBI_NADH, "name": "NADH"},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_go_0016279_rejects_unmodified_protein_pair() -> None:
    cls = ProteinLysineNMethyltransferaseActivity()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_PROTEIN_LYSINE, "name": "protein lysine"},
            {"chebi_id": CHEBI_SAM, "name": "SAM"},
        ],
        right=[
            {"chebi_id": CHEBI_PROTEIN_LYSINE, "name": "protein lysine"},
            {"chebi_id": CHEBI_SAH, "name": "SAH"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_go_0009982_rejects_identity_reaction() -> None:
    cls = PseudouridineSynthaseActivity()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_URIDINE_IN_RNA, "name": "uridine in RNA"}],
        right=[{"chebi_id": CHEBI_URIDINE_IN_RNA, "name": "uridine in RNA"}],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_go_0016417_rejects_thioester_hydrolysis() -> None:
    cls = SAcyltransferaseActivity()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_ACETYL_COA, "name": "acetyl-CoA"}],
        right=[
            {"chebi_id": CHEBI_COA, "name": "CoA"},
            {"chebi_id": CHEBI_THIOACETATE, "name": "thioacetate"},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_go_0031559_rejects_uncyclized_product() -> None:
    cls = OxidosqualeneCyclaseActivity()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_OXIDOSQUALENE, "name": "oxidosqualene"}],
        right=[{"chebi_id": CHEBI_OXIDOSQUALENE, "name": "oxidosqualene"}],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_go_0016423_rejects_missing_modified_trna_product() -> None:
    cls = TRNAGuanineMethyltransferaseActivity()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_GUANOSINE_IN_TRNA, "name": "guanosine in tRNA"},
            {"chebi_id": CHEBI_SAM, "name": "SAM"},
        ],
        right=[
            {"chebi_id": CHEBI_GUANOSINE_IN_TRNA, "name": "guanosine in tRNA"},
            {"chebi_id": CHEBI_SAH, "name": "SAH"},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is False


@pytest.mark.parametrize(
    "class_name",
    [
        "OxidoreductaseActingOnTheCHOHGroupOfDonorsQuinoneOrSimilarCompoundAsAcceptor",
        "ProteinLysineNMethyltransferaseActivity",
        "LysineNMethyltransferaseActivity",
        "HistoneMethyltransferaseActivity",
        "HistoneModifyingActivity",
        "PseudouridineSynthaseActivity",
        "SAcyltransferaseActivity",
        "OxidosqualeneCyclaseActivity",
        "TRNAGuanineMethyltransferaseActivity",
    ],
)
def test_batch25_eval_smoke(class_name: str) -> None:
    metrics = evaluate_reaction_class(class_name, cache_dir="cache", go_only=True)
    assert metrics.true_positives >= 0
