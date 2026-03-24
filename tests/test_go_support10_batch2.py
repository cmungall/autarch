"""Focused tests for the next support-qualified GO classes."""

from __future__ import annotations

from typing_extensions import TypedDict

from autarch.datamodel import Participant, PolymerType, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_CTP,
    CHEBI_DIPHOSPHATE,
    CHEBI_GTP,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_UTP,
)
from autarch.ontology.acyl_coa_hydrolase import AcylCoAHydrolase
from autarch.ontology.alcohol_dehydrogenase_nad_p import AlcoholDehydrogenaseNADP
from autarch.ontology.aldehyde_dehydrogenase_nad import AldehydeDehydrogenaseNAD
from autarch.ontology.amino_acid_racemase import AminoAcidRacemase
from autarch.ontology.amino_acid_transmembrane_transporter_activity import (
    AminoAcidTransmembraneTransporterActivity,
)
from autarch.ontology.cytidylyltransferase import Cytidylyltransferase
from autarch.ontology.guanylyltransferase import Guanylyltransferase
from autarch.ontology.uridylyltransferase import Uridylyltransferase

CHEBI_ACETYL_COA = "CHEBI:57288"
CHEBI_COA = "CHEBI:57287"
CHEBI_METHIONINE = "CHEBI:57844"
CHEBI_TRIGLYCERIDE = "CHEBI:17855"
CHEBI_ACETATE = "CHEBI:30089"


class ParticipantInput(TypedDict, total=False):
    chebi_id: str
    smiles: str
    name: str
    count: int
    location: str
    polymer_type: PolymerType
    polymer_index: str


def _reaction(left: list[ParticipantInput], right: list[ParticipantInput]) -> Reaction:
    return Reaction(
        left_participants=[Participant(**participant) for participant in left],
        right_participants=[Participant(**participant) for participant in right],
    )


def _atp_coupled_transport(left_substrate: ParticipantInput, right_substrate: ParticipantInput) -> Reaction:
    return _reaction(
        left=[
            left_substrate,
            {"chebi_id": "CHEBI:15422", "smiles": "P", "name": "ATP"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            right_substrate,
            {"chebi_id": CHEBI_ADP, "smiles": "P", "name": "ADP"},
            {"chebi_id": CHEBI_DIPHOSPHATE, "smiles": "OP(=O)(O)OP(=O)(O)O", "name": "diphosphate"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )


def test_go_0047661_positive_amino_acid_racemase() -> None:
    cls = AminoAcidRacemase()
    reaction = _reaction(
        left=[{"smiles": "C[C@H]([NH3+])C(=O)[O-]", "name": "L-alanine"}],
        right=[{"smiles": "C[C@@H]([NH3+])C(=O)[O-]", "name": "D-alanine"}],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_go_0015171_positive_amino_acid_transmembrane_transporter_activity() -> None:
    cls = AminoAcidTransmembraneTransporterActivity()
    reaction = _atp_coupled_transport(
        {
            "chebi_id": CHEBI_METHIONINE,
            "smiles": "CSCC[C@H]([NH3+])C(=O)[O-]",
            "name": "L-methionine",
            "location": "out",
        },
        {
            "chebi_id": CHEBI_METHIONINE,
            "smiles": "CSCC[C@H]([NH3+])C(=O)[O-]",
            "name": "L-methionine",
            "location": "in",
        },
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_go_0015171_rejects_carbohydrate_transport() -> None:
    cls = AminoAcidTransmembraneTransporterActivity()
    reaction = _atp_coupled_transport(
        {
            "chebi_id": "CHEBI:53455",
            "smiles": "OC[C@H]1O[C@H](O)[C@@H](O)[C@H](O)C1",
            "name": "D-xylopyranose",
            "location": "out",
        },
        {
            "chebi_id": "CHEBI:53455",
            "smiles": "OC[C@H]1O[C@H](O)[C@@H](O)[C@H](O)C1",
            "name": "D-xylopyranose",
            "location": "in",
        },
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_go_0004029_positive_aldehyde_dehydrogenase_nad() -> None:
    cls = AldehydeDehydrogenaseNAD()
    reaction = _reaction(
        left=[
            {"smiles": "CC=O", "name": "acetaldehyde"},
            {"smiles": "OP(=O)([O-])[O-]", "name": "phosphate"},
            {"chebi_id": CHEBI_NAD_PLUS, "name": "NAD+"},
        ],
        right=[
            {"smiles": "CC(=O)OP(=O)([O-])[O-]", "name": "acetyl phosphate"},
            {"chebi_id": CHEBI_NADH, "name": "NADH"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_go_0004029_rejects_nadp_branch() -> None:
    cls = AldehydeDehydrogenaseNAD()
    reaction = _reaction(
        left=[
            {"smiles": "CC=O", "name": "acetaldehyde"},
            {"smiles": "OP(=O)([O-])[O-]", "name": "phosphate"},
            {"chebi_id": CHEBI_NADP_PLUS, "name": "NADP+"},
        ],
        right=[
            {"smiles": "CC(=O)OP(=O)([O-])[O-]", "name": "acetyl phosphate"},
            {"chebi_id": CHEBI_NADPH, "name": "NADPH"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_go_0008106_positive_alcohol_dehydrogenase_nadp() -> None:
    cls = AlcoholDehydrogenaseNADP()
    reaction = _reaction(
        left=[
            {"smiles": "CCO", "name": "ethanol"},
            {"chebi_id": CHEBI_NADP_PLUS, "name": "NADP+"},
        ],
        right=[
            {"smiles": "CC=O", "name": "acetaldehyde"},
            {"chebi_id": CHEBI_NADPH, "name": "NADPH"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_go_0008106_rejects_nad_branch() -> None:
    cls = AlcoholDehydrogenaseNADP()
    reaction = _reaction(
        left=[
            {"smiles": "CCO", "name": "ethanol"},
            {"chebi_id": CHEBI_NAD_PLUS, "name": "NAD+"},
        ],
        right=[
            {"smiles": "CC=O", "name": "acetaldehyde"},
            {"chebi_id": CHEBI_NADH, "name": "NADH"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_go_0016289_positive_acyl_coa_hydrolase() -> None:
    cls = AcylCoAHydrolase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_ACETYL_COA, "smiles": "CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc23)C(O)C1OP(=O)(O)O", "name": "acetyl-CoA"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"chebi_id": CHEBI_ACETATE, "smiles": "CC(=O)[O-]", "name": "acetate"},
            {"chebi_id": CHEBI_COA, "smiles": "SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc23)C(O)C1OP(=O)(O)O", "name": "CoA"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def _polymer_extension_reaction(donor_chebi: str) -> Reaction:
    return _reaction(
        left=[
            {"name": "RNA", "polymer_type": PolymerType.RNA, "polymer_index": "n"},
            {"chebi_id": donor_chebi, "name": "nucleotide donor"},
        ],
        right=[
            {"name": "RNA", "polymer_type": PolymerType.RNA, "polymer_index": "n+1"},
            {"chebi_id": CHEBI_DIPHOSPHATE, "name": "diphosphate"},
        ],
    )


def test_go_0070567_positive_cytidylyltransferase() -> None:
    cls = Cytidylyltransferase()
    assert cls.check_membership_impl(_polymer_extension_reaction(CHEBI_CTP)).is_member is True
    assert cls.check_membership_impl(_polymer_extension_reaction(CHEBI_GTP)).is_member is False


def test_go_0070568_positive_guanylyltransferase() -> None:
    cls = Guanylyltransferase()
    assert cls.check_membership_impl(_polymer_extension_reaction(CHEBI_GTP)).is_member is True
    assert cls.check_membership_impl(_polymer_extension_reaction(CHEBI_CTP)).is_member is False


def test_go_0070569_positive_uridylyltransferase() -> None:
    cls = Uridylyltransferase()
    assert cls.check_membership_impl(_polymer_extension_reaction(CHEBI_UTP)).is_member is True
    assert cls.check_membership_impl(_polymer_extension_reaction(CHEBI_CTP)).is_member is False
