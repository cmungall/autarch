"""Focused tests for batch 12 specific classifiers."""

from __future__ import annotations

from typing_extensions import TypedDict

from autarch.datamodel import Participant, Reaction
from autarch.molecules import (
    CHEBI_AMP,
    CHEBI_ATP,
    CHEBI_COA,
    CHEBI_DIPHOSPHATE,
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NAD_PLUS,
    CHEBI_NH4,
    CHEBI_O2,
    CHEBI_SAH,
    CHEBI_SAM,
)
from autarch.ontology.cholate_coa_ligase import CholateCoALigase
from autarch.ontology.delta3_delta2_enoyl_coa_isomerase import Delta3Delta2EnoylCoAIsomerase
from autarch.ontology.glycine_oxidase import GlycineOxidase
from autarch.ontology.l_histidine_n_alpha_methyltransferase import (
    LHistidineNAlphaMethyltransferase,
)
from autarch.ontology.three_hydroxyacyl_coa_dehydratase import (
    ThreeHydroxyacylCoADehydratase,
)

CHEBI_CROTONYL_COA = "CHEBI:57330"
CHEBI_GLYCINE = "CHEBI:57305"
CHEBI_GLYOXYLATE = "CHEBI:36655"
CHEBI_CHOLATE = "CHEBI:29747"
CHEBI_CHOLOYL_COA = "CHEBI:57373"
CHEBI_BENZOATE = "CHEBI:16150"
CHEBI_BENZOYL_COA = "CHEBI:57369"
CHEBI_DIMETHYL_HISTIDINE = "CHEBI:57610"
CHEBI_HERCYNINE = "CHEBI:15781"
CHEBI_CATECHOL = "CHEBI:17847"
CHEBI_GUAIACOL = "CHEBI:18199"


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


def test_go_0018812_positive_three_hydroxyacyl_coa_dehydratase() -> None:
    cls = ThreeHydroxyacylCoADehydratase()
    reaction = _reaction(
        left=[
            {
                "smiles": "[1*][C@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-]",
                "name": "generic 3-hydroxyacyl-CoA",
            },
        ],
        right=[
            {
                "smiles": "[1*]/C=C/C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-]",
                "name": "generic (2E)-enoyl-CoA",
            },
            {"chebi_id": CHEBI_H2O, "smiles": "[H]O[H]", "name": "water"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "dehydratase" in result.explanation.lower() or "enoyl" in result.explanation.lower()


def test_go_0018812_rejects_enoyl_reduction() -> None:
    cls = ThreeHydroxyacylCoADehydratase()
    reaction = _reaction(
        left=[
            {
                "chebi_id": CHEBI_CROTONYL_COA,
                "smiles": "CCCCCCCCC/C=C/C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-]",
                "name": "trans-dodec-2-enoyl-CoA",
            },
            {"chebi_id": CHEBI_NADH, "name": "NADH"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
        right=[
            {
                "smiles": "CCCCCCCCCCC[C@H]2CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-]",
                "name": "reduced acyl-CoA analog",
            },
            {"chebi_id": CHEBI_NAD_PLUS, "name": "NAD+"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_go_0018812_rejects_hmg_coa_dehydration() -> None:
    cls = ThreeHydroxyacylCoADehydratase()
    reaction = _reaction(
        left=[
            {
                "smiles": "CC(C)(COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-])[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)C[C@@](C)(O)CC(=O)[O-]",
                "name": "(3S)-3-hydroxy-3-methylglutaryl-CoA",
            },
        ],
        right=[
            {
                "smiles": "C/C(=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-])CC(=O)[O-]",
                "name": "3-methyl-(2E)-glutaconyl-CoA",
            },
            {"chebi_id": CHEBI_H2O, "smiles": "[H]O[H]", "name": "water"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_go_0043799_positive_glycine_oxidase() -> None:
    cls = GlycineOxidase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_GLYCINE, "smiles": "[NH3+]CC(=O)[O-]", "name": "glycine"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
            {"chebi_id": CHEBI_H2O, "smiles": "[H]O[H]", "name": "water"},
        ],
        right=[
            {"chebi_id": CHEBI_GLYOXYLATE, "smiles": "[H]C(=O)C(=O)[O-]", "name": "glyoxylate"},
            {"chebi_id": CHEBI_H2O2, "smiles": "[H]OO[H]", "name": "hydrogen peroxide"},
            {"chebi_id": CHEBI_NH4, "smiles": "[H][N+]([H])([H])[H]", "name": "ammonium"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "glycine oxidase" in result.explanation.lower() or "small amino acid" in result.explanation.lower()


def test_go_0043799_rejects_large_amino_acid_oxidation() -> None:
    cls = GlycineOxidase()
    reaction = _reaction(
        left=[
            {"smiles": "[NH3+]CCCC[C@H](N)C(=O)[O-]", "name": "lysine-like substrate"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
            {"chebi_id": CHEBI_H2O, "smiles": "[H]O[H]", "name": "water"},
        ],
        right=[
            {"smiles": "O=CCCCC(=O)[O-]", "name": "oxo acid"},
            {"chebi_id": CHEBI_H2O2, "smiles": "[H]OO[H]", "name": "hydrogen peroxide"},
            {"chebi_id": CHEBI_NH4, "smiles": "[H][N+]([H])([H])[H]", "name": "ammonium"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_go_0052706_positive_histidine_n_alpha_methyltransferase() -> None:
    cls = LHistidineNAlphaMethyltransferase()
    reaction = _reaction(
        left=[
            {
                "chebi_id": CHEBI_DIMETHYL_HISTIDINE,
                "smiles": "C[NH+](C)[C@@H](CC1=CNC=N1)C(=O)[O-]",
                "name": "dimethylhistidine",
            },
            {
                "chebi_id": CHEBI_SAM,
                "smiles": "C[S+](CC[C@H]([NH3+])C(=O)[O-])C[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1O",
                "name": "SAM",
            },
        ],
        right=[
            {
                "chebi_id": CHEBI_HERCYNINE,
                "smiles": "C[N+](C)(C)[C@@H](CC1=CNC=N1)C(=O)[O-]",
                "name": "hercynine",
            },
            {
                "chebi_id": CHEBI_SAH,
                "smiles": "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](CSCC[C@H]([NH3+])C(=O)[O-])[C@@H](O)[C@H]1O",
                "name": "SAH",
            },
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "histidine" in result.explanation.lower() or "methyltransferase" in result.explanation.lower()


def test_go_0052706_rejects_non_histidine_methyltransferase() -> None:
    cls = LHistidineNAlphaMethyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_SAM, "name": "SAM"},
            {"chebi_id": CHEBI_CATECHOL, "smiles": "Oc1ccccc1O", "name": "catechol"},
        ],
        right=[
            {"chebi_id": CHEBI_SAH, "name": "SAH"},
            {"chebi_id": CHEBI_GUAIACOL, "smiles": "COc1ccccc1O", "name": "guaiacol"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_go_0052706_rejects_carnosine_methylation() -> None:
    cls = LHistidineNAlphaMethyltransferase()
    reaction = _reaction(
        left=[
            {
                "smiles": "[NH3+]CCC(=O)N[C@@H](CC1=CNC=N1)C(=O)[O-]",
                "name": "carnosine",
            },
            {
                "chebi_id": CHEBI_SAM,
                "smiles": "C[S+](CC[C@H]([NH3+])C(=O)[O-])C[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1O",
                "name": "SAM",
            },
        ],
        right=[
            {
                "smiles": "CN1C=NC=C1C[C@H](NC(=O)CC[NH3+])C(=O)[O-]",
                "name": "anserine",
            },
            {
                "chebi_id": CHEBI_SAH,
                "smiles": "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](CSCC[C@H]([NH3+])C(=O)[O-])[C@@H](O)[C@H]1O",
                "name": "SAH",
            },
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_go_0004165_positive_delta3_delta2_enoyl_coa_isomerase() -> None:
    cls = Delta3Delta2EnoylCoAIsomerase()
    reaction = _reaction(
        left=[
            {
                "smiles": "CCCCCCCC/C=C\\CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-]",
                "name": "(3Z)-dodecenoyl-CoA",
            },
        ],
        right=[
            {
                "chebi_id": CHEBI_CROTONYL_COA,
                "smiles": "CCCCCCCCC/C=C/C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-]",
                "name": "(2E)-dodecenoyl-CoA",
            },
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "isomerase" in result.explanation.lower() or "delta(3)" in result.explanation.lower()


def test_go_0004165_rejects_enoyl_reduction() -> None:
    cls = Delta3Delta2EnoylCoAIsomerase()
    reaction = _reaction(
        left=[
            {
                "chebi_id": CHEBI_CROTONYL_COA,
                "smiles": "CCCCCCCCC/C=C/C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-]",
                "name": "(2E)-dodecenoyl-CoA",
            },
            {"chebi_id": CHEBI_NADH, "name": "NADH"},
        ],
        right=[
            {
                "smiles": "CCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-]",
                "name": "dodecanoyl-CoA analog",
            },
            {"chebi_id": CHEBI_NAD_PLUS, "name": "NAD+"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_go_0004165_rejects_acp_enoyl_isomerization() -> None:
    cls = Delta3Delta2EnoylCoAIsomerase()
    reaction = _reaction(
        left=[
            {
                "smiles": "*N[C@@H](COP(=O)([O-])OCC(C)(C)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)/C=C/CCCCCCC)C(*)=O",
                "name": "(2E)-decenoyl-[ACP]",
            },
        ],
        right=[
            {
                "smiles": "*N[C@@H](COP(=O)([O-])OCC(C)(C)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)C/C=C\\CCCCCC)C(*)=O",
                "name": "(3Z)-decenoyl-[ACP]",
            },
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_go_0047747_positive_cholate_coa_ligase() -> None:
    cls = CholateCoALigase()
    reaction = _reaction(
        left=[
            {
                "chebi_id": CHEBI_CHOLATE,
                "smiles": "[H][C@@]12C[C@H](O)CC[C@]1(C)[C@@]1([H])C[C@H](O)[C@@]3(C)[C@@]([H])(CC[C@]3([H])[C@H](C)CCC(=O)[O-])[C@]1([H])[C@H](O)C2",
                "name": "cholate",
            },
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
            {"chebi_id": CHEBI_COA, "name": "CoA"},
        ],
        right=[
            {
                "chebi_id": CHEBI_CHOLOYL_COA,
                "smiles": "[H][C@@]12C[C@H](O)CC[C@]1(C)[C@@]1([H])C[C@H](O)[C@@]3(C)[C@@]([H])(CC[C@]3([H])[C@]([H])(C)CCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]3O[C@@H](N4C=NC5=C4N=CN=C5N)[C@H](O)[C@@H]3OP(=O)([O-])[O-])[C@]1([H])[C@H](O)C2",
                "name": "choloyl-CoA",
            },
            {"chebi_id": CHEBI_AMP, "name": "AMP"},
            {"chebi_id": CHEBI_DIPHOSPHATE, "smiles": "O=P([O-])([O-])OP(=O)([O-])O", "name": "diphosphate"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "cholate" in result.explanation.lower() or "bile acid" in result.explanation.lower()


def test_go_0047747_rejects_non_bile_acid_coa_ligase() -> None:
    cls = CholateCoALigase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_BENZOATE, "smiles": "O=C([O-])C1=CC=CC=C1", "name": "benzoate"},
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
            {"chebi_id": CHEBI_COA, "name": "CoA"},
        ],
        right=[
            {"chebi_id": CHEBI_BENZOYL_COA, "smiles": "CC(C)(COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-])[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)C1=CC=CC=C1", "name": "benzoyl-CoA"},
            {"chebi_id": CHEBI_AMP, "name": "AMP"},
            {"chebi_id": CHEBI_DIPHOSPHATE, "smiles": "O=P([O-])([O-])OP(=O)([O-])O", "name": "diphosphate"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False
