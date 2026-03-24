"""Focused tests for batch 13 specific classifiers."""

from __future__ import annotations

from typing_extensions import TypedDict

from autarch.datamodel import Participant, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_ATP,
    CHEBI_DIPHOSPHATE,
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_H_PLUS,
    CHEBI_O2,
    CHEBI_PHOSPHATE,
    CHEBI_SAH,
    CHEBI_SAM,
)
from autarch.ontology.carnosine_synthase import CarnosineSynthase
from autarch.ontology.methyl_ethyl_malonyl_coa_decarboxylase import (
    MethylEthylMalonylCoADecarboxylase,
)
from autarch.ontology.n1_acetylpolyamine_oxidase import (
    N1AcetylpolyamineOxidase,
)
from autarch.ontology.pinene_synthase import PineneSynthase
from autarch.ontology.tocopherol_c_methyltransferase import TocopherolCMethyltransferase

CHEBI_BETA_ALANINE = "CHEBI:57966"
CHEBI_HISTIDINE = "CHEBI:57595"
CHEBI_CARNOSINE = "CHEBI:57485"
CHEBI_GLYCINE = "CHEBI:57305"
CHEBI_BETA_ALANYL_GLYCINE = "CHEBI:137478"

CHEBI_METHYLMALONYL_COA = "CHEBI:57326"
CHEBI_ETHYLMALONYL_COA = "CHEBI:60909"
CHEBI_PROPIONYL_COA = "CHEBI:57392"
CHEBI_BUTYRYL_COA = "CHEBI:57371"
CHEBI_CO2 = "CHEBI:16526"
CHEBI_MALONYL_COA = "CHEBI:15525"
CHEBI_ACETYL_COA = "CHEBI:15351"

CHEBI_GPP = "CHEBI:58057"
CHEBI_ALPHA_PINENE_MINUS = "CHEBI:28660"
CHEBI_LIMONENE = "CHEBI:15384"

CHEBI_GAMMA_TOCOPHEROL = "CHEBI:18185"
CHEBI_ALPHA_TOCOPHEROL = "CHEBI:18145"
CHEBI_CATECHOL = "CHEBI:17847"
CHEBI_GUAIACOL = "CHEBI:18199"

CHEBI_N1_ACETYLSPERMINE = "CHEBI:58101"
CHEBI_N1_ACETYLSPERMIDINE = "CHEBI:58324"
CHEBI_DIACETYLSPERMINE = "CHEBI:58550"
CHEBI_3_ACETAMIDOPROPANAL = "CHEBI:30322"
CHEBI_SPERMIDINE = "CHEBI:57834"
CHEBI_PUTRESCINE = "CHEBI:326268"
CHEBI_SPERMINE = "CHEBI:45725"
CHEBI_APABUTANAL = "CHEBI:58869"


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


def test_go_0047730_positive_carnosine_synthase() -> None:
    cls = CarnosineSynthase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_BETA_ALANINE, "smiles": "[NH3+]CCC(=O)[O-]", "name": "beta-alanine"},
            {
                "chebi_id": CHEBI_HISTIDINE,
                "smiles": "[NH3+][C@@H](CC1=CNC=N1)C(=O)[O-]",
                "name": "L-histidine",
            },
            {
                "chebi_id": CHEBI_ATP,
                "smiles": "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]1O",
                "name": "ATP",
            },
        ],
        right=[
            {
                "chebi_id": CHEBI_CARNOSINE,
                "smiles": "[NH3+]CCC(=O)N[C@@H](CC1=CNC=N1)C(=O)[O-]",
                "name": "carnosine",
            },
            {"chebi_id": CHEBI_ADP, "name": "ADP"},
            {"chebi_id": CHEBI_PHOSPHATE, "smiles": "O=P([O-])([O-])O", "name": "phosphate"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "carnosine synthase" in result.explanation.lower() or "beta-alanine" in result.explanation.lower()


def test_go_0047730_rejects_other_dipeptide_ligation() -> None:
    cls = CarnosineSynthase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_BETA_ALANINE, "smiles": "[NH3+]CCC(=O)[O-]", "name": "beta-alanine"},
            {"chebi_id": CHEBI_GLYCINE, "smiles": "[NH3+]CC(=O)[O-]", "name": "glycine"},
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
        ],
        right=[
            {
                "chebi_id": CHEBI_BETA_ALANYL_GLYCINE,
                "smiles": "[NH3+]CCC(=O)NCC(=O)[O-]",
                "name": "beta-alanyl-glycine",
            },
            {"chebi_id": CHEBI_ADP, "name": "ADP"},
            {"chebi_id": CHEBI_PHOSPHATE, "name": "phosphate"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_go_0004492_positive_methyl_ethyl_malonyl_coa_decarboxylase() -> None:
    cls = MethylEthylMalonylCoADecarboxylase()
    reaction = _reaction(
        left=[
            {
                "chebi_id": CHEBI_METHYLMALONYL_COA,
                "smiles": "C[C@H](C(=O)[O-])C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-]",
                "name": "methylmalonyl-CoA",
            },
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
        right=[
            {
                "chebi_id": CHEBI_PROPIONYL_COA,
                "smiles": "CCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-]",
                "name": "propionyl-CoA",
            },
            {"chebi_id": CHEBI_CO2, "smiles": "O=C=O", "name": "carbon dioxide"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "decarboxylase" in result.explanation.lower() or "methylmalonyl" in result.explanation.lower()


def test_go_0004492_rejects_plain_malonyl_coa_decarboxylation() -> None:
    cls = MethylEthylMalonylCoADecarboxylase()
    reaction = _reaction(
        left=[
            {
                "chebi_id": CHEBI_MALONYL_COA,
                "smiles": "[O-]C(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-]",
                "name": "malonyl-CoA",
            },
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
        right=[
            {
                "chebi_id": CHEBI_ACETYL_COA,
                "smiles": "CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-]",
                "name": "acetyl-CoA",
            },
            {"chebi_id": CHEBI_CO2, "smiles": "O=C=O", "name": "carbon dioxide"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_go_0050550_positive_pinene_synthase() -> None:
    cls = PineneSynthase()
    reaction = _reaction(
        left=[
            {
                "chebi_id": CHEBI_GPP,
                "smiles": "CC(C)=CCC/C(C)=C/COP(=O)([O-])OP(=O)([O-])[O-]",
                "name": "geranyl diphosphate",
            },
        ],
        right=[
            {
                "chebi_id": CHEBI_ALPHA_PINENE_MINUS,
                "smiles": "CC1=CC[C@H]2C[C@@H]1C2(C)C",
                "name": "alpha-pinene",
            },
            {"chebi_id": CHEBI_DIPHOSPHATE, "smiles": "O=P([O-])([O-])OP(=O)([O-])O", "name": "diphosphate"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "pinene synthase" in result.explanation.lower() or "pinene" in result.explanation.lower()


def test_go_0050550_rejects_non_pinene_monoterpene_cyclization() -> None:
    cls = PineneSynthase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_GPP, "name": "geranyl diphosphate"},
        ],
        right=[
            {
                "chebi_id": CHEBI_LIMONENE,
                "smiles": "CC1=CCC(CC1)C(=C)C",
                "name": "limonene",
            },
            {"chebi_id": CHEBI_DIPHOSPHATE, "name": "diphosphate"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_go_0050342_positive_tocopherol_c_methyltransferase() -> None:
    cls = TocopherolCMethyltransferase()
    reaction = _reaction(
        left=[
            {
                "chebi_id": CHEBI_GAMMA_TOCOPHEROL,
                "smiles": "CC1=C(O)\\C=C2\\CC[C@@](C)(CCC[C@H](C)CCC[C@H](C)CCCC(C)C)O\\C2=C\\1C",
                "name": "gamma-tocopherol",
            },
            {"chebi_id": CHEBI_SAM, "name": "SAM"},
        ],
        right=[
            {
                "chebi_id": CHEBI_ALPHA_TOCOPHEROL,
                "smiles": "CC1=C2\\CC[C@@](C)(CCC[C@H](C)CCC[C@H](C)CCCC(C)C)O\\C2=C(C)\\C(C)=C\\1O",
                "name": "alpha-tocopherol",
            },
            {"chebi_id": CHEBI_SAH, "name": "SAH"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "tocopherol" in result.explanation.lower() or "methyltransferase" in result.explanation.lower()


def test_go_0050342_rejects_generic_o_methyltransferase() -> None:
    cls = TocopherolCMethyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_CATECHOL, "smiles": "Oc1ccccc1O", "name": "catechol"},
            {"chebi_id": CHEBI_SAM, "name": "SAM"},
        ],
        right=[
            {"chebi_id": CHEBI_GUAIACOL, "smiles": "COc1ccccc1O", "name": "guaiacol"},
            {"chebi_id": CHEBI_SAH, "name": "SAH"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_go_0052903_positive_n_acetylpolyamine_oxidase() -> None:
    cls = N1AcetylpolyamineOxidase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_N1_ACETYLSPERMINE, "smiles": "CC(=O)NCCC[NH2+]CCCC[NH2+]CCC[NH3+]", "name": "N(1)-acetylspermine"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
            {"chebi_id": CHEBI_H2O, "smiles": "[H]O[H]", "name": "water"},
        ],
        right=[
            {"chebi_id": CHEBI_3_ACETAMIDOPROPANAL, "smiles": "[H]C(=O)CCNC(C)=O", "name": "3-acetamidopropanal"},
            {"chebi_id": CHEBI_SPERMIDINE, "smiles": "[NH3+]CCCC[NH2+]CCC[NH3+]", "name": "spermidine"},
            {"chebi_id": CHEBI_H2O2, "smiles": "[H]OO[H]", "name": "hydrogen peroxide"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "acetylpolyamine oxidase" in result.explanation.lower() or "acetamidopropanal" in result.explanation.lower()


def test_go_0052903_rejects_spermine_oxidase_branch() -> None:
    cls = N1AcetylpolyamineOxidase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_SPERMINE, "smiles": "[NH3+]CCC[NH2+]CCCC[NH2+]CCC[NH3+]", "name": "spermine"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
            {"chebi_id": CHEBI_H2O, "smiles": "[H]O[H]", "name": "water"},
        ],
        right=[
            {"chebi_id": CHEBI_APABUTANAL, "smiles": "[NH3+]CCCNCCCC=O", "name": "N-(3-aminopropyl)-4-aminobutanal"},
            {"chebi_id": CHEBI_PUTRESCINE, "smiles": "[NH3+]CCCC[NH3+]", "name": "putrescine"},
            {"chebi_id": CHEBI_H2O2, "smiles": "[H]OO[H]", "name": "hydrogen peroxide"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False
