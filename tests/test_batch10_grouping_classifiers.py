"""Focused tests for batch 10 grouping classifiers."""

from __future__ import annotations

from typing_extensions import TypedDict

from autarch.datamodel import Participant, Reaction
from autarch.molecules import (
    CHEBI_AMP,
    CHEBI_ATP,
    CHEBI_CMP,
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_H_PLUS,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_O2,
)
from autarch.ontology.diphosphotransferase import Diphosphotransferase
from autarch.ontology.oxidoreductase_acting_on_a_sulfur_group_of_donors_oxygen_as_acceptor import (
    OxidoreductaseActingOnASulfurGroupOfDonorsOxygenAsAcceptor,
)
from autarch.ontology.oxidoreductase_acting_on_the_aldehyde_or_oxo_group_of_donors_oxygen_as_acceptor import (
    OxidoreductaseActingOnTheAldehydeOrOxoGroupOfDonorsOxygenAsAcceptor,
)
from autarch.ontology.oxidoreductase_acting_on_the_ch_ch_group_of_donors_iron_sulfur_protein_as_acceptor import (
    OxidoreductaseActingOnTheCHCHGroupOfDonorsIronSulfurProteinAsAcceptor,
)
from autarch.ontology.sialyltransferase import Sialyltransferase

CHEBI_CMP_NEUAC = "CHEBI:57812"
CHEBI_OXIDIZED_4FE4S = "CHEBI:33722"
CHEBI_REDUCED_4FE4S = "CHEBI:33723"


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


def test_ec_137_positive_iron_sulfur_ch_ch_redox() -> None:
    cls = OxidoreductaseActingOnTheCHCHGroupOfDonorsIronSulfurProteinAsAcceptor()
    reaction = _reaction(
        left=[
            {
                "smiles": "[H][C@]1(CC2=N/C(=C\\C3=C(CCC(=O)[O-])C(C)=C(/C=C4\\NC(=O)C(C)=C4C=C)N3)C(CCC(=O)[O-])=C2C)NC(=O)C(C=C)=C1C",
                "name": "15,16-dihydrobiliverdin",
            },
            {
                "chebi_id": CHEBI_OXIDIZED_4FE4S,
                "smiles": "[SH]12[Fe]3[SH]4[Fe]1[SH]1[Fe+]2[SH]3[Fe+]41",
                "name": "oxidized [4Fe-4S] cluster",
                "count": 2,
            },
        ],
        right=[
            {
                "smiles": "C=CC1=C(C)/C(=C/C2=N/C(=C\\C3=C(CCC(=O)[O-])C(C)=C(/C=C4\\NC(=O)C(C)=C4C=C)N3)C(CCC(=O)[O-])=C2C)NC1=O",
                "name": "biliverdin-like product",
            },
            {
                "chebi_id": CHEBI_REDUCED_4FE4S,
                "smiles": "[SH]12[Fe]3[SH]4[Fe]1[SH]1[Fe]2[SH]3[Fe+]41",
                "name": "reduced [4Fe-4S] cluster",
                "count": 2,
            },
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron", "count": 2},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "iron-sulfur" in result.explanation.lower()


def test_ec_137_rejects_nicotinamide_ch_ch_redox() -> None:
    cls = OxidoreductaseActingOnTheCHCHGroupOfDonorsIronSulfurProteinAsAcceptor()
    reaction = _reaction(
        left=[
            {"smiles": "O=C1NC(=O)C=CC1", "name": "dihydroorotate-like substrate"},
            {"chebi_id": CHEBI_NADP_PLUS, "name": "NADP+"},
        ],
        right=[
            {"smiles": "O=C1NC(=O)C=CC1", "name": "orotate-like product"},
            {"chebi_id": CHEBI_NADPH, "name": "NADPH"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_123_positive_aldehyde_oxidase() -> None:
    cls = OxidoreductaseActingOnTheAldehydeOrOxoGroupOfDonorsOxygenAsAcceptor()
    reaction = _reaction(
        left=[
            {"smiles": "CC=O", "name": "acetaldehyde"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "CC(=O)[O-]", "name": "acetate"},
            {"chebi_id": CHEBI_H2O2, "smiles": "OO", "name": "hydrogen peroxide"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "oxygen" in result.explanation.lower()


def test_ec_123_rejects_alcohol_oxidase() -> None:
    cls = OxidoreductaseActingOnTheAldehydeOrOxoGroupOfDonorsOxygenAsAcceptor()
    reaction = _reaction(
        left=[
            {"smiles": "CCO", "name": "ethanol"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
        ],
        right=[
            {"smiles": "CC=O", "name": "acetaldehyde"},
            {"chebi_id": CHEBI_H2O2, "smiles": "OO", "name": "hydrogen peroxide"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_sialyltransferase_positive_cmp_neuac_transfer() -> None:
    cls = Sialyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_CMP_NEUAC, "name": "CMP-Neu5Ac"},
            {"name": "beta-galactoside acceptor"},
        ],
        right=[
            {"chebi_id": CHEBI_CMP, "name": "CMP"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
            {"name": "sialylated product"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "sialyl" in result.explanation.lower()


def test_sialyltransferase_rejects_cmp_neuac_hydrolysis() -> None:
    cls = Sialyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_CMP_NEUAC, "name": "CMP-Neu5Ac"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"chebi_id": CHEBI_CMP, "name": "CMP"},
            {"name": "free sialic acid"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_276_positive_diphosphotransferase() -> None:
    cls = Diphosphotransferase()
    reaction = _reaction(
        left=[
            {
                "smiles": "O=P([O-])([O-])OC[C@H]1OC(O)[C@H](O)[C@@H]1O",
                "name": "ribose 5-phosphate",
            },
            {
                "chebi_id": CHEBI_ATP,
                "smiles": "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]1O",
                "name": "ATP",
            },
        ],
        right=[
            {
                "smiles": "O=P([O-])([O-])OC[C@H]1O[C@H](OP(=O)([O-])OP(=O)([O-])[O-])[C@H](O)[C@@H]1O",
                "name": "ribose diphosphate product",
            },
            {
                "chebi_id": CHEBI_AMP,
                "smiles": "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(=O)([O-])[O-])[C@@H](O)[C@H]1O",
                "name": "AMP",
            },
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "diphospho" in result.explanation.lower() or "pyrophosph" in result.explanation.lower()


def test_ec_276_rejects_simple_kinase() -> None:
    cls = Diphosphotransferase()
    reaction = _reaction(
        left=[
            {"smiles": "OCC1OC(O)C(O)C(O)C1O", "name": "glucose"},
            {
                "chebi_id": CHEBI_ATP,
                "smiles": "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]1O",
                "name": "ATP",
            },
        ],
        right=[
            {"smiles": "O=P([O-])([O-])OCC1OC(O)C(O)C(O)C1O", "name": "glucose phosphate"},
            {"name": "ADP"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_183_positive_sulfite_oxidase() -> None:
    cls = OxidoreductaseActingOnASulfurGroupOfDonorsOxygenAsAcceptor()
    reaction = _reaction(
        left=[
            {"smiles": "O=S([O-])[O-]", "name": "sulfite"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "O=S(=O)([O-])[O-]", "name": "sulfate"},
            {"chebi_id": CHEBI_H2O2, "smiles": "OO", "name": "hydrogen peroxide"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "sulfur" in result.explanation.lower()


def test_ec_183_rejects_nad_linked_sulfur_redox() -> None:
    cls = OxidoreductaseActingOnASulfurGroupOfDonorsOxygenAsAcceptor()
    reaction = _reaction(
        left=[
            {"smiles": "O=S([O-])[O-]", "name": "sulfite"},
            {"chebi_id": CHEBI_NADP_PLUS, "name": "NADP+"},
        ],
        right=[
            {"smiles": "O=S(=O)([O-])[O-]", "name": "sulfate"},
            {"chebi_id": CHEBI_NADPH, "name": "NADPH"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False
