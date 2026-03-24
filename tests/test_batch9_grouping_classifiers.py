"""Focused tests for batch 9 grouping classifiers."""

from __future__ import annotations

from typing_extensions import TypedDict

from autarch.datamodel import Participant, Reaction
from autarch.molecules import (
    CHEBI_DIPHOSPHATE,
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_H_PLUS,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_O2,
    CHEBI_PHOSPHATE,
)
from autarch.ontology.acyltransferase_acyl_groups_converted_into_alkyl_on_transfer import (
    AcyltransferaseAcylGroupsConvertedIntoAlkylOnTransfer,
)
from autarch.ontology.carboxyl_or_carbamoyltransferase import (
    CarboxylOrCarbamoyltransferase,
)
from autarch.ontology.diphosphoric_monoester_hydrolase import (
    DiphosphoricMonoesterHydrolase,
)
from autarch.ontology.oxidoreductase_acting_on_ch_or_ch2_groups_disulfide_as_acceptor import (
    OxidoreductaseActingOnCHOrCH2GroupsDisulfideAsAcceptor,
)
from autarch.ontology.oxidoreductase_acting_on_diphenols_and_related_substances_as_donors_oxygen_as_acceptor import (
    OxidoreductaseActingOnDiphenolsAndRelatedSubstancesAsDonorsOxygenAsAcceptor,
)


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


def test_ec_1174_positive_disulfide_acceptor_ch_redox() -> None:
    cls = OxidoreductaseActingOnCHOrCH2GroupsDisulfideAsAcceptor()
    reaction = _reaction(
        left=[
            {"smiles": "CSSC", "name": "dimethyl disulfide"},
            {
                "smiles": "NC1=NC=NC2=C1N=CN2[C@H]1C[C@H](O)[C@@H](COP(=O)([O-])OP(=O)([O-])[O-])O1",
                "name": "dADP",
            },
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "CS", "name": "methanethiol", "count": 2},
            {
                "smiles": "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]1O",
                "name": "ADP",
            },
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "disulfide" in result.explanation.lower()


def test_ec_1174_rejects_sulfur_group_redox_branch() -> None:
    cls = OxidoreductaseActingOnCHOrCH2GroupsDisulfideAsAcceptor()
    reaction = _reaction(
        left=[
            {"smiles": "CSSC", "name": "dimethyl disulfide"},
            {"smiles": "CS", "name": "methanethiol", "count": 2},
        ],
        right=[
            {"smiles": "CS", "name": "methanethiol", "count": 2},
            {"smiles": "CSSC", "name": "dimethyl disulfide"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_1103_positive_diphenol_oxidase() -> None:
    cls = OxidoreductaseActingOnDiphenolsAndRelatedSubstancesAsDonorsOxygenAsAcceptor()
    reaction = _reaction(
        left=[
            {"smiles": "Oc1ccc(O)cc1", "name": "hydroquinone"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
        ],
        right=[
            {"smiles": "O=C1C=CC(=O)C=C1", "name": "benzoquinone"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "diphenol" in result.explanation.lower()


def test_ec_1103_rejects_alcohol_oxidase() -> None:
    cls = OxidoreductaseActingOnDiphenolsAndRelatedSubstancesAsDonorsOxygenAsAcceptor()
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


def test_ec_213_positive_carbamoyl_transfer() -> None:
    cls = CarboxylOrCarbamoyltransferase()
    reaction = _reaction(
        left=[
            {"smiles": "NC(=O)OP(=O)(O)O", "name": "carbamoyl phosphate"},
            {"smiles": "NCCCC[C@H](N)C(=O)O", "name": "ornithine"},
        ],
        right=[
            {"smiles": "NC(=O)NCCCC[C@H](N)C(=O)O", "name": "citrulline"},
            {"chebi_id": CHEBI_PHOSPHATE, "smiles": "OP(=O)(O)O", "name": "phosphate"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "carbamoyl" in result.explanation.lower()


def test_ec_213_rejects_simple_acyl_transfer() -> None:
    cls = CarboxylOrCarbamoyltransferase()
    reaction = _reaction(
        left=[
            {"smiles": "CC(=O)SCC", "name": "acetyl thioester"},
            {"smiles": "NCC(=O)O", "name": "glycine"},
        ],
        right=[
            {"smiles": "CC(=O)NCC(=O)O", "name": "N-acetylglycine"},
            {"smiles": "SCC", "name": "thiol carrier"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_317_positive_diphosphoester_hydrolysis() -> None:
    cls = DiphosphoricMonoesterHydrolase()
    reaction = _reaction(
        left=[
            {"smiles": "C=CCOP(=O)(O)OP(=O)(O)O", "name": "allyl diphosphate"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "C=CCO", "name": "allyl alcohol"},
            {
                "chebi_id": CHEBI_DIPHOSPHATE,
                "smiles": "OP(=O)(O)OP(=O)(O)O",
                "name": "diphosphate",
            },
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "diphosphate" in result.explanation.lower()


def test_ec_317_rejects_pyrophosphatase() -> None:
    cls = DiphosphoricMonoesterHydrolase()
    reaction = _reaction(
        left=[
            {
                "chebi_id": CHEBI_DIPHOSPHATE,
                "smiles": "OP(=O)(O)OP(=O)(O)O",
                "name": "diphosphate",
            },
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"chebi_id": CHEBI_PHOSPHATE, "smiles": "OP(=O)(O)O", "name": "phosphate", "count": 2},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_233_positive_acyl_to_alkyl_transfer() -> None:
    cls = AcyltransferaseAcylGroupsConvertedIntoAlkylOnTransfer()
    reaction = _reaction(
        left=[
            {"smiles": "CC(=O)SCC", "name": "acetyl thioester"},
            {"smiles": "O=C(O)CC(=O)C(=O)O", "name": "oxaloacetate"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "SCC", "name": "thiol carrier"},
            {"smiles": "O=C(O)CC(O)(CC(=O)O)C(=O)O", "name": "citrate"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "acyltransferase" in result.explanation.lower()


def test_ec_233_rejects_standard_acyltransferase() -> None:
    cls = AcyltransferaseAcylGroupsConvertedIntoAlkylOnTransfer()
    reaction = _reaction(
        left=[
            {"smiles": "CC(=O)SCC", "name": "acetyl thioester"},
            {"smiles": "NCC(=O)O", "name": "glycine"},
        ],
        right=[
            {"smiles": "CC(=O)NCC(=O)O", "name": "N-acetylglycine"},
            {"smiles": "SCC", "name": "thiol carrier"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_1103_rejects_nadph_quinone_redox() -> None:
    cls = OxidoreductaseActingOnDiphenolsAndRelatedSubstancesAsDonorsOxygenAsAcceptor()
    reaction = _reaction(
        left=[
            {"smiles": "O=C1C=CC(=O)C=C1", "name": "benzoquinone"},
            {"chebi_id": CHEBI_NADPH, "name": "NADPH"},
        ],
        right=[
            {"smiles": "Oc1ccc(O)cc1", "name": "hydroquinone"},
            {"chebi_id": CHEBI_NADP_PLUS, "name": "NADP(+)"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False
