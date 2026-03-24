"""Focused tests for the next batch of EC level-3 grouping classifiers."""

from typing_extensions import TypedDict

from autarch.datamodel import Participant, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_ATP,
    CHEBI_CO2,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_O2,
)
from autarch.ontology.two_oxoglutarate_dependent_dioxygenase import (
    TwoOxoglutarateDependentDioxygenase,
)
from autarch.ontology.phosphotransferase_carboxyl_group_as_acceptor import (
    PhosphotransferaseCarboxylGroupAsAcceptor,
)
from autarch.ontology.phosphotransferase_nitrogenous_group_as_acceptor import (
    PhosphotransferaseNitrogenousGroupAsAcceptor,
)


class ParticipantInput(TypedDict, total=False):
    """Subset of participant fields used in classifier tests."""

    chebi_id: str
    smiles: str
    name: str
    count: int


def _reaction(
    left: list[ParticipantInput],
    right: list[ParticipantInput],
) -> Reaction:
    return Reaction(
        left_participants=[Participant(**participant) for participant in left],
        right_participants=[Participant(**participant) for participant in right],
    )


def test_ec_11411_positive_two_oxoglutarate_dependent_dioxygenase() -> None:
    cls = TwoOxoglutarateDependentDioxygenase()
    reaction = _reaction(
        left=[
            {"smiles": "O=C([O-])CCCC(=O)[O-]", "name": "glutarate"},
            {"chebi_id": "CHEBI:16810", "smiles": "O=C([O-])CCC(=O)C(=O)[O-]", "name": "2-oxoglutarate"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
        ],
        right=[
            {"smiles": "O=C([O-])CC[C@H](O)C(=O)[O-]", "name": "(S)-2-hydroxyglutarate"},
            {"chebi_id": "CHEBI:30031", "smiles": "O=C([O-])CCC(=O)[O-]", "name": "succinate"},
            {"chebi_id": CHEBI_CO2, "smiles": "O=C=O", "name": "carbon dioxide"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "2-oxoglutarate" in result.explanation.lower()


def test_ec_11411_rejects_nadph_monooxygenase() -> None:
    cls = TwoOxoglutarateDependentDioxygenase()
    reaction = _reaction(
        left=[
            {"smiles": "c1ccccc1", "name": "benzene"},
            {"chebi_id": CHEBI_NADPH, "name": "NADPH"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
        ],
        right=[
            {"smiles": "Oc1ccccc1", "name": "phenol"},
            {"chebi_id": CHEBI_NADP_PLUS, "name": "NADP+"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_272_positive_carboxyl_group_acceptor() -> None:
    cls = PhosphotransferaseCarboxylGroupAsAcceptor()
    reaction = _reaction(
        left=[
            {"smiles": "CC(=O)[O-]", "name": "acetate"},
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
        ],
        right=[
            {"smiles": "CC(=O)OP(=O)([O-])[O-]", "name": "acetyl phosphate"},
            {"chebi_id": CHEBI_ADP, "name": "ADP"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "carboxyl" in result.explanation.lower()


def test_ec_272_rejects_alcohol_group_acceptor() -> None:
    cls = PhosphotransferaseCarboxylGroupAsAcceptor()
    reaction = _reaction(
        left=[
            {"smiles": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O", "name": "glucose"},
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
        ],
        right=[
            {"smiles": "O=P([O-])([O-])OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O", "name": "glucose 6-phosphate"},
            {"chebi_id": CHEBI_ADP, "name": "ADP"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_273_positive_nitrogenous_group_acceptor() -> None:
    cls = PhosphotransferaseNitrogenousGroupAsAcceptor()
    reaction = _reaction(
        left=[
            {"smiles": "CN(CC(=O)[O-])C(N)=[NH2+]", "name": "creatine"},
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
        ],
        right=[
            {"smiles": "CN(CC(=O)[O-])C(=[NH2+])NP(=O)([O-])[O-]", "name": "N-phosphocreatine"},
            {"chebi_id": CHEBI_ADP, "name": "ADP"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "nitrogenous" in result.explanation.lower()


def test_ec_273_rejects_carboxyl_group_acceptor() -> None:
    cls = PhosphotransferaseNitrogenousGroupAsAcceptor()
    reaction = _reaction(
        left=[
            {"smiles": "CC(=O)[O-]", "name": "acetate"},
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
        ],
        right=[
            {"smiles": "CC(=O)OP(=O)([O-])[O-]", "name": "acetyl phosphate"},
            {"chebi_id": CHEBI_ADP, "name": "ADP"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False
