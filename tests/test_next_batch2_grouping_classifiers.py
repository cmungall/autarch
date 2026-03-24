"""Focused tests for the next batch of EC level-3 grouping classifiers."""

from typing_extensions import TypedDict

from autarch.datamodel import Participant, PolymerType, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_ATP,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADP_PLUS,
    CHEBI_NADPH,
    CHEBI_NAD_PLUS,
    CHEBI_NH3,
    CHEBI_NH4,
    CHEBI_PHOSPHATE,
)
from autarch.ontology.aminoacyltransferase import Aminoacyltransferase
from autarch.ontology.hydrolase_acting_on_acid_anhydrides_in_phosphorus_containing_anhydrides import (
    HydrolaseActingOnAcidAnhydridesInPhosphorusContainingAnhydrides,
)
from autarch.ontology.oxidoreductase_acting_on_other_nitrogenous_compounds_as_donors_with_nad_or_nadp_as_acceptor import (
    OxidoreductaseActingOnOtherNitrogenousCompoundsAsDonorsWithNADOrNADPAsAcceptor,
)


class ParticipantInput(TypedDict, total=False):
    """Subset of participant fields used in classifier tests."""

    chebi_id: str
    smiles: str
    name: str
    count: int
    polymer_type: PolymerType


ATP_SMILES = (
    "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])"
    "[C@@H](O)[C@H]1O"
)
ADP_SMILES = (
    "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]1O"
)
AMP_SMILES = "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(=O)([O-])[O-])[C@@H](O)[C@H]1O"
PHOSPHATE_SMILES = "O=P([O-])([O-])[O-]"
DIPHOSPHATE_SMILES = "O=P([O-])([O-])OP(=O)([O-])[O-]"


def _reaction(
    left: list[ParticipantInput],
    right: list[ParticipantInput],
) -> Reaction:
    return Reaction(
        left_participants=[Participant(**participant) for participant in left],
        right_participants=[Participant(**participant) for participant in right],
    )


def test_ec_171_positive_nad_linked_nitrogenous_redox() -> None:
    cls = OxidoreductaseActingOnOtherNitrogenousCompoundsAsDonorsWithNADOrNADPAsAcceptor()
    reaction = _reaction(
        left=[
            {"chebi_id": "CHEBI:16301", "smiles": "[O-][N]=O", "name": "nitrite"},
            {"chebi_id": CHEBI_NAD_PLUS, "name": "NAD(+)"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {
                "chebi_id": "CHEBI:17632",
                "smiles": "O=[N+]([O-])[O-]",
                "name": "nitrate",
            },
            {"chebi_id": CHEBI_NADH, "name": "NADH"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "nad" in result.explanation.lower()


def test_ec_171_rejects_ch_nh_redox() -> None:
    cls = OxidoreductaseActingOnOtherNitrogenousCompoundsAsDonorsWithNADOrNADPAsAcceptor()
    reaction = _reaction(
        left=[
            {
                "smiles": "[H][C@@]1(C(=O)[O-])CCCC[NH2+]1",
                "name": "L-pipecolate",
            },
            {"chebi_id": CHEBI_NADP_PLUS, "name": "NADP(+)"},
        ],
        right=[
            {
                "smiles": "O=C([O-])C1=[NH+]CCCC1",
                "name": "Delta(1)-piperideine-2-carboxylate",
            },
            {"chebi_id": CHEBI_NADPH, "name": "NADPH"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_361_positive_atp_hydrolysis() -> None:
    cls = HydrolaseActingOnAcidAnhydridesInPhosphorusContainingAnhydrides()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_ATP, "smiles": ATP_SMILES, "name": "ATP"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"chebi_id": CHEBI_ADP, "smiles": ADP_SMILES, "name": "ADP"},
            {"chebi_id": CHEBI_PHOSPHATE, "smiles": PHOSPHATE_SMILES, "name": "phosphate"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "anhydride" in result.explanation.lower()


def test_ec_361_positive_acyl_phosphate_hydrolysis() -> None:
    cls = HydrolaseActingOnAcidAnhydridesInPhosphorusContainingAnhydrides()
    reaction = _reaction(
        left=[
            {"smiles": "CC(=O)OP(=O)([O-])[O-]", "name": "acetyl phosphate"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "CC(=O)[O-]", "name": "acetate"},
            {"chebi_id": CHEBI_PHOSPHATE, "smiles": PHOSPHATE_SMILES, "name": "phosphate"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_361_rejects_phosphomonoester_hydrolysis() -> None:
    cls = HydrolaseActingOnAcidAnhydridesInPhosphorusContainingAnhydrides()
    reaction = _reaction(
        left=[
            {
                "smiles": "O=P([O-])([O-])OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
                "name": "glucose 6-phosphate",
            },
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O", "name": "glucose"},
            {"chebi_id": CHEBI_PHOSPHATE, "smiles": PHOSPHATE_SMILES, "name": "phosphate"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_232_positive_aminoacyl_trna_transfer() -> None:
    cls = Aminoacyltransferase()
    reaction = _reaction(
        left=[
            {
                "chebi_id": "CHEBI:78597",
                "name": "an N-terminal L-alpha-aminoacyl-[protein]",
                "polymer_type": PolymerType.PROTEIN,
                "smiles": "*C(=O)[C@H](*)[NH3+]",
            },
            {
                "chebi_id": "CHEBI:78513",
                "name": "L-arginyl-tRNA(Arg)",
                "polymer_type": PolymerType.TRNA,
                "smiles": "*P(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OC(=O)[C@@H]([NH3+])CCCNC(N)=[NH2+]",
            },
        ],
        right=[
            {
                "chebi_id": "CHEBI:83562",
                "name": "an N-terminal L-arginyl-L-aminoacyl-[protein]",
                "polymer_type": PolymerType.PROTEIN,
                "smiles": "*C(=O)[C@H](*)NC(=O)[C@@H]([NH3+])CCCNC(N)=[NH2+]",
            },
            {
                "chebi_id": "CHEBI:78442",
                "name": "tRNA(Arg)",
                "polymer_type": PolymerType.TRNA,
                "smiles": "*P(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C(N)N=CN=C32)[C@H](O)[C@@H]1O",
            },
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "aminoacyl" in result.explanation.lower()


def test_ec_232_positive_glutamyl_transfer() -> None:
    cls = Aminoacyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": "CHEBI:57416", "smiles": "C[C@@H]([NH3+])C(=O)[O-]", "name": "D-alanine"},
            {"chebi_id": "CHEBI:58359", "smiles": "NC(=O)CC[C@H]([NH3+])C(=O)[O-]", "name": "L-glutamine"},
        ],
        right=[
            {
                "chebi_id": "CHEBI:57915",
                "smiles": "C[C@@H](NC(=O)CC[C@H]([NH3+])C(=O)[O-])C(=O)[O-]",
                "name": "gamma-L-glutamyl-D-alanine",
            },
            {"chebi_id": CHEBI_NH4, "smiles": "[H][N+]([H])([H])[H]", "name": "ammonium"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_232_rejects_atp_dependent_amide_ligase() -> None:
    cls = Aminoacyltransferase()
    reaction = _reaction(
        left=[
            {"smiles": "[NH3+]CCC[C@H]([NH3+])C(=O)[O-]", "name": "L-lysine"},
            {"smiles": "O=C([O-])CC[C@H]([NH3+])C(=O)[O-]", "name": "L-glutamate"},
            {"chebi_id": CHEBI_ATP, "smiles": ATP_SMILES, "name": "ATP"},
            {"chebi_id": CHEBI_NH3, "smiles": "N", "name": "ammonia"},
        ],
        right=[
            {"smiles": "[NH3+]CCC[C@H](NC(=O)CC[C@H]([NH3+])C(=O)[O-])C(=O)[O-]", "name": "glutamyl-lysine"},
            {"chebi_id": CHEBI_ADP, "smiles": ADP_SMILES, "name": "ADP"},
            {"chebi_id": CHEBI_PHOSPHATE, "smiles": PHOSPHATE_SMILES, "name": "phosphate"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False
