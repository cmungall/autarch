"""Focused tests for cache-supported exact classifiers."""

from __future__ import annotations

from typing_extensions import TypedDict

from autarch.datamodel import Participant, PolymerType, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_ATP,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
)
from autarch.ontology.aldose_1_dehydrogenase_nad_p import Aldose1DehydrogenaseNADP
from autarch.ontology.d_aminoacyl_trna_deacylase import DAminoacylTRNADeacylase
from autarch.ontology.dihydrocarveol_dehydrogenase import DihydrocarveolDehydrogenase
from autarch.ontology.nucleoside_diphosphate_kinase import NucleosideDiphosphateKinase
from autarch.ontology.purine_nucleoside_phosphorylase import PurineNucleosidePhosphorylase

CHEBI_UDP = "CHEBI:58223"
CHEBI_UTP = "CHEBI:46398"
CHEBI_AMP = "CHEBI:16027"
CHEBI_GUANOSINE = "CHEBI:16750"
CHEBI_GUANINE = "CHEBI:16235"
CHEBI_PHOSPHATE = "CHEBI:16838"
CHEBI_RIBOSE_1_PHOSPHATE = "CHEBI:57720"
CHEBI_URIDINE = "CHEBI:17568"
CHEBI_DIHYDROCARVEOL = "CHEBI:50215"
CHEBI_DIHYDROCARVONE = "CHEBI:23733"
CHEBI_D_GLUCOSE = "CHEBI:4167"
CHEBI_D_GLUCONO_15_LACTONE = "CHEBI:16217"
CHEBI_D_TYROSYL_TRNA = "CHEBI:78723"
CHEBI_D_TYROSINE = "CHEBI:58570"
CHEBI_TRNA = "CHEBI:78442"


class ParticipantInput(TypedDict, total=False):
    chebi_id: str
    smiles: str
    inchi: str
    name: str
    count: int
    polymer_type: PolymerType


def _reaction(left: list[ParticipantInput], right: list[ParticipantInput]) -> Reaction:
    return Reaction(
        left_participants=[Participant(**participant) for participant in left],
        right_participants=[Participant(**participant) for participant in right],
    )


def test_go_0004550_positive_nucleoside_diphosphate_kinase() -> None:
    cls = NucleosideDiphosphateKinase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_UDP, "name": "UDP"},
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
        ],
        right=[
            {"chebi_id": CHEBI_UTP, "name": "UTP"},
            {"chebi_id": CHEBI_ADP, "name": "ADP"},
        ],
    )

    result = cls.check_membership_impl(reaction)
    assert result.is_member is True
    assert "nucleoside diphosphate kinase" in result.explanation.lower()


def test_go_0004550_rejects_wrong_nucleotide_branch() -> None:
    cls = NucleosideDiphosphateKinase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_UDP, "name": "UDP"},
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
        ],
        right=[
            {"chebi_id": CHEBI_AMP, "name": "AMP"},
            {"chebi_id": CHEBI_ADP, "name": "ADP"},
        ],
    )

    assert cls.check_membership_impl(reaction).is_member is False


def test_go_0004731_positive_purine_nucleoside_phosphorylase() -> None:
    cls = PurineNucleosidePhosphorylase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_GUANOSINE, "name": "guanosine"},
            {"chebi_id": CHEBI_PHOSPHATE, "name": "phosphate"},
        ],
        right=[
            {"chebi_id": CHEBI_GUANINE, "name": "guanine"},
            {"chebi_id": CHEBI_RIBOSE_1_PHOSPHATE, "name": "alpha-D-ribose 1-phosphate"},
        ],
    )

    result = cls.check_membership_impl(reaction)
    assert result.is_member is True
    assert "purine-nucleoside phosphorylase" in result.explanation.lower()


def test_go_0004731_rejects_pyrimidine_nucleoside() -> None:
    cls = PurineNucleosidePhosphorylase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_URIDINE, "name": "uridine"},
            {"chebi_id": CHEBI_PHOSPHATE, "name": "phosphate"},
        ],
        right=[
            {"chebi_id": CHEBI_RIBOSE_1_PHOSPHATE, "name": "alpha-D-ribose 1-phosphate"},
            {"chebi_id": CHEBI_GUANINE, "name": "guanine"},
        ],
    )

    assert cls.check_membership_impl(reaction).is_member is False


def test_ec_111296_positive_dihydrocarveol_dehydrogenase() -> None:
    cls = DihydrocarveolDehydrogenase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_DIHYDROCARVEOL, "name": "dihydrocarveol"},
            {"chebi_id": CHEBI_NAD_PLUS, "name": "NAD+"},
        ],
        right=[
            {"chebi_id": CHEBI_DIHYDROCARVONE, "name": "dihydrocarvone"},
            {"chebi_id": CHEBI_NADH, "name": "NADH"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)
    assert result.is_member is True
    assert "dihydrocarveol dehydrogenase" in result.explanation.lower()


def test_ec_111296_rejects_wrong_redox_cofactor() -> None:
    cls = DihydrocarveolDehydrogenase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_DIHYDROCARVEOL, "name": "dihydrocarveol"},
            {"chebi_id": CHEBI_NADP_PLUS, "name": "NADP+"},
        ],
        right=[
            {"chebi_id": CHEBI_DIHYDROCARVONE, "name": "dihydrocarvone"},
            {"chebi_id": CHEBI_NADPH, "name": "NADPH"},
        ],
    )

    assert cls.check_membership_impl(reaction).is_member is False


def test_ec_111359_positive_aldose_1_dehydrogenase_nad_p() -> None:
    cls = Aldose1DehydrogenaseNADP()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_D_GLUCOSE, "name": "D-glucose"},
            {"chebi_id": CHEBI_NADP_PLUS, "name": "NADP+"},
        ],
        right=[
            {"chebi_id": CHEBI_D_GLUCONO_15_LACTONE, "name": "D-glucono-1,5-lactone"},
            {"chebi_id": CHEBI_NADPH, "name": "NADPH"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)
    assert result.is_member is True
    assert "aldose 1-dehydrogenase" in result.explanation.lower()


def test_ec_111359_rejects_wrong_oxidized_product() -> None:
    cls = Aldose1DehydrogenaseNADP()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_D_GLUCOSE, "name": "D-glucose"},
            {"chebi_id": CHEBI_NADP_PLUS, "name": "NADP+"},
        ],
        right=[
            {"chebi_id": CHEBI_RIBOSE_1_PHOSPHATE, "name": "alpha-D-ribose 1-phosphate"},
            {"chebi_id": CHEBI_NADPH, "name": "NADPH"},
        ],
    )

    assert cls.check_membership_impl(reaction).is_member is False


def test_go_0051499_positive_d_aminoacyl_trna_deacylase() -> None:
    cls = DAminoacylTRNADeacylase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_D_TYROSYL_TRNA, "name": "D-tyrosyl-tRNA", "polymer_type": PolymerType.TRNA},
            {"chebi_id": CHEBI_H2O, "name": "water"},
        ],
        right=[
            {"chebi_id": CHEBI_D_TYROSINE, "name": "D-tyrosine"},
            {"chebi_id": CHEBI_TRNA, "name": "tRNA", "polymer_type": PolymerType.TRNA},
        ],
    )

    result = cls.check_membership_impl(reaction)
    assert result.is_member is True
    assert "d-aminoacyl-trna deacylase" in result.explanation.lower()


def test_go_0051499_accepts_reverse_orientation() -> None:
    cls = DAminoacylTRNADeacylase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_D_TYROSINE, "name": "D-tyrosine"},
            {"chebi_id": CHEBI_TRNA, "name": "tRNA", "polymer_type": PolymerType.TRNA},
        ],
        right=[
            {"chebi_id": CHEBI_D_TYROSYL_TRNA, "name": "D-tyrosyl-tRNA", "polymer_type": PolymerType.TRNA},
            {"chebi_id": CHEBI_H2O, "name": "water"},
        ],
    )

    assert cls.check_membership_impl(reaction).is_member is True


def test_go_0051499_rejects_missing_water() -> None:
    cls = DAminoacylTRNADeacylase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_D_TYROSYL_TRNA, "name": "D-tyrosyl-tRNA", "polymer_type": PolymerType.TRNA},
        ],
        right=[
            {"chebi_id": CHEBI_D_TYROSINE, "name": "D-tyrosine"},
            {"chebi_id": CHEBI_TRNA, "name": "tRNA", "polymer_type": PolymerType.TRNA},
        ],
    )

    assert cls.check_membership_impl(reaction).is_member is False
