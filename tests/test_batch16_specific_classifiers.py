"""Focused tests for batch 16 specific classifiers."""

from __future__ import annotations

from typing_extensions import TypedDict

from autarch.datamodel import Participant, Reaction
from autarch.molecules import CHEBI_CO2, CHEBI_H2O, CHEBI_H_PLUS, CHEBI_SAH, CHEBI_SAM
from autarch.ontology.rs1_benzyl_1234_tetrahydroisoquinoline_n_methyltransferase import (
    RS1Benzyl1234TetrahydroisoquinolineNMethyltransferase,
)
from autarch.ontology.glycerol_dehydratase import GlycerolDehydratase
from autarch.ontology.five_methyltetrahydropteroyltriglutamate_homocysteine_s_methyltransferase import (
    FiveMethyltetrahydropteroyltriglutamateHomocysteineSMethyltransferase,
)
from autarch.ontology.prephenate_dehydratase import PrephenateDehydratase
from autarch.ontology.trimethylsulfonium_tetrahydrofolate_n_methyltransferase import (
    TrimethylsulfoniumTetrahydrofolateNMethyltransferase,
)

CHEBI_BENZYL_TETRAHYDROISOQUINOLINE = "CHEBI:57902"
CHEBI_N_METHYL_BENZYL_TETRAHYDROISOQUINOLINE = "CHEBI:57598"
CHEBI_METHYLTHF_TRIGLU = "CHEBI:58207"
CHEBI_THF_TRIGLU = "CHEBI:58140"
CHEBI_HOMOCYSTEINE = "CHEBI:58199"
CHEBI_METHIONINE = "CHEBI:57844"
CHEBI_TRIMETHYLSULFONIUM = "CHEBI:17434"
CHEBI_DIMETHYL_SULFIDE = "CHEBI:17437"
CHEBI_THF = "CHEBI:57453"
CHEBI_METHYL_THF = "CHEBI:18608"
CHEBI_GLYCEROL = "CHEBI:17522"
CHEBI_HYDROXYPROPANAL = "CHEBI:17871"
CHEBI_GLYCERALDEHYDE = "CHEBI:17167"
CHEBI_PREPHENATE = "CHEBI:29934"
CHEBI_PHENYLPYRUVATE = "CHEBI:18005"
CHEBI_HYDROXYPHENYLPYRUVATE = "CHEBI:35290"


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


def test_go_0030776_positive_benzyl_tetrahydroisoquinoline_n_methyltransferase() -> None:
    cls = RS1Benzyl1234TetrahydroisoquinolineNMethyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_BENZYL_TETRAHYDROISOQUINOLINE, "name": "(RS)-1-benzyl-1,2,3,4-tetrahydroisoquinoline"},
            {"chebi_id": CHEBI_SAM, "name": "SAM"},
        ],
        right=[
            {"chebi_id": CHEBI_N_METHYL_BENZYL_TETRAHYDROISOQUINOLINE, "name": "N-methyl benzyl tetrahydroisoquinoline"},
            {"chebi_id": CHEBI_SAH, "name": "SAH"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "isoquinoline" in result.explanation.lower()


def test_go_0030776_rejects_unmethylated_isoquinoline_product() -> None:
    cls = RS1Benzyl1234TetrahydroisoquinolineNMethyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_BENZYL_TETRAHYDROISOQUINOLINE, "name": "(RS)-1-benzyl-1,2,3,4-tetrahydroisoquinoline"},
            {"chebi_id": CHEBI_SAM, "name": "SAM"},
        ],
        right=[
            {"chebi_id": CHEBI_BENZYL_TETRAHYDROISOQUINOLINE, "name": "(RS)-1-benzyl-1,2,3,4-tetrahydroisoquinoline"},
            {"chebi_id": CHEBI_SAH, "name": "SAH"},
        ],
    )

    assert cls.check_membership_impl(reaction).is_member is False


def test_go_0003871_positive_methyltetrahydropteroyltriglutamate_homocysteine_s_methyltransferase() -> None:
    cls = FiveMethyltetrahydropteroyltriglutamateHomocysteineSMethyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_METHYLTHF_TRIGLU, "name": "5-methyltetrahydropteroyltri-L-glutamate"},
            {"chebi_id": CHEBI_HOMOCYSTEINE, "name": "L-homocysteine"},
        ],
        right=[
            {"chebi_id": CHEBI_THF_TRIGLU, "name": "tetrahydropteroyltri-L-glutamate"},
            {"chebi_id": CHEBI_METHIONINE, "name": "L-methionine"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "homocysteine" in result.explanation.lower() or "methyltransferase" in result.explanation.lower()


def test_go_0003871_rejects_missing_methionine_product() -> None:
    cls = FiveMethyltetrahydropteroyltriglutamateHomocysteineSMethyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_METHYLTHF_TRIGLU, "name": "5-methyltetrahydropteroyltri-L-glutamate"},
            {"chebi_id": CHEBI_HOMOCYSTEINE, "name": "L-homocysteine"},
        ],
        right=[
            {"chebi_id": CHEBI_THF_TRIGLU, "name": "tetrahydropteroyltri-L-glutamate"},
            {"chebi_id": CHEBI_HOMOCYSTEINE, "name": "L-homocysteine"},
        ],
    )

    assert cls.check_membership_impl(reaction).is_member is False


def test_go_0047147_positive_trimethylsulfonium_tetrahydrofolate_n_methyltransferase() -> None:
    cls = TrimethylsulfoniumTetrahydrofolateNMethyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_TRIMETHYLSULFONIUM, "name": "trimethylsulfonium"},
            {"chebi_id": CHEBI_THF, "name": "tetrahydrofolate"},
        ],
        right=[
            {"chebi_id": CHEBI_DIMETHYL_SULFIDE, "name": "dimethyl sulfide"},
            {"chebi_id": CHEBI_METHYL_THF, "name": "5-methyltetrahydrofolate"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "trimethylsulfonium" in result.explanation.lower() or "tetrahydrofolate" in result.explanation.lower()


def test_go_0047147_rejects_dimethylsulfide_only_conversion() -> None:
    cls = TrimethylsulfoniumTetrahydrofolateNMethyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_TRIMETHYLSULFONIUM, "name": "trimethylsulfonium"},
            {"chebi_id": CHEBI_THF, "name": "tetrahydrofolate"},
        ],
        right=[
            {"chebi_id": CHEBI_DIMETHYL_SULFIDE, "name": "dimethyl sulfide"},
            {"chebi_id": CHEBI_THF, "name": "tetrahydrofolate"},
        ],
    )

    assert cls.check_membership_impl(reaction).is_member is False


def test_go_0046405_positive_glycerol_dehydratase() -> None:
    cls = GlycerolDehydratase()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_GLYCEROL, "name": "glycerol"}],
        right=[
            {"chebi_id": CHEBI_HYDROXYPROPANAL, "name": "3-hydroxypropanal"},
            {"chebi_id": CHEBI_H2O, "name": "water"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "glycerol dehydratase" in result.explanation.lower()


def test_go_0046405_rejects_glyceraldehyde_formation() -> None:
    cls = GlycerolDehydratase()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_GLYCEROL, "name": "glycerol"}],
        right=[
            {"chebi_id": CHEBI_GLYCERALDEHYDE, "name": "glyceraldehyde"},
            {"chebi_id": CHEBI_H2O, "name": "water"},
        ],
    )

    assert cls.check_membership_impl(reaction).is_member is False


def test_go_0004664_positive_prephenate_dehydratase() -> None:
    cls = PrephenateDehydratase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_PREPHENATE, "name": "prephenate"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
        right=[
            {"chebi_id": CHEBI_PHENYLPYRUVATE, "name": "3-phenylpyruvate"},
            {"chebi_id": CHEBI_CO2, "name": "carbon dioxide"},
            {"chebi_id": CHEBI_H2O, "name": "water"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "prephenate dehydratase" in result.explanation.lower()


def test_go_0004664_rejects_hydroxylated_ketoacid_product() -> None:
    cls = PrephenateDehydratase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_PREPHENATE, "name": "prephenate"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
        right=[
            {"chebi_id": CHEBI_HYDROXYPHENYLPYRUVATE, "name": "hydroxyphenylpyruvate"},
            {"chebi_id": CHEBI_CO2, "name": "carbon dioxide"},
            {"chebi_id": CHEBI_H2O, "name": "water"},
        ],
    )

    assert cls.check_membership_impl(reaction).is_member is False
