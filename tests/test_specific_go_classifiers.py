"""Tests for GO-specific classifiers with strict concept alignment."""

from autarch.datamodel import Participant, Reaction
from autarch.ontology.aldose_1_epimerase import Aldose1Epimerase
from autarch.ontology.pectate_lyase import PectateLyase
from autarch.ontology.pectinesterase import Pectinesterase
from autarch.ontology.ribulose_bisphosphate_carboxylase import (
    RibuloseBisphosphateCarboxylase,
)
from autarch.ontology.squalene_monooxygenase import SqualeneMonooxygenase


def _reaction(left: list[dict], right: list[dict]) -> Reaction:
    return Reaction(
        left_participants=[Participant(**p) for p in left],
        right_participants=[Participant(**p) for p in right],
    )


def test_squalene_monooxygenase_positive() -> None:
    cls = SqualeneMonooxygenase()
    reaction = _reaction(
        left=[
            {"chebi_id": "CHEBI:15440"},  # squalene
            {"chebi_id": "CHEBI:57618"},  # reduced [NADPH--hemoprotein reductase]
            {"chebi_id": "CHEBI:15379"},  # O2
        ],
        right=[
            {"chebi_id": "CHEBI:15441"},  # (S)-2,3-epoxysqualene
            {"chebi_id": "CHEBI:58210"},  # oxidized [NADPH--hemoprotein reductase]
            {"chebi_id": "CHEBI:15377"},  # H2O
            {"chebi_id": "CHEBI:15378"},  # H+
        ],
    )

    result = cls.check_membership_impl(reaction)
    assert result.is_member is True
    assert "squalene monooxygenase" in result.explanation.lower()


def test_squalene_monooxygenase_rejects_generic_monooxygenase() -> None:
    cls = SqualeneMonooxygenase()
    reaction = _reaction(
        left=[
            {"chebi_id": "CHEBI:15379"},  # O2
            {"chebi_id": "CHEBI:57783"},  # NADPH
            {"chebi_id": "CHEBI:17230"},  # generic organic substrate
        ],
        right=[
            {"chebi_id": "CHEBI:57857"},  # NADP+
            {"chebi_id": "CHEBI:15377"},  # H2O
            {"chebi_id": "CHEBI:28885"},  # generic oxygenated product
        ],
    )

    result = cls.check_membership_impl(reaction)
    assert result.is_member is False


def test_rubisco_matches_rhea_23124_direction() -> None:
    cls = RibuloseBisphosphateCarboxylase()
    reaction = _reaction(
        left=[
            {"chebi_id": "CHEBI:58272"},  # 3-phosphoglycerate
            {"chebi_id": "CHEBI:58272"},
            {"chebi_id": "CHEBI:15378"},  # H+
            {"chebi_id": "CHEBI:15378"},
        ],
        right=[
            {"chebi_id": "CHEBI:57870"},  # ribulose 1,5-bisphosphate
            {"chebi_id": "CHEBI:16526"},  # CO2
            {"chebi_id": "CHEBI:15377"},  # H2O
        ],
    )

    result = cls.check_membership_impl(reaction)
    assert result.is_member is True
    assert "rubisco" in result.explanation.lower()


def test_rubisco_rejects_non_rubisco_carboxylation() -> None:
    cls = RibuloseBisphosphateCarboxylase()
    reaction = _reaction(
        left=[
            {"chebi_id": "CHEBI:32816"},  # pyruvate
            {"chebi_id": "CHEBI:16526"},  # CO2
            {"chebi_id": "CHEBI:15422"},  # ATP
        ],
        right=[
            {"chebi_id": "CHEBI:16452"},  # oxaloacetate
            {"chebi_id": "CHEBI:16761"},  # ADP
            {"chebi_id": "CHEBI:18367"},  # phosphate
        ],
    )

    result = cls.check_membership_impl(reaction)
    assert result.is_member is False


def test_aldose_1_epimerase_positive_glucose() -> None:
    cls = Aldose1Epimerase()
    reaction = _reaction(
        left=[{"chebi_id": "CHEBI:15444"}],  # alpha-D-glucose
        right=[{"chebi_id": "CHEBI:15903"}],  # beta-D-glucose
    )

    result = cls.check_membership_impl(reaction)
    assert result.is_member is True
    assert "aldose 1-epimerase" in result.explanation.lower()


def test_aldose_1_epimerase_rejects_label_only_fallback() -> None:
    cls = Aldose1Epimerase()
    reaction = _reaction(
        left=[{"name": "alpha-D-galactose"}],
        right=[{"name": "beta-D-galactose"}],
    )
    reaction.label = "alpha-D-galactose = beta-D-galactose"

    result = cls.check_membership_impl(reaction)
    assert result.is_member is False


def test_aldose_1_epimerase_rejects_non_anomer_pair() -> None:
    cls = Aldose1Epimerase()
    reaction = _reaction(
        left=[{"chebi_id": "CHEBI:17925"}],  # D-glucose (generic)
        right=[{"chebi_id": "CHEBI:17634"}],  # D-fructose
    )

    result = cls.check_membership_impl(reaction)
    assert result.is_member is False


def test_pectinesterase_positive() -> None:
    cls = Pectinesterase()
    reaction = _reaction(
        left=[
            {
                "chebi_id": "CHEBI:140522",
                "name": "[(1->4)-alpha-D-galacturonosyl methyl ester]",
            },
            {"chebi_id": "CHEBI:15377"},  # water
        ],
        right=[
            {"chebi_id": "CHEBI:140523", "name": "[(1->4)-alpha-D-galacturonosyl]"},
            {"chebi_id": "CHEBI:17790"},  # methanol
            {"chebi_id": "CHEBI:15378"},  # H+
        ],
    )

    result = cls.check_membership_impl(reaction)
    assert result.is_member is True
    assert "pectinesterase" in result.explanation.lower()


def test_pectinesterase_rejects_missing_methanol() -> None:
    cls = Pectinesterase()
    reaction = _reaction(
        left=[
            {
                "chebi_id": "CHEBI:140522",
                "name": "[(1->4)-alpha-D-galacturonosyl methyl ester]",
            },
            {"chebi_id": "CHEBI:15377"},
        ],
        right=[
            {"chebi_id": "CHEBI:140523", "name": "[(1->4)-alpha-D-galacturonosyl]"},
            {"chebi_id": "CHEBI:15378"},
        ],
    )

    result = cls.check_membership_impl(reaction)
    assert result.is_member is False


def test_pectate_lyase_positive() -> None:
    cls = PectateLyase()
    reaction = _reaction(
        left=[
            {"name": "pectate polymer"},
        ],
        right=[
            {"name": "4-deoxy-alpha-D-gluc-4-enuronosyl oligomer"},
            {"name": "shortened pectate"},
        ],
    )
    reaction.label = "pectate polymer = 4-deoxy-alpha-D-gluc-4-enuronosyl oligomer + shortened pectate"

    result = cls.check_membership_impl(reaction)
    assert result.is_member is True


def test_pectate_lyase_rejects_hydrolysis() -> None:
    cls = PectateLyase()
    reaction = _reaction(
        left=[
            {"name": "pectate polymer"},
            {"chebi_id": "CHEBI:15377"},
        ],
        right=[
            {"name": "4-deoxy-alpha-D-gluc-4-enuronosyl oligomer"},
        ],
    )
    reaction.label = "pectate polymer + H2O = 4-deoxy-alpha-D-gluc-4-enuronosyl oligomer"

    result = cls.check_membership_impl(reaction)
    assert result.is_member is False
