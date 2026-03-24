"""Focused tests for the next GO support>=10 batch."""

from __future__ import annotations

from typing_extensions import TypedDict

from autarch.datamodel import Participant, PolymerType, Reaction
from autarch.evaluation import evaluate_reaction_class
from autarch.ontology.dicarboxylic_acid_transmembrane_transporter_activity import (
    DicarboxylicAcidTransmembraneTransporterActivity,
)
from autarch.ontology.fatty_acid_ligase_activity import FattyAcidLigaseActivity
from autarch.ontology.neutral_l_amino_acid_transmembrane_transporter_activity import (
    NeutralLAminoAcidTransmembraneTransporterActivity,
)
from autarch.ontology.organic_acid_sodium_symporter_activity import (
    OrganicAcidSodiumSymporterActivity,
)
from autarch.ontology.p_type_ion_transporter_activity import PTypeIonTransporterActivity
from autarch.ontology.p_type_transmembrane_transporter_activity import (
    PTypeTransmembraneTransporterActivity,
)
from autarch.ontology.protein_kinase_activity import ProteinKinaseActivity
from autarch.ontology.protein_serine_threonine_kinase_activity import (
    ProteinSerineThreonineKinaseActivity,
)
from autarch.ontology.sesquiterpene_synthase_activity import SesquiterpeneSynthaseActivity
from autarch.ontology.transition_metal_ion_transmembrane_transporter_activity import (
    TransitionMetalIonTransmembraneTransporterActivity,
)

CHEBI_ADP = "CHEBI:456216"
CHEBI_ATP = "CHEBI:30616"
CHEBI_COA = "CHEBI:57287"
CHEBI_DIPHOSPHATE = "CHEBI:33019"
CHEBI_FARNESYL_DIPHOSPHATE = "CHEBI:175763"
CHEBI_FERROUS = "CHEBI:29033"
CHEBI_H2O = "CHEBI:15377"
CHEBI_METHIONINE = "CHEBI:57844"
CHEBI_SODIUM = "CHEBI:29101"
CHEBI_SUCCINATE = "CHEBI:15741"


class ParticipantInput(TypedDict, total=False):
    chebi_id: str
    smiles: str
    name: str
    location: str
    polymer_type: PolymerType


def _reaction(left: list[ParticipantInput], right: list[ParticipantInput]) -> Reaction:
    return Reaction(
        left_participants=[Participant(**participant) for participant in left],
        right_participants=[Participant(**participant) for participant in right],
    )


def test_transition_metal_ion_transmembrane_transporter_activity_positive() -> None:
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_FERROUS, "name": "Fe2+", "location": "out"}],
        right=[{"chebi_id": CHEBI_FERROUS, "name": "Fe2+", "location": "in"}],
    )
    assert TransitionMetalIonTransmembraneTransporterActivity().check_membership_impl(reaction).is_member is True


def test_p_type_wrappers_positive() -> None:
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_SODIUM, "name": "sodium ion", "location": "out"},
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
            {"chebi_id": CHEBI_H2O, "name": "water"},
        ],
        right=[
            {"chebi_id": CHEBI_SODIUM, "name": "sodium ion", "location": "in"},
            {"chebi_id": CHEBI_ADP, "name": "ADP"},
        ],
    )
    assert PTypeTransmembraneTransporterActivity().check_membership_impl(reaction).is_member is True
    assert PTypeIonTransporterActivity().check_membership_impl(reaction).is_member is True


def test_neutral_l_amino_acid_transmembrane_transporter_activity_positive() -> None:
    reaction = _reaction(
        left=[
            {
                "chebi_id": CHEBI_METHIONINE,
                "smiles": "CSCC[C@H]([NH3+])C(=O)[O-]",
                "name": "L-methionine",
                "location": "out",
            }
        ],
        right=[
            {
                "chebi_id": CHEBI_METHIONINE,
                "smiles": "CSCC[C@H]([NH3+])C(=O)[O-]",
                "name": "L-methionine",
                "location": "in",
            }
        ],
    )
    assert (
        NeutralLAminoAcidTransmembraneTransporterActivity().check_membership_impl(reaction).is_member
        is True
    )


def test_organic_acid_sodium_symporter_activity_positive() -> None:
    reaction = _reaction(
        left=[
            {
                "chebi_id": CHEBI_SUCCINATE,
                "smiles": "O=C([O-])CCC(=O)[O-]",
                "name": "succinate",
                "location": "out",
            },
            {"chebi_id": CHEBI_SODIUM, "name": "sodium ion", "location": "out"},
        ],
        right=[
            {
                "chebi_id": CHEBI_SUCCINATE,
                "smiles": "O=C([O-])CCC(=O)[O-]",
                "name": "succinate",
                "location": "in",
            },
            {"chebi_id": CHEBI_SODIUM, "name": "sodium ion", "location": "in"},
        ],
    )
    assert OrganicAcidSodiumSymporterActivity().check_membership_impl(reaction).is_member is True


def test_dicarboxylic_acid_transmembrane_transporter_activity_positive() -> None:
    reaction = _reaction(
        left=[
            {
                "chebi_id": CHEBI_SUCCINATE,
                "smiles": "O=C([O-])CCC(=O)[O-]",
                "name": "succinate",
                "location": "out",
            }
        ],
        right=[
            {
                "chebi_id": CHEBI_SUCCINATE,
                "smiles": "O=C([O-])CCC(=O)[O-]",
                "name": "succinate",
                "location": "in",
            }
        ],
    )
    assert DicarboxylicAcidTransmembraneTransporterActivity().check_membership_impl(reaction).is_member is True


def test_sesquiterpene_synthase_activity_positive() -> None:
    reaction = _reaction(
        left=[
            {
                "chebi_id": CHEBI_FARNESYL_DIPHOSPHATE,
                "name": "farnesyl diphosphate",
            }
        ],
        right=[
            {"chebi_id": "CHEBI:49290", "name": "gamma-humulene"},
            {"chebi_id": CHEBI_DIPHOSPHATE, "name": "diphosphate"},
        ],
    )
    assert SesquiterpeneSynthaseActivity().check_membership_impl(reaction).is_member is True


def test_fatty_acid_ligase_activity_positive() -> None:
    reaction = _reaction(
        left=[
            {"smiles": "CCCCCCCC(=O)[O-]", "name": "octanoate"},
            {"chebi_id": CHEBI_COA, "name": "CoA"},
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
        ],
        right=[
            {
                "smiles": "CCCCCCCC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(C(O)C1O)n1cnc2c(N)ncnc12",
                "name": "octanoyl-CoA",
            },
            {"chebi_id": "CHEBI:456215", "name": "AMP"},
            {"chebi_id": CHEBI_DIPHOSPHATE, "name": "diphosphate"},
        ],
    )
    assert FattyAcidLigaseActivity().check_membership_impl(reaction).is_member is True


def test_protein_kinase_wrappers_have_positive_support() -> None:
    protein_metrics = evaluate_reaction_class(
        ProteinKinaseActivity.__name__,
        cache_dir="cache",
        go_only=True,
    )
    ser_thr_metrics = evaluate_reaction_class(
        ProteinSerineThreonineKinaseActivity.__name__,
        cache_dir="cache",
        go_only=True,
    )
    assert protein_metrics.true_positives > 0
    assert ser_thr_metrics.true_positives > 0
