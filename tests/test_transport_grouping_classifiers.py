"""Focused tests for transport grouping classifiers."""

from __future__ import annotations

from typing_extensions import TypedDict

from autarch.datamodel import Participant, Reaction
from autarch.molecules import CHEBI_ADP, CHEBI_ATP, CHEBI_H2O, CHEBI_H_PLUS
from autarch.ontology.abc_type_carbohydrate_transporter import (
    ABCTypeCarbohydrateTransporter,
)
from autarch.ontology.translocation_of_amino_acids_and_peptides_linked_to_the_hydrolysis_of_a_nucleoside_triphosphate import (
    TranslocationOfAminoAcidsAndPeptidesLinkedToTheHydrolysisOfANucleosideTriphosphate,
)
from autarch.ontology.translocation_of_inorganic_anions_and_their_chelates_linked_to_the_hydrolysis_of_a_nucleoside_triphosphate import (
    TranslocationOfInorganicAnionsAndTheirChelatesLinkedToTheHydrolysisOfANucleosideTriphosphate,
)
from autarch.ontology.translocation_of_inorganic_cations_linked_to_the_hydrolysis_of_a_nucleoside_triphosphate import (
    TranslocationOfInorganicCationsLinkedToTheHydrolysisOfANucleosideTriphosphate,
)

CHEBI_POLYPHOSPHATE = "CHEBI:16838"
CHEBI_MAGNESIUM_2 = "CHEBI:18420"
CHEBI_SULFATE = "CHEBI:16189"
CHEBI_D_XYLOSE = "CHEBI:53455"
CHEBI_L_METHIONINE = "CHEBI:57844"
CHEBI_TAURINE = "CHEBI:507393"
CHEBI_GLYCINE_BETAINE = "CHEBI:17750"
CHEBI_GENERIC_COBALAMIN_CARRIER = "CHEBI:140785"


class ParticipantInput(TypedDict, total=False):
    """Subset of participant fields used in transport tests."""

    chebi_id: str
    smiles: str
    name: str
    location: str


def _reaction(left: list[ParticipantInput], right: list[ParticipantInput]) -> Reaction:
    return Reaction(
        left_participants=[Participant(**participant) for participant in left],
        right_participants=[Participant(**participant) for participant in right],
    )


def _atp_coupled_transport(
    substrate_left: ParticipantInput,
    substrate_right: ParticipantInput,
) -> Reaction:
    return _reaction(
        left=[
            substrate_left,
            {"chebi_id": CHEBI_ATP, "smiles": "P", "name": "ATP"},
            {"chebi_id": CHEBI_H2O, "smiles": "[H]O[H]", "name": "water"},
        ],
        right=[
            substrate_right,
            {"chebi_id": CHEBI_ADP, "smiles": "P", "name": "ADP"},
            {"chebi_id": CHEBI_POLYPHOSPHATE, "name": "polyphosphate"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )


def test_ec_732_positive_inorganic_anion_translocation() -> None:
    cls = TranslocationOfInorganicAnionsAndTheirChelatesLinkedToTheHydrolysisOfANucleosideTriphosphate()
    reaction = _atp_coupled_transport(
        {
            "chebi_id": CHEBI_SULFATE,
            "smiles": "[O-]S([O-])(=O)=O",
            "name": "sulfate",
            "location": "out",
        },
        {
            "chebi_id": CHEBI_SULFATE,
            "smiles": "[O-]S([O-])(=O)=O",
            "name": "sulfate",
            "location": "in",
        },
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "anions" in result.explanation.lower()


def test_ec_732_rejects_cation_transport() -> None:
    cls = TranslocationOfInorganicAnionsAndTheirChelatesLinkedToTheHydrolysisOfANucleosideTriphosphate()
    reaction = _atp_coupled_transport(
        {
            "chebi_id": CHEBI_MAGNESIUM_2,
            "smiles": "[Mg+2]",
            "name": "magnesium(2+)",
            "location": "out",
        },
        {
            "chebi_id": CHEBI_MAGNESIUM_2,
            "smiles": "[Mg+2]",
            "name": "magnesium(2+)",
            "location": "in",
        },
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_722_positive_inorganic_cation_translocation() -> None:
    cls = TranslocationOfInorganicCationsLinkedToTheHydrolysisOfANucleosideTriphosphate()
    reaction = _atp_coupled_transport(
        {
            "chebi_id": CHEBI_MAGNESIUM_2,
            "smiles": "[Mg+2]",
            "name": "magnesium(2+)",
            "location": "out",
        },
        {
            "chebi_id": CHEBI_MAGNESIUM_2,
            "smiles": "[Mg+2]",
            "name": "magnesium(2+)",
            "location": "in",
        },
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "cations" in result.explanation.lower()


def test_ec_722_rejects_taurine_transport() -> None:
    cls = TranslocationOfInorganicCationsLinkedToTheHydrolysisOfANucleosideTriphosphate()
    reaction = _atp_coupled_transport(
        {
            "chebi_id": CHEBI_TAURINE,
            "smiles": "[NH3+]CCS(=O)(=O)[O-]",
            "name": "taurine",
            "location": "out",
        },
        {
            "chebi_id": CHEBI_TAURINE,
            "smiles": "[NH3+]CCS(=O)(=O)[O-]",
            "name": "taurine",
            "location": "in",
        },
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_742_positive_amino_acid_peptide_translocation() -> None:
    cls = TranslocationOfAminoAcidsAndPeptidesLinkedToTheHydrolysisOfANucleosideTriphosphate()
    reaction = _atp_coupled_transport(
        {
            "chebi_id": CHEBI_L_METHIONINE,
            "smiles": "CSCC[C@H]([NH3+])C(=O)[O-]",
            "name": "L-methionine zwitterion",
            "location": "out",
        },
        {
            "chebi_id": CHEBI_L_METHIONINE,
            "smiles": "CSCC[C@H]([NH3+])C(=O)[O-]",
            "name": "L-methionine zwitterion",
            "location": "in",
        },
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "amino acids and peptides" in result.explanation.lower()


def test_ec_742_rejects_carbohydrate_transport() -> None:
    cls = TranslocationOfAminoAcidsAndPeptidesLinkedToTheHydrolysisOfANucleosideTriphosphate()
    reaction = _atp_coupled_transport(
        {
            "chebi_id": CHEBI_D_XYLOSE,
            "smiles": "OC[C@H]1O[C@H](O)[C@@H](O)[C@H](O)C1",
            "name": "D-xylopyranose",
            "location": "out",
        },
        {
            "chebi_id": CHEBI_D_XYLOSE,
            "smiles": "OC[C@H]1O[C@H](O)[C@@H](O)[C@H](O)C1",
            "name": "D-xylopyranose",
            "location": "in",
        },
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_742_rejects_amino_acid_derivative_transport() -> None:
    cls = TranslocationOfAminoAcidsAndPeptidesLinkedToTheHydrolysisOfANucleosideTriphosphate()
    reaction = _atp_coupled_transport(
        {
            "chebi_id": CHEBI_GLYCINE_BETAINE,
            "smiles": "C[N+](C)(C)CC(=O)[O-]",
            "name": "glycine betaine",
            "location": "out",
        },
        {
            "chebi_id": CHEBI_GLYCINE_BETAINE,
            "smiles": "C[N+](C)(C)CC(=O)[O-]",
            "name": "glycine betaine",
            "location": "in",
        },
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_go_0043211_positive_abc_type_carbohydrate_transport() -> None:
    cls = ABCTypeCarbohydrateTransporter()
    reaction = _atp_coupled_transport(
        {
            "chebi_id": CHEBI_D_XYLOSE,
            "smiles": "OC[C@H]1O[C@H](O)[C@@H](O)[C@H](O)C1",
            "name": "D-xylopyranose",
            "location": "out",
        },
        {
            "chebi_id": CHEBI_D_XYLOSE,
            "smiles": "OC[C@H]1O[C@H](O)[C@@H](O)[C@H](O)C1",
            "name": "D-xylopyranose",
            "location": "in",
        },
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "carbohydrate" in result.explanation.lower()


def test_go_0043211_rejects_amino_acid_transport() -> None:
    cls = ABCTypeCarbohydrateTransporter()
    reaction = _atp_coupled_transport(
        {
            "chebi_id": CHEBI_L_METHIONINE,
            "smiles": "CSCC[C@H]([NH3+])C(=O)[O-]",
            "name": "L-methionine zwitterion",
            "location": "out",
        },
        {
            "chebi_id": CHEBI_L_METHIONINE,
            "smiles": "CSCC[C@H]([NH3+])C(=O)[O-]",
            "name": "L-methionine zwitterion",
            "location": "in",
        },
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_go_0043211_ignores_generic_transport_carrier_placeholder() -> None:
    cls = ABCTypeCarbohydrateTransporter()
    reaction = _reaction(
        left=[
            {
                "chebi_id": CHEBI_GENERIC_COBALAMIN_CARRIER,
                "name": "generic carrier scaffold",
                "location": "out",
            },
            {"chebi_id": CHEBI_ATP, "smiles": "P", "name": "ATP"},
            {"chebi_id": CHEBI_H2O, "smiles": "[H]O[H]", "name": "water"},
            {"name": "an R-cob(III)alamin", "location": "out"},
        ],
        right=[
            {
                "chebi_id": CHEBI_GENERIC_COBALAMIN_CARRIER,
                "name": "generic carrier scaffold",
                "location": "in",
            },
            {"chebi_id": CHEBI_ADP, "smiles": "P", "name": "ADP"},
            {"chebi_id": CHEBI_POLYPHOSPHATE, "name": "polyphosphate"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
            {"name": "an R-cob(III)alamin", "location": "in"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False
