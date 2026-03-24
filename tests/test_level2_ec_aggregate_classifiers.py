"""Focused tests for level-2 EC aggregate classifiers."""

from __future__ import annotations

from typing_extensions import TypedDict

from autarch.datamodel import Participant, Reaction
from autarch.molecules import (
    CHEBI_ATP,
    CHEBI_DIPHOSPHATE,
    CHEBI_FAD,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_O2,
)
from autarch.ontology.carbon_carbon_lyases import CarbonCarbonLyases
from autarch.ontology.carbon_oxygen_lyases import CarbonOxygenLyases
from autarch.ontology.acting_on_ester_bonds import ActingOnEsterBonds
from autarch.ontology.acting_on_paired_donors_with_incorporation_or_reduction_of_molecular_oxygen import (
    ActingOnPairedDonorsWithIncorporationOrReductionOfMolecularOxygen,
)
from autarch.ontology.transferring_alkyl_or_aryl_groups_other_than_methyl_groups import (
    TransferringAlkylOrArylGroupsOtherThanMethylGroups,
)
from autarch.ontology.transferring_nitrogenous_groups import (
    TransferringNitrogenousGroups,
)
from autarch.ontology.transferring_one_carbon_groups import (
    TransferringOneCarbonGroups,
)
from autarch.ontology.transferring_sulfur_containing_groups import (
    TransferringSulfurContainingGroups,
)

CHEBI_SAM = "CHEBI:59789"
CHEBI_SAH = "CHEBI:57856"
CHEBI_GIBBERELLIN_A4 = "CHEBI:73251"
CHEBI_GIBBERELLIN_A4_METHYL_ESTER = "CHEBI:73252"
CHEBI_REDUCED_ETF = "CHEBI:58307"
CHEBI_TRIPHOSPHATE = "CHEBI:18036"
CHEBI_COB_II_ALAMIN = "CHEBI:16304"
CHEBI_ADENOSYLCOBALAMIN = "CHEBI:18408"
CHEBI_TYROSINE = "CHEBI:58315"
CHEBI_2_OXOGLUTARATE = "CHEBI:16810"
CHEBI_TYROSINE_KETO_ACID = "CHEBI:36242"
CHEBI_GLUTAMATE = "CHEBI:29985"
CHEBI_THIOSULFATE = "CHEBI:32594"
CHEBI_HYDROGEN_CYANIDE = "CHEBI:18432"
CHEBI_THIOCYANATE = "CHEBI:30690"
CHEBI_SULFITE = "CHEBI:18498"
CHEBI_CITRATE = "CHEBI:30769"
CHEBI_OXALOACETATE = "CHEBI:16452"
CHEBI_ACETATE = "CHEBI:30089"
CHEBI_FARNESYL_DIPHOSPHATE = "CHEBI:175763"
CHEBI_GAMMA_HUMULENE = "CHEBI:49290"
CHEBI_R_LIMONENE = "CHEBI:15382"
CHEBI_LIMONENE_12_EPOXIDE = "CHEBI:16431"


class ParticipantInput(TypedDict, total=False):
    chebi_id: str
    smiles: str
    count: int
    location: str


def _reaction(left: list[ParticipantInput], right: list[ParticipantInput]) -> Reaction:
    return Reaction(
        left_participants=[Participant(**participant) for participant in left],
        right_participants=[Participant(**participant) for participant in right],
    )


def test_ec_21_positive_transferase_transferring_one_carbon_groups() -> None:
    cls = TransferringOneCarbonGroups()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_GIBBERELLIN_A4}, {"chebi_id": CHEBI_SAM}],
        right=[{"chebi_id": CHEBI_GIBBERELLIN_A4_METHYL_ESTER}, {"chebi_id": CHEBI_SAH}],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_25_positive_transferase_transferring_alkyl_or_aryl_groups() -> None:
    cls = TransferringAlkylOrArylGroupsOtherThanMethylGroups()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_COB_II_ALAMIN},
            {"chebi_id": CHEBI_COB_II_ALAMIN},
            {"chebi_id": CHEBI_REDUCED_ETF},
            {"chebi_id": CHEBI_ATP},
            {"chebi_id": CHEBI_ATP},
        ],
        right=[
            {"chebi_id": CHEBI_ADENOSYLCOBALAMIN},
            {"chebi_id": CHEBI_ADENOSYLCOBALAMIN},
            {"chebi_id": CHEBI_TRIPHOSPHATE},
            {"chebi_id": CHEBI_TRIPHOSPHATE},
            {"chebi_id": CHEBI_FAD},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_26_positive_transferase_transferring_nitrogenous_groups() -> None:
    cls = TransferringNitrogenousGroups()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_TYROSINE},
            {"chebi_id": CHEBI_2_OXOGLUTARATE},
        ],
        right=[
            {"chebi_id": CHEBI_TYROSINE_KETO_ACID},
            {"chebi_id": CHEBI_GLUTAMATE},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_28_positive_transferase_transferring_sulfur_containing_groups() -> None:
    cls = TransferringSulfurContainingGroups()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_THIOSULFATE, "smiles": "[H]SS(=O)(=O)[O-]"},
            {"chebi_id": CHEBI_HYDROGEN_CYANIDE, "smiles": "C#N"},
        ],
        right=[
            {"chebi_id": CHEBI_THIOCYANATE, "smiles": "N#C[S-]"},
            {"chebi_id": CHEBI_SULFITE, "smiles": "O=S([O-])[O-]"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "count": 2},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_31_positive_hydrolase_acting_on_ester_bonds() -> None:
    cls = ActingOnEsterBonds()
    reaction = _reaction(
        left=[{"smiles": "CC(=O)OC"}, {"chebi_id": CHEBI_H2O, "smiles": "O"}],
        right=[{"chebi_id": CHEBI_ACETATE, "smiles": "CC(=O)[O-]"}, {"smiles": "CO"}],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_41_positive_carbon_carbon_lyase() -> None:
    cls = CarbonCarbonLyases()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_CITRATE, "smiles": "O=C([O-])CC(O)(CC(=O)[O-])C(=O)[O-]"}],
        right=[
            {"chebi_id": CHEBI_OXALOACETATE, "smiles": "O=C([O-])CC(=O)C(=O)[O-]"},
            {"chebi_id": CHEBI_ACETATE, "smiles": "CC(=O)[O-]"},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_42_positive_carbon_oxygen_lyase() -> None:
    cls = CarbonOxygenLyases()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_FARNESYL_DIPHOSPHATE}],
        right=[{"chebi_id": CHEBI_GAMMA_HUMULENE}, {"chebi_id": CHEBI_DIPHOSPHATE}],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_114_positive_level2_oxygenase_parent() -> None:
    cls = ActingOnPairedDonorsWithIncorporationOrReductionOfMolecularOxygen()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_R_LIMONENE},
            {"chebi_id": CHEBI_NADPH},
            {"chebi_id": CHEBI_O2},
            {"chebi_id": CHEBI_H_PLUS},
        ],
        right=[
            {"chebi_id": CHEBI_LIMONENE_12_EPOXIDE},
            {"chebi_id": CHEBI_NADP_PLUS},
            {"chebi_id": CHEBI_H2O},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True
