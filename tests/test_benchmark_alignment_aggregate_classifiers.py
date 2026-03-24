"""Focused tests for benchmark-aligned aggregate classifiers."""

from __future__ import annotations

from typing_extensions import TypedDict

from autarch.datamodel import Participant, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_ATP,
    CHEBI_DIPHOSPHATE,
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_H_PLUS,
    CHEBI_NADP_PLUS,
    CHEBI_NADPH,
    CHEBI_O2,
)
from autarch.ontology.hydrolase_acting_on_acid_anhydrides import (
    HydrolaseActingOnAcidAnhydrides,
)
from autarch.ontology.hydrolase_acting_on_ester_bonds import (
    HydrolaseActingOnEsterBonds,
)
from autarch.ontology.carbon_carbon_lyase import CarbonCarbonLyase
from autarch.ontology.carbon_oxygen_lyase import CarbonOxygenLyase
from autarch.ontology.hydro_lyase import HydroLyase
from autarch.ontology.oxidoreductase_acting_on_peroxide_as_acceptor import (
    OxidoreductaseActingOnPeroxideAsAcceptor,
)
from autarch.ontology.oxidoreductase_acting_on_paired_donors_with_incorporation_or_reduction_of_molecular_oxygen import (
    OxidoreductaseActingOnPairedDonorsWithIncorporationOrReductionOfMolecularOxygen,
)
from autarch.ontology.oxidoreductase_acting_on_the_aldehyde_or_oxo_group_of_donors import (
    OxidoreductaseActingOnTheAldehydeOrOxoGroupOfDonors,
)
from autarch.ontology.transferase_transferring_one_carbon_groups import (
    TransferaseTransferringOneCarbonGroups,
)

CHEBI_R_LIMONENE = "CHEBI:15382"
CHEBI_LIMONENE_12_EPOXIDE = "CHEBI:16431"
CHEBI_ACETALDEHYDE = "CHEBI:15343"
CHEBI_ACETATE = "CHEBI:30089"
CHEBI_GLYCEROL = "CHEBI:17522"
CHEBI_3_HYDROXYPROPANAL = "CHEBI:17871"
CHEBI_PHENOL = "CHEBI:15882"
CHEBI_PHOSPHATE = "CHEBI:18367"
CHEBI_GIBBERELLIN_A4 = "CHEBI:73251"
CHEBI_GIBBERELLIN_A4_METHYL_ESTER = "CHEBI:73252"
CHEBI_SAM = "CHEBI:59789"
CHEBI_SAH = "CHEBI:57856"
CHEBI_CITRATE = "CHEBI:30769"
CHEBI_OXALOACETATE = "CHEBI:16452"
CHEBI_ACETATE_CO2 = "CHEBI:30089"
CHEBI_FARNESYL_DIPHOSPHATE = "CHEBI:175763"
CHEBI_GAMMA_HUMULENE = "CHEBI:49290"


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


def test_go_aligned_oxygenase_parent_wrapper_matches() -> None:
    cls = OxidoreductaseActingOnPairedDonorsWithIncorporationOrReductionOfMolecularOxygen()
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


def test_go_aligned_transferase_one_carbon_wrapper_matches() -> None:
    cls = TransferaseTransferringOneCarbonGroups()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_GIBBERELLIN_A4}, {"chebi_id": CHEBI_SAM}],
        right=[{"chebi_id": CHEBI_GIBBERELLIN_A4_METHYL_ESTER}, {"chebi_id": CHEBI_SAH}],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_go_aligned_hydrolase_ester_wrapper_matches() -> None:
    cls = HydrolaseActingOnEsterBonds()
    reaction = _reaction(
        left=[{"smiles": "CC(=O)OC"}, {"chebi_id": CHEBI_H2O, "smiles": "O"}],
        right=[{"chebi_id": CHEBI_ACETATE, "smiles": "CC(=O)[O-]"}, {"smiles": "CO"}],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_go_aligned_carbon_oxygen_lyase_wrapper_matches() -> None:
    cls = CarbonOxygenLyase()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_FARNESYL_DIPHOSPHATE}],
        right=[{"chebi_id": CHEBI_GAMMA_HUMULENE}, {"chebi_id": CHEBI_DIPHOSPHATE}],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_go_aligned_carbon_carbon_lyase_wrapper_matches() -> None:
    cls = CarbonCarbonLyase()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_CITRATE, "smiles": "O=C([O-])CC(O)(CC(=O)[O-])C(=O)[O-]"}],
        right=[
            {"chebi_id": CHEBI_OXALOACETATE, "smiles": "O=C([O-])CC(=O)C(=O)[O-]"},
            {"chebi_id": CHEBI_ACETATE_CO2, "smiles": "CC(=O)[O-]"},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_hydrolase_acting_on_acid_anhydrides_positive() -> None:
    cls = HydrolaseActingOnAcidAnhydrides()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_ATP}, {"chebi_id": CHEBI_H2O, "smiles": "O"}],
        right=[{"chebi_id": CHEBI_ADP}, {"chebi_id": CHEBI_PHOSPHATE}],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_hydro_lyase_positive() -> None:
    cls = HydroLyase()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_GLYCEROL, "smiles": "OCC(O)CO"}],
        right=[
            {"chebi_id": CHEBI_3_HYDROXYPROPANAL, "smiles": "OCCC=O"},
            {"chebi_id": CHEBI_H2O, "smiles": "O"},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_aldehyde_oxo_donor_parent_positive() -> None:
    cls = OxidoreductaseActingOnTheAldehydeOrOxoGroupOfDonors()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_ACETALDEHYDE, "smiles": "CC=O"},
            {"chebi_id": CHEBI_NADP_PLUS},
            {"chebi_id": CHEBI_H2O, "smiles": "O"},
        ],
        right=[
            {"chebi_id": CHEBI_ACETATE, "smiles": "CC(=O)[O-]"},
            {"chebi_id": CHEBI_NADPH},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "count": 2},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_peroxide_acceptor_parent_positive() -> None:
    cls = OxidoreductaseActingOnPeroxideAsAcceptor()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_PHENOL},
            {"chebi_id": CHEBI_H2O2},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]"},
        ],
        right=[
            {"smiles": "c1ccccc1[O]"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "count": 2},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True
