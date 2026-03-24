"""Focused tests for completing the remaining explicit EC level-2 classes."""

from __future__ import annotations

import pytest

from autarch.datamodel import Participant, PolymerType, Reaction
from autarch.molecules import (
    CHEBI_AMP,
    CHEBI_ATP,
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_H_PLUS,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NH3,
    CHEBI_O2,
)
from autarch.ontology.acting_on_carbon_carbon_bonds import ActingOnCarbonCarbonBonds
from autarch.ontology.acting_on_ether_bonds import ActingOnEtherBonds
from autarch.ontology.acting_on_halogen_in_donors import ActingOnHalogenInDonors
from autarch.ontology.acting_on_nadh_or_nadph import ActingOnNADHOrNADPH
from autarch.ontology.acting_on_superoxide_as_acceptor import ActingOnSuperoxideAsAcceptor
from autarch.ontology.atp_independent_chelatases import ATPIndependentChelatases
from autarch.ontology.carbon_nitrogen_lyases import CarbonNitrogenLyases
from autarch.ontology.cis_trans_isomerase import CisTransIsomerase
from autarch.ontology.enzymes_using_h2_as_reductant import EnzymesUsingH2AsReductant
from autarch.ontology.forming_carbon_oxygen_bonds import FormingCarbonOxygenBonds
from autarch.ontology.forming_nitrogen_nitrogen_bonds import FormingNitrogenNitrogenBonds
from autarch.ontology.nitrogen_oxygen_lyases import NitrogenOxygenLyases
from autarch.ontology.other_enzymes_using_o2_as_oxidant import OtherEnzymesUsingO2AsOxidant
from autarch.ontology.translocation_of_hydrons import TranslocationOfHydrons

CHEBI_SUPEROXIDE = "CHEBI:18421"
CHEBI_EPOXIDE = "CHEBI:35762"
CHEBI_DIOL = "CHEBI:23824"
CHEBI_FUMARYLACETOACETATE = "CHEBI:18034"
CHEBI_PEPTIDYL_PROLINE_CIS = "CHEBI:83155"
CHEBI_PEPTIDYL_PROLINE_TRANS = "CHEBI:83154"
CHEBI_DIPHOSPHATE = "CHEBI:33019"
CHEBI_TYROSINE = "CHEBI:58315"
CHEBI_TRNA = "CHEBI:64908"
CHEBI_CYTOCHROME_C_REDUCED = "CHEBI:29034"
CHEBI_CYTOCHROME_C_OXIDIZED = "CHEBI:29033"


def _reaction(
    left: list[Participant],
    right: list[Participant],
    label: str | None = None,
) -> Reaction:
    return Reaction(left_participants=left, right_participants=right, label=label or "")


def test_ec_16_positive_acting_on_nadh_or_nadph() -> None:
    reaction = _reaction(
        left=[
            Participant(smiles="O=C1C=CC(=O)C([O-])=C1", name="2-hydroxy-1,4-benzoquinone"),
            Participant(chebi_id=CHEBI_NADPH),
            Participant(chebi_id=CHEBI_H_PLUS),
        ],
        right=[
            Participant(smiles="OC1=CC(O)=C(O)C=C1", name="benzene-1,2,4-triol"),
            Participant(chebi_id=CHEBI_NADP_PLUS),
        ],
        label="2-hydroxy-1,4-benzoquinone + NADPH + H+ = benzene-1,2,4-triol + NADP+",
    )
    assert ActingOnNADHOrNADPH().check_membership(reaction).is_member is True


def test_ec_115_positive_acting_on_superoxide_as_acceptor() -> None:
    reaction = _reaction(
        left=[
            Participant(chebi_id=CHEBI_SUPEROXIDE),
            Participant(chebi_id=CHEBI_H_PLUS),
        ],
        right=[
            Participant(chebi_id=CHEBI_O2),
            Participant(chebi_id=CHEBI_H2O2),
        ],
    )
    assert ActingOnSuperoxideAsAcceptor().check_membership(reaction).is_member is True


def test_ec_33_positive_acting_on_ether_bonds() -> None:
    reaction = _reaction(
        left=[
            Participant(chebi_id=CHEBI_EPOXIDE, smiles="CC1CO1"),
            Participant(chebi_id=CHEBI_H2O, smiles="O"),
        ],
        right=[Participant(chebi_id=CHEBI_DIOL, smiles="CC(O)CO")],
    )
    assert ActingOnEtherBonds().check_membership(reaction).is_member is True


def test_ec_37_positive_acting_on_carbon_carbon_bonds() -> None:
    reaction = _reaction(
        left=[
            Participant(chebi_id=CHEBI_FUMARYLACETOACETATE),
            Participant(chebi_id=CHEBI_H2O),
        ],
        right=[Participant(name="fumarate"), Participant(name="acetoacetate")],
        label="fumarylacetoacetate + water = fumarate + acetoacetate",
    )
    assert ActingOnCarbonCarbonBonds().check_membership(reaction).is_member is True


def test_ec_43_positive_carbon_nitrogen_lyases() -> None:
    reaction = _reaction(
        left=[Participant(chebi_id=CHEBI_TYROSINE)],
        right=[Participant(name="trans-cinnamate"), Participant(chebi_id=CHEBI_NH3)],
        label="L-phenylalanine = trans-cinnamate + ammonia",
    )
    assert CarbonNitrogenLyases().check_membership(reaction).is_member is True


def test_ec_52_positive_cis_trans_isomerases() -> None:
    reaction = _reaction(
        left=[Participant(chebi_id=CHEBI_PEPTIDYL_PROLINE_CIS)],
        right=[Participant(chebi_id=CHEBI_PEPTIDYL_PROLINE_TRANS)],
    )
    assert CisTransIsomerase().check_membership(reaction).is_member is True


def test_ec_61_positive_forming_carbon_oxygen_bonds() -> None:
    reaction = _reaction(
        left=[
            Participant(polymer_type=PolymerType.TRNA, chebi_id=CHEBI_TRNA),
            Participant(chebi_id=CHEBI_TYROSINE),
            Participant(chebi_id=CHEBI_ATP),
        ],
        right=[
            Participant(polymer_type=PolymerType.TRNA, chebi_id=CHEBI_TRNA),
            Participant(chebi_id=CHEBI_AMP),
            Participant(chebi_id=CHEBI_DIPHOSPHATE),
        ],
    )
    assert FormingCarbonOxygenBonds().check_membership(reaction).is_member is True


def test_ec_71_positive_translocation_of_hydrons() -> None:
    reaction = _reaction(
        left=[
            Participant(chebi_id=CHEBI_CYTOCHROME_C_REDUCED),
            Participant(chebi_id=CHEBI_O2),
            Participant(chebi_id=CHEBI_H_PLUS),
        ],
        right=[
            Participant(chebi_id=CHEBI_CYTOCHROME_C_OXIDIZED),
            Participant(chebi_id=CHEBI_H2O),
        ],
    )
    assert TranslocationOfHydrons().check_membership(reaction).is_member is True


@pytest.mark.parametrize(
    ("cls",),
    [
        (ActingOnHalogenInDonors,),
        (EnzymesUsingH2AsReductant,),
        (OtherEnzymesUsingO2AsOxidant,),
        (NitrogenOxygenLyases,),
        (ATPIndependentChelatases,),
        (FormingNitrogenNitrogenBonds,),
    ],
)
def test_empty_level2_shells_stay_honest(cls: type) -> None:
    reaction = _reaction(
        left=[Participant(chebi_id=CHEBI_ATP)],
        right=[Participant(chebi_id=CHEBI_AMP)],
    )
    classifier = cls()
    assert classifier.supports_evaluation(reaction) is False
    assert classifier.check_membership(reaction).is_member is False
