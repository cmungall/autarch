"""Focused tests for batch 7 grouping classifiers."""

from __future__ import annotations

from typing_extensions import TypedDict

from autarch.datamodel import Participant, Reaction
from autarch.molecules import (
    CHEBI_CMP,
    CHEBI_FAD,
    CHEBI_FADH2,
    CHEBI_FMN,
    CHEBI_FMNH2,
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_H_PLUS,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_O2,
)
from autarch.ontology.oxidoreductase_acting_on_paired_donors_with_incorporation_or_reduction_of_molecular_oxygen_reduced_flavin_or_flavoprotein_as_one_donor_and_incorporation_of_one_atom_of_oxygen import (
    OxidoreductaseActingOnPairedDonorsWithIncorporationOrReductionOfMolecularOxygenReducedFlavinOrFlavoproteinAsOneDonorAndIncorporationOfOneAtomOfOxygen,
)
from autarch.ontology.oxidoreductase_acting_on_paired_donors_with_oxidation_of_a_pair_of_donors_resulting_in_the_reduction_of_molecular_oxygen_to_two_molecules_of_water import (
    OxidoreductaseActingOnPairedDonorsWithOxidationOfAPairOfDonorsResultingInTheReductionOfMolecularOxygenToTwoMoleculesOfWater,
)
from autarch.ontology.oxidoreductase_acting_on_the_ch_ch_group_of_donors_oxygen_as_acceptor import (
    OxidoreductaseActingOnTheCHCHGroupOfDonorsOxygenAsAcceptor,
)
from autarch.ontology.phosphotransferase_for_other_substituted_phosphate_groups import (
    PhosphotransferaseForOtherSubstitutedPhosphateGroups,
)
from autarch.ontology.racemase_and_epimerase_acting_on_carbohydrates_and_derivatives import (
    RacemaseAndEpimeraseActingOnCarbohydratesAndDerivatives,
)


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


def test_ec_11414_positive_flavin_linked_monooxygenation() -> None:
    cls = (
        OxidoreductaseActingOnPairedDonorsWithIncorporationOrReductionOfMolecularOxygenReducedFlavinOrFlavoproteinAsOneDonorAndIncorporationOfOneAtomOfOxygen()
    )
    reaction = _reaction(
        left=[
            {"smiles": "O=C([O-])Cc1ccc(O)cc1", "name": "4-hydroxyphenylacetate"},
            {"chebi_id": CHEBI_FMNH2, "name": "FMNH2"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
        ],
        right=[
            {"smiles": "O=C([O-])Cc1cc(O)c(O)cc1", "name": "3,4-dihydroxyphenylacetate"},
            {"chebi_id": CHEBI_FMN, "name": "FMN"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "one atom of oxygen" in result.explanation.lower()


def test_ec_11414_rejects_two_water_desaturase_branch() -> None:
    cls = (
        OxidoreductaseActingOnPairedDonorsWithIncorporationOrReductionOfMolecularOxygenReducedFlavinOrFlavoproteinAsOneDonorAndIncorporationOfOneAtomOfOxygen()
    )
    reaction = _reaction(
        left=[
            {"smiles": "CCCCCCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OCC1OC(C(O)C1OP(=O)([O-])[O-])n1cnc2c(N)ncnc12", "name": "octadecanoyl-CoA"},
            {"chebi_id": CHEBI_FADH2, "name": "FADH2"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron", "count": 2},
        ],
        right=[
            {"smiles": "CCCCCCCC/C=C\\CCCCCCCC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OCC1OC(C(O)C1OP(=O)([O-])[O-])n1cnc2c(N)ncnc12", "name": "octadecenoyl-CoA"},
            {"chebi_id": CHEBI_FAD, "name": "FAD"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water", "count": 2},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_11419_positive_desaturase_like_oxygen_reduction() -> None:
    cls = (
        OxidoreductaseActingOnPairedDonorsWithOxidationOfAPairOfDonorsResultingInTheReductionOfMolecularOxygenToTwoMoleculesOfWater()
    )
    reaction = _reaction(
        left=[
            {"smiles": "CCCCCCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OCC1OC(C(O)C1OP(=O)([O-])[O-])n1cnc2c(N)ncnc12", "name": "octadecanoyl-CoA"},
            {"chebi_id": CHEBI_FADH2, "name": "FADH2"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron", "count": 2},
        ],
        right=[
            {"smiles": "CCCCCCCC/C=C\\CCCCCCCC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OCC1OC(C(O)C1OP(=O)([O-])[O-])n1cnc2c(N)ncnc12", "name": "octadecenoyl-CoA"},
            {"chebi_id": CHEBI_FAD, "name": "FAD"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water", "count": 2},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "two molecules of water" in result.explanation.lower()


def test_ec_11419_rejects_flavin_monooxygenation() -> None:
    cls = (
        OxidoreductaseActingOnPairedDonorsWithOxidationOfAPairOfDonorsResultingInTheReductionOfMolecularOxygenToTwoMoleculesOfWater()
    )
    reaction = _reaction(
        left=[
            {"smiles": "c1ccc2ccccc2s1(=O)=O", "name": "dibenzothiophene sulfone"},
            {"chebi_id": CHEBI_FMNH2, "name": "FMNH2"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
        ],
        right=[
            {"smiles": "O=S([O-])c1ccccc1c1ccccc1O", "name": "2-hydroxybiphenyl-2-sulfinate"},
            {"chebi_id": CHEBI_FMN, "name": "FMN"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_513_positive_carbohydrate_epimerization() -> None:
    cls = RacemaseAndEpimeraseActingOnCarbohydratesAndDerivatives()
    reaction = _reaction(
        left=[
            {
                "smiles": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
                "inchi": "InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6+/m1/s1",
                "name": "alpha-D-glucose",
            }
        ],
        right=[
            {
                "smiles": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
                "inchi": "InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6-/m1/s1",
                "name": "beta-D-glucose",
            }
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "carbohydrate" in result.explanation.lower()


def test_ec_513_rejects_amino_acid_racemization() -> None:
    cls = RacemaseAndEpimeraseActingOnCarbohydratesAndDerivatives()
    reaction = _reaction(
        left=[
            {
                "smiles": "N[C@@H](C)C(=O)O",
                "inchi": "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1",
                "name": "L-alanine",
            }
        ],
        right=[
            {
                "smiles": "N[C@H](C)C(=O)O",
                "inchi": "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2+/m1/s1",
                "name": "D-alanine",
            }
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_278_positive_substituted_phosphate_transfer() -> None:
    cls = PhosphotransferaseForOtherSubstitutedPhosphateGroups()
    reaction = _reaction(
        left=[
            {
                "smiles": "NC1=NC(=O)N([C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])OCC[NH3+])[C@@H](O)[C@H]2O)C=C1",
                "name": "CDP-ethanolamine",
            },
            {"smiles": "[NH3+][C@@H](CO)C(=O)[O-]", "name": "L-serine"},
        ],
        right=[
            {"smiles": "[NH3+]CCOP(=O)([O-])OC[C@H]([NH3+])C(=O)[O-]", "name": "L-serine phosphoethanolamine"},
            {"chebi_id": CHEBI_CMP, "name": "CMP"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "substituted phosphate" in result.explanation.lower()


def test_ec_278_rejects_simple_kinase() -> None:
    cls = PhosphotransferaseForOtherSubstitutedPhosphateGroups()
    reaction = _reaction(
        left=[
            {"smiles": "C[N+](C)(C)CCO", "name": "choline"},
            {"smiles": "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]1O", "name": "ATP"},
        ],
        right=[
            {"smiles": "C[N+](C)(C)CCOP(=O)([O-])[O-]", "name": "phosphocholine"},
            {"smiles": "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]1O", "name": "ADP"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_133_positive_oxygen_acceptor_desaturation() -> None:
    cls = OxidoreductaseActingOnTheCHCHGroupOfDonorsOxygenAsAcceptor()
    reaction = _reaction(
        left=[
            {"smiles": "O=C1CCNC(=O)N1", "name": "5,6-dihydrouracil"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
        ],
        right=[
            {"smiles": "O=C1C=CNC(=O)N1", "name": "uracil"},
            {"chebi_id": CHEBI_H2O2, "smiles": "OO", "name": "hydrogen peroxide"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "oxygen as acceptor" in result.explanation.lower()


def test_ec_133_rejects_nadph_monooxygenation() -> None:
    cls = OxidoreductaseActingOnTheCHCHGroupOfDonorsOxygenAsAcceptor()
    reaction = _reaction(
        left=[
            {"smiles": "CC1=CCC(CC1)C(=C)C", "name": "limonene"},
            {"chebi_id": CHEBI_NADPH, "name": "NADPH"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
        right=[
            {"smiles": "CC1=CCC(CC1O)C(=C)C", "name": "limonenol"},
            {"chebi_id": CHEBI_NADP_PLUS, "name": "NADP(+)"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False
