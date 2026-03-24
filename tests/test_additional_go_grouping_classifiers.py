"""Focused tests for additional GO-aligned grouping classifiers."""

from typing_extensions import TypedDict

from autarch.datamodel import Participant, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_ATP,
    CHEBI_CMP,
    CHEBI_FMN,
    CHEBI_FMNH2,
    CHEBI_DIPHOSPHATE,
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADP_PLUS,
    CHEBI_NADPH,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
    CHEBI_PHOSPHATE,
    CHEBI_UDP,
)
from autarch.ontology.acid_amino_acid_ligase import AcidAminoAcidLigase
from autarch.ontology.acyltransferase_transferring_groups_other_than_amino_acyl_groups import (
    AcyltransferaseTransferringGroupsOtherThanAminoAcylGroups,
)
from autarch.ontology.hexosyltransferase import Hexosyltransferase
from autarch.ontology.hydrolase_acting_on_carbon_nitrogen_but_not_peptide_bonds_in_linear_amides import (
    HydrolaseActingOnCarbonNitrogenButNotPeptideBondsInLinearAmides,
)
from autarch.ontology.intramolecular_oxidoreductase_interconverting_aldoses_and_ketoses import (
    IntramolecularOxidoreductaseInterconvertingAldosesAndKetoses,
)
from autarch.ontology.oxidoreductase_acting_on_paired_donors_with_incorporation_or_reduction_of_molecular_oxygen_nadh_or_nadph_as_one_donor_and_incorporation_of_one_atom_of_oxygen import (
    OxidoreductaseActingOnPairedDonorsWithIncorporationOrReductionOfMolecularOxygenNADHOrNADPHAsOneDonorAndIncorporationOfOneAtomOfOxygen,
)
from autarch.ontology.oxidoreductase_acting_on_single_donors_with_incorporation_of_molecular_oxygen_incorporation_of_two_atoms_of_oxygen import (
    OxidoreductaseActingOnSingleDonorsWithIncorporationOfMolecularOxygenIncorporationOfTwoAtomsOfOxygen,
)
from autarch.ontology.oxidoreductase_acting_on_the_ch_nh_group_of_donors_nad_or_nadp_as_acceptor import (
    OxidoreductaseActingOnTheCHNHGroupOfDonorsNADOrNADPAsAcceptor,
)
from autarch.ontology.pentosyltransferase import Pentosyltransferase
from autarch.ontology.thiolester_hydrolase import ThiolesterHydrolase


class ParticipantInput(TypedDict, total=False):
    """Subset of participant fields used in classifier tests."""

    chebi_id: str
    smiles: str
    name: str
    count: int


def _reaction(
    left: list[ParticipantInput],
    right: list[ParticipantInput],
) -> Reaction:
    return Reaction(
        left_participants=[Participant(**participant) for participant in left],
        right_participants=[Participant(**participant) for participant in right],
    )


def test_ec_11413_positive_nadph_monooxygenase() -> None:
    cls = (
        OxidoreductaseActingOnPairedDonorsWithIncorporationOrReductionOfMolecularOxygenNADHOrNADPHAsOneDonorAndIncorporationOfOneAtomOfOxygen()
    )
    reaction = _reaction(
        left=[
            {"smiles": "Cc1ccccc1", "name": "toluene"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
            {"chebi_id": CHEBI_NADPH, "name": "NADPH"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
        right=[
            {"smiles": "OCc1ccccc1", "name": "benzyl alcohol"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
            {"chebi_id": CHEBI_NADP_PLUS, "name": "NADP(+)"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "one atom of oxygen" in result.explanation.lower()


def test_ec_11413_rejects_dioxygenase_without_nadph() -> None:
    cls = (
        OxidoreductaseActingOnPairedDonorsWithIncorporationOrReductionOfMolecularOxygenNADHOrNADPHAsOneDonorAndIncorporationOfOneAtomOfOxygen()
    )
    reaction = _reaction(
        left=[
            {"smiles": "Oc1ccccc1O", "name": "catechol"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
        ],
        right=[
            {"smiles": "O=C(O)C=CC=CC(=O)O", "name": "muconic acid"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_11413_accepts_reduced_flavin_monooxygenase_proxy() -> None:
    cls = (
        OxidoreductaseActingOnPairedDonorsWithIncorporationOrReductionOfMolecularOxygenNADHOrNADPHAsOneDonorAndIncorporationOfOneAtomOfOxygen()
    )
    reaction = _reaction(
        left=[
            {"smiles": "O=C([O-])/C=C/C1=CC=CC=C1", "name": "trans-cinnamate"},
            {
                "chebi_id": "CHEBI:57618",
                "smiles": "CC1=C(C)C=C2C(=C1)NC1=C(NC(=O)NC1=O)N2C[C@H](O)[C@H](O)[C@H](O)COP(=O)([O-])[O-]",
                "name": "reduced hemoprotein reductase proxy",
            },
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
        ],
        right=[
            {"smiles": "O=C([O-])/C=C/C1=CC=C(O)C=C1", "name": "p-coumarate"},
            {
                "chebi_id": "CHEBI:58210",
                "smiles": "CC1=C(C)C=C2C(=C1)/N=C1/C(=O)[N-]C(=O)N=C1N2C[C@H](O)[C@H](O)[C@H](O)COP(=O)([O-])[O-]",
                "name": "oxidized hemoprotein reductase proxy",
            },
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_11413_accepts_methane_monooxygenase() -> None:
    cls = (
        OxidoreductaseActingOnPairedDonorsWithIncorporationOrReductionOfMolecularOxygenNADHOrNADPHAsOneDonorAndIncorporationOfOneAtomOfOxygen()
    )
    reaction = _reaction(
        left=[
            {"smiles": "C", "name": "methane"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
            {"chebi_id": CHEBI_NADH, "name": "NADH"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
        right=[
            {"smiles": "CO", "name": "methanol"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
            {"chebi_id": CHEBI_NAD_PLUS, "name": "NAD(+)"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_11413_rejects_nadh_dioxygenase_without_water() -> None:
    cls = (
        OxidoreductaseActingOnPairedDonorsWithIncorporationOrReductionOfMolecularOxygenNADHOrNADPHAsOneDonorAndIncorporationOfOneAtomOfOxygen()
    )
    reaction = _reaction(
        left=[
            {"smiles": "c1ccccc1", "name": "benzene"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
            {"chebi_id": CHEBI_NADH, "name": "NADH"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
        right=[
            {"smiles": "O[C@@H]1C=CC=C[C@@H]1O", "name": "cis-dihydrodiol"},
            {"chebi_id": CHEBI_NAD_PLUS, "name": "NAD(+)"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_11311_positive_dioxygenase() -> None:
    cls = (
        OxidoreductaseActingOnSingleDonorsWithIncorporationOfMolecularOxygenIncorporationOfTwoAtomsOfOxygen()
    )
    reaction = _reaction(
        left=[
            {"smiles": "Oc1ccccc1O", "name": "catechol"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
        ],
        right=[
            {"smiles": "O=C(O)C=CC=CC(=O)O", "name": "muconic acid"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "two atoms of oxygen" in result.explanation.lower()


def test_ec_11311_rejects_nadph_monooxygenase() -> None:
    cls = (
        OxidoreductaseActingOnSingleDonorsWithIncorporationOfMolecularOxygenIncorporationOfTwoAtomsOfOxygen()
    )
    reaction = _reaction(
        left=[
            {"smiles": "Cc1ccccc1", "name": "toluene"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
            {"chebi_id": CHEBI_NADPH, "name": "NADPH"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
        right=[
            {"smiles": "OCc1ccccc1", "name": "benzyl alcohol"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
            {"chebi_id": CHEBI_NADP_PLUS, "name": "NADP(+)"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_151_positive_scaffold_preserving_n_redox() -> None:
    cls = OxidoreductaseActingOnTheCHNHGroupOfDonorsNADOrNADPAsAcceptor()
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

    assert result.is_member is True
    assert "ch-nh" in result.explanation.lower()


def test_ec_151_positive_flavin_reductase_mode() -> None:
    cls = OxidoreductaseActingOnTheCHNHGroupOfDonorsNADOrNADPAsAcceptor()
    reaction = _reaction(
        left=[
            {
                "chebi_id": CHEBI_FMNH2,
                "smiles": "CC1=C(C)C=C2C(=C1)NC1=C(NC(=O)NC1=O)N2C[C@H](O)[C@H](O)[C@H](O)COP(=O)([O-])[O-]",
                "name": "FMNH2",
            },
            {"chebi_id": CHEBI_NAD_PLUS, "name": "NAD(+)"},
        ],
        right=[
            {
                "chebi_id": CHEBI_FMN,
                "smiles": "CC1=C(C)C=C2C(=C1)/N=C1/C(=O)[N-]C(=O)N=C1N2C[C@H](O)[C@H](O)[C@H](O)COP(=O)([O-])[O-]",
                "name": "FMN",
            },
            {"chebi_id": CHEBI_NADH, "name": "NADH"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron", "count": 2},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_151_positive_water_assisted_opine_cleavage() -> None:
    cls = OxidoreductaseActingOnTheCHNHGroupOfDonorsNADOrNADPAsAcceptor()
    reaction = _reaction(
        left=[
            {
                "smiles": "C[C@@H]([NH2+]CC(=O)[O-])C(=O)[O-]",
                "name": "N-(carboxymethyl)-D-alanine",
            },
            {"chebi_id": CHEBI_NAD_PLUS, "name": "NAD(+)"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "[NH3+]CC(=O)[O-]", "name": "glycine"},
            {"smiles": "CC(=O)C(=O)[O-]", "name": "pyruvate"},
            {"chebi_id": CHEBI_NADH, "name": "NADH"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_151_rejects_ch_oh_dehydrogenase_with_incidental_nitrogen() -> None:
    cls = OxidoreductaseActingOnTheCHNHGroupOfDonorsNADOrNADPAsAcceptor()
    reaction = _reaction(
        left=[
            {"smiles": "OC1CC(N)CC1C(=O)O", "name": "4-hydroxy-L-proline"},
            {"chebi_id": CHEBI_NAD_PLUS, "name": "NAD(+)"},
        ],
        right=[
            {"smiles": "O=C1CC(N)CC1C(=O)O", "name": "4-oxo-L-proline"},
            {"chebi_id": CHEBI_NADH, "name": "NADH"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_151_rejects_o2_dependent_amine_oxidase() -> None:
    cls = OxidoreductaseActingOnTheCHNHGroupOfDonorsNADOrNADPAsAcceptor()
    reaction = _reaction(
        left=[
            {"smiles": "CNCC(=O)[O-]", "name": "sarcosine"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "[NH3+]CC(=O)[O-]", "name": "glycine"},
            {"smiles": "C=O", "name": "formaldehyde"},
            {"chebi_id": CHEBI_H2O2, "smiles": "OO", "name": "hydrogen peroxide"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_531_positive_cyclic_glucose6p_to_fructose6p() -> None:
    cls = IntramolecularOxidoreductaseInterconvertingAldosesAndKetoses()
    reaction = _reaction(
        left=[
            {
                "smiles": "O=P([O-])([O-])OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
                "name": "alpha-D-glucose 6-phosphate",
            }
        ],
        right=[
            {
                "smiles": "O=P([O-])([O-])OC[C@H]1O[C@](O)(CO)[C@@H](O)[C@@H]1O",
                "name": "beta-D-fructose 6-phosphate",
            }
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "aldose-ketose" in result.explanation.lower()


def test_ec_531_positive_triose_phosphate_isomerase() -> None:
    cls = IntramolecularOxidoreductaseInterconvertingAldosesAndKetoses()
    reaction = _reaction(
        left=[
            {
                "smiles": "[H]C(=O)[C@H](O)COP(=O)([O-])[O-]",
                "name": "D-glyceraldehyde 3-phosphate",
            }
        ],
        right=[
            {
                "smiles": "O=C(CO)COP(=O)([O-])[O-]",
                "name": "dihydroxyacetone phosphate",
            }
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_531_positive_modified_ribose_isomerization() -> None:
    cls = IntramolecularOxidoreductaseInterconvertingAldosesAndKetoses()
    reaction = _reaction(
        left=[
            {
                "smiles": "NC(=O)C1=C(/N=C/N[C@@H]2O[C@H](COP(=O)([O-])[O-])[C@@H](O)[C@H]2O)N([C@@H]2O[C@H](COP(=O)([O-])[O-])[C@@H](O)[C@H]2O)C=N1",
                "name": "ribose-containing imidazole precursor",
            }
        ],
        right=[
            {
                "smiles": "NC(=O)C1=C(NC=NCC(=O)[C@H](O)[C@H](O)COP(=O)([O-])[O-])N([C@@H]2O[C@H](COP(=O)([O-])[O-])[C@@H](O)[C@H]2O)C=N1",
                "name": "ribulose-containing imidazole product",
            }
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_531_rejects_aldose_epimerization() -> None:
    cls = IntramolecularOxidoreductaseInterconvertingAldosesAndKetoses()
    reaction = _reaction(
        left=[
            {
                "smiles": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
                "name": "alpha-D-glucose",
            }
        ],
        right=[
            {
                "smiles": "OC[C@H]1O[C@@H](O)[C@@H](O)[C@H](O)[C@@H]1O",
                "name": "alpha-D-mannose",
            }
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_531_rejects_isomerase_with_external_cofactor() -> None:
    cls = IntramolecularOxidoreductaseInterconvertingAldosesAndKetoses()
    reaction = _reaction(
        left=[
            {
                "smiles": "[H]C(=O)[C@H](O)COP(=O)([O-])[O-]",
                "name": "D-glyceraldehyde 3-phosphate",
            },
            {"chebi_id": CHEBI_NAD_PLUS, "name": "NAD(+)"},
        ],
        right=[
            {
                "smiles": "O=C(CO)COP(=O)([O-])[O-]",
                "name": "dihydroxyacetone phosphate",
            },
            {"chebi_id": CHEBI_NADH, "name": "NADH"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_351_positive_simple_amidase() -> None:
    cls = HydrolaseActingOnCarbonNitrogenButNotPeptideBondsInLinearAmides()
    reaction = _reaction(
        left=[
            {"smiles": "CCCCC(N)=O", "name": "pentanamide"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "CCCCC(=O)[O-]", "name": "pentanoate"},
            {"smiles": "[NH4+]", "name": "ammonium"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "amide" in result.explanation.lower()


def test_ec_351_positive_side_chain_deamidation() -> None:
    cls = HydrolaseActingOnCarbonNitrogenButNotPeptideBondsInLinearAmides()
    reaction = _reaction(
        left=[
            {"smiles": "N[C@@H](CC(N)=O)C(=O)O", "name": "L-asparagine"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "N[C@@H](CC(=O)O)C(=O)O", "name": "L-aspartate"},
            {"smiles": "[NH4+]", "name": "ammonium"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_351_positive_deacetylation() -> None:
    cls = HydrolaseActingOnCarbonNitrogenButNotPeptideBondsInLinearAmides()
    reaction = _reaction(
        left=[
            {
                "smiles": "CC(=O)N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O",
                "name": "N-acetylglucosamine",
            },
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {
                "smiles": "[NH3+][C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O",
                "name": "glucosamine",
            },
            {"smiles": "CC(=O)[O-]", "name": "acetate"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_351_rejects_peptide_hydrolysis() -> None:
    cls = HydrolaseActingOnCarbonNitrogenButNotPeptideBondsInLinearAmides()
    reaction = _reaction(
        left=[
            {"smiles": "N[C@@H](C)C(=O)N[C@@H](C)C(=O)O", "name": "alanylalanine"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "N[C@@H](C)C(=O)O", "name": "alanine", "count": 2},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_351_rejects_ester_hydrolysis() -> None:
    cls = HydrolaseActingOnCarbonNitrogenButNotPeptideBondsInLinearAmides()
    reaction = _reaction(
        left=[
            {"smiles": "CC(=O)OC", "name": "methyl acetate"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "CC(=O)[O-]", "name": "acetate"},
            {"smiles": "CO", "name": "methanol"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_242_positive_nucleoside_phosphorolysis() -> None:
    cls = Pentosyltransferase()
    reaction = _reaction(
        left=[
            {
                "smiles": "NC1=NC2=C(N=CN2[C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)C(=O)N1",
                "name": "guanosine",
            },
            {"smiles": "O=P([O-])([O-])O", "name": "phosphate"},
        ],
        right=[
            {
                "smiles": "O=P([O-])([O-])O[C@H]1O[C@H](CO)[C@@H](O)[C@H]1O",
                "name": "alpha-D-ribose 1-phosphate",
            },
            {
                "smiles": "NC1=NC2=C(N=CN2)C(=O)N1",
                "name": "guanine",
            },
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "pentose" in result.explanation.lower()


def test_ec_242_positive_prpp_exchange() -> None:
    cls = Pentosyltransferase()
    reaction = _reaction(
        left=[
            {
                "smiles": "O=C1C=C(C(=O)[O-])N([C@@H]2O[C@H](COP(=O)([O-])[O-])[C@@H](O)[C@H]2O)C(=O)N1",
                "name": "orotidine 5'-phosphate",
            },
            {"smiles": "O=P([O-])([O-])OP(=O)([O-])O", "name": "diphosphate"},
        ],
        right=[
            {"smiles": "O=C1C=C(C(=O)[O-])NC(=O)N1", "name": "orotate"},
            {
                "smiles": "O=P([O-])([O-])OC[C@H]1O[C@H](OP(=O)([O-])OP(=O)([O-])[O-])[C@H](O)[C@@H]1O",
                "name": "PRPP",
            },
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_242_positive_acceptor_exchange() -> None:
    cls = Pentosyltransferase()
    reaction = _reaction(
        left=[
            {"smiles": "CC1=C/C2=C(\\C=C/1C)NC=N2", "name": "5,6-dimethylbenzimidazole"},
            {
                "smiles": "O=C([O-])C1=CC=C[N+]([C@@H]2O[C@H](COP(=O)([O-])[O-])[C@@H](O)[C@H]2O)=C1",
                "name": "nicotinate beta-D-ribonucleotide",
            },
        ],
        right=[
            {
                "smiles": "CC1=CC2=C(C=C1C)N([C@H]1O[C@H](COP(=O)([O-])[O-])[C@@H](O)[C@H]1O)C=N2",
                "name": "alpha-ribazole 5'-phosphate",
            },
            {"smiles": "O=C([O-])C1=CC=CN=C1", "name": "nicotinate"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_242_rejects_hexose_glycosyltransferase() -> None:
    cls = Pentosyltransferase()
    reaction = _reaction(
        left=[
            {
                "smiles": "O=P([O-])([O-])OP(=O)([O-])OC[C@H]1O[C@@H](n2ccc(N)nc2=O)[C@H](O)[C@@H]1O",
                "name": "UDP-glucose proxy",
            },
            {"smiles": "CO", "name": "methanol"},
        ],
        right=[
            {"smiles": "CO[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O", "name": "methyl glucoside"},
            {"chebi_id": CHEBI_UDP, "name": "UDP"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_242_rejects_nad_dependent_adp_ribosylation_proxy() -> None:
    cls = Pentosyltransferase()
    reaction = _reaction(
        left=[
            {
                "chebi_id": CHEBI_NAD_PLUS,
                "name": "NAD(+)",
            },
            {"smiles": "NCCS", "name": "simple amine acceptor"},
        ],
        right=[
            {"smiles": "NCCN([C@H]1O[C@H](COP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]1O)S", "name": "ADP-ribosylated adduct"},
            {"smiles": "NC1=CC=CC=C1", "name": "nicotinamide proxy"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_231_positive_thioester_to_amine_transfer() -> None:
    cls = AcyltransferaseTransferringGroupsOtherThanAminoAcylGroups()
    reaction = _reaction(
        left=[
            {"smiles": "CC(=O)SC", "name": "methyl thioacetate"},
            {"smiles": "NCCO", "name": "ethanolamine"},
        ],
        right=[
            {"smiles": "CC(=O)NCCO", "name": "N-acetylethanolamine"},
            {"smiles": "CS", "name": "methanethiol"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "acyltransferase" in result.explanation.lower()


def test_ec_231_positive_acyl_phosphate_transfer() -> None:
    cls = AcyltransferaseTransferringGroupsOtherThanAminoAcylGroups()
    reaction = _reaction(
        left=[
            {"smiles": "CC(=O)OP(=O)([O-])[O-]", "name": "acetyl phosphate"},
            {"smiles": "NCCCC[C@H](N)C(=O)O", "name": "lysine"},
        ],
        right=[
            {"smiles": "CC(=O)NCCCC[C@H](N)C(=O)O", "name": "N-acetyllysine"},
            {"chebi_id": CHEBI_PHOSPHATE, "name": "phosphate"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_231_positive_acyl_sugar_transfer() -> None:
    cls = AcyltransferaseTransferringGroupsOtherThanAminoAcylGroups()
    reaction = _reaction(
        left=[
            {
                "smiles": "CC(=O)O[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O",
                "name": "acetyl glucose",
            },
            {"smiles": "C[N+](C)(C)CCO", "name": "choline"},
        ],
        right=[
            {"smiles": "CC(=O)OCC[N+](C)(C)C", "name": "acetylcholine"},
            {"smiles": "OC[C@H]1O[C@H](O)[C@@H](O)[C@H](O)[C@H]1O", "name": "glucose"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_231_rejects_hydrolysis() -> None:
    cls = AcyltransferaseTransferringGroupsOtherThanAminoAcylGroups()
    reaction = _reaction(
        left=[
            {"smiles": "CC(=O)SC", "name": "methyl thioacetate"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "CC(=O)[O-]", "name": "acetate"},
            {"smiles": "CS", "name": "methanethiol"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_231_rejects_trans_thioester_transfer() -> None:
    cls = AcyltransferaseTransferringGroupsOtherThanAminoAcylGroups()
    reaction = _reaction(
        left=[
            {"smiles": "CC(=O)SC", "name": "methyl thioacetate"},
            {"smiles": "CCC(=O)[O-]", "name": "propionate"},
        ],
        right=[
            {"smiles": "CCC(=O)SC", "name": "methyl thiopropionate"},
            {"smiles": "CC(=O)[O-]", "name": "acetate"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_632_positive_acid_amino_acid_ligase() -> None:
    cls = AcidAminoAcidLigase()
    reaction = _reaction(
        left=[
            {"smiles": "NCC(=O)O", "name": "glycine"},
            {"smiles": "NCC(=O)O", "name": "glycine"},
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
        ],
        right=[
            {"smiles": "NCC(=O)NCC(=O)O", "name": "glycylglycine"},
            {"chebi_id": CHEBI_ADP, "name": "ADP"},
            {"chebi_id": CHEBI_PHOSPHATE, "name": "phosphate"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "acid-amino acid ligase" in result.explanation.lower()


def test_ec_632_rejects_ammonia_ligase() -> None:
    cls = AcidAminoAcidLigase()
    reaction = _reaction(
        left=[
            {"smiles": "N[C@@H](CCC(=O)O)C(=O)O", "name": "glutamate"},
            {"smiles": "N", "name": "ammonia"},
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
        ],
        right=[
            {"smiles": "N[C@@H](CCC(N)=O)C(=O)O", "name": "glutamine"},
            {"chebi_id": CHEBI_ADP, "name": "ADP"},
            {"chebi_id": CHEBI_PHOSPHATE, "name": "phosphate"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_632_accepts_acid_amine_ligation_when_amino_acid_is_acyl_donor() -> None:
    cls = AcidAminoAcidLigase()
    reaction = _reaction(
        left=[
            {"smiles": "[NH3+]CCC1=CNC=N1", "name": "histamine"},
            {
                "smiles": "[NH3+][C@@H](CCC(=O)[O-])C(=O)[O-]",
                "name": "L-glutamate",
            },
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
        ],
        right=[
            {
                "smiles": "[NH3+][C@@H](CCC(=O)NCCC1=CNC=N1)C(=O)[O-]",
                "name": "gamma-glutamylhistamine",
            },
            {"chebi_id": CHEBI_ADP, "name": "ADP"},
            {"chebi_id": CHEBI_PHOSPHATE, "name": "phosphate"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_632_accepts_ctp_coupled_acid_amino_acid_ligation() -> None:
    cls = AcidAminoAcidLigase()
    reaction = _reaction(
        left=[
            {
                "smiles": "CC(C)(COP(=O)([O-])[O-])[C@@H](O)C(=O)NCCC(=O)[O-]",
                "name": "phosphopantothenate",
            },
            {
                "smiles": "[NH3+][C@@H](CS)C(=O)[O-]",
                "name": "L-cysteine",
            },
            {
                "chebi_id": "CHEBI:37563",
                "smiles": "NC1=NC(=O)N([C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]2O)C=C1",
                "name": "CTP",
            },
        ],
        right=[
            {
                "smiles": "CC(C)(COP(=O)([O-])[O-])[C@@H](O)C(=O)NCCC(=O)N[C@@H](CS)C(=O)[O-]",
                "name": "phosphopantothenoylcysteine",
            },
            {
                "chebi_id": CHEBI_CMP,
                "smiles": "NC1=NC(=O)N([C@@H]2O[C@H](COP(=O)([O-])[O-])[C@@H](O)[C@H]2O)C=C1",
                "name": "CMP",
            },
            {"chebi_id": CHEBI_DIPHOSPHATE, "name": "diphosphate"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_312_positive_thiolester_hydrolase() -> None:
    cls = ThiolesterHydrolase()
    reaction = _reaction(
        left=[
            {"smiles": "CC(=O)SC", "name": "methyl thioacetate"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "CC(=O)O", "name": "acetic acid"},
            {"smiles": "CS", "name": "methanethiol"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "thiolester hydrolase" in result.explanation.lower()


def test_ec_312_rejects_oxygen_ester_hydrolysis() -> None:
    cls = ThiolesterHydrolase()
    reaction = _reaction(
        left=[
            {"smiles": "CC(=O)OC", "name": "methyl acetate"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "CC(=O)O", "name": "acetic acid"},
            {"smiles": "CO", "name": "methanol"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_241_positive_hexosyltransferase() -> None:
    cls = Hexosyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": "CHEBI:18066", "name": "UDP-D-glucose"},
            {"smiles": "Oc1ccccc1", "name": "phenol"},
        ],
        right=[
            {"chebi_id": CHEBI_UDP, "name": "UDP"},
            {
                "smiles": "OC[C@H]1O[C@@H](Oc2ccccc2)[C@H](O)[C@@H](O)[C@H]1O",
                "name": "phenyl glucoside",
            },
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "hexosyltransferase" in result.explanation.lower()


def test_ec_241_rejects_non_hexose_nucleotide_sugar() -> None:
    cls = Hexosyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": "CHEBI:57812", "name": "CMP-sialic acid"},
            {"smiles": "Oc1ccccc1", "name": "phenol"},
        ],
        right=[
            {"chebi_id": CHEBI_CMP, "name": "CMP"},
            {"smiles": "Oc1ccccc1", "name": "acceptor-derived product"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_241_accepts_hexose_phosphorylase_branch() -> None:
    cls = Hexosyltransferase()
    reaction = _reaction(
        left=[
            {"smiles": "O[C@H]1[C@H](O)O[C@@H](CO)[C@@H](O)[C@H]1O[C@H]1O[C@@H](CO)[C@@H](O)[C@H](O)[C@H]1O", "name": "cellobiose"},
            {"chebi_id": CHEBI_PHOSPHATE, "name": "phosphate"},
        ],
        right=[
            {
                "smiles": "O[P](=O)([O-])O[C@H]1[C@H](O)O[C@@H](CO)[C@@H](O)[C@H]1O",
                "name": "glucose 1-phosphate",
            },
            {"smiles": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O", "name": "glucose"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_241_rejects_phosphoglycosyl_transferase() -> None:
    cls = Hexosyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": "CHEBI:18307", "name": "UDP-galactose"},
            {"smiles": "CCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCOP(=O)([O-])[O-]", "name": "undecaprenyl phosphate"},
        ],
        right=[
            {"smiles": "CCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCOP(=O)([O-])OP(=O)([O-])O[C@H]1O[C@@H](CO)[C@H](O)[C@@H](O)[C@H]1O", "name": "galactosyl undecaprenyl diphosphate"},
            {"chebi_id": "CHEBI:57865", "name": "UMP"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False
