"""Focused tests for the next five grouping classifiers."""

from __future__ import annotations

from typing_extensions import TypedDict

from autarch.datamodel import Participant, Reaction
from autarch.molecules import (
    CHEBI_FAD,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
)
from autarch.ontology.hydrolase_acting_on_carbon_nitrogen_but_not_peptide_bonds_in_cyclic_amides import (
    HydrolaseActingOnCarbonNitrogenButNotPeptideBondsInCyclicAmides,
)
from autarch.ontology.oxidoreductase_acting_on_paired_donors_with_incorporation_or_reduction_of_molecular_oxygen_nadh_or_nadph_as_one_donor_and_incorporation_of_two_atoms_of_oxygen_into_one_donor import (
    OxidoreductaseActingOnPairedDonorsWithIncorporationOrReductionOfMolecularOxygenNADHOrNADPHAsOneDonorAndIncorporationOfTwoAtomsOfOxygenIntoOneDonor,
)
from autarch.ontology.oxidoreductase_acting_on_paired_donors_with_incorporation_or_reduction_of_molecular_oxygen_reduced_iron_sulfur_protein_as_one_donor_and_incorporation_of_one_atom_of_oxygen import (
    OxidoreductaseActingOnPairedDonorsWithIncorporationOrReductionOfMolecularOxygenReducedIronSulfurProteinAsOneDonorAndIncorporationOfOneAtomOfOxygen,
)
from autarch.ontology.oxidoreductase_acting_on_the_ch_ch_group_of_donors_with_a_flavin_as_acceptor import (
    OxidoreductaseActingOnTheCHCHGroupOfDonorsWithAFlavinAsAcceptor,
)
from autarch.ontology.oxo_acid_lyase import OxoAcidLyase

CHEBI_ADP = "CHEBI:456216"
CHEBI_ATP = "CHEBI:30616"
CHEBI_PHOSPHATE = "CHEBI:16838"
CHEBI_IRON_II = "CHEBI:29033"
CHEBI_IRON_III = "CHEBI:29034"
CHEBI_REDUCED_ETF = "CHEBI:58307"

ATP_SMILES = (
    "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])"
    "[C@@H](O)[C@H]1O"
)
ADP_SMILES = (
    "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]1O"
)
FAD_SMILES = (
    "CC1=CC2=C(C=C1C)N(C[C@H](O)[C@H](O)[C@H](O)COP(=O)([O-])OP(=O)([O-])"
    "OC[C@H]1O[C@@H](N3C=NC4=C3N=CN=C4N)[C@H](O)[C@@H]1O)C1=NC(=O)[N-]C(=O)C1=N2"
)
REDUCED_ETF_SMILES = (
    "CC1=CC2=C(C=C1C)N(C[C@H](O)[C@H](O)[C@H](O)COP(=O)([O-])OP(=O)([O-])"
    "OC[C@H]1O[C@@H](N3C=NC4=C3N=CN=C4N)[C@H](O)[C@@H]1O)C1=C(N2)C(=O)NC(=O)N1"
)


class ParticipantInput(TypedDict, total=False):
    """Subset of participant fields used in classifier tests."""

    chebi_id: str
    smiles: str
    name: str
    count: int


def _reaction(left: list[ParticipantInput], right: list[ParticipantInput]) -> Reaction:
    return Reaction(
        left_participants=[Participant(**participant) for participant in left],
        right_participants=[Participant(**participant) for participant in right],
    )


def test_ec_352_positive_cyclic_amide_ring_opening() -> None:
    cls = HydrolaseActingOnCarbonNitrogenButNotPeptideBondsInCyclicAmides()
    reaction = _reaction(
        left=[
            {"smiles": "O=C1CCNC(=O)N1", "name": "5,6-dihydrouracil"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "NC(=O)NCCC(=O)[O-]", "name": "N-carbamoyl-beta-alaninate"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "cyclic amide" in result.explanation.lower()


def test_ec_352_positive_atp_dependent_lactam_hydrolysis() -> None:
    cls = HydrolaseActingOnCarbonNitrogenButNotPeptideBondsInCyclicAmides()
    reaction = _reaction(
        left=[
            {"smiles": "O=C1CC[C@@H](C(=O)[O-])N1", "name": "5-oxo-L-proline"},
            {"chebi_id": CHEBI_ATP, "smiles": ATP_SMILES, "name": "ATP"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water", "count": 2},
        ],
        right=[
            {"smiles": "[NH3+][C@@H](CCC(=O)[O-])C(=O)[O-]", "name": "L-glutamate"},
            {"chebi_id": CHEBI_ADP, "smiles": ADP_SMILES, "name": "ADP"},
            {"chebi_id": CHEBI_PHOSPHATE, "smiles": "O=P([O-])([O-])O", "name": "polyphosphate"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_352_rejects_cyclic_amidine_hydrolysis() -> None:
    cls = HydrolaseActingOnCarbonNitrogenButNotPeptideBondsInCyclicAmides()
    reaction = _reaction(
        left=[
            {"smiles": "CN1C(=N)N=C(C1=O)N", "name": "creatinine"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "CN(CC(=O)[O-])C(=N)N", "name": "creatine"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_138_positive_flavin_dependent_desaturation() -> None:
    cls = OxidoreductaseActingOnTheCHCHGroupOfDonorsWithAFlavinAsAcceptor()
    reaction = _reaction(
        left=[
            {
                "smiles": "CC(C)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])"
                "OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-]",
                "name": "isovaleryl-CoA",
            },
            {"chebi_id": CHEBI_FAD, "smiles": FAD_SMILES, "name": "FAD"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
        right=[
            {
                "smiles": "CC(C)=CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])"
                "OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-]",
                "name": "3-methylbut-2-enoyl-CoA",
            },
            {"chebi_id": CHEBI_REDUCED_ETF, "smiles": REDUCED_ETF_SMILES, "name": "reduced ETF"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "flavin" in result.explanation.lower()


def test_ec_138_rejects_nad_dependent_ch_ch_redox() -> None:
    cls = OxidoreductaseActingOnTheCHCHGroupOfDonorsWithAFlavinAsAcceptor()
    reaction = _reaction(
        left=[
            {
                "smiles": "CC(C)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])"
                "OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-]",
                "name": "isovaleryl-CoA",
            },
            {"chebi_id": CHEBI_NAD_PLUS, "name": "NAD(+)"},
        ],
        right=[
            {
                "smiles": "CC(C)=CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])"
                "OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-]",
                "name": "3-methylbut-2-enoyl-CoA",
            },
            {"chebi_id": CHEBI_NADH, "name": "NADH"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_11412_positive_nad_linked_dioxygenation() -> None:
    cls = (
        OxidoreductaseActingOnPairedDonorsWithIncorporationOrReductionOfMolecularOxygenNADHOrNADPHAsOneDonorAndIncorporationOfTwoAtomsOfOxygenIntoOneDonor()
    )
    reaction = _reaction(
        left=[
            {"smiles": "O=C([O-])C1=CC=CC=C1", "name": "benzoate"},
            {"chebi_id": CHEBI_NADH, "name": "NADH"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
        right=[
            {"smiles": "O=C([O-])[C@@]1(O)C=CC=C[C@@H]1O", "name": "cis-dihydrodiol"},
            {"chebi_id": CHEBI_NAD_PLUS, "name": "NAD(+)"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "two atoms of oxygen" in result.explanation.lower()


def test_ec_11412_rejects_ferredoxin_branch() -> None:
    cls = (
        OxidoreductaseActingOnPairedDonorsWithIncorporationOrReductionOfMolecularOxygenNADHOrNADPHAsOneDonorAndIncorporationOfTwoAtomsOfOxygenIntoOneDonor()
    )
    reaction = _reaction(
        left=[
            {"smiles": "CCCCCCCC", "name": "octane"},
            {"chebi_id": CHEBI_IRON_II, "smiles": "[Fe+2]", "name": "iron(2+)", "count": 2},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron", "count": 2},
        ],
        right=[
            {"chebi_id": CHEBI_IRON_III, "smiles": "[Fe+3]", "name": "iron(3+)", "count": 2},
            {"smiles": "CCCCCCCCO", "name": "octan-1-ol"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_413_positive_oxo_acid_lyase() -> None:
    cls = OxoAcidLyase()
    reaction = _reaction(
        left=[
            {"smiles": "O=C([O-])CC(O)(CC(=O)[O-])C(=O)[O-]", "name": "citrate"},
        ],
        right=[
            {"smiles": "O=C([O-])CC(=O)C(=O)[O-]", "name": "oxaloacetate"},
            {"smiles": "CC(=O)[O-]", "name": "acetate"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "oxo-acid lyase" in result.explanation.lower()


def test_ec_413_rejects_decarboxylase() -> None:
    cls = OxoAcidLyase()
    reaction = _reaction(
        left=[
            {"smiles": "CC(=O)C(=O)[O-]", "name": "pyruvate"},
        ],
        right=[
            {"smiles": "CC=O", "name": "acetaldehyde"},
            {"smiles": "O=C=O", "name": "carbon dioxide"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_11415_positive_iron_sulfur_monooxygenase() -> None:
    cls = (
        OxidoreductaseActingOnPairedDonorsWithIncorporationOrReductionOfMolecularOxygenReducedIronSulfurProteinAsOneDonorAndIncorporationOfOneAtomOfOxygen()
    )
    reaction = _reaction(
        left=[
            {"smiles": "CCCCCCCC", "name": "octane"},
            {"chebi_id": CHEBI_IRON_II, "smiles": "[Fe+2]", "name": "iron(2+)", "count": 2},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron", "count": 2},
        ],
        right=[
            {"chebi_id": CHEBI_IRON_III, "smiles": "[Fe+3]", "name": "iron(3+)", "count": 2},
            {"smiles": "CCCCCCCCO", "name": "octan-1-ol"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "iron-sulfur" in result.explanation.lower()


def test_ec_11415_rejects_nadph_monooxygenase() -> None:
    cls = (
        OxidoreductaseActingOnPairedDonorsWithIncorporationOrReductionOfMolecularOxygenReducedIronSulfurProteinAsOneDonorAndIncorporationOfOneAtomOfOxygen()
    )
    reaction = _reaction(
        left=[
            {"smiles": "CCCCCCCC", "name": "octane"},
            {"chebi_id": CHEBI_NADPH, "name": "NADPH"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
        right=[
            {"smiles": "CCCCCCCCO", "name": "octan-1-ol"},
            {"chebi_id": CHEBI_NADP_PLUS, "name": "NADP(+)"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False
