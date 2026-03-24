"""Focused tests for additional level-3 grouping classifiers."""

from typing_extensions import TypedDict

from autarch.datamodel import Participant, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_AMP,
    CHEBI_ATP,
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_NH4,
    CHEBI_O2,
    CHEBI_PHOSPHATE,
)
from autarch.ontology.acid_ammonia_or_amide_ligase import AcidAmmoniaOrAmideLigase
from autarch.ontology.coa_transferase import CoATransferase
from autarch.ontology.hydrolase_acting_on_carbon_nitrogen_but_not_peptide_bonds_in_linear_amidines import (
    HydrolaseActingOnCarbonNitrogenButNotPeptideBondsInLinearAmidines,
)
from autarch.ontology.oxidoreductase_acting_on_a_sulfur_group_of_donors_nad_p_as_acceptor import (
    OxidoreductaseActingOnASulfurGroupOfDonorsNADPAsAcceptor,
)
from autarch.ontology.oxidoreductase_acting_on_the_ch_nh2_group_of_donors_oxygen_as_acceptor import (
    OxidoreductaseActingOnTheCHNH2GroupOfDonorsOxygenAsAcceptor,
)
from autarch.ontology.oxidoreductase_acting_on_the_ch_nh_group_of_donors_oxygen_as_acceptor import (
    OxidoreductaseActingOnTheCHNHGroupOfDonorsOxygenAsAcceptor,
)
from autarch.ontology.phosphotransferase_phosphate_group_as_acceptor import (
    PhosphotransferasePhosphateGroupAsAcceptor,
)
from autarch.ontology.racemase_and_epimerase_acting_on_amino_acids_and_derivatives import (
    RacemaseAndEpimeraseActingOnAminoAcidsAndDerivatives,
)
from autarch.ontology.sulfurtransferase import Sulfurtransferase
from autarch.ontology.sulfuric_ester_hydrolase import SulfuricEsterHydrolase


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


def test_ec_283_positive_coa_transfer() -> None:
    cls = CoATransferase()
    reaction = _reaction(
        left=[
            {"smiles": "CCC(=O)[O-]", "name": "propanoate"},
            {
                "smiles": "CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-]",
                "name": "acetyl-CoA",
            },
        ],
        right=[
            {
                "smiles": "CCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-]",
                "name": "propanoyl-CoA",
            },
            {"smiles": "CC(=O)[O-]", "name": "acetate"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "coa" in result.explanation.lower()


def test_ec_283_rejects_simple_acyl_transfer_without_new_thioester() -> None:
    cls = CoATransferase()
    reaction = _reaction(
        left=[
            {
                "smiles": "CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-]",
                "name": "acetyl-CoA",
            },
            {"smiles": "OCC[N+](C)(C)C", "name": "choline"},
        ],
        right=[
            {"smiles": "CC(=O)OCC[N+](C)(C)C", "name": "acetylcholine"},
            {
                "smiles": "SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-]",
                "name": "CoA",
            },
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_316_positive_sulfuric_ester_hydrolase() -> None:
    cls = SulfuricEsterHydrolase()
    reaction = _reaction(
        left=[
            {"smiles": "C[C@@H](OS(=O)(=O)[O-])C(=O)[O-]", "name": "(R)-2-O-sulfolactate"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "C[C@@H](O)C(=O)[O-]", "name": "(R)-lactate"},
            {"smiles": "O=S(=O)([O-])[O-]", "name": "sulfate"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "sulfuric ester" in result.explanation.lower()


def test_ec_316_rejects_sulfamate_hydrolysis() -> None:
    cls = SulfuricEsterHydrolase()
    reaction = _reaction(
        left=[
            {"smiles": "NS(=O)(=O)OC1CCCCC1", "name": "cyclohexylsulfamate"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "NC1CCCCC1", "name": "cyclohexylamine"},
            {"smiles": "O=S(=O)([O-])[O-]", "name": "sulfate"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_274_positive_nucleoside_monophosphate_kinase() -> None:
    cls = PhosphotransferasePhosphateGroupAsAcceptor()
    reaction = _reaction(
        left=[
            {
                "smiles": "NC1=NC2=C(N=CN2[C@@H]2O[C@H](COP(=O)([O-])[O-])[C@@H](O)[C@H]2O)C(=O)N1",
                "name": "GMP",
            },
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
        ],
        right=[
            {
                "smiles": "NC1=NC2=C(N=CN2[C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]2O)C(=O)N1",
                "name": "GDP",
            },
            {"chebi_id": CHEBI_ADP, "name": "ADP"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "phosphate acceptor" in result.explanation.lower()


def test_ec_274_positive_bisphosphate_phosphokinase() -> None:
    cls = PhosphotransferasePhosphateGroupAsAcceptor()
    reaction = _reaction(
        left=[
            {
                "smiles": "O=P([O-])([O-])OC[C@H]1O[C@H](OP(=O)([O-])[O-])[C@H](O)[C@@H]1O",
                "name": "alpha-D-ribose 1,5-bisphosphate",
            },
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
        ],
        right=[
            {
                "smiles": "O=P([O-])([O-])OC[C@H]1O[C@H](OP(=O)([O-])OP(=O)([O-])[O-])[C@H](O)[C@@H]1O",
                "name": "phosphoribosyl diphosphate",
            },
            {"chebi_id": CHEBI_ADP, "name": "ADP"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_274_positive_amp_kinase_branch() -> None:
    cls = PhosphotransferasePhosphateGroupAsAcceptor()
    reaction = _reaction(
        left=[
            {
                "chebi_id": CHEBI_AMP,
                "smiles": "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(=O)([O-])[O-])[C@@H](O)[C@H]1O",
                "name": "AMP",
            },
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
        ],
        right=[
            {
                "chebi_id": CHEBI_ADP,
                "smiles": "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]1O",
                "name": "ADP",
                "count": 2,
            }
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_274_rejects_simple_kinase() -> None:
    cls = PhosphotransferasePhosphateGroupAsAcceptor()
    reaction = _reaction(
        left=[
            {"smiles": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O", "name": "glucose"},
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
        ],
        right=[
            {
                "smiles": "O=P([O-])([O-])OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
                "name": "glucose 6-phosphate",
            },
            {"chebi_id": CHEBI_ADP, "name": "ADP"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_153_positive_amine_oxidase() -> None:
    cls = OxidoreductaseActingOnTheCHNHGroupOfDonorsOxygenAsAcceptor()
    reaction = _reaction(
        left=[
            {"smiles": "C[NH2+]CCCC[C@H]([NH3+])C(=O)[O-]", "name": "N(6)-methyl-L-lysine"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "[NH3+]CCCC[C@H]([NH3+])C(=O)[O-]", "name": "L-lysine"},
            {"smiles": "C=O", "name": "formaldehyde"},
            {"chebi_id": CHEBI_H2O2, "smiles": "OO", "name": "hydrogen peroxide"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "oxygen acceptor" in result.explanation.lower()


def test_ec_153_positive_polyamine_oxidase() -> None:
    cls = OxidoreductaseActingOnTheCHNHGroupOfDonorsOxygenAsAcceptor()
    reaction = _reaction(
        left=[
            {"smiles": "[NH3+]CCCC[NH2+]CCC[NH3+]", "name": "spermidine"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "[H]C(=O)CC[NH3+]", "name": "3-aminopropanal"},
            {"smiles": "[NH3+]CCCC[NH3+]", "name": "putrescine"},
            {"chebi_id": CHEBI_H2O2, "smiles": "OO", "name": "hydrogen peroxide"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_153_rejects_ch_oh_oxidase() -> None:
    cls = OxidoreductaseActingOnTheCHNHGroupOfDonorsOxygenAsAcceptor()
    reaction = _reaction(
        left=[
            {"smiles": "CCO", "name": "ethanol"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
        ],
        right=[
            {"smiles": "CC=O", "name": "acetaldehyde"},
            {"chebi_id": CHEBI_H2O2, "smiles": "OO", "name": "hydrogen peroxide"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_153_accepts_hydroxynicotine_oxidative_ring_opening() -> None:
    cls = OxidoreductaseActingOnTheCHNHGroupOfDonorsOxygenAsAcceptor()
    reaction = _reaction(
        left=[
            {
                "smiles": "[H][C@]1(C2=CN=C(O)C=C2)CCC[NH+]1C",
                "name": "(R)-6-hydroxynicotine",
            },
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {
                "smiles": "C[NH2+]CCCC(=O)C1=CN=C(O)C=C1",
                "name": "6-hydroxypseudooxynicotine",
            },
            {"chebi_id": CHEBI_H2O2, "smiles": "OO", "name": "hydrogen peroxide"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_153_rejects_sarcosine_oxidative_deamination() -> None:
    cls = OxidoreductaseActingOnTheCHNHGroupOfDonorsOxygenAsAcceptor()
    reaction = _reaction(
        left=[
            {"smiles": "C[NH2+]CC(=O)[O-]", "name": "sarcosine"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "C[NH3+]", "name": "methylamine"},
            {"smiles": "[H]C(=O)C(=O)[O-]", "name": "glyoxylate"},
            {"chebi_id": CHEBI_H2O2, "smiles": "OO", "name": "hydrogen peroxide"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_353_positive_arginine_deiminase() -> None:
    cls = HydrolaseActingOnCarbonNitrogenButNotPeptideBondsInLinearAmidines()
    reaction = _reaction(
        left=[
            {"smiles": "NC(=[NH2+])NCCC[C@H]([NH3+])C(=O)[O-]", "name": "L-arginine"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "NC(=O)NCCC[C@H]([NH3+])C(=O)[O-]", "name": "L-citrulline"},
            {"smiles": "[NH4+]", "name": "ammonium"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "amidine" in result.explanation.lower()


def test_ec_353_positive_guanidinoacetase() -> None:
    cls = HydrolaseActingOnCarbonNitrogenButNotPeptideBondsInLinearAmidines()
    reaction = _reaction(
        left=[
            {"smiles": "NC(=[NH2+])NCC(=O)[O-]", "name": "guanidinoacetate"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "NC(N)=O", "name": "urea"},
            {"smiles": "[NH3+]CC(=O)[O-]", "name": "glycine"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_353_rejects_linear_amide_hydrolysis() -> None:
    cls = HydrolaseActingOnCarbonNitrogenButNotPeptideBondsInLinearAmidines()
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

    assert result.is_member is False


def test_ec_353_positive_formimidoyl_hydrolysis() -> None:
    cls = HydrolaseActingOnCarbonNitrogenButNotPeptideBondsInLinearAmidines()
    reaction = _reaction(
        left=[
            {"smiles": "[H]C(=[NH2+])N[C@@H](CCC(=O)[O-])C(=O)[O-]", "name": "N-formimidoyl-L-glutamate"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "[H]C(N)=O", "name": "formamide"},
            {"smiles": "[NH3+][C@@H](CCC(=O)[O-])C(=O)[O-]", "name": "L-glutamate"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_353_positive_succinyl_arginine_hydrolase() -> None:
    cls = HydrolaseActingOnCarbonNitrogenButNotPeptideBondsInLinearAmidines()
    reaction = _reaction(
        left=[
            {
                "smiles": "NC(=[NH2+])NCCC[C@H](NC(=O)CCC(=O)[O-])C(=O)[O-]",
                "name": "N(2)-succinyl-L-arginine",
            },
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water", "count": 2},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron", "count": 2},
        ],
        right=[
            {
                "smiles": "[NH3+]CCC[C@H](NC(=O)CCC(=O)[O-])C(=O)[O-]",
                "name": "N(2)-succinyl-L-ornithine",
            },
            {"chebi_id": CHEBI_NH4, "smiles": "[NH4+]", "name": "ammonium", "count": 2},
            {"smiles": "O=C=O", "name": "carbon dioxide"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_143_positive_amino_acid_oxidase() -> None:
    cls = OxidoreductaseActingOnTheCHNH2GroupOfDonorsOxygenAsAcceptor()
    reaction = _reaction(
        left=[
            {"smiles": "[NH3+]CC(=O)[O-]", "name": "glycine"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "[H]C(=O)C(=O)[O-]", "name": "glyoxylate"},
            {"chebi_id": CHEBI_H2O2, "smiles": "OO", "name": "hydrogen peroxide"},
            {"chebi_id": CHEBI_NH4, "smiles": "[NH4+]", "name": "ammonium"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "ch-nh2" in result.explanation.lower()


def test_ec_143_positive_imine_forming_amino_acid_oxidase() -> None:
    cls = OxidoreductaseActingOnTheCHNH2GroupOfDonorsOxygenAsAcceptor()
    reaction = _reaction(
        left=[
            {"smiles": "[NH3+][C@@H](CC(=O)[O-])C(=O)[O-]", "name": "L-aspartate"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
        ],
        right=[
            {"smiles": "[NH2+]=C(CC(=O)[O-])C(=O)[O-]", "name": "iminoaspartate"},
            {"chebi_id": CHEBI_H2O2, "smiles": "OO", "name": "hydrogen peroxide"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_143_positive_sarcosine_oxidase_branch() -> None:
    cls = OxidoreductaseActingOnTheCHNH2GroupOfDonorsOxygenAsAcceptor()
    reaction = _reaction(
        left=[
            {"smiles": "C[NH2+]CC(=O)[O-]", "name": "sarcosine"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "C[NH3+]", "name": "methylamine"},
            {"smiles": "[H]C(=O)C(=O)[O-]", "name": "glyoxylate"},
            {"chebi_id": CHEBI_H2O2, "smiles": "OO", "name": "hydrogen peroxide"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_143_rejects_ch_nh_oxidase_boundary() -> None:
    cls = OxidoreductaseActingOnTheCHNH2GroupOfDonorsOxygenAsAcceptor()
    reaction = _reaction(
        left=[
            {"smiles": "C[NH2+]CCCC[C@H]([NH3+])C(=O)[O-]", "name": "N(6)-methyl-L-lysine"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "[NH3+]CCCC[C@H]([NH3+])C(=O)[O-]", "name": "L-lysine"},
            {"smiles": "C=O", "name": "formaldehyde"},
            {"chebi_id": CHEBI_H2O2, "smiles": "OO", "name": "hydrogen peroxide"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_511_positive_amino_acid_racemase() -> None:
    cls = RacemaseAndEpimeraseActingOnAminoAcidsAndDerivatives()
    reaction = _reaction(
        left=[
            {"smiles": "[NH3+][C@@H](CO)C(=O)[O-]", "name": "L-serine"},
        ],
        right=[
            {"smiles": "[NH3+][C@H](CO)C(=O)[O-]", "name": "D-serine"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "stereo" in result.explanation.lower()


def test_ec_511_positive_peptidyl_amino_acid_racemase() -> None:
    cls = RacemaseAndEpimeraseActingOnAminoAcidsAndDerivatives()
    reaction = _reaction(
        left=[
            {"smiles": "*N[C@@H](CO)C(*)=O", "name": "L-seryl-[protein]"},
        ],
        right=[
            {"smiles": "*N[C@H](CO)C(*)=O", "name": "D-seryl-[protein]"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_511_rejects_sugar_epimerase() -> None:
    cls = RacemaseAndEpimeraseActingOnAminoAcidsAndDerivatives()
    reaction = _reaction(
        left=[
            {
                "smiles": "O=P([O-])([O-])OP(=O)([O-])OCC1OC(n2ccc(=O)[nH]c2=O)C(O)C1O",
                "name": "UDP-glucose",
            },
        ],
        right=[
            {
                "smiles": "O=P([O-])([O-])OP(=O)([O-])OCC1OC(n2ccc(=O)[nH]c2=O)C(O)C1O",
                "name": "UDP-galactose",
            },
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_181_positive_glutathione_reductase_branch() -> None:
    cls = OxidoreductaseActingOnASulfurGroupOfDonorsNADPAsAcceptor()
    reaction = _reaction(
        left=[
            {"smiles": "[NH3+][C@@H](CCC(=O)N[C@@H](CS)C(=O)NCC(=O)[O-])C(=O)[O-]", "name": "glutathione", "count": 2},
            {"chebi_id": CHEBI_NADP_PLUS, "smiles": "NC(=O)C1=CC=C[N+]([C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])OC[C@H]3O[C@@H](N4C=NC5=C4N=CN=C5N)[C@H](OP(=O)([O-])[O-])[C@@H]3O)[C@@H](O)[C@H]2O)=C1", "name": "NADP+"},
        ],
        right=[
            {"smiles": "[NH3+][C@@H](CCC(=O)N[C@@H](CSSC[C@H](NC(=O)CC[C@H]([NH3+])C(=O)[O-])C(=O)NCC(=O)[O-])C(=O)NCC(=O)[O-])C(=O)[O-]", "name": "glutathione disulfide"},
            {"chebi_id": CHEBI_NADPH, "smiles": "NC(=O)C1=CN([C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])OC[C@H]3O[C@@H](N4C=NC5=C4N=CN=C5N)[C@H](OP(=O)([O-])[O-])[C@@H]3O)[C@@H](O)[C@H]2O)C=CC1", "name": "NADPH"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "sulfur" in result.explanation.lower()


def test_ec_181_positive_sulfide_oxidation() -> None:
    cls = OxidoreductaseActingOnASulfurGroupOfDonorsNADPAsAcceptor()
    reaction = _reaction(
        left=[
            {"smiles": "[H][S-]", "name": "hydrogen sulfide"},
            {"chebi_id": CHEBI_NADP_PLUS, "smiles": "NC(=O)C1=CC=C[N+]([C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])OC[C@H]3O[C@@H](N4C=NC5=C4N=CN=C5N)[C@H](OP(=O)([O-])[O-])[C@@H]3O)[C@@H](O)[C@H]2O)=C1", "name": "NADP+", "count": 3},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water", "count": 3},
        ],
        right=[
            {"smiles": "O=S([O-])[O-]", "name": "sulfite"},
            {"chebi_id": CHEBI_NADPH, "smiles": "NC(=O)C1=CN([C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])OC[C@H]3O[C@@H](N4C=NC5=C4N=CN=C5N)[C@H](OP(=O)([O-])[O-])[C@@H]3O)[C@@H](O)[C@H]2O)C=CC1", "name": "NADPH", "count": 3},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron", "count": 4},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_181_rejects_nad_dependent_alcohol_dehydrogenase() -> None:
    cls = OxidoreductaseActingOnASulfurGroupOfDonorsNADPAsAcceptor()
    reaction = _reaction(
        left=[
            {"smiles": "CCO", "name": "ethanol"},
            {"chebi_id": CHEBI_NAD_PLUS, "smiles": "NC(=O)C1=CC=C[N+]([C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])OC[C@H]3O[C@@H](N4C=NC5=C4N=CN=C5N)[C@H](O)[C@@H]3O)[C@@H](O)[C@H]2O)=C1", "name": "NAD+"},
        ],
        right=[
            {"smiles": "CC=O", "name": "acetaldehyde"},
            {"chebi_id": CHEBI_NADH, "smiles": "NC(=O)C1=CN([C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])OC[C@H]3O[C@@H](N4C=NC5=C4N=CN=C5N)[C@H](O)[C@@H]3O)[C@@H](O)[C@H]2O)C=CC1", "name": "NADH"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_281_positive_thiocyanate_forming_sulfurtransferase() -> None:
    cls = Sulfurtransferase()
    reaction = _reaction(
        left=[
            {"smiles": "[H]SS(=O)(=O)[O-]", "name": "thiosulfate"},
            {"smiles": "C#N", "name": "hydrogen cyanide"},
        ],
        right=[
            {"smiles": "N#C[S-]", "name": "thiocyanate"},
            {"smiles": "O=S([O-])[O-]", "name": "sulfite"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron", "count": 2},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "sulfur" in result.explanation.lower()


def test_ec_281_positive_mercaptopyruvate_sulfurtransferase() -> None:
    cls = Sulfurtransferase()
    reaction = _reaction(
        left=[
            {"smiles": "O=C([O-])C(=O)CS", "name": "3-mercaptopyruvate"},
            {"smiles": "*N[C@@H](CS)C(*)=O", "name": "[thioredoxin]-dithiol", "count": 2},
        ],
        right=[
            {"smiles": "*N[C@@H](CSSC[C@H](N*)C(*)=O)C(*)=O", "name": "[thioredoxin]-disulfide"},
            {"smiles": "[H][S-]", "name": "hydrogen sulfide"},
            {"smiles": "CC(=O)C(=O)[O-]", "name": "pyruvate"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_281_rejects_sulfur_group_nadp_oxidoreductase() -> None:
    cls = Sulfurtransferase()
    reaction = _reaction(
        left=[
            {"smiles": "[NH3+][C@@H](CCC(=O)N[C@@H](CS)C(=O)NCC(=O)[O-])C(=O)[O-]", "name": "glutathione", "count": 2},
            {"chebi_id": CHEBI_NADP_PLUS, "smiles": "NC(=O)C1=CC=C[N+]([C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])OC[C@H]3O[C@@H](N4C=NC5=C4N=CN=C5N)[C@H](OP(=O)([O-])[O-])[C@@H]3O)[C@@H](O)[C@H]2O)=C1", "name": "NADP+"},
        ],
        right=[
            {"smiles": "[NH3+][C@@H](CCC(=O)N[C@@H](CSSC[C@H](NC(=O)CC[C@H]([NH3+])C(=O)[O-])C(=O)NCC(=O)[O-])C(=O)NCC(=O)[O-])C(=O)[O-]", "name": "glutathione disulfide"},
            {"chebi_id": CHEBI_NADPH, "smiles": "NC(=O)C1=CN([C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])OC[C@H]3O[C@@H](N4C=NC5=C4N=CN=C5N)[C@H](OP(=O)([O-])[O-])[C@@H]3O)[C@@H](O)[C@H]2O)C=CC1", "name": "NADPH"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_631_positive_ammonia_ligase() -> None:
    cls = AcidAmmoniaOrAmideLigase()
    reaction = _reaction(
        left=[
            {"smiles": "[NH3+][C@@H](CC(=O)[O-])C(=O)[O-]", "name": "L-aspartate"},
            {"chebi_id": CHEBI_NH4, "smiles": "[H][N+]([H])([H])[H]", "name": "ammonium"},
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
        ],
        right=[
            {"smiles": "NC(=O)C[C@H]([NH3+])C(=O)[O-]", "name": "L-asparagine"},
            {"chebi_id": CHEBI_AMP, "name": "AMP"},
            {"chebi_id": "CHEBI:33019", "name": "diphosphate"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "acid-ammonia" in result.explanation.lower()


def test_ec_631_positive_amine_acceptor_ligase() -> None:
    cls = AcidAmmoniaOrAmideLigase()
    reaction = _reaction(
        left=[
            {"smiles": "[NH3+]CCCC[NH3+]", "name": "putrescine"},
            {
                "smiles": "[NH3+][C@@H](CCC(=O)[O-])C(=O)[O-]",
                "name": "L-glutamate",
            },
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
        ],
        right=[
            {
                "smiles": "[NH3+]CCCCNC(=O)CC[C@H]([NH3+])C(=O)[O-]",
                "name": "gamma-L-glutamylputrescine",
            },
            {"chebi_id": CHEBI_ADP, "name": "ADP"},
            {"chebi_id": CHEBI_PHOSPHATE, "name": "phosphate"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_631_positive_lipoate_protein_ligase_cache_backed() -> None:
    cls = AcidAmmoniaOrAmideLigase()
    reaction = _reaction(
        left=[
            {"chebi_id": "CHEBI:29969", "name": "L-lysyl-[lipoyl-carrier protein]"},
            {"chebi_id": "CHEBI:83088", "name": "(R)-lipoate"},
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
        ],
        right=[
            {"chebi_id": "CHEBI:83099", "name": "N(6)-[(R)-lipoyl]-L-lysyl-[lipoyl-carrier protein]"},
            {"chebi_id": CHEBI_AMP, "name": "AMP"},
            {"chebi_id": "CHEBI:33019", "name": "diphosphate"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_631_rejects_acid_amino_acid_ligase() -> None:
    cls = AcidAmmoniaOrAmideLigase()
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

    assert result.is_member is False


def test_ec_631_rejects_bicarbonate_carbamoylation() -> None:
    cls = AcidAmmoniaOrAmideLigase()
    reaction = _reaction(
        left=[
            {"smiles": "O=C([O-])O", "name": "hydrogencarbonate"},
            {"chebi_id": CHEBI_NH4, "smiles": "[H][N+]([H])([H])[H]", "name": "ammonium"},
            {"chebi_id": CHEBI_ATP, "name": "ATP", "count": 2},
        ],
        right=[
            {"smiles": "NC(=O)OP(=O)([O-])[O-]", "name": "carbamoyl phosphate"},
            {"chebi_id": CHEBI_ADP, "name": "ADP", "count": 2},
            {"chebi_id": CHEBI_PHOSPHATE, "name": "phosphate"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron", "count": 2},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False
