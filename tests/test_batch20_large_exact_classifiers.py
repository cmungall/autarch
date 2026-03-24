"""Focused tests for large exact-class batch 20."""

from __future__ import annotations

from typing_extensions import TypedDict

from autarch.datamodel import Participant, PolymerType, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS, CHEBI_NADH, CHEBI_NADPH, CHEBI_NADP_PLUS, CHEBI_NAD_PLUS, CHEBI_O2
from autarch.ontology.atp_diphosphatase import ATPDiphosphatase
from autarch.ontology.homomethionine_n_monooxygenase import HomomethionineNMonooxygenase
from autarch.ontology.l_glutamate_gamma_semialdehyde_dehydrogenase import LGlutamateGammaSemialdehydeDehydrogenase
from autarch.ontology.methoxylated_aromatic_compound_corrinoid_protein_co_methyltransferase import MethoxylatedAromaticCompoundCorrinoidProteinCoMethyltransferase
from autarch.ontology.methylxanthine_n1_demethylase import MethylxanthineN1Demethylase
from autarch.ontology.monoterpene_epsilon_lactone_hydrolase import MonoterpeneEpsilonLactoneHydrolase
from autarch.ontology.n_terminal_amino_acid_n_alpha_acetyltransferase_nat_a import NTerminalAminoAcidNAlphaAcetyltransferaseNatA
from autarch.ontology.n_terminal_methionine_n_alpha_acetyltransferase_nat_c import NTerminalMethionineNAlphaAcetyltransferaseNatC
from autarch.ontology.n_terminal_methionine_n_alpha_acetyltransferase_nat_e import NTerminalMethionineNAlphaAcetyltransferaseNatE
from autarch.ontology.nucleotide_diphosphatase import NucleotideDiphosphatase
from autarch.ontology.two_three_cyclic_nucleotide_two_phosphodiesterase import TwoThreeCyclicNucleotideTwoPhosphodiesterase
from autarch.ontology.vitamin_d3_24_hydroxylase import VitaminD324Hydroxylase

CHEBI_ATP = "CHEBI:30616"
CHEBI_AMP = "CHEBI:456215"
CHEBI_DIPHOSPHATE = "CHEBI:33019"
CHEBI_DUTP = "CHEBI:61555"
CHEBI_DUMP = "CHEBI:246422"
CHEBI_COA = "CHEBI:57287"
CHEBI_ACETYL_COA = "CHEBI:57288"
CHEBI_NAT_A_SUBSTRATE = "CHEBI:64739"
CHEBI_NAT_A_PRODUCT = "CHEBI:133375"
CHEBI_NAT_C_SUBSTRATE = "CHEBI:133377"
CHEBI_NAT_C_PRODUCT = "CHEBI:133378"
CHEBI_NAT_E_SUBSTRATE = "CHEBI:133407"
CHEBI_NAT_E_PRODUCT = "CHEBI:133406"
CHEBI_CORRINOID_COI = "CHEBI:85033"
CHEBI_CORRINOID_METHYL = "CHEBI:85035"
CHEBI_GUAIACOL = "CHEBI:28591"
CHEBI_CATECHOL = "CHEBI:18135"
CHEBI_CALCITRIOL = "CHEBI:17823"
CHEBI_CALCITETROL = "CHEBI:47799"
CHEBI_ADRENOXIN_REDUCED = "CHEBI:33738"
CHEBI_ADRENOXIN_OXIDIZED = "CHEBI:33737"
CHEBI_CAFFEINE = "CHEBI:27732"
CHEBI_THEOBROMINE = "CHEBI:28946"
CHEBI_FORMALDEHYDE = "CHEBI:16842"
CHEBI_P5C = "CHEBI:17388"
CHEBI_GLUTAMATE = "CHEBI:29985"
CHEBI_CYCLIC_GMP = "CHEBI:60837"
CHEBI_THREE_GMP = "CHEBI:60732"
CHEBI_LACTONE = "CHEBI:50238"
CHEBI_HYDROXY_ACID = "CHEBI:64224"
CHEBI_HOMOMETHIONINE = "CHEBI:134632"
CHEBI_HOMOMETHIONINE_OXIME = "CHEBI:134682"
CHEBI_REDUCED_HEMOREDUCTASE = "CHEBI:57618"
CHEBI_OXIDIZED_HEMOREDUCTASE = "CHEBI:58210"
CHEBI_CO2 = "CHEBI:16526"


class ParticipantInput(TypedDict, total=False):
    chebi_id: str
    smiles: str
    inchi: str
    name: str
    count: int
    polymer_type: PolymerType


def _reaction(left: list[ParticipantInput], right: list[ParticipantInput]) -> Reaction:
    return Reaction(
        left_participants=[Participant(**participant) for participant in left],
        right_participants=[Participant(**participant) for participant in right],
    )


def test_ec_3619_positive_nucleotide_diphosphatase() -> None:
    cls = NucleotideDiphosphatase()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_DUTP, "name": "dUTP"}, {"chebi_id": CHEBI_H2O, "name": "water"}],
        right=[{"chebi_id": CHEBI_DUMP, "name": "dUMP"}, {"chebi_id": CHEBI_DIPHOSPHATE, "name": "diphosphate"}, {"chebi_id": CHEBI_H_PLUS, "name": "hydron"}],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_3619_rejects_wrong_product_branch() -> None:
    cls = NucleotideDiphosphatase()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_DUTP}, {"chebi_id": CHEBI_H2O}],
        right=[{"chebi_id": CHEBI_AMP}, {"chebi_id": CHEBI_DIPHOSPHATE}, {"chebi_id": CHEBI_H_PLUS}],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_ec_23618_positive_atp_diphosphatase() -> None:
    cls = ATPDiphosphatase()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_ATP}, {"chebi_id": CHEBI_H2O}],
        right=[{"chebi_id": CHEBI_AMP}, {"chebi_id": CHEBI_DIPHOSPHATE}, {"chebi_id": CHEBI_H_PLUS}],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_23618_rejects_non_triphosphate_branch() -> None:
    cls = ATPDiphosphatase()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_AMP}, {"chebi_id": CHEBI_H2O}],
        right=[{"chebi_id": CHEBI_AMP}, {"chebi_id": CHEBI_DIPHOSPHATE}, {"chebi_id": CHEBI_H_PLUS}],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_ec_231255_positive_nat_a() -> None:
    cls = NTerminalAminoAcidNAlphaAcetyltransferaseNatA()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_NAT_A_SUBSTRATE, "polymer_type": PolymerType.PROTEIN}, {"chebi_id": CHEBI_ACETYL_COA}],
        right=[{"chebi_id": CHEBI_NAT_A_PRODUCT, "polymer_type": PolymerType.PROTEIN}, {"chebi_id": CHEBI_COA}, {"chebi_id": CHEBI_H_PLUS}],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_231255_rejects_wrong_protein_pair() -> None:
    cls = NTerminalAminoAcidNAlphaAcetyltransferaseNatA()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_NAT_A_SUBSTRATE, "polymer_type": PolymerType.PROTEIN}, {"chebi_id": CHEBI_ACETYL_COA}],
        right=[{"chebi_id": CHEBI_NAT_C_PRODUCT, "polymer_type": PolymerType.PROTEIN}, {"chebi_id": CHEBI_COA}, {"chebi_id": CHEBI_H_PLUS}],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_ec_231256_positive_nat_c() -> None:
    cls = NTerminalMethionineNAlphaAcetyltransferaseNatC()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_NAT_C_SUBSTRATE, "polymer_type": PolymerType.PROTEIN}, {"chebi_id": CHEBI_ACETYL_COA}],
        right=[{"chebi_id": CHEBI_NAT_C_PRODUCT, "polymer_type": PolymerType.PROTEIN}, {"chebi_id": CHEBI_COA}, {"chebi_id": CHEBI_H_PLUS}],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_231256_accepts_reverse_orientation() -> None:
    cls = NTerminalMethionineNAlphaAcetyltransferaseNatC()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_NAT_C_PRODUCT, "polymer_type": PolymerType.PROTEIN}, {"chebi_id": CHEBI_COA}, {"chebi_id": CHEBI_H_PLUS}],
        right=[{"chebi_id": CHEBI_NAT_C_SUBSTRATE, "polymer_type": PolymerType.PROTEIN}, {"chebi_id": CHEBI_ACETYL_COA}],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_231258_positive_nat_e() -> None:
    cls = NTerminalMethionineNAlphaAcetyltransferaseNatE()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_NAT_E_SUBSTRATE, "polymer_type": PolymerType.PROTEIN}, {"chebi_id": CHEBI_ACETYL_COA}],
        right=[{"chebi_id": CHEBI_NAT_E_PRODUCT, "polymer_type": PolymerType.PROTEIN}, {"chebi_id": CHEBI_COA}, {"chebi_id": CHEBI_H_PLUS}],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_211382_positive_corrinoid_co_methyltransferase() -> None:
    cls = MethoxylatedAromaticCompoundCorrinoidProteinCoMethyltransferase()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_GUAIACOL}, {"chebi_id": CHEBI_CORRINOID_COI, "polymer_type": PolymerType.PROTEIN}, {"chebi_id": CHEBI_H_PLUS}],
        right=[{"chebi_id": CHEBI_CATECHOL}, {"chebi_id": CHEBI_CORRINOID_METHYL, "polymer_type": PolymerType.PROTEIN}],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_211382_accepts_reverse_orientation() -> None:
    cls = MethoxylatedAromaticCompoundCorrinoidProteinCoMethyltransferase()
    reaction = _reaction(
        left=[{"chebi_id": "CHEBI:52678"}, {"chebi_id": CHEBI_CORRINOID_METHYL, "polymer_type": PolymerType.PROTEIN}],
        right=[{"chebi_id": "CHEBI:59114"}, {"chebi_id": CHEBI_CORRINOID_COI, "polymer_type": PolymerType.PROTEIN}, {"chebi_id": CHEBI_H_PLUS}],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_1141516_positive_vitamin_d3_24_hydroxylase() -> None:
    cls = VitaminD324Hydroxylase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_CALCITRIOL},
            {"chebi_id": CHEBI_ADRENOXIN_REDUCED, "count": 2},
            {"chebi_id": CHEBI_O2},
            {"chebi_id": CHEBI_H_PLUS, "count": 2},
        ],
        right=[
            {"chebi_id": CHEBI_CALCITETROL},
            {"chebi_id": CHEBI_ADRENOXIN_OXIDIZED, "count": 2},
            {"chebi_id": CHEBI_H2O},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_11413178_positive_methylxanthine_n1_demethylase() -> None:
    cls = MethylxanthineN1Demethylase()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_CAFFEINE}, {"chebi_id": CHEBI_NADPH}, {"chebi_id": CHEBI_O2}, {"chebi_id": CHEBI_H_PLUS}],
        right=[{"chebi_id": CHEBI_THEOBROMINE}, {"chebi_id": CHEBI_FORMALDEHYDE}, {"chebi_id": CHEBI_NADP_PLUS}, {"chebi_id": CHEBI_H2O}],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_12188_positive_l_glutamate_gamma_semialdehyde_dehydrogenase() -> None:
    cls = LGlutamateGammaSemialdehydeDehydrogenase()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_P5C}, {"chebi_id": CHEBI_NAD_PLUS}, {"chebi_id": CHEBI_H2O, "count": 2}],
        right=[{"chebi_id": CHEBI_GLUTAMATE}, {"chebi_id": CHEBI_NADH}, {"chebi_id": CHEBI_H_PLUS}],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_31416_positive_two_three_cyclic_nucleotide_two_phosphodiesterase() -> None:
    cls = TwoThreeCyclicNucleotideTwoPhosphodiesterase()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_CYCLIC_GMP}, {"chebi_id": CHEBI_H2O}],
        right=[{"chebi_id": CHEBI_THREE_GMP}, {"chebi_id": CHEBI_H_PLUS}],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_31183_positive_monoterpene_epsilon_lactone_hydrolase() -> None:
    cls = MonoterpeneEpsilonLactoneHydrolase()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_LACTONE}, {"chebi_id": CHEBI_H2O}],
        right=[{"chebi_id": CHEBI_HYDROXY_ACID}, {"chebi_id": CHEBI_H_PLUS}],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_1141442_positive_homomethionine_n_monooxygenase() -> None:
    cls = HomomethionineNMonooxygenase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_HOMOMETHIONINE},
            {"chebi_id": CHEBI_REDUCED_HEMOREDUCTASE, "count": 2},
            {"chebi_id": CHEBI_O2, "count": 2},
        ],
        right=[
            {"chebi_id": CHEBI_HOMOMETHIONINE_OXIME},
            {"chebi_id": CHEBI_OXIDIZED_HEMOREDUCTASE, "count": 2},
            {"chebi_id": CHEBI_CO2},
            {"chebi_id": CHEBI_H2O, "count": 3},
            {"chebi_id": CHEBI_H_PLUS, "count": 2},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True
