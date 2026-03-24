"""Focused tests for supported batch 21 classifiers."""

from __future__ import annotations

from typing_extensions import TypedDict

from autarch.datamodel import Participant, PolymerType, Reaction
from autarch.molecules import CHEBI_ADP, CHEBI_ATP, CHEBI_H2O, CHEBI_H_PLUS, CHEBI_NADH, CHEBI_NADPH, CHEBI_NADP_PLUS, CHEBI_NAD_PLUS, CHEBI_O2
from autarch.ontology.deoxynucleoside_kinase import DeoxynucleosideKinase
from autarch.ontology.delta24_sterol_reductase import Delta24SterolReductase
from autarch.ontology.gdp_4_dehydro_d_rhamnose_reductase import GDP4DehydroDRhamnoseReductase
from autarch.ontology.methylxanthine_n3_demethylase import MethylxanthineN3Demethylase
from autarch.ontology.n_terminal_methionine_n_alpha_acetyltransferase_nat_b import NTerminalMethionineNAlphaAcetyltransferaseNatB
from autarch.ontology.peptide_o_fucosyltransferase import PeptideOFucosyltransferase
from autarch.ontology.polyamine_oxidase import PolyamineOxidase
from autarch.ontology.r_amidase import RAmidase
from autarch.ontology.ribonucleoside_diphosphate_reductase_thioredoxin_disulfide_as_acceptor import RibonucleosideDiphosphateReductaseThioredoxinDisulfideAsAcceptor
from autarch.ontology.rrna_small_subunit_pseudouridine_methyltransferase_nep1 import RRNASmallSubunitPseudouridineMethyltransferaseNep1
from autarch.ontology.six_seven_dihydropteridine_reductase import SixSevenDihydropteridineReductase
from autarch.ontology.trna_cytidine_5_methyltransferase import TRNACytidine5Methyltransferase

CHEBI_DGDP = "CHEBI:58595"
CHEBI_GDP = "CHEBI:58189"
CHEBI_THIOREDOXIN_DISULFIDE = "CHEBI:50058"
CHEBI_THIOREDOXIN_DITHIOL = "CHEBI:29950"
CHEBI_N1_ACETYLSPERMINE = "CHEBI:58101"
CHEBI_SPERMINE = "CHEBI:45725"
CHEBI_SPERMIDINE = "CHEBI:57834"
CHEBI_ACETAMIDOPROPANAL = "CHEBI:30322"
CHEBI_AMINOPROPANAL = "CHEBI:133427"
CHEBI_H2O2 = "CHEBI:16240"
CHEBI_GDP_6_DEOXY_TALOSE = "CHEBI:57638"
CHEBI_GDP_4_DEHYDRO_RHAMNOSE = "CHEBI:57964"
CHEBI_THYMIDINE = "CHEBI:17748"
CHEBI_DTMP = "CHEBI:63528"
CHEBI_ADENOSINE = "CHEBI:16335"
CHEBI_LATHOSTEROL = "CHEBI:134072"
CHEBI_CHOLESTA_7_24_DIENOL = "CHEBI:147457"
CHEBI_TETRAHYDROPTERIDINE = "CHEBI:28889"
CHEBI_DIHYDROPTERIDINE = "CHEBI:30156"
CHEBI_SERYL_PROTEIN = "CHEBI:29999"
CHEBI_FUCOSYL_SERYL_PROTEIN = "CHEBI:189632"
CHEBI_GDP_BETA_L_FUCOSE = "CHEBI:57273"
CHEBI_METHIONYL_ASPARTYL_PROTEIN = "CHEBI:133045"
CHEBI_ACETYL_METHIONYL_ASPARTYL_PROTEIN = "CHEBI:133063"
CHEBI_ACETYL_COA = "CHEBI:57288"
CHEBI_COA = "CHEBI:57287"
CHEBI_CYTIDINE_IN_TRNA = "CHEBI:82748"
CHEBI_METHYLCYTIDINE_IN_TRNA = "CHEBI:74483"
CHEBI_PSEUDOURIDINE_IN_RRNA = "CHEBI:65314"
CHEBI_METHYLPSEUDOURIDINE_IN_RRNA = "CHEBI:74890"
CHEBI_SAM = "CHEBI:59789"
CHEBI_SAH = "CHEBI:57856"
CHEBI_R_PIPERAZINE_CARBOXAMIDE = "CHEBI:58916"
CHEBI_R_PIPERAZINE_CARBOXYLATE = "CHEBI:58917"
CHEBI_AMMONIUM = "CHEBI:28938"
CHEBI_THEOBROMINE = "CHEBI:28946"
CHEBI_SEVEN_METHYLXANTHINE = "CHEBI:48991"
CHEBI_FORMALDEHYDE = "CHEBI:16842"


class ParticipantInput(TypedDict, total=False):
    chebi_id: str
    name: str
    count: int
    polymer_type: PolymerType


def _reaction(left: list[ParticipantInput], right: list[ParticipantInput]) -> Reaction:
    return Reaction(
        left_participants=[Participant(**participant) for participant in left],
        right_participants=[Participant(**participant) for participant in right],
    )


def test_ec_11741_positive_ribonucleoside_diphosphate_reductase() -> None:
    cls = RibonucleosideDiphosphateReductaseThioredoxinDisulfideAsAcceptor()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_DGDP},
            {"chebi_id": CHEBI_THIOREDOXIN_DISULFIDE},
            {"chebi_id": CHEBI_H2O},
        ],
        right=[
            {"chebi_id": CHEBI_THIOREDOXIN_DITHIOL, "count": 2},
            {"chebi_id": CHEBI_GDP},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_11741_accepts_generic_cached_thioredoxin_shape() -> None:
    cls = RibonucleosideDiphosphateReductaseThioredoxinDisulfideAsAcceptor()
    reaction = _reaction(
        left=[
            {"chebi_id": "CHEBI:73316"},
            {"chebi_id": CHEBI_THIOREDOXIN_DISULFIDE},
            {"chebi_id": CHEBI_H2O},
        ],
        right=[
            {"chebi_id": "CHEBI:57930"},
            {"chebi_id": "CHEBI:29950"},
            {"chebi_id": "CHEBI:29950"},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_go_0046592_positive_polyamine_oxidase() -> None:
    cls = PolyamineOxidase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_N1_ACETYLSPERMINE},
            {"chebi_id": CHEBI_O2},
            {"chebi_id": CHEBI_H2O},
        ],
        right=[
            {"chebi_id": CHEBI_ACETAMIDOPROPANAL},
            {"chebi_id": CHEBI_SPERMIDINE},
            {"chebi_id": CHEBI_H2O2},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_go_0046592_rejects_wrong_polyamine_branch() -> None:
    cls = PolyamineOxidase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_SPERMINE},
            {"chebi_id": CHEBI_O2},
            {"chebi_id": CHEBI_H2O},
        ],
        right=[
            {"chebi_id": CHEBI_ACETAMIDOPROPANAL},
            {"chebi_id": CHEBI_SPERMIDINE},
            {"chebi_id": CHEBI_H2O2},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_go_0042356_positive_gdp_4_dehydro_d_rhamnose_reductase() -> None:
    cls = GDP4DehydroDRhamnoseReductase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_GDP_6_DEOXY_TALOSE},
            {"chebi_id": CHEBI_NADP_PLUS},
        ],
        right=[
            {"chebi_id": CHEBI_GDP_4_DEHYDRO_RHAMNOSE},
            {"chebi_id": CHEBI_NADPH},
            {"chebi_id": CHEBI_H_PLUS},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_go_0019136_positive_deoxynucleoside_kinase() -> None:
    cls = DeoxynucleosideKinase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_THYMIDINE},
            {"chebi_id": CHEBI_ATP},
        ],
        right=[
            {"chebi_id": CHEBI_DTMP},
            {"chebi_id": CHEBI_ADP},
            {"chebi_id": CHEBI_H_PLUS},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_go_0019136_rejects_ribonucleoside_kinase() -> None:
    cls = DeoxynucleosideKinase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_ADENOSINE},
            {"chebi_id": CHEBI_ATP},
        ],
        right=[
            {"chebi_id": CHEBI_DTMP},
            {"chebi_id": CHEBI_ADP},
            {"chebi_id": CHEBI_H_PLUS},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_go_0050614_positive_delta24_sterol_reductase() -> None:
    cls = Delta24SterolReductase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_LATHOSTEROL},
            {"chebi_id": CHEBI_NADP_PLUS},
        ],
        right=[
            {"chebi_id": CHEBI_CHOLESTA_7_24_DIENOL},
            {"chebi_id": CHEBI_NADPH},
            {"chebi_id": CHEBI_H_PLUS},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_go_0004155_positive_six_seven_dihydropteridine_reductase() -> None:
    cls = SixSevenDihydropteridineReductase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_TETRAHYDROPTERIDINE},
            {"chebi_id": CHEBI_NAD_PLUS},
        ],
        right=[
            {"chebi_id": CHEBI_DIHYDROPTERIDINE},
            {"chebi_id": CHEBI_NADH},
            {"chebi_id": CHEBI_H_PLUS},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_go_0046922_positive_peptide_o_fucosyltransferase() -> None:
    cls = PeptideOFucosyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_SERYL_PROTEIN, "polymer_type": PolymerType.PROTEIN},
            {"chebi_id": CHEBI_GDP_BETA_L_FUCOSE},
        ],
        right=[
            {"chebi_id": CHEBI_FUCOSYL_SERYL_PROTEIN, "polymer_type": PolymerType.PROTEIN},
            {"chebi_id": CHEBI_GDP},
            {"chebi_id": CHEBI_H_PLUS},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_231254_positive_nat_b() -> None:
    cls = NTerminalMethionineNAlphaAcetyltransferaseNatB()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_METHIONYL_ASPARTYL_PROTEIN, "polymer_type": PolymerType.PROTEIN},
            {"chebi_id": CHEBI_ACETYL_COA},
        ],
        right=[
            {"chebi_id": CHEBI_ACETYL_METHIONYL_ASPARTYL_PROTEIN, "polymer_type": PolymerType.PROTEIN},
            {"chebi_id": CHEBI_COA},
            {"chebi_id": CHEBI_H_PLUS},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_231254_accepts_reverse_orientation() -> None:
    cls = NTerminalMethionineNAlphaAcetyltransferaseNatB()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_ACETYL_METHIONYL_ASPARTYL_PROTEIN, "polymer_type": PolymerType.PROTEIN},
            {"chebi_id": CHEBI_COA},
            {"chebi_id": CHEBI_H_PLUS},
        ],
        right=[
            {"chebi_id": CHEBI_METHIONYL_ASPARTYL_PROTEIN, "polymer_type": PolymerType.PROTEIN},
            {"chebi_id": CHEBI_ACETYL_COA},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_go_0016428_positive_trna_cytidine_5_methyltransferase() -> None:
    cls = TRNACytidine5Methyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_CYTIDINE_IN_TRNA},
            {"chebi_id": CHEBI_SAM},
        ],
        right=[
            {"chebi_id": CHEBI_METHYLCYTIDINE_IN_TRNA},
            {"chebi_id": CHEBI_SAH},
            {"chebi_id": CHEBI_H_PLUS},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_211260_positive_nep1() -> None:
    cls = RRNASmallSubunitPseudouridineMethyltransferaseNep1()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_PSEUDOURIDINE_IN_RRNA},
            {"chebi_id": CHEBI_SAM},
        ],
        right=[
            {"chebi_id": CHEBI_METHYLPSEUDOURIDINE_IN_RRNA},
            {"chebi_id": CHEBI_SAH},
            {"chebi_id": CHEBI_H_PLUS},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_351100_positive_r_amidase() -> None:
    cls = RAmidase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_R_PIPERAZINE_CARBOXAMIDE},
            {"chebi_id": CHEBI_H2O},
        ],
        right=[
            {"chebi_id": CHEBI_R_PIPERAZINE_CARBOXYLATE},
            {"chebi_id": CHEBI_AMMONIUM},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_11413179_positive_methylxanthine_n3_demethylase() -> None:
    cls = MethylxanthineN3Demethylase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_THEOBROMINE},
            {"chebi_id": CHEBI_NADH},
            {"chebi_id": CHEBI_O2},
            {"chebi_id": CHEBI_H_PLUS},
        ],
        right=[
            {"chebi_id": CHEBI_SEVEN_METHYLXANTHINE},
            {"chebi_id": CHEBI_FORMALDEHYDE},
            {"chebi_id": CHEBI_NAD_PLUS},
            {"chebi_id": CHEBI_H2O},
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True
