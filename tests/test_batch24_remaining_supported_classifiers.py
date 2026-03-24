"""Focused tests for the remaining support-qualified exact classifiers."""

from __future__ import annotations

from autarch.datamodel import Participant, Reaction
from autarch.molecules import CHEBI_H_PLUS, CHEBI_NADH, CHEBI_NADPH, CHEBI_NADP_PLUS, CHEBI_NAD_PLUS
from autarch.ontology.dihydrolipoyl_dehydrogenase_nadh import DihydrolipoylDehydrogenaseNADH
from autarch.ontology.gamma_humulene_synthase import GammaHumuleneSynthase
from autarch.ontology.trna_dihydrouridine_1617_synthase_nad_p import TRNADihydrouridine1617SynthaseNADP
from autarch.ontology.trna_dihydrouridine_20a20b_synthase_nad_p import TRNADihydrouridine20A20BSynthaseNADP

CHEBI_FARNESYL_DIPHOSPHATE = "CHEBI:175763"
CHEBI_DIPHOSPHATE = "CHEBI:33019"
CHEBI_GAMMA_HUMULENE = "CHEBI:49290"
CHEBI_TRNA_DIHYDROURIDINE = "CHEBI:74443"
CHEBI_TRNA_URIDINE = "CHEBI:65315"
CHEBI_DIHYDROLIPOYL_PROTEIN = "CHEBI:83100"
CHEBI_LIPOYL_PROTEIN = "CHEBI:83099"


def test_ec_42356_positive_gamma_humulene_synthase() -> None:
    cls = GammaHumuleneSynthase()
    reaction = Reaction(
        left_participants=[Participant(chebi_id=CHEBI_FARNESYL_DIPHOSPHATE)],
        right_participants=[
            Participant(chebi_id=CHEBI_GAMMA_HUMULENE),
            Participant(chebi_id=CHEBI_DIPHOSPHATE),
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_42356_rejects_unsupported_sesquiterpene_branch() -> None:
    cls = GammaHumuleneSynthase()
    reaction = Reaction(
        left_participants=[Participant(chebi_id=CHEBI_FARNESYL_DIPHOSPHATE)],
        right_participants=[
            Participant(chebi_id="CHEBI:3992"),
            Participant(chebi_id=CHEBI_DIPHOSPHATE),
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_ec_13188_positive_trna_dihydrouridine_1617_synthase_nadp() -> None:
    cls = TRNADihydrouridine1617SynthaseNADP()
    reaction = Reaction(
        left_participants=[
            Participant(chebi_id=CHEBI_TRNA_DIHYDROURIDINE),
            Participant(chebi_id=CHEBI_NADP_PLUS),
        ],
        right_participants=[
            Participant(chebi_id=CHEBI_TRNA_URIDINE),
            Participant(chebi_id=CHEBI_NADPH),
            Participant(chebi_id=CHEBI_H_PLUS),
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_ec_13190_positive_trna_dihydrouridine_20a20b_synthase_nad() -> None:
    cls = TRNADihydrouridine20A20BSynthaseNADP()
    reaction = Reaction(
        left_participants=[
            Participant(chebi_id=CHEBI_TRNA_DIHYDROURIDINE),
            Participant(chebi_id=CHEBI_NAD_PLUS),
        ],
        right_participants=[
            Participant(chebi_id=CHEBI_TRNA_URIDINE),
            Participant(chebi_id=CHEBI_NADH),
            Participant(chebi_id=CHEBI_H_PLUS),
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True


def test_go_0004148_positive_dihydrolipoyl_dehydrogenase_nadh() -> None:
    cls = DihydrolipoylDehydrogenaseNADH()
    reaction = Reaction(
        left_participants=[
            Participant(chebi_id=CHEBI_DIHYDROLIPOYL_PROTEIN),
            Participant(chebi_id=CHEBI_NAD_PLUS),
        ],
        right_participants=[
            Participant(chebi_id=CHEBI_LIPOYL_PROTEIN),
            Participant(chebi_id=CHEBI_NADH),
            Participant(chebi_id=CHEBI_H_PLUS),
        ],
    )
    assert cls.check_membership_impl(reaction).is_member is True
