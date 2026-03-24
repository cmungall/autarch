"""Focused tests for batch 17 specific classifiers."""

from __future__ import annotations

from typing_extensions import TypedDict

from autarch.datamodel import Participant, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_AMP,
    CHEBI_ATP,
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_H_PLUS,
    CHEBI_O2,
    CHEBI_SAH,
    CHEBI_SAM,
)
from autarch.ontology.adp_dependent_nadh_or_nadph_hydrate_dehydratase import (
    ADPDependentNADHOrNADPHHydrateDehydratase,
)
from autarch.ontology.atp_dependent_nadh_or_nadph_hydrate_dehydratase import (
    ATPDependentNADHOrNADPHHydrateDehydratase,
)
from autarch.ontology.caffeine_synthase import CaffeineSynthase
from autarch.ontology.phosphoethanolamine_n_methyltransferase import (
    PhosphoethanolamineNMethyltransferase,
)
from autarch.ontology.polyamine_oxidase_propane_1_3_diamine_forming import (
    PolyamineOxidasePropane13DiamineForming,
)

CHEBI_7_METHYLXANTHINE = "CHEBI:48991"
CHEBI_THEOBROMINE = "CHEBI:28946"
CHEBI_17_DIMETHYLXANTHINE = "CHEBI:25858"
CHEBI_CAFFEINE = "CHEBI:27732"

CHEBI_PHOSPHOETHANOLAMINE = "CHEBI:58190"
CHEBI_N_METHYL_PHOSPHOETHANOLAMINE = "CHEBI:57781"
CHEBI_NN_DIMETHYL_PHOSPHOETHANOLAMINE = "CHEBI:58641"
CHEBI_PHOSPHOCHOLINE = "CHEBI:295975"

CHEBI_SPERMIDINE = "CHEBI:57834"
CHEBI_SPERMINE = "CHEBI:45725"
CHEBI_N1_ACETYLSPERMINE = "CHEBI:58101"
CHEBI_PROPANE_13_DIAMINE = "CHEBI:57484"
CHEBI_4_AMINOBUTANAL = "CHEBI:58264"
CHEBI_3_AMINOPROPYL_4_AMINOBUTANAL = "CHEBI:58869"
CHEBI_3_ACETAMIDOPROPYL_4_AMINOBUTANAL = "CHEBI:58858"

CHEBI_NADHX = "CHEBI:64074"
CHEBI_NADPHX = "CHEBI:64076"
CHEBI_NADH = "CHEBI:57945"
CHEBI_NADPH = "CHEBI:57783"
CHEBI_POLYPHOSPHATE = "CHEBI:16838"


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


def test_go_0102741_positive_caffeine_synthase() -> None:
    cls = CaffeineSynthase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_THEOBROMINE, "name": "theobromine"},
            {"chebi_id": CHEBI_SAM, "name": "SAM"},
        ],
        right=[
            {"chebi_id": CHEBI_CAFFEINE, "name": "caffeine"},
            {"chebi_id": CHEBI_SAH, "name": "SAH"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)
    assert result.is_member is True
    assert "caffeine synthase" in result.explanation.lower() or "caffeine" in result.explanation.lower()


def test_go_0102741_rejects_non_caffeine_xanthine_branch() -> None:
    cls = CaffeineSynthase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_7_METHYLXANTHINE, "name": "7-methylxanthine"},
            {"chebi_id": CHEBI_SAM, "name": "SAM"},
        ],
        right=[
            {"chebi_id": CHEBI_17_DIMETHYLXANTHINE, "name": "1,7-dimethylxanthine"},
            {"chebi_id": CHEBI_SAH, "name": "SAH"},
        ],
    )

    assert cls.check_membership_impl(reaction).is_member is False


def test_go_0000234_positive_phosphoethanolamine_n_methyltransferase() -> None:
    cls = PhosphoethanolamineNMethyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_PHOSPHOETHANOLAMINE, "name": "phosphoethanolamine"},
            {"chebi_id": CHEBI_SAM, "name": "SAM"},
        ],
        right=[
            {"chebi_id": CHEBI_N_METHYL_PHOSPHOETHANOLAMINE, "name": "N-methylethanolamine phosphate"},
            {"chebi_id": CHEBI_SAH, "name": "SAH"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)
    assert result.is_member is True
    assert "phosphoethanolamine" in result.explanation.lower()


def test_go_0000234_rejects_wrong_phosphobase_product() -> None:
    cls = PhosphoethanolamineNMethyltransferase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_PHOSPHOETHANOLAMINE, "name": "phosphoethanolamine"},
            {"chebi_id": CHEBI_SAM, "name": "SAM"},
        ],
        right=[
            {"chebi_id": CHEBI_PHOSPHOCHOLINE, "name": "phosphocholine"},
            {"chebi_id": CHEBI_SAH, "name": "SAH"},
        ],
    )

    assert cls.check_membership_impl(reaction).is_member is False


def test_go_0052900_positive_spermine_oxidase_propane_1_3_diamine_forming() -> None:
    cls = PolyamineOxidasePropane13DiamineForming()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_SPERMINE, "name": "spermine"},
            {"chebi_id": CHEBI_O2, "name": "dioxygen"},
            {"chebi_id": CHEBI_H2O, "name": "water"},
        ],
        right=[
            {"chebi_id": CHEBI_3_AMINOPROPYL_4_AMINOBUTANAL, "name": "N-(3-aminopropyl)-4-aminobutanal"},
            {"chebi_id": CHEBI_PROPANE_13_DIAMINE, "name": "propane-1,3-diamine"},
            {"chebi_id": CHEBI_H2O2, "name": "hydrogen peroxide"},
        ],
    )

    result = cls.check_membership_impl(reaction)
    assert result.is_member is True
    assert "propane-1,3-diamine" in result.explanation.lower() or "spermine oxidase" in result.explanation.lower()


def test_go_0052900_rejects_branch_without_propane_1_3_diamine() -> None:
    cls = PolyamineOxidasePropane13DiamineForming()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_SPERMIDINE, "name": "spermidine"},
            {"chebi_id": CHEBI_O2, "name": "dioxygen"},
            {"chebi_id": CHEBI_H2O, "name": "water"},
        ],
        right=[
            {"chebi_id": CHEBI_4_AMINOBUTANAL, "name": "4-aminobutanal"},
            {"chebi_id": CHEBI_SPERMIDINE, "name": "spermidine"},
            {"chebi_id": CHEBI_H2O2, "name": "hydrogen peroxide"},
        ],
    )

    result = cls.check_membership_impl(reaction)
    assert result.is_member is False


def test_go_0047453_positive_atp_dependent_nadh_or_nadph_hydrate_dehydratase() -> None:
    cls = ATPDependentNADHOrNADPHHydrateDehydratase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_NADPHX, "name": "NADPHX"},
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
        ],
        right=[
            {"chebi_id": CHEBI_ADP, "name": "ADP"},
            {"chebi_id": CHEBI_POLYPHOSPHATE, "name": "phosphate"},
            {"chebi_id": CHEBI_NADPH, "name": "NADPH"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)
    assert result.is_member is True
    assert "nadhx" in result.explanation.lower() or "dehydratase" in result.explanation.lower()


def test_go_0047453_rejects_adp_dependent_branch() -> None:
    cls = ATPDependentNADHOrNADPHHydrateDehydratase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_NADHX, "name": "NADHX"},
            {"chebi_id": CHEBI_ADP, "name": "ADP"},
        ],
        right=[
            {"chebi_id": CHEBI_AMP, "name": "AMP"},
            {"chebi_id": CHEBI_POLYPHOSPHATE, "name": "phosphate"},
            {"chebi_id": CHEBI_NADH, "name": "NADH"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
    )

    assert cls.check_membership_impl(reaction).is_member is False


def test_go_0052855_positive_adp_dependent_nadh_or_nadph_hydrate_dehydratase() -> None:
    cls = ADPDependentNADHOrNADPHHydrateDehydratase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_NADHX, "name": "NADHX"},
            {"chebi_id": CHEBI_ADP, "name": "ADP"},
        ],
        right=[
            {"chebi_id": CHEBI_AMP, "name": "AMP"},
            {"chebi_id": CHEBI_POLYPHOSPHATE, "name": "phosphate"},
            {"chebi_id": CHEBI_NADH, "name": "NADH"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)
    assert result.is_member is True
    assert "dehydratase" in result.explanation.lower()


def test_go_0052855_rejects_atp_dependent_branch() -> None:
    cls = ADPDependentNADHOrNADPHHydrateDehydratase()
    reaction = _reaction(
        left=[
            {"chebi_id": CHEBI_NADPHX, "name": "NADPHX"},
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
        ],
        right=[
            {"chebi_id": CHEBI_ADP, "name": "ADP"},
            {"chebi_id": CHEBI_POLYPHOSPHATE, "name": "phosphate"},
            {"chebi_id": CHEBI_NADPH, "name": "NADPH"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
    )

    assert cls.check_membership_impl(reaction).is_member is False
