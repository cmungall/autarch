"""Focused tests for the next cache-backed GO support batch."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from typing_extensions import TypedDict

from autarch.datamodel import Participant, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS, CHEBI_O2
from autarch.ontology.amine_n_methyltransferase_activity import AmineNMethyltransferaseActivity
from autarch.ontology.aureusidin_synthase_activity import AureusidinSynthaseActivity
from autarch.ontology.feruloyl_coa_hydratase_lyase_activity import FeruloylCoAHydrataseLyaseActivity
from autarch.ontology.glyceollin_synthase_activity import GlyceollinSynthaseActivity
from autarch.ontology.limonene_12_monooxygenase_nadh_or_nadph_activity import Limonene12MonooxygenaseNADHOrNADPHActivity
from autarch.ontology.methyl_co_iii_methylamine_specific_corrinoid_protein_coenzyme_m_methyltransferase_activity import (
    MethylCoIIIMethylamineSpecificCorrinoidProteinCoenzymeMMethyltransferaseActivity,
)

CHEBI_LIMONENE_S = "CHEBI:15383"
CHEBI_LIMONENE_EPOXIDE = "CHEBI:16431"
CHEBI_METHYL_CORRINOID = "CHEBI:85035"
CHEBI_CORRINOID_CO_I = "CHEBI:85033"
CHEBI_COENZYME_M = "CHEBI:58319"
CHEBI_METHYL_COM = "CHEBI:58286"
CHEBI_AUREUSIDIN_SUBSTRATE = "CHEBI:77622"
CHEBI_AUREUSIDIN_PRODUCT = "CHEBI:66905"
CHEBI_GLYCEOLLIN_SUBSTRATE = "CHEBI:50036"
CHEBI_GLYCEOLLIN_PRODUCT = "CHEBI:16470"
CHEBI_REDUCED_HEMOPROTEIN_REDUCTASE = "CHEBI:57618"
CHEBI_OXIDIZED_HEMOPROTEIN_REDUCTASE = "CHEBI:58210"
CHEBI_CAFFEOYL_COA = "CHEBI:87136"
CHEBI_DIHYDROXYBENZALDEHYDE = "CHEBI:50205"
CHEBI_ACETYL_COA = "CHEBI:57288"
CHEBI_GENERIC_PRIMARY_AMINE = "CHEBI:65296"
CHEBI_GENERIC_METHYLATED_AMINE = "CHEBI:131823"
CHEBI_SAM_ZWITTERION = "CHEBI:59789"
CHEBI_SAH_ZWITTERION = "CHEBI:57856"


class ParticipantInput(TypedDict, total=False):
    chebi_id: str
    smiles: str
    inchi: str
    name: str
    count: int


def _load_rhea_reaction(rhea_id: str) -> Reaction:
    chebi_to_smiles = json.loads(Path("cache/chebi_smiles.json").read_text())
    for line in Path("cache/rhea_reactions.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        term = json.loads(line)
        if term["rhea_id"] == rhea_id:
            reaction = term["reaction"]
            left_participants = []
            right_participants = []
            for participant in reaction["left_participants"]:
                chebi_id = participant.get("chebi_id")
                smiles = participant.get("smiles") or (chebi_to_smiles.get(chebi_id) if chebi_id else None)
                left_participants.append(Participant(
                    chebi_id=chebi_id,
                    smiles=smiles,
                    inchi=participant.get("inchi"),
                    location=participant.get("location"),
                    name=participant.get("name"),
                    polymer_index=participant.get("polymer_index"),
                    polymer_type=participant.get("polymer_type"),
                ))
            for participant in reaction["right_participants"]:
                chebi_id = participant.get("chebi_id")
                smiles = participant.get("smiles") or (chebi_to_smiles.get(chebi_id) if chebi_id else None)
                right_participants.append(Participant(
                    chebi_id=chebi_id,
                    smiles=smiles,
                    inchi=participant.get("inchi"),
                    location=participant.get("location"),
                    name=participant.get("name"),
                    polymer_index=participant.get("polymer_index"),
                    polymer_type=participant.get("polymer_type"),
                ))
            return Reaction(left_participants=left_participants, right_participants=right_participants, label=term.get("label", ""))
    raise AssertionError(f"Missing cached RHEA reaction: {rhea_id}")


def _reaction(left: list[ParticipantInput], right: list[ParticipantInput]) -> Reaction:
    return Reaction(
        left_participants=[Participant(**participant) for participant in left],
        right_participants=[Participant(**participant) for participant in right],
    )


@pytest.mark.parametrize(
    ("classifier", "rhea_id"),
    [
        (Limonene12MonooxygenaseNADHOrNADPHActivity, "RHEA:26085"),
        (MethylCoIIIMethylamineSpecificCorrinoidProteinCoenzymeMMethyltransferaseActivity, "RHEA:18773"),
        (AureusidinSynthaseActivity, "RHEA:34195"),
        (GlyceollinSynthaseActivity, "RHEA:22844"),
        (FeruloylCoAHydrataseLyaseActivity, "RHEA:36307"),
        (AmineNMethyltransferaseActivity, "RHEA:23136"),
    ],
)
def test_batch28_positive_cached_examples(classifier, rhea_id) -> None:
    reaction = _load_rhea_reaction(rhea_id)
    assert classifier().check_membership(reaction).is_member is True


def test_limonene_12_monooxygenase_requires_nicotinamide_pair() -> None:
    cls = Limonene12MonooxygenaseNADHOrNADPHActivity()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_LIMONENE_S}, {"chebi_id": CHEBI_O2}, {"chebi_id": CHEBI_H_PLUS}],
        right=[{"chebi_id": CHEBI_LIMONENE_EPOXIDE}, {"chebi_id": CHEBI_H2O}],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_corrinoid_methyltransferase_requires_methyl_com() -> None:
    cls = MethylCoIIIMethylamineSpecificCorrinoidProteinCoenzymeMMethyltransferaseActivity()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_METHYL_CORRINOID}, {"chebi_id": CHEBI_COENZYME_M}],
        right=[{"chebi_id": CHEBI_CORRINOID_CO_I}, {"chebi_id": CHEBI_H_PLUS}],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_aureusidin_synthase_requires_oxygen() -> None:
    cls = AureusidinSynthaseActivity()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_AUREUSIDIN_SUBSTRATE}],
        right=[{"chebi_id": CHEBI_AUREUSIDIN_PRODUCT}, {"chebi_id": CHEBI_H2O}],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_glyceollin_synthase_requires_reductase() -> None:
    cls = GlyceollinSynthaseActivity()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_GLYCEOLLIN_SUBSTRATE}, {"chebi_id": CHEBI_O2}],
        right=[{"chebi_id": CHEBI_GLYCEOLLIN_PRODUCT}, {"chebi_id": CHEBI_H2O}, {"chebi_id": CHEBI_H_PLUS}],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_feruloyl_coa_hydratase_lyase_requires_water() -> None:
    cls = FeruloylCoAHydrataseLyaseActivity()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_CAFFEOYL_COA}],
        right=[{"chebi_id": CHEBI_DIHYDROXYBENZALDEHYDE}, {"chebi_id": CHEBI_ACETYL_COA}],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_amine_n_methyltransferase_requires_sah_product() -> None:
    cls = AmineNMethyltransferaseActivity()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_GENERIC_PRIMARY_AMINE}, {"chebi_id": CHEBI_SAM_ZWITTERION}],
        right=[{"chebi_id": CHEBI_GENERIC_METHYLATED_AMINE}, {"chebi_id": CHEBI_H_PLUS}],
    )
    assert cls.check_membership_impl(reaction).is_member is False
