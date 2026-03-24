"""Focused tests for the next support-2 GO batch."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from typing_extensions import TypedDict

from autarch.datamodel import Participant, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS, CHEBI_UDP
from autarch.ontology.acetylajmaline_esterase_activity import AcetylajmalineEsteraseActivity
from autarch.ontology.d_arabinitol_dehydrogenase_nadp_activity import DArabinitolDehydrogenaseNADPActivity
from autarch.ontology.fifteen_oxoprostaglandin_13_reductase_nad_or_nadp_activity import FifteenOxoprostaglandin13ReductaseNADOrNADPActivity
from autarch.ontology.glc2_man9_glcnac2_oligosaccharide_glucosidase_activity import Glc2Man9GlcNAc2OligosaccharideGlucosidaseActivity
from autarch.ontology.malonate_semialdehyde_dehydrogenase_acetylating_activity import MalonateSemialdehydeDehydrogenaseAcetylatingActivity
from autarch.ontology.protein_o_acetylglucosaminyltransferase_activity import ProteinOAcetylglucosaminyltransferaseActivity
from autarch.ontology.serine_trna_ligase_activity import SerineTRNALigaseActivity
from autarch.ontology.short_chain_fatty_acyl_coa_dehydrogenase_activity import ShortChainFattyAcylCoADehydrogenaseActivity

CHEBI_SERINE = "CHEBI:33384"
CHEBI_TRNA = "CHEBI:78442"
CHEBI_SERYL_TRNA = "CHEBI:78533"
CHEBI_OXIDIZED_ETF = "CHEBI:57692"
CHEBI_BUTYRYL_COA = "CHEBI:57371"
CHEBI_CROTONOYL_COA = "CHEBI:57332"
CHEBI_MALONATE_SEMIALDEHYDE = "CHEBI:33190"
CHEBI_D_ARABINITOL = "CHEBI:18333"
CHEBI_RIBULOSE = "CHEBI:17173"
CHEBI_ACETYLAJMALINE = "CHEBI:58679"
CHEBI_AJMALINE = "CHEBI:58567"
CHEBI_ACETATE = "CHEBI:30089"
CHEBI_UDP_GLCNAC = "CHEBI:57705"
CHEBI_SERYL_PROTEIN = "CHEBI:29999"
CHEBI_GLCNAC_SERYL_PROTEIN = "CHEBI:90838"
CHEBI_GLC2MAN9 = "CHEBI:59082"
CHEBI_GLCMAN9 = "CHEBI:59080"
CHEBI_BETA_D_GLUCOSE = "CHEBI:15903"


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
            return Reaction(
                left_participants=[
                    Participant(
                        chebi_id=participant.get("chebi_id"),
                        smiles=participant.get("smiles") or chebi_to_smiles.get(participant.get("chebi_id")),
                        inchi=participant.get("inchi"),
                        location=participant.get("location"),
                        name=participant.get("name"),
                        polymer_index=participant.get("polymer_index"),
                        polymer_type=participant.get("polymer_type"),
                        count=participant.get("count") or 1,
                    )
                    for participant in reaction["left_participants"]
                ],
                right_participants=[
                    Participant(
                        chebi_id=participant.get("chebi_id"),
                        smiles=participant.get("smiles") or chebi_to_smiles.get(participant.get("chebi_id")),
                        inchi=participant.get("inchi"),
                        location=participant.get("location"),
                        name=participant.get("name"),
                        polymer_index=participant.get("polymer_index"),
                        polymer_type=participant.get("polymer_type"),
                        count=participant.get("count") or 1,
                    )
                    for participant in reaction["right_participants"]
                ],
                label=term.get("label", ""),
            )
    raise AssertionError(f"Missing cached RHEA reaction: {rhea_id}")



def _reaction(left: list[ParticipantInput], right: list[ParticipantInput]) -> Reaction:
    return Reaction(
        left_participants=[Participant(**participant) for participant in left],
        right_participants=[Participant(**participant) for participant in right],
    )


@pytest.mark.parametrize(
    ("classifier", "rhea_id"),
    [
        (SerineTRNALigaseActivity, "RHEA:12292"),
        (ShortChainFattyAcylCoADehydrogenaseActivity, "RHEA:24004"),
        (MalonateSemialdehydeDehydrogenaseAcetylatingActivity, "RHEA:22988"),
        (DArabinitolDehydrogenaseNADPActivity, "RHEA:11868"),
        (AcetylajmalineEsteraseActivity, "RHEA:22124"),
        (FifteenOxoprostaglandin13ReductaseNADOrNADPActivity, "RHEA:11912"),
        (ProteinOAcetylglucosaminyltransferaseActivity, "RHEA:48908"),
        (Glc2Man9GlcNAc2OligosaccharideGlucosidaseActivity, "RHEA:55996"),
    ],
)
def test_batch30_positive_cached_examples(classifier, rhea_id) -> None:
    reaction = _load_rhea_reaction(rhea_id)
    assert classifier().check_membership(reaction).is_member is True


def test_serine_trna_ligase_requires_atp() -> None:
    cls = SerineTRNALigaseActivity()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_TRNA}, {"chebi_id": CHEBI_SERINE}],
        right=[{"chebi_id": CHEBI_SERYL_TRNA}],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_short_chain_fatty_acyl_coa_dehydrogenase_requires_etf() -> None:
    cls = ShortChainFattyAcylCoADehydrogenaseActivity()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_BUTYRYL_COA}, {"chebi_id": CHEBI_H_PLUS}],
        right=[{"chebi_id": CHEBI_CROTONOYL_COA}],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_d_arabinitol_dehydrogenase_requires_nadp() -> None:
    cls = DArabinitolDehydrogenaseNADPActivity()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_D_ARABINITOL}],
        right=[{"chebi_id": CHEBI_RIBULOSE}],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_protein_o_acetylglucosaminyltransferase_requires_udp_glcnac() -> None:
    cls = ProteinOAcetylglucosaminyltransferaseActivity()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_SERYL_PROTEIN}],
        right=[{"chebi_id": CHEBI_GLCNAC_SERYL_PROTEIN}, {"chebi_id": CHEBI_UDP}, {"chebi_id": CHEBI_H_PLUS}],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_glc2_man9_glcnac2_oligosaccharide_glucosidase_requires_glucose_release() -> None:
    cls = Glc2Man9GlcNAc2OligosaccharideGlucosidaseActivity()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_GLC2MAN9}, {"chebi_id": CHEBI_H2O}],
        right=[{"chebi_id": CHEBI_GLCMAN9}],
    )
    assert cls.check_membership_impl(reaction).is_member is False
