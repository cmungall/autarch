"""Focused tests for the next cache-backed GO support batch."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from typing_extensions import TypedDict

from autarch.datamodel import Participant, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS, CHEBI_NADPH, CHEBI_NADP_PLUS, CHEBI_O2, CHEBI_PHOSPHATE
from autarch.ontology.coniferyl_aldehyde_dehydrogenase_nad_or_nadp_activity import ConiferylAldehydeDehydrogenaseNADOrNADPActivity
from autarch.ontology.estradiol_17_beta_dehydrogenase_nad_or_nadp_activity import Estradiol17BetaDehydrogenaseNADOrNADPActivity
from autarch.ontology.l_aminoadipate_semialdehyde_dehydrogenase_nad_or_nadp_activity import LAminoadipateSemialdehydeDehydrogenaseNADOrNADPActivity
from autarch.ontology.l_aspartate_dehydrogenase_nad_or_nadp_activity import LAspartateDehydrogenaseNADOrNADPActivity
from autarch.ontology.leucyl_trna_protein_transferase_activity import LeucylTRNAProteinTransferaseActivity
from autarch.ontology.mannosyl_oligosaccharide_12_alpha_mannosidase_activity import MannosylOligosaccharide12AlphaMannosidaseActivity
from autarch.ontology.nadh_or_nadph_oxidase_h2o2_forming_activity import NADHOrNADPHOxidaseH2O2FormingActivity
from autarch.ontology.nitroquinoline_n_oxide_reductase_nadh_or_nadph_activity import NitroquinolineNOxideReductaseNADHOrNADPHActivity
from autarch.ontology.protein_histidine_phosphatase_activity import ProteinHistidinePhosphataseActivity
from autarch.ontology.short_chain_2_methyl_fatty_acyl_coa_dehydrogenase_activity import ShortChain2MethylFattyAcylCoADehydrogenaseActivity
from autarch.ontology.succinyl_coa_3_oxo_acid_coa_transferase_activity import SuccinylCoA3OxoAcidCoATransferaseActivity
from autarch.ontology.three_hydroxyphenylacetate_6_hydroxylase_activity import ThreeHydroxyphenylacetate6HydroxylaseActivity

CHEBI_OXIDIZED_ETF = "CHEBI:57692"
CHEBI_REDUCED_ETF = "CHEBI:58307"
CHEBI_TWO_METHYLBUTANOYL_COA = "CHEBI:57336"
CHEBI_TWO_METHYLBUTENOYL_COA = "CHEBI:57337"
CHEBI_L_ALLYSINE = "CHEBI:58321"
CHEBI_L_2_AMINOADIPATE = "CHEBI:58672"
CHEBI_ESTRADIOL = "CHEBI:16469"
CHEBI_ESTRONE = "CHEBI:17263"
CHEBI_SUCCINYL_COA = "CHEBI:57292"
CHEBI_ACETOACETATE = "CHEBI:13705"
CHEBI_ACETOACETYL_COA = "CHEBI:57286"
CHEBI_SUCCINATE = "CHEBI:30031"
CHEBI_L_ASPARTATE = "CHEBI:29991"
CHEBI_OXALOACETATE = "CHEBI:16452"
CHEBI_AMMONIUM = "CHEBI:28938"
CHEBI_3_HYDROXYPHENYLACETATE = "CHEBI:58149"
CHEBI_HOMOGENTISATE = "CHEBI:16169"
CHEBI_CONIFERYL_ALDEHYDE = "CHEBI:16547"
CHEBI_FERULATE = "CHEBI:29749"
CHEBI_HISTIDYL_PROTEIN = "CHEBI:29979"
CHEBI_NTELE_PHOSPHOHISTIDINE = "CHEBI:83586"
CHEBI_LEUCYL_TRNA = "CHEBI:78494"
CHEBI_TRNA_LEU = "CHEBI:78442"
CHEBI_LYSYL_PROTEIN = "CHEBI:65249"
CHEBI_LEUCYL_LYSYL_PROTEIN = "CHEBI:133043"
CHEBI_NITROQUINOLINE_N_OXIDE = "CHEBI:16907"
CHEBI_HYDROXYAMINO_QUINOLINE_N_OXIDE = "CHEBI:28469"
CHEBI_MAN9_GLYCAN = "CHEBI:139493"
CHEBI_MAN5_GLYCAN = "CHEBI:59087"
CHEBI_BETA_D_MANNOSE = "CHEBI:28563"


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
                    count=participant.get("count") or 1,
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
                    count=participant.get("count") or 1,
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
        (ShortChain2MethylFattyAcylCoADehydrogenaseActivity, "RHEA:43780"),
        (LAminoadipateSemialdehydeDehydrogenaseNADOrNADPActivity, "RHEA:12304"),
        (Estradiol17BetaDehydrogenaseNADOrNADPActivity, "RHEA:24612"),
        (SuccinylCoA3OxoAcidCoATransferaseActivity, "RHEA:25480"),
        (NADHOrNADPHOxidaseH2O2FormingActivity, "RHEA:11260"),
        (LAspartateDehydrogenaseNADOrNADPActivity, "RHEA:11784"),
        (ThreeHydroxyphenylacetate6HydroxylaseActivity, "RHEA:22204"),
        (ConiferylAldehydeDehydrogenaseNADOrNADPActivity, "RHEA:23964"),
        (ProteinHistidinePhosphataseActivity, "RHEA:47960"),
        (LeucylTRNAProteinTransferaseActivity, "RHEA:12340"),
        (NitroquinolineNOxideReductaseNADHOrNADPHActivity, "RHEA:22476"),
        (MannosylOligosaccharide12AlphaMannosidaseActivity, "RHEA:56008"),
    ],
)
def test_batch29_positive_cached_examples(classifier, rhea_id) -> None:
    reaction = _load_rhea_reaction(rhea_id)
    assert classifier().check_membership(reaction).is_member is True


def test_short_chain_2_methyl_fatty_acyl_coa_dehydrogenase_requires_etf_pair() -> None:
    cls = ShortChain2MethylFattyAcylCoADehydrogenaseActivity()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_TWO_METHYLBUTANOYL_COA}],
        right=[{"chebi_id": CHEBI_TWO_METHYLBUTENOYL_COA}],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_l_aminoadipate_semialdehyde_dehydrogenase_requires_water() -> None:
    cls = LAminoadipateSemialdehydeDehydrogenaseNADOrNADPActivity()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_L_ALLYSINE}, {"chebi_id": CHEBI_NADP_PLUS}],
        right=[{"chebi_id": CHEBI_L_2_AMINOADIPATE}, {"chebi_id": CHEBI_NADPH}, {"chebi_id": CHEBI_H_PLUS, "count": 2}],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_succinyl_coa_3_oxo_acid_coa_transferase_requires_succinate() -> None:
    cls = SuccinylCoA3OxoAcidCoATransferaseActivity()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_ACETOACETATE}, {"chebi_id": CHEBI_SUCCINYL_COA}],
        right=[{"chebi_id": CHEBI_ACETOACETYL_COA}],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_nadh_or_nadph_oxidase_h2o2_forming_requires_hydrogen_peroxide() -> None:
    cls = NADHOrNADPHOxidaseH2O2FormingActivity()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_NADPH}, {"chebi_id": CHEBI_O2}, {"chebi_id": CHEBI_H_PLUS}],
        right=[{"chebi_id": CHEBI_NADP_PLUS}],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_three_hydroxyphenylacetate_6_hydroxylase_requires_oxygen() -> None:
    cls = ThreeHydroxyphenylacetate6HydroxylaseActivity()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_3_HYDROXYPHENYLACETATE}, {"chebi_id": CHEBI_NADPH}, {"chebi_id": CHEBI_H_PLUS}],
        right=[{"chebi_id": CHEBI_HOMOGENTISATE}, {"chebi_id": CHEBI_NADP_PLUS}, {"chebi_id": CHEBI_H2O}],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_protein_histidine_phosphatase_requires_water() -> None:
    cls = ProteinHistidinePhosphataseActivity()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_NTELE_PHOSPHOHISTIDINE}],
        right=[{"chebi_id": CHEBI_HISTIDYL_PROTEIN}, {"chebi_id": CHEBI_PHOSPHATE}],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_leucyl_trna_protein_transferase_requires_trna_product() -> None:
    cls = LeucylTRNAProteinTransferaseActivity()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_LEUCYL_TRNA}, {"chebi_id": CHEBI_LYSYL_PROTEIN}],
        right=[{"chebi_id": CHEBI_LEUCYL_LYSYL_PROTEIN}, {"chebi_id": CHEBI_H_PLUS}],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_nitroquinoline_n_oxide_reductase_requires_two_nicotinamide_equivalents() -> None:
    cls = NitroquinolineNOxideReductaseNADHOrNADPHActivity()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_NITROQUINOLINE_N_OXIDE}, {"chebi_id": CHEBI_NADPH}, {"chebi_id": CHEBI_H_PLUS}, {"chebi_id": CHEBI_H_PLUS}],
        right=[{"chebi_id": CHEBI_HYDROXYAMINO_QUINOLINE_N_OXIDE}, {"chebi_id": CHEBI_NADP_PLUS}, {"chebi_id": CHEBI_H2O}],
    )
    assert cls.check_membership_impl(reaction).is_member is False


def test_mannosyl_oligosaccharide_12_alpha_mannosidase_requires_mannose_release() -> None:
    cls = MannosylOligosaccharide12AlphaMannosidaseActivity()
    reaction = _reaction(
        left=[{"chebi_id": CHEBI_MAN9_GLYCAN}, {"chebi_id": CHEBI_H2O, "count": 4}],
        right=[{"chebi_id": CHEBI_MAN5_GLYCAN}],
    )
    assert cls.check_membership_impl(reaction).is_member is False
