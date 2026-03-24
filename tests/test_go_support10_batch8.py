from pathlib import Path
import json

import pytest

from autarch.datamodel import Participant, Reaction
from autarch.evaluation import evaluate_reaction_class
from autarch.ontology.acyl_coa_dehydrogenase_activity import (
    AcylCoADehydrogenaseActivity,
)
from autarch.ontology.amino_acid_dehydrogenase_nad_or_nadp import (
    AminoAcidDehydrogenaseNADOrNADP,
)
from autarch.ontology.beta_glucosidase_activity import BetaGlucosidaseActivity
from autarch.ontology.catalytic_activity_acting_on_a_glycoprotein import (
    CatalyticActivityActingOnAGlycoprotein,
)
from autarch.ontology.deacetylase_activity import DeacetylaseActivity
from autarch.ontology.demethylase_activity import DemethylaseActivity
from autarch.ontology.l_amino_acid_n_acetyltransferase_activity import (
    LAminoAcidNAcetyltransferaseActivity,
)
from autarch.ontology.protein_methyltransferase_activity import (
    ProteinMethyltransferaseActivity,
)
from autarch.ontology.steroid_dehydrogenase_activity import (
    SteroidDehydrogenaseActivity,
)
from autarch.ontology.steroid_dehydrogenase_activity_acting_on_the_ch_oh_group_of_donors_nad_or_nadp_as_acceptor import (
    SteroidDehydrogenaseActivityActingOnTheCHOHGroupOfDonorsNADOrNADPAsAcceptor,
)
from autarch.ontology.steroid_hydroxylase_activity import (
    SteroidHydroxylaseActivity,
)


def _load_rhea_reaction(rhea_id: str):
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
                left_participants.append(
                    Participant(
                        chebi_id=chebi_id,
                        smiles=smiles,
                        inchi=participant.get("inchi"),
                        location=participant.get("location"),
                        name=participant.get("name"),
                        polymer_index=participant.get("polymer_index"),
                        polymer_type=participant.get("polymer_type"),
                    )
                )
            for participant in reaction["right_participants"]:
                chebi_id = participant.get("chebi_id")
                smiles = participant.get("smiles") or (chebi_to_smiles.get(chebi_id) if chebi_id else None)
                right_participants.append(
                    Participant(
                        chebi_id=chebi_id,
                        smiles=smiles,
                        inchi=participant.get("inchi"),
                        location=participant.get("location"),
                        name=participant.get("name"),
                        polymer_index=participant.get("polymer_index"),
                        polymer_type=participant.get("polymer_type"),
                    )
                )
            return Reaction(
                left_participants=left_participants,
                right_participants=right_participants,
                label=term.get("label", ""),
            )
    raise AssertionError(f"Missing cached RHEA reaction: {rhea_id}")


@pytest.mark.parametrize(
    ("classifier", "rhea_id"),
    [
        (ProteinMethyltransferaseActivity, "RHEA:10024"),
        (SteroidHydroxylaseActivity, "RHEA:15629"),
        (SteroidDehydrogenaseActivity, "RHEA:14929"),
        (
            SteroidDehydrogenaseActivityActingOnTheCHOHGroupOfDonorsNADOrNADPAsAcceptor,
            "RHEA:14929",
        ),
        (DeacetylaseActivity, "RHEA:15941"),
        (DemethylaseActivity, "RHEA:13021"),
        (AminoAcidDehydrogenaseNADOrNADP, "RHEA:11156"),
        (BetaGlucosidaseActivity, "RHEA:11956"),
        (CatalyticActivityActingOnAGlycoprotein, "RHEA:11456"),
        (AcylCoADehydrogenaseActivity, "RHEA:24004"),
        (LAminoAcidNAcetyltransferaseActivity, "RHEA:24292"),
    ],
)
def test_go_support10_batch8_positive_examples(classifier, rhea_id):
    reaction = _load_rhea_reaction(rhea_id)
    assert classifier().check_membership(reaction).is_member is True


@pytest.mark.parametrize(
    "class_name",
    [
        "ProteinMethyltransferaseActivity",
        "SteroidHydroxylaseActivity",
        "SteroidDehydrogenaseActivity",
        "SteroidDehydrogenaseActivityActingOnTheCHOHGroupOfDonorsNADOrNADPAsAcceptor",
        "DeacetylaseActivity",
        "DemethylaseActivity",
        "AminoAcidDehydrogenaseNADOrNADP",
        "BetaGlucosidaseActivity",
        "CatalyticActivityActingOnAGlycoprotein",
        "AcylCoADehydrogenaseActivity",
        "LAminoAcidNAcetyltransferaseActivity",
    ],
)
def test_go_support10_batch8_eval_smoke(class_name):
    metrics = evaluate_reaction_class(class_name, cache_dir="cache", go_only=True)
    assert metrics.true_positives > 0
