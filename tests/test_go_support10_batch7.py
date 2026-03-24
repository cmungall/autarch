from pathlib import Path

import pytest

from autarch.datamodel import RheaTerm
from autarch.evaluation import evaluate_reaction_class
from autarch.ontology.acyl_coa_desaturase_activity import AcylCoADesaturaseActivity
from autarch.ontology.glucosidase_activity import GlucosidaseActivity
from autarch.ontology.glycerophospholipase_activity import GlycerophospholipaseActivity
from autarch.ontology.inositol_phosphate_kinase_activity import (
    InositolPhosphateKinaseActivity,
)
from autarch.ontology.inositol_phosphate_phosphatase_activity import (
    InositolPhosphatePhosphataseActivity,
)
from autarch.ontology.lipid_kinase_activity import LipidKinaseActivity
from autarch.ontology.phosphatidylinositol_phosphate_phosphatase_activity import (
    PhosphatidylinositolPhosphatePhosphataseActivity,
)
from autarch.ontology.phospholipase_activity import PhospholipaseActivity
from autarch.ontology.prenyl_diphosphate_synthase_activity import (
    PrenylDiphosphateSynthaseActivity,
)


def _load_rhea_reaction(rhea_id: str):
    for line in Path("cache/rhea_reactions.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        term = RheaTerm.model_validate_json(line)
        if term.rhea_id == rhea_id:
            assert term.reaction is not None
            return term.reaction
    raise AssertionError(f"Missing cached RHEA reaction: {rhea_id}")


@pytest.mark.parametrize(
    ("classifier", "rhea_id"),
    [
        (PrenylDiphosphateSynthaseActivity, "RHEA:22408"),
        (InositolPhosphatePhosphataseActivity, "RHEA:19797"),
        (PhosphatidylinositolPhosphatePhosphataseActivity, "RHEA:17193"),
        (LipidKinaseActivity, "RHEA:12709"),
        (InositolPhosphateKinaseActivity, "RHEA:11020"),
        (GlycerophospholipaseActivity, "RHEA:10604"),
        (GlucosidaseActivity, "RHEA:11956"),
        (AcylCoADesaturaseActivity, "RHEA:19721"),
    ],
)
def test_go_support10_batch7_positive_examples(classifier, rhea_id):
    reaction = _load_rhea_reaction(rhea_id)
    assert classifier().check_membership_impl(reaction).is_member is True


def test_phospholipase_activity_wraps_glycerophospholipase():
    reaction = _load_rhea_reaction("RHEA:10604")
    assert PhospholipaseActivity().check_membership_impl(reaction).is_member is True


@pytest.mark.parametrize(
    "class_name",
    [
        "PrenylDiphosphateSynthaseActivity",
        "InositolPhosphatePhosphataseActivity",
        "PhosphatidylinositolPhosphatePhosphataseActivity",
        "LipidKinaseActivity",
        "InositolPhosphateKinaseActivity",
        "GlycerophospholipaseActivity",
        "GlucosidaseActivity",
        "AcylCoADesaturaseActivity",
    ],
)
def test_go_support10_batch7_eval_smoke(class_name):
    metrics = evaluate_reaction_class(class_name, cache_dir="cache", go_only=True)
    assert metrics.true_positives > 0
