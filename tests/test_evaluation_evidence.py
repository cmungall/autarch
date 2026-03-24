"""Tests for class-specific evaluation evidence handling."""

from __future__ import annotations

import json

import pytest

from autarch.datamodel import Participant, PolymerType, Reaction
from autarch.evaluation import ec_matches_prefix, evaluate_reaction_class
from autarch.ontology.atp_hydrolysis import ATPHydrolysis
from autarch.ontology.carbonate_dehydratase import CarbonateDehydratase
from autarch.ontology.glycosyltransferase import Glycosyltransferase
from autarch.ontology.hydrolase import Hydrolase
from autarch.ontology.methyltransferase import Methyltransferase
from autarch.ontology.polysialic_acid_o_acetyltransferase import (
    PolysialicAcidOAcetyltransferase,
)
from autarch.ontology.primary_active_transmembrane_transporter import (
    PrimaryActiveTransmembraneTransporter,
)
from autarch.ontology.ribulose_bisphosphate_carboxylase import (
    RibuloseBisphosphateCarboxylase,
)
from autarch.ontology.sialyltransferase import Sialyltransferase
from autarch.ontology.sumo_transferase import SUMOTransferase
from autarch.ontology.superoxide_dismutase import SuperoxideDismutase
from autarch.ontology.ubiquitin_protein_ligase import UbiquitinProteinLigase

CHEBI_CMP_NEUAC = "CHEBI:57812"
CHEBI_CMP = "CHEBI:60377"
CHEBI_H_PLUS = "CHEBI:15378"


def test_sialyltransferase_supports_identifier_only_evaluation() -> None:
    reaction = Reaction(
        left_participants=[
            Participant(chebi_id=CHEBI_CMP_NEUAC, name="CMP-Neu5Ac", smiles="CC"),
            Participant(name="generic acceptor"),
        ],
        right_participants=[
            Participant(chebi_id=CHEBI_CMP, name="CMP", smiles="C"),
            Participant(chebi_id=CHEBI_H_PLUS, name="hydron", smiles="[H+]"),
            Participant(name="sialylated product"),
        ],
    )

    assert Sialyltransferase.supports_evaluation(reaction) is True
    assert Hydrolase.supports_evaluation(reaction) is False


@pytest.mark.parametrize(
    "reaction_class",
    [
        ATPHydrolysis,
        CarbonateDehydratase,
        Glycosyltransferase,
        Methyltransferase,
        PrimaryActiveTransmembraneTransporter,
        RibuloseBisphosphateCarboxylase,
        SUMOTransferase,
        SuperoxideDismutase,
        UbiquitinProteinLigase,
    ],
)
def test_identifier_driven_classes_support_identifier_only_evaluation(reaction_class) -> None:
    reaction = Reaction(
        left_participants=[
            Participant(chebi_id="CHEBI:15422", name="ATP"),
            Participant(name="generic substrate"),
        ],
        right_participants=[
            Participant(chebi_id="CHEBI:16761", name="ADP"),
            Participant(name="generic product"),
        ],
    )

    assert reaction_class.supports_evaluation(reaction) is True


def test_polysialic_o_acetyltransferase_supports_polymer_identifier_evidence() -> None:
    reaction = Reaction(
        left_participants=[
            Participant(
                name="sialic acid polymer",
                polymer_type=PolymerType.SIALIC_ACID_POLYMER,
                polymer_index="n",
            ),
            Participant(chebi_id="CHEBI:57288", name="acetyl-CoA"),
        ],
        right_participants=[
            Participant(
                name="acetylated sialic acid polymer",
                polymer_type=PolymerType.SIALIC_ACID_POLYMER,
                polymer_index="n",
            ),
            Participant(chebi_id="CHEBI:57287", name="CoA"),
        ],
    )

    assert PolysialicAcidOAcetyltransferase.supports_evaluation(reaction) is True


def test_evaluate_reaction_class_uses_class_specific_evidence(tmp_path) -> None:
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()

    go_terms = [
        {
            "go_id": "GO:0008373",
            "label": "sialyltransferase activity",
            "definition": "Catalysis of the transfer of sialic acid to an acceptor molecule.",
            "rhea_ids": [],
            "ec_numbers": ["2.4.3.-"],
            "ancestors": ["GO:0008373", "GO:0016757"],
        },
        {
            "go_id": "GO:0001665",
            "label": "alpha-N-acetylgalactosaminide alpha-2,6-sialyltransferase activity",
            "definition": "Catalysis of transfer from CMP-N-acetylneuraminate to an acceptor.",
            "rhea_ids": ["RHEA:11136"],
            "ec_numbers": ["2.4.3.3"],
            "ancestors": ["GO:0001665", "GO:0008373", "GO:0016757"],
        },
    ]
    with (cache_dir / "go_terms.jsonl").open("w") as stream:
        for term in go_terms:
            stream.write(json.dumps(term) + "\n")

    reaction_record = {
        "rhea_id": "RHEA:11136",
        "label": "CMP-Neu5Ac + generic acceptor = CMP + sialylated product + H(+)",
        "reaction": {
            "left_participants": [
                {"chebi_id": CHEBI_CMP_NEUAC, "name": "CMP-Neu5Ac", "smiles": "CC"},
                {"name": "generic acceptor"},
            ],
            "right_participants": [
                {"chebi_id": CHEBI_CMP, "name": "CMP", "smiles": "C"},
                {"chebi_id": CHEBI_H_PLUS, "name": "hydron", "smiles": "[H+]"},
                {"name": "sialylated product"},
            ],
        },
        "ec_numbers": ["2.4.3.3"],
        "go_terms": ["GO:0001665"],
        "direction": "bidirectional",
        "parent_rhea_id": None,
    }
    (cache_dir / "rhea_reactions.jsonl").write_text(json.dumps(reaction_record) + "\n")
    (cache_dir / "chebi_smiles.json").write_text("{}")

    metrics = evaluate_reaction_class(
        "Sialyltransferase",
        cache_dir=str(cache_dir),
        go_only=True,
    )

    assert metrics.true_positives == 1
    assert metrics.false_positives == 0
    assert metrics.false_negatives == 0


def test_ec_matches_prefix_handles_exact_ec_numbers() -> None:
    assert ec_matches_prefix("1.1.1.296", "1.1.1.296") is True
    assert ec_matches_prefix("2.7.4.6", "2.7.4.6") is True
    assert ec_matches_prefix("2.7.4.6", "2.7.4.-") is True
    assert ec_matches_prefix("2.7.4.6", "2.7.5.-") is False
