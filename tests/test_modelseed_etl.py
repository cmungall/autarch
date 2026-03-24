"""Tests for ModelSEED ETL and ChEBI bridge utilities."""

import json
from pathlib import Path

from autarch.etl.modelseed_etl import (
    ModelSeedETL,
    cache_modelseed_dataset,
    load_modelseed_to_chebi_mappings,
    summarize_modelseed_cache,
)


COMPOUND_HEADER = [
    "id",
    "abbreviation",
    "name",
    "formula",
    "mass",
    "source",
    "inchikey",
    "charge",
    "is_core",
    "is_obsolete",
    "linked_compound",
    "is_cofactor",
    "deltag",
    "deltagerr",
    "pka",
    "pkb",
    "abstract_compound",
    "comprised_of",
    "aliases",
    "smiles",
    "notes",
]

REACTION_HEADER = [
    "id",
    "abbreviation",
    "name",
    "code",
    "stoichiometry",
    "is_transport",
    "equation",
    "definition",
    "reversibility",
    "direction",
    "abstract_reaction",
    "pathways",
    "aliases",
    "ec_numbers",
    "deltag",
    "deltagerr",
    "compound_ids",
    "status",
    "is_obsolete",
    "linked_reaction",
    "notes",
    "source",
]


def _write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    lines = ["\t".join(header)]
    lines.extend("\t".join(row) for row in rows)
    path.write_text("\n".join(lines) + "\n")


def _write_modelseed_test_files(base_dir: Path) -> None:
    aliases_dir = base_dir / "Aliases"
    aliases_dir.mkdir(parents=True, exist_ok=True)

    _write_tsv(
        base_dir / "compounds.tsv",
        COMPOUND_HEADER,
        [
            [
                "cpd00001",
                "h2o",
                "water",
                "H2O",
                "18.0",
                "Primary Database",
                "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
                "0",
                "1",
                "0",
                "null",
                "0",
                "0",
                "0",
                "null",
                "null",
                "null",
                "null",
                "Name: water; H2O|KEGG: C00001",
                "O",
                "",
            ],
            [
                "cpd00002",
                "ala",
                "alanine",
                "C3H7NO2",
                "89.0",
                "Primary Database",
                "QNAYBMKLOCPYGJ-REOHCLBHSA-N",
                "0",
                "1",
                "0",
                "null",
                "0",
                "0",
                "0",
                "null",
                "null",
                "null",
                "null",
                "Name: alanine; L-alanine",
                "C[C@H](N)C(=O)O",
                "",
            ],
            [
                "cpd00003",
                "h",
                "proton",
                "H",
                "1.0",
                "Primary Database",
                "",
                "1",
                "1",
                "0",
                "null",
                "0",
                "0",
                "0",
                "null",
                "null",
                "null",
                "null",
                "Name: proton; H+",
                "",
                "",
            ],
            [
                "cpd00004",
                "nh3",
                "ammonia",
                "NH3",
                "17.0",
                "Primary Database",
                "",
                "0",
                "1",
                "0",
                "null",
                "0",
                "0",
                "0",
                "null",
                "null",
                "null",
                "null",
                "Name: ammonia",
                "",
                "",
            ],
        ],
    )

    _write_tsv(
        base_dir / "reactions.tsv",
        REACTION_HEADER,
        [
            [
                "rxn00001",
                "H2Ot",
                "water transport",
                "(1) cpd00001[0] <=> (1) cpd00001[1]",
                '-1:cpd00001:0:0:"Water";1:cpd00001:1:0:"Water"',
                "1",
                "(1) cpd00001[0] <=> (1) cpd00001[1]",
                "(1) Water[0] <=> (1) Water[1]",
                "=",
                "=",
                "null",
                "null",
                "Name: water transport",
                "",
                "0",
                "0",
                "cpd00001",
                "OK",
                "0",
                "null",
                "",
                "Primary Database",
            ]
        ],
    )

    _write_tsv(
        aliases_dir / "Unique_ModelSEED_Reaction_Aliases.txt",
        ["ModelSEED ID", "External ID", "Source"],
        [
            ["rxn00001", "10000", "rhea"],
            ["rxn00001", "10004", "rhea"],
            ["rxn00001", "R00001", "KEGG"],
        ],
    )

    _write_tsv(
        aliases_dir / "Unique_ModelSEED_Reaction_ECs.txt",
        ["ModelSEED ID", "External ID", "Source"],
        [["rxn00001", "1.1.1.1", "Enzyme Class"]],
    )


def _write_chebi_lookup_cache(cache_dir: Path) -> None:
    (cache_dir / "smiles_to_chebi.json").write_text(
        json.dumps(
            {
                "O": "CHEBI:15377",
                "CC(N)C(=O)O": "CHEBI:16449",
            }
        )
    )
    (cache_dir / "chebi_names.json").write_text(
        json.dumps(
            {
                "CHEBI:15377": "water",
                "CHEBI:16449": "alanine",
                "CHEBI:16134": "ammonia",
            }
        )
    )


def test_modelseed_compound_bridge_maps_by_multiple_methods(tmp_path: Path):
    """Bridge compounds via exact SMILES, no-stereo SMILES, and name fallback."""
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    _write_chebi_lookup_cache(cache_dir)

    modelseed_dir = cache_dir / "modelseed"
    _write_modelseed_test_files(modelseed_dir)

    etl = ModelSeedETL(cache_dir=modelseed_dir)
    compounds = etl.load_compounds()
    mappings = etl.build_chebi_bridge(compounds, cache_dir=cache_dir)

    assert mappings["cpd00001"].chebi_id == "CHEBI:15377"
    assert mappings["cpd00001"].mapping_method == "exact_smiles"

    assert mappings["cpd00002"].chebi_id == "CHEBI:16449"
    assert mappings["cpd00002"].mapping_method == "no_stereo_smiles"

    assert mappings["cpd00003"].chebi_id == "CHEBI:15378"
    assert mappings["cpd00003"].mapping_method == "common_name"

    assert mappings["cpd00004"].chebi_id == "CHEBI:16134"
    assert mappings["cpd00004"].mapping_method == "name_cache"


def test_modelseed_reaction_records_include_rhea_ec_and_autarch_reaction(tmp_path: Path):
    """Build ModelSEED reaction records enriched with Rhea IDs and ChEBI-mapped participants."""
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    _write_chebi_lookup_cache(cache_dir)

    modelseed_dir = cache_dir / "modelseed"
    _write_modelseed_test_files(modelseed_dir)

    etl = ModelSeedETL(cache_dir=modelseed_dir)
    compounds = etl.load_compounds()
    mappings = etl.build_chebi_bridge(compounds, cache_dir=cache_dir)
    reactions = etl.load_reactions(compounds, mappings)

    record = reactions["rxn00001"]
    assert record.rhea_ids == ["RHEA:10000", "RHEA:10004"]
    assert record.ec_numbers == ["1.1.1.1"]
    assert record.is_transport is True
    assert record.fully_mapped_to_chebi is True
    assert record.compound_ids == ["cpd00001"]

    assert record.reaction is not None
    assert record.reaction.is_transport_reaction() is True
    assert record.reaction.left_participants[0].chebi_id == "CHEBI:15377"
    assert record.reaction.left_participants[0].location == "0"
    assert record.reaction.right_participants[0].location == "1"


def test_modelseed_fractional_stoichiometry_is_preserved_without_crashing(tmp_path: Path):
    """Fractional ModelSEED coefficients should be preserved as stoichiometry strings."""
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    _write_chebi_lookup_cache(cache_dir)

    modelseed_dir = cache_dir / "modelseed"
    _write_modelseed_test_files(modelseed_dir)
    _write_tsv(
        modelseed_dir / "reactions.tsv",
        REACTION_HEADER,
        [
            [
                "rxn00002",
                "half-water",
                "fractional water reaction",
                "",
                '-0.5:cpd00001:0:0:"Water";1:cpd00004:0:0:"ammonia"',
                "0",
                "",
                "fractional water reaction",
                "=",
                "=",
                "null",
                "null",
                "",
                "",
                "0",
                "0",
                "cpd00001;cpd00004",
                "OK",
                "0",
                "null",
                "",
                "Primary Database",
            ]
        ],
    )
    _write_tsv(
        modelseed_dir / "Aliases" / "Unique_ModelSEED_Reaction_Aliases.txt",
        ["ModelSEED ID", "External ID", "Source"],
        [],
    )
    _write_tsv(
        modelseed_dir / "Aliases" / "Unique_ModelSEED_Reaction_ECs.txt",
        ["ModelSEED ID", "External ID", "Source"],
        [],
    )

    etl = ModelSeedETL(cache_dir=modelseed_dir)
    compounds = etl.load_compounds()
    mappings = etl.build_chebi_bridge(compounds, cache_dir=cache_dir)
    reactions = etl.load_reactions(compounds, mappings)

    record = reactions["rxn00002"]
    assert record.reaction is not None
    assert record.reaction.left_participants[0].count == 1
    assert record.reaction.left_participants[0].stoichiometry == "0.5"
    assert record.reaction.right_participants[0].stoichiometry is None


def test_cache_modelseed_dataset_writes_jsonl_outputs(tmp_path: Path, monkeypatch):
    """Cache helper writes compound, mapping, and reaction artifacts."""
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    _write_chebi_lookup_cache(cache_dir)

    modelseed_dir = cache_dir / "modelseed"
    _write_modelseed_test_files(modelseed_dir)

    def _skip_download(self):
        return None

    monkeypatch.setattr(ModelSeedETL, "download_required_files", _skip_download)

    summary = cache_modelseed_dataset(cache_dir=cache_dir)

    assert summary["total_compounds"] == 4
    assert summary["mapped_compounds"] == 4
    assert summary["total_reactions"] == 1
    assert summary["reactions_with_rhea"] == 1
    assert summary["fully_mapped_reactions"] == 1

    compounds_file = cache_dir / "modelseed_compounds.jsonl"
    mappings_file = cache_dir / "modelseed_to_chebi.jsonl"
    reactions_file = cache_dir / "modelseed_reactions.jsonl"

    assert compounds_file.exists()
    assert mappings_file.exists()
    assert reactions_file.exists()

    loaded_mappings = load_modelseed_to_chebi_mappings(mappings_file)
    assert loaded_mappings["cpd00001"].chebi_id == "CHEBI:15377"

    reaction_lines = reactions_file.read_text().strip().splitlines()
    assert len(reaction_lines) == 1
    reaction_payload = json.loads(reaction_lines[0])
    assert reaction_payload["modelseed_id"] == "rxn00001"
    assert reaction_payload["rhea_ids"] == ["RHEA:10000", "RHEA:10004"]


def test_summarize_modelseed_cache_reports_mapping_coverage(tmp_path: Path, monkeypatch):
    """Summarize cached ModelSEED coverage from local cache artifacts."""
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    _write_chebi_lookup_cache(cache_dir)

    modelseed_dir = cache_dir / "modelseed"
    _write_modelseed_test_files(modelseed_dir)

    def _skip_download(self):
        return None

    monkeypatch.setattr(ModelSeedETL, "download_required_files", _skip_download)

    cache_modelseed_dataset(cache_dir=cache_dir)
    summary = summarize_modelseed_cache(cache_dir)

    assert summary["total_compounds"] == 4
    assert summary["mapped_compounds"] == 4
    assert summary["total_reactions"] == 1
    assert summary["reactions_with_rhea"] == 1
    assert summary["cached_reactions"] == 1
    assert summary["cached_reactions_with_rhea"] == 1
    assert summary["cached_fully_mapped_reactions"] == 1

    assert summary["compound_mapping_methods"]["exact_smiles"] == 1
    assert summary["compound_mapping_methods"]["no_stereo_smiles"] == 1
    assert summary["compound_mapping_methods"]["common_name"] == 1
    assert summary["compound_mapping_methods"]["name_cache"] == 1

    assert summary["reaction_alias_sources"]["rhea"]["distinct_ids"] == 1
    assert summary["reaction_alias_sources"]["rhea"]["rows"] == 2
    assert summary["reaction_alias_sources"]["KEGG"]["distinct_ids"] == 1
    assert summary["ec_coverage"]["distinct_reactions"] == 1
    assert summary["compound_alias_sources"]["KEGG"]["distinct_ids"] == 1
