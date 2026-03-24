"""ModelSEED ETL and bridge utilities.

This module downloads a small set of canonical ModelSEED biochemistry files,
maps ModelSEED compounds into the local ChEBI-centric identifier space, and
converts ModelSEED reactions into ``autarch.datamodel.Reaction`` objects.

The goal is pragmatic:
1. Keep ModelSEED as an auxiliary corpus for out-of-distribution evaluation.
2. Reuse existing ChEBI lookup caches and RDKit normalization.
3. Emit artifacts that fit the existing cache-oriented workflow.
"""

from __future__ import annotations

import csv
import json
import logging
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, TypedDict

import requests
from pydantic import BaseModel, Field
from rdkit import Chem
from rdkit.Chem.inchi import MolToInchiKey

from autarch.datamodel import Participant, Reaction
from autarch.etl.chebi_normalization import normalize_chebi_by_name
from autarch.etl.chebi_smiles import canonicalize_smiles, load_smiles_to_chebi_cache

logger = logging.getLogger(__name__)

MODELSEED_RAW_BASE = (
    "https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/"
)

MODELSEED_FILES = {
    "compounds": "compounds.tsv",
    "reactions": "reactions.tsv",
    "reaction_aliases": "Aliases/Unique_ModelSEED_Reaction_Aliases.txt",
    "reaction_ecs": "Aliases/Unique_ModelSEED_Reaction_ECs.txt",
}


class AliasSourceSummary(TypedDict):
    """Row and distinct-ID counts for one alias source."""

    rows: int
    distinct_ids: int


class ReactionECCoverageSummary(TypedDict):
    """Coverage summary for the cached reaction-EC file."""

    rows: int
    distinct_reactions: int


class ModelSeedCacheSummary(TypedDict):
    """Typed summary of cached ModelSEED bridge coverage."""

    total_compounds: int
    mapped_compounds: int
    compound_mapping_methods: dict[str, int]
    compound_alias_sources: dict[str, AliasSourceSummary]
    total_reactions: int
    reactions_with_rhea: int
    cached_reactions: int
    cached_reactions_with_rhea: int
    cached_fully_mapped_reactions: int
    reaction_alias_sources: dict[str, AliasSourceSummary]
    ec_coverage: ReactionECCoverageSummary


class ModelSeedCompound(BaseModel):
    """ModelSEED compound record."""

    modelseed_id: str
    abbreviation: str = ""
    name: str = ""
    formula: Optional[str] = None
    mass: Optional[float] = None
    inchikey: Optional[str] = None
    charge: Optional[int] = None
    aliases: dict[str, list[str]] = Field(default_factory=dict)
    smiles: Optional[str] = None
    notes: str = ""


class ModelSeedCompoundMapping(BaseModel):
    """Bridge from ModelSEED compound IDs to ChEBI IDs."""

    modelseed_id: str
    modelseed_name: str = ""
    chebi_id: Optional[str] = None
    mapping_method: Optional[str] = None


class ModelSeedStoichiometryParticipant(BaseModel):
    """Participant parsed from ModelSEED stoichiometry."""

    modelseed_id: str
    name: str = ""
    count: int = 1
    stoichiometry: Optional[str] = None
    location: Optional[str] = None


class ModelSeedReactionRecord(BaseModel):
    """ModelSEED reaction record enriched for Autarch."""

    modelseed_id: str
    abbreviation: str = ""
    name: str = ""
    definition: str = ""
    equation: str = ""
    stoichiometry: str = ""
    is_transport: bool = False
    status: str = ""
    ec_numbers: list[str] = Field(default_factory=list)
    rhea_ids: list[str] = Field(default_factory=list)
    compound_ids: list[str] = Field(default_factory=list)
    reaction: Optional[Reaction] = None
    fully_mapped_to_chebi: bool = False


def _normalize_null(value: Optional[str]) -> Optional[str]:
    """Normalize ModelSEED null-like string values."""
    if value is None:
        return None
    cleaned = value.strip()
    if not cleaned or cleaned.lower() in {"null", "none", "nan"}:
        return None
    return cleaned


def _parse_float(value: Optional[str]) -> Optional[float]:
    """Parse optional float fields."""
    cleaned = _normalize_null(value)
    if cleaned is None:
        return None
    return float(cleaned)


def _parse_int(value: Optional[str]) -> Optional[int]:
    """Parse optional integer fields."""
    cleaned = _normalize_null(value)
    if cleaned is None:
        return None
    return int(cleaned)


def _parse_bool_flag(value: Optional[str]) -> bool:
    """Parse ModelSEED boolean-ish fields."""
    cleaned = (_normalize_null(value) or "").lower()
    return cleaned in {"1", "true", "t", "yes"}


def _unique(values: Iterable[str]) -> list[str]:
    """Return values preserving order and removing duplicates."""
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        if value and value not in seen:
            seen.add(value)
            ordered.append(value)
    return ordered


def parse_modelseed_alias_blob(text: Optional[str]) -> dict[str, list[str]]:
    """Parse the pipe-delimited alias column used in ModelSEED TSV files.

    Example input:
        ``Name: water; H2O|KEGG: C00001|BiGG: h2o``
    """
    cleaned = _normalize_null(text)
    if cleaned is None:
        return {}

    aliases: dict[str, list[str]] = {}
    for group in cleaned.split("|"):
        if ":" not in group:
            continue
        source, values = group.split(":", 1)
        entries = [value.strip() for value in values.split(";") if value.strip()]
        if entries:
            aliases[source.strip()] = entries
    return aliases


def parse_modelseed_stoichiometry(
    text: Optional[str],
) -> tuple[list[ModelSeedStoichiometryParticipant], list[ModelSeedStoichiometryParticipant]]:
    """Parse ModelSEED stoichiometry into left/right participants."""
    cleaned = _normalize_null(text)
    if cleaned is None:
        return [], []

    left: list[ModelSeedStoichiometryParticipant] = []
    right: list[ModelSeedStoichiometryParticipant] = []

    for raw_part in cleaned.split(";"):
        part = raw_part.strip()
        if not part:
            continue

        coefficient, compound_id, location, _community, name = part.split(":", 4)
        coefficient_value = float(coefficient)
        count = abs(int(coefficient_value)) if coefficient_value.is_integer() else 1
        stoichiometry = (
            coefficient.strip().lstrip("+-")
            if not coefficient_value.is_integer()
            else None
        )
        participant = ModelSeedStoichiometryParticipant(
            modelseed_id=compound_id.strip(),
            name=name.strip().strip('"'),
            count=count,
            stoichiometry=stoichiometry,
            location=location.strip() if location.strip() else None,
        )

        if coefficient_value < 0:
            left.append(participant)
        else:
            right.append(participant)

    return left, right


def _build_secondary_chebi_indexes(
    smiles_to_chebi: dict[str, str],
) -> tuple[dict[str, str], dict[str, str]]:
    """Build relaxed ChEBI structure indexes from the canonical SMILES cache."""
    inchikey_to_chebi: dict[str, str] = {}
    no_stereo_to_chebi: dict[str, str] = {}

    for canonical_smiles, chebi_id in smiles_to_chebi.items():
        mol = Chem.MolFromSmiles(canonical_smiles)
        if mol is None:
            continue

        inchikey = MolToInchiKey(mol)
        no_stereo = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)

        inchikey_to_chebi.setdefault(inchikey, chebi_id)
        no_stereo_to_chebi.setdefault(no_stereo, chebi_id)

    return inchikey_to_chebi, no_stereo_to_chebi


def _build_name_lookup(chebi_to_name: dict[str, str]) -> dict[str, str]:
    """Build lowercase exact-name lookup from the cached ChEBI label map."""
    lookup: dict[str, str] = {}
    for chebi_id, name in chebi_to_name.items():
        normalized = name.strip().lower()
        if normalized:
            lookup.setdefault(normalized, chebi_id)
    return lookup


def _candidate_names(compound: ModelSeedCompound) -> list[str]:
    """Return candidate names for compound-based fallback mapping."""
    name_aliases = compound.aliases.get("Name", []) + compound.aliases.get("name", [])
    candidates = [compound.name, *name_aliases]
    return _unique([name.strip() for name in candidates if name and name.strip()])


@dataclass
class ModelSeedETL:
    """ModelSEED ETL rooted at a local cache directory."""

    cache_dir: Path = Path("cache/modelseed")
    force_download: bool = False

    def __post_init__(self) -> None:
        """Ensure cache directory exists."""
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def download_file(self, relative_path: str) -> Path:
        """Download a ModelSEED raw file if needed."""
        local_path = self.cache_dir / relative_path
        local_path.parent.mkdir(parents=True, exist_ok=True)

        if local_path.exists() and not self.force_download:
            logger.info("Using cached ModelSEED file: %s", local_path)
            return local_path

        url = f"{MODELSEED_RAW_BASE}{relative_path}"
        logger.info("Downloading ModelSEED file from %s", url)

        response = requests.get(url, stream=True)
        response.raise_for_status()

        with open(local_path, "wb") as handle:
            for chunk in response.iter_content(chunk_size=8192):
                handle.write(chunk)

        return local_path

    def download_required_files(self) -> dict[str, Path]:
        """Download the minimal set of files needed for the bridge."""
        return {
            key: self.download_file(relative_path)
            for key, relative_path in MODELSEED_FILES.items()
        }

    def _path_for(self, file_key: str) -> Path:
        """Return local path for a required file, downloading if missing."""
        relative_path = MODELSEED_FILES[file_key]
        path = self.cache_dir / relative_path
        if not path.exists():
            path = self.download_file(relative_path)
        return path

    def load_compounds(self) -> dict[str, ModelSeedCompound]:
        """Load ModelSEED compounds from TSV."""
        compounds: dict[str, ModelSeedCompound] = {}
        with open(self._path_for("compounds"), newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                record = ModelSeedCompound(
                    modelseed_id=row["id"],
                    abbreviation=row.get("abbreviation", "") or "",
                    name=row.get("name", "") or "",
                    formula=_normalize_null(row.get("formula")),
                    mass=_parse_float(row.get("mass")),
                    inchikey=_normalize_null(row.get("inchikey")),
                    charge=_parse_int(row.get("charge")),
                    aliases=parse_modelseed_alias_blob(row.get("aliases")),
                    smiles=_normalize_null(row.get("smiles")),
                    notes=row.get("notes", "") or "",
                )
                compounds[record.modelseed_id] = record
        return compounds

    def build_chebi_bridge(
        self,
        compounds: dict[str, ModelSeedCompound],
        cache_dir: str | Path = "cache",
    ) -> dict[str, ModelSeedCompoundMapping]:
        """Map ModelSEED compounds into ChEBI using local structure/name caches."""
        smiles_to_chebi, chebi_to_name = load_smiles_to_chebi_cache(str(cache_dir))
        inchikey_to_chebi, no_stereo_to_chebi = _build_secondary_chebi_indexes(
            smiles_to_chebi
        )
        name_to_chebi = _build_name_lookup(chebi_to_name)

        mappings: dict[str, ModelSeedCompoundMapping] = {}
        for compound in compounds.values():
            chebi_id: Optional[str] = None
            method: Optional[str] = None

            canonical = canonicalize_smiles(compound.smiles) if compound.smiles else None
            if canonical and canonical in smiles_to_chebi:
                chebi_id = smiles_to_chebi[canonical]
                method = "exact_smiles"

            if chebi_id is None and compound.inchikey:
                chebi_id = inchikey_to_chebi.get(compound.inchikey)
                if chebi_id:
                    method = "inchikey"

            if chebi_id is None and compound.smiles:
                mol = Chem.MolFromSmiles(compound.smiles)
                if mol is not None:
                    no_stereo = Chem.MolToSmiles(
                        mol, canonical=True, isomericSmiles=False
                    )
                    chebi_id = no_stereo_to_chebi.get(no_stereo)
                    if chebi_id:
                        method = "no_stereo_smiles"

            if chebi_id is None:
                for candidate in _candidate_names(compound):
                    normalized = candidate.lower().strip()
                    chebi_id = name_to_chebi.get(normalized)
                    if chebi_id:
                        method = "name_cache"
                        break

                    chebi_id = normalize_chebi_by_name(candidate)
                    if chebi_id:
                        method = "common_name"
                        break

            mappings[compound.modelseed_id] = ModelSeedCompoundMapping(
                modelseed_id=compound.modelseed_id,
                modelseed_name=compound.name,
                chebi_id=chebi_id,
                mapping_method=method,
            )

        return mappings

    def _load_reaction_aliases_by_source(self, source_name: str) -> dict[str, list[str]]:
        """Load reaction aliases keyed by ModelSEED ID for a specific source."""
        aliases: dict[str, list[str]] = {}
        with open(self._path_for("reaction_aliases"), newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                if row["Source"].strip().lower() != source_name.lower():
                    continue
                modelseed_id = row["ModelSEED ID"].strip()
                external_id = row["External ID"].strip()
                aliases.setdefault(modelseed_id, []).append(external_id)
        return {
            modelseed_id: _unique(values)
            for modelseed_id, values in aliases.items()
        }

    def _load_reaction_ecs(self) -> dict[str, list[str]]:
        """Load reaction EC numbers keyed by ModelSEED ID."""
        ec_map: dict[str, list[str]] = {}
        with open(self._path_for("reaction_ecs"), newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                modelseed_id = row["ModelSEED ID"].strip()
                ec_number = row["External ID"].strip()
                ec_map.setdefault(modelseed_id, []).append(ec_number)
        return {
            modelseed_id: _unique(values)
            for modelseed_id, values in ec_map.items()
        }

    def _build_reaction_object(
        self,
        stoichiometry: str,
        compounds: dict[str, ModelSeedCompound],
        mappings: dict[str, ModelSeedCompoundMapping],
        include_locations: bool,
        label: str,
    ) -> Reaction:
        """Convert parsed ModelSEED stoichiometry into an Autarch Reaction."""
        left_specs, right_specs = parse_modelseed_stoichiometry(stoichiometry)

        def _participant_from_spec(
            spec: ModelSeedStoichiometryParticipant,
        ) -> Participant:
            compound = compounds.get(spec.modelseed_id)
            mapping = mappings.get(spec.modelseed_id)
            return Participant(
                chebi_id=mapping.chebi_id if mapping else None,
                smiles=compound.smiles if compound else None,
                name=(compound.name if compound and compound.name else spec.name) or None,
                formula=compound.formula if compound else None,
                charge=compound.charge if compound else None,
                count=spec.count,
                stoichiometry=spec.stoichiometry,
                location=spec.location if include_locations else None,
            )

        return Reaction(
            left_participants=[_participant_from_spec(spec) for spec in left_specs],
            right_participants=[_participant_from_spec(spec) for spec in right_specs],
            label=label,
        )

    def load_reactions(
        self,
        compounds: dict[str, ModelSeedCompound],
        mappings: dict[str, ModelSeedCompoundMapping],
        limit: Optional[int] = None,
    ) -> dict[str, ModelSeedReactionRecord]:
        """Load ModelSEED reactions enriched with Rhea/EC links and Reaction objects."""
        rhea_aliases = self._load_reaction_aliases_by_source("rhea")
        reaction_ecs = self._load_reaction_ecs()

        reactions: dict[str, ModelSeedReactionRecord] = {}
        with open(self._path_for("reactions"), newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for index, row in enumerate(reader):
                if limit is not None and index >= limit:
                    break

                modelseed_id = row["id"]
                rhea_ids = [
                    f"RHEA:{external_id}"
                    for external_id in rhea_aliases.get(modelseed_id, [])
                ]
                ec_numbers = _unique(
                    [
                        *reaction_ecs.get(modelseed_id, []),
                        *[
                            value.strip()
                            for value in re.split(r"[;|]", row.get("ec_numbers", ""))
                            if value.strip() and value.strip().lower() != "null"
                        ],
                    ]
                )
                compound_ids = [
                    value.strip()
                    for value in (row.get("compound_ids", "") or "").split(";")
                    if value.strip() and value.strip().lower() != "null"
                ]

                definition = row.get("definition", "") or row.get("name", "") or ""
                reaction = self._build_reaction_object(
                    stoichiometry=row.get("stoichiometry", "") or "",
                    compounds=compounds,
                    mappings=mappings,
                    include_locations=_parse_bool_flag(row.get("is_transport")),
                    label=definition,
                )
                fully_mapped_to_chebi = all(
                    participant.chebi_id is not None
                    for participant in reaction.all_participants()
                )

                reactions[modelseed_id] = ModelSeedReactionRecord(
                    modelseed_id=modelseed_id,
                    abbreviation=row.get("abbreviation", "") or "",
                    name=row.get("name", "") or "",
                    definition=definition,
                    equation=row.get("equation", "") or "",
                    stoichiometry=row.get("stoichiometry", "") or "",
                    is_transport=_parse_bool_flag(row.get("is_transport")),
                    status=row.get("status", "") or "",
                    ec_numbers=ec_numbers,
                    rhea_ids=rhea_ids,
                    compound_ids=compound_ids,
                    reaction=reaction,
                    fully_mapped_to_chebi=fully_mapped_to_chebi,
                )

        return reactions


def serialize_modelseed_compounds(
    compounds: dict[str, ModelSeedCompound], output_file: str | Path
) -> None:
    """Serialize ModelSEED compounds to JSONL."""
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as handle:
        for compound in compounds.values():
            json.dump(compound.model_dump(), handle)
            handle.write("\n")


def serialize_modelseed_mappings(
    mappings: dict[str, ModelSeedCompoundMapping], output_file: str | Path
) -> None:
    """Serialize ModelSEED→ChEBI mappings to JSONL."""
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as handle:
        for mapping in mappings.values():
            json.dump(mapping.model_dump(), handle)
            handle.write("\n")


def serialize_modelseed_reactions(
    reactions: dict[str, ModelSeedReactionRecord], output_file: str | Path
) -> None:
    """Serialize ModelSEED reaction records to JSONL."""
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as handle:
        for reaction in reactions.values():
            json.dump(reaction.model_dump(), handle)
            handle.write("\n")


def load_modelseed_to_chebi_mappings(
    input_file: str | Path,
) -> dict[str, ModelSeedCompoundMapping]:
    """Load serialized ModelSEED→ChEBI mappings from JSONL."""
    mappings: dict[str, ModelSeedCompoundMapping] = {}
    with open(input_file) as handle:
        for line in handle:
            if not line.strip():
                continue
            mapping = ModelSeedCompoundMapping(**json.loads(line))
            mappings[mapping.modelseed_id] = mapping
    return mappings


def load_modelseed_reaction_records(
    input_file: str | Path,
) -> dict[str, ModelSeedReactionRecord]:
    """Load serialized ModelSEED reaction records from JSONL."""
    reactions: dict[str, ModelSeedReactionRecord] = {}
    with open(input_file) as handle:
        for line in handle:
            if not line.strip():
                continue
            reaction = ModelSeedReactionRecord(**json.loads(line))
            reactions[reaction.modelseed_id] = reaction
    return reactions


def _summarize_alias_sources(
    alias_map: dict[str, dict[str, list[str]]],
) -> dict[str, AliasSourceSummary]:
    """Summarize alias coverage as row counts and distinct ID counts."""
    rows: Counter[str] = Counter()
    distinct: Counter[str] = Counter()

    for source_map in alias_map.values():
        for source, values in source_map.items():
            if not values:
                continue
            rows[source] += len(values)
            distinct[source] += 1

    return {
        source: {"rows": rows[source], "distinct_ids": distinct[source]}
        for source in sorted(rows)
    }


def _summarize_reaction_alias_file(
    input_file: str | Path,
) -> dict[str, AliasSourceSummary]:
    """Summarize ModelSEED reaction aliases from the cached alias TSV."""
    rows: Counter[str] = Counter()
    distinct_ids: dict[str, set[str]] = {}

    with open(input_file, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            source = row["Source"].strip()
            modelseed_id = row["ModelSEED ID"].strip()
            rows[source] += 1
            distinct_ids.setdefault(source, set()).add(modelseed_id)

    return {
        source: {"rows": rows[source], "distinct_ids": len(distinct_ids[source])}
        for source in sorted(rows)
    }


def _summarize_reaction_ec_file(input_file: str | Path) -> ReactionECCoverageSummary:
    """Summarize ModelSEED EC coverage from the cached EC TSV."""
    rows = 0
    reaction_ids: set[str] = set()

    with open(input_file, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows += 1
            reaction_ids.add(row["ModelSEED ID"].strip())

    return {"rows": rows, "distinct_reactions": len(reaction_ids)}


def _count_tsv_rows(input_file: str | Path) -> int:
    """Count data rows in a TSV file, excluding the header."""
    with open(input_file, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            next(reader)
        except StopIteration:
            return 0
        return sum(1 for _ in reader)


def summarize_modelseed_cache(
    cache_dir: str | Path = "cache",
) -> ModelSeedCacheSummary:
    """Summarize cached ModelSEED bridge and alias coverage."""
    cache_path = Path(cache_dir)
    modelseed_dir = cache_path / "modelseed"
    mappings_path = cache_path / "modelseed_to_chebi.jsonl"
    reactions_path = cache_path / "modelseed_reactions.jsonl"

    required_paths = [
        modelseed_dir / MODELSEED_FILES["compounds"],
        modelseed_dir / MODELSEED_FILES["reaction_aliases"],
        modelseed_dir / MODELSEED_FILES["reaction_ecs"],
        mappings_path,
        reactions_path,
    ]
    missing = [path for path in required_paths if not path.exists()]
    if missing:
        missing_str = ", ".join(str(path) for path in missing[:3])
        raise FileNotFoundError(
            f"ModelSEED cache files not found ({missing_str}). "
            f"Run 'autarch cache-modelseed --cache-dir {cache_path}' first."
        )

    etl = ModelSeedETL(cache_dir=modelseed_dir)
    compounds = etl.load_compounds()
    mappings = load_modelseed_to_chebi_mappings(mappings_path)
    cached_reactions = load_modelseed_reaction_records(reactions_path)

    compound_mapping_methods = Counter(
        mapping.mapping_method or "unmapped" for mapping in mappings.values()
    )
    compound_alias_sources = _summarize_alias_sources(
        {compound.modelseed_id: compound.aliases for compound in compounds.values()}
    )
    reaction_alias_sources = _summarize_reaction_alias_file(
        modelseed_dir / MODELSEED_FILES["reaction_aliases"]
    )
    ec_coverage = _summarize_reaction_ec_file(
        modelseed_dir / MODELSEED_FILES["reaction_ecs"]
    )
    total_reactions = _count_tsv_rows(modelseed_dir / MODELSEED_FILES["reactions"])
    rhea_alias_summary = reaction_alias_sources.get("rhea")

    return {
        "total_compounds": len(compounds),
        "mapped_compounds": sum(1 for mapping in mappings.values() if mapping.chebi_id),
        "compound_mapping_methods": dict(sorted(compound_mapping_methods.items())),
        "compound_alias_sources": compound_alias_sources,
        "total_reactions": total_reactions,
        "reactions_with_rhea": (
            rhea_alias_summary["distinct_ids"] if rhea_alias_summary is not None else 0
        ),
        "cached_reactions": len(cached_reactions),
        "cached_reactions_with_rhea": sum(
            1 for reaction in cached_reactions.values() if reaction.rhea_ids
        ),
        "cached_fully_mapped_reactions": sum(
            1 for reaction in cached_reactions.values() if reaction.fully_mapped_to_chebi
        ),
        "reaction_alias_sources": reaction_alias_sources,
        "ec_coverage": ec_coverage,
    }


def cache_modelseed_dataset(
    cache_dir: str | Path = "cache",
    limit: Optional[int] = None,
    force_download: bool = False,
) -> dict[str, int]:
    """Cache ModelSEED compounds, reactions, and the ModelSEED→ChEBI bridge."""
    cache_path = Path(cache_dir)
    cache_path.mkdir(parents=True, exist_ok=True)

    etl = ModelSeedETL(cache_dir=cache_path / "modelseed", force_download=force_download)
    etl.download_required_files()

    compounds = etl.load_compounds()
    mappings = etl.build_chebi_bridge(compounds, cache_dir=cache_path)
    reactions = etl.load_reactions(compounds, mappings, limit=limit)

    serialize_modelseed_compounds(compounds, cache_path / "modelseed_compounds.jsonl")
    serialize_modelseed_mappings(mappings, cache_path / "modelseed_to_chebi.jsonl")
    serialize_modelseed_reactions(reactions, cache_path / "modelseed_reactions.jsonl")

    return {
        "total_compounds": len(compounds),
        "mapped_compounds": sum(1 for mapping in mappings.values() if mapping.chebi_id),
        "total_reactions": len(reactions),
        "reactions_with_rhea": sum(1 for reaction in reactions.values() if reaction.rhea_ids),
        "fully_mapped_reactions": sum(
            1 for reaction in reactions.values() if reaction.fully_mapped_to_chebi
        ),
    }
