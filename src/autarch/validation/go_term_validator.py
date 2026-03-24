"""Strict GO term validation for reaction classifier naming consistency.

This validator enforces concept alignment across:
- GO label from GO_ID, or EC label when GO_ID is absent
- classifier class name
- ontology module filename
- classifier class docstring (first line)

It also performs semantic consistency checks when GO metadata are available:
- EC mapping consistency between classifier EC metadata and GO xrefs
- GO hierarchy consistency with classifier inheritance
- Detection of mapping-driven source-code exclusions that can overfit evals
- Basic GO-definition vs source-logic conflict checks

Normalization allows:
- CamelCase vs snake_case differences
- omission of the non-informative token "activity"
"""

from dataclasses import dataclass
from inspect import getdoc
import inspect
import json
from pathlib import Path
import re
from typing import Any, Optional, Type, cast

from oaklib import get_adapter  # type: ignore[import-untyped]

from autarch.ontology.reaction import ReactionClass

DEFAULT_BASELINE_PATH = Path(__file__).with_name("go_term_validator_baseline.json")
FAILURE_FLAG_FIELDS = (
    "file_matches_class",
    "go_matches_class",
    "description_matches_class",
    "ec_matches_go",
    "inheritance_matches_go",
    "source_matches_go",
    "definition_matches_go",
)


NON_INFORMATIVE_TOKENS = {"activity"}
NUMERIC_TOKEN_NORMALIZATION = {
    "zero": "0",
    "one": "1",
    "two": "2",
    "three": "3",
    "four": "4",
    "five": "5",
    "six": "6",
    "seven": "7",
    "eight": "8",
    "nine": "9",
    "ten": "10",
    "eleven": "11",
    "twelve": "12",
    "thirteen": "13",
    "fourteen": "14",
    "fifteen": "15",
    "sixteen": "16",
    "seventeen": "17",
    "eighteen": "18",
    "nineteen": "19",
    "twenty": "20",
    "sulphur": "sulfur",
}
GO_MAPPING_ANTI_PATTERN_RE = re.compile(r"NON_[0-9]{4,}_[A-Z0-9_]*EXCLUSION", re.IGNORECASE)
GO_MAPPING_TEXT_ANTI_PATTERNS = (
    "outside go:",
    "current evaluation graph",
    "mapped outside go:",
)
TEXT_MINING_ANTI_PATTERNS = (
    (re.compile(r"\breaction\.label\b"), "uses reaction.label text parsing"),
    (re.compile(r"\bparticipant\.name\b"), "uses participant.name text parsing"),
    (re.compile(r"\b[A-Z0-9_]*LABEL[A-Z0-9_]*\b"), "defines label-based keyword lists"),
    (re.compile(r"\b[A-Z0-9_]*NAME_[A-Z0-9_]*\b"), "defines name-based keyword lists"),
)


@dataclass
class GoTermMetadata:
    """GO metadata required for semantic validation."""

    go_id: str
    label: str
    definition: str
    namespace: str
    ec_numbers: list[str]
    ancestors: list[str]


def _split_camel_case(value: str) -> str:
    """Split CamelCase boundaries to spaces."""
    value = re.sub(r"([A-Z]+)([A-Z][a-z])", r"\1 \2", value)
    return re.sub(r"([a-z0-9])([A-Z])", r"\1 \2", value)


def _normalize_phrase(value: str) -> str:
    """Normalize a phrase for strict concept matching."""
    if not value:
        return ""

    expanded = _split_camel_case(value).replace("_", " ").replace("-", " ")
    expanded = re.sub(r"NAD\(P\)\+", "NADOrNADP", expanded, flags=re.IGNORECASE)
    expanded = re.sub(r"NAD\(P\)H", "NADHOrNADPH", expanded, flags=re.IGNORECASE)
    tokens = re.findall(r"[A-Za-z0-9]+", expanded.lower())
    tokens = [NUMERIC_TOKEN_NORMALIZATION.get(token, token) for token in tokens]
    tokens = [token for token in tokens if token not in NON_INFORMATIVE_TOKENS]
    return "".join(tokens)


def _normalize_ec_label_fragment(value: str) -> str:
    """Normalize EC label fragments into a concept phrase."""
    cleaned = value.strip()
    if ". " in cleaned:
        cleaned = cleaned.split(". ", 1)[0].strip()
    lowered = cleaned.lower()
    if lowered.startswith("catalysing the "):
        cleaned = cleaned[len("Catalysing the "):]
        lowered = cleaned.lower()
    if lowered.startswith("catalyzing the "):
        cleaned = cleaned[len("Catalyzing the "):]
        lowered = cleaned.lower()
    if lowered.startswith("translocation ") and not lowered.startswith("translocation of "):
        cleaned = "translocation of " + cleaned[len("translocation "):]
    return cleaned


def _extract_description_phrase(cls: Type[ReactionClass]) -> str:
    """Extract the class description phrase from the leading docstring header."""
    doc = getdoc(cls) or ""
    if not doc:
        return ""

    lines = [line.strip() for line in doc.strip().splitlines() if line.strip()]
    if not lines:
        return ""

    header_lines = [lines[0]]
    for line in lines[1:]:
        if header_lines[-1].endswith((".", "!", "?", ":", ";")):
            break
        if not re.match(r"^[a-z(\[]", line):
            break
        header_lines.append(line)

    first_line = " ".join(header_lines).strip()
    cleaned = re.sub(r"^(classifier|class)\s+for\s+", "", first_line, flags=re.IGNORECASE)
    cleaned = re.sub(r"\s+reaction\s+classification\.?$", "", cleaned, flags=re.IGNORECASE)
    cleaned = re.sub(r"\s+reactions?\.?$", "", cleaned, flags=re.IGNORECASE)
    return cleaned.strip()


def _module_stem_from_class(cls: Type[ReactionClass]) -> str:
    """Get module stem (filename without .py) from class module."""
    return cls.__module__.split(".")[-1]


def _safe_get_class_source(cls: Type[ReactionClass]) -> str:
    """Get class source code, returning empty string if unavailable."""
    try:
        return inspect.getsource(cls)
    except Exception:
        return ""


def _normalize_ec(ec: str) -> str:
    """Normalize an EC string to 4 components with '-' wildcards."""
    cleaned = ec.strip().replace("EC:", "")
    parts = cleaned.split(".")
    while len(parts) < 4:
        parts.append("-")
    return ".".join(parts[:4])


def _ec_specificity(ec: str) -> int:
    """Count non-wildcard levels in an EC code."""
    return sum(1 for part in _normalize_ec(ec).split(".") if part != "-")


def _ec_matches_prefix(ec_number: str, ec_prefix: str) -> bool:
    """Check if an EC number matches a wildcard prefix pattern."""
    ec_number_parts = _normalize_ec(ec_number).split(".")
    ec_prefix_parts = _normalize_ec(ec_prefix).split(".")
    for number_part, prefix_part in zip(ec_number_parts, ec_prefix_parts):
        if prefix_part == "-":
            return True
        if number_part != prefix_part:
            return False
    return True


def _ec_query_id(ec_prefix: str) -> str:
    """Convert an EC prefix into the most specific queryable EC CURIE."""
    normalized = _normalize_ec(ec_prefix)
    parts = normalized.split(".")
    while parts and parts[-1] == "-":
        parts.pop()
    return f"EC:{'.'.join(parts)}" if parts else "EC:"


def _ec_parent_id(ec_id: str) -> Optional[str]:
    """Return the parent EC CURIE for a queryable EC CURIE."""
    cleaned = ec_id.replace("EC:", "")
    parts = cleaned.split(".")
    if len(parts) <= 1:
        return None
    return f"EC:{'.'.join(parts[:-1])}"


def _extract_ec_xrefs(values: list[str]) -> list[str]:
    """Extract EC numbers from GO xref strings."""
    ec_values: list[str] = []
    for value in values:
        if value.startswith("EC:"):
            ec_values.append(value.replace("EC:", ""))
            continue
        # Example: <http://purl.uniprot.org/enzyme/1.1.1.1>
        match = re.search(r"/enzyme/([0-9\-]+\.[0-9\-]+\.[0-9\-]+\.[0-9\-]+)", value)
        if match:
            ec_values.append(match.group(1))
    return sorted(set(ec_values))


def _get_all_reaction_classes() -> list[Type[ReactionClass]]:
    """Discover all concrete ReactionClass subclasses."""
    # Import side effects from ontology __init__ register all classes.
    import autarch.ontology  # noqa: F401

    def get_all_subclasses(base: type[Any]) -> set[Type[ReactionClass]]:
        subclasses: set[Type[ReactionClass]] = set()
        for subclass in base.__subclasses__():
            if not issubclass(subclass, ReactionClass):
                continue
            reaction_subclass = cast(Type[ReactionClass], subclass)
            subclasses.add(reaction_subclass)
            subclasses.update(get_all_subclasses(reaction_subclass))
        return subclasses

    all_subclasses = get_all_subclasses(ReactionClass)
    concrete_classes = [
        cls
        for cls in all_subclasses
        if not getattr(cls, "__abstractmethods__", None)
        if not cls.__dict__.get("EXCLUDE_FROM_DISCOVERY", False)
    ]
    return sorted(concrete_classes, key=lambda cls: cls.__name__)


@dataclass
class ValidationResult:
    """Result of validating a single reaction class."""

    class_name: str
    go_id: str
    go_name: str
    module_stem: str
    description: str
    go_definition: str
    go_namespace: str
    go_ec_numbers: list[str]
    normalized_class: str
    normalized_go: str
    normalized_module: str
    normalized_description: str
    file_matches_class: bool
    go_matches_class: bool
    description_matches_class: bool
    ec_matches_go: bool
    inheritance_matches_go: bool
    source_matches_go: bool
    definition_matches_go: bool
    matches: bool
    explanation: str

    @property
    def failure_flags(self) -> list[str]:
        """Return the names of failed validation checks."""
        return [
            field_name
            for field_name in FAILURE_FLAG_FIELDS
            if not getattr(self, field_name)
        ]


class GoTermValidator:
    """Validator for strict GO/class/file/description concept alignment."""

    def __init__(
        self,
        adapter: Optional[Any] = None,
        ec_adapter: Optional[Any] = None,
        go_terms_cache_path: str = "cache/go_terms.jsonl",
        go_term_map: Optional[dict[str, dict[str, Any] | GoTermMetadata]] = None,
    ):
        """Initialize with GO ontology adapter and optional GO metadata cache."""
        self.adapter = adapter if adapter is not None else get_adapter("sqlite:obo:go")
        self.ec_adapter = (
            ec_adapter
            if ec_adapter is not None
            else adapter
            if adapter is not None
            else get_adapter("sqlite:obo:ec")
        )
        if go_term_map is not None:
            self.go_term_map = {
                go_id: self._coerce_go_term_metadata(go_id, value)
                for go_id, value in go_term_map.items()
            }
        else:
            self.go_term_map = self._load_go_term_cache(go_terms_cache_path)

    def _coerce_go_term_metadata(
        self, go_id: str, value: dict[str, Any] | GoTermMetadata
    ) -> GoTermMetadata:
        """Coerce cached GO metadata into a typed object."""
        if isinstance(value, GoTermMetadata):
            return value
        return GoTermMetadata(
            go_id=go_id,
            label=value.get("label", ""),
            definition=value.get("definition", ""),
            namespace=value.get("namespace", ""),
            ec_numbers=list(value.get("ec_numbers", [])),
            ancestors=list(value.get("ancestors", [])),
        )

    def _load_go_term_cache(self, cache_path: str) -> dict[str, GoTermMetadata]:
        """Load GO metadata from cache/go_terms.jsonl when available."""
        path = Path(cache_path)
        if not path.exists():
            return {}
        metadata_map: dict[str, GoTermMetadata] = {}
        with open(path) as handle:
            for line in handle:
                raw = json.loads(line)
                go_id = raw.get("go_id")
                if not go_id:
                    continue
                metadata_map[go_id] = GoTermMetadata(
                    go_id=go_id,
                    label=raw.get("label", ""),
                    definition=raw.get("definition", ""),
                    namespace=raw.get("namespace", ""),
                    ec_numbers=list(raw.get("ec_numbers", [])),
                    ancestors=list(raw.get("ancestors", [])),
                )
        return metadata_map

    def _get_go_term_metadata(self, go_id: str) -> Optional[GoTermMetadata]:
        """Get GO metadata from cache or adapter."""
        if go_id in self.go_term_map:
            return self.go_term_map[go_id]

        # Fallback: query adapter directly for metadata.
        try:
            label = self.adapter.label(go_id) or ""
        except Exception:
            return None

        definition = ""
        namespace = ""
        ec_numbers: list[str] = []
        ancestors: list[str] = []

        try:
            metadata = self.adapter.entity_metadata_map(go_id)
            definitions = metadata.get("IAO:0000115", [])
            if isinstance(definitions, list) and definitions:
                definition = str(definitions[0])
            elif isinstance(definitions, str):
                definition = definitions

            namespaces = metadata.get("oio:hasOBONamespace", [])
            if isinstance(namespaces, list) and namespaces:
                namespace = str(namespaces[0])
            elif isinstance(namespaces, str):
                namespace = namespaces

            dbxrefs = metadata.get("oio:hasDbXref", [])
            if isinstance(dbxrefs, str):
                dbxrefs = [dbxrefs]
            ec_numbers = _extract_ec_xrefs([x for x in dbxrefs if isinstance(x, str)])
        except Exception:
            pass

        try:
            ancestors = list(self.adapter.ancestors(go_id, predicates=["rdfs:subClassOf"]))
        except Exception:
            ancestors = []

        term = GoTermMetadata(
            go_id=go_id,
            label=label,
            definition=definition,
            namespace=namespace,
            ec_numbers=ec_numbers,
            ancestors=ancestors,
        )
        self.go_term_map[go_id] = term
        return term

    def get_go_term_label(self, go_id: str) -> Optional[str]:
        """Get GO term label using OAK API."""
        metadata = self._get_go_term_metadata(go_id)
        if metadata and metadata.label:
            return metadata.label
        try:
            return self.adapter.label(go_id)
        except Exception as e:
            print(f"Error getting GO term info for {go_id}: {e}")
            return None

    def get_ec_term_label(self, ec_id: str) -> Optional[str]:
        """Get EC term label using the EC adapter."""
        try:
            label = self.ec_adapter.label(ec_id)
            if not label:
                return None

            normalized_label = _normalize_ec_label_fragment(label)
            lowered_label = normalized_label.lower()
            if not lowered_label.startswith(("linked to ", "with ", "using ", "in ")) and lowered_label != "miscellaneous":
                return normalized_label

            parent_id = _ec_parent_id(ec_id)
            if not parent_id:
                return normalized_label

            parent_label = self.ec_adapter.label(parent_id)
            if not parent_label:
                return normalized_label

            normalized_parent = _normalize_ec_label_fragment(parent_label)
            return f"{normalized_parent} {normalized_label.lower()}"
        except Exception as e:
            print(f"Error getting EC term info for {ec_id}: {e}")
            return None

    def _check_ec_consistency(
        self,
        class_ec_prefix: str,
        go_ec_numbers: list[str],
    ) -> tuple[bool, str]:
        """Validate class EC metadata against GO EC xrefs."""
        if not class_ec_prefix or not go_ec_numbers:
            return True, ""

        normalized_class_ec = _normalize_ec(class_ec_prefix)
        normalized_go_ecs = [_normalize_ec(ec) for ec in go_ec_numbers]

        matches_any = any(
            _ec_matches_prefix(go_ec, normalized_class_ec)
            for go_ec in normalized_go_ecs
        )
        if not matches_any:
            return (
                False,
                f"class EC '{class_ec_prefix}' not compatible with GO EC xrefs {go_ec_numbers}",
            )

        # If GO has exactly one fully specified EC xref, reject broader classifier EC prefixes.
        if len(normalized_go_ecs) == 1:
            go_ec = normalized_go_ecs[0]
            if _ec_specificity(go_ec) == 4 and _ec_specificity(normalized_class_ec) < 4:
                return (
                    False,
                    f"class EC '{class_ec_prefix}' is broader than exact GO EC '{go_ec}'",
                )

        return True, ""

    def _check_inheritance_consistency(
        self,
        cls: Type[ReactionClass],
        go_ancestors: list[str],
    ) -> tuple[bool, str]:
        """Check that direct classifier parents are compatible in GO hierarchy."""
        if not go_ancestors:
            return True, ""

        parent_go_ids = []
        for base in cls.__bases__:
            if not issubclass(base, ReactionClass) or base is ReactionClass:
                continue
            parent_go_id = getattr(base, "GO_ID", None)
            if parent_go_id:
                parent_go_ids.append(parent_go_id)

        missing = [go_id for go_id in parent_go_ids if go_id not in set(go_ancestors)]
        if missing:
            return (
                False,
                f"parent GO_ID(s) {missing} not found in child GO ancestors",
            )
        return True, ""

    def _check_source_consistency(self, source_text: str) -> tuple[bool, str]:
        """Detect source anti-patterns (mapping-driven logic and text-mining heuristics)."""
        if not source_text:
            return True, ""

        issues = []
        if GO_MAPPING_ANTI_PATTERN_RE.search(source_text):
            issues.append("contains GO-ID keyed exclusion list (mapping-driven logic)")

        lower_source = source_text.lower()
        for pattern in GO_MAPPING_TEXT_ANTI_PATTERNS:
            if pattern in lower_source:
                issues.append(f"contains mapping-driven phrase '{pattern}'")

        for regex_pattern, message in TEXT_MINING_ANTI_PATTERNS:
            if regex_pattern.search(source_text):
                issues.append(message)

        if issues:
            return False, "; ".join(issues)
        return True, ""

    def _check_definition_consistency(
        self,
        go_definition: str,
        source_text: str,
    ) -> tuple[bool, str]:
        """Check for direct conflicts between GO definition semantics and source logic."""
        if not go_definition or not source_text:
            return True, ""

        definition_lower = go_definition.lower()
        source_lower = source_text.lower()
        issues = []

        # GO: "... aldehyde or ketone ..." should not be gated to narrow subcontexts.
        if "aldehyde or ketone" in definition_lower and "sugar_polyol_context" in source_lower:
            issues.append("ketone branch is context-gated despite GO definition allowing generic aldehyde/ketone")

        if issues:
            return False, "; ".join(issues)
        return True, ""

    def validate_class(self, cls: Type[ReactionClass]) -> ValidationResult:
        """Validate a single reaction class."""
        class_name = cls.__name__
        go_id = cls.GO_ID
        class_ec_prefix = getattr(cls, "EC_NUMBER_PREFIX", "")
        module_stem = _module_stem_from_class(cls)
        description = _extract_description_phrase(cls)
        source_text = _safe_get_class_source(cls)

        normalized_class = _normalize_phrase(class_name)
        normalized_module = _normalize_phrase(module_stem)
        normalized_description = _normalize_phrase(description)

        # Checks that do not require GO lookup
        file_matches_class = normalized_module == normalized_class
        description_matches_class = normalized_description == normalized_class

        if not go_id:
            if not class_ec_prefix:
                issues = ["missing GO_ID", "missing EC_NUMBER_PREFIX"]
                if not file_matches_class:
                    issues.append("filename does not match class name")
                if not description_matches_class:
                    issues.append("description does not match class name")
                return ValidationResult(
                    class_name=class_name,
                    go_id="None",
                    go_name="No GO ID specified",
                    module_stem=module_stem,
                    description=description,
                    go_definition="",
                    go_namespace="",
                    go_ec_numbers=[],
                    normalized_class=normalized_class,
                    normalized_go="",
                    normalized_module=normalized_module,
                    normalized_description=normalized_description,
                    file_matches_class=file_matches_class,
                    go_matches_class=False,
                    description_matches_class=description_matches_class,
                    ec_matches_go=True,
                    inheritance_matches_go=True,
                    source_matches_go=True,
                    definition_matches_go=True,
                    matches=False,
                    explanation="❌ " + "; ".join(issues),
                )

            ec_query_id = _ec_query_id(class_ec_prefix)
            ec_label = self.get_ec_term_label(ec_query_id) if ec_query_id else None
            normalized_ec = _normalize_phrase(_normalize_ec_label_fragment(ec_label or ""))
            label_matches_class = normalized_ec == normalized_class
            source_matches_go, source_issue = self._check_source_consistency(source_text)
            issues = []
            if not file_matches_class:
                issues.append("filename does not match class name")
            if not label_matches_class:
                issues.append(f"EC label '{ec_label or 'Could not retrieve'}' != class '{class_name}'")
            if not description_matches_class:
                issues.append("description does not match class name")
            if not source_matches_go:
                issues.append(source_issue)
            return ValidationResult(
                class_name=class_name,
                go_id=ec_query_id or "None",
                go_name=ec_label or "Could not retrieve",
                module_stem=module_stem,
                description=description,
                go_definition="",
                go_namespace="ec",
                go_ec_numbers=[_normalize_ec(class_ec_prefix)],
                normalized_class=normalized_class,
                normalized_go=normalized_ec,
                normalized_module=normalized_module,
                normalized_description=normalized_description,
                file_matches_class=file_matches_class,
                go_matches_class=label_matches_class,
                description_matches_class=description_matches_class,
                ec_matches_go=True,
                inheritance_matches_go=True,
                source_matches_go=source_matches_go,
                definition_matches_go=True,
                matches=len(issues) == 0,
                explanation="✅ Strict concept alignment" if not issues else "❌ " + "; ".join(issues),
            )

        go_metadata = self._get_go_term_metadata(go_id)
        go_name = go_metadata.label if go_metadata and go_metadata.label else self.get_go_term_label(go_id)
        if not go_name:
            issues = [f"could not retrieve GO term {go_id}"]
            if not file_matches_class:
                issues.append("filename does not match class name")
            if not description_matches_class:
                issues.append("description does not match class name")
            return ValidationResult(
                class_name=class_name,
                go_id=go_id,
                go_name="Could not retrieve",
                module_stem=module_stem,
                description=description,
                go_definition="",
                go_namespace="",
                go_ec_numbers=[],
                normalized_class=normalized_class,
                normalized_go="",
                normalized_module=normalized_module,
                normalized_description=normalized_description,
                file_matches_class=file_matches_class,
                go_matches_class=False,
                description_matches_class=description_matches_class,
                ec_matches_go=True,
                inheritance_matches_go=True,
                source_matches_go=True,
                definition_matches_go=True,
                matches=False,
                explanation="❌ " + "; ".join(issues),
            )

        go_definition = go_metadata.definition if go_metadata else ""
        go_namespace = go_metadata.namespace if go_metadata else ""
        go_ec_numbers = go_metadata.ec_numbers if go_metadata else []
        go_ancestors = go_metadata.ancestors if go_metadata else []

        normalized_go = _normalize_phrase(go_name)
        go_matches_class = normalized_go == normalized_class
        ec_matches_go, ec_issue = self._check_ec_consistency(class_ec_prefix, go_ec_numbers)
        inheritance_matches_go, inheritance_issue = self._check_inheritance_consistency(
            cls, go_ancestors
        )
        source_matches_go, source_issue = self._check_source_consistency(source_text)
        definition_matches_go, definition_issue = self._check_definition_consistency(
            go_definition, source_text
        )

        issues = []
        if not file_matches_class:
            issues.append(
                f"filename '{module_stem}' != class '{class_name}'"
            )
        if not go_matches_class:
            issues.append(
                f"GO label '{go_name}' != class '{class_name}'"
            )
        if not description_matches_class:
            issues.append(
                f"description '{description}' != class '{class_name}'"
            )
        if not ec_matches_go:
            issues.append(ec_issue)
        if not inheritance_matches_go:
            issues.append(inheritance_issue)
        if not source_matches_go:
            issues.append(source_issue)
        if not definition_matches_go:
            issues.append(definition_issue)

        matches = len(issues) == 0
        explanation = "✅ Strict concept alignment" if matches else "❌ " + "; ".join(issues)

        return ValidationResult(
            class_name=class_name,
            go_id=go_id,
            go_name=go_name,
            module_stem=module_stem,
            description=description,
            go_definition=go_definition,
            go_namespace=go_namespace,
            go_ec_numbers=go_ec_numbers,
            normalized_class=normalized_class,
            normalized_go=normalized_go,
            normalized_module=normalized_module,
            normalized_description=normalized_description,
            file_matches_class=file_matches_class,
            go_matches_class=go_matches_class,
            description_matches_class=description_matches_class,
            ec_matches_go=ec_matches_go,
            inheritance_matches_go=inheritance_matches_go,
            source_matches_go=source_matches_go,
            definition_matches_go=definition_matches_go,
            matches=matches,
            explanation=explanation,
        )

    def validate_all_classes(self) -> list[ValidationResult]:
        """Validate all concrete reaction classes."""
        classes = _get_all_reaction_classes()
        return [self.validate_class(cls) for cls in classes]

    def print_validation_report(self, results: list[ValidationResult]) -> bool:
        """Print a formatted validation report."""
        print("🔍 Validating strict GO/class/file/description + semantic consistency...\n")

        issues = [result for result in results if not result.matches]

        for result in results:
            status = "✅" if result.matches else "❌"
            print(f"{status} {result.class_name} -> {result.go_id}")
            print(f"   GO name: {result.go_name}")
            if result.go_namespace:
                print(f"   GO namespace: {result.go_namespace}")
            if result.go_ec_numbers:
                print(f"   GO EC xrefs: {', '.join(result.go_ec_numbers)}")
            print(f"   filename: {result.module_stem}.py")
            print(f"   description: {result.description or '(missing)'}")
            print(f"   {result.explanation}")
            print()

        print(f"📊 Summary: Checked {len(results)} classes")
        print(f"✅ Strictly aligned: {len(results) - len(issues)}")
        print(f"❌ Inconsistencies: {len(issues)}")

        return len(issues) == 0

    def validate(self) -> bool:
        """Run full validation and return success status."""
        results = self.validate_all_classes()
        return self.print_validation_report(results)

    def load_baseline(
        self,
        baseline_path: Path | str = DEFAULT_BASELINE_PATH,
    ) -> dict[str, list[str]]:
        """Load a baseline of known validation failures."""
        path = Path(baseline_path)
        if not path.exists():
            return {}
        with open(path) as handle:
            raw = json.load(handle)
        if not isinstance(raw, dict):
            return {}
        baseline: dict[str, list[str]] = {}
        for class_name, flags in raw.items():
            if isinstance(class_name, str) and isinstance(flags, list):
                baseline[class_name] = [
                    flag for flag in flags if flag in FAILURE_FLAG_FIELDS
                ]
        return baseline

    def build_baseline_map(self, results: list[ValidationResult]) -> dict[str, list[str]]:
        """Build a baseline map from validation results."""
        return {
            result.class_name: result.failure_flags
            for result in results
            if result.failure_flags
        }

    def write_baseline(
        self,
        baseline_path: Path | str = DEFAULT_BASELINE_PATH,
    ) -> dict[str, list[str]]:
        """Write the current failure baseline to disk."""
        results = self.validate_all_classes()
        baseline = self.build_baseline_map(results)
        path = Path(baseline_path)
        with open(path, "w") as handle:
            json.dump(baseline, handle, indent=2, sort_keys=True)
            handle.write("\n")
        return baseline

    def compare_against_baseline(
        self,
        results: list[ValidationResult],
        baseline: dict[str, list[str]],
    ) -> tuple[list[ValidationResult], list[str]]:
        """Find unexpected failures and resolved baseline entries."""
        unexpected: list[ValidationResult] = []
        current = self.build_baseline_map(results)
        for result in results:
            baseline_flags = set(baseline.get(result.class_name, []))
            current_flags = set(result.failure_flags)
            if current_flags - baseline_flags:
                unexpected.append(result)

        resolved = [
            class_name
            for class_name, flags in baseline.items()
            if class_name not in current or set(current[class_name]) < set(flags)
        ]
        return unexpected, sorted(resolved)

    def print_regression_report(
        self,
        results: list[ValidationResult],
        baseline: dict[str, list[str]],
    ) -> bool:
        """Print validation report relative to a known baseline."""
        print("🔍 Validating GO/class/file/description + semantic consistency against baseline...\n")

        unexpected, resolved = self.compare_against_baseline(results, baseline)

        if unexpected:
            print("❌ New validation regressions:\n")
            for result in unexpected:
                baseline_flags = set(baseline.get(result.class_name, []))
                new_flags = [
                    flag for flag in result.failure_flags if flag not in baseline_flags
                ]
                print(f"❌ {result.class_name} -> {result.go_id}")
                print(f"   GO name: {result.go_name}")
                print(f"   filename: {result.module_stem}.py")
                print(f"   description: {result.description or '(missing)'}")
                print(f"   New failures: {', '.join(new_flags)}")
                print(f"   {result.explanation}")
                print()
        else:
            print("✅ No new validation regressions.\n")

        if resolved:
            print("✅ Improved relative to baseline:")
            for class_name in resolved:
                print(f"   - {class_name}")
            print()

        print(f"📊 Summary: Checked {len(results)} classes")
        print(f"✅ Baseline regressions: {0 if not unexpected else len(unexpected)}")
        print(f"✅ Resolved relative to baseline: {len(resolved)}")
        print(f"📌 Baseline known issues: {len(baseline)}")

        return len(unexpected) == 0

    def validate_against_baseline(
        self,
        baseline_path: Path | str = DEFAULT_BASELINE_PATH,
    ) -> bool:
        """Run validation and fail only on new issues relative to baseline."""
        results = self.validate_all_classes()
        baseline = self.load_baseline(baseline_path)
        return self.print_regression_report(results, baseline)
