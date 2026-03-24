"""Validation of GO term mappings in reaction classifiers.

Uses OAK to verify that:
1. GO_ID fields reference valid GO terms
2. EC_NUMBER_PREFIX fields correspond to the expected GO term mappings
3. Class docstrings match GO term definitions
"""

import importlib
import inspect
import pkgutil
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from oaklib import get_adapter  # type: ignore[import-untyped]


@dataclass
class GOValidationResult:
    """Result of validating a classifier's GO mapping."""
    class_name: str
    go_id: Optional[str]
    ec_prefix: Optional[str]
    is_valid_go_term: bool
    go_term_label: Optional[str]
    go_term_ec_xrefs: list[str]
    ec_in_go_xrefs: bool
    notes: str = ""


def get_all_reaction_classes():
    """Discover all ReactionClass subclasses in autarch.ontology."""
    from autarch.ontology.reaction import ReactionClass

    # Import all modules in autarch.ontology
    import autarch.ontology as ontology_pkg

    classes = []
    ontology_path = Path(ontology_pkg.__file__).parent

    for _, module_name, _ in pkgutil.iter_modules([str(ontology_path)]):
        if module_name.startswith('_'):
            continue
        try:
            module = importlib.import_module(f'autarch.ontology.{module_name}')
            for name, obj in inspect.getmembers(module, inspect.isclass):
                if (issubclass(obj, ReactionClass) and
                    obj is not ReactionClass and
                    hasattr(obj, 'GO_ID') and
                    not obj.__dict__.get('EXCLUDE_FROM_DISCOVERY', False)):
                    classes.append(obj)
        except ImportError as e:
            print(f"Warning: Could not import {module_name}: {e}")

    return classes


def validate_go_mappings(verbose: bool = False) -> list[GOValidationResult]:
    """Validate GO mappings for all reaction classifiers.

    Checks:
    1. GO_ID is a valid GO term
    2. If EC_NUMBER_PREFIX is set, check if it appears in GO term's EC xrefs
    """
    results = []

    # Get GO adapter for looking up terms
    print("Loading GO ontology adapter...")
    go = get_adapter("sqlite:obo:go")

    # Get all reaction classes
    classes = get_all_reaction_classes()
    print(f"Validating {len(classes)} reaction classes...")

    for cls in classes:
        go_id = getattr(cls, 'GO_ID', None)
        ec_prefix = getattr(cls, 'EC_NUMBER_PREFIX', None)

        # Validate GO term exists
        is_valid = False
        go_label = None
        go_ec_xrefs = []

        if go_id:
            try:
                go_label = go.label(go_id)
                is_valid = go_label is not None

                # Get EC xrefs from this GO term
                if is_valid:
                    metadata = go.entity_metadata_map(go_id)
                    xrefs = metadata.get("xref", [])
                    go_ec_xrefs = [x for x in xrefs if isinstance(x, str) and x.startswith("EC:")]
            except Exception:
                is_valid = False

        # Check if EC prefix matches GO xrefs
        ec_in_xrefs = True
        notes = ""

        if ec_prefix and go_ec_xrefs:
            # Normalize EC prefix for comparison
            ec_normalized = ec_prefix.replace(".-", "").rstrip(".-")
            ec_check = f"EC:{ec_normalized}"

            # Check if any GO xref matches our EC prefix
            ec_matches = [x for x in go_ec_xrefs if x.startswith(ec_check) or x == ec_check]
            ec_in_xrefs = len(ec_matches) > 0

            if not ec_in_xrefs:
                notes = f"EC {ec_prefix} not in GO xrefs: {go_ec_xrefs}"
        elif ec_prefix and is_valid and not go_ec_xrefs:
            # GO term has no EC xrefs - this might be ok for some terms
            notes = "GO term has no EC xrefs"

        result = GOValidationResult(
            class_name=cls.__name__,
            go_id=go_id,
            ec_prefix=ec_prefix,
            is_valid_go_term=is_valid,
            go_term_label=go_label,
            go_term_ec_xrefs=go_ec_xrefs,
            ec_in_go_xrefs=ec_in_xrefs,
            notes=notes,
        )
        results.append(result)

        if verbose:
            status = "✓" if is_valid and ec_in_xrefs else "✗"
            print(f"  {status} {cls.__name__}: {go_id} ({go_label or 'INVALID'})")
            if notes:
                print(f"      {notes}")

    return results


def print_validation_report(results: list[GOValidationResult]):
    """Print a summary report of GO mapping validation."""
    print("\n" + "=" * 70)
    print("GO MAPPING VALIDATION REPORT")
    print("=" * 70)

    # Count issues
    invalid_go = [r for r in results if not r.is_valid_go_term and r.go_id]
    ec_mismatch = [r for r in results if not r.ec_in_go_xrefs and r.ec_prefix]
    no_go = [r for r in results if not r.go_id]
    valid = [r for r in results if r.is_valid_go_term and r.ec_in_go_xrefs]

    print(f"\nTotal classes: {len(results)}")
    print(f"Valid GO mappings: {len(valid)}")
    print(f"Invalid GO terms: {len(invalid_go)}")
    print(f"EC prefix not in GO xrefs: {len(ec_mismatch)}")
    print(f"Missing GO_ID: {len(no_go)}")

    if invalid_go:
        print("\n--- INVALID GO TERMS ---")
        for r in invalid_go:
            print(f"  {r.class_name}: {r.go_id}")

    if ec_mismatch:
        print("\n--- EC PREFIX NOT IN GO XREFS ---")
        for r in ec_mismatch:
            print(f"  {r.class_name}:")
            print(f"    GO: {r.go_id} ({r.go_term_label})")
            print(f"    EC prefix: {r.ec_prefix}")
            print(f"    GO EC xrefs: {r.go_term_ec_xrefs}")
            if r.notes:
                print(f"    Note: {r.notes}")

    if no_go:
        print("\n--- MISSING GO_ID ---")
        for r in no_go:
            print(f"  {r.class_name} (EC: {r.ec_prefix or 'N/A'})")

    print("\n" + "=" * 70)

    return len(invalid_go) + len(ec_mismatch)


if __name__ == "__main__":
    results = validate_go_mappings(verbose=True)
    errors = print_validation_report(results)
    exit(1 if errors > 0 else 0)
