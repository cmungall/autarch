"""CLI interface for autarch."""

import json
import logging
from pathlib import Path
from typing import Optional, List

import typer
from rich.console import Console
from rich.table import Table
from typing_extensions import Annotated

from autarch.classifier import ReactionClassifier
from autarch.io import load_reaction_from_yaml
from autarch.display import display_classification_results, display_available_classes
from autarch.etl.chebi_etl import fetch_go_enzyme_mappings, serialize_go_terms
from autarch.etl.rhea_etl import serialize_rhea_reactions
from autarch.etl.chebi_rdkit_cache import save_cache as save_rdkit_cache
from autarch.evaluation import evaluate_reaction_class
from autarch.validation import (
    ClassSummarizer,
    GoTermValidator,
    print_go_cache_validation_report,
    validate_go_cache,
)
from autarch.export import export_html
from autarch.classification_research import (
    all_evidence_valid,
    apply_research_record,
    artifact_paths,
    build_research_prompt,
    collect_reaction_class_metadata,
    extract_research_record,
    load_research_record,
    parse_key_value_options,
    run_deep_research,
    validate_research_record,
    write_research_artifacts,
    write_validation_report,
)

# Set up logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)

# Cache directory
CACHE_DIR = Path("cache")

app = typer.Typer(help="autarch: Agentic Reaction Classification")
console = Console()


@app.command()
def classify(
    rhea_ids: Annotated[
        Optional[List[str]],
        typer.Argument(
            help=(
                "RHEA IDs to classify (e.g., 10000 10001 or RHEA:10000). "
                "For backward compatibility, a single YAML file path is also accepted."
            ),
        ),
    ] = None,
    reaction_class: Annotated[
        Optional[str],
        typer.Option(
            "--class",
            "-c",
            help="Specific reaction class to test (e.g., 'Hydrolase'). If not specified, tests all classes.",
        ),
    ] = None,
    yaml_file: Annotated[
        Optional[Path],
        typer.Option(
            "--yaml",
            "-y",
            help="Path to YAML file containing reaction definition (alternative to RHEA IDs)",
            exists=True,
            readable=True,
        ),
    ] = None,
    cache_dir: Annotated[
        str, typer.Option("--cache-dir", "-d", help="Cache directory path")
    ] = "cache",
    output_format: Annotated[
        str, typer.Option("--format", "-f", help="Output format: table, json, or yaml")
    ] = "table",
    show_explanation: Annotated[
        bool,
        typer.Option(
            "--explanation",
            "-e",
            "--verbose",
            "-v",
            help="Show classification explanations",
        ),
    ] = False,
    show_positive_only: Annotated[
        bool, typer.Option("--positive-only", "-p", help="Only show positive classifications")
    ] = False,
) -> None:
    """Classify chemical reactions by RHEA ID.

    The most common usage is to pass one or more RHEA IDs directly.
    Reactions are loaded from the RHEA cache (run 'autarch cache-rhea' first).

    Examples:
        autarch classify 10000                    # Classify single reaction
        autarch classify 10000 10001 10002        # Classify multiple reactions
        autarch classify RHEA:10000               # RHEA: prefix is optional
        autarch classify 10000 -c Hydrolase       # Test specific class only
        autarch classify 10000 -e                 # Show explanations
        autarch classify 10000 --positive-only    # Only show matches
        autarch classify --yaml reaction.yaml     # From YAML file (advanced)
    """
    try:
        # Initialize classifier
        classifier = ReactionClassifier()

        # Validate reaction class if specified
        if reaction_class and reaction_class not in classifier.reaction_classes:
            typer.echo(f"Unknown reaction class '{reaction_class}'")
            typer.echo(
                f"Available classes: {', '.join(sorted(classifier.reaction_classes.keys()))}"
            )
            raise typer.Exit(1)

        reactions_to_classify = []

        # Backward compatibility: allow `classify reaction.yaml` as positional input.
        if not yaml_file and rhea_ids and len(rhea_ids) == 1:
            candidate = Path(rhea_ids[0])
            if (
                candidate.exists()
                and candidate.is_file()
                and candidate.suffix.lower() in {".yaml", ".yml"}
            ):
                yaml_file = candidate
                rhea_ids = None

        if yaml_file and rhea_ids:
            console.print("[red]Error: Provide either RHEA IDs or --yaml, not both[/red]")
            raise typer.Exit(1)

        # Load from YAML if specified
        if yaml_file:
            reaction = load_reaction_from_yaml(yaml_file)
            reactions_to_classify.append((yaml_file.name, reaction))
        elif rhea_ids:
            # Load from RHEA cache
            cache_path = Path(cache_dir)
            rhea_cache = cache_path / "rhea_reactions.jsonl"

            if not rhea_cache.exists():
                console.print(f"[red]Error: RHEA cache not found at {rhea_cache}[/red]")
                console.print("[yellow]Run 'autarch cache-rhea' first to build the cache[/yellow]")
                raise typer.Exit(1)

            # Load all cached reactions
            rhea_reactions = {}
            with open(rhea_cache) as f:
                for line in f:
                    reaction_data = json.loads(line)
                    rhea_reactions[reaction_data["rhea_id"]] = reaction_data

            # Clean and look up each RHEA ID
            for rid in rhea_ids:
                # Normalize ID (remove RHEA: prefix if present)
                clean_id = rid.replace("RHEA:", "").strip()
                rhea_id = f"RHEA:{clean_id}"

                if rhea_id not in rhea_reactions:
                    console.print(f"[yellow]Warning: {rhea_id} not found in cache[/yellow]")
                    continue

                # Convert cached data to Reaction object
                from autarch.datamodel import Reaction, Participant

                reaction_data = rhea_reactions[rhea_id]
                reaction_dict = reaction_data.get("reaction", {})

                if not reaction_dict:
                    console.print(f"[yellow]Warning: {rhea_id} has no reaction data[/yellow]")
                    continue

                # Build Reaction object
                left_participants = [
                    Participant(**p) for p in reaction_dict.get("left_participants", [])
                ]
                right_participants = [
                    Participant(**p) for p in reaction_dict.get("right_participants", [])
                ]

                reaction = Reaction(
                    left_participants=left_participants,
                    right_participants=right_participants,
                    label=reaction_data.get("label", ""),
                )

                reactions_to_classify.append((rhea_id, reaction))
        else:
            console.print("[red]Error: Please provide RHEA IDs or use --yaml[/red]")
            console.print("\nUsage examples:")
            console.print("  autarch classify 10000")
            console.print("  autarch classify 10000 10001 10002")
            console.print("  autarch classify --yaml reaction.yaml")
            raise typer.Exit(1)

        if not reactions_to_classify:
            console.print("[red]Error: No valid reactions to classify[/red]")
            raise typer.Exit(1)

        # Classify each reaction
        for reaction_id, reaction in reactions_to_classify:
            if len(reactions_to_classify) > 1:
                console.print(f"\n[bold cyan]═══ {reaction_id} ═══[/bold cyan]")
                if reaction.label:
                    console.print(f"[dim]{reaction.label}[/dim]")

            # Perform classification
            if reaction_class:
                result = classifier.classify_with_class(
                    reaction, classifier.reaction_classes[reaction_class]
                )
                results = {reaction_class: result}
            else:
                results = classifier.classify(reaction)

            # Filter to positive only if requested
            if show_positive_only:
                results = {k: v for k, v in results.items() if v.is_member}
                if not results:
                    console.print("[dim]No positive classifications[/dim]")
                    continue

            # Display results
            title = f"Classification Results for {reaction_id}"
            display_classification_results(
                results, format=output_format, verbose=show_explanation, title=title
            )

    except typer.Exit:
        raise
    except FileNotFoundError as e:
        console.print(f"[red]Error: File not found - {e}[/red]")
        raise typer.Exit(1)
    except Exception as e:
        console.print(f"[red]Error: {e}[/red]")
        raise typer.Exit(1)


@app.command()
def list_classes() -> None:
    """List all available reaction classes."""
    classifier = ReactionClassifier()
    display_available_classes(classifier.reaction_classes)


@app.command()
def cache_go(
    output_format: Annotated[
        str, typer.Option("--format", "-f", help="Output format: jsonl or csv")
    ] = "jsonl",
    go_terms: Annotated[
        Optional[str],
        typer.Option(
            "--terms",
            "-t",
            help="Comma-separated list of GO terms to fetch (e.g., GO:0016787,GO:0016740)",
        ),
    ] = None,
    cache_dir: Annotated[
        str, typer.Option("--cache-dir", "-d", help="Cache directory path")
    ] = "cache",
) -> None:
    """Cache GO enzyme classifications to the cache directory.

    Examples:
        autarch cache-go
        autarch cache-go --terms GO:0016787,GO:0016740
        autarch cache-go --format csv
    """
    try:
        # Create cache directory if it doesn't exist
        cache_path = Path(cache_dir)
        cache_path.mkdir(exist_ok=True)

        # Parse GO terms if provided
        go_term_list = None
        if go_terms:
            go_term_list = [term.strip() for term in go_terms.split(",")]
            console.print(
                f"[cyan]Fetching {len(go_term_list)} specific GO terms...[/cyan]"
            )
        else:
            console.print(
                "[cyan]Fetching all molecular function GO terms (this may take a while)...[/cyan]"
            )

        # Fetch GO terms
        go_mappings = fetch_go_enzyme_mappings(go_term_list)

        # Save to cache
        output_file = cache_path / f"go_terms.{output_format}"
        serialize_go_terms(go_mappings, str(output_file), output_format)

        console.print(
            f"[green]✓ Cached {len(go_mappings)} GO terms to {output_file}[/green]"
        )

        # Show summary statistics
        total_rhea = sum(len(term.rhea_ids) for term in go_mappings.values())
        total_ec = sum(len(term.ec_numbers) for term in go_mappings.values())
        console.print(f"[dim]Total RHEA mappings: {total_rhea}[/dim]")
        console.print(f"[dim]Total EC mappings: {total_ec}[/dim]")

    except Exception as e:
        console.print(f"[red]Error caching GO terms: {e}[/red]")
        raise typer.Exit(1)


@app.command()
def cache_rhea(
    output_format: Annotated[
        str, typer.Option("--format", "-f", help="Output format: jsonl or csv")
    ] = "jsonl",
    rhea_ids: Annotated[
        Optional[str],
        typer.Option(
            "--ids",
            "-i",
            help="Comma-separated list of RHEA IDs to fetch (e.g., RHEA:18777,RHEA:18778)",
        ),
    ] = None,
    limit: Annotated[
        Optional[int],
        typer.Option(
            "--limit", "-l", help="Limit number of reactions to fetch (for testing)"
        ),
    ] = None,
    cache_dir: Annotated[
        str, typer.Option("--cache-dir", "-d", help="Cache directory path")
    ] = "cache",
    only_with_go: Annotated[
        bool,
        typer.Option(
            "--only-with-go", help="Only fetch RHEA reactions that have GO annotations"
        ),
    ] = False,
    only_with_ec: Annotated[
        bool,
        typer.Option(
            "--only-with-ec", help="Only fetch RHEA reactions that have EC numbers"
        ),
    ] = False,
    test_mode: Annotated[
        bool,
        typer.Option(
            "--test",
            "-t",
            help="Test mode: saves to test cache, doesn't overwrite main cache",
        ),
    ] = False,
    include_participants: Annotated[
        bool,
        typer.Option(
            "--participants/--no-participants",
            help="Include reaction participants (substrates/products)",
        ),
    ] = True,
) -> None:
    """Cache RHEA reactions to the cache directory.

    By default, fetches ALL reactions from RHEA (~12,000 reactions with EC).
    Use --only-with-ec to limit to reactions with EC numbers.
    Use --test mode for testing without overwriting the main cache.

    Examples:
        autarch cache-rhea                        # Fetch all reactions
        autarch cache-rhea --only-with-ec         # Only reactions with EC numbers
        autarch cache-rhea --ids RHEA:18777,RHEA:18778  # Specific reactions
        autarch cache-rhea --limit 100 --test     # Test with 100 reactions
    """
    try:
        # Create cache directory if it doesn't exist
        cache_path = Path(cache_dir)
        cache_path.mkdir(exist_ok=True)

        # Parse RHEA IDs if provided
        rhea_id_list = None
        if rhea_ids:
            rhea_id_list = [id.strip() for id in rhea_ids.split(",")]
            console.print(
                f"[cyan]Fetching {len(rhea_id_list)} specific RHEA reactions...[/cyan]"
            )
        else:
            if limit:
                console.print(
                    f"[cyan]Fetching {limit} RHEA reactions (test mode)...[/cyan]"
                )
            elif only_with_ec:
                console.print(
                    "[cyan]Fetching all RHEA reactions with EC numbers (~7,500 reactions)...[/cyan]"
                )
                console.print("[dim]This may take several minutes...[/dim]")
            else:
                console.print(
                    "[cyan]Fetching ALL RHEA reactions (~17,000 reactions)...[/cyan]"
                )
                console.print(
                    "[yellow]⚠️  This will take 15-30 minutes. Use --only-with-ec or --limit for faster results.[/yellow]"
                )

        # If --only-with-go, get RHEA IDs from GO cache
        if only_with_go:
            go_cache_file = cache_path / "go_terms.jsonl"
            if not go_cache_file.exists():
                console.print(
                    f"[red]Error: GO cache not found at {go_cache_file}[/red]"
                )
                console.print("[yellow]Run 'autarch cache-go' first[/yellow]")
                raise typer.Exit(1)

            # Collect all RHEA IDs with GO mappings
            rhea_ids_from_go = set()
            with open(go_cache_file) as f:
                for line in f:
                    term_data = json.loads(line)
                    rhea_ids_from_go.update(term_data.get("rhea_ids", []))

            console.print(
                f"[cyan]Found {len(rhea_ids_from_go)} RHEA IDs with GO mappings[/cyan]"
            )

            # Convert to numeric and compute master IDs
            # RHEA uses hierarchical IDs: master (divisible by 4) + 1/2/3 for directional
            master_ids = set()
            for rid in rhea_ids_from_go:
                if rid.startswith("RHEA:"):
                    numeric_id = int(rid[5:])
                    # Round down to nearest multiple of 4 to get master ID
                    master_id = (numeric_id // 4) * 4
                    master_ids.add(master_id)
                    # Also add all directional variants
                    for offset in range(4):
                        master_ids.add(master_id + offset)

            console.print(
                f"[cyan]Expanded to {len(master_ids)} IDs (including directional variants)[/cyan]"
            )

            # Convert to list format expected by ETL
            rhea_id_list = [str(mid) for mid in master_ids]

        # Fetch RHEA reactions using TSV ETL
        console.print("[green]Using TSV-based ETL (fast, reliable)[/green]")
        from autarch.etl.rhea_tsv_etl import RheaTSVETL

        etl = RheaTSVETL(cache_dir=cache_path / "rhea_tsv")
        reactions = etl.load_reactions_batch(rhea_ids=rhea_id_list, limit=limit)

        # Pass GO cache path for GO term mappings
        go_cache_path = cache_path / "go_terms.jsonl"
        rhea_list = etl.create_rhea_terms(
            reactions,
            include_go=True,
            go_cache_path=go_cache_path,
            master_only=True  # Only store master (direction-neutral) IDs
        )

        # Convert list to dict for compatibility
        rhea_mappings = {term.rhea_id.replace("RHEA:", ""): term for term in rhea_list}

        # Save to cache
        if test_mode:
            # Save to test cache to avoid overwriting main cache
            output_file = cache_path / f"rhea_reactions_test.{output_format}"
        else:
            output_file = cache_path / f"rhea_reactions.{output_format}"
        serialize_rhea_reactions(rhea_mappings, str(output_file), output_format)

        # Save the RDKit cache (InChI, canonical SMILES, chiral centers)
        save_rdkit_cache()
        console.print("[dim]Saved RDKit molecular data cache[/dim]")

        console.print(
            f"[green]✓ Cached {len(rhea_mappings)} RHEA reactions to {output_file}[/green]"
        )

        # Show summary statistics
        with_inputs = sum(
            1
            for r in rhea_mappings.values()
            if r.reaction and r.reaction.left_participants
        )
        with_outputs = sum(
            1
            for r in rhea_mappings.values()
            if r.reaction and r.reaction.right_participants
        )
        console.print(f"[dim]Reactions with inputs: {with_inputs}[/dim]")
        console.print(f"[dim]Reactions with outputs: {with_outputs}[/dim]")

    except Exception as e:
        console.print(f"[red]Error caching RHEA reactions: {e}[/red]")
        raise typer.Exit(1)


@app.command()
def cache_chebi(
    cache_dir: Annotated[
        str, typer.Option("--cache-dir", "-d", help="Cache directory path")
    ] = "cache",
) -> None:
    """Cache CHEBI SMILES for all CHEBI IDs in cached RHEA reactions.

    This command reads the cached RHEA reactions and fetches SMILES
    strings for all CHEBI IDs found in the reactions.

    Examples:
        autarch cache-chebi
        autarch cache-chebi --cache-dir test-cache
    """
    try:
        from autarch.etl.chebi_smiles import fetch_chebi_smiles_from_rhea

        cache_path = Path(cache_dir)
        rhea_cache = cache_path / "rhea_reactions.jsonl"

        if not rhea_cache.exists():
            # Check for JSON format too
            rhea_cache = cache_path / "rhea_reactions.json"
            if not rhea_cache.exists():
                console.print(f"[red]Error: RHEA cache not found at {cache_path}[/red]")
                console.print("[yellow]Run 'autarch cache-rhea' first[/yellow]")
                raise typer.Exit(1)

        console.print("[cyan]Fetching SMILES for ChEBI IDs in RHEA reactions...[/cyan]")
        console.print("[dim]This may take several minutes for large datasets...[/dim]")

        # Fetch SMILES for all ChEBI IDs in RHEA
        chebi_to_smiles = fetch_chebi_smiles_from_rhea(cache_dir)

        # Save to cache
        output_file = cache_path / "chebi_smiles.json"
        with open(output_file, "w") as f:
            json.dump(chebi_to_smiles, f, indent=2)

        console.print(
            f"[green]✓ Cached SMILES for {len(chebi_to_smiles)} CHEBI IDs to {output_file}[/green]"
        )

    except Exception as e:
        console.print(f"[red]Error caching CHEBI SMILES: {e}[/red]")
        raise typer.Exit(1)


@app.command()
def cache_smiles_lookup(
    cache_dir: Annotated[
        str, typer.Option("--cache-dir", "-d", help="Cache directory path")
    ] = "cache",
) -> None:
    """Build SMILES→ChEBI reverse lookup from ChEBI database.

    This command builds a cache that maps canonical SMILES strings to
    ChEBI IDs, enabling lookup of ChEBI identifiers from reaction SMILES.

    Requires the ChEBI SQLite database at ~/.data/oaklib/chebi.db
    (install with: runoak cache-xref-sqlite -p chebi)

    Examples:
        autarch cache-smiles-lookup
        autarch cache-smiles-lookup --cache-dir test-cache
    """
    from autarch.etl.chebi_smiles import cache_smiles_to_chebi

    cache_path = Path(cache_dir)
    cache_path.mkdir(exist_ok=True)

    console.print("[cyan]Building SMILES→ChEBI reverse lookup...[/cyan]")
    console.print("[dim]This reads the ChEBI database and canonicalizes all SMILES[/dim]")

    cache_smiles_to_chebi(cache_dir)

    console.print("[green]✓ SMILES lookup cache built successfully[/green]")


@app.command()
def cache_modelseed(
    cache_dir: Annotated[
        str, typer.Option("--cache-dir", "-d", help="Cache directory path")
    ] = "cache",
    limit: Annotated[
        Optional[int],
        typer.Option(
            "--limit",
            "-l",
            help="Limit number of ModelSEED reactions to transform (for testing)",
        ),
    ] = None,
    force_download: Annotated[
        bool,
        typer.Option(
            "--force",
            help="Redownload raw ModelSEED files even if cached locally",
        ),
    ] = False,
) -> None:
    """Cache ModelSEED reactions and a ModelSEED→ChEBI bridge.

    Uses the local SMILES→ChEBI cache built by ``autarch cache-smiles-lookup``.

    Examples:
        autarch cache-modelseed
        autarch cache-modelseed --limit 100
        autarch cache-modelseed --cache-dir test-cache --force
    """
    try:
        from autarch.etl.modelseed_etl import cache_modelseed_dataset

        console.print(
            "[cyan]Caching ModelSEED compounds, reactions, and ChEBI bridge...[/cyan]"
        )
        console.print(
            "[dim]This uses local ChEBI structure caches for identifier mapping[/dim]"
        )

        summary = cache_modelseed_dataset(
            cache_dir=cache_dir,
            limit=limit,
            force_download=force_download,
        )

        console.print(
            f"[green]✓ Cached {summary['total_compounds']} ModelSEED compounds[/green]"
        )
        console.print(f"[dim]Mapped to CHEBI: {summary['mapped_compounds']}[/dim]")
        console.print(f"[dim]Cached reactions: {summary['total_reactions']}[/dim]")
        console.print(
            f"[dim]Reactions with RHEA aliases: {summary['reactions_with_rhea']}[/dim]"
        )
        console.print(
            f"[dim]Fully mapped reactions: {summary['fully_mapped_reactions']}[/dim]"
        )
    except FileNotFoundError as e:
        console.print(f"[red]Error: {e}[/red]")
        console.print("[yellow]Run 'autarch cache-smiles-lookup' first[/yellow]")
        raise typer.Exit(1)
    except Exception as e:
        console.print(f"[red]Error caching ModelSEED bridge: {e}[/red]")
        raise typer.Exit(1)


@app.command()
def summarize_modelseed(
    cache_dir: Annotated[
        str, typer.Option("--cache-dir", "-d", help="Cache directory path")
    ] = "cache",
    output_format: Annotated[
        str, typer.Option("--format", "-f", help="Output format: table or json")
    ] = "table",
) -> None:
    """Summarize cached ModelSEED mapping coverage."""
    try:
        from autarch.etl.modelseed_etl import summarize_modelseed_cache

        summary = summarize_modelseed_cache(cache_dir)

        if output_format == "json":
            typer.echo(json.dumps(summary, indent=2))
            return

        if output_format != "table":
            console.print(f"[red]Error: Unsupported format '{output_format}'[/red]")
            raise typer.Exit(1)

        overview = Table(title="ModelSEED Coverage Summary")
        overview.add_column("Metric", style="cyan")
        overview.add_column("Value", style="green", justify="right")
        overview.add_row("Total compounds", str(summary["total_compounds"]))
        overview.add_row("Mapped compounds to CHEBI", str(summary["mapped_compounds"]))
        overview.add_row("Total reactions", str(summary["total_reactions"]))
        overview.add_row("Reactions with RHEA aliases", str(summary["reactions_with_rhea"]))
        overview.add_row("EC rows", str(summary["ec_coverage"]["rows"]))
        overview.add_row(
            "Distinct reactions with EC",
            str(summary["ec_coverage"]["distinct_reactions"]),
        )
        overview.add_row("Cached transformed reactions", str(summary["cached_reactions"]))
        overview.add_row(
            "Cached reactions with RHEA",
            str(summary["cached_reactions_with_rhea"]),
        )
        overview.add_row(
            "Cached fully mapped reactions",
            str(summary["cached_fully_mapped_reactions"]),
        )
        console.print(overview)

        method_table = Table(title="Compound Mapping Methods")
        method_table.add_column("Method", style="cyan")
        method_table.add_column("Count", style="green", justify="right")
        for method, count in summary["compound_mapping_methods"].items():
            method_table.add_row(method, str(count))
        console.print(method_table)

        reaction_table = Table(title="Reaction Alias Coverage")
        reaction_table.add_column("Source", style="cyan")
        reaction_table.add_column("Rows", style="green", justify="right")
        reaction_table.add_column("Distinct Reactions", style="yellow", justify="right")
        for source, counts in summary["reaction_alias_sources"].items():
            reaction_table.add_row(
                source,
                str(counts["rows"]),
                str(counts["distinct_ids"]),
            )
        console.print(reaction_table)

        compound_table = Table(title="Compound Alias Coverage")
        compound_table.add_column("Source", style="cyan")
        compound_table.add_column("Rows", style="green", justify="right")
        compound_table.add_column("Distinct Compounds", style="yellow", justify="right")
        for source, counts in summary["compound_alias_sources"].items():
            compound_table.add_row(
                source,
                str(counts["rows"]),
                str(counts["distinct_ids"]),
            )
        console.print(compound_table)
    except FileNotFoundError as e:
        console.print(f"[red]Error: {e}[/red]")
        raise typer.Exit(1)
    except typer.Exit:
        raise
    except Exception as e:
        console.print(f"[red]Error summarizing ModelSEED cache: {e}[/red]")
        raise typer.Exit(1)


@app.command()
def benchmark_modelseed(
    cache_dir: Annotated[
        str, typer.Option("--cache-dir", "-d", help="Cache directory path")
    ] = "cache",
    output_dir: Annotated[
        Optional[str],
        typer.Option(
            "--output-dir",
            "-o",
            help="Directory for benchmark artifacts (default: <cache-dir>/modelseed_benchmark)",
        ),
    ] = None,
    limit: Annotated[
        Optional[int],
        typer.Option(
            "--limit",
            "-l",
            help="Limit selected benchmark reactions after filtering",
        ),
    ] = None,
    include_rhea: Annotated[
        bool,
        typer.Option(
            "--include-rhea",
            help="Keep ModelSEED reactions that already have RHEA aliases",
        ),
    ] = False,
    include_transport: Annotated[
        bool,
        typer.Option(
            "--include-transport",
            help="Keep transport reactions in the benchmark slice",
        ),
    ] = False,
    all_statuses: Annotated[
        bool,
        typer.Option(
            "--all-statuses",
            help="Keep non-OK ModelSEED status codes instead of restricting to status=OK",
        ),
    ] = False,
    allow_partial_mapping: Annotated[
        bool,
        typer.Option(
            "--allow-partial-mapping",
            help="Keep reactions that are not fully mapped into ChEBI",
        ),
    ] = False,
    require_ec: Annotated[
        bool,
        typer.Option(
            "--require-ec",
            help="Restrict the benchmark slice to reactions with at least one EC number",
        ),
    ] = False,
    output_format: Annotated[
        str, typer.Option("--format", "-f", help="Output format: table or json")
    ] = "table",
) -> None:
    """Materialize a filtered ModelSEED benchmark slice and cached predictions."""
    try:
        from autarch.modelseed_benchmark import build_modelseed_benchmark

        if output_format != "json":
            console.print("[cyan]Building ModelSEED benchmark slice...[/cyan]")
            console.print(
                "[dim]Default slice: non-RHEA, status=OK, non-transport, fully mapped[/dim]"
            )

        summary = build_modelseed_benchmark(
            cache_dir=cache_dir,
            output_dir=output_dir,
            limit=limit,
            include_rhea=include_rhea,
            include_transport=include_transport,
            all_statuses=all_statuses,
            allow_partial_mapping=allow_partial_mapping,
            require_ec=require_ec,
        )

        if output_format == "json":
            typer.echo(json.dumps(summary, indent=2))
            return

        if output_format != "table":
            console.print(f"[red]Error: Unsupported format '{output_format}'[/red]")
            raise typer.Exit(1)

        overview = Table(title="ModelSEED Benchmark Summary")
        overview.add_column("Metric", style="cyan")
        overview.add_column("Value", style="green", justify="right")
        overview.add_row("Input cached reactions", str(summary["input_reactions"]))
        overview.add_row("Selected benchmark reactions", str(summary["selected_reactions"]))
        overview.add_row("Selected with EC", str(summary["selected_with_ec"]))
        overview.add_row("Selected without EC", str(summary["selected_without_ec"]))
        overview.add_row("Positive predictions", str(summary["positive_reactions"]))
        overview.add_row("Positive rate", str(summary["positive_rate"]))
        overview.add_row("Classifier count", str(summary["classifiers_run"]))
        overview.add_row("Runtime (s)", str(summary["duration_seconds"]))
        console.print(overview)

        skip_table = Table(title="Filtered Out")
        skip_table.add_column("Reason", style="cyan")
        skip_table.add_column("Count", style="green", justify="right")
        for reason, count in summary["skipped_counts"].items():
            skip_table.add_row(reason, str(count))
        console.print(skip_table)

        class_table = Table(title="Top Positive Classes")
        class_table.add_column("Class", style="cyan")
        class_table.add_column("Count", style="green", justify="right")
        for class_name, count in summary["top_positive_classes"][:15]:
            class_table.add_row(class_name, str(count))
        console.print(class_table)

        artifact_table = Table(title="Artifacts")
        artifact_table.add_column("Artifact", style="cyan")
        artifact_table.add_column("Path", style="green")
        for label, path in summary["artifacts"].items():
            artifact_table.add_row(label, str(path))
        console.print(artifact_table)
    except FileNotFoundError as e:
        console.print(f"[red]Error: {e}[/red]")
        console.print("[yellow]Run 'autarch cache-smiles-lookup' first[/yellow]")
        raise typer.Exit(1)
    except typer.Exit:
        raise
    except Exception as e:
        console.print(f"[red]Error building ModelSEED benchmark: {e}[/red]")
        raise typer.Exit(1)


@app.command()
def fetch_rhea(
    rhea_ids: Annotated[
        List[str], typer.Argument(help="RHEA IDs to fetch (e.g., 10000 10001 10002)")
    ],
    output_format: Annotated[
        str, typer.Option("--format", "-f", help="Output format")
    ] = "json",
    output_file: Annotated[
        Optional[str],
        typer.Option("--output", "-o", help="Output file path (default: stdout)"),
    ] = None,
    verbose: Annotated[
        bool, typer.Option("--verbose", "-v", help="Show progress for each reaction")
    ] = False,
    cache_dir: Annotated[
        str, typer.Option("--cache-dir", "-d", help="Cache directory path")
    ] = "cache",
) -> None:
    """Fetch specific RHEA reactions using TSV data.

    This command fetches reaction data from RHEA TSV files.

    Examples:
        autarch fetch-rhea 10000 10001 10002
        autarch fetch-rhea 10000 --format json --output reaction.json
        autarch fetch-rhea 10192 10260 -v
    """
    try:
        from autarch.etl.rhea_tsv_etl import RheaTSVETL
        from pathlib import Path
        import json

        console.print(
            f"[cyan]Fetching {len(rhea_ids)} RHEA reactions from TSV data...[/cyan]"
        )

        etl = RheaTSVETL(cache_dir=Path(cache_dir) / "rhea_tsv")

        # Clean IDs (remove RHEA: prefix if present)
        clean_ids = [rid.replace("RHEA:", "") for rid in rhea_ids]

        reactions = etl.load_reactions_batch(rhea_ids=clean_ids)
        rhea_terms = etl.create_rhea_terms(reactions)

        if verbose:
            for term in rhea_terms:
                console.print(f"[dim]Fetched {term.rhea_id}[/dim]")
                if term.reaction:
                    console.print(
                        f"  [dim]Left: {len(term.reaction.left_participants)} compounds[/dim]"
                    )
                    console.print(
                        f"  [dim]Right: {len(term.reaction.right_participants)} compounds[/dim]"
                    )

        console.print(
            f"[green]✓ Successfully fetched {len(rhea_terms)} reactions[/green]"
        )

        # Format output
        if output_format == "json":
            output_data = json.dumps(
                [rhea_term.model_dump() for rhea_term in rhea_terms],
                indent=2,
                default=str,
            )
        elif output_format == "jsonl":
            output_data = "\n".join(
                json.dumps(rhea_term.model_dump(), default=str)
                for rhea_term in rhea_terms
            )
        else:
            # Simple text format
            lines = []
            for rhea_term in rhea_terms:
                lines.append(f"{rhea_term.rhea_id}: {rhea_term.label}")
                if rhea_term.reaction:
                    left_names = [
                        p.name or p.smiles[:20] + "..." if p.smiles else "?"
                        for p in rhea_term.reaction.left_participants
                    ]
                    right_names = [
                        p.name or p.smiles[:20] + "..." if p.smiles else "?"
                        for p in rhea_term.reaction.right_participants
                    ]
                    if left_names:
                        lines.append(f"  Left: {', '.join(left_names)}")
                    if right_names:
                        lines.append(f"  Right: {', '.join(right_names)}")
            output_data = "\n".join(lines)

        # Output results
        if output_file:
            with open(output_file, "w") as f:
                f.write(output_data)
            console.print(f"[green]✓ Saved to {output_file}[/green]")
        else:
            console.print(output_data)

    except Exception as e:
        console.print(f"[red]Error fetching RHEA reactions: {e}[/red]")
        raise typer.Exit(1)


@app.command()
def eval(
    reaction_class: Annotated[
        str, typer.Argument(help="Reaction class to evaluate (e.g., 'Hydrolase')")
    ],
    go_term: Annotated[
        Optional[str],
        typer.Option(
            "--go",
            "-g",
            help="GO term ID for the reaction class (e.g., GO:0016787 for hydrolase)",
        ),
    ] = None,
    ec_prefix: Annotated[
        Optional[str],
        typer.Option(
            "--ec",
            "-e",
            help="EC number prefix (e.g., 3.1.1 for esterases)",
        ),
    ] = None,
    verbose: Annotated[
        bool,
        typer.Option("--verbose", "-v", help="Show detailed results for each reaction"),
    ] = False,
    cache_dir: Annotated[
        str, typer.Option("--cache-dir", "-d", help="Cache directory path")
    ] = "cache",
    go_only: Annotated[
        bool,
        typer.Option(
            "--go-only/--include-ec",
            help="Only evaluate reactions with GO mappings (default). "
            "Use --include-ec to also use EC numbers as ground truth.",
        ),
    ] = True,
) -> None:
    """Evaluate a reaction class against cached GO and RHEA data.

    This command loads cached GO and RHEA data and evaluates how well
    a reaction class classifier matches the expected RHEA reactions
    based on GO classifications.

    By default, only reactions with GO term mappings are evaluated (--go-only).
    This ensures proper ground truth. Use --include-ec to also use EC numbers.

    Examples:
        autarch eval Hydrolase --go GO:0016787
        autarch eval Hydrolase --verbose
        autarch eval Hydrolase --include-ec  # Include EC-based ground truth
    """
    try:
        console.print(f"[cyan]Evaluating {reaction_class} classification...[/cyan]")

        # Try to get GO_ID from class if not provided
        if not go_term:
            from autarch.classifier import ReactionClassifier

            classifier = ReactionClassifier()
            if reaction_class in classifier.reaction_classes:
                reaction_class_obj = classifier.reaction_classes[reaction_class]
                if hasattr(reaction_class_obj, "GO_ID"):
                    go_term = reaction_class_obj.GO_ID
                    console.print(f"[dim]Using GO_ID from class: {go_term}[/dim]")

        # Run evaluation
        metrics = evaluate_reaction_class(
            reaction_class=reaction_class,
            go_term=go_term,
            ec_prefix=ec_prefix,
            cache_dir=cache_dir,
            verbose=verbose,
            go_only=go_only,
        )

        # Display evaluation metrics
        console.print("\n[green]Evaluation complete![/green]")

        # Create metrics table
        table = Table(title=f"Classification Metrics for {reaction_class}")
        table.add_column("Metric", style="cyan", width=25)
        table.add_column("Value", style="green", justify="right")

        # Confusion matrix values
        table.add_row("True Positives (TP)", str(metrics.true_positives))
        table.add_row("False Positives (FP)", str(metrics.false_positives))
        table.add_row("True Negatives (TN)", str(metrics.true_negatives))
        table.add_row("False Negatives (FN)", str(metrics.false_negatives))
        table.add_row("", "")  # Separator

        # Performance metrics
        table.add_row("Precision", f"{metrics.precision:.3f}")
        table.add_row("Recall (Sensitivity)", f"{metrics.recall:.3f}")
        table.add_row("Specificity", f"{metrics.specificity:.3f}")
        table.add_row("Accuracy", f"{metrics.accuracy:.3f}")
        table.add_row("F1 Score", f"{metrics.f1_score:.3f}")
        table.add_row("MCC", f"{metrics.matthews_correlation_coefficient:.3f}")
        table.add_row("", "")  # Separator

        # Summary
        table.add_row("Total Evaluated", str(metrics.total))

        console.print(table)

        # Always show false positives and false negatives for debugging
        # Load reaction data for labels
        cache_path = Path(cache_dir)
        rhea_cache = cache_path / "rhea_reactions.jsonl"
        rhea_reactions = {}
        if rhea_cache.exists():
            with open(rhea_cache) as f:
                for line in f:
                    reaction_data = json.loads(line)
                    rhea_reactions[reaction_data["rhea_id"]] = reaction_data

        # Show false positives (should be negative but classified as positive)
        if metrics.fp_reactions:
            n_show = min(10, len(metrics.fp_reactions))
            console.print(
                f"\n[red]Top {n_show} False Positives (predicted {reaction_class} but not actually):[/red]"
            )
            for rhea_id in metrics.fp_reactions[:n_show]:
                label = rhea_reactions.get(rhea_id, {}).get("label", "No label")
                explanation = (
                    metrics.fp_explanations.get(rhea_id, "No explanation available")
                    if metrics.fp_explanations
                    else "No explanation"
                )
                console.print(f"  • {rhea_id}: {label}")
                console.print(f"    [dim]→ {explanation}[/dim]")

        # Show false negatives (should be positive but classified as negative)
        if metrics.fn_reactions:
            n_show = min(10, len(metrics.fn_reactions))
            console.print(
                f"\n[yellow]Top {n_show} False Negatives (actually {reaction_class} but not predicted):[/yellow]"
            )
            for rhea_id in metrics.fn_reactions[:n_show]:
                label = rhea_reactions.get(rhea_id, {}).get("label", "No label")
                explanation = (
                    metrics.fn_explanations.get(rhea_id, "No explanation available")
                    if metrics.fn_explanations
                    else "No explanation"
                )
                console.print(f"  • {rhea_id}: {label}")
                console.print(f"    [dim]→ {explanation}[/dim]")

        # Show true positives only in verbose mode
        if verbose and metrics.tp_reactions:
            n_show = min(10, len(metrics.tp_reactions))
            console.print(
                f"\n[green]Top {n_show} True Positives (correctly identified):[/green]"
            )
            for rhea_id in metrics.tp_reactions[:n_show]:
                label = rhea_reactions.get(rhea_id, {}).get("label", "No label")
                console.print(f"  • {rhea_id}: {label}")

    except Exception as e:
        console.print(f"[red]Error during evaluation: {e}[/red]")
        raise typer.Exit(1)


@app.command()
def members(
    reaction_class: Annotated[
        str,
        typer.Argument(help="Reaction class name (e.g., 'Kinase', 'Hydrolase')"),
    ],
    cache_dir: Annotated[
        str, typer.Option("--cache-dir", "-d", help="Cache directory path")
    ] = "cache",
    show_ec: Annotated[
        bool, typer.Option("--show-ec", "-e", help="Show EC numbers")
    ] = True,
    show_go: Annotated[
        bool, typer.Option("--show-go", "-g", help="Show GO terms")
    ] = True,
    limit: Annotated[
        Optional[int], typer.Option("--limit", "-l", help="Limit number of results")
    ] = None,
    verbose: Annotated[
        bool, typer.Option("--verbose", "-v", help="Show reaction details")
    ] = False,
) -> None:
    """Show asserted RHEA reactions for a reaction class.
    
    This command shows which RHEA reactions are expected to belong to a 
    reaction class based on GO term and EC number associations in the 
    cached data. This helps debug classifier performance by showing what
    the ground truth data expects.
    
    Examples:
        autarch members Kinase
        autarch members Hydrolase --limit 10
        autarch members Oxidoreductase --verbose
        autarch members TransferaseTransferringPhosphorusContainingGroups --show-ec
    """
    try:
        # Initialize classifier to get the reaction class
        classifier = ReactionClassifier()
        
        if reaction_class not in classifier.reaction_classes:
            console.print(f"[red]Error: Unknown reaction class '{reaction_class}'[/red]")
            console.print(f"[dim]Available classes: {', '.join(classifier.reaction_classes.keys())}[/dim]")
            raise typer.Exit(1)
        
        reaction_class_obj = classifier.reaction_classes[reaction_class]
        
        # Get GO term and EC prefix from class
        go_term = getattr(reaction_class_obj, "GO_ID", None)
        ec_prefix = getattr(reaction_class_obj, "EC_NUMBER_PREFIX", None)
        
        console.print(f"[cyan]Finding RHEA reactions for {reaction_class}[/cyan]")
        if go_term:
            console.print(f"[dim]GO term: {go_term}[/dim]")
        if ec_prefix:
            console.print(f"[dim]EC prefix: {ec_prefix}[/dim]")
        
        # Load cached data
        cache_path = Path(cache_dir)
        
        # Load GO terms
        go_terms = {}
        go_cache = cache_path / "go_terms.jsonl"
        if go_cache.exists():
            with open(go_cache) as f:
                for line in f:
                    go_data = json.loads(line)
                    go_terms[go_data["go_id"]] = go_data
        
        # Load RHEA reactions
        rhea_reactions = {}
        rhea_cache = cache_path / "rhea_reactions.jsonl"
        if rhea_cache.exists():
            with open(rhea_cache) as f:
                for line in f:
                    reaction_data = json.loads(line)
                    rhea_reactions[reaction_data["rhea_id"]] = reaction_data
        
        # Find positive RHEA IDs (same logic as evaluation)
        positive_rhea_ids = set()
        
        # Check EC numbers if provided
        if ec_prefix:
            for rhea_id, rhea_data in rhea_reactions.items():
                ec_numbers = rhea_data.get("ec_numbers", [])
                for ec in ec_numbers:
                    if ec.startswith(ec_prefix):
                        positive_rhea_ids.add(rhea_id)
                        break
        
        # Check GO terms if provided
        if go_term and go_term in go_terms:
            # Include RHEA IDs from the term itself
            term_data = go_terms[go_term]
            positive_rhea_ids.update(term_data.get("rhea_ids", []))
            
            # Find all RHEA reactions associated with descendants
            for term in go_terms.values():
                if go_term in term.get("ancestors", []) and term["go_id"] != go_term:
                    positive_rhea_ids.update(term.get("rhea_ids", []))
            
            # Also check RHEA reactions directly
            for rhea_id, rhea_data in rhea_reactions.items():
                for rhea_go_term in rhea_data.get("go_terms", []):
                    if rhea_go_term in go_terms:
                        term_ancestors = go_terms[rhea_go_term].get("ancestors", [])
                        if go_term in term_ancestors or rhea_go_term == go_term:
                            positive_rhea_ids.add(rhea_id)
        
        # Sort RHEA IDs
        sorted_rhea_ids = sorted(positive_rhea_ids)
        
        # Apply limit if specified
        if limit:
            sorted_rhea_ids = sorted_rhea_ids[:limit]
        
        # Display results
        console.print(f"\n[green]Found {len(positive_rhea_ids)} RHEA reactions[/green]")
        if limit and len(positive_rhea_ids) > limit:
            console.print(f"[dim]Showing first {limit}[/dim]")
        
        # Create table
        table = Table(title=f"RHEA Reactions for {reaction_class}")
        table.add_column("RHEA ID", style="cyan")
        table.add_column("Reaction", style="white")
        if show_ec:
            table.add_column("EC", style="yellow")
        if show_go:
            table.add_column("GO Terms", style="magenta")
        
        for rhea_id in sorted_rhea_ids:
            rhea_data = rhea_reactions.get(rhea_id, {})
            
            # Get reaction label
            label = rhea_data.get("label", "No label")
            if len(label) > 60 and not verbose:
                label = label[:57] + "..."
            
            # Get EC numbers
            ec_numbers = rhea_data.get("ec_numbers", [])
            ec_str = ", ".join(ec_numbers) if ec_numbers else "-"
            
            # Get GO terms
            go_term_list = rhea_data.get("go_terms", [])
            go_str = ", ".join(go_term_list[:3]) if go_term_list else "-"
            if len(go_term_list) > 3:
                go_str += f" (+{len(go_term_list)-3})"
            
            # Build row
            row = [rhea_id, label]
            if show_ec:
                row.append(ec_str)
            if show_go:
                row.append(go_str)
            
            table.add_row(*row)
        
        console.print(table)
        
        # Show detailed info in verbose mode
        if verbose:
            console.print("\n[cyan]Detailed reaction information:[/cyan]")
            for rhea_id in sorted_rhea_ids[:5]:  # Show first 5 in detail
                rhea_data = rhea_reactions.get(rhea_id, {})
                console.print(f"\n[bold]{rhea_id}[/bold]: {rhea_data.get('label', 'No label')}")
                
                # Show participants
                reaction = rhea_data.get("reaction", {})
                if reaction:
                    left = reaction.get("left_participants", [])
                    right = reaction.get("right_participants", [])
                    
                    left_str = " + ".join([p.get("name", p.get("chebi_id", "?")) for p in left])
                    right_str = " + ".join([p.get("name", p.get("chebi_id", "?")) for p in right])
                    console.print(f"  [dim]{left_str} → {right_str}[/dim]")
                
                # Show EC and GO
                if rhea_data.get("ec_numbers"):
                    console.print(f"  EC: {', '.join(rhea_data['ec_numbers'])}")
                if rhea_data.get("go_terms"):
                    console.print(f"  GO: {', '.join(rhea_data['go_terms'][:5])}")
        
    except Exception as e:
        console.print(f"[red]Error: {e}[/red]")
        raise typer.Exit(1)


@app.command()
def eval_all(
    verbose: Annotated[
        bool,
        typer.Option("--verbose", "-v", help="Show detailed results for each class"),
    ] = False,
    cache_dir: Annotated[
        str, typer.Option("--cache-dir", "-d", help="Cache directory path")
    ] = "cache",
    output_dir: Annotated[
        Optional[Path],
        typer.Option(
            "--output-dir",
            "-o",
            help="Directory to save detailed evaluation reports (CSV files, FP/FN reports, etc.)",
        ),
    ] = None,
    visualize: Annotated[
        bool,
        typer.Option(
            "--visualize",
            "-viz",
            help="Generate visualizations after evaluation",
        ),
    ] = False,
    show_plots: Annotated[
        bool,
        typer.Option(
            "--show-plots",
            "-s",
            help="Show interactive plots in browser (requires --visualize)",
        ),
    ] = False,
    go_only: Annotated[
        bool,
        typer.Option(
            "--go-only/--include-ec",
            help="Only evaluate reactions with GO mappings (default). "
            "Use --include-ec to also use EC numbers as ground truth.",
        ),
    ] = True,
) -> None:
    """Evaluate all reaction classes and provide macro/micro statistics.

    This command evaluates all available reaction classifiers against cached
    GO and RHEA data, providing both per-class metrics and overall macro/micro
    statistics. Optionally generates comprehensive visualizations.

    By default, only reactions with GO term mappings are evaluated (--go-only).
    This ensures proper ground truth. Use --include-ec to also use EC numbers.

    Examples:
        autarch eval-all
        autarch eval-all --verbose
        autarch eval-all --output-dir results/
        autarch eval-all -o results/ --verbose
        autarch eval-all --include-ec  # Include EC-based ground truth
    """
    try:
        import pandas as pd

        console.print("[cyan]Evaluating all reaction classes...[/cyan]")

        # Initialize classifier to get all classes
        classifier = ReactionClassifier()

        # Store results for each class
        all_metrics = []
        class_results = {}

        for class_name, reaction_class_obj in classifier.reaction_classes.items():
            console.print(f"\n[dim]Evaluating {class_name}...[/dim]")

            # Get GO term from class if available
            go_term = getattr(reaction_class_obj, "GO_ID", None)

            if not go_term:
                console.print(
                    f"[yellow]Warning: No GO_ID found for {class_name}, skipping[/yellow]"
                )
                continue

            try:
                # Run evaluation for this class
                metrics = evaluate_reaction_class(
                    reaction_class=class_name,
                    go_term=go_term,
                    cache_dir=cache_dir,
                    verbose=False,  # Don't show verbose output for individual classes
                    go_only=go_only,
                )

                class_results[class_name] = metrics

                # Store metrics for aggregation
                all_metrics.append(
                    {
                        "class": class_name,
                        "go_term": go_term,
                        "tp": metrics.true_positives,
                        "fp": metrics.false_positives,
                        "tn": metrics.true_negatives,
                        "fn": metrics.false_negatives,
                        "precision": metrics.precision,
                        "recall": metrics.recall,
                        "specificity": metrics.specificity,
                        "accuracy": metrics.accuracy,
                        "f1_score": metrics.f1_score,
                        "mcc": metrics.matthews_correlation_coefficient,
                        "total": metrics.total,
                    }
                )

            except Exception as e:
                console.print(f"[red]Error evaluating {class_name}: {e}[/red]")
                continue

        if not all_metrics:
            console.print("[red]No classes could be evaluated[/red]")
            raise typer.Exit(1)

        # Convert to DataFrame for easy aggregation
        df = pd.DataFrame(all_metrics)

        # Calculate micro-averaged metrics (pool all predictions)
        micro_tp = df["tp"].sum()
        micro_fp = df["fp"].sum()
        micro_tn = df["tn"].sum()
        micro_fn = df["fn"].sum()

        micro_precision = (
            micro_tp / (micro_tp + micro_fp) if (micro_tp + micro_fp) > 0 else 0
        )
        micro_recall = (
            micro_tp / (micro_tp + micro_fn) if (micro_tp + micro_fn) > 0 else 0
        )
        micro_f1 = (
            2 * (micro_precision * micro_recall) / (micro_precision + micro_recall)
            if (micro_precision + micro_recall) > 0
            else 0
        )
        micro_accuracy = (
            (micro_tp + micro_tn) / (micro_tp + micro_tn + micro_fp + micro_fn)
            if (micro_tp + micro_tn + micro_fp + micro_fn) > 0
            else 0
        )

        # Calculate macro-averaged metrics (average of per-class metrics)
        macro_precision = df["precision"].mean()
        macro_recall = df["recall"].mean()
        macro_f1 = df["f1_score"].mean()
        macro_accuracy = df["accuracy"].mean()
        macro_mcc = df["mcc"].mean()
        class_total_min = int(df["total"].min()) if len(df) > 0 else 0
        class_total_max = int(df["total"].max()) if len(df) > 0 else 0
        benchmark_mode_label = "GO-only" if go_only else "GO+EC merged"

        # Save detailed reports if output directory is specified
        if output_dir:
            from datetime import datetime

            # Create output directory
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            console.print(f"\n[cyan]Saving detailed reports to {output_dir}/[/cyan]")

            # Load RHEA reactions for labels
            cache_path = Path(cache_dir)
            rhea_cache = cache_path / "rhea_reactions.jsonl"
            rhea_reactions = {}
            if rhea_cache.exists():
                with open(rhea_cache) as f:
                    for line in f:
                        reaction_data = json.loads(line)
                        rhea_reactions[reaction_data["rhea_id"]] = reaction_data
            go_linked_count = sum(
                1 for reaction_data in rhea_reactions.values() if reaction_data.get("go_terms")
            )
            annotated_count = sum(
                1
                for reaction_data in rhea_reactions.values()
                if reaction_data.get("go_terms") or reaction_data.get("ec_numbers")
            )

            # 1. Save main results CSV with all metrics
            results_csv = output_dir / "evaluation_results.csv"
            df.to_csv(results_csv, index=False)
            console.print(f"  • Saved evaluation metrics to {results_csv}")

            # 2. Save detailed predictions CSV (one row per reaction outcome)
            detailed_results = []
            for class_name, metrics in class_results.items():
                go_term = next(
                    (m["go_term"] for m in all_metrics if m["class"] == class_name),
                    None,
                )

                # Add true positives
                if metrics.tp_reactions:
                    for rhea_id in metrics.tp_reactions:
                        detailed_results.append(
                            {
                                "class": class_name,
                                "go_term": go_term,
                                "rhea_id": rhea_id,
                                "label": rhea_reactions.get(rhea_id, {}).get(
                                    "label", ""
                                ),
                                "prediction": "positive",
                                "actual": "positive",
                                "outcome": "TP",
                                "explanation": metrics.tp_explanations.get(rhea_id, "")
                                if metrics.tp_explanations
                                else "",
                            }
                        )

                # Add false positives
                if metrics.fp_reactions:
                    for rhea_id in metrics.fp_reactions:
                        detailed_results.append(
                            {
                                "class": class_name,
                                "go_term": go_term,
                                "rhea_id": rhea_id,
                                "label": rhea_reactions.get(rhea_id, {}).get(
                                    "label", ""
                                ),
                                "prediction": "positive",
                                "actual": "negative",
                                "outcome": "FP",
                                "explanation": metrics.fp_explanations.get(rhea_id, "")
                                if metrics.fp_explanations
                                else "",
                            }
                        )

                # Add true negatives
                if metrics.tn_reactions:
                    for rhea_id in metrics.tn_reactions:
                        detailed_results.append(
                            {
                                "class": class_name,
                                "go_term": go_term,
                                "rhea_id": rhea_id,
                                "label": rhea_reactions.get(rhea_id, {}).get(
                                    "label", ""
                                ),
                                "prediction": "negative",
                                "actual": "negative",
                                "outcome": "TN",
                                "explanation": metrics.tn_explanations.get(rhea_id, "")
                                if metrics.tn_explanations
                                else "",
                            }
                        )

                # Add false negatives
                if metrics.fn_reactions:
                    for rhea_id in metrics.fn_reactions:
                        detailed_results.append(
                            {
                                "class": class_name,
                                "go_term": go_term,
                                "rhea_id": rhea_id,
                                "label": rhea_reactions.get(rhea_id, {}).get(
                                    "label", ""
                                ),
                                "prediction": "negative",
                                "actual": "positive",
                                "outcome": "FN",
                                "explanation": metrics.fn_explanations.get(rhea_id, "")
                                if metrics.fn_explanations
                                else "",
                            }
                        )

            if detailed_results:
                detailed_csv = output_dir / "detailed_predictions.csv"
                detailed_df = pd.DataFrame(detailed_results)
                detailed_df.to_csv(detailed_csv, index=False)
                console.print(f"  • Saved detailed predictions to {detailed_csv}")

            # 3. Save FP/FN reports for each class
            for class_name, metrics in class_results.items():
                class_dir = output_dir / class_name.lower()
                class_dir.mkdir(exist_ok=True)

                # Save false positives
                if metrics.fp_reactions:
                    fp_file = class_dir / "false_positives.txt"
                    with open(fp_file, "w") as f:
                        f.write(f"False Positives for {class_name}\n")
                        f.write(f"{'=' * 50}\n")
                        f.write(f"Total: {len(metrics.fp_reactions)}\n\n")
                        for rhea_id in metrics.fp_reactions:
                            label = rhea_reactions.get(rhea_id, {}).get(
                                "label", "No label"
                            )
                            explanation = (
                                metrics.fp_explanations.get(rhea_id, "No explanation")
                                if metrics.fp_explanations
                                else "No explanation"
                            )
                            f.write(f"{rhea_id}: {label}\n")
                            f.write(f"  → {explanation}\n\n")
                    console.print(
                        f"  • Saved {class_name} false positives to {fp_file}"
                    )

                # Save false negatives
                if metrics.fn_reactions:
                    fn_file = class_dir / "false_negatives.txt"
                    with open(fn_file, "w") as f:
                        f.write(f"False Negatives for {class_name}\n")
                        f.write(f"{'=' * 50}\n")
                        f.write(f"Total: {len(metrics.fn_reactions)}\n\n")
                        for rhea_id in metrics.fn_reactions:
                            label = rhea_reactions.get(rhea_id, {}).get(
                                "label", "No label"
                            )
                            explanation = (
                                metrics.fn_explanations.get(rhea_id, "No explanation")
                                if metrics.fn_explanations
                                else "No explanation"
                            )
                            f.write(f"{rhea_id}: {label}\n")
                            f.write(f"  → {explanation}\n\n")
                    console.print(
                        f"  • Saved {class_name} false negatives to {fn_file}"
                    )

            # 4. Save summary statistics
            summary_file = output_dir / "summary_statistics.txt"
            with open(summary_file, "w") as f:
                f.write("Evaluation Summary Report\n")
                f.write(f"Generated: {datetime.now().isoformat()}\n")
                f.write(f"{'=' * 60}\n\n")

                f.write(f"Benchmark Mode: {benchmark_mode_label}\n")
                f.write(f"Classes Evaluated: {len(all_metrics)}\n")
                benchmark_candidate_count = go_linked_count if go_only else annotated_count
                f.write(
                    f"Benchmark Candidate Reactions: {benchmark_candidate_count}\n"
                )
                f.write(
                    f"Class-specific Evaluated Reactions: {class_total_min}-{class_total_max}\n\n"
                )

                f.write("MACRO-AVERAGED METRICS\n")
                f.write("-" * 30 + "\n")
                f.write(f"  Precision: {macro_precision:.3f}\n")
                f.write(f"  Recall: {macro_recall:.3f}\n")
                f.write(f"  F1 Score: {macro_f1:.3f}\n")
                f.write(f"  Accuracy: {macro_accuracy:.3f}\n")
                f.write(f"  MCC: {macro_mcc:.3f}\n\n")

                f.write("MICRO-AVERAGED METRICS\n")
                f.write("-" * 30 + "\n")
                f.write(f"  Precision: {micro_precision:.3f}\n")
                f.write(f"  Recall: {micro_recall:.3f}\n")
                f.write(f"  F1 Score: {micro_f1:.3f}\n")
                f.write(f"  Accuracy: {micro_accuracy:.3f}\n\n")

                f.write("CONFUSION MATRIX TOTALS\n")
                f.write("-" * 30 + "\n")
                f.write(f"  True Positives: {micro_tp}\n")
                f.write(f"  False Positives: {micro_fp}\n")
                f.write(f"  True Negatives: {micro_tn}\n")
                f.write(f"  False Negatives: {micro_fn}\n\n")

                f.write("PER-CLASS SUMMARY\n")
                f.write("-" * 30 + "\n")
                for row in all_metrics:
                    f.write(f"\n{row['class']} ({row['go_term']}):\n")
                    f.write(f"  Precision: {row['precision']:.3f}\n")
                    f.write(f"  Recall: {row['recall']:.3f}\n")
                    f.write(f"  F1 Score: {row['f1_score']:.3f}\n")
                    f.write(f"  MCC: {row['mcc']:.3f}\n")
                    f.write(
                        f"  TP={row['tp']}, FP={row['fp']}, TN={row['tn']}, FN={row['fn']}\n"
                    )

            console.print(f"  • Saved summary statistics to {summary_file}")

            from autarch.complexity import (
                analyze_all_classifiers,
                save_complexity_report,
            )
            from autarch.benchmark_divergence import save_benchmark_divergence_report
            from autarch.curation_objective import save_curation_objective_report

            complexity_report = analyze_all_classifiers()
            save_complexity_report(complexity_report, output_dir)
            console.print(f"  • Saved complexity report to {output_dir}")

            save_curation_objective_report(df, complexity_report, output_dir)
            console.print(f"  • Saved curation objective report to {output_dir}")

            divergence_csv, divergence_summary = save_benchmark_divergence_report(
                output_dir,
                cache_dir,
            )
            console.print(f"  • Saved benchmark divergence report to {divergence_csv}")
            console.print(f"  • Saved benchmark divergence summary to {divergence_summary}")

            console.print(f"\n[green]All reports saved to {output_dir}/[/green]")

        # Display results
        console.print(
            f"\n[green]Evaluation complete for {len(all_metrics)} classes![/green]\n"
        )

        # Per-class results table
        if verbose:
            class_table = Table(title="Per-Class Performance")
            class_table.add_column("Class", style="cyan")
            class_table.add_column("GO Term", style="dim")
            class_table.add_column("Precision", justify="right")
            class_table.add_column("Recall", justify="right")
            class_table.add_column("F1", justify="right")
            class_table.add_column("MCC", justify="right")
            class_table.add_column("Total", justify="right")

            for row in all_metrics:
                class_table.add_row(
                    row["class"],
                    row["go_term"],
                    f"{row['precision']:.3f}",
                    f"{row['recall']:.3f}",
                    f"{row['f1_score']:.3f}",
                    f"{row['mcc']:.3f}",
                    str(row["total"]),
                )

            console.print(class_table)
            console.print()

        # Summary statistics table
        summary_table = Table(title="Overall Statistics")
        summary_table.add_column("Metric Type", style="cyan", width=20)
        summary_table.add_column("Metric", style="yellow", width=25)
        summary_table.add_column("Value", style="green", justify="right")

        # Macro-averaged metrics
        summary_table.add_row("Macro-averaged", "Precision", f"{macro_precision:.3f}")
        summary_table.add_row("Macro-averaged", "Recall", f"{macro_recall:.3f}")
        summary_table.add_row("Macro-averaged", "F1 Score", f"{macro_f1:.3f}")
        summary_table.add_row("Macro-averaged", "Accuracy", f"{macro_accuracy:.3f}")
        summary_table.add_row("Macro-averaged", "MCC", f"{macro_mcc:.3f}")
        summary_table.add_row("", "", "")  # Separator

        # Micro-averaged metrics
        summary_table.add_row("Micro-averaged", "Precision", f"{micro_precision:.3f}")
        summary_table.add_row("Micro-averaged", "Recall", f"{micro_recall:.3f}")
        summary_table.add_row("Micro-averaged", "F1 Score", f"{micro_f1:.3f}")
        summary_table.add_row("Micro-averaged", "Accuracy", f"{micro_accuracy:.3f}")
        summary_table.add_row("", "", "")  # Separator

        # Totals
        summary_table.add_row("Totals", "Benchmark Mode", benchmark_mode_label)
        summary_table.add_row("Totals", "Classes Evaluated", str(len(all_metrics)))
        summary_table.add_row(
            "Totals",
            "Evaluated Reactions / Class",
            f"{class_total_min}-{class_total_max}",
        )
        summary_table.add_row("Totals", "Total TP", str(micro_tp))
        summary_table.add_row("Totals", "Total FP", str(micro_fp))
        summary_table.add_row("Totals", "Total TN", str(micro_tn))
        summary_table.add_row("Totals", "Total FN", str(micro_fn))

        console.print(summary_table)
        
        # Generate visualizations if requested
        if visualize:
            if not output_dir:
                console.print(
                    "[red]Error: --output-dir is required when using --visualize[/red]"
                )
                raise typer.Exit(1)
            
            console.print("\n[cyan]Generating visualizations...[/cyan]")
            
            # Import visualization module
            from autarch import visualization
            import matplotlib
            
            # Set backend for non-interactive mode
            if not show_plots:
                matplotlib.use('Agg')
            
            # Set up visualization directory
            viz_output_dir = output_dir / "visualizations"
            viz_output_dir.mkdir(exist_ok=True, parents=True)
            
            # Setup style
            visualization.setup_style()
            
            # Load the data we just saved
            eval_df = pd.read_csv(output_dir / "evaluation_results.csv")
            if (output_dir / "detailed_predictions.csv").exists():
                pd.read_csv(output_dir / "detailed_predictions.csv")
            
            console.print("[dim]Creating visualizations:[/dim]")
            
            # Generate all visualizations
            console.print("  • Performance bar charts...")
            visualization.create_performance_barchart(eval_df, viz_output_dir)
            
            console.print("  • Interactive heatmap...")
            visualization.create_interactive_heatmap(eval_df, viz_output_dir)
            
            console.print("  • Confusion matrix grid...")
            visualization.create_confusion_matrix_grid(eval_df, viz_output_dir)
            
            console.print("  • Radar charts...")
            visualization.create_radar_chart(eval_df, viz_output_dir)
            
            console.print("  • EC class analysis...")
            visualization.create_ec_class_analysis(eval_df, viz_output_dir)
            
            console.print("  • Performance distributions...")
            visualization.create_performance_distribution(eval_df, viz_output_dir)
            
            console.print("  • 3D performance scatter...")
            visualization.create_interactive_3d_scatter(eval_df, viz_output_dir)
            
            console.print("  • Sunburst hierarchy chart...")
            visualization.create_sunburst_chart(eval_df, viz_output_dir)

            # Complexity analysis
            console.print("  • Complexity analysis...")
            from autarch.complexity import analyze_all_classifiers, save_complexity_report
            complexity_report = analyze_all_classifiers()
            viz_output_dir_str = str(viz_output_dir)
            visualization.create_complexity_chart(complexity_report, viz_output_dir_str)
            visualization.create_complexity_vs_performance(
                eval_df, complexity_report, viz_output_dir_str
            )
            visualization.create_complexity_vs_performance_static(
                eval_df,
                complexity_report,
                viz_output_dir / "complexity_vs_f1.png",
            )
            save_complexity_report(complexity_report, viz_output_dir)

            from autarch.rhea_text_embeddings import save_rhea_browser_html

            console.print("  • RHEA reaction browser...")
            try:
                save_rhea_browser_html(
                    viz_output_dir / "rhea_browser.html",
                    Path(cache_dir),
                    results_dir=output_dir,
                )
            except FileNotFoundError:
                console.print(
                    "    [yellow]Skipped: RHEA cache not available for reaction browser[/yellow]"
                )

            console.print("  • Comprehensive report...")
            visualization.generate_summary_report(eval_df, viz_output_dir, complexity_report)

            console.print(f"\n[green] Visualizations saved to {viz_output_dir}[/green]")
            
            # Open interactive plots if requested
            if show_plots:
                import webbrowser
                console.print("[cyan]Opening interactive visualizations in browser...[/cyan]")
                webbrowser.open(str(viz_output_dir / "interactive_heatmap.html"))

    except Exception as e:
        console.print(f"[red]Error during evaluation: {e}[/red]")
        raise typer.Exit(1)


@app.command()
def visualize(
    results_dir: Annotated[
        Path,
        typer.Option(
            "--results-dir",
            "-r",
            help="Directory containing evaluation results",
        ),
    ] = Path("results"),
    output_dir: Annotated[
        Optional[Path],
        typer.Option(
            "--output-dir",
            "-o",
            help="Directory to save visualizations (default: results/visualizations)",
        ),
    ] = None,
    show_plots: Annotated[
        bool,
        typer.Option(
            "--show/--no-show",
            "-s/-ns",
            help="Show interactive plots in browser",
        ),
    ] = False,
) -> None:
    """Generate comprehensive visualizations from evaluation results.
    
    Creates beautiful charts including:
    - Performance bar charts and heatmaps
    - Confusion matrices grid
    - 3D performance scatter plots
    - Radar charts for top performers
    - EC class analysis
    - Interactive HTML visualizations
    
    Examples:
        autarch visualize
        autarch visualize --results-dir custom-results/
        autarch visualize --show  # Opens interactive plots
        autarch visualize -o output/charts/
    """
    try:
        import matplotlib
        
        # Set backend for non-interactive mode
        if not show_plots:
            matplotlib.use('Agg')
        
        console.print("[cyan]Generating visualizations from evaluation results...[/cyan]")
        
        # Check if results directory exists
        if not results_dir.exists():
            console.print(f"[red]Error: Results directory '{results_dir}' does not exist[/red]")
            console.print("[yellow]Please run 'autarch eval-all --output-dir results' first[/yellow]")
            raise typer.Exit(1)
        
        # Check for required files
        eval_results_file = results_dir / "evaluation_results.csv"
        if not eval_results_file.exists():
            console.print(f"[red]Error: Required file '{eval_results_file}' not found[/red]")
            console.print("[yellow]Please run 'autarch eval-all --output-dir results' first[/yellow]")
            raise typer.Exit(1)
        
        # Set up output directory
        if output_dir is None:
            output_dir = results_dir / "visualizations"
        
        # Create temporary module to override paths
        import sys
        import types
        from typing import Any
        
        # Create a temporary module with our settings
        temp_module: Any = types.ModuleType('viz_config')
        temp_module.RESULTS_DIR = str(results_dir)
        temp_module.OUTPUT_DIR = str(output_dir)
        temp_module.SHOW_PLOTS = show_plots
        sys.modules['viz_config'] = temp_module
        
        # Modify the visualization module to use our settings
        from autarch import visualization
        
        # Override the main function to use our paths
        
        def custom_main():
            visualization.setup_style()
            
            # Create output directory
            from pathlib import Path
            output_path = Path(str(output_dir))
            output_path.mkdir(exist_ok=True, parents=True)
            
            # Load data with custom path
            console.print(f"Loading evaluation results from {results_dir}...")
            eval_df, detailed_df = visualization.load_data(str(results_dir))
            
            console.print(f"Loaded {len(eval_df)} classifier results")
            console.print("\nGenerating visualizations...\n")
            
            # Generate all visualizations
            console.print("1. Creating performance bar charts...")
            visualization.create_performance_barchart(eval_df, output_path)
            
            console.print("2. Creating interactive heatmap...")
            visualization.create_interactive_heatmap(eval_df, output_path)
            
            console.print("3. Creating confusion matrix grid...")
            visualization.create_confusion_matrix_grid(eval_df, output_path)
            
            console.print("4. Creating radar charts for top performers...")
            visualization.create_radar_chart(eval_df, output_path)
            
            console.print("5. Analyzing by EC class...")
            visualization.create_ec_class_analysis(eval_df, output_path)
            
            console.print("6. Creating performance distribution plots...")
            visualization.create_performance_distribution(eval_df, output_path)
            
            console.print("7. Creating 3D performance scatter plot...")
            visualization.create_interactive_3d_scatter(eval_df, output_path)
            
            console.print("8. Creating sunburst hierarchy chart...")
            visualization.create_sunburst_chart(eval_df, output_path)
            
            console.print("9. Creating complexity analysis...")
            from autarch.complexity import analyze_all_classifiers
            complexity_report = analyze_all_classifiers()
            visualization.create_complexity_chart(complexity_report, str(output_path))
            visualization.create_complexity_vs_performance(
                eval_df,
                complexity_report,
                output_path,
            )
            visualization.create_complexity_vs_performance_static(
                eval_df,
                complexity_report,
                output_path / "complexity_vs_f1.png",
            )

            from autarch.rhea_text_embeddings import save_rhea_browser_html

            console.print("10. Creating RHEA reaction browser...")
            try:
                save_rhea_browser_html(
                    output_path / "rhea_browser.html",
                    Path("cache"),
                    results_dir=results_dir,
                )
            except FileNotFoundError:
                console.print("   Skipped: RHEA cache not available for reaction browser")

            console.print("11. Generating comprehensive report...")
            visualization.generate_summary_report(
                eval_df,
                output_path,
                complexity_report,
            )
            
            console.print(f"\n✨ All visualizations saved to {output_path}")
            console.print("\nInteractive visualizations:")
            console.print(f"  - {output_path}/interactive_heatmap.html")
            console.print(f"  - {output_path}/radar_charts.html")
            console.print(f"  - {output_path}/3d_performance_scatter.html")
            console.print(f"  - {output_path}/sunburst_performance.html")
            console.print(f"  - {output_path}/complexity_vs_performance.html")
            console.print(f"  - {output_path}/rhea_browser.html")
            
            console.print("\nStatic visualizations:")
            console.print(f"  - {output_path}/performance_metrics.png")
            console.print(f"  - {output_path}/confusion_matrices.png")
            console.print(f"  - {output_path}/ec_class_analysis.png")
            console.print(f"  - {output_path}/performance_distributions.png")
            console.print(f"  - {output_path}/complexity_vs_f1.png")
            
            console.print(f"\nReport: {output_path}/comprehensive_report.txt")
            
            # Open interactive plots if requested
            if show_plots:
                import webbrowser
                console.print("\n[cyan]Opening interactive visualizations in browser...[/cyan]")
                webbrowser.open(str(output_path / "interactive_heatmap.html"))
        
        # Run the custom main
        custom_main()
        
        console.print("\n[green]✓ Visualization generation complete![/green]")
        
    except Exception as e:
        console.print(f"[red]Error generating visualizations: {e}[/red]")
        raise typer.Exit(1)


@app.command()
def validate_go_names():
    """Validate GO term naming consistency for all reaction classifiers."""
    try:
        validator = GoTermValidator()
        success = validator.validate_against_baseline()
        
        if not success:
            raise typer.Exit(1)
            
    except Exception as e:
        console.print(f"[red]Error during validation: {e}[/red]")
        raise typer.Exit(1)


@app.command("validate-go-cache")
def validate_go_cache_cmd(
    cache_dir: Annotated[
        str, typer.Option("--cache-dir", "-d", help="Cache directory path")
    ] = "cache",
):
    """Validate cached GO ancestor closure against the local GO sqlite."""
    try:
        report = validate_go_cache(Path(cache_dir) / "go_terms.jsonl")
        success = print_go_cache_validation_report(report)
        if not success:
            raise typer.Exit(1)
    except Exception as e:
        console.print(f"[red]Error during GO cache validation: {e}[/red]")
        raise typer.Exit(1)


@app.command()
def validate_go_names_strict():
    """Run strict GO/class/file/description and semantic validation without baseline filtering."""
    try:
        validator = GoTermValidator()
        success = validator.validate()

        if not success:
            raise typer.Exit(1)

    except Exception as e:
        console.print(f"[red]Error during strict validation: {e}[/red]")
        raise typer.Exit(1)


@app.command()
def summarize_classes():
    """Generate comprehensive summary table of all reaction classes with GO, EC, and RHEA mappings."""
    try:
        summarizer = ClassSummarizer()
        summarizer.generate_summary()

    except Exception as e:
        console.print(f"[red]Error generating class summary: {e}[/red]")
        raise typer.Exit(1)


@app.command()
def export(
    results_dir: Annotated[
        Path,
        typer.Option(
            "--results-dir",
            "-r",
            help="Directory containing evaluation results (evaluation_results.csv)",
        ),
    ] = Path("eval-results"),
    output_dir: Annotated[
        Optional[Path],
        typer.Option(
            "--output-dir",
            "-o",
            help="Directory to write HTML files (default: results_dir/html)",
        ),
    ] = None,
    open_browser: Annotated[
        bool,
        typer.Option(
            "--open/--no-open",
            help="Open index page in browser after export",
        ),
    ] = False,
) -> None:
    """Export evaluation results as standalone HTML pages.

    Generates an interactive HTML report with:
    - Index page showing all reaction classes sorted by F1 score
    - Individual class pages with source code, metrics, and TP/FP/TN/FN drilldown

    Examples:
        autarch export
        autarch export --results-dir results/
        autarch export -r eval-results -o docs/report --open
    """
    try:
        results_path = Path(results_dir)

        # Check for required files
        if not results_path.exists():
            console.print(f"[red]Error: Results directory '{results_path}' not found[/red]")
            console.print("[yellow]Run 'autarch eval-all --output-dir results' first[/yellow]")
            raise typer.Exit(1)

        if not (results_path / "evaluation_results.csv").exists():
            console.print(f"[red]Error: evaluation_results.csv not found in {results_path}[/red]")
            console.print("[yellow]Run 'autarch eval-all --output-dir results' first[/yellow]")
            raise typer.Exit(1)

        console.print(f"[cyan]Exporting HTML report from {results_path}...[/cyan]")

        # Run export
        output_path = export_html(
            results_dir=results_path,
            output_dir=output_dir,
            verbose=True,
        )

        console.print(f"\n[green]✓ HTML report exported to {output_path}[/green]")
        console.print(f"  Index: {output_path / 'index.html'}")

        # Open in browser if requested
        if open_browser:
            import webbrowser
            console.print("[cyan]Opening in browser...[/cyan]")
            webbrowser.open(str(output_path / "index.html"))

    except Exception as e:
        console.print(f"[red]Error exporting HTML: {e}[/red]")
        raise typer.Exit(1)


@app.command("research-class")
def research_class(
    reaction_class: Annotated[
        str,
        typer.Argument(help="Reaction class name to research (e.g., 'Kinase')"),
    ],
    provider: Annotated[
        Optional[str],
        typer.Option(
            "--provider",
            "-p",
            help="Deep research provider to use (e.g., openai, falcon, perplexity)",
        ),
    ] = None,
    model: Annotated[
        Optional[str],
        typer.Option("--model", "-m", help="Provider-specific model override"),
    ] = None,
    provider_param: Annotated[
        Optional[List[str]],
        typer.Option(
            "--param",
            help="Provider parameter as key=value; may be repeated",
        ),
    ] = None,
    output_dir: Annotated[
        Path,
        typer.Option(
            "--output-dir",
            "-o",
            help="Directory for raw reports and extracted evidence",
        ),
    ] = Path("docs/research/classification"),
    cache_dir: Annotated[
        str,
        typer.Option(
            "--cache-dir",
            "-d",
            help="GO cache directory used to enrich the research prompt",
        ),
    ] = "cache",
    reference_cache_dir: Annotated[
        Path,
        typer.Option(
            "--reference-cache-dir",
            help="Cache directory for fetched references during validation",
        ),
    ] = Path("references_cache"),
    evidence_count: Annotated[
        int,
        typer.Option(
            "--evidence-count",
            help="Maximum number of exact literature snippets to request",
            min=2,
            max=5,
        ),
    ] = 3,
    apply_docstring: Annotated[
        bool,
        typer.Option(
            "--apply-docstring",
            help="Apply the validated evidence block to the classifier docstring",
        ),
    ] = False,
) -> None:
    """Run deep research for one classifier and extract exact evidence snippets."""
    try:
        metadata = collect_reaction_class_metadata(reaction_class, cache_dir=cache_dir)
        prompt = build_research_prompt(metadata, evidence_count=evidence_count)
        params = parse_key_value_options(provider_param)

        console.print(f"[cyan]Researching {metadata.class_name}...[/cyan]")
        result = run_deep_research(
            prompt,
            provider=provider,
            model=model,
            provider_params=params,
        )

        record = extract_research_record(result.markdown, metadata)
        paths = write_research_artifacts(
            metadata,
            prompt,
            result.markdown,
            record,
            output_dir,
        )

        validation_results = validate_research_record(
            record,
            cache_dir=reference_cache_dir,
        )
        write_validation_report(record, validation_results, paths.validation_path)

        console.print(f"[green]✓ Research artifacts written to {paths.class_dir}[/green]")
        console.print(f"[dim]Query: {paths.query_path}[/dim]")
        console.print(f"[dim]Report: {paths.report_path}[/dim]")
        console.print(f"[dim]Evidence: {paths.record_path}[/dim]")
        console.print(f"[dim]Validation: {paths.validation_path}[/dim]")
        console.print(f"[dim]Provider: {getattr(result, 'provider', 'unknown')}[/dim]")
        if getattr(result, "model", None):
            console.print(f"[dim]Model: {result.model}[/dim]")

        if not all_evidence_valid(validation_results):
            console.print("\n[red]Evidence validation failed:[/red]")
            for item in validation_results:
                if item.is_valid:
                    continue
                console.print(f"  • [{item.reference_id}] {item.message}")
            raise typer.Exit(1)

        console.print("[green]✓ All evidence snippets validated successfully[/green]")

        if apply_docstring:
            updated_path = apply_research_record(record)
            console.print(f"[green]✓ Updated classifier docstring in {updated_path}[/green]")

    except ValueError as e:
        console.print(f"[red]Error: {e}[/red]")
        raise typer.Exit(1)
    except RuntimeError as e:
        console.print(f"[red]Error: {e}[/red]")
        raise typer.Exit(1)
    except Exception as e:
        console.print(f"[red]Error running deep research: {e}[/red]")
        raise typer.Exit(1)


@app.command("validate-class-research")
def validate_class_research(
    reaction_class: Annotated[
        str,
        typer.Argument(help="Reaction class name whose saved evidence should be validated"),
    ],
    output_dir: Annotated[
        Path,
        typer.Option(
            "--output-dir",
            "-o",
            help="Directory containing saved classifier research artifacts",
        ),
    ] = Path("docs/research/classification"),
    cache_dir: Annotated[
        str,
        typer.Option(
            "--cache-dir",
            "-d",
            help="GO cache directory used to locate the classifier metadata",
        ),
    ] = "cache",
    reference_cache_dir: Annotated[
        Path,
        typer.Option(
            "--reference-cache-dir",
            help="Cache directory for fetched references during validation",
        ),
    ] = Path("references_cache"),
) -> None:
    """Validate an existing saved classifier evidence file."""
    try:
        metadata = collect_reaction_class_metadata(reaction_class, cache_dir=cache_dir)
        paths = artifact_paths(metadata, output_dir)
        if not paths.record_path.exists():
            console.print(f"[red]Error: No saved evidence file at {paths.record_path}[/red]")
            raise typer.Exit(1)

        record = load_research_record(paths.record_path)
        validation_results = validate_research_record(
            record,
            cache_dir=reference_cache_dir,
        )
        write_validation_report(record, validation_results, paths.validation_path)

        if not all_evidence_valid(validation_results):
            console.print("[red]Saved evidence has validation issues:[/red]")
            for item in validation_results:
                if item.is_valid:
                    continue
                console.print(f"  • [{item.reference_id}] {item.message}")
            raise typer.Exit(1)

        console.print("[green]✓ All evidence snippets validated successfully[/green]")
        console.print(f"[dim]Validation report: {paths.validation_path}[/dim]")

    except ValueError as e:
        console.print(f"[red]Error: {e}[/red]")
        raise typer.Exit(1)
    except RuntimeError as e:
        console.print(f"[red]Error: {e}[/red]")
        raise typer.Exit(1)
    except Exception as e:
        console.print(f"[red]Error validating saved research: {e}[/red]")
        raise typer.Exit(1)


@app.command("apply-class-research")
def apply_class_research(
    reaction_class: Annotated[
        str,
        typer.Argument(help="Reaction class name whose saved evidence should update the docstring"),
    ],
    output_dir: Annotated[
        Path,
        typer.Option(
            "--output-dir",
            "-o",
            help="Directory containing saved classifier research artifacts",
        ),
    ] = Path("docs/research/classification"),
    cache_dir: Annotated[
        str,
        typer.Option(
            "--cache-dir",
            "-d",
            help="GO cache directory used to locate the classifier metadata",
        ),
    ] = "cache",
    reference_cache_dir: Annotated[
        Path,
        typer.Option(
            "--reference-cache-dir",
            help="Cache directory for fetched references during validation",
        ),
    ] = Path("references_cache"),
) -> None:
    """Apply saved, validated classifier evidence to the Python docstring."""
    try:
        metadata = collect_reaction_class_metadata(reaction_class, cache_dir=cache_dir)
        paths = artifact_paths(metadata, output_dir)
        if not paths.record_path.exists():
            console.print(f"[red]Error: No saved evidence file at {paths.record_path}[/red]")
            raise typer.Exit(1)

        record = load_research_record(paths.record_path)
        validation_results = validate_research_record(
            record,
            cache_dir=reference_cache_dir,
        )
        write_validation_report(record, validation_results, paths.validation_path)

        if not all_evidence_valid(validation_results):
            console.print("[red]Refusing to apply unvalidated evidence:[/red]")
            for item in validation_results:
                if item.is_valid:
                    continue
                console.print(f"  • [{item.reference_id}] {item.message}")
            raise typer.Exit(1)

        updated_path = apply_research_record(record)
        console.print(f"[green]✓ Updated classifier docstring in {updated_path}[/green]")

    except ValueError as e:
        console.print(f"[red]Error: {e}[/red]")
        raise typer.Exit(1)
    except RuntimeError as e:
        console.print(f"[red]Error: {e}[/red]")
        raise typer.Exit(1)
    except Exception as e:
        console.print(f"[red]Error applying saved research: {e}[/red]")
        raise typer.Exit(1)


def main():
    """Main entry point for the CLI."""
    app()


if __name__ == "__main__":
    main()
