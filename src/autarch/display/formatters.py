"""Formatters for displaying classification results."""

import json
from typing import Dict

import yaml
from rich.console import Console
from rich.table import Table

from autarch.datamodel import ClassificationResult

console = Console()


def display_classification_results(
    results: Dict[str, ClassificationResult],
    format: str = "table",
    verbose: bool = False,
    title: str = "Classification Results",
) -> None:
    """Display classification results in the specified format.

    Args:
        results: Dictionary mapping class names to classification results
        format: Output format - "table", "json", or "yaml"
        verbose: Whether to show detailed explanations
        title: Title for the output (used in table format)
    """
    if format == "json":
        output = format_results_as_json(results)
        print(output)
    elif format == "yaml":
        output = format_results_as_yaml(results)
        print(output)
    else:  # table
        table = format_results_as_table(results, verbose, title)
        console.print(table)

        # Show matching classes summary
        matching_classes = [
            name for name, result in results.items() if result.is_member
        ]
        if matching_classes:
            console.print(
                f"\n[green]✓ Reaction matches: {', '.join(matching_classes)}[/green]"
            )
        else:
            console.print("\n[yellow]No matching reaction classes found[/yellow]")


def format_results_as_json(results: Dict[str, ClassificationResult]) -> str:
    """Format classification results as JSON.

    Args:
        results: Dictionary mapping class names to classification results

    Returns:
        JSON string representation
    """
    output = {
        class_name: {"is_member": result.is_member, "explanation": result.explanation}
        for class_name, result in results.items()
    }
    return json.dumps(output, indent=2)


def format_results_as_yaml(results: Dict[str, ClassificationResult]) -> str:
    """Format classification results as YAML.

    Args:
        results: Dictionary mapping class names to classification results

    Returns:
        YAML string representation
    """
    output = {
        class_name: {"is_member": result.is_member, "explanation": result.explanation}
        for class_name, result in results.items()
    }
    return yaml.dump(output, default_flow_style=False)


def format_results_as_table(
    results: Dict[str, ClassificationResult],
    verbose: bool = False,
    title: str = "Classification Results",
) -> Table:
    """Format classification results as a Rich table.

    Args:
        results: Dictionary mapping class names to classification results
        verbose: Whether to show detailed explanations
        title: Title for the table

    Returns:
        Rich Table object
    """
    table = Table(title=title)
    table.add_column("Reaction Class", style="cyan")
    table.add_column("Is Member", style="bold")
    if verbose:
        table.add_column("Explanation", style="dim")

    for class_name, result in results.items():
        is_member_str = "✓ Yes" if result.is_member else "✗ No"
        is_member_style = "green" if result.is_member else "red"

        if verbose:
            table.add_row(
                class_name,
                f"[{is_member_style}]{is_member_str}[/{is_member_style}]",
                result.explanation,
            )
        else:
            table.add_row(
                class_name, f"[{is_member_style}]{is_member_str}[/{is_member_style}]"
            )

    return table


def display_available_classes(reaction_classes: Dict[str, type]) -> None:
    """Display available reaction classes.

    Args:
        reaction_classes: Dictionary mapping class names to class types
    """
    table = Table(title="Available Reaction Classes")
    table.add_column("Class Name", style="cyan")
    table.add_column("Module", style="dim")

    for class_name, class_type in reaction_classes.items():
        module_name = class_type.__module__
        table.add_row(class_name, module_name)

    console.print(table)
    console.print(f"\n[dim]Total classes: {len(reaction_classes)}[/dim]")
