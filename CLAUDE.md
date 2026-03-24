# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

The goal of this project is to use AI and coding tools to generative declarative programs/rules
for classifying specific reactions into reactions classes (e.g. hydrolase, oxidoreductase) based on
the chemical structure of the reactants and products.

A good place to start is

 * src/autarch/ontology/kinase.py

## Best Practices

### Keep it simple and declarative

Use the pattern matching framework where appropriate

Avoid hardcoding ad-hoc exceptions and rules. It is better to be explainable and
parsimonious than to cover every edge case.

## Workflows (via justfile)

### Cache Management
The system requires cached data from GO and RHEA databases for evaluation.

**Production cache** (in `cache/`):
- `just cache-go-all` - Cache all GO molecular function terms (~10k terms)
- `just cache-rhea` - Cache RHEA reactions with GO mappings (~4800 master IDs)
- `just cache-chebi` - Cache ChEBI SMILES for reactions

**Test cache** (in `test-cache/`):
- `just cache-go-test` - Cache test GO terms (hydrolase & oxidoreductase)
- `just cache-rhea-limited` - Cache 100 RHEA reactions for testing
- `just cache-all-test` - Cache all test data

### Evaluation
- `just eval-all` - Evaluate all 62 reaction classes with summary statistics
- `just eval-all-stats` - Evaluate with detailed output to `eval-results/`
- `just eval-hydrolase` - Evaluate Hydrolase classifier only
- `uv run autarch eval Kinase` - Evaluate specific reaction class

### Testing and Quality
- `just test` - Run all tests (pytest, mypy, format checks)
- `just pytest` - Run Python tests only
- `just doctest` - Run doctests in source files
- `just mypy` - Run type checking
- `just format` - Run ruff linting/formatting checks

### Validation
- `just validate-go-names` - Validate GO term name consistency
- `just summarize-classes` - Generate summary table of all reaction classes

### CLI Commands
- `uv run autarch classify reaction.yaml` - Classify a reaction from YAML
- `uv run autarch list-classes` - List available reaction classes
- `uv run autarch cache-go --terms GO:0016787` - Cache specific GO terms
- `uv run autarch cache-rhea` - Cache all RHEA reactions
- `uv run autarch eval Hydrolase` - Evaluate classifier performance

## Cache Directory Structure

```
cache/
├── go_terms.jsonl           # GO molecular function terms with ancestors/RHEA mappings
├── rhea_reactions.jsonl     # RHEA reactions (master IDs only, with GO mappings)
├── chebi_smiles.json        # ChEBI ID → SMILES mapping
└── rhea_tsv/                # Raw TSV files from RHEA FTP
    ├── rhea-reaction-smiles.tsv
    ├── rhea2ec.tsv
    ├── chebiId_name.tsv
    └── rhea-reactions.txt.gz

eval-results/                # Evaluation output (after `just eval-all-stats`)
├── evaluation_results.csv   # Per-class metrics
├── detailed_predictions.csv # All predictions with explanations
├── summary_statistics.txt   # Overall metrics
└── <classname>/             # Per-class false positives/negatives
    ├── false_positives.txt
    └── false_negatives.txt
```

### RHEA ID Convention
- **Master IDs only**: The cache stores direction-neutral RHEA IDs (divisible by 4)
- GO terms map to master IDs (e.g., `RHEA:10000` not `RHEA:10001`)
- Each RheaTerm includes `go_terms` list directly for efficient evaluation

## High-Level Architecture

### Core Classification System
The project implements a plugin-based reaction classification system:

1. **ReactionClass (ABC)** (`src/autarch/ontology/reaction.py`)
   - Abstract base class for all reaction classifiers
   - Subclasses implement `check_membership_impl()` to determine if a reaction belongs to their class

2. **Concrete Classifiers** (`src/autarch/ontology/`)
   - `Hydrolase` - Classifies hydrolysis reactions (EC 3.x.x.x)
   - `Oxidoreductase` - Classifies oxidation-reduction reactions (EC 1.x.x.x)
   - `HydrolaseSulfurNitrogen` - Specialized hydrolase for S-N bonds
   - Each uses RDKit SMARTS patterns to detect characteristic reaction patterns

3. **ReactionClassifier** (`src/autarch/classifier.py`)
   - Main orchestrator that discovers all ReactionClass subclasses
   - Runs reactions through all available classifiers
   - Returns classification results with explanations

### Data Flow
1. **Input**: Reactions defined as YAML with SMILES strings for reactants/products
2. **Processing**: 
   - SMILES → RDKit Mol objects
   - Pattern matching using SMARTS
   - Bond difference analysis between reactants and products
3. **Output**: Classification results with boolean membership and explanations

### ETL and Evaluation Pipeline
- **GO ETL** (`etl/etl.py`): Fetches GO molecular function terms from OAK/Ubergraph
- **RHEA TSV ETL** (`etl/rhea_tsv_etl.py`): Loads reactions from RHEA FTP TSV files
  - Uses master IDs only (direction-neutral, divisible by 4)
  - Includes GO term mappings directly in each RheaTerm
  - Handles polymer reactions with (n) stoichiometry notation
- **Evaluation** (`evaluation.py`): Compares classifier predictions against GO ground truth
  - Only evaluates reactions with GO mappings (proper ground truth)
  - Uses GO term hierarchy for positive/negative classification
- **Caching**: JSONL serialization for offline processing

### Data Models (`datamodel.py`)
- `Stoi`: Chemical species with SMILES and stoichiometry
- `Reaction`: Left and right participants (reactants/products)
- `ClassificationResult`: Boolean result with explanation
- `GoTerm`: GO term with RHEA/EC mappings
- `RheaTerm`: RHEA reaction with metadata

## Important Development Practices

### Test-Driven Development
- Write tests BEFORE implementing features
- Don't create mock tests to bypass failures
- Keep trying until functionality works correctly
- Extensive use of doctests for documentation and testing

### Code Style
- Use docstrings for all functions and classes
- AVOID try/except blocks unless interfacing with external systems
- For deterministic code, exceptions indicate logic errors that need fixing
- Use type hints throughout

### Testing Strategy
- Unit tests for individual classifiers
- Integration tests with real RHEA data
- Parametrized tests for multiple test cases
- Doctests for usage examples

## Technology Stack
- **Python 3.10+** with `uv` for dependency management
- **RDKit** for cheminformatics (SMILES, SMARTS, molecule manipulation)
- **Typer** for CLI interface
- **Pydantic** for data validation
- **OAKlib** for ontology operations
- **pytest** for testing
- **mypy** for type checking
- **ruff** for linting/formatting
- **just** as command runner

## Key Configuration Files
- `pyproject.toml` - Project configuration and dependencies
- `justfile` - Main command recipes
- `project.justfile` - Project-specific recipes (cache, eval commands)
- `uv.lock` - Locked dependency versions

## Development Workflow
1. Use `uv add` for new dependencies
2. Run commands through `just` or `uv run`
3. Dynamic versioning from git tags
4. Documentation at https://cmungall.github.io/autarch
