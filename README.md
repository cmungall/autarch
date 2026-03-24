
# Autarch

**Automated Agentic Reaction Classification Hierarchy**: an ontology- and pattern-driven system for classifying enzymatic reactions.

## Project Overview

Autarch builds upon the concepts introduced in [C3P (CHEBI Classification Programs)](https://github.com/chemkg/c3p) but extends the focus from compound classification to **reaction classification**. While C3P uses LLMs to generate classification programs for chemical compounds based on CHEBI classes, Autarch specializes in understanding and categorizing enzymatic reactions using:

- **Gene Ontology (GO) molecular function terms** for enzyme classification
- **RHEA reaction database** for comprehensive reaction coverage
- **RDKit cheminformatics** for reaction pattern detection
- **SMARTS pattern matching** for identifying reaction mechanisms

### Key Goals

1. **Automated Enzyme Function Prediction**: Classify reactions according to EC (Enzyme Commission) numbers and GO molecular function terms
2. **Reaction Mechanism Understanding**: Detect bond changes, functional group transformations, and reaction centers
3. **Ontology-Driven Classification**: Leverage existing biological ontologies (GO, RHEA, EC) for standardized classification
4. **Extensible Plugin Architecture**: Easy addition of new reaction classifiers through a simple abstract base class pattern
5. **Evaluation Framework**: Compare classifier predictions against curated databases for accuracy assessment

### Current Capabilities

- **80+ reaction classifiers** spanning major EC/GO functional groups.
- **Declarative pattern DSL** for explainable reaction-class rules.
- **Polymer-aware reactions** via `(n)`, `(n+1)`, `(n-1)` notation support.
- **Multi-format input**: YAML reactions and cached RHEA reactions.

## Quick Start

```bash
# Install with uv
uv pip install autarch

# Classify a reaction from YAML
uv run autarch classify --yaml reaction.yaml

# List available reaction classes
uv run autarch list-classes

# Cache data for evaluation (required first time)
just cache-go-all      # Cache GO terms (~10k terms)
just cache-rhea        # Cache RHEA reactions with GO mappings
just cache-chebi       # Cache ChEBI SMILES

# Evaluate all classifiers
just eval-all          # Evaluate all configured classes (82 in latest snapshot)
just eval-all-stats    # With detailed output to eval-results/

# Evaluate specific classifier
uv run autarch eval Hydrolase
```

## Documentation Website

[https://ai4curation.github.io/autarch](https://ai4curation.github.io/autarch)

## Cache and Evaluation Data

The system uses cached data from GO and RHEA databases:

```
cache/
├── go_terms.jsonl           # GO molecular function terms with RHEA mappings
├── rhea_reactions.jsonl     # RHEA reactions with GO mappings
├── chebi_smiles.json        # ChEBI ID → SMILES mapping
└── rhea_tsv/                # Raw TSV files from RHEA FTP

eval-results/                # Evaluation output
├── evaluation_results.csv   # Per-class metrics (precision, recall, F1)
├── detailed_predictions.csv # All predictions with explanations
└── <classname>/             # Per-class error analysis
```

**Key design decisions:**
- Cache entries are direction-neutral and include GO mappings for efficient evaluation.
- Each cached RHEA reaction includes GO terms directly for lookup/evaluation.
- Evaluation uses reactions with GO mappings as ground truth.

## Architecture

Autarch implements a plugin-based classification system where each reaction type is defined as a Python class inheriting from `ReactionClass`:

```python
from autarch.ontology import ReactionClass

class MyEnzymeClass(ReactionClass):
    """Custom enzyme reaction classifier"""
    
    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        # Implement pattern matching logic using RDKit
        # Return ClassificationResult with boolean and explanation
        pass
```

The system automatically discovers all subclasses and applies them to input reactions, providing multi-class classification with detailed explanations.

## Use Cases

- **Drug Discovery**: Understanding enzyme-substrate interactions and metabolic pathways
- **Synthetic Biology**: Designing enzymatic cascades and predicting reaction outcomes
- **Bioinformatics**: Annotating enzyme functions from reaction data
- **Chemical Education**: Teaching reaction mechanisms through computational analysis

## Development

```bash
# Clone the repository
git clone https://github.com/ai4curation/autarch
cd autarch

# Install with uv
uv sync

# Run tests
just test

# Run specific classifier tests
uv run pytest tests/test_classifier.py::test_hydrolase_classification

# Type checking
just mypy

# Format code
just format
```

## Contributing

We welcome contributions! Areas of particular interest:

1. **New Reaction Classifiers**: Implement classifiers for transferases, ligases, lyases, and isomerases
2. **Pattern Libraries**: Expand SMARTS patterns for better reaction coverage
3. **Integration**: Connect with additional reaction databases (KEGG, MetaCyc, etc.)
4. **ML Integration**: Incorporate machine learning models for pattern discovery

## Credits

This project builds on ideas from [C3P](https://github.com/chemkg/c3p) and uses the [monarch-project-copier](https://github.com/monarch-initiative/monarch-project-copier) template.
