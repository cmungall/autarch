## Add your own just recipes here. This is imported by the main justfile.

# ============== Cache management ==============

# Cache GO terms (hydrolase and oxidoreductase for testing)
[group('cache')]
cache-go-test:
    mkdir -p test-cache
    uv run autarch cache-go --terms GO:0016787,GO:0016491 --cache-dir test-cache

# Cache all GO molecular function terms
[group('cache')]
cache-go-all:
    uv run autarch cache-go

# Cache all RHEA reactions (uses TSV from RHEA FTP)
[group('cache')]
cache-rhea:
    uv run autarch cache-rhea

# Cache specific RHEA reactions for testing
[group('cache')]
cache-rhea-test:
    mkdir -p test-cache
    uv run autarch cache-rhea --ids RHEA:18777,RHEA:18778,RHEA:18779,RHEA:18780 --cache-dir test-cache

# Cache limited RHEA reactions for testing
[group('cache')]
cache-rhea-limited:
    mkdir -p test-cache
    uv run autarch cache-rhea --limit 100 --cache-dir test-cache

# Cache limited RHEA reactions for testing with more reactions
[group('cache')]
cache-rhea-test-large:
    mkdir -p test-cache
    uv run autarch cache-rhea --limit 1000 --cache-dir test-cache

# Cache CHEBI SMILES for production
[group('cache')]
cache-chebi:
    uv run autarch cache-chebi

# Cache ModelSEED bridge artifacts
[group('cache')]
cache-modelseed:
    uv run autarch cache-modelseed

# Summarize cached ModelSEED mapping coverage
[group('analysis')]
summarize-modelseed:
    uv run autarch summarize-modelseed

# Materialize the default ModelSEED OOD benchmark slice and predictions
[group('analysis')]
benchmark-modelseed:
    uv run autarch benchmark-modelseed

# Materialize the EC-backed ModelSEED OOD benchmark slice and predictions
[group('analysis')]
benchmark-modelseed-ec:
    uv run autarch benchmark-modelseed --require-ec

# Cache CHEBI SMILES for testing
[group('cache')]
cache-chebi-test:
    mkdir -p test-cache
    uv run autarch cache-chebi --cache-dir test-cache

# Cache all data (limited for testing)
[group('cache')]
cache-all-test: cache-go-test cache-rhea-limited cache-chebi-test

# Clean cache directory
[group('cache')]
cache-clean:
    rm -rf cache/

# ============== Evaluation ==============

# Evaluate Hydrolase class against test cached data
[group('evaluation')]
eval-hydrolase-test:
    uv run autarch eval Hydrolase --cache-dir test-cache

# Evaluate Oxidoreductase class against test cached data
[group('evaluation')]
eval-oxidoreductase-test:
    uv run autarch eval Oxidoreductase --cache-dir test-cache

# Evaluate Hydrolase class against production cache
[group('evaluation')]
eval-hydrolase:
    uv run autarch eval Hydrolase

# Evaluate Oxidoreductase class against production cache
[group('evaluation')]
eval-oxidoreductase:
    uv run autarch eval Oxidoreductase

# Run all test evaluations
[group('evaluation')]
eval-all-test: eval-hydrolase-test eval-oxidoreductase-test

# Run all production evaluations with macro/micro statistics
[group('evaluation')]
eval-all:
    uv run autarch eval-all

eval-all-stats:
    uv run autarch eval-all --visualize -s -o eval-results

# Run all production evaluations with verbose output
[group('evaluation')]
eval-all-verbose:
    uv run autarch eval-all --verbose

# ============== Validation ==============

# Validate GO term mappings in all reaction classes using OAK
[group('validation')]
validate-go-mappings:
    uv run python -m autarch.validation.go_mapping_validator

# Validate GO term name consistency for all reaction classes
[group('validation')]
validate-go-names:
    uv run autarch validate-go-names

[group('validation')]
validate-go-cache:
    uv run autarch validate-go-cache

[group('validation')]
validate-go-names-strict:
    uv run autarch validate-go-names-strict

# Generate comprehensive summary table of all reaction classes
[group('validation')]
summarize-classes:
    uv run autarch summarize-classes

# ============== Research ==============

# Optional env vars:
#   ARCTURUS_RESEARCH_PROVIDER=openai
#   ARCTURUS_RESEARCH_MODEL=<provider-model-name>
[group('research')]
research-class CLASS:
    uv run --python 3.12 autarch research-class {{CLASS}} ${ARCTURUS_RESEARCH_PROVIDER:+--provider $ARCTURUS_RESEARCH_PROVIDER} ${ARCTURUS_RESEARCH_MODEL:+--model $ARCTURUS_RESEARCH_MODEL}

[group('research')]
research-class-apply CLASS:
    uv run --python 3.12 autarch research-class {{CLASS}} --apply-docstring ${ARCTURUS_RESEARCH_PROVIDER:+--provider $ARCTURUS_RESEARCH_PROVIDER} ${ARCTURUS_RESEARCH_MODEL:+--model $ARCTURUS_RESEARCH_MODEL}

[group('research')]
research-validate CLASS:
    uv run --python 3.12 autarch validate-class-research {{CLASS}}

[group('research')]
research-apply CLASS:
    uv run --python 3.12 autarch apply-class-research {{CLASS}}

# ============== Complexity Analysis ==============

# Analyze classifier complexity (cyclomatic complexity, ChEBI IDs, etc.)
[group('analysis')]
complexity:
    uv run python -c "from autarch.complexity import analyze_all_classifiers, print_complexity_report; print_complexity_report(analyze_all_classifiers())"

# Analyze classifier complexity with verbose output
[group('analysis')]
complexity-verbose:
    uv run python -c "from autarch.complexity import analyze_all_classifiers, print_complexity_report; print_complexity_report(analyze_all_classifiers(), verbose=True)"

# Save complexity report to eval-results/
[group('analysis')]
complexity-report:
    uv run python -c "from autarch.complexity import analyze_all_classifiers, print_complexity_report, save_complexity_report; from pathlib import Path; r = analyze_all_classifiers(); print_complexity_report(r); save_complexity_report(r, Path('eval-results'))"

# Save curator-facing objective report to eval-results/
[group('analysis')]
curation-report:
    uv run python -c "from pathlib import Path; import pandas as pd; from autarch.complexity import analyze_all_classifiers; from autarch.curation_objective import save_curation_objective_report; eval_df = pd.read_csv('eval-results/evaluation_results.csv'); save_curation_objective_report(eval_df, analyze_all_classifiers(), Path('eval-results'))"

# Save the manuscript RHEA/GO/EC overlap figure
[group('analysis')]
paper-figure-overlap:
    uv run python -c "from pathlib import Path; from autarch.visualization import create_resource_overlap_figure; create_resource_overlap_figure(Path('docs/manuscript/figures/resource_overlap.png'))"

# Save the manuscript classifier architecture figure
[group('analysis')]
paper-figure-architecture:
    uv run python -c "from pathlib import Path; from autarch.visualization import create_classifier_architecture_figure; create_classifier_architecture_figure(Path('docs/manuscript/figures/classifier_architecture.png'))"

# Save the manuscript complexity-vs-F1 figure
[group('analysis')]
paper-figure-complexity:
    uv run python -c "from pathlib import Path; import pandas as pd; from autarch.complexity import analyze_all_classifiers; from autarch.visualization import create_complexity_vs_performance_static; eval_df = pd.read_csv('eval-results/evaluation_results.csv'); create_complexity_vs_performance_static(eval_df, analyze_all_classifiers(), Path('docs/manuscript/figures/complexity_vs_f1.png'))"

# Save the manuscript system overview figure
[group('analysis')]
paper-figure-overview:
    uv run python -c "from pathlib import Path; from autarch.visualization import create_system_overview_figure; create_system_overview_figure(Path('docs/manuscript/figures/system_overview.png'))"

# Save the manuscript F1-vs-support figure
[group('analysis')]
paper-figure-support:
    uv run python -c "from pathlib import Path; import pandas as pd; from autarch.visualization import create_f1_vs_support_static; eval_df = pd.read_csv('eval-results/evaluation_results.csv'); create_f1_vs_support_static(eval_df, Path('docs/manuscript/figures/f1_vs_support.png'))"

# Save the manuscript concrete reaction-example figure
[group('analysis')]
paper-figure-conflicts:
    uv run python -c "from pathlib import Path; from autarch.visualization import create_scope_conflict_examples_figure; create_scope_conflict_examples_figure(Path('docs/manuscript/figures/scope_conflict_examples.png'))"

# Regenerate all main manuscript figures
[group('analysis')]
paper-figures: paper-figure-overlap paper-figure-architecture paper-figure-overview paper-figure-support paper-figure-complexity paper-figure-conflicts

# ============== Documentation ==============

# Run MkDocs on the repo-specific local port.
[group('documentation')]
serve-docs:
    @echo "http://{{mkdocs-host}}:{{mkdocs-port}}/autarch/"
    uv run mkdocs serve -a {{mkdocs-host}}:{{mkdocs-port}}

# Generate HTML documentation for all reaction classes
[group('documentation')]
html:
    uv run autarch export --results-dir eval-results --output-dir docs/report --open
