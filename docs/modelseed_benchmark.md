# ModelSEED OOD Benchmark

This note summarizes how ModelSEED fits into `autarch` as an
out-of-distribution benchmark source.

## Why ModelSEED

ModelSEED is useful here as an auxiliary reaction corpus rather than as a
replacement for the current `RHEA + GO + ChEBI` benchmark stack.

What it adds:

- Broad reaction coverage outside the current GO-grounded RHEA benchmark
- Reaction-level mappings to Rhea and EC
- Transport and status metadata that make it possible to define cleaner
  evaluation slices

What it does not add:

- GO-based reaction classification labels comparable to the current RHEA
  evaluation set
- Direct ChEBI-first compound normalization

## Current Overlap Snapshot

These counts were computed from the local cache on May 30, 2026.

Reaction-level overlap:

- Total ModelSEED reactions: `43,774`
- ModelSEED reactions with any Rhea alias: `11,623`
- ModelSEED reactions without a Rhea alias: `32,151`

Rhea master-ID overlap after normalizing ModelSEED directional aliases through
`cache/rhea_tsv/rhea-directions.tsv`:

- Distinct Rhea master IDs represented by ModelSEED: `7,637`
- Full local Rhea master universe: `17,783`
- Overlap: `7,473`

Relative to the repo's current cached Rhea benchmark in
`cache/rhea_reactions.jsonl`:

- Cached Rhea masters: `5,336`
- Overlap with ModelSEED-normalized Rhea masters: `4,983`

## Default Benchmark Slice

The `benchmark-modelseed` command materializes a conservative OOD slice from
cached ModelSEED reactions. The default filters are:

- exclude reactions with Rhea aliases
- require `status=OK`
- exclude transport reactions
- require full ChEBI mapping
- exclude reactions with fractional stoichiometry

Command:

```bash
uv run autarch benchmark-modelseed
```

Artifacts are written to `cache/modelseed_benchmark/`:

- `candidates.jsonl`
- `predictions.jsonl`
- `summary.json`

The EC-backed variant is written separately to `cache/modelseed_benchmark_ec/`:

```bash
uv run autarch benchmark-modelseed --require-ec
```

## Default Slice Results

Using the default filters on May 30, 2026:

- Input ModelSEED reactions: `43,774`
- Selected benchmark reactions: `9,477`
- Selected reactions with ECs: `6,448`
- Selected reactions without ECs: `3,029`
- Classifiers run: `488`
- Positive classifier predictions: `7,907`
- Positive rate: `83.4%`

Filtered out before benchmarking:

- Rhea-linked: `11,623`
- Non-`OK` status: `10,522`
- Partial ChEBI mapping: `6,935`
- Transport: `5,198`
- Fractional stoichiometry: `19`

Restricting the same benchmark to EC-backed reactions only:

- Selected reactions: `6,448`
- Classifiers run: `488`
- Positive predictions: `5,343`
- Positive rate: `82.9%`

## What Looks Interesting

The main result is not coverage. It is the hit rate.

Even on the conservative default slice, the classifier stack predicts at least
one class for most OOD ModelSEED reactions. The dominant positive classes are
generic:

- `Oxidoreductase` / `OtherOxidoreductase`: `3,347`
- `Hydrolase`: `1,498`
- `Transferase`: `1,286`
- `Lyase` / `OtherLyase`: about `930`
- `Isomerase` / `OtherIsomerase`: `447`

This makes ModelSEED a useful stress test for classifier permissiveness.

Some predictions look plausible. For example, `rxn00019`
(`2-Nitropropane:oxygen 2-oxidoreductase`) is classified in the
oxidoreductase family.

Some predictions look suspicious. For example:

- `rxn00021` (`benzaldehyde-lyase`) is classified as `Isomerase`
- `rxn00047` (a uronate lyase) is also classified as `Isomerase`

Those cases suggest that OOD evaluation against ModelSEED is likely to be more
useful for finding overgeneralization than for producing a new source of
ground-truth labels.

## Practical Interpretation

ModelSEED fits best as:

- an OOD benchmark corpus
- a coverage-expansion corpus for error analysis
- a source of EC-backed reactions outside the current GO-grounded RHEA set

It does not currently fit as:

- a replacement for the RHEA benchmark
- a GO-labeled evaluation source
- a canonical identifier space for the project

## Recommended Next Steps

- Rank the most suspicious OOD positives by class and EC mismatch
- Compare predicted classes against EC prefixes where available
- Track benchmark precision manually on a small curated slice before using
  aggregate hit rates as a quality signal
