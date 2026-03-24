# Lean Spec Scaffold

This directory holds a first-pass Lean 4 scaffold for the symbolic reaction DSL.

The current intent is narrow:

- Treat Lean as a reference semantics for the symbolic subset of reaction
  matching.
- Keep RDKit, SMILES parsing, and dataset ETL in Python.
- Start with the `Kinase` pattern fragment before expanding to other classes.

## Files

- `Autarch/IR.lean`: symbolic IR mirroring `src/autarch/lean_spec.py`
- `Autarch/Match.lean`: executable matcher with current greedy semantics
- `Autarch/Kinase.lean`: kinase patterns and simple proof targets
- `Main.lean`: tiny entry point for `lake run`

## Status

- This is a scaffold, not a validated Lean build.
- The semantics currently mirror the Python optional-first greedy matcher,
  because that is what the repository executes today.
- The Python bridge is testable now via `tests/test_lean_spec.py`.

## Once Lean Is Installed

```bash
cd lean-spec
lake build
lake run
```

The Python bridge that emits JSON IR lives in `src/autarch/lean_spec.py`.
