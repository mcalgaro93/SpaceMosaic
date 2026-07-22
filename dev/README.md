# Development Scripts

This folder is for ad hoc method development and tuning experiments.

## Scope

- Fast iteration scripts for evolving algorithms (for example getPatches tuning).
- Benchmark runs and exploratory diagnostics.
- Not part of package unit testing or release validation.

## Suggested Layout

- getPatches_tuning_sandbox.R: main ad hoc tuning script.
- runs/: optional outputs from local runs (CSV summaries, plots, notes).

## Usage Notes

- Run scripts from the repository root.
- Keep scripts reproducible enough to re-run, but lightweight for iteration speed.
- When logic stabilizes and should be regression-tested, promote it to formal tests.
