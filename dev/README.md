# Development Scripts

This folder is for ad hoc method development and tuning experiments.

## Scope

- Fast iteration scripts for evolving algorithms (for example getPatches tuning).
- Benchmark runs and exploratory diagnostics.
- Not part of package unit testing or release validation.

## Suggested Layout

- getPatches_tuning_sandbox.R: main ad hoc tuning script.
- tune_stem_cell_niche_patches.R: parameter sweep and biological/geometric
  diagnostics for long crypt-villus patches in vignette Scenario 2. Performance
  is summarized directly by two biological metrics;
  near-optimal runs are in the top 10% for both. An interpretable decision tree
  identifies parameter rules associated with those results. Independent
  configurations run in parallel on macOS/Linux; set `SPACEMOSAIC_WORKERS` to
  control the number of worker processes.
- moran_toy_example.R: minimal reproducible `patchDE()` to `moranTest()`
  workflow with a spatially structured gene and a noise control. It prints the
  gene-by-patch results, writes `dev/runs/moran_toy_example.pdf`, and performs
  lightweight smoke checks. Set `SPACEMOSAIC_MORAN_PLOT` to use a different
  plot path; formal regression coverage remains in `tests/testthat/`.
- limmaDE_toy_report.Rmd: concise report explaining the research
  question, rationale, simulation and real within-patch comparison of
  `hastyDE()` and `limmaDE()`. The rendered HTML is written to
  `dev/runs/limmaDE_toy_report.html`.
- runs/: optional outputs from local runs (CSV summaries, plots, notes).

## Usage Notes

- Run scripts from the repository root.
- Keep scripts reproducible enough to re-run, but lightweight for iteration speed.
- When logic stabilizes and should be regression-tested, promote it to formal tests.
