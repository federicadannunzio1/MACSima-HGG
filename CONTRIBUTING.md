# Contributing

## Workflow (GitHub Flow)

1. Start from an up-to-date `main`: `git switch main && git pull`.
2. Create a short-lived branch for one topic: `git switch -c <name>/<topic>`.
3. Commit small, self-contained changes (one theme per commit, message in the imperative:
   "Read clinical metadata from data/clinical.tsv").
4. Push the branch and open a **draft** pull request against `main`.
5. The repository owner reviews and merges (squash). Nobody commits directly to `main`.

A pull request changes one thing. Refactoring and documentation PRs must leave every
result unchanged; analysis changes go in their own PRs and state which numbers change
and why.

## Privacy

The data come from a clinical cohort of minors.

- Identify patients only by `patient_id` (e.g. `N17_135`). Never write names, dates of birth
  or other identifiers in code, comments, commit messages, issues or pull requests.
- Clinical values and sample metadata live in `data/` (ignored by git), never in the code.
- Files under `example/` are synthetic.
- Before pushing, check the diff: `git diff origin/main...HEAD`.

## Checking that results do not change

A refactoring is correct only if the pipeline produces the same outputs as `main`.

1. Run `main` once on the full data and keep its `output/` as the reference.
2. Run your branch on the same data into another folder (`HGG_BASE_DIR`).
3. Compare the two folders table by table (a comparison script will be added with the
   test suite). Expected: identical tables, or differences within the run-to-run noise
   of `main` itself. Any other difference is a bug of the refactoring.

Note: Gaussian-mixture thresholds (step 06) depend on package versions, so the reference
and the branch must run in the same R environment.

## Code style

- R, tidyverse style; paths only through `BASE_DIR`, `DATA_DIR`, `OUT_DIR` from
  `code/utils/constants.R`; no `setwd()`.
- Code, comments and documentation in English.
