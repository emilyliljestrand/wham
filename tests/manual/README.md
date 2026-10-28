# Historical Refit Checks

Run from the repository root in separate Rscript processes:

```text
Rscript tests/manual/refit-releases.R generate v1.0.0
Rscript tests/manual/refit-releases.R refit v1.0.0
Rscript tests/manual/refit-releases.R refit-defaults v1.0.0
Rscript tests/manual/refit-releases.R run v1.0.1 v1.0.2 v1.0.3 v1.0.4 v1.0.5 v1.0.6
Rscript tests/manual/refit-releases.R summarize
```

Set `WHAM_CURRENT_LIBRARY` to the library directory containing the current WHAM
installation if it is not on R's library path. The refit stages load the installed
native model and the workspace versions of `R/migrate_wham_input.R` and
`R/refit_wham.R`; the installed R functions are not modified.

`generate` requires Git, the local release tag, TMB, dplyr, tidyr, and a working R
C++ toolchain. It extracts the release, compiles an unchanged build copy of its
model source, and sources its historical R functions. It does not install an old
WHAM package or replace the current one. The toolchain is current, not a
recreation of the release's original R/TMB environment.

The harness evaluates the input-construction expressions from each tag's own
example-1 test, retaining that release's selectivity and likelihood settings.
The four original example-1 configurations cover recruitment-only and full
numbers-at-age random effects, each with multinomial and logistic-normal
age compositions. Generated likelihoods are checked against the release's
archived test results. Five additional cases cover SCAA, recruitment AR1 over
years, abundance AR1 over ages, abundance AR1 over years, and two-dimensional
abundance AR1. These additional cases have no archived reference optimum.
`refit` checks the migrated likelihood at the saved fixed-parameter estimates
before optimization, as well as after refitting. It disables retrospectives, standard errors, OSA,
reference points, and posterior sampling to isolate the migration and likelihood
check. `refit-defaults` leaves the saved call unchanged to exercise default
behavior. Warnings are retained separately; NLL agreement does not imply a
well-conditioned Hessian or successful uncertainty calculations.

Generation also evaluates each historical DLL at the saved fixed and latent
parameter values and saves a `_reference.rds` report for biological comparisons.
Some early fits cache reports that are inconsistent with their parameter lists,
so these checks use the historical DLL as the reference rather than assuming all
cached report fields describe the same parameter state. `reference <tag>` can
create these reports for existing fixtures without refitting. The original saved
NLL remains the benchmark for the marginal-likelihood checks.

Outputs are written under `oldfitobjects/release_refits/<tag>/`, configurable
with `WHAM_RELEASE_RESULTS`: provenance, original fits, fresh refits, comparison
text, and CSV summaries. No refit results are reused. The NLL tolerance is `1e-6`
absolute. The regression tests in `test_prepare_wham_input_from_fit.R` also check
the v1.0.0 through v1.0.6 fixtures when present; they skip absent fixtures.
The `run` stage launches separate generation and refit processes, saves their
logs, and stops on the first failure. Generation is repeatable but not cached;
use the `refit` stage alone to rerun migration against already-generated fits.

Other tags may be supplied, but a successful v1.0.0 run alone does not establish
coverage of every release or model configuration. Environmental effects,
selectivity or mortality random effects, weight-based mortality, steepness,
stored projections, and estimated early mortality baselines with age-specific
offsets remain unsupported by the legacy migration.
For the early packed age-composition parameter layout, supported likelihood codes
are 1 (multinomial), 2 (Dirichlet-multinomial), and 7 (logistic-normal with zeros
missing). Historical code 7 becomes current code 5. Stored likelihood components
determine whether to retain the original effective sample sizes (weighted
variance) or use unit effective sample sizes (unweighted variance). The originals
are retained in `input$legacy_age_comp_Neff`. Missing or ambiguous evidence and
other early likelihood codes are explicitly rejected.
The later matrix parameter layout retains its likelihood codes and sample-size
weights independently of whether observation indices need rebuilding.

`summarize` checks all six requested releases by default and writes
`oldfitobjects/release_refits/release_summary.csv`. It requires all nine cases per
release, four archived NLL matches, matching migrated and refitted NLLs, optimizer
convergence, and no migration errors. Warnings are counted separately.

Verified with WHAM 2.1.0.9012 and R 4.4.3: all 54 cases across v1.0.1 through
v1.0.6 pass the `1e-6` absolute NLL tolerance before and after refitting. The
largest refit difference was `1.65e-8`. Five v1.0.1 cases retain singular-Hessian
warnings; v1.0.2 through v1.0.6 produced no refit warnings in this matrix.