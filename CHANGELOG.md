# Changelog

All notable changes to `ddd` will be documented in this file.
Format: [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
Versioning: [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3.0] - 2026-10-10

### Added
- **Two-step GMM combination across comparison cohorts** (`gmm`): for each
  (g,t), the estimator combines every identified not-yet-treated cohort (plus
  the never-enabled group) using inverse-variance weights built from the
  estimators' influence functions. `control(notyet) gmm` uses the full set of
  not-yet-treated cohorts.
- **Cluster bootstrap**: `cluster()` now resamples clusters (via
  `bsample, cluster() idcluster()`) instead of units for `reps()`.

### Changed
- `_ddd_cell` can output its influence function (`savif()`); `_ddd_estall`
  collects them and performs the GMM combination.

### Verified
- Multi-period (true ATT = 1.5): `control(notyet) gmm` overall ≈ 1.515;
  `control(never)` ≈ 1.513; cluster bootstrap (40 reps, 60 clusters) SE ≈ 0.075.

## [0.2.0] - 2026-10-10

### Added
- **Multinomial-logit generalized propensity score** (`psmethod(mnl)`, default),
  estimated for the four (S,Q) cells; automatic fallback to a linear-probability
  working model (`psmethod(lpm)`) if `mlogit` fails to converge (iterate cap 30).
- **Not-yet-treated comparison cohorts** (`control(notyet)`): for each (g,t) the
  estimator uses the earliest not-yet-enabled cohort (> t), falling back to the
  never-enabled group when none is identifiable. `control(never)` remains default.
- **Sup-t uniform confidence bands** (`cband`): unit-level bootstrap over cells,
  multiplier-style sup-t critical value; band bounds posted in `e(band)` and the
  critical value in `e(cband_crit)`.
- **Multi-period example** (`examples/ddd_simdata_multi.do`) with staggered enabling.

### Changed
- `_ddd_estall` now threads `control()` and `psmethod()` and selects the comparison
  cohort per (g,t).
- Bootstrap stores per-cell draws and reports a band matrix.

### Verified
- Two-period: ATT ≈ 1.97 (true 2) for dr / ra / ipw / mnl / lpm.
- Multi-period (true ATT = 1.5): control(notyet) overall ≈ 1.516, control(never) ≈ 1.513.
- `cband` critical values > 1.96 (e.g. 2.37 for the 5-cell multi-period design).

## [0.1.0] - 2026-10-09

### Added
- **Core**
  - `ddd` main command: parses and validates a triple-differences design,
    builds the internal time index, maps `enable()` to cohort ids, reshapes to
    unit level, estimates all `ATT(g,t)`, aggregates, and posts `e()` results.
  - `_ddd_cell`: RA / IPW / DR DDD estimator for one cell, implementing the
    three-comparison DR form (Eq. 3.5) of Ortiz-Villavicencio & Sant'Anna (2025).
  - `_ddd_estall`: loops over cohorts and post periods.
  - `_ddd_display`: formatted output, decoupled from estimation.
- **Inference**
  - influence-function (analytic) standard errors for each `ATT(g,t)`;
  - optional unit-level bootstrap (`reps()`) for the overall ATT (SE + percentile CI).
- **Design support**
  - two-period and multi-period designs;
  - never-enabled comparison group;
  - time-invariant covariates;
  - eligibility partition (`eligible()`).
- **Interface**
  - `e(att)`, `e(se_att)`, `e(ci_lower)`, `e(ci_upper)`, `e(ncells)`;
  - `e(att_gt)`, `e(es)` matrices;
  - event-study plot (suppress with `nograph`).
- **Packaging / docs**
  - `ddd.sthlp`, `ddd.pkg`, `stata.toc`, `README.md`, `CHANGELOG.md`;
  - `examples/ddd_simdata.do`, `examples/ddd_example.do`.

### Notes
- v0.1 generalized propensity score uses a linear-probability working model
  with trimming; a multinomial-logit version is planned for v0.2.
- No cluster bootstrap or uniform bands yet.

### Verified
- **Cross-validated against R `triplediff` 0.2.4**: on the same simulated data,
  the DR-DDD point estimate matches exactly (ATT = 2.0317).
- Smoke-tested on **Stata 19.0 MP** (see `test_run.do`). On the simulated
  two-period design with true `ATT(2,2) = 2`:
  - `method(dr)` → ATT = 2.0317 (SE 0.1056);
  - `method(ra)` with `reps(50)` → ATT = 2.0317, bootstrap SE = 0.0870,
    95% CI = [1.8401, 2.1630];
  - `method(ipw)` → ATT = 2.0317 (SE 0.1189);
  - `covariates(eligible) method(dr)` → ATT = 2.0317.

### Known limitations
- Requires a never-enabled comparison group.
- Aggregation weights are cohort-size based (simpler than the paper's).
