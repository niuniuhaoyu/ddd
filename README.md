# ddd — Triple Differences (DDD) estimators for Stata

**Version 0.3.0**

Implements the regression-adjustment (RA), inverse-probability-weighting (IPW),
and doubly-robust (DR) estimators of group-time average treatment effects
`ATT(g,t)` in triple-differences designs, following
**Ortiz-Villavicencio & Sant'Anna (2025), arXiv:2505.09942**.

## Why this package

The paper shows that two extremely common DDD practices are generally **invalid**:

- **three-way fixed effects** regressions, and
- **"DDD = difference of two DiDs"**

both fail when covariates are needed for identification. The correct DR DDD
estimator is a function of **three** DiD comparisons — against `(g,0)`, `(cmp,1)`
and `(cmp,0)` — where `g` is the treated cohort, `cmp` the comparison cohort,
`1/0` the eligibility partition. In staggered designs, pooling all
not-yet-treated units as controls also biases the estimates.

## Install (local)

```stata
adopath + "D:\CFPS\pkgs\ddd"
```

## Syntax

```stata
ddd depvar [if] [in], unit(varname) enable(varname) eligible(varname) time(varname) ///
    [covariates(varlist) method(ra|ipw|dr) psmethod(mnl|lpm) control(never|notyet) gmm ///
     aggregation(none|cohort|overall|es) reps(#) seed(#) level(#) cband cluster(varname) nograph]
```

- `enable()` — first period the unit's **group** enables the policy (0 = never);
- `eligible()` — 0/1 **population partition** that qualifies for treatment;
- A unit is treated when `time >= enable` **and** `eligible == 1`.

## Quick start

```stata
do examples/ddd_simdata.do
ddd y, unit(id) enable(enable) eligible(eligible) time(time) method(dr)
ddd y, unit(id) enable(enable) eligible(eligible) time(time) method(ra) reps(100) seed(1)
```

## Stored results

- `e(att)`, `e(se_att)`, `e(ci_lower)`, `e(ci_upper)`, `e(ncells)`
- `e(cband_crit)` — sup-t critical value (when `cband`)
- `e(att_gt)` — matrix `g, t, ATT, SE, n`
- `e(es)` — matrix `e, ES(e), SE, .`
- `e(band)` — matrix `g, t, ATT, se_boot, uniform_lo, uniform_hi`

## Package structure

| File | Role |
| --- | --- |
| `ddd.ado` | main command: parse → validate → prep → estimate → aggregate → bootstrap → display → `ereturn` |
| `_ddd_cell.ado` | one DDD cell (RA / IPW / DR), Eq. (3.5) with **three** comparisons |
| `_ddd_estall.ado` | loops over cohorts/times, returns the `ATT(g,t)` matrix |
| `_ddd_display.ado` | formatted output (printing fully decoupled from estimation) |
| `ddd.sthlp` | help file |
| `examples/` | simulated data + end-to-end example |

## Status / roadmap

- **v0.3.0 (current)**: adds **two-step GMM combination across comparison
  cohorts** (`gmm`) and **cluster bootstrap** (`cluster()`), on top of v0.2's
  multinomial-logit propensity score, `control(notyet)`, and sup-t bands.
- **v0.4 (planned)**: CS-style aggregation weights; covariates in the `gmm`
  cross-covariance.

## Caveats

- `covariates()` must be **time-invariant** within unit.
- v0.1 aggregation weights cohorts by treated-cell size; this is simpler than
  the paper's optimal scheme.
- The analytic SEs hold nuisance functions fixed at their estimates, so they
  understate uncertainty in the covariate case; use `reps()` for the overall ATT.

## References

Ortiz-Villavicencio, M. and P. H. C. Sant'Anna (2025). *Better Understanding
Triple Differences Estimators.* arXiv:2505.09942.

## Author

Haoyu Niu
