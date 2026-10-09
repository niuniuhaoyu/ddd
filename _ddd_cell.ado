*! _ddd_cell.ado
*! One DDD cell estimator for a given treated cohort g and comparison cohort cmp,
*! implementing the regression-adjustment (RA), inverse-probability-weighting (IPW),
*! and doubly-robust (DR) estimators of Ortiz-Villavicencio & Sant'Anna (2025).
*!
*! The DR estimator is Eq. (3.5) of the paper: a function of THREE DiD comparisons
*!   (g,0), (cmp,1), (cmp,0)  -- treated cell is (g,1).
*!
*! Returns:
*!   r(att)          DDD point estimate
*!   r(se_analytic)  influence-function SE (nuisances held fixed at estimates)
*!   r(n)            number of observations used
*!   r(ok)           1 if all four cells non-empty, 0 otherwise
program define _ddd_cell, rclass
    version 16
    syntax varname, G(integer) CMP(integer) SVAR(varname) QVAR(varname) ///
        [COVariates(varlist) METhod(string) PSMethod(string) TRIM(real 0.001) ///
         SAVif(name)]

    local dy "`varlist'"
    local method = lower("`method'")
    if "`method'" == "" local method "dr"
    local psmethod = lower("`psmethod'")
    if "`psmethod'" == "" local psmethod "mnl"
    if !inlist("`method'","ra","ipw","dr") {
        di as error "_ddd_cell: method() must be ra, ipw, or dr"
        exit 198
    }

    tempvar c1 c2 c3 c4
    quietly gen byte `c1' = (`svar'==`g'   & `qvar'==1 & !missing(`dy'))
    quietly gen byte `c2' = (`svar'==`g'   & `qvar'==0 & !missing(`dy'))
    quietly gen byte `c3' = (`svar'==`cmp' & `qvar'==1 & !missing(`dy'))
    quietly gen byte `c4' = (`svar'==`cmp' & `qvar'==0 & !missing(`dy'))

    foreach c in 1 2 3 4 {
        quietly count if `c`c''
        if r(N) == 0 {
            return scalar ok = 0
            return scalar att = .
            return scalar se_analytic = .
            return scalar n = .
            exit
        }
    }
    return scalar ok = 1

    *------------------------------------------------------------------
    * 1. Outcome regression m^(c)(X) = E[DeltaY | cell c, X]  (c = 2,3,4)
    *    RA and DR use them; IPW sets m = 0.
    *------------------------------------------------------------------
    tempvar m2 m3 m4
    foreach c in 2 3 4 {
        quietly gen double `m`c'' = 0
    }
    if "`method'" != "ipw" {
        if "`covariates'" == "" {
            foreach c in 2 3 4 {
                quietly summarize `dy' if `c`c'', meanonly
                quietly replace `m`c'' = r(mean)
            }
        }
        else {
            tempvar cons
            quietly gen double `cons' = 1
            foreach c in 2 3 4 {
                quietly regress `dy' `covariates' if `c`c''
                matrix __ddd_b = e(b)
                matrix colnames __ddd_b = `covariates' `cons'
                matrix score `m`c'' = __ddd_b, replace
            }
        }
    }

    *------------------------------------------------------------------
    * 2. Generalized propensity score p^(c)(X) = P[cell=c | X]
    *    psmethod(mnl): multinomial logit (default, per the paper)
    *    psmethod(lpm): linear-probability working model (fallback)
    *    RA sets p to the unconditional cell shares (weights cancel out).
    *------------------------------------------------------------------
    tempvar p1 p2 p3 p4
    foreach c in 1 2 3 4 {
        quietly summarize `c`c'', meanonly
        quietly gen double `p`c'' = r(mean)
    }
    if "`method'" != "ra" {
        if "`covariates'" != "" {
            local okps = 0

            if "`psmethod'" == "mnl" {
                tempvar cellid
                quietly gen byte `cellid' = .
                quietly replace `cellid' = 1 if `c1'
                quietly replace `cellid' = 2 if `c2'
                quietly replace `cellid' = 3 if `c3'
                quietly replace `cellid' = 4 if `c4'
                capture quietly mlogit `cellid' `covariates' if !missing(`cellid'), base(1) iterate(30)
                if _rc == 0 {
                    capture quietly predict double `p1' `p2' `p3' `p4' if !missing(`cellid')
                    if _rc == 0 local okps = 1
                }
            }

            if `okps' == 0 {
                * fallback: linear-probability working model
                tempvar cons2
                quietly gen double `cons2' = 1
                foreach c in 1 2 3 4 {
                    quietly regress `c`c'' `covariates' if `c1'|`c2'|`c3'|`c4'
                    matrix __ddd_b = e(b)
                    matrix colnames __ddd_b = `covariates' `cons2'
                    matrix score `p`c'' = __ddd_b, replace
                }
            }

            foreach c in 1 2 3 4 {
                quietly replace `p`c'' = max(`p`c'', `trim')
                quietly replace `p`c'' = min(`p`c'', 1-`trim')
            }
        }
    }

    *------------------------------------------------------------------
    * 3. Weights: w_trt (treated) and w_comp^(c) (comparison cells)
    *------------------------------------------------------------------
    tempvar wtrt
    quietly summarize `c1', meanonly
    quietly gen double `wtrt' = `c1'/r(mean)

    tempvar w2 w3 w4
    foreach c in 2 3 4 {
        tempvar num
        quietly gen double `num' = `c`c'' * `p1' / `p`c''
        quietly summarize `num', meanonly
        quietly gen double `w`c'' = `num'/r(mean)
    }

    *------------------------------------------------------------------
    * 4. DR DDD (Eq. 3.5): three terms, signs + + -
    *------------------------------------------------------------------
    tempvar e1 e2 e3 ifv
    quietly gen double `e1' = (`wtrt' - `w2')*(`dy' - `m2')
    quietly gen double `e2' = (`wtrt' - `w3')*(`dy' - `m3')
    quietly gen double `e3' = (`wtrt' - `w4')*(`dy' - `m4')

    quietly summarize `e1', meanonly
    local a1 = r(mean)
    quietly summarize `e2', meanonly
    local a2 = r(mean)
    quietly summarize `e3', meanonly
    local a3 = r(mean)
    local att = `a1' + `a2' - `a3'

    quietly gen double `ifv' = `e1' + `e2' - `e3' - `att'
    quietly summarize `ifv'
    local sd = r(sd)
    quietly count if !missing(`ifv')
    local n = r(N)
    local se = cond(`n'>1, `sd'/sqrt(`n'), .)

    if "`savif'" != "" {
        quietly gen double `savif' = `ifv'
    }

    return scalar att = `att'
    return scalar se_analytic = `se'
    return scalar n = `n'
    return scalar n_trt = .
end
