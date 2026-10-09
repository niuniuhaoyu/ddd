*! ddd: Double (Triple) Difference-in-Differences estimators
*! Version 0.1.0   2026-10-09
*!
*! Reference:
*!   Ortiz-Villavicencio, M. and P. H. C. Sant'Anna (2025). "Better Understanding
*!   Triple Differences Estimators." arXiv:2505.09942.
*!
*! Implements regression-adjustment (RA), inverse-probability-weighting (IPW),
*! and doubly-robust (DR) estimators of ATT(g,t) in Triple Differences designs.
*!
*! Syntax:
*!   ddd depvar [if] [in], unit(varname) enable(varname) eligible(varname)
*!       time(varname) [covariates(varlist) method(ra|ipw|dr) aggregation(...)
*!       reps(#) seed(#) level(#) cluster(varname) nograph]
*!
*! enable()  : first period the unit's group ENABLES the policy (0 = never)
*! eligible(): 1 if the unit's population partition qualifies for treatment
*! unit()    : panel unit identifier
program define ddd, eclass
    version 16
    syntax varlist(min=1 max=1 numeric) [if] [in], ///
        UNIT(varname numeric) ENAble(varname numeric) ELIGible(varname numeric) ///
        TIME(varname numeric) ///
        [COVariates(varlist) METhod(string) PSMethod(string) COntrol(string) ///
         AGGregation(string) REPS(integer 0) SEED(integer 12345) LEVel(real 95) ///
         CLuster(varname numeric) NOGraph CBand GMM]
    preserve
    capture noisily _ddd_body `0'
    local __ddd_rc = _rc
    capture restore
    if `__ddd_rc' exit `__ddd_rc'
end

program define _ddd_body, eclass
    version 16
    syntax varlist(min=1 max=1 numeric) [if] [in], ///
        UNIT(varname numeric) ENAble(varname numeric) ELIGible(varname numeric) ///
        TIME(varname numeric) ///
        [COVariates(varlist) METhod(string) PSMethod(string) COntrol(string) ///
         AGGregation(string) REPS(integer 0) SEED(integer 12345) LEVel(real 95) ///
         CLuster(varname numeric) NOGraph CBand GMM]

    local outcome `varlist'
    local method = lower("`method'")
    if "`method'" == "" local method "dr"
    if !inlist("`method'","ra","ipw","dr") {
        di as error "ddd: method() must be ra, ipw, or dr"
        exit 198
    }
    if "`aggregation'" == "" local aggregation "cohort"
    if !inlist("`aggregation'","none","cohort","overall","es") {
        di as error "ddd: aggregation() must be none, cohort, overall, or es"
        exit 198
    }
    if `level' <= 0 | `level' >= 100 {
        di as error "ddd: level() must be between 0 and 100"
        exit 198
    }
    local covopt ""
    if "`covariates'" != "" local covopt "covariates(`covariates')"
    if "`control'" == "" local control "never"
    if "`psmethod'" == "" local psmethod "mnl"
    local gmmopt ""
    if "`gmm'" != "" local gmmopt "gmm"

    marksample touse
    quietly count if `touse'
    if r(N) == 0 {
        di as error "ddd: no observations"
        exit 2000
    }

    *===================================================================
    * PREP
    *===================================================================
    quietly keep if `touse'

    tempvar tidx sidx qvar uid
    * -- contiguous internal time index --
    quietly levelsof `time', local(tvals) clean
    local T : word count `tvals'
    if `T' < 2 {
        di as error "ddd: at least two time periods are required"
        exit 198
    }
    quietly gen long `tidx' = .
    local k = 0
    foreach v of local tvals {
        local ++k
        quietly replace `tidx' = `k' if `time' == `v'
    }
    quietly count if missing(`tidx')
    if r(N) {
        di as error "ddd: time() has missing/unmapped values"
        exit 459
    }
    * -- enable -> internal index (0 = never enabled) --
    quietly gen int `sidx' = 0
    foreach v of local tvals {
        local kk : list posof "`v'" in tvals
        quietly replace `sidx' = `kk' if `enable' == `v'
    }
    quietly count if `enable' != 0 & `sidx' == 0
    if r(N) {
        di as error "ddd: enable() values must be either 0 (never) or a value of time()"
        exit 459
    }
    * -- eligibility --
    quietly gen byte `qvar' = (`eligible' != 0)
    * -- unit id --
    quietly egen `uid' = group(`unit')
    * -- covariates must be time-invariant --
    if "`covariates'" != "" {
        foreach x of local covariates {
            quietly bysort `uid': gen byte __ddd_chg = (`x' != `x'[1])
            quietly count if __ddd_chg
            if r(N) {
                di as error "ddd: covariate `x' varies within unit; DDD covariates must be time-invariant"
                exit 459
            }
            quietly drop __ddd_chg
        }
    }
    * -- need a never-enabled comparison group --
    quietly count if `sidx' == 0
    if r(N) == 0 {
        di as error "ddd: no never-enabled units; v0.1 requires a never-enabled comparison group"
        exit 459
    }
    * -- uniqueness of (unit, time) --
    capture isid `unit' `tidx'
    if _rc {
        di as error "ddd: duplicate (unit,time) observations"
        exit 459
    }

    * -- keep only what is needed, then reshape to unit level (wide) --
    keep `uid' `sidx' `qvar' `tidx' `outcome' `covariates' `cluster'
    local iis "`uid' `sidx' `qvar'"
    if "`covariates'" != "" local iis "`iis' `covariates'"
    if "`cluster'"   != "" local iis "`iis' `cluster'"
    quietly reshape wide `outcome', i(`iis') j(`tidx')

    * -- treated cohorts --
    quietly levelsof `sidx' if `sidx' >= 2 & `qvar' == 1, local(cohorts) clean
    if "`cohorts'" == "" {
        di as error "ddd: no treated (enabled & eligible) units found"
        exit 459
    }
    * -- cohort sizes (for aggregation weights) --
    foreach g of local cohorts {
        quietly count if `sidx' == `g' & `qvar' == 1
        local w_`g' = r(N)
    }

    *===================================================================
    * POINT ESTIMATE
    *===================================================================
    quietly _ddd_estall, depvar(`outcome') svar(`sidx') qvar(`qvar') nt(`T') ///
        cohorts(`cohorts') `covopt' method(`method') control(`control') ///
        psmethod(`psmethod') `gmmopt'
    matrix ATT = r(att_gt)
    local ncells = r(ncells)
    if `ncells' == 0 {
        di as error "ddd: no identified ATT(g,t) cells"
        exit 2000
    }

    * -- weighted aggregation helper: overall + event-study --
    local att_over = 0
    local var_over = 0
    local wsum = 0
    forvalues r = 1/`=rowsof(ATT)' {
        local g = ATT[`r',1]
        local wg = `w_`g''
        local att_over = `att_over' + `wg'*ATT[`r',3]
        local wsum = `wsum' + `wg'
    }
    local att_over = `att_over'/`wsum'
    * event study ES(e)
    local emax = 0
    local emin = 0
    forvalues r = 1/`=rowsof(ATT)' {
        local e = ATT[`r',2] - ATT[`r',1]
        if `e' > `emax' local emax = `e'
        if `e' < `emin' local emin = `e'
    }
    local ne = `emax' - `emin' + 1
    if `ne' < 1 local ne = 1
    matrix ES = J(`ne',4,.)
    local row = 0
    forvalues e = `emin'/`emax' {
        local num = 0
        local den = 0
        local vsum = 0
        forvalues r = 1/`=rowsof(ATT)' {
            if (ATT[`r',2]-ATT[`r',1]) == `e' {
                local g = ATT[`r',1]
                local wg = `w_`g''
                local num = `num' + `wg'*ATT[`r',3]
                local den = `den' + `wg'
                local vsum = `vsum' + (`wg'^2)*(ATT[`r',4]^2)
            }
        }
        if `den' > 0 {
            local ++row
            matrix ES[`row',1] = `e'
            matrix ES[`row',2] = `num'/`den'
            matrix ES[`row',3] = sqrt(`vsum')/`den'
            matrix ES[`row',4] = .
        }
    }
    if `row' == 0 {
        matrix ES = J(1,4,.)
    }
    else if `row' < `ne' {
        matrix ES = ES[1..`row',1..4]
    }

    *===================================================================
    * BOOTSTRAP (optional, unit-level): SE/CI + sup-t uniform band
    *===================================================================
    if "`cband'" != "" & `reps' == 0 {
        di as error "ddd: cband requires reps() > 0"
        exit 198
    }
    local ncells_p = rowsof(ATT)
    local se_over = .
    local ci_lo = .
    local ci_hi = .
    local nboot = 0
    local crit = .
    matrix BAND = J(`ncells_p',6,.)
    if `reps' > 0 {
        tempfile basedata
        quietly save `basedata'
        tempname BOV
        matrix `BOV' = J(`reps',1,.)
        matrix BOOTGT = J(`reps',`ncells_p',.)
        set seed `seed'
        forvalues b = 1/`reps' {
            quietly use `basedata', clear
            if "`cluster'" != "" {
                quietly bsample, cluster(`cluster') idcluster(__ddd_bc)
            }
            else {
                quietly bsample
            }
            * recompute cohort sizes on the resample
            foreach g of local cohorts {
                quietly count if `sidx' == `g' & `qvar' == 1
                local wb_`g' = r(N)
            }
            capture quietly _ddd_estall, depvar(`outcome') svar(`sidx') qvar(`qvar') ///
                nt(`T') cohorts(`cohorts') `covopt' method(`method') ///
                control(`control') psmethod(`psmethod') `gmmopt'
            if _rc == 0 {
                matrix __B = r(att_gt)
                if rowsof(__B) > 0 {
                    local num = 0
                    local den = 0
                    forvalues r = 1/`=rowsof(__B)' {
                        local gg = __B[`r',1]
                        local tt = __B[`r',2]
                        local wg = `wb_`gg''
                        local num = `num' + `wg'*__B[`r',3]
                        local den = `den' + `wg'
                        forvalues c = 1/`ncells_p' {
                            if (ATT[`c',1]==`gg') & (ATT[`c',2]==`tt') {
                                matrix BOOTGT[`b',`c'] = __B[`r',3]
                            }
                        }
                    }
                    if `den' > 0 {
                        matrix `BOV'[`b',1] = `num'/`den'
                        local ++nboot
                    }
                }
            }
        }
        if `nboot' > 1 {
            * -- overall SE + percentile CI --
            preserve
            clear
            svmat double `BOV', names(bov)
            quietly summarize bov1
            local se_over = r(sd)
            _pctile bov1, p(2.5 97.5)
            local ci_lo = r(r1)
            local ci_hi = r(r2)
            restore
            * -- per-cell bootstrap SE + sup-t critical value --
            matrix SDV = J(`ncells_p',1,.)
            preserve
            clear
            svmat double BOOTGT, names(bg)
            tempvar tmax
            quietly gen double `tmax' = 0
            forvalues c = 1/`ncells_p' {
                quietly summarize bg`c'
                matrix SDV[`c',1] = r(sd)
                local sdc = r(sd)
                local mpc = ATT[`c',3]
                if `sdc' < . & `sdc' > 0 {
                    quietly replace `tmax' = max(`tmax', abs((bg`c' - `mpc')/`sdc')) if !missing(bg`c')
                }
            }
            local alpha = (100-`level')/100
            _pctile `tmax', p(`=100*(1-`alpha')')
            local crit = r(r1)
            restore
            * -- assemble band matrix (g, t, att, se_boot, lo, hi) --
            forvalues c = 1/`ncells_p' {
                local a = ATT[`c',3]
                local s = SDV[`c',1]
                matrix BAND[`c',1] = ATT[`c',1]
                matrix BAND[`c',2] = ATT[`c',2]
                matrix BAND[`c',3] = `a'
                matrix BAND[`c',4] = `s'
                matrix BAND[`c',5] = `a' - `crit'*`s'
                matrix BAND[`c',6] = `a' + `crit'*`s'
            }
        }
    }

    *===================================================================
    * DISPLAY
    *===================================================================
    local dispopt "seover(`se_over') cilo(`ci_lo') cihi(`ci_hi')"
    local graphopt ""
    if "`nograph'" != "" local graphopt "nograph"
    _ddd_display, att(ATT) es(ES) overall(`att_over') `dispopt' ///
        method(`method') depvar(`outcome') cohorts(`cohorts') `graphopt' ///
        aggregation(`aggregation')

    *===================================================================
    * ereturn
    *===================================================================
    ereturn clear
    ereturn scalar att = `att_over'
    ereturn scalar se_att = `se_over'
    ereturn scalar ci_lower = `ci_lo'
    ereturn scalar ci_upper = `ci_hi'
    ereturn scalar ncells = `ncells'
    ereturn scalar nboot = `nboot'
    ereturn scalar cband_crit = `crit'
    ereturn local cmd "ddd"
    ereturn local depvar "`outcome'"
    ereturn local method "`method'"
    ereturn local aggregation "`aggregation'"
    ereturn matrix att_gt = ATT
    ereturn matrix es = ES
    ereturn matrix band = BAND
end
