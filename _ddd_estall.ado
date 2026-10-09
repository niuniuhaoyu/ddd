*! _ddd_estall.ado
*! Estimate all ATT(g,t) cells under the DDD design.
*!
*! Inputs (current dataset must be UNIT-LEVEL / wide: one row per unit,
*! with outcome periods stored as `<depvar>1 ... <depvar>T`, plus
*! `svar()` = first-enabled-period index (0 = never enabled) and
*! `qvar()` = eligible partition (0/1)).
*!
*! control(never):  compare against never-enabled units (index 0).
*! control(notyet): compare against not-yet-enabled cohorts (> t).
*! gmm:             combine all identified comparison cohorts by two-step
*!                  (inverse-variance) GMM built from influence functions.
*!
*! Returns:
*!   r(att_gt)   matrix with columns: g, t, att, se, n
*!   r(ncells)   number of identified cells
program define _ddd_estall, rclass
    version 16
    syntax, DEPVAR(name) SVAR(varname) QVAR(varname) NT(integer) ///
        COHORTS(numlist) [COVariates(varlist) METhod(string) PSMethod(string) ///
        CONTrol(string) CMP(integer 0) GMM]

    local covopt ""
    if "`covariates'" != "" local covopt "covariates(`covariates')"
    local psopt ""
    if "`psmethod'" != "" local psopt "psmethod(`psmethod')"
    if "`control'" == "" local control "never"

    tempname pf
    tempfile res
    postfile `pf' double g double t double att double se double nobs using `res', replace

    foreach g of local cohorts {
        quietly count if `svar'==`g' & `qvar'==1
        if r(N) > 0 {
            local gprev = `g' - 1
            forvalues t = `g'/`nt' {
                * -- candidate comparison cohorts --
                if "`control'" == "never" {
                    local cmplist "0"
                }
                else {
                    local cmplist ""
                    foreach c of local cohorts {
                        if `c' > `t' local cmplist "`cmplist' `c'"
                    }
                    local cmplist "`cmplist' 0"
                }

                tempvar dy
                quietly gen double `dy' = `depvar'`t' - `depvar'`gprev'

                * -- evaluate each comparison cohort, keep influence functions --
                local j = 0
                foreach cmp of local cmplist {
                    tempvar psi
                    quietly _ddd_cell `dy', g(`g') cmp(`cmp') svar(`svar') qvar(`qvar') ///
                        `covopt' `psopt' method(`method') savif(`psi')
                    if r(ok) == 1 {
                        local ++j
                        local tau_`j' = r(att)
                        local psi_`j' "`psi'"
                    }
                }

                if `j' > 0 {
                    quietly count
                    local nobs = r(N)
                    local att = .
                    local se  = .
                    if `j' == 1 | "`gmm'" == "" {
                        local att = `tau_1'
                        quietly summarize `psi_1'
                        local se = r(sd)/sqrt(r(N))
                    }
                    else {
                        * -- two-step (inverse-variance) GMM combination --
                        local psilist ""
                        forvalues k = 1/`j' {
                            local psilist "`psilist' `psi_`k''"
                        }
                        mkmat `psilist', matrix(PSI)
                        matrix OMEGA = (PSI' * PSI) / (`nobs'^2)
                        matrix ONES  = J(`j',1,1)
                        matrix OMINV = invsym(OMEGA)
                        matrix W    = OMINV * ONES
                        matrix WN   = ONES' * W
                        matrix W    = W / WN[1,1]
                        matrix TAU  = J(1,`j',.)
                        forvalues k = 1/`j' {
                            matrix TAU[1,`k'] = `tau_`k''
                        }
                        matrix ATTm = TAU * W
                        local att = ATTm[1,1]
                        matrix DEN  = ONES' * OMINV * ONES
                        local se = sqrt(1/DEN[1,1])
                    }
                    post `pf' (`g') (`t') (`att') (`se') (`nobs')
                }
                quietly drop `dy'
            }
        }
    }
    postclose `pf'

    preserve
    quietly use `res', clear
    quietly count
    local ncells = r(N)
    if `ncells' > 0 {
        mkmat g t att se nobs, matrix(__ddd_M)
    }
    else {
        matrix __ddd_M = J(1,5,.)
    }
    restore

    return matrix att_gt = __ddd_M
    return scalar ncells = `ncells'
end
