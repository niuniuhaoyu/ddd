*! _ddd_display.ado
*! Formatted output for ddd. Keeps all printing separate from estimation.
program define _ddd_display
    version 16
    syntax, ATT(name) ES(name) OVERall(real) ///
        [SEOver(real 0) CIlo(real 0) CIhi(real 0) METhod(string) DEPVar(string) ///
         COHorts(numlist) NOGraph AGGregation(string)]

    di ""
    di as txt "{hline 78}"
    di as txt "Triple Differences (DDD) estimation"
    di as txt "  Ortiz-Villavicencio & Sant'Anna (2025), arXiv:2505.09942"
    di as txt "{hline 78}"
    di as txt "Outcome:        " as res "`depvar'"
    di as txt "Estimator:      " as res "`method'"
    di as txt "Cohorts (g):    " as res "`cohorts'"
    di as txt "Aggregation:    " as res "`aggregation'"
    di ""

    di as txt "ATT(g,t):"
    di as txt "        g        t         ATT          SE"
    local nr = rowsof(`att')
    forvalues r = 1/`nr' {
        di as txt %9.0f `att'[`r',1] "  " %8.0f `att'[`r',2] "  " ///
            %11.4f `att'[`r',3] "  " %11.4f `att'[`r',4]
    }

    if rowsof(`es') > 0 {
        di ""
        di as txt "Event study ES(e):"
        di as txt "        e         ATT          SE        CI low       CI high"
        local ne = rowsof(`es')
        forvalues r = 1/`ne' {
            local elo = `es'[`r',2] - 1.96*`es'[`r',3]
            local ehi = `es'[`r',2] + 1.96*`es'[`r',3]
            di as txt %9.0f `es'[`r',1] "  " %11.4f `es'[`r',2] "  " ///
                %11.4f `es'[`r',3] "  " %11.4f `elo' "  " %11.4f `ehi'
        }
    }

    di ""
    di as txt "{hline 78}"
    di as txt "Overall ATT:    " as res %11.4f `overall'
    if `seover' < . {
        di as txt "  Bootstrap SE: " as res %11.4f `seover'
    }
    if `cilo' < . {
        di as txt "  95% CI:       [" as res %9.4f `cilo' as txt ", " ///
            as res %9.4f `cihi' as txt "]"
    }
    di as txt "{hline 78}"

    if "`nograph'" == "" {
        if rowsof(`es') > 0 {
            preserve
            clear
            local ne = rowsof(`es')
            quietly set obs `ne'
            quietly gen double __e = .
            quietly gen double __att = .
            quietly gen double __lo = .
            quietly gen double __hi = .
            forvalues r = 1/`ne' {
                quietly replace __e   = `es'[`r',1] in `r'
                quietly replace __att = `es'[`r',2] in `r'
                quietly replace __lo  = `es'[`r',2] - 1.96*`es'[`r',3] in `r'
                quietly replace __hi  = `es'[`r',2] + 1.96*`es'[`r',3] in `r'
            }
            quietly twoway ///
                (rarea __lo __hi __e, color(gs13%50) lwidth(none)) ///
                (connected __att __e, lcolor(navy) mcolor(navy) msymbol(O)), ///
                title("DDD event study") ///
                xtitle("Elapsed time (e = t - g)") ytitle("ATT(g,g+e)") ///
                yline(0, lpattern(dash) lcolor(gs8)) legend(off) name(ddd_es, replace)
            restore
        }
    }
end
