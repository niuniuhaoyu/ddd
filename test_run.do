clear all
set more off
version 16
local PKG "C:/Users/35223/Desktop/默认项目/程序包/ddd"
adopath + "`PKG'"

do "`PKG'/examples/ddd_simdata.do"

display as txt _n "=== 1. dr, no covariates ==="
capture noisily ddd y, unit(id) enable(enable) eligible(eligible) time(time) method(dr)
display as txt "rc=" _rc "  e(att)=" e(att)

display as txt _n "=== 2. dr + covariates, psmethod(mnl) ==="
capture noisily ddd y, unit(id) enable(enable) eligible(eligible) time(time) ///
    covariates(x) method(dr) psmethod(mnl)
display as txt "rc=" _rc "  e(att)=" e(att)

display as txt _n "=== 3. dr + covariates, psmethod(lpm) fallback ==="
capture noisily ddd y, unit(id) enable(enable) eligible(eligible) time(time) ///
    covariates(x) method(dr) psmethod(lpm)
display as txt "rc=" _rc "  e(att)=" e(att)

display as txt _n "=== 4. control(notyet) on 2-period (falls back to never) ==="
capture noisily ddd y, unit(id) enable(enable) eligible(eligible) time(time) ///
    method(dr) control(notyet)
display as txt "rc=" _rc "  e(att)=" e(att)

*---------------------------------------------------------------------------
display as txt _n "=== 5. MULTI-PERIOD, control(notyet) ==="
do "`PKG'/examples/ddd_simdata_multi.do"
capture noisily ddd y, unit(id) enable(enable) eligible(eligible) time(time) ///
    method(dr) control(notyet)
display as txt "rc=" _rc "  e(att)=" e(att)

display as txt _n "=== 6. MULTI-PERIOD, control(never) ==="
capture noisily ddd y, unit(id) enable(enable) eligible(eligible) time(time) ///
    method(dr) control(never)
display as txt "rc=" _rc "  e(att)=" e(att)

display as txt _n "=== 7. cband (sup-t uniform band), 2-period ==="
capture noisily ddd y, unit(id) enable(enable) eligible(eligible) time(time) ///
    method(dr) reps(60) seed(2) cband
display as txt "rc=" _rc "  crit=" e(cband_crit)
capture matrix list e(band)

display as txt _n "=== 8. cband, MULTI-PERIOD (mnl covariates not used) ==="
do "`PKG'/examples/ddd_simdata_multi.do"
capture noisily ddd y, unit(id) enable(enable) eligible(eligible) time(time) ///
    method(dr) control(notyet) reps(60) seed(3) cband
display as txt "rc=" _rc "  crit=" e(cband_crit)
capture matrix list e(band)

display as txt _n "=== 9. MULTI-PERIOD gmm combination ==="
do "`PKG'/examples/ddd_simdata_multi.do"
capture noisily ddd y, unit(id) enable(enable) eligible(eligible) time(time) ///
    method(dr) control(notyet) gmm
display as txt "rc=" _rc "  e(att)=" e(att)

display as txt _n "=== 10. cluster bootstrap ==="
capture noisily ddd y, unit(id) enable(enable) eligible(eligible) time(time) ///
    method(dr) control(notyet) reps(40) seed(4) cluster(prov)
display as txt "rc=" _rc "  e(att)=" e(att) "  se=" e(se_att)

display as txt _n "done. (true ATT = 2 for two-period, 1.5 for multi-period)"
