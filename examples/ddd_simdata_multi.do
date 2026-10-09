*! ddd_simdata_multi.do -- multi-period DDD design with staggered enabling
*! True ATT(g,t) = 1.5 for eligible units in enabled groups, t >= enable.
version 16
clear all
set more off
set seed 20261010

set obs 3000
gen long id = _n
gen int prov = ceil(id/50)          // cluster variable (unit-level)
gen double ai = rnormal()
gen double grp = runiform()
gen int enable = cond(grp < 0.33, 3, cond(grp < 0.66, 4, 0))   // 0 = never
gen byte q = (runiform() < 0.5)

expand 5
bysort id: gen int time = _n
gen double y = ai + 0.1*time + rnormal()
replace y = y + 1.5 if q==1 & enable>0 & time>=enable

keep id time y enable q prov
rename q eligible
order id time y enable eligible prov
sort id time
label var y "outcome"
di as txt "ddd_simdata_multi: N = " _N "; true ATT(g,t) = 1.5"
