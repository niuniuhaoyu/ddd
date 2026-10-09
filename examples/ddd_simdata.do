*! ddd_simdata.do -- simulate a two-period DDD design with a known ATT(2,2)
*! Usage:  do ddd_simdata.do
*! Creates in memory: id, time, y, enable, eligible
version 16
clear all
set more off
set seed 20261009

set obs 4000
gen long id = _n
gen double ai = rnormal()
gen byte grp = (runiform() < 0.5)      // 1 = group enables the policy at t = 2
gen byte q   = (runiform() < 0.5)      // eligible population partition
gen double x = rnormal()               // time-invariant covariate

gen double y1 = ai + rnormal()
gen double y2 = ai + 0.5 + rnormal()          // common time trend
replace  y2   = y2 + 2.0 if grp==1 & q==1     // TRUE ATT(2,2) = 2

reshape long y, i(id) j(t)
gen int enable   = cond(grp==1, 2, 0)
gen int eligible = q
gen int time     = t
order id time y enable eligible x
sort id time
label var y        "outcome"
label var enable   "first period the group enables the policy (0 = never)"
label var eligible "population partition eligible for treatment (0/1)"
label var x        "time-invariant covariate"

di as txt "ddd_simdata: simulated N = " _N " obs; true ATT(2,2) = 2"
