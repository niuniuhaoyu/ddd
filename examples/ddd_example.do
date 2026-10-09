*! ddd_example.do -- minimal end-to-end example for the ddd package
*! Usage:  do ddd_example.do
version 16
clear all
set more off

* --- make the package findable (edit the path if you move the folder) ---
capture adopath + "."
capture adopath + "D:\CFPS\pkgs\ddd"

* --- simulate ---
do ddd_simdata

* --- DR (default), no covariates ---
ddd y, unit(id) enable(enable) eligible(eligible) time(time) method(dr)

* --- RA with unit-level bootstrap for the overall ATT ---
ddd y, unit(id) enable(enable) eligible(eligible) time(time) method(ra) reps(100) seed(1)

di as txt _n "Reference: true ATT(2,2) = 2."
di as txt "The DIFFERENCE OF TWO DiDs (three-way FE) is also unbiased here"
di as txt "because there are no covariates; the paper's point is that it is"
di as txt "BIASED once covariates matter (see method(dr) with covariates())."
