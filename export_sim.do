clear all
set more off
adopath + "C:/Users/35223/Desktop/默认项目/程序包/ddd"
do "C:/Users/35223/Desktop/默认项目/程序包/ddd/examples/ddd_simdata.do"
keep id time y enable eligible
export delimited using "C:/Users/35223/Desktop/默认项目/程序包/ddd/ddd_sim.csv", replace
display "exported ddd_sim.csv"
