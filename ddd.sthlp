{smcl}
{* *! version 0.1.0  2026-10-09}{...}
{vieweralsosee "help ddd" "ddd"}{...}
{viewerjumpto "Syntax" "ddd##syntax"}{...}
{viewerjumpto "Description" "ddd##description"}{...}
{title:Title}

{phang2}
{bf:ddd} -- Triple Differences (DDD) estimators
{hline}

{pstd}
Ortiz-Villavicencio & Sant'Anna (2025), arXiv:2505.09942

{title:Syntax}

{phang2}
{cmd:ddd} {it:depvar} {ifin} {c -(}, {cmd:unit(}{it:varname}{cmd:)} {cmd:enable(}{it:varname}{cmd:)} {cmd:eligible(}{it:varname}{cmd:)} {cmd:time(}{it:varname}{cmd:)} {break}
    {cmd:,} {cmd:covariates(}{it:varlist}{cmd:)} {cmd:method(}{cmd:ra}{cmd:|}{cmd:ipw}{cmd:|}{cmd:dr}{cmd:)} {cmd:aggregation(}{cmd:none}{cmd:|}{cmd:cohort}{cmd:|}{cmd:overall}{cmd:|}{cmd:es}{cmd:)} {break}
        {cmd:reps(}{it:#}{cmd:)} {cmd:seed(}{it:#}{cmd:)} {cmd:level(}{it:#}{cmd:)} {cmd:cluster(}{it:varname}{cmd:)} {cmd:nograph}

{title:Description}

{pstd}
{cmd:ddd} estimates group-time average treatment effects {it:ATT(g,t)} in a
triple-differences (DDD) design, using the regression-adjustment (RA),
inverse-probability-weighting (IPW), and doubly-robust (DR) estimators of
Ortiz-Villavicencio and Sant'Anna (2025).

{pstd}
A DDD design requires a unit to satisfy two criteria to be effectively treated:
(a) it belongs to a group that has enabled the policy by period {it:t}, and
(b) it belongs to a population partition that qualifies for the treatment. The
command therefore requires {cmd:enable()} (the first period each group enables
the policy; 0 = never) and {cmd:eligible()} (0/1 partition indicator).

{pstd}
The DR estimator is a function of {bf:three} DiD comparisons -- the treated
cell against (g,0), (cmp,1) and (cmp,0) -- not two. Common implementations
("three-way fixed effects", or "difference of two DiDs") are generally invalid
when covariates matter; see the paper.

{title:Options}

{phang}
{cmd:unit(}{it:varname}{cmd:)} specifies the panel unit identifier.

{phang}
{cmd:enable(}{it:varname}{cmd:)} specifies the first period the unit's group
enables the policy; 0 (or a value not in {cmd:time()}) means never enabled.

{phang}
{cmd:eligible(}{it:varname}{cmd:)} is the 0/1 population-partition indicator.

{phang}
{cmd:time(}{it:varname}{cmd:)} is the time period.

{phang}
{cmd:covariates(}{it:varlist}{cmd:)} adds time-invariant covariates X.

{phang}
{cmd:method(ra|ipw|dr)} selects the estimator; default is {cmd:dr}.

{phang}
{cmd:aggregation(none|cohort|overall|es)}; v0.1 reports ATT(g,t) and ES(e).

{phang}
{cmd:reps(}{it:#}{cmd:)} number of unit-level bootstrap replications for the
overall ATT; default 0 (no bootstrap).

{phang}
{cmd:cluster(}{it:varname}{cmd:)} reserved for cluster bootstrap (v0.2).

{title:Stored results}

{pstd}
{cmd:ddd} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(att)}}overall ATT{p_end}
{synopt:{cmd:e(se_att)}}bootstrap SE of the overall ATT{p_end}
{synopt:{cmd:e(ci_lower)}}{cmd:e(ci_upper)}percentile CI{p_end}
{synopt:{cmd:e(ncells)}}number of identified ATT(g,t) cells{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(att_gt)}}g, t, ATT, SE, n{p_end}
{synopt:{cmd:e(es)}}e, ES(e), SE, .{p_end}

{title:Example}

{phang2}
{cmd:ddd y, unit(id) enable(enable) eligible(eligible) time(t) method(dr) reps(50)}

{title:Author}

{pstd}
Haoyu Niu

{title:References}

{pstd}
Ortiz-Villavicencio, M. and P. H. C. Sant'Anna (2025). Better Understanding
Triple Differences Estimators. arXiv:2505.09942.

{pstd}
Callaway, B. and P. H. C. Sant'Anna (2021). Difference-in-Differences with
Multiple Time Periods. {it:Journal of Econometrics} 225(2): 200--230.
