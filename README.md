# Reassessing the Effect of Bank Branching Deregulation on Income Inequality

![Language](https://img.shields.io/badge/language-Stata-blue)
![Methods](https://img.shields.io/badge/methods-Goodman--Bacon%20%7C%20CSA%20%7C%20BJS%20%7C%20Honest%20DiD-orange)
![Topic](https://img.shields.io/badge/topic-Applied%20Econometrics-green)
![Status](https://img.shields.io/badge/status-Complete-brightgreen)

A modern-econometrics reassessment of Beck, Levine, and Levkov (2010, *Journal of Finance*), one of the most-cited papers on bank deregulation and inequality. Re-runs the original staggered two-way fixed-effects design with four post-2020 estimators — Goodman-Bacon decomposition, Callaway–Sant'Anna, Borusyak–Jaravel–Spiess, and Rambachan–Roth Honest DiD — to test whether the original finding survives once heterogeneous treatment effects and parallel-trends violations are taken seriously.

---

## Research Question

Beck, Levine, and Levkov (2010) used a staggered difference-in-differences design to argue that the U.S. intrastate bank branching deregulations of 1976–1994 substantially reduced state-level income inequality. Since 2010, the econometrics literature has shown that staggered TWFE estimators can be severely biased when treatment effects are heterogeneous across cohorts or over time. This project asks:

> Does the Beck–Levine–Levkov finding survive when re-estimated with modern, heterogeneity-robust DiD estimators — and if not, can the bias be characterized?

## Headline Findings

- **Goodman-Bacon decomposition fails.** Most of the weight in the original TWFE estimate comes from "forbidden" comparisons: later-treated states being compared to already-treated states, or treated states being compared to the 12 states that deregulated before the 1976 sample window. The clean treated-vs.-not-yet-treated comparisons receive disproportionately small weight.
- **Heterogeneity-robust estimators reverse the sign.** Both Callaway–Sant'Anna and Borusyak–Jaravel–Spiess produce statistically insignificant or even *positive* (inequality-worsening) point estimates on the log Gini, in stark contrast to the original −0.022 (p < 0.01).
- **BJS suffers from long-run parallel-trends violations.** The BJS event study shows large, sometimes statistically significant pre-trends, suggesting its strong long-pre-period identifying assumption is not satisfied on this panel.
- **CSA sensitivity tells a more nuanced story.** Restricting the CSA estimation window to 1975–1993 (excluding the federally-mandated Riegle-Neal deregulations) returns a **negative and weakly significant** estimate (−0.015, p < 0.10). The negative effect from Beck et al. partially re-emerges when control groups are restricted to states comparable in time.
- **Honest DiD bounds.** Under Rambachan–Roth sensitivity bounds, the CSA confidence interval crosses zero even at modest values of `M̄ = 0.25`, meaning small parallel-trends violations are sufficient to overturn the estimate.

**Bottom line:** the original finding does not cleanly survive modernization, but the negative-effect hypothesis is not cleanly rejected either. The result is more attenuated and more sensitive than Beck et al. report. This is a more nuanced conclusion than Baker, Larcker, and Wang (2022), who flagged this paper as a clear case where CSA overturns the original estimate.

## Methodology

| Stage | Approach |
|---|---|
| Replication | Re-estimates Beck et al. (2010) Table II, Panel A using the original TWFE specification on log Gini, logistic Gini, log Theil, log 90/10, and log 75/25. |
| Diagnostic 1 | **Goodman-Bacon (2021) decomposition** of the TWFE estimator into its constituent 2×2 DiD comparisons, with deregulation cohorts binned into five-year windows. |
| Estimator 1 | **Callaway–Sant'Anna (2021)** group-time ATTs with not-yet-treated controls; aggregated to overall ATT and event-study coefficients across all five inequality measures. |
| Estimator 2 | **Borusyak–Jaravel–Spiess (2024)** imputation estimator with 10 pre-trend lags and 15 post-treatment horizons. |
| Sensitivity 1 | TWFE re-estimated dropping the 12 always-treated states; CSA re-estimated on restricted windows (1975–1993 and 1975–1990) to test sensitivity to control-group composition. |
| Sensitivity 2 | **Rambachan–Roth (2023) Honest DiD** with the relative-magnitudes (Δ^RM) restriction over `M̄ ∈ {0, 0.25, 0.5, …, 2.0}`. |
| Percentile-by-percentile | CSA re-estimation across 19 percentiles (5th to 95th, by 5s) of the state-level income distribution, replicating Beck et al. Figure 2 under the modern estimator. |

## Data

The analysis operates on `macro_workfile.dta`, the state-year aggregate panel Beck et al. constructed from CPS microdata and used for their Table II. **This file is not redistributed here** — it must be rebuilt from Beck et al.'s original replication package, which includes the CPS aggregation scripts. See "Running the Analysis" below for the step-by-step procedure.

The final panel covers 48 states plus D.C. over 1976–2007, with the following key variables:

- `gini`, `theil`, `p10`, `p25`, `p75`, `p90` — inequality measures and percentile incomes
- `branch_reform` — year of intrastate branching deregulation (the staggered treatment)
- `_intra` — post-deregulation indicator
- `gsp_pc_growth`, `prop_blacks`, `prop_dropouts`, `prop_female_headed`, `unemploymentrate` — state-year controls

Underlying microdata comes from the Current Population Survey (CPS), restricted to individuals aged 25–54 with non-negative personal income, aggregated to the state-year level by Beck et al.'s preprocessing scripts.

## Repository Contents

```
.
├── DuiDavidssonRDModernizationPaper.pdf   ← full paper with motivation, decomposition derivation, results
├── DuiDavidssonRDModernizationCode.do     ← complete Stata replication script
└── README.md
```

## Running the Analysis

### 1. Install Stata packages

The script depends on several user-written Stata packages:

```stata
ssc install bacondecomp        // Goodman-Bacon decomposition
ssc install csdid              // Callaway-Sant'Anna
ssc install did_imputation     // Borusyak-Jaravel-Spiess
ssc install event_plot         // event-study plotting
ssc install honestdid          // Rambachan-Roth Honest DiD
ssc install estout             // regression tables
ssc install reghdfe            // high-dimensional fixed effects (dependency)
```

### 2. Obtain the data

This project does **not** redistribute Beck, Levine, and Levkov's data. To reproduce, you'll need to build their `macro_workfile.dta` yourself from their original replication materials:

1. **Download the Beck et al. replication package** from the *Journal of Finance* replication archive for Beck, Levine, and Levkov (2010) or from Ross Levine's website. This package contains the original `.do` files, the deregulation-year crosswalk, and the README.
2. **Follow the README in the Beck et al. replication package.** Their preprocessing scripts download CPS microdata for individuals aged 25–54 with non-negative personal income over 1976–2007 from IPUMS-CPS or the Census Bureau, apply the sample restrictions described in the original paper (drop South Dakota and Delaware, etc.), and aggregate the microdata into state-year inequality measures (Gini, Theil, and percentile ratios) to produce `macro_workfile.dta`. Their scripts must be run in the order specified in their README — the CPS aggregation step has to complete before `macro_workfile.dta` exists.
3. **Place the resulting `macro_workfile.dta` in your working directory.** Once it's there, this project's `.do` file picks up from that point: no CPS handling, no microdata aggregation — the modernization script operates on the same state-year panel Beck et al. used for their Table II.

#### ⚠️ IPUMS variable renames since 2009

Beck et al.'s `micro_workfile.do` was written in 2009 and hard-codes CPS variable names from the IPUMS extract format used at that time. **Several of those variables have since been renamed in IPUMS-CPS extracts**, which means the original `infix` block and subsequent `label var`/`replace` statements will fail (or, worse, silently produce wrong results) on a modern extract. Before running Beck et al.'s preprocessing scripts, you'll need to update the variable names.

Based on a current IPUMS-CPS extract (ASEC, 1977–2007, with the variables Beck et al. specify), the following renames apply:

| Beck et al. (2009) name | Current IPUMS name | Notes |
|---|---|---|
| `hhwt` | `asecwth` | ASEC household weight; IPUMS now names supplement weights explicitly. |
| `perwt` | `asecwt` | ASEC person weight, same logic. |
| `uhrswork` | `uhrsworkly` | "Usual hours worked per week, *last year*" — IPUMS added the `ly` suffix to distinguish from contemporaneous hours. |
| `educrec` | `educ` | IPUMS consolidated the education recode variables. |
| `hhwt04` | *(no direct replacement)* | 2004 supplement weight used for the ASEC expansion. The modern `asecwth` already incorporates the relevant adjustments, but verify against the IPUMS [user note on ASEC weights](https://cps.ipums.org/cps/asec_weights.shtml) for your sample years. |
| `perwt04` | *(no direct replacement)* | Same as above for person weights; verify against `asecwt`. |

Variables unchanged: `year`, `serial`, `gq`, `statefip`, `hhincome`, `relate`, `age`, `sex`, `race`, `hispan`, `educ99`, `higrade`, `empstat`, `classwly`, `wkswork1`, `inctot`, `incwage`, `incbus`.

In practice, this means editing Beck et al.'s `micro_workfile.do` to: (1) update the variable list in the IPUMS extract request, (2) update the `infix` block to use the new names and column positions from your fixed-width extract, (3) `rename` the new names back to Beck et al.'s names *before* running their downstream scripts (cleanest option — leaves `macro_workfile.do` untouched), or alternatively edit the downstream scripts to use the new names throughout.

If you encounter additional renames not listed here, you can verify any variable's current name on the [IPUMS-CPS variable browser](https://cps.ipums.org/cps-action/variables/group).

### 3. Update the path and run

In `DuiDavidssonRDModernizationCode.do`, edit the `cd` command on line 14 to point to your working directory:

```stata
cd "/path/to/your/working/directory";
```

Then run the full script in Stata. It will produce:

- TWFE replication estimates (Section 1)
- Goodman-Bacon decomposition plot (Section 2)
- CSA overall ATTs and event studies, plus the percentile-by-percentile re-estimation (Section 3)
- BJS overall ATTs and event studies (Section 4)
- A combined three-panel event-study comparison (Section 5)
- Honest DiD sensitivity bounds (Section 6)
- TWFE robustness dropping always-treated states (Section 7)
- CSA window-restricted sensitivity estimates (Section 8)

All figures are exported as `.png` and `.gph` files; all regression tables are stored as Stata `estimates` and printed to the log via `estout`.

## Skills Demonstrated

- **Modern causal inference**: full pipeline of post-2020 DiD diagnostics — Goodman-Bacon decomposition, Callaway–Sant'Anna with not-yet-treated controls, Borusyak–Jaravel–Spiess imputation, and Rambachan–Roth Honest DiD sensitivity bounds.
- **Critical replication**: not just re-running the original specification, but systematically stress-testing it under multiple identification assumptions and reporting where the original finding holds and where it breaks.
- **Stata proficiency**: substantial use of `csdid`, `did_imputation`, `bacondecomp`, `honestdid`, and `event_plot`, with custom event-study plot composition (`graph combine`) and percentile loops for figure replication.
- **Identification reasoning**: explicit treatment of how different estimators imply different parallel-trends assumptions, and how the choice of control group composition (never-treated vs. not-yet-treated, with or without window restriction) affects what each estimator is identifying.
- **Honest reporting**: a mixed-evidence finding presented as mixed evidence, rather than forcing a clean narrative in either direction. The paper engages directly with Baker, Larcker, and Wang (2022)'s claim that this paper is a clear-cut CSA reversal and adds nuance to that conclusion.

## Context

This was an independent research project completed for ECON 562 (Research Design and Policy Evaluation in Economics) at the University of British Columbia. The replication strategy, diagnostic design, sensitivity tests, and interpretation are my own work; the underlying data and original TWFE specification are from Beck, Levine, and Levkov (2010).

## Citation

```bibtex
@unpublished{Davidsson2026BankDereg,
 author = {Davidsson, Dui},
 title = {Reassessing the Effect of Bank Branching Deregulation on Income Inequality},
 year = {2026},
 note = {Course project, PhD Econometrics, University of British Columbia},
}
```

## Author

**Dui Davidsson** — [GitHub](https://github.com/DavidssonDui)

---

*Full citations and methodological references are in the [paper](DuiDavidssonRDModernizationPaper.pdf).*
