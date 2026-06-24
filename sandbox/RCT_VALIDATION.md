# RCT subgroup-discovery: validation summary

Evidence that the `EffectXshift(rct = TRUE)` workflow discovers effect-modified
subgroups in randomized trials with honest, post-selection inference. All
simulations use a randomized binary treatment (alpha = 0.5), ranger + GLM
nuisances, 2 folds (5 in the demo), and the held-out CV-TMLE region estimates.
Scripts live alongside this file in `sandbox/`.

## 1. Does it invent subgroups? (type-I / false discovery)  — `rct_null_typeI.R`

Null DGP: **homogeneous** treatment effect (tau = 1 everywhere), 5 covariates, no
effect modification. 200 replicates, n = 1000.

| null setting | proposes a subgroup | **false-positive rate** (V−V^c CI excludes 0) |
|---|---|---|
| simple baseline | 0% | **0.0%** |
| nonlinear baseline (learner tempted) | 71% | **6.5%** (target ~5%) |

Even when a flexible learner is tempted to split on noise 71% of the time, the
sample-split / CV-TMLE inference holds the false-discovery rate near nominal.
This is the key contrast with naive RCT subgroup analysis.

## 2. Does it find real subgroups, with valid CIs? (power / coverage)  — `rct_breadth.R`

Alternative DGP: treatment effect 2 in the modified subgroup, 0 elsewhere.
True V−V^c contrast = 2. 40 replicates per setting.

| setting | detection | contrast est (truth 2) | 95% CI coverage |
|---|---|---|---|
| binary modifier, n=500 | 100% | 1.92 | 0.88 |
| binary modifier, n=1000 | 100% | 1.98 | 0.95 |
| binary modifier, n=2000 | 100% | 1.99 | 1.00 |
| continuous (threshold) modifier, n=2000 | 100% | 1.98 | 0.98 |

100% modifier detection; bias → 0 with n; coverage at/above nominal for
n ≥ 1000 (slightly under at n=500, as expected). Earlier focused runs
(`rct_validation.R`, `rct_incps_validation.R`) give the per-region breakdown:
ate and incps both ~0.92/0.98/0.94 coverage for V / V^c / contrast.

## 3. Single-arm prognostic risk mode (`target = "risk"`)  — `rct_risk_validation.R`

Binary-outcome DGP, true high-risk region {W1 = 1} (risk 0.721) vs {W1 = 0}
(risk 0.279). 40 replicates, n = 1500.

- Discovers W1 as the risk driver: **100%**
- Held-out region-risk estimates: V = 0.727 (truth 0.721), V^c = 0.275 (truth 0.279)
- Held-out V-risk CI coverage: **90%**; V risk > V^c risk: 100%

The risk mode reports a held-out *prediction* (region rate), not a causal effect.

## 4. Why the RCT path is well-behaved (unlike the continuous-mixture path)

The doubly-robust remainder is `(g_n − g)(Q_n − Q)`. In an RCT the propensity
`g` is **known** (randomization), so the remainder vanishes regardless of the
outcome model — the CV-TMLE SE is valid at moderate n. (The continuous
mixed-exposure path must estimate a conditional density `g(A|W)`, which is why
its CIs are anti-conservative at small n; see the `mixture_*` scripts.)

## Reproduce

```r
source("sandbox/rct_demo.R")              # end-to-end demo on one trial
Rscript sandbox/rct_null_typeI.R          # type-I (add NULL_HARD=TRUE for the hard null)
Rscript sandbox/rct_breadth.R             # power / coverage across n + modifier type
Rscript sandbox/rct_risk_validation.R     # prognostic risk mode
```
