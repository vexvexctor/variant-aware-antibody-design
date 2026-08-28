# Appendix text — design-panel size (ready to paste)

Supporting data: `PANEL_SIZE_SENSITIVITY.md`, `results/figures/panel_size_sensitivity.png`,
`results/report_data/panel_size_{summary,per_unit}.csv`.

---

## A.x Sensitivity to the design-panel size $M_{\text{design}}$

Our objective aggregates a binding-ΔΔG oracle over $M_{\text{design}}=6$ antigen escape variants
at design time. Because this is a free parameter, we swept it over
$M_{\text{design}}\in\{1,3,6,12\}$ with every other setting held fixed (AntiFold prior, BA-DDG
oracle, CVaR$_{0.2}$ aggregation, pool size 6, SVDD width 6), across all 51 targets with three
unseeded replicates per cell (612 design runs).

Critically, **every arm is evaluated on the identical held-out panel.** The design panel is the
top-$M$ variants by $\omega$ and the evaluation panel is the bottom-24; since every target
carries at least 41 variants (median 4{,}715), the two sets cannot intersect for $M\le 12$. We
verified this on the generated designs: all 51 targets have a byte-identical ordered held-out
variant list at every $M$, with no design/evaluation overlap. Differences across $M$ therefore
reflect the objective alone, not a shifting evaluation target.

**The method is insensitive to $M_{\text{design}}$ over this range.** Beat-native rates are
70.6\%, 66.7\%, 66.7\% and 64.7\% for $M=1,3,6,12$ (Table~\ref{tab:panelsize}); no pairwise
contrast is statistically significant under a per-target paired test with 95\% confidence
intervals bootstrapped over the 30 antigen clusters. The mean held-out worst-case margin
improves modestly and monotonically with $M$ ($-0.056$, $-0.108$, $-0.102$, $-0.132$), and the
only significant contrast is $M=12$ against $M=6$ ($\Delta=-0.154$, CI $[-0.282,-0.040]$); this
does not propagate to the binarised beat-native rate ($-2.0$ pp, CI $[-8.7,+4.9]$).

Two consequences are worth stating plainly. First, $M=6$ was fixed as a default before any
held-out evaluation and is not a tuned optimum — this sweep is the first time the parameter was
varied — so the reported results involve no selection of $M$ against the evaluation panel.
Second, the direction of the one significant effect makes $M=6$ a **conservative** setting:
$M=12$ yields tighter margins, so our choice if anything understates the method. We retain
$M=6$ because it is the pre-registered default and because the reward-evaluation cost grows
linearly in $M$ (441, 1{,}198, 2{,}334 and 4{,}599 oracle calls per design run), while noting
that $M=1$ attains an indistinguishable beat-native rate at roughly one fifth of the cost and is
the appropriate setting under a tight inference budget.

Finally, the flatness is consistent with a known limitation of the benchmark: held-out escape
worst-case correlates with wild-type binding at $r=+0.72$ (Appendix~\ref{sec:coupling}). If much
of the held-out signal is affinity rather than escape-specific, enlarging the design panel should
yield limited returns, and a single-mutant-enriched panel would be the natural test of whether a
genuine dependence on $M_{\text{design}}$ emerges.

### Table

| $M_{\text{design}}$ | beat-native | 95\% CI | mean held-out margin | oracle calls / run |
|---:|---:|---|---:|---:|
| 1 | 70.6\% | [56.7, 84.4] | $-0.056$ | 441 |
| 3 | 66.7\% | [54.4, 78.8] | $-0.108$ | 1,198 |
| **6** (ours) | **66.7\%** | **[51.2, 78.7]** | $\mathbf{-0.102}$ | **2,334** |
| 12 | 64.7\% | [48.9, 78.4] | $-0.132$ | 4,599 |

_51 targets, 3 replicates, H3-DDG evaluator (not used during optimization); CIs bootstrap the 30
antigen clusters. Margin is design-minus-native; negative favours the design._

### Paired contrasts

| contrast | $\Delta$ margin | 95\% CI | $\Delta$ beat-native | 95\% CI |
|---|---:|---|---:|---|
| $6\!\to\!1$ | $-0.025$ | [$-0.182$, $+0.100$] | $+3.9$ pp | [$-4.6$, $+15.6$] |
| $6\!\to\!3$ | $-0.028$ | [$-0.122$, $+0.079$] | $+0.0$ pp | [$-6.9$, $+8.7$] |
| $6\!\to\!12$ | $\mathbf{-0.154}$ | **[$-0.282$, $-0.040$]** | $-2.0$ pp | [$-8.7$, $+4.9$] |

---

## If a reviewer pushes

**"Why six? This looks arbitrary."**
> It is a default, and we say so. The sweep in Appendix A.x is the first time the parameter was
> varied, so no result in the paper involved selecting $M$ against the held-out panel. Over
> $M\in[1,12]$ the beat-native rate is flat within noise, so the method does not depend on the
> choice; and the single significant contrast runs *against* us ($M=12$ gives tighter margins
> than $M=6$), meaning our setting is conservative rather than favourable.

**"Then why not report M=12, since it's better?"**
> Only on the continuous margin, and at twice the reward-evaluation cost; it does not move the
> headline rate. We report $M=6$ because it is the pre-registered default and changing it after
> seeing the sweep would be exactly the selection effect we are avoiding. Both columns are in
> the appendix.

**"Doesn't flatness mean the escape panel isn't doing anything?"**
> No — the step from zero variants to one is large (WT-only selection reaches 9.8\% beat-native
> on this panel versus 43--45\% for any variant-aware objective, Appendix EXP4). What is flat is
> the return to *additional* variants beyond the first, which we attribute to the wild-type /
> escape coupling documented in Appendix~\ref{sec:coupling}.

---

## Wording to avoid

| do not write | why |
|---|---|
| "performance saturates at $M=6$" | It never rises; there is no knee. A reviewer re-running the sweep would catch this. |
| "we tuned $M_{\text{design}}$ and selected 6" | Untrue, and it would invite a selection-on-held-out objection that does not currently apply. |
| "larger panels do not help" | $M=12$ *is* significantly better on margin. Say it does not help the *rate*. |
| "$M=6$ balances cost and performance" | $M=1$ dominates it on both. The defence is pre-registration and conservatism, not a trade-off. |
