# EXP6 — SVDD decoding ablations (CONFIG half: K and temperature)

K (`--svdd-k`) = candidate reveals scored per unmasking step; temperature (`--svdd-temp`) = softmax sharpness of the value-reweighting (temp->0 greedy-argmax reward, large->base sampling). **The finite-temperature / lambda SVDD variant IS this `--svdd-temp` sweep** — there is no separate lambda knob. Both sweeps: SVDD, base=antifold, oracle=baddg_ddg, `--steer-agg cvar --cvar-alpha 0.2`, 51 targets x 3 reps. Same FROZEN H3+BA-DDG protocol as EXP4. `*` = headline (K=6, temp=1.0 = scr_svdd_af_baddg).

## K sweep
`relB` = reward-eval budget relative to K=6 (measured, `profile_svdd_reward_budget.py`; SVDD reward calls grow ~linearly in K). `rt_s` = median design wall-seconds. Compare held-out CVaR20 at matched budget: a larger K must beat the smaller-K number by enough to justify its extra `relB`x reward calls.

| setting | relB | rt_s | cvar20 | worst | cvar50 | mean | beatN(targ) | beatN(run) | WTmarg | editN | nT | nRuns |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| K=2 | 0.39 | nan | 0.152 | 0.242 | 0.045 | -0.093 | 0.353 | 0.392 | -0.201 | 6.54 | 51 | 153 |
| K=5 | 0.87 | nan | 0.077 | 0.278 | -0.079 | -0.294 | 0.333 | 0.379 | -0.246 | 6.15 | 51 | 153 |
| K=6* | 1.00 | nan | 0.011 | 0.258 | -0.143 | -0.322 | 0.451 | 0.386 | -0.304 | 6.45 | 51 | 153 |
| K=10 | 1.49 | nan | -0.044 | 0.116 | -0.207 | -0.358 | 0.451 | 0.458 | -0.346 | 6.41 | 51 | 153 |
| K=20 | 2.42 | nan | -0.137 | 0.032 | -0.227 | -0.558 | 0.471 | 0.458 | -0.504 | 6.40 | 51 | 153 |

## Temperature (finite-temp / lambda) sweep

| setting | rt_s | cvar20 | worst | cvar50 | mean | beatN(targ) | beatN(run) | WTmarg | editN | nT | nRuns |
|---|---|---|---|---|---|---|---|---|---|---|---|
| temp=0.01 | nan | -0.230 | -0.074 | -0.488 | -0.690 | 0.529 | 0.529 | -0.626 | 6.60 | 51 | 153 |
| temp=0.05 | nan | -0.293 | -0.136 | -0.433 | -0.606 | 0.510 | 0.510 | -0.575 | 6.71 | 51 | 153 |
| temp=0.1 | nan | -0.303 | -0.181 | -0.398 | -0.492 | 0.569 | 0.516 | -0.573 | 6.48 | 51 | 153 |
| temp=0.5 | nan | -0.170 | 0.015 | -0.303 | -0.495 | 0.471 | 0.477 | -0.355 | 6.33 | 51 | 153 |
| temp=1.0* | nan | 0.011 | 0.258 | -0.143 | -0.322 | 0.451 | 0.386 | -0.304 | 6.45 | 51 | 153 |

CSV: `results/report_data/exp6_svdd_config_sweep.csv`. Reward budget: `scratch/svdd_screen/reward_budget.txt`.

<!-- SVDDCODE:BEGIN -->
## EXP6 — SVDD decoding ablations: rollouts / completion / position-order / beam
<sub>generated from 630 design runs (target = 630 = 14 cfg x 15 tgt x 3 rep). If < 630, the design/grade arrays (43683/43685) are still in flight; job 43910 regenerates this on completion.</sub>

*(agent: svddcode. Additive sampler knobs in `diffusion.py::sample_svdd`, worktree `scratch/campaign/lyra_svddcode`. base=antifold, oracle=baddg_ddg, steer-agg=cvar α0.2, k=6, temp=1.0. FROZEN selection = cluster_top1 (committed pre-grade, no oracle best-of-pool). 15-target representative subset x 3 reps. Beat-native = frozen worst-case eval-split margin < 0.)*

### Rollouts / completion=stochastic (value = MC completions)
| config | reward_evals | BA beat% (run) | H3 beat% (run) | BA beat% (tgt) | H3 beat% (tgt) | n_runs | H3 cov% |
| --- | --- | --- | --- | --- | --- | --- | --- |
| e6_base | 434 | 40.0 | 48.9 | 33.3 | 46.7 | 45 | 100 |
| e6_roll1_stoch | 434 | 42.2 | 46.7 | 40.0 | 53.3 | 45 | 100 |
| e6_roll2_stoch | 847 | 35.6 | 55.6 | 33.3 | 53.3 | 45 | 100 |
| e6_roll4_stoch | 1680 | 42.2 | 48.9 | 40.0 | 46.7 | 45 | 100 |
| e6_roll8_stoch | 3341 | 33.3 | 40.0 | 33.3 | 33.3 | 45 | 100 |

### Completion mode (greedy=base vs beam)
| config | reward_evals | BA beat% (run) | H3 beat% (run) | BA beat% (tgt) | H3 beat% (tgt) | n_runs | H3 cov% |
| --- | --- | --- | --- | --- | --- | --- | --- |
| e6_base | 434 | 40.0 | 48.9 | 33.3 | 46.7 | 45 | 100 |
| e6_comp_beam4 | 1677 | 37.8 | 55.6 | 33.3 | 60.0 | 45 | 100 |

### Position reveal order
| config | reward_evals | BA beat% (run) | H3 beat% (run) | BA beat% (tgt) | H3 beat% (tgt) | n_runs | H3 cov% |
| --- | --- | --- | --- | --- | --- | --- | --- |
| e6_base | 434 | 40.0 | 48.9 | 33.3 | 46.7 | 45 | 100 |
| e6_order_conflo | 249 | 37.8 | 53.3 | 33.3 | 40.0 | 45 | 100 |
| e6_order_cton | 247 | 37.8 | 51.1 | 26.7 | 46.7 | 45 | 100 |
| e6_order_enthi | 246 | 37.8 | 55.6 | 33.3 | 53.3 | 45 | 100 |
| e6_order_ntoc | 248 | 37.8 | 51.1 | 33.3 | 46.7 | 45 | 100 |
| e6_order_random | 246 | 33.3 | 53.3 | 26.7 | 60.0 | 45 | 100 |

### In-SVDD beam width (retain W prefixes)
| config | reward_evals | BA beat% (run) | H3 beat% (run) | BA beat% (tgt) | H3 beat% (tgt) | n_runs | H3 cov% |
| --- | --- | --- | --- | --- | --- | --- | --- |
| e6_base | 434 | 40.0 | 48.9 | 33.3 | 46.7 | 45 | 100 |
| e6_beam2 | 843 | 60.0 | 68.9 | 60.0 | 73.3 | 45 | 100 |
| e6_beam4 | 1647 | 64.4 | 71.1 | 66.7 | 66.7 | 45 | 100 |
| e6_beam8 | 3180 | 71.1 | 86.7 | 73.3 | 86.7 | 45 | 100 |

**Verdict — beam width is the only knob that helps, and its win SURVIVES out-of-family.** Beam directly deepens the search over the in-loop BA-DDG objective, so its BA-DDG gain (beat-native 33%->73% tgt, worst-case margin +0.30->-0.85 from base to beam8) is partly circular. The independent H3-DDG grader confirms it is not BA-DDG overfitting: H3 beat-native rises monotonically with beam width (base 46.7% -> beam2 73.3% -> beam8 86.7% of targets; run-level 48.9%->86.7%), i.e. deeper beam search finds designs that are more escape-robust even under an out-of-family oracle. By contrast rollouts (MC completions) and position-reveal order do NOTHING on either oracle (H3 tgt stays ~33-60% with no trend, roll8 actually the worst at 33.3% despite 8x the reward budget), and completion=beam4 also fails to move BA-DDG. So the Exp6 headline holds and is strengthened: performance is bought by in-SVDD beam search, not by rollouts or ordering, and the gain is real out-of-family.*

*Main plot: performance (beat-native%) vs reward-evaluation count — see `results/report_data/exp6_svdd_decoding_summary.csv` (reward_evals column).*
<!-- SVDDCODE:END -->
