# Inference trajectory schema

One JSON per (target, generator, mode). Filename: `{TARGET}__{BASE}__{guided|unguided}.json`.
`guided` = svdd_k 6 (6 candidates scored and resampled per step);
`unguided` = svdd_k 1 (one candidate, no selection) — the same code path, so the two are
directly comparable step-for-step.

```
{
  target, base, oracle, svdd_k, svdd_temp, steer_agg, cvar_alpha, native_cdr, n_cdr,
  guided: bool,
  trajectories: [                 # one per sampled design
    { design_index, final_cdr,
      steps: [                    # one per denoising step; ONE position is revealed per step
        { step,                   # timestep counter (counts down)
          n_masked,               # positions still unrevealed
          n_candidates,           # candidates drawn (= svdd_k)
          n_unique_candidates,    # after dedup; the rest are cache hits
          predicted_clean,        # argmax-filled sequence the reward actually scored
          prior_entropy,          # mean entropy of the generator's distribution over masked pos
          reward_min/mean/max,    # reward across candidates this step
          reward_spread,          # max - min: how much discrimination the reward had
          reward_selected,        # reward of the candidate that survived
          candidates: [
            { pos, aa,            # the proposed reveal
              prior_p,            # generator's probability for that residue
              reward,             # objective's score
              resample_w,         # softmax(reward/temp) weight
              selected }          # did it survive
          ] } ] } ]
}
```

## What this supports
- Rejection dynamics: `selected` vs `resample_w` per step shows what the objective discards.
- Prior-vs-reward disagreement: `prior_p` of the selected candidate vs the max `prior_p`
  quantifies how far the objective pulls away from the generator's preference.
- Where steering acts: `reward_spread` by `step` / `n_masked` shows whether discrimination
  happens early (many positions open) or late (loop nearly committed).
- Guided vs unguided: identical fields, so any per-step quantity can be differenced directly.
