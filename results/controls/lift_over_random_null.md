## Lift over each arm's OWN edit-distance-matched random control

Beat-native alone confounds design quality with how far the design moved from the native CDR: more edits score worse under the grader regardless of who proposed them. Each arm is therefore compared to a random control that mutates the SAME number of positions per design. LIFT = arm rate minus its own null rate, on the targets where both are graded; it is the part of the score not explained by edit distance.

Restricted to targets where ALL FOUR cells are graded (base, ours, and both nulls), so the two lifts are computed on the SAME targets and are comparable. Per-arm lift on differing subsets is not.

| generator | arm | beat-native | its random null | LIFT pp | Δ LIFT (ours-base) | n |
|---|---|---|---|---|---|---|
| AntiFold | `base_antifold` | 50.1% | 29.8% | **+20.4** | | 190 |
| AntiFold | `ours_antifold` | 64.3% | 21.9% | **+42.4** | **+22.0** ** [+12, +33] | 190 |
| AbMPNN | `base_abmpnn_ar` | 77.1% | 34.0% | **+43.1** | | 190 |
| AbMPNN | `ours_abmpnn_live` | 83.4% | 20.4% | **+62.9** | **+19.8** ** [+9, +30] | 190 |
| ProteinMPNN | `base_proteinmpnn_ar` | 53.5% | 7.9% | **+45.6** | | 190 |
| ProteinMPNN | `ours_proteinmpnn_live` | 62.6% | 7.2% | **+55.4** | **+9.8** ** [+3, +18] | 190 |

## Old-51 subset vs the new targets (is the expansion changing the story?)
