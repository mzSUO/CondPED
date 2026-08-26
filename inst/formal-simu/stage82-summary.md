# Stage 8.2 Summary — Simulation I formal 500 reps

- scenarios: I-0 null, I-1 single, I-2 two_linked_trait_specific (r2 = 0.3, spve = 0.020)
- R = 500 per scenario; 1500/1500 replicates ok; retries: 0 (recovered 0)
- wall time 14024 s; peak RSS 104 MB; workers = 4

## Primary endpoints (freeze 1.53)

| scenario | full recovery | secondary recovery | attribution exact | direction | cond effect RMSE | marker FDP |
|---|---|---|---|---|---|---|
| I-0 null | — | — | — | — | — | — (see Type I below) |
| I-1 single | 0.858 | — | 0.951 | 0.974 | 0.031 | 0.071 |
| I-2 two-linked | 0.778 | 0.778 | 0.814 | 0.820 | 0.061 (beta RMSE; see note) | 0.051 |

## Secondary diagnostics

| scenario | signal_count_exact | extra-signal rate | novel-locus extra share |
|---|---|---|---|
| I-1 | 0.949 | 0.106 | 0.697 |
| I-2 | 0.649 | 0.486 | (see extra_provenance.csv) |

## Oracle-causal-set resolver diagnostic (I-2, reps 1-50): recovery = 0.98 (49/50)

## Extra-signal provenance (I-1 + I-2)

- novel_locus_fp: 350; split_region_fp: 109; within_locus_fp: 43

## Type I (I-0 null)

- omnibus Type I (marker-level BH-selected null share / any-FP rate): 0.049

Note: the frozen evaluator defines conditional_effect_* only for
single-signal modes; for I-2 the table reports beta RMSE (effect
estimation error on matched signals) instead. I-1 marker FDP (0.071) is
slightly above the nominal 0.05 — BH controls E[FDP] under
independence/PRDS; mildly correlated test statistics can exceed it.
