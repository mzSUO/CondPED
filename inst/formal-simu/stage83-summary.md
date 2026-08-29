# Stage 8.3 Summary — Simulation II formal 500 reps

- II-A: A1 trait_specific / A2 two_trait_concordant / A3 two_trait_antagonistic / A4 broad_concordant (signal_oracle)
- II-B: R1 highly_representable (rho = 0.05) / R3 strongly_nonredundant (signal_trait_oracle)
- R = 500 per scenario; 3000/3000 replicates ok; retries: 398
- wall time 6024 s; peak RSS 107 MB; workers = 4

## II-A: CondPED vs ASSET trait-set comparison (freeze 1.31)

CondPED (frozen evaluator):

| scenario | exact | precision | recall | Jaccard | direction |
|---|---|---|---|---|---|
| trait_specific | 0.940 | 0.971 | 0.996 | 0.967 | 0.940 |
| two_trait_concordant | 0.926 | 0.979 | 0.989 | 0.969 | 0.950 |
| two_trait_antagonistic | 0.926 | 0.979 | 0.993 | 0.972 | 0.926 |
| broad_concordant | 0.784 | 1.000 | 0.938 | 0.938 | 1.000 |

ASSET (subset reference; no rho/Rep per freeze 1.35):

| scenario | exact | precision | recall | Jaccard | direction ok |
|---|---|---|---|---|---|
| trait_specific | 0.552 | 0.774 | 1.000 | 0.774 | 0.552 |
| two_trait_concordant | 0.136 | 0.960 | 0.571 | 0.556 | 0.918 |
| two_trait_antagonistic | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| broad_concordant | 0.000 | 1.000 | 0.406 | 0.406 | 1.000 |

## II-B: representation structure (freeze 1.34)

| scenario | eta bias | eta RMSE | rho MAE | rho map MAE | beta cov95 | thr-side acc | min Rep card | Rep exact | Irr recovery |
|---|---|---|---|---|---|---|---|---|---|
| R1 | 0.0003 | 0.0403 | 0.0347 | 0.0347 | 0.954 | 0.955 | 0.722 | 0.722 | 0.722 |
| R3 | 0.0003 | 0.0403 | 0.0552 | 0.0552 | 0.945 | 1.000 | 0.994 | 0.994 | 0.994 |

ASSET breadth reference: R1 = 1.04, R3 = 3.13 (CondPED min Rep: 0.72 vs 0.99; truth: 1 vs 4).

Note: the frozen evaluator defines beta_coverage (effect estimates);
eta coverage95 is not an evaluator output and is therefore not reported.

## Fairness / truth-map checks (freeze 1.16 / 1.18)

- R1: full trait set 500/500, concordant direction 500/500, truth min Rep in [1, 1], mean realized PVE 0.0288
- R3: full trait set 500/500, concordant direction 0/500, truth min Rep in [4, 4], mean realized PVE 0.0081
- The simulator's acceptance loop regenerates any effect vector not
  satisfying the architecture conditions; violations above must be 0.

## Effect-scale matching check (freeze 1.17)

- mean truth Q-form v*beta' Sigma_P_ref^{-1} beta: R1 = 0.0800, R3 = 0.0800 (matched by construction);
- mean realized marginal signal PVE: R1 = 0.0288, R3 = 0.0081 (architecture-driven difference, expected).

## Sign-pattern distribution (freeze v2 honest reporting)

- R1: ++++ n=500 (100.0%)
- R3: +--+ n=269 (53.8%); +-+- n=231 (46.2%)
