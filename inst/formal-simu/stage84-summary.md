# Stage 8.4 Summary — Simulation III formal 500 reps

- E1 single_highly_representable (++++, rho=0.05) / E2 single_nonredundant / E3 linked_pseudo_multitrait (r2=0.3, spve=0.020)
- R = 500 per scenario; 1500/1500 replicates ok; retries: 0
- wall time 13695 s; peak RSS 108 MB; workers = 4

## E1 vs E2: representation architecture (three paired pipelines)

| scenario | locus det | overall recovery | cand exact | direction | rep card | Rep exact | Irr | rho MAE | beta cov95 |
|---|---|---|---|---|---|---|---|---|---|
| E1 | 0.994 | 0.340 | 0.459 | 0.996 | 0.855 | 0.855 | 0.855 | 0.0356 | 0.958 |
| E2 | 1.000 | 0.578 | 0.580 | 0.858 | 0.578 | 0.578 | 0.578 | 0.1234 | 0.962 |

- eta bias E1 0.0004 / E2 0.0003; eta RMSE E1 0.0367 / E2 0.0376; threshold-side acc 0.975 / 0.926
- Q-form strength: E1 0.0800 vs E2 0.0800 (freeze 1.17 matched); beta RMSE 0.0451 / 0.0440
- sign patterns: E1: ++++ n=500 (100.0%); E2: +--+ n=226 (45.2%); +-+- n=273 (54.6%); ++-- n=1 (0.2%)
- pipeline breadth (descriptive): lead infl E1 -2.93 / E2 -0.84; resolved infl E1 -2.93 / E2 -0.84

## E3: conditional effect contamination (freeze 1.27, MAIN evidence)

- empirical r2 between causals: mean 0.302
- |contamination| marginal: 0.1369
- |contamination| resolved-signal: 0.0594
- |contamination| oracle-conditional: 0.0455
- reduction marginal -> conditional: 0.0913 (67%)
- reduction marginal -> resolved: 0.0838 (61%)

Trait-breadth inflation is descriptive only (freeze 1.27 downgrade):
- lead breadth inflation 0.110; resolved breadth inflation 0.184; pseudo-multitrait rate 0.079

Note: causal1 acts on Trait1 only (true Trait2 effect = 0); contamination
= |estimated Trait2 effect of causal1| under each analysis mode.
