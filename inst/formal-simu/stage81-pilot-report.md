# Stage 8.1 Pilot Report

- grid: 12 scenarios x R = 20 reps (freeze doc 05 模拟20260825.md)
- launcher: bash run_sim.sh stage81 inst/formal-simu/run-stage81.R
- main grid wall time: 2117 s; total (incl. retry/resume/probe): 2153 s
- peak RSS (master + workers): 107 MB
- replicates: 240 ok / 240 total

## Pilot checks (freeze 1.62)

| # | check | result | evidence |
|---|-------|--------|----------|
| 1 | 1. realized LD (target 0.3) | PASS | median realized mean r2 = 0.300 over 40 linked reps (range 0.300-0.300) |
| 2 | 2. signal PVE (target 0.02, per-active-trait) | PASS | fixed-dirs scenarios: median realized per-active-trait PVE = 0.0200 (n=180 signals); scale-match scenarios by design: median realized total PVE = 0.0157 (Q-form matched) |
| 3 | 3. Sigma_G^bg PSD | PASS | min eigenvalue over 240 reps = 1.605e-01 |
| 4 | 4. truth rho map | PASS | 240 reps with representation_map, 0 with non-finite/negative rho |
| 5 | 5. Rep/Irr truth (R1 min Rep=1, R3 min Rep=4) | PASS | R1/R3 reps checked: 40, truth min-Rep size correct: 40 |
| 6 | 6. ASSET correlation matrix | PASS | paired reps: 60; ASSET ran: 53 (both adapters ok); no-locus (ASSET not run): 7; adapter failures: 0 |
| 7 | 7. runtime | PASS | main grid wall = 2117s, total incl. retry/resume/probe = 2153s (240 reps, 4 workers) |
| 8 | 8. signal matching | PASS | 240 evaluated reps; representative-to-causal r2 out-of-range violations: 0 |
| 9 | 9. missing-output handling | PASS | resume skipped 240 ok reps; deleted probe regenerated: TRUE; content identical: TRUE |
| 10 | 10. numerical failure rate == 0 | PASS | 240/240 replicates ok; retries needed: 0 (recovered: 0) |

## Per-scenario completion

- linked_pseudo_multitrait: 20/20 ok; single_highly_representable: 20/20 ok; single_nonredundant: 20/20 ok; null: 20/20 ok; single_multi_trait: 20/20 ok; two_linked_trait_specific: 20/20 ok; broad_concordant: 20/20 ok; highly_representable: 20/20 ok; strongly_nonredundant: 20/20 ok; trait_specific: 20/20 ok; two_trait_antagonistic: 20/20 ok; two_trait_concordant: 20/20 ok
