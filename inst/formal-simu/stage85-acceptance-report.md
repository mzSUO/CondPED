# Stage 8.5 Acceptance Report — Formal Simulation Freeze

Aggregation over Stage 8.2 (Simulation I), 8.3 (Simulation II) and 8.4
(Simulation III) formal 500-rep results. All proportion estimates carry
Wilson 95% CIs; continuous metrics carry Monte Carlo SE.

Freeze-doc section references (current v2 document): Freeze Criteria
§1.63; metric tiers §1.67/§1.68/§1.69; figures §1.71/§1.72/§1.73;
pairing §1.64; E3 evidence §1.27; sign reporting §1.16 (v2).

## Headline numbers

- Stage 8.2: 1500/1500 ok | Stage 8.3: 3000/3000 ok | Stage 8.4: 1500/1500 ok; 0 numerical failure overall

### Tier1

| stage | scenario | metric | estimate [Wilson 95% CI] | n |
|---|---|---|---|---|
| 8.2 | I-0 | type1_omnibus (BH 0.05 calibration) | 0.049 [0.033, 0.072] | 500 |
| 8.2 | I-1 | locus_detection | 0.858 [0.825, 0.886] | 500 |
| 8.2 | I-1 | signal_recovery | 0.636 [0.590, 0.680] | 429 |
| 8.2 | I-1 | attribution_exact | 0.951 [0.926, 0.968] | 429 |
| 8.2 | I-1 | direction_recovery (pattern) | 0.974 [0.955, 0.986] | 429 |
| 8.2 | I-1 | beta_rmse (effect error) | 0.042 [0.040, 0.044] | 429 |
| 8.2 | I-2 | locus_detection | 0.872 [0.840, 0.898] | 500 |
| 8.2 | I-2 | signal_recovery | 0.000 [0.000, 0.009] | 436 |
| 8.2 | I-2 | secondary_recovery (evaluator) | 0.830 [0.795, 0.860] | 500 |
| 8.2 | I-2 | full_recovery (strict matching) | 0.778 [0.740, 0.812] | 500 |
| 8.2 | I-2 | secondary_recovery (strict matching) | 0.778 [0.740, 0.812] | 500 |
| 8.2 | I-2 | attribution_exact | 0.814 [0.775, 0.848] | 436 |
| 8.2 | I-2 | direction_recovery (pattern) | 0.820 [0.781, 0.853] | 436 |
| 8.2 | I-2 | beta_rmse (effect error) | 0.061 [0.059, 0.064] | 436 |
| 8.3 | A1 | attribution_exact (CondPED) | 0.940 [0.916, 0.958] | 500 |
| 8.3 | A1 | direction_recovery (pattern) | 0.940 [0.916, 0.958] | 500 |
| 8.3 | A1 | rho_MAE | 0.009 [0.005, 0.013] | 498 |
| 8.3 | A1 | rep_cardinality_recovery | 0.018 [0.009, 0.034] | 500 |
| 8.3 | A1 | rep_exact_family_recovery | 0.018 [0.009, 0.034] | 500 |
| 8.3 | A1 | irr_recovery | 0.018 [0.009, 0.034] | 500 |
| 8.3 | A2 | attribution_exact (CondPED) | 0.926 [0.900, 0.946] | 500 |
| 8.3 | A2 | direction_recovery (pattern) | 0.950 [0.927, 0.966] | 500 |
| 8.3 | A2 | rho_MAE | 0.056 [0.052, 0.060] | 500 |
| 8.3 | A2 | rep_cardinality_recovery | 0.638 [0.595, 0.679] | 500 |
| 8.3 | A2 | rep_exact_family_recovery | 0.614 [0.571, 0.656] | 500 |
| 8.3 | A2 | irr_recovery | 0.614 [0.571, 0.656] | 500 |
| 8.3 | A3 | attribution_exact (CondPED) | 0.926 [0.900, 0.946] | 500 |
| 8.3 | A3 | direction_recovery (pattern) | 0.926 [0.900, 0.946] | 500 |
| 8.3 | A3 | rho_MAE | 0.033 [0.029, 0.037] | 500 |
| 8.3 | A3 | rep_cardinality_recovery | 0.986 [0.971, 0.993] | 500 |
| 8.3 | A3 | rep_exact_family_recovery | 0.986 [0.971, 0.993] | 500 |
| 8.3 | A3 | irr_recovery | 0.986 [0.971, 0.993] | 500 |
| 8.3 | A4 | attribution_exact (CondPED) | 0.784 [0.746, 0.818] | 500 |
| 8.3 | A4 | direction_recovery (pattern) | 1.000 [0.992, 1.000] | 500 |
| 8.3 | A4 | rho_MAE | 0.061 [0.059, 0.064] | 500 |
| 8.3 | A4 | rep_cardinality_recovery | 0.776 [0.737, 0.810] | 500 |
| 8.3 | A4 | rep_exact_family_recovery | 0.428 [0.385, 0.472] | 500 |
| 8.3 | A4 | irr_recovery | 0.082 [0.061, 0.109] | 500 |
| 8.3 | R1 | attribution_exact (CondPED) | 1.000 [0.992, 1.000] | 500 |
| 8.3 | R1 | direction_recovery (pattern) | 0.838 [0.803, 0.868] | 500 |
| 8.3 | R1 | rho_MAE | 0.035 [0.033, 0.036] | 500 |
| 8.3 | R1 | rep_cardinality_recovery | 0.722 [0.681, 0.759] | 500 |
| 8.3 | R1 | rep_exact_family_recovery | 0.722 [0.681, 0.759] | 500 |
| 8.3 | R1 | irr_recovery | 0.722 [0.681, 0.759] | 500 |
| 8.3 | R3 | attribution_exact (CondPED) | 1.000 [0.992, 1.000] | 500 |
| 8.3 | R3 | direction_recovery (pattern) | 1.000 [0.992, 1.000] | 500 |
| 8.3 | R3 | rho_MAE | 0.055 [0.053, 0.057] | 500 |
| 8.3 | R3 | rep_cardinality_recovery | 0.994 [0.983, 0.998] | 500 |
| 8.3 | R3 | rep_exact_family_recovery | 0.994 [0.983, 0.998] | 500 |
| 8.3 | R3 | irr_recovery | 0.994 [0.983, 0.998] | 500 |
| 8.3 | A1 | attribution_jaccard (CondPED) | 0.967 [0.955, 0.979] | 500 |
| 8.3 | A2 | attribution_jaccard (CondPED) | 0.969 [0.958, 0.979] | 500 |
| 8.3 | A3 | attribution_jaccard (CondPED) | 0.972 [0.963, 0.981] | 500 |
| 8.3 | A4 | attribution_jaccard (CondPED) | 0.938 [0.926, 0.949] | 500 |
| 8.4 | E1 | locus_detection | 0.994 [0.983, 0.998] | 500 |
| 8.4 | E1 | attribution_exact | 0.459 [0.415, 0.503] | 497 |
| 8.4 | E1 | direction_recovery (pattern) | 0.996 [0.985, 0.999] | 497 |
| 8.4 | E1 | rho_MAE | 0.036 [0.034, 0.037] | 497 |
| 8.4 | E1 | rep_exact_family_recovery | 0.855 [0.821, 0.883] | 497 |
| 8.4 | E1 | beta_rmse (effect error) | 0.045 [0.043, 0.047] | 497 |
| 8.4 | E2 | locus_detection | 1.000 [0.992, 1.000] | 500 |
| 8.4 | E2 | attribution_exact | 0.580 [0.536, 0.622] | 500 |
| 8.4 | E2 | direction_recovery (pattern) | 0.858 [0.825, 0.886] | 500 |
| 8.4 | E2 | rho_MAE | 0.123 [0.113, 0.134] | 500 |
| 8.4 | E2 | rep_exact_family_recovery | 0.578 [0.534, 0.621] | 500 |
| 8.4 | E2 | beta_rmse (effect error) | 0.044 [0.042, 0.046] | 500 |
| 8.4 | E3 | locus_detection | 0.874 [0.842, 0.900] | 500 |
| 8.4 | E3 | attribution_exact | 0.805 [0.766, 0.840] | 437 |
| 8.4 | E3 | direction_recovery (pattern) | 0.809 [0.769, 0.843] | 437 |
| 8.4 | E3 | rho_MAE | 0.014 [0.010, 0.018] | 437 |
| 8.4 | E3 | rep_exact_family_recovery | 0.101 [0.076, 0.132] | 437 |
| 8.4 | E3 | beta_rmse (effect error) | 0.062 [0.060, 0.065] | 437 |

### Tier2

| stage | scenario | metric | estimate [Wilson 95% CI] | n |
|---|---|---|---|---|
| 8.2 | I-0 | marker_fdp (pooled) | 0.906 [0.758, 0.968] |   32 |
| 8.2 | I-1 | marker_fdp (pooled) | 0.071 [0.051, 0.098] |  463 |
| 8.2 | I-2 | marker_fdp (pooled) | 0.051 [0.045, 0.056] | 6126 |

### Tier3

| stage | scenario | metric | estimate [Wilson 95% CI] | n |
|---|---|---|---|---|
| 8.2 | I-0 | signal_count_exact | 0.000 [0.000, 0.562] |   3 |
| 8.2 | I-0 | extra_signal_rate | 0.056 [0.039, 0.080] | 500 |
| 8.2 | I-1 | signal_count_exact | 0.949 [0.924, 0.966] | 429 |
| 8.2 | I-1 | extra_signal_rate | 0.055 [0.038, 0.078] | 500 |
| 8.2 | I-2 | signal_count_exact | 0.649 [0.603, 0.692] | 436 |
| 8.2 | I-2 | extra_signal_rate | 0.228 [0.194, 0.267] | 500 |
| 8.4 | E1 | lead_breadth_inflation_rate (descriptive) | 0.000 [0.000, 0.008] | 497 |
| 8.4 | E1 | pseudo_multitrait_rate (descriptive) | 0.000 [0.000, 0.008] | 497 |
| 8.4 | E2 | lead_breadth_inflation_rate (descriptive) | 0.000 [0.000, 0.008] | 500 |
| 8.4 | E2 | pseudo_multitrait_rate (descriptive) | 0.006 [0.002, 0.017] | 500 |
| 8.4 | E3 | lead_breadth_inflation_rate (descriptive) | 0.105 [0.080, 0.137] | 438 |
| 8.4 | E3 | pseudo_multitrait_rate (descriptive) | 0.079 [0.057, 0.108] | 438 |

## Freeze Criteria audit (§1.63)

| criterion | result | evidence |
|---|---|---|
| 1.63.1 Statistical validity (Sigma_G^bg PSD) | PASS | Stage 8.1 pilot: min eigen(Sigma_G^bg) = 0.1605 over 240 reps; all Stage 8.2-8.4 reps passed the simulator's PSD acceptance (3000+1500+1500 ok) |
| 1.63.2 Numerical validity (0 rank-deficient, 0 numerical failure) | PASS | Stage 8.2 1500/1500 + Stage 8.3 3000/3000 + Stage 8.4 1500/1500 ok; 0 rank_deficient; 0 numerical failure in all manifests |
| 1.63.3 Truth validity (R1/R3 complete truth map) | PASS | Stage 8.3 fairness checks: R1/R3 full trait set 500/500 each; truth min Rep in [1,1] (R1) and [4,4] (R3); repflip note: all R1 failures are threshold_cross, none structural |
| 1.63.4 LD validity (realized r2 in target regime) | PASS | realized causal r2 = 0.302 (E3, n=500) and 0.300 (Stage 8.1 pilot, 40 linked reps) vs target 0.3 |
| 1.63.5 Effect validity (realized PVE near target) | PASS | Stage 8.1 pilot: fixed-dirs scenarios median realized per-active-trait PVE = 0.0200 (target 0.02); scale-match scenarios Q-form matched (8.3: 0.0800 vs 0.0800; 8.4: 0.0800 vs 0.0800) |
| 1.63.6 Pipeline validity (strict pairing) | PASS | seed family master_seed=20260826 verified on probes (8.2: 25/25, 8.3: 25/25, 8.4: 25/25 identical to .seed_for_rep); paired three-pipeline structure verified (8.2: 25/25, 8.4: 25/25) |

## Consistency checks

- seed family (master_seed = 20260826, .seed_for_rep deterministic): stage82 25/25; stage83 25/25; stage84 25/25
- paired three-pipeline structure (Lead+ASSET / Resolved+ASSET / CondPED on the same replicate): stage82 25/25; stage84 25/25

## Figure data tables

- Figure 2 (§1.71): figure2_e3_contamination.csv — conditional effect
  contamination gradient (marginal 0.137 -> resolved 0.059 -> conditional
  0.046); trait-breadth descriptive only (§1.27).
- Figure 3 (§1.72): figure3_r1_vs_r3.csv — R1 vs R3 full metrics; sign
  distributions honestly reported (R1: ++++ 100%; R3: +--+ 53.8% /
  +-+- 46.2%) per freeze v2.
- Figure 4 (§1.73): figure4_pipelines.csv — E1/E2/E3 three-pipeline
  comparison (lead vs resolved breadth, ASSET jaccard/direction,
  CondPED attribution/Rep).
