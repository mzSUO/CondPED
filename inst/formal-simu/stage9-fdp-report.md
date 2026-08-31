# Stage 9 batch 11: FDR internal calibration (freeze §1.52)

- null / single / mixed, R=200 each, 600/600 ok

| scenario | n | marker FDP [Wilson 95%] | type1 | extra rate |
|---|---|---|---|---|
| mixed_multisignal | 200 | 0.044 [0.040, 0.049] | 0.093 | 0.412 |
| null | 200 | 0.867 [0.621, 0.963] | 0.049 | 0.065 |
| single_multi_trait | 200 | 0.079 [0.048, 0.126] | 0.050 | 0.052 |

- 正式 claim（§1.52）：marker FDP ≈ 0.05（含 CI）
- locus/signal FDP 不作为 primary error-control claim（marker BH 不自动产生 signal-level 控制）
