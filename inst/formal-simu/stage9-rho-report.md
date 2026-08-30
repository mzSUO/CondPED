# Stage 9 batch 4: rho sensitivity (freeze §1.59.1)

- R1 target_loss in {0.02, 0.08}, R=200/cell, 400/400 ok

| rho | rep_exact | rep_card | rho MAE | thr-side acc | cand exact |
|---|---|---|---|---|---|
| 0.02 | 0.930 | 0.930 | 0.0313 | 0.988 | 1.000 |
| 0.08 | 0.465 | 0.465 | 0.0384 | 0.916 | 1.000 |

- reference: main rho=0.05 (Stage 8.3, R=500): rep_exact 0.722, rho MAE 0.035
- 核对点：R1/R3 结论方向不随 rho 点翻转
