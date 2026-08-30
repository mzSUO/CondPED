# Stage 9 batch 3: LD sensitivity (freeze §1.57) + E3 contamination gradient (§1.27)

- target_r2 in {0.1, 0.5, 0.7} (R=200/cell) + main 0.3 (Stage 8.4, R=500); 600/600 ok

| target r2 | empirical r2 | |contam| marginal | |contam| conditional | reduction |
|---|---|---|---|---|
| 0.1 | 0.102 | 0.0770 | 0.0415 | 0.0355 |
| 0.3 | 0.301 | 0.1369 | 0.0455 | 0.0913 |
| 0.5 | 0.494 | 0.1790 | 0.0521 | 0.1268 |
| 0.7 | 0.606 | 0.1820 | 0.0678 | 0.1141 |

- **monotone in r2: YES**
- 核对点（§1.27）：marginal contamination 随 r2 单调上升，conditional 显著更低
