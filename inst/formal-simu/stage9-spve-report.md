# Stage 9 batch 1: secondary-signal PVE sensitivity (freeze §1.56)

- spve grid {0.0075, 0.015, 0.0175}, R=200/cell, paired, 600/600 ok

| spve | locus det | secondary power | overall recovery | cand exact | direction | count_exact | extra rate |
|---|---|---|---|---|---|---|---|
| 0.0075 | 0.805 | 0.480 | 0.000 | 0.699 | 0.727 | 0.391 | 0.153 |
| 0.0150 | 0.790 | 0.675 | 0.000 | 0.741 | 0.747 | 0.639 | 0.205 |
| 0.0175 | 0.880 | 0.820 | 0.000 | 0.753 | 0.764 | 0.642 | 0.234 |

- reference: main spve=0.020 at R=500 (Stage 8.2): secondary 0.778, full 0.778, extra 0.486
- 结论核对点：secondary power 应随 spve 单调上升
