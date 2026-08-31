# Stage 9 Summary — Supplementary Sensitivity Series

All batches: deterministic seeds (master 20260826), per-rep RDS
checkpoints, workers = 4, one automatic retry, 0 numerical failure
hard standard. Section numbers refer to the current freeze document
(05 模拟20260825.md, constraint fix v2).

## stage9-spve

# Stage 9 batch 1: secondary-signal PVE sensitivity (freeze §1.56)

- spve grid {0.0075, 0.015, 0.0175}, R=200/cell, paired, 600/600 ok

| spve | locus det | secondary power | overall recovery | cand exact | direction | count_exact | extra rate |
|---|---|---|---|---|---|---|---|
| 0.0075 | 0.805 | 0.480 | 0.000 | 0.699 | 0.727 | 0.391 | 0.153 |
| 0.0150 | 0.790 | 0.675 | 0.000 | 0.741 | 0.747 | 0.639 | 0.205 |
| 0.0175 | 0.880 | 0.820 | 0.000 | 0.753 | 0.764 | 0.642 | 0.234 |

- reference: main spve=0.020 at R=500 (Stage 8.2): secondary 0.778, full 0.778, extra 0.486
- 结论核对点：secondary power 应随 spve 单调上升

## 补充：严格匹配口径的 recovery（与 Stage 8.2 一致）

| spve | full recovery | secondary recovery |
|---|---|---|---|
| 0.0075 | 0.470 | 0.470 |
| 0.0150 | 0.665 | 0.665 |
| 0.0175 | 0.770 | 0.770 |

## stage9-n

# Stage 9 batch 2: sample-size sensitivity (freeze §1.55)

- n in {500, 750, 1250}, R=200/cell, paired, 600/600 ok

| n | locus det | secondary power | overall recovery | cand exact | direction | count_exact | extra rate | mean rep runtime (s) |
|---|---|---|---|---|---|---|---|---|
| 500 | 0.285 | 0.200 | 0.000 | 0.439 | 0.439 | 0.351 | 0.068 | 16.1 |
| 750 | 0.605 | 0.535 | 0.000 | 0.607 | 0.620 | 0.521 | 0.132 | 25.5 |
| 1250 | 0.980 | 0.950 | 0.000 | 0.901 | 0.903 | 0.781 | 0.281 | 46.2 |

- reference: main n=1000 at R=500 (Stage 8.2): full/secondary 0.778
- 核对点：recovery 随 n 单调不降；runtime 随 n 增长可接受

- 严格匹配口径（Stage 8.2 约定）：full recovery 0.530, secondary 0.530

### 严格匹配口径按 n 分组

| n | full recovery | secondary recovery |
|---|---|---|---|
| 500 | 0.155 | 0.155 |
| 750 | 0.490 | 0.490 |
| 1250 | 0.945 | 0.945 |

## stage9-r2

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

## stage9-rho

# Stage 9 batch 4: rho sensitivity (freeze §1.59.1)

- R1 target_loss in {0.02, 0.08}, R=200/cell, 400/400 ok

| rho | rep_exact | rep_card | rho MAE | thr-side acc | cand exact |
|---|---|---|---|---|---|
| 0.02 | 0.930 | 0.930 | 0.0313 | 0.988 | 1.000 |
| 0.08 | 0.465 | 0.465 | 0.0384 | 0.916 | 1.000 |

- reference: main rho=0.05 (Stage 8.3, R=500): rep_exact 0.722, rho MAE 0.035
- 核对点：R1/R3 结论方向不随 rho 点翻转

## stage9-tau

# Stage 9 batch 5: tau sensitivity (freeze §1.59)

- R1 (rho=0.05) at tau in {0.05, 0.20}, R=200/cell, 200/200 ok

| tau | rep_exact | rep_card | rho MAE | thr-side acc |
|---|---|---|---|---|
| 0.20 | 0.995 | 0.995 | 0.0340 | 0.999 |

- reference: main tau=0.10 (Stage 8.3, R=500): rep_exact 0.722
- 核对点：tau=0.20 更宽 → rep_exact 应升（边界效应方向性）
- **tau=0.05 格结构性不可行**：R1 要求 truth rho < tau，rho_main=0.05 与
  tau=0.05 重合于边界，acceptance rule 全部拒绝（200/200 failed，
  本目录保留失败 rep 为证据）——tau 敏感性在 tau <= rho_main 处
  不是性能下降而是架构不存在

## stage9-boundary

# Stage 9 batch 6: boundary sensitivity (freeze §1.60)

- rho* in {0.08, 0.10, 0.12}, tau = 0.10, R=200/cell, 200/600 ok

## generation feasibility

| rho* | n | n_ok | feasible rate | fail codes |
|---|---|---|---|---|
| 0.08 | 200 | 200 | 1.000 |  |
| 0.10 | 200 | 0 | 0.000 | unstable |
| 0.12 | 200 | 0 | 0.000 | unstable |

## P(rho_hat crosses tau) on feasible cells

| rho* | n | |rho*-tau| | P(cross) |
|---|---|---|---|
| 0.08 | 200 | 0.02 | 0.535 |

- 参考点（Stage 8.3 主场景，rho*=0.05, tau=0.10, n=500）：P(cross) = 0.278
- 注意：rho* >= tau 的 cell 若 generation-infeasible，属于结构性事实：
  R1 架构按定义要求 truth rho < tau，边界外侧的真值不可生成——
  这正是 boundary sensitivity 要报告的机制边界

## stage9-q3

# Stage 9 batch 7: K_Q = 3 scalability (freeze §1.61)

- three_linked_trait_specific, R=200, 200/200 ok

- locus detection 0.995; overall recovery 0.000; count_exact 0.528; extra rate 0.348
- mean signals estimated 4.12; mean rep runtime 35.4s

- 定位（§1.61）：K_Q=3 只声明确算稳定与统计有效，不进主文
- 0 numerical failure = 计算稳定；locus detection 与 recovery 为统计有效性的 sanity 证据

## stage9-i3e4

# Stage 9 batch 8: I-3 two_heterogeneous + E4 mixed_multisignal (freeze §1.24 / §1.46)

- R=200 each, paired, 400/400 ok

| scenario | locus det | secondary power | overall recovery | cand exact | direction | count_exact | extra rate |
|---|---|---|---|---|---|---|---|
| mixed_multisignal | 0.990 | 0.875 | 0.222 | 0.715 | 0.816 | 0.667 | 0.457 |
| two_heterogeneous | 0.990 | 0.910 | 0.220 | 0.740 | 0.798 | 0.641 | 0.376 |

- 定位：I-3 证明 locus 内 multi-trait 与 trait-specific signal 共存；E4 只作展示
- E4 的 extra rate 预期偏高（novel-locus FDR 传播，Stage 7.5/8.5 已定性）

## stage9-i4

# Stage 9 batch 9: Sigma_P = I4 sensitivity (freeze §1.58)

- two_linked under correlation=independent, R=200, 200/200 ok

- locus detection 0.945; secondary power 0.750; overall recovery 0.000; cand exact 0.505; direction 0.521
- count_exact 0.407; extra rate 0.346; pseudo-multitrait rate 0.257

- 对照主场景（block cor, Stage 8.2）：pseudo rate 应在 I4 下明显下降；
  若不降，说明 pseudo-mixing 有 trait-correlation 之外的来源

- 严格匹配口径（Stage 8.2 约定）：full recovery 0.470, secondary 0.470

## stage9-r2arch

# Stage 9 batch 10: R2 partially_representable (freeze §1.14)

- R=200, 200/200 ok

- rep_exact 0.910; rep_card 0.930; rho MAE 0.0345; thr-side acc 0.990; Irr 0.910; cand exact 1.000

- 定位（§1.14/§1.32）：R2 放 supplementary，作为 R1/R3 主对比的中间架构

## stage9-fdp

# Stage 9 batch 11: FDR internal calibration (freeze §1.52)

- null / single / mixed, R=200 each, 600/600 ok

| scenario | n | marker FDP [Wilson 95%] | type1 | extra rate |
|---|---|---|---|---|
| mixed_multisignal | 200 | 0.044 [0.040, 0.049] | 0.093 | 0.412 |
| null | 200 | 0.867 [0.621, 0.963] | 0.049 | 0.065 |
| single_multi_trait | 200 | 0.079 [0.048, 0.126] | 0.050 | 0.052 |

- 正式 claim（§1.52）：marker FDP ≈ 0.05（含 CI）
- locus/signal FDP 不作为 primary error-control claim（marker BH 不自动产生 signal-level 控制）

## stage9-iib-i4

# Stage 9 batch 12: II-B-I4 strict same-sign comparison (freeze §1.58 ext.)

- Sigma_P = I4; R1/R3 with strict ++++ rejection, R=500 each, 1000/1000 ok

| scenario | n | rep_exact | rep_card | rho MAE | thr acc | Irr | cand exact | direction | eta bias | eta RMSE | mean attempts |
|---|---|---|---|---|---|---|---|---|---|---|---|
| highly_representable | 500 | 0.740 | 0.740 | 0.0222 | 0.958 | 0.740 | 1.000 | 0.588 | 0.0000 | 0.0444 | 8.4 |
| strongly_nonredundant | 500 | 0.790 | 0.790 | 0.0607 | 0.986 | 0.790 | 1.000 | 1.000 | -0.0010 | 0.0457 | 7.9 |

- empirical acceptance rate: highly_representable 0.119; strongly_nonredundant 0.127 (probe expected ~0.13; material deviation must be flagged)
- wall 10078 s
