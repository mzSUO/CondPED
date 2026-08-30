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
