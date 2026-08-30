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
