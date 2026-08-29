# Stage 8.3 II-B rep_exact 诊断说明（read-only 分析）

> 纯分析，不改任何代码与结果。数据：`output/stage83/repflip_diagnostics.csv`（1000 reps 逐条分解）与 `output/stage83/r1_singleton_rho_hat.png`。

## 1. R1 rep_exact = 0.722 的失败分解

500 reps 中 361 ok / 139 fail。失败类别（每条失败 rep 逐一归类）：

| 类别 | 数量 | 占比 |
|---|---|---|
| threshold_cross（singleton rho_hat ≥ tau=0.10 而 truth rho < tau） | **139** | **100%** |
| wrong_singleton_trait（cardinality 对但认错 trait） | 0 | 0% |
| cardinality_error_other | 0 | 0% |
| family_mismatch_other | 0 | 0% |

**结论：R1 的全部 rep_exact 失败都是 threshold misclassification**——即
truth singleton（{Trait1}，truth rho 固定 = 0.05）的估计值 ρ̂ 越过
tau=0.10 边界，导致该 rep 的 min Rep cardinality 被判为 ≥2。没有任何
cardinality 估计本身的错误（在 ρ̂ < tau 的 rep 中，cardinality 全部正确）。

singleton ρ̂ 分布（n=500，图见 output/stage83/r1_singleton_rho_hat.png）：

- median = 0.0714，q90 = 0.133，q99 = 0.200
- P(ρ̂ ≥ tau) = 0.278 —— 与 rep_exact 失败率（0.278）**逐一对应**（139/139）
- ρ̂ 相对 truth（0.05）存在约 +0.02 的上偏，叠加 SD ≈ 0.04–0.05 的估计噪声，使距 tau 仅 0.05 的真值有 ~28% 的概率被误判越界

R3 对照：497 ok / 3 fail（family_mismatch_other），无 threshold 问题——
R3 的 Rep=全集判定不涉及 singleton 边界，0.994 的恢复率与机制一致。

## 2. 与 Stage 7.8 预实验的配置差异（R1/R3）

| 配置项 | Stage 7.8 预实验 | Stage 8.3 正式 | 影响 |
|---|---|---|---|
| analysis mode | 自定义 full discovery pipeline（scan→loci→resolve→pipeline_condped_full） | signal_trait_oracle（truth-defined signal） | **关键差异**：7.8 的 rep_exact 混入了 discovery/resolution 层的失败；8.3 隔离纯 representation 层 |
| reps | 50 | 500 | 7.8 的 Monte Carlo 95% CI 约 ±0.14（R3）/±0.06（R1） |
| master_seed | 20260819 | 20260826 | 不同抽样 |
| R1 符号约束 | 无约束（混合符号 ~42%） | concordant-only（++++ 100%） | 改变 ρ̂ 分布形态 |
| R3 符号约束 | 无约束 | 无约束（自然可行锥 +--+ / +-+-） | 一致 |
| target_loss (R1) | 0.05 | 0.05 | 同 |
| n / m / p / locus_pve / tau / block cor | 同 | 同 | — |
| effect-scale matching（Q-form） | scale_match 开启 | scale_match 开启（实测 0.0800 vs 0.0800） | 同 |

### rep_exact 变化解释

- **R3: 0.56 → 0.994**：主要不是方法变化，而是**测量口径**。7.8 经
  discovery+resolution 全管线测量，信号解析失败的 rep 表现为 Rep 失败；
  8.3 在 truth-defined signal 下隔离 representation 层，R3 的真实
  representation 恢复率为 0.994。叠加 n=50 的宽置信区间，两个数字不矛盾。
- **R1: 0.94 → 0.722**：正式口径下的失败机制已由第 1 节定量锁定为
  **singleton ρ̂ 越界**（P(ρ̂≥tau)=0.278 与失败率一一对应）。0.94 是
  n=50、不同 seed 体系、经 discovery 口径（未评估 rep 不计入分母）且
  符号无约束的试点估计；其点估计偏高主要可由小样本波动与口径差异
  解释，不需要也不应据此调整 tau 或 scenario。R1 的 0.722 是 tau=0.10
  与 truth rho=0.05 相距过近这一**实验设计属性**的直接体现（与
  Stage 7.2/7.6 记录的 rho-near-tau boundary 现象一致），不是实现缺陷。

## 3. 对论文叙事的含义

- 报告中应将 R1 的 rep_exact=0.722 表述为 **threshold-side
  misclassification rate**（tau 边界效应），并同时报告 thr-side
  accuracy（0.955）与 rho MAE（0.035）——后两者显示 ρ̂ 本身估计良好，
  失败集中在离 tau 最近的决策边界上。
- 若希望 R1 rep_exact 更高，正确做法是增大 truth rho 与 tau 的 margin
  （属 sensitivity 设计选择），而不是修改 tau 或评价器。

## 附：复现

```bash
Rscript inst/formal-simu/analyze-stage83-repflip.R
# 输出: output/stage83/repflip_diagnostics.csv, r1_singleton_rho_hat.png
```
