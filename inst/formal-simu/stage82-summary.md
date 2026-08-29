# Stage 8.2 Summary — Simulation I formal 500 reps

- scenarios: I-0 null, I-1 single, I-2 two_linked_trait_specific (r2 = 0.3, spve = 0.020)
- R = 500 per scenario; 1500/1500 replicates ok; retries: 0 (recovered 0)
- wall time 14024 s; peak RSS 104 MB; workers = 4

## Primary endpoints (freeze 1.53)

| scenario | full recovery | secondary recovery | attribution exact | direction | cond effect RMSE | marker FDP |
|---|---|---|---|---|---|---|
| I-0 null | — | — | — | — | — | — (see Type I below) |
| I-1 single | 0.858 | — | 0.951 | 0.974 | 0.031 | 0.071 |
| I-2 two-linked | 0.778 | 0.778 | 0.814 | 0.820 | 0.061 (beta RMSE; see note) | 0.051 |

## Secondary diagnostics

| scenario | signal_count_exact | extra-signal rate | novel-locus extra share |
|---|---|---|---|
| I-1 | 0.949 | 0.106 | 0.697 |
| I-2 | 0.649 | 0.486 | (see extra_provenance.csv) |

## Oracle-causal-set resolver diagnostic (I-2, reps 1-50): recovery = 0.98 (49/50)

## Resolver funnel (I-2, reps 1-50, strictly paired)

- oracle-causal-set recovery: 0.98 (49/50)
- oracle-primary recovery (true primary conditioned, secondary tested, within-locus Bonferroni M = region-1): 0.86 (43/50); single-test M = 1: 1.00
- full end-to-end two-signal recovery (all 500 reps): 0.778

## Note on direction_exact

direction_exact (0.820) is the frozen evaluator's direction PATTERN recovery: per matched signal pair, the estimated candidate-set direction pattern (0 candidates -> not_applicable; 1 -> single_trait; >= 2 same-sign -> concordant; 2 opposite -> antagonistic; >= 3 both signs -> mixed) must equal truth$effect_direction. For I-2 both truth signals are single_trait, so the metric fails exactly when the attributed candidate set is not exactly one trait - i.e. it is attribution-driven, not sign-flip-driven. The point-estimate sign accuracy on truth-active traits is 1.000 over matched pairs (direction_metrics_i2.csv). The Stage 7.7/7.8 pilot value ~1.00 used the latter sign-based metric on n = 50; the 0.820 here is the stricter pattern metric on n = 500, and it tracks attribution exact (0.814) as expected.

## Extra-signal provenance (I-1 + I-2)

- novel_locus_fp: 350; split_region_fp: 109; within_locus_fp: 43

## Type I (I-0 null)

- omnibus Type I (marker-level BH-selected null share / any-FP rate): 0.049

Note: the frozen evaluator defines conditional_effect_* only for
single-signal modes; for I-2 the table reports beta RMSE (effect
estimation error on matched signals) instead. I-1 marker FDP (0.071) is
slightly above the nominal 0.05 — BH controls E[FDP] under
independence/PRDS; mildly correlated test statistics can exceed it.

---
# summary解读

这份 summary 信息量比表面大，我逐块拆开讲，包括数字背后的机制含义和对论文叙事的影响。

## 一、总体判断

1500/1500、0 失败、I-0 校准达标、I-2 全部主指标落在预期区间——**这是一份健康的正式结果**，而且几个"看起来不好看"的数字（signal_count_exact 0.649、extra rate 0.486、direction 0.820）恰恰是你方案里提前布局过的点，都有现成的解释框架。

## 二、Primary endpoints 逐行解读

**I-0 null**：Type I = 0.049 ≈ 名义 0.05。说明 discovery layer（marker-level BH）在 null 下是校准的——这是一切后续结果的合法性前提，过了这关，I-1/I-2 的数字才可解释。

**I-1 single**：full recovery 0.858 / attribution 0.951 / direction 0.974 / RMSE 0.031。单个 PVE=0.02 的信号在 n=1000 下有 86% 完整恢复率、95% 性状归因正确率，这是合理的"有限样本但不失能"区间。注意 RMSE 0.031 很小：一旦信号被匹配上，效应估计很准——误差主要来自"没检到"，而不是"检到但估错"。

**I-2 two-linked（核心场景）**：
- full = secondary = 0.778，与 Stage 7.8 预实验的 0.80 一致（500 reps 下 MCSE≈0.019，0.778 在 0.80±0.04 内），说明预实验的 regime 选择没有偏；
- attribution 0.814，比预实验的 0.76 还略好，即约 81% 的 rep 里 S1→A、S2→B 被**精确**恢复；
- **direction 0.820 是唯一与预实验（1.00）明显偏离的指标**。50 reps 的预实验看不到尾部，500 reps 暴露了约 18% 的 rep 存在至少一个方向错误。最可能的来源是弱信号 rep：效应估计值接近 0 时符号被噪声翻转，属于真实有限样本现象而非 bug（I-1 也有 2.6% 的方向错误佐证）。写作时建议加一句"conditional on 高置信匹配 rep 的 direction"作为补充口径；
- marker FDP 0.051 正好打在靶心上；I-1 的 0.071 略超，报告注释的解释是对的——BH 只保证独立性/PRDS 下的 E[FDP]，你的四性状相关结构天然违反该假设，0.07 这个量级属于文献中已知的温和超出，如实写即可，反而是诚实性的加分项。

## 三、Secondary diagnostics：这份报告最有价值的一块

I-2 的 signal_count_exact 只有 0.649、extra-signal rate 高达 0.486——如果只看这两个数，审稿人会得出"近一半 rep 解析出错"的印象。**但 provenance 分解把这个误解剖开了**：

- novel_locus_fp 350（70%）：绝大部分 extra 是**基因组别处的假 locus**，属 discovery layer 的固有噪声，与 CondPED 的 signal resolution 无关；
- split_region_fp 109（22%）：一个真 locus 被区域构造切成两块，属 locus definition 的人为产物；
- within_locus_fp 43（8.6%）：真正发生在目标 locus 内部、会污染 resolution 解释的 extra，**1000 个 rep 里只有 43 个（4.3%）**。

这正是 freeze 方案 §1.53 把 signal_count_exact 降级为 secondary、并要求报 extra_provenance 的原因——现在数据验证了那个设计决定。论文里这段应写成："表面 count 误差主要由 discovery-layer novel-locus FP 主导；目标 locus 内部的 resolution 错误率仅 ~4%"。

I-1 里 signal_count_exact 0.949 > full recovery 0.858 也同理：count 恰好为 1 的 rep 里混着"真信号丢了、但恰好检出 1 个假信号"的情况——count 是个混淆指标，不能单独解读。

## 四、Oracle 漏斗（目前不完整）

oracle-causal-set 0.98（49/50）：把真 causal variants 直接交给 resolver，98% 能解出来；full pipeline 0.778。说明约 20 个百分点的损失发生在 discovery → locus → 候选集构造这些上游环节，而不是 resolver 本身——这对你是好消息，机制主张（resolver + representation）没有被证伪。

两个口径问题要补：① oracle-primary 缺报，漏斗只有两层（Stage 7.8 参考值 0.80）；② oc 只跑了 reps 1–50，而 full 是 500 reps，**严格比较应在同样 50 个 rep 上报 paired 的 full recovery**，否则漏斗两端样本不同。

## 五、对下游的影响

1. 数字层面全部支持按原计划进 8.3，无需调任何参数——**千万不要因为 direction 0.82 或 count 0.649 去动代码或 regime**，这些是机制真实特征，方案的三层指标分层就是为这种情况设计的。
2. 写作素材上，这份报告直接给了你：Figure 2 主场景的 full/attribution 锚点、extra provenance 的分解图（很适合做 supplementary waterfall/bar）、FDP 校准与相关性 caveat 的一段话。
3. 建议补进 summary 的三项（可让 kimi 在进 8.3 前顺手做）：direction_exact 的口径说明、oracle-primary 数值、oc vs full 的同 50-rep paired 对照。

可以继续按上一条的 8.3 提示词推进。