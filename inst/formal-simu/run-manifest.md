# CondPED Formal Simulation Run Manifest

> 补建于 Stage 8.2 完成后（2026-08-27）。参数口径以 freeze 方案
> `inst/formal-simu/05 模拟20260825.md` 为准。

## Freeze 参数（以 05 模拟20260825.md 为准）

| 参数 | 值 |
|---|---|
| scenario（主场景） | two_linked_trait_specific |
| n | 1000 |
| m | 4 |
| p | 1000 |
| h2 | 0.5 |
| locus_pve | 0.02 |
| secondary_signal_pve | 0.020 |
| target_r2 | 0.3 |
| rho_main | 0.05 |
| tau (tolerance) | 0.10 |
| correlation | block |
| q (K_Q) | 2 |
| reps（正式） | 500 |

## 运行记录

### Stage 8.1 pilot（2026-08-26）

- commit：`6d479b9`
- 范围：12 个主场景（I-0/I-1/I-2、II-A A1–A4、II-B R1/R3、E1/E2/E3）× R=20 = 240 reps
- workers = 4；主 grid wall = 2117s（总计 2153s）；峰值 RSS ≈ 104MB
- 结果：240/240 ok，0 numerical failure；§1.62 十项 pilot 检查全部 PASS
- 报告：stage81-pilot-report.md（本目录）；原始结果 output/stage81/（不进 git）
- 启动：`bash run_sim.sh stage81 inst/formal-simu/run-stage81.R`

### Stage 8.2 Simulation I 正式 500 reps（2026-08-27）

- commit：`c761ce7`
- 范围：I-0 null / I-1 single / I-2 two_linked_trait_specific × R=500 = 1500 reps
- workers = 4；wall = 14024s；峰值 RSS ≈ 104MB
- 结果：1500/1500 ok，0 numerical failure，失败重试 0 触发
- 关键结果：I-0 Type I 0.049；I-1 full 0.858 / attribution 0.951 /
  direction 0.974；I-2 full=secondary 0.778 / attribution 0.814 /
  direction(pattern) 0.820 / marker FDP 0.051 / signal_count_exact 0.649 /
  oracle-causal-set 0.98 (49/50)
- 报告：stage82-summary.md（本目录）；原始结果 output/stage82/（不进 git）
- 启动：`bash run_sim.sh stage82 inst/formal-simu/run-stage82.R`

## 基础设施说明

- `run_sim.sh` 已在 Stage 8.1 移植为本机可用的 `--user` scope 版本
  （commit 6d479b9）；所有正式批量运行必须经此入口（AGENTS.md 附加条款 1）。
- 所有正式模拟脚本在 `inst/formal-simu/` 根目录，结果在
  `inst/formal-simu/output/<阶段名>/`（gitignored）。
