# CondPED state

## 当前状态（2026-08-26）

- 分支 dev/hierarchical；统计定义已冻结（Methods > Interface Contract > Runbook > 代码）。
- 工作规范：见 AGENTS.md（已合并 SERVER_CODEX_RULES.md 全部规则 + 三条附加条款：
  批量模拟必须经 run_sim.sh；禁止修改 production code / simulation registry /
  evaluator；Stage 8 起正式模拟脚本放 inst/formal-simu/ 根目录、结果放
  inst/formal-simu/output/<阶段名>/（已加入 .gitignore，只 commit 脚本与汇总报告）；
  freeze 参数以 inst/formal-simu/05 模拟20260825.md 为准）。
- 测试：devtools::test() 1513 PASS / 0 FAIL / 0 WARN / 1 SKIP；
  devtools::check 0 errors / 0 warnings。

## 已完成的关键工作

- Stage 7.3.5：scan Q 跨平台数值一致性——`.safe_inverse` 加 scale-aware
  zero floor（intercept 共线 marker → rank 0），测试 contract 更新。
- Stage 7.4：same-locus causal geometry 修正（q≥2 causals 改为 region 中心
  1kb cluster，原为 49kb 两端 vs 5kb window）；secondary PVE 校准
  0.0075→0.0175；fit_mt_null REML 溢出 guard；真实 ASSET 验证通过。
- Stage 7.5：full→oracle 失败归因审计；signal_count_exact 被 novel-locus
  FDR 系统性压低的机制确认。
- Stage 7.6：方法/模拟/评价器三层归因；oracle 阶梯证明 resolver 无损，
  无 C 类方法问题。
- Stage 7.7：operating-regime 校准（PVE/n/LD/rho/q 网格，各 50 reps）。
- Stage 7.8：operating-regime audit，推荐 formal baseline：
  two_linked_trait_specific, n=1000, locus_pve=0.02, spve=0.020, r2=0.3,
  rho=0.05, q=2（详见 inst/validation/output/stage78/stage78_final_freeze.csv）。
- Rcpp 迁移：scan_mt_omnibus per-SNP GLS block 迁至 RcppArmadillo
  （src/gls_blocks.cpp）；R 参考保留为 .gls_block_components_r()；
  金标准测试 tests/testthat/test-gls-block-cpp.R（154 个）；
  benchmark n=800/p=100k：~208s→26.9s（×8），数值逐位一致；
  Stage 7.8 复跑结果与迁移前 bitwise 一致。

## Stage 8.2（2026-08-27，已完成）

- `inst/formal-simu/run-stage82.R`：Simulation I 正式 500 reps——I-0 null /
  I-1 single / I-2 two_linked_trait_specific（r2=0.3，spve=0.020），paired
  三 pipeline 严格配对，deterministic seed（master 20260826），per-rep RDS
  checkpoint，每 50 reps 落 partial RDS + 进度打印，失败自动重试一次
  （0 失败未触发）。经 `bash run_sim.sh stage82 ...` 启动，4 workers，
  wall 14024s，峰值 RSS 104MB，**1500/1500 ok，0 numerical failure**。
- `inst/formal-simu/analyze-stage82.R`：freeze §1.53 指标汇总。
  关键结果：I-0 Type I 0.049（BH 达标）；I-1 full recovery 0.858、
  attribution 0.951、direction 0.974；I-2 full/secondary recovery 0.778、
  attribution 0.814、direction 0.820、marker FDP 0.051、
  signal_count_exact 0.649、oracle-causal-set 0.98（49/50）；
  extra provenance：novel_locus 350 / split_region 109 / within_locus 43。
- 汇总：inst/formal-simu/output/stage82/summary.md（提交副本
  inst/formal-simu/stage82-summary.md）；output/ 不进 git。
- 注：evaluator 的 conditional_effect_* 仅 single-signal 模式有定义，I-2
  报告 beta RMSE；I-1 marker FDP 0.071 略高于 0.05（相关检验下 BH 已知
  现象，已在报告中注明）。
- run-manifest.md 至今不存在（Stage 8.1/8.2 均以 freeze 文档为准）。

## Stage 8.1 pilot（2026-08-26，已完成）

- `inst/formal-simu/run-stage81.R`：12 个主场景（I-0/I-1/I-2、II-A A1–A4、
  II-B R1/R3、E1/E2/E3）× R=20 = 240 reps，经 `bash run_sim.sh stage81 ...`
  启动；三 pipeline 严格配对（runner paired 模式，deterministic seed）；
  per-rep RDS checkpoint；失败自动重试一次（本次 0 失败，retry 未触发）。
- run_sim.sh 已修复为本机可用的 `--user` scope 版本（原 sudo/--wait/Nice
  与本机 systemd 不兼容，且 sudo 路径会产生 root -owned 输出与 PATH 缺失）。
- pilot 结果：240/240 ok；10 项 freeze 检查（§1.62）全部 PASS
  （realized LD 0.300；fixed-dirs 场景 per-active-trait realized PVE 0.0200；
  Sigma_G^bg PSD；truth rho map / Rep/Irr truth 正常；ASSET 适配器
  53/53 运行全 ok，7 个 no-locus rep 属 discovery 结果；resume/probe
  再生一致；0 numerical failure）。主 grid wall 2117s，峰值 RSS 107MB。
- 报告：inst/formal-simu/stage81-pilot-report.md（提交）；
  原始结果在 inst/formal-simu/output/stage81/（gitignored，不提交）。
- 满足 freeze 硬性条件（0 numerical failure），可进入正式 500-rep。
