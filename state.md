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

## 下一步

- Stage 8：正式 500-replicate simulation（按 05 模拟20260825.md 冻结参数，
  经 run_sim.sh 启动）。
