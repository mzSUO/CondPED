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

## Stage 8.4（Simulation III，2026-08-30 完成）

- 场景：E1 single_highly_representable（++++, rho=0.05, conforming loop）、
  E2 single_nonredundant（runner 自然可行锥）、E3 linked_pseudo_multitrait
  （r2=0.3, spve=0.020），paired 三 pipeline，R=500 × 3 = 1500 reps，
  **1500/1500 ok，0 numerical failure**，wall 13695s，峰值 RSS 108MB。
- E1 vs E2：cand_exact 0.459/0.580，rep_exact 0.855/0.578，eta bias
  ~0.0003-0.0004，rho MAE 0.036/0.123，beta cov95 0.958/0.962，
  Q-form 强度匹配 0.0800/0.0800（§1.17）；符号分布 E1 全 ++++，
  E2 +--+ 45.2% / +-+- 54.6% / ++-- 0.2%（freeze v2 如实报告）。
- **E3 主证据（freeze 1.27 conditional effect contamination）**：
  真实 r2=0.302 下，causal1 对 Trait2 的 |污染| marginal 0.137 →
  resolved 0.059 → oracle-conditional 0.046；reduction 67%（→conditional）
  / 61%（→resolved）。trait-breadth inflation 仅作描述（lead 0.110 /
  resolved 0.184 / pseudo rate 0.079）。
- 产物：output/stage84/（不进 git）；提交副本 stage84-summary.md。

## Stage 8.3（Simulation II，2026-08-29 完成）

- 场景：II-A（A1–A4, signal_oracle）+ II-B（R1/R3, signal_trait_oracle），
  R=500 × 6 = 3000 reps，全部 ok，0 numerical failure。
- **freeze constraint fix v2**（commit `e856d26`）：可行性探测发现 R1/R3
  可行符号锥在冻结 Sigma_P_ref 下不相交（R1: ++++ 0.87；R3: +-+- 0.525 /
  +--+ 0.47，交集为空）→ freeze 文档 §1.13/1.16/1.33/1.42/1.43 改为
  "强度匹配 + 符号分布如实报告"；严格 ++++ 对比交由 I4 supplementary。
  I4 探测：两架构 100% 接受、++++ 各 13%（可行，待用户决定是否正式跑）。
- R1：500 个 concordant reps（可行锥主导模式，rejection 生成，
  n_sim_attempts 记录）；R3：runner 无约束 500 reps（接受率 ~100%，
  实际符号分布 +--+ 53.8% / +-+- 46.2%）。
- 关键结果：II-A CondPED exact 0.784–0.940 vs ASSET exact
  0.000–1.000（antagonistic ASSET 1.00 占优，broad ASSET 0.00）；
  II-B R1 rep_exact 0.722 / R3 0.994；eta bias ~0.0003、rho MAE
  0.035/0.055；Q-form 强度匹配 0.0800 vs 0.0800（§1.17 成立）。
- 事故记录：retry 的递归 glob 误删 R3 rejection-failed 留档（已修
  run-stage83.R 排除规则）；留档重建为 ARCHIVE.csv + 3 个确定性证据 RDS。
- 产物：output/stage83/（不进 git）；提交副本 stage83-summary.md。

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
- 补救补建（2026-08-27）：新增 inst/formal-simu/run-manifest.md（freeze
  参数 + Stage 8.1/8.2 运行记录）；stage82-summary.md 补两段文档——
  direction_exact 定义（pattern-of-attributed-set；0.820 vs 试点 ~1.00
  的差为指标口径差异，点估计符号准确率实为 1.000）与 I-2 oracle 漏斗
  （oracle-causal-set 0.98 → oracle-primary M49 0.86 / M1 1.00 →
  full 0.778）。不改任何结果数据。

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

## Stage 8.5（汇总与冻结验收，2026-08-30 完成）

- `inst/formal-simu/analyze-stage85.R`：汇总 stage82/83/84 全部 500-rep
  结果，所有比例类估计附 Wilson 95% CI（n=500），连续指标附 MC SE；
  主表 87 行（main_results_freeze.csv，提交副本 stage85-main-results.csv）。
- Figure 数据表：F2 = E3 contamination 梯度（marginal 0.137 → resolved
  0.059 → conditional 0.046，freeze §1.71/§1.27）；F3 = R1 vs R3 全指标
  + 符号分布（§1.72）；F4 = E1/E2/E3 三 pipeline 对比（§1.73）。
  指标分层按当前文档 §1.67/§1.68/§1.69 执行。
- **Freeze Criteria（当前文档 §1.63.1–1.63.6）全部 PASS**：
  PSD（pilot min eigen 0.16 + acceptance 全过）；0 rank-deficient /
  0 numerical failure（1500+3000+1500 全 ok）；R1/R3 truth map 完整；
  realized r2 0.302/0.300 ≈ 0.3；realized PVE 0.0200 ≈ 0.02、Q-form
  匹配 0.0800/0.0800；pipeline 严格配对（seed 25/25×3、结构 25/25×2）。
- 一致性核对：master_seed=20260826 体系三阶段一致（.seed_for_rep 逐点
  验证 75/75）；paired 三 pipeline 结构一致（50/50）。
- 报告：acceptance_report.md（提交副本 stage85-acceptance-report.md）；
  output/stage85/ 不进 git。章节号引用均按当前 v2 文档核实（§1.63
  criteria、§1.67–1.69 tiers、§1.71–1.73 figures、§1.64 pairing、
  §1.27 E3、§1.16 v2 符号报告）。

## 下一步（可选，待用户决定）

- II-B-I4 supplementary（Sigma_P=I4 严格 ++++ R1/R3 对比，探测已证可行）。
- Sensitivity 系列（n 梯度、spve 梯度、LD 梯度、rho 边界、q=3）。
