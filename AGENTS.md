# AGENTS.md — CondPED

本文件是所有代码代理在本仓库工作时的必读规范。内容与 `SERVER_CODEX_RULES.md`
完全一致并合并了其全部规则；两者冲突时以本文件为准（本文件即其超集）。

---

## 1. 项目概览

CondPED（Conditional Projection for Decomposing Genetic Effects in
Multi-Trait GWAS）是一个 R 包：对已发现的多性状关联 locus 做两层条件解析——
variant 维度上区分同一 locus 内的多个独立信号（conditional signal
resolution），trait 维度上做 trait attribution 与条件可表示性分解
（eta / rho / Rep / Irr）。比较器为 ASSET（可选依赖）。

- 当前版本 v1.0（DESCRIPTION 0.1.0），分支 `dev/hierarchical`。
- 统计定义已冻结，优先级：Methods > Interface Contract > Runbook > 代码。
- 全基因组扫描热点路径（`scan_mt_omnibus()` 的 per-SNP GLS block）已迁移至
  Rcpp/RcppArmadillo（`src/gls_blocks.cpp`），R 参考实现保留为
  `CondPED:::.gls_block_components_r()`，金标准测试见
  `tests/testthat/test-gls-block-cpp.R`。

## 2. 环境与常用命令

```bash
# 共享服务器上所有 R 进程默认单线程 BLAS（外层并行时禁止嵌套并行）：
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

# 加载 / 测试 / 文档 / 检查
Rscript -e 'devtools::load_all(".", quiet=TRUE)'
Rscript -e 'devtools::test()'                    # 全量；可用 filter="xxx"
Rscript -e 'devtools::test(filter="gls-block-cpp")'
Rscript -e 'Rcpp::compileAttributes(); devtools::document()'   # 改 C++ 或 roxygen 后
Rscript -e 'devtools::check(document=FALSE, manual=FALSE)'

# ASSET 依赖检查（需要 ASSET 的任务先跑）：
Rscript -e 'cat("ASSET:", requireNamespace("ASSET", quietly=TRUE), "\n")'

# 正式批量模拟（强制入口，见 §4 附加规则 1）：
bash run_sim.sh <标签> <脚本>
```

注意：本服务器 NFS 存在时钟偏差，改 C++ 后若出现 undefined symbol /
stale object，先 `rm -f src/*.o src/*.so && touch src/*.cpp` 再重建。
devtools::document() 重建时若遇 `pandoc` 缺失，profvis 报告用
`selfcontained = FALSE` 保存。

## 3. 仓库结构

```
R/                  # 包源码（统计方法 + runner + evaluator）
src/                # Rcpp/Armadillo kernel（gls_blocks.cpp + RcppExports）
tests/testthat/     # testthat 第 3 版
inst/validation/    # 诊断 / 校准 / pilot 脚本（output/ 为其输出目录，gitignored）
inst/formal-simu/   # 正式模拟设计文档与正式批量脚本
man/                # roxygen2 生成，勿手改
run_sim.sh          # 正式批量模拟的唯一启动入口
SERVER_CODEX_RULES.md  # 服务器执行规则（已并入本文件 §4）
```

当前关键事实（Stage 7.3.5–7.8 与 Rcpp 迁移后）：

- 全套测试 1513 PASS / 0 FAIL / 0 WARN / 1 SKIP（SKIP 为 ASSET 已安装的
  预期分支）。
- 正式 500-rep 推荐 baseline：`two_linked_trait_specific`，n=1000，m=4，
  p=1000，locus_pve=0.02，secondary_signal_pve=0.020，target_r2=0.3，
  rho=0.05，tau=0.10，block cor，q=2。**一切以
  `inst/formal-simu/05 模拟20260825.md` 为准**（见 §4 附加规则 3）。

---

## 4. 服务器与项目执行规则（SERVER_CODEX_RULES.md 全文合并 + 附加条款）

### 附加条款（本文件新增，优先级同最高）

1. **正式批量模拟必须通过 `bash run_sim.sh <标签> <脚本>` 启动**，
   禁止直接 `Rscript` 跑批量任务。
2. **任何情况下不修改 production code、simulation registry、
   evaluator**；只能新增 `inst/validation/` 与 `inst/formal-simu/` 下的
   脚本和 `output/` 下的结果。（历史 stage 的代码修改均为当时用户明确
   授权，不构成先例。）
3. **正式模拟的 freeze 方案文档是 `inst/formal-simu/05 模拟20260825.md`，
   所有参数以它为准。**

### 4.1. 适用范围与优先级

> 适用范围：共享 Linux 服务器上的 CondPED 及其他科研计算项目。
> 目标：高吞吐、可复现、低风险、低 token 消耗；充分利用可安全使用/已分配
> 的计算资源，同时避免 OOM、进程失控、误伤其他项目或用户。
> 优先级：安全与数据完整性 > 可复现性 > 正确性 > 吞吐 > 便利性。

### 4.2. 每次任务开始必须先做

在任何修改、运行或后台计算前，先确认：

```bash
pwd
git status --short
git branch --show-current
git log -1 --oneline
nproc
free -h
uptime
```

若当前目录不是目标项目目录，停止，不执行后续命令。

CondPED 默认工作目录：`/data2/smz/CondPED`
CondPED 默认运行输出：`/data2/smz/CondPED/inst/validation/output/`
除非用户明确指定其他路径。

### 4.3. 每条 Shell 命令执行前的 Pre-command Gate

每次执行命令前快速检查以下 8 项（内部确认即可）：

1. 目录：命令是否只作用于当前项目或明确允许的输出目录？
2. 范围：是否可能递归影响父目录、其他项目、其他用户？
3. Git：是否会覆盖未提交修改或 tracked 文件？
4. CPU：并行度是否合理，是否存在 nested parallelism / oversubscription？
5. 内存：预计峰值内存是否可能触发 OOM 或大量 swap？
6. 进程：停止/重启命令是否只针对本任务明确 PID？
7. 输出：结果是否写入正确、可 resume、不会覆盖旧结果的位置？
8. 可复现：seed、scenario_id、rep_id、参数是否固定并可追踪？

若任一项不确定：先检查，再执行；不要猜。

### 4.4. 共享服务器安全边界

允许：读取/修改当前项目目录；写入当前项目明确的 output/log/tmp 目录；
查看本用户自己的进程与系统总体资源状态；启动/停止本任务创建的进程。

除非用户明确授权，禁止：

```text
sudo
su
rm -rf
chmod -R
chown -R
killall
pkill -u
git clean -fdx
git reset --hard
```

也禁止：修改其他用户文件；删除其他项目文件；kill 不属于本任务的 PID；
修改系统级 R / Python / shared library；修改全局 shell 配置；在项目目录
之外做代码修改；扫描与当前任务无关的其他用户目录；使用不受控的递归
`find /`、`grep -R /` 等全系统扫描。

需要访问项目外路径时：停止并先向用户报告原因。

### 4.5. CPU 使用原则：追求吞吐，不做资源失控

- 大量独立 replicate 优先 replicate-level parallelism，使用现有 runner 的
  `workers = N`；不要为了并行重新实现 RNG 或任务调度。
- 外层已有多个 worker 时，内部 BLAS/OpenMP 默认单线程（见 §2 的环境变量），
  避免 N workers × M BLAS threads 的 oversubscription。
- 不无条件使用全部 CPU。先检查 `nproc` / `free -h` / `uptime` /
  `ps -u "$USER" -o pid,ppid,%cpu,%mem,rss,etime,cmd --sort=-%cpu | head -25`。
  优先做短 benchmark（workers = 2, 4, 8, 16），必要时再测 24/32。
  选择标准：单位墙钟时间完成的 replicate 数最多，且内存安全。
  不是 worker 数越多越好。

### 4.6. 内存管理原则：绝不以 OOM 换速度

- 启动前记录 `free -h`；新 workload 先用少量真实 replicate 估计单 worker
  峰值 RSS。粗略安全条件：workers × peak_RSS_per_worker < 可安全使用内存，
  并始终保留明显的系统/其他任务余量。
- 运行中定期检查 `free -h` 与
  `ps -u "$USER" -o pid,ppid,%cpu,%mem,rss,etime,cmd --sort=-%mem | head -20`。
  出现以下任一情况时降低 workers：available memory 持续快速下降；swap
  持续增长；fork failure；OOM；worker RSS 异常增大；增加 workers 后总吞吐
  反而下降。
- 避免 master process 同时聚合大量大对象；每个 replicate 完成后写盘、
  释放大对象、再进入下一个任务。

### 4.7. 进程管理

- 查看本用户进程优先：`ps -u "$USER" -o pid,ppid,%cpu,%mem,rss,etime,cmd
  --sort=-%cpu | head -30`。
- 只允许停止当前任务明确启动、PID 已确认、命令行已核对的进程。
- 优先 `kill <PID>`；必要时最后才 `kill -9 <PID>`。
- 禁止 `killall R` / `pkill R` / `pkill -u "$USER"`（可能误杀其他项目）。

### 4.8. 后台任务规则

长时间任务必须：有明确日志；有明确输出目录；可 resume；能根据 PID 定位；
不依赖当前 SSH 会话的交互输入。推荐至少保存：command、start time、PID、
workers、seed/master_seed、scenario/grid、output path、log path。
不要重复启动同一任务；启动前先检查本用户进程和目标 output。

### 4.9. 输出、临时文件与磁盘

- 大规模 simulation output 不进入 Git；CondPED 默认忽略
  `inst/validation/output/`；运行输出必须和 source code / formal
  configuration / validation scripts 分离。
- 长任务优先 atomic save：write temp → successful close → rename。
- Resume 默认：status=ok → skip；temp/incomplete → rerun；failed → 保存
  失败信息；是否重跑 failed 由显式参数决定。

### 4.10. Git 原则

修改前：`git status --short`、`git branch --show-current`。
提交前：`git diff --check`、`git diff --cached --check`、
`git status --short`、`git diff --cached --stat`。

服务器默认：不自动 commit；不自动 push；不自动 merge；不自动 rebase；
不自动 reset hard。运行结果目录不得 stage。

### 4.11. R / Python 环境原则

未经允许禁止系统级安装或全局环境修改（sudo apt/yum install、pip 到系统
Python、全局 conda 修改）。缺包时：报告包名；优先项目/用户级环境；等用户
决定。CondPED 需要 ASSET 的任务先检查
`Rscript -e 'cat("ASSET:", requireNamespace("ASSET", quietly=TRUE), "\n")'`。

### 4.12. Token / 对话成本控制

优先执行而不是长篇解释。不重复用户已给出的背景；不重复完整日志；不把
大文件整段打印到对话；使用 grep / sed -n / head / tail 精准读取；只读与
当前问题相关的代码片段；不一次性展开整个仓库；不反复描述计划；不生成
无用 Markdown 报告；过程中只汇报 blocker / 关键发现；最终一次性汇总。
命令输出过长时优先 `... | tail -50` / `head -50` / `grep -n` /
`sed -n '120,220p'`。

### 4.13. 修改代码的最小化原则

优先最小修改、局部修复、已有函数复用。禁止为了"顺手优化"而大规模重构、
更换接口、重命名大量字段、新写一套已有统计量、改动与当前任务无关模块。
每个 bug fix 都应尽可能附带最小 regression test。

### 4.14. CondPED 项目专属冻结规则

除非用户明确批准，以下核心定义视为冻结：

```text
omnibus statistic
BH
Holm trait attribution
conditional score statistic
within-locus Bonferroni
local conditional threshold
eta
rho
Rep
Irr
tau
ASSET implementation
signal matching definition
Sigma_P_ref trait-space role
```

诊断发现性能问题时：先定位层级；先做 oracle / controlled diagnostic；
区分 simulation scenario / discovery / locus construction / signal
resolution / attribution / rho-Rep；报告证据；不自动改核心算法。

### 4.15. CondPED simulation 特别原则

RNG：每个 replicate 只能由 master_seed、canonical scenario_id、rep_id
决定。必须保持：grid reorder invariant、workers invariant、
sequential == parallel、resume invariant。禁止以 worker id、grid row、
execution order、wall-clock time 决定 seed。

Sim III paired comparison：同一 replicate 的 lead_asset /
resolved_asset / condped_full 必须共享 same Y / same G / same truth /
same seed / same simulation identity，不得分别重新模拟。

Runner 职责：grid、seed、simulate、analyze、evaluate、save、failure、
resume、parallel。不得在 runner 内重新实现 truth / eta / rho / Rep /
Irr / ASSET statistic / signal matching。

### 4.16. 每次命令执行前的极简自检格式

内部按 `PATH → GIT → SCOPE → CPU → MEM → PID → OUTPUT → RNG` 顺序确认。
若全部安全，直接执行；若任一项有风险，只输出一句
`BLOCKED: <具体风险 + 需要用户决定的事项>`，不要输出长篇安全说明。

### 4.17. 长计算任务的执行模式

开始前：environment check → small smoke test → worker benchmark →
choose safe workers → launch resumable run → monitor memory/CPU →
summarize。不要一上来直接 500 reps、开最大 workers 或重写 runner。

### 4.18. 完成一个任务后的最小汇报

默认只报告：STATUS / CHANGED / RAN / RESULT / RESOURCE / TEST /
BLOCKER / NEXT。避免重复粘贴完整日志。

### 4.19. 出现异常时

- OOM / swap 暴涨：立即停止提交新 worker；等已运行 worker 安全退出；
  降低 workers；保留失败日志；不删除其他任务。
- 输出异常：先检查 seed / scenario_id / rep_id / input identity /
  resume status，再怀疑统计方法。
- Git 冲突：停止自动处理并报告。
- 路径不确定：停止，不执行破坏性命令。

### 4.20. 一句话原则

> 只碰自己的项目；先确认再执行；独立 replicate 用并行；外层并行时
> BLAS 单线程；CPU 追求吞吐、内存留余量；只杀自己的 PID；结果可
> resume；核心统计不因结果难看而擅自修改；尽量少读、少说、少耗 token。
