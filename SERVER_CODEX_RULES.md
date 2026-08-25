本轮所有 Shell 命令、代码修改、进程管理和计算资源使用均必须遵守该文件。
每条命令执行前执行其中的 Pre-command Gate。
如本轮任务指令与规则冲突，先报告冲突，不得自行绕过。

# 1. Server Codex Execution Rules
## 1.1. CondPED 项目服务器执行原则（兼容其他科研项目）

> **适用范围**：共享 Linux 服务器上的 CondPED 及其他科研计算项目。  
> **目标**：高吞吐、可复现、低风险、低 token 消耗；充分利用**可安全使用/已分配**的计算资源，同时避免 OOM、进程失控、误伤其他项目或用户。  
> **优先级**：安全与数据完整性 > 可复现性 > 正确性 > 吞吐 > 便利性。

---

## 1.2. 每次任务开始必须先做

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

若当前目录不是目标项目目录，**停止，不执行后续命令**。

CondPED 默认工作目录：

```text
/data2/smz/CondPED
```

CondPED 默认运行输出：

```text
/data2/smz/CondPED/inst/validation/output/
```

除非用户明确指定其他路径。

---

# 2. 每条 Shell 命令执行前的 Pre-command Gate

**每次执行命令前都必须快速检查以下 8 项。无需长篇输出，只需内部确认。**

1. **目录**：命令是否只作用于当前项目或明确允许的输出目录？
2. **范围**：是否可能递归影响父目录、其他项目、其他用户？
3. **Git**：是否会覆盖未提交修改或 tracked 文件？
4. **CPU**：并行度是否合理，是否存在 nested parallelism / oversubscription？
5. **内存**：预计峰值内存是否可能触发 OOM 或大量 swap？
6. **进程**：停止/重启命令是否只针对本任务明确 PID？
7. **输出**：结果是否写入正确、可 resume、不会覆盖旧结果的位置？
8. **可复现**：seed、scenario_id、rep_id、参数是否固定并可追踪？

若任一项不确定：**先检查，再执行；不要猜。**

---

# 3. 共享服务器安全边界

## 3.1. 允许范围

默认只能：

- 读取/修改当前项目目录；
- 写入当前项目明确的 output/log/tmp 目录；
- 查看本用户自己的进程；
- 查看系统总体资源状态；
- 启动/停止本任务创建的进程。

## 3.2. 禁止事项

除非用户明确授权，否则禁止：

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

也禁止：

- 修改其他用户文件；
- 删除其他项目文件；
- kill 不属于本任务的 PID；
- 修改系统级 R / Python / shared library；
- 修改全局 shell 配置；
- 在项目目录之外做代码修改；
- 扫描与当前任务无关的其他用户目录；
- 使用不受控的递归 `find /`、`grep -R /` 等全系统扫描。

需要访问项目外路径时：**停止并先向用户报告原因。**

---

# 4. CPU 使用原则：追求吞吐，不做资源失控

## 4.1. 并行优先级

对于大量独立 replicate：

> **优先 replicate-level parallelism。**

CondPED simulation 的 replicate 独立，优先使用现有 runner 的：

```r
workers = N
```

不要为了并行重新实现一套 RNG 或任务调度。

## 4.2. 防止 nested parallelism

外层已有多个 worker 时，内部 BLAS/OpenMP 默认设为单线程：

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
```

避免：

```text
N workers × 每个 worker M 个 BLAS threads
```

导致 CPU oversubscription、内存膨胀和吞吐下降。

## 4.3. worker 选择

没有调度器时，不无条件使用全部 CPU。

先检查：

```bash
nproc
free -h
uptime
ps -u "$USER" -o pid,ppid,%cpu,%mem,rss,etime,cmd --sort=-%cpu | head -25
```

优先做短 benchmark，例如：

```text
workers = 2, 4, 8, 16
```

必要时再测试 24/32。

选择标准：

> **单位墙钟时间完成的 replicate 数最多，且内存安全。**

不是 worker 数越多越好。

---

# 5. 内存管理原则：绝不以 OOM 换速度

## 5.1. 启动前

记录：

```bash
free -h
```

若是新 workload，先用少量真实 replicate 估计单 worker 峰值 RSS。

粗略安全条件：

```text
workers × peak_RSS_per_worker < 可安全使用内存
```

始终保留明显的系统/其他任务余量。

## 5.2. 运行中

定期检查：

```bash
free -h
ps -u "$USER" -o pid,ppid,%cpu,%mem,rss,etime,cmd --sort=-%mem | head -20
```

出现以下任一情况时降低 workers：

- available memory 持续快速下降；
- swap 持续增长；
- fork failure；
- OOM；
- worker RSS 异常增大；
- 增加 workers 后总吞吐反而下降。

## 5.3. 大对象原则

避免 master process 同时聚合大量 `Y/G/fit/K/large covariance/all replicate objects`。

每个 replicate 完成后：

1. 写盘；
2. 释放大对象；
3. 再进入下一个任务。

---

# 6. 进程管理

查看本用户进程优先：

```bash
ps -u "$USER" -o pid,ppid,%cpu,%mem,rss,etime,cmd --sort=-%cpu | head -30
```

停止任务时，只允许停止：

- 当前任务明确启动；
- PID 已确认；
- 命令行已核对；

的进程。

优先：

```bash
kill <PID>
```

必要时最后才使用：

```bash
kill -9 <PID>
```

禁止：

```bash
killall R
pkill R
pkill -u "$USER"
```

因为可能误杀本用户其他项目。

---

# 7. 后台任务规则

长时间任务必须：

- 有明确日志；
- 有明确输出目录；
- 可 resume；
- 能根据 PID 定位；
- 不依赖当前 SSH 会话的交互输入。

推荐至少保存：

```text
command
start time
PID
workers
seed/master_seed
scenario/grid
output path
log path
```

不要重复启动同一任务。启动前先检查本用户进程和目标 output。

---

# 8. 输出、临时文件与磁盘

## 8.1. 运行结果

大规模 simulation output 不进入 Git。

CondPED 默认忽略：

```text
inst/validation/output/
```

运行输出必须和 source code / formal configuration / validation scripts 分离。

## 8.2. Atomic save

长任务优先：

```text
write temp → successful close → rename to final file
```

防止中断后留下伪完整结果。

## 8.3. Resume

默认：

- `status=ok` → skip；
- temp/incomplete → rerun；
- failed → 保存失败信息；
- 是否重跑 failed 由显式参数决定。

---

# 9. Git 原则

修改前：

```bash
git status --short
git branch --show-current
```

提交前：

```bash
git diff --check
git diff --cached --check
git status --short
git diff --cached --stat
```

服务器默认：

- 不自动 commit；
- 不自动 push；
- 不自动 merge；
- 不自动 rebase；
- 不自动 reset hard。

运行结果目录不得 stage。

---

# 10. R / Python 环境原则

未经允许禁止系统级安装或全局环境修改，例如：

```text
sudo apt install
sudo yum install
pip install 到系统 Python
全局 conda 修改
```

缺包时：

1. 报告包名；
2. 优先项目/用户级环境；
3. 等用户决定。

CondPED 需要 ASSET 的任务先检查：

```bash
Rscript -e 'cat("ASSET:", requireNamespace("ASSET", quietly=TRUE), "\n")'
```

---

# 11. Token / 对话成本控制

Codex 应优先**执行而不是长篇解释**。

默认规则：

- 不重复用户已给出的背景；
- 不重复完整日志；
- 不把大文件整段打印到对话；
- 使用 `grep`, `sed -n`, `head`, `tail` 精准读取；
- 只读与当前问题相关的代码片段；
- 不一次性展开整个仓库；
- 不反复描述计划；
- 不生成无用 Markdown 报告；
- 过程中只汇报 blocker / 关键发现；
- 最终一次性汇总。

命令输出过长时优先：

```bash
... | tail -50
... | head -50
grep -n "pattern" file
sed -n '120,220p' file
```

避免 `cat` 巨大日志、打印整个大型对象、dump 整个 repository。

---

# 12. 修改代码的最小化原则

优先：

> **最小修改、局部修复、已有函数复用。**

禁止为了“顺手优化”而：

- 大规模重构；
- 更换接口；
- 重命名大量字段；
- 新写一套已有统计量；
- 改动与当前任务无关模块。

每个 bug fix 都应尽可能附带最小 regression test。

---

# 13. CondPED 项目专属冻结规则

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

诊断发现性能问题时：

1. 先定位层级；
2. 先做 oracle / controlled diagnostic；
3. 区分 simulation scenario / discovery / locus construction / signal resolution / attribution / rho-Rep；
4. 报告证据；
5. 不自动改核心算法。

---

# 14. CondPED simulation 特别原则

## 14.1. RNG

每个 replicate 只能由：

```text
master_seed
canonical scenario_id
rep_id
```

决定。

必须保持：

```text
grid reorder invariant
workers invariant
sequential == parallel
resume invariant
```

禁止以 worker id、grid row、execution order、wall-clock time 决定 seed。

## 14.2. Sim III paired comparison

同一 replicate 的：

```text
lead_asset
resolved_asset
condped_full
```

必须共享：

```text
same Y
same G
same truth
same seed
same simulation identity
```

不得分别重新模拟。

## 14.3. Runner 职责

Runner 只负责：

```text
grid
seed
simulate
analyze
evaluate
save
failure
resume
parallel
```

不得在 runner 内重新实现 truth / eta / rho / Rep / Irr / ASSET statistic / signal matching。

---

# 15. 每次命令执行前的极简自检格式

为了节省 token，不需要每次向用户展开 8 项检查。

Codex 在内部按如下顺序确认：

```text
PATH → GIT → SCOPE → CPU → MEM → PID → OUTPUT → RNG
```

若全部安全，直接执行。

若任一项有风险，只输出一句：

```text
BLOCKED: <具体风险 + 需要用户决定的事项>
```

不要输出长篇安全说明。

---

# 16. 长计算任务的执行模式

开始前：

```text
1. environment check
2. small smoke test
3. worker benchmark
4. choose safe workers
5. launch resumable run
6. monitor memory/CPU
7. summarize
```

不要一上来直接 500 reps、开最大 workers 或重写 runner。

---

# 17. 完成一个任务后的最小汇报

默认只报告：

```text
STATUS
CHANGED
RAN
RESULT
RESOURCE
TEST
BLOCKER
NEXT
```

示例：

```text
STATUS: PASS
CHANGED: 2 files
RAN: 100 reps, workers=16
RESULT: secondary recovery 0.82
RESOURCE: 42 min, peak RSS 19 GB
TEST: 1338 PASS
BLOCKER: none
NEXT: freeze scenario
```

避免重复粘贴完整日志。

---

# 18. 出现异常时

## 18.1. OOM / swap 暴涨

立即：

1. 停止继续提交新 worker；
2. 等已运行 worker 安全退出；
3. 降低 workers；
4. 保留失败日志；
5. 不删除其他任务。

## 18.2. 输出异常

先检查：

```text
seed
scenario_id
rep_id
input identity
resume status
```

再怀疑统计方法。

## 18.3. Git 冲突

停止自动处理并报告。

## 18.4. 路径不确定

停止，不执行破坏性命令。

---

# 19. 本文件的执行要求

**每次开始新任务时先读取本文件。**

**每次执行 Shell 命令前按第 1 节的 Pre-command Gate 检查。**

**每次启动长计算前按第 15 节执行。**

若用户临时指令与本文件冲突：用户明确的新指令优先，但不得突破共享服务器安全、文件权限和数据完整性底线。

---

## 19.1. 一句话原则

> **只碰自己的项目；先确认再执行；独立 replicate 用并行；外层并行时 BLAS 单线程；CPU 追求吞吐、内存留余量；只杀自己的 PID；结果可 resume；核心统计不因结果难看而擅自修改；尽量少读、少说、少耗 token。**
