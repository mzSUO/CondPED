# CondPED 函数接口参考

- 包：**CondPED** v0.1.0 —— Conditional Projection for Decomposing Genetic Effects in Multi-Trait GWAS
- 核心流程：三层工作流 —— Layer 1 边际 OLS 扫描 → Layer 2 条件投影 → Layer 3 成对双向 MR → 关联模式分类
- 统计口径：NAMESPACE 共导出 **32 个函数**；包内另定义约 **20 个内部辅助函数**（未导出）
- 依赖：`mvtnorm`（Imports）；LMM 分支可选依赖 `lme4` / `lmerTest`

## 模块与源文件对照

| 模块 | 源文件 | 导出函数数 | 内部函数数 |
|---|---|---|---|
| 主流程入口 | `R/condped.R` | 1 | 1 |
| Layer 1 边际扫描 | `R/gwas_marginal.R` | 2 | 1 |
| Layer 2 条件投影 | `R/conditional.R` | 4 | 3 |
| Layer 3 MR | `R/mr_v2.R` | 9 | 13（含 `bidirectional_mr`，见文末备注） |
| 模式分类 | `R/classify_v2.R` | 3 | 5 |
| 数据模拟 | `R/simulate_data.R` | 13 | 0 |

---

## 1. 主流程入口（`R/condped.R`）

### `condped()` — 三层工作流主函数

```r
condped(geno, pheno, covariates = NULL, alpha1 = 0.05, alpha2 = NULL,
        min_maf = 0.01, ld_matrix = NULL, traits = NULL,
        output_dir = "./condped_output", n_boot = 200L, verbose = TRUE)
```

| 参数 | 默认值 | 类型 / 含义 |
|---|---|---|
| `geno` | — | 数值基因型矩阵，个体 × SNP |
| `pheno` | — | 数值表型矩阵 / data.frame，个体 × 性状（首列非数值时视为 ID 列丢弃） |
| `covariates` | `NULL` | 可选协变量 data.frame；所有层分析前先对表型做残差化 |
| `alpha1` | `0.05` | Layer 1 family-wise 显著性水平（SNP-性状阈值由 `layer1_marginal_scan()` 内部 Bonferroni 计算） |
| `alpha2` | `NULL` | Layer 2 阈值；`NULL` 时取 `0.05 / (Layer-2 位点数 × 性状数)`，无 Layer-2 位点时取 0.05 |
| `min_maf` | `0.01` | Layer 1 最小 MAF |
| `ld_matrix` | `NULL` | 可选 SNP 相关矩阵 R，用于 MR 的 LD 剪枝 |
| `traits` | `NULL` | 可选性状名；默认取 pheno 列名，无列名时生成 `Trait1...Traitm` |
| `output_dir` | `"./condped_output"` | 输出目录（不存在则递归创建） |
| `n_boot` | `200L` | MR bootstrap 重复次数 |
| `verbose` | `TRUE` | 是否打印进度 |

- **输入校验**：个体数一致、全数值、无缺失、至少 2 个性状、SNP 列名自动补齐且不重复、traits 长度与唯一性检查（不满足即 `stop`）。
- **返回**：`invisible` list，同时写盘 `condped_result.rds` 和 `classification.csv`：
  - `layer1` — `layer1_marginal_scan()` 原始结果
  - `layer2` — `compute_conditional_phenotype(Y)` 结果（无 Layer-2 位点时为 `NULL`）
  - `cond_effects` — `fit_conditional_model()` 结果 data.frame（列：`trait, locus, pval_cond, sig_cond`）
  - `layer3` — `run_all_trait_pairs_mr()` 成对 MR 结果列表（无 Layer-2 位点时为空 list）
  - `classification` — 最终分类 data.frame（列：`SNPID, trait1, trait2, pattern, direction`），固定以 `include_not_detected = FALSE, require_steiger = FALSE, reject_heterogeneity = FALSE, require_robust_direction = FALSE` 调用 `classify_multitrait_patterns()`
  - `params` — 分析参数 list（`alpha1, alpha1_used, alpha2, alpha_mr, min_maf, traits, n_boot`）
  - 边界情况：Layer 1 无显著位点时提前返回，`layer2/cond_effects/layer3` 为 `NULL`、classification 为空 data.frame

### `.remove_covariates()` — 内部函数（未导出）

```r
.remove_covariates(Y, covariates)
```

| 参数 | 含义 |
|---|---|
| `Y` | 数值矩阵，个体 × 性状 |
| `covariates` | 协变量 data.frame / 矩阵，经 `model.matrix(~ .)` 构造设计矩阵（含截距） |

- **返回**：与 `Y` 同维度同列名的数值残差矩阵（逐列 `stats::lm.fit(design, y)$residuals`）。
- **功能**：从所有性状中扣除协变量效应。无缺失值处理（缺失检查由 `condped()` 前置完成）。

---

## 2. Layer 1 边际扫描（`R/gwas_marginal.R`）

### `layer1_marginal_scan()` — Layer 1 主函数

```r
layer1_marginal_scan(X, Y, Z = NULL, env = NULL, alpha = 0.05,
                     min_maf = 0.01, verbose = TRUE)
```

| 参数 | 默认值 | 类型 / 含义 |
|---|---|---|
| `X` | — | n×p 数值基因型矩阵（需列名 = SNP ID；支持 RIL -1/1、F2 -1/0/1 编码） |
| `Y` | — | n×m 数值表型矩阵（需列名 = 性状名，个体顺序与 X 一致） |
| `Z` | `NULL` | n×q 额外固定协变量矩阵（不含截距）；若同时给 `env`，环境哑变量追加其后 |
| `env` | `NULL` | 长度 n 的因子/字符向量，环境指示；自动哑变量编码（丢弃第一水平） |
| `alpha` | `0.05` | 家族显著性水平，阈值 = `alpha / (p × m)` |
| `min_maf` | `0.01` | MAF 过滤阈值，低于者扫描前移除 |
| `verbose` | `TRUE` | 是否打印进度 |

- **返回**：命名 list：
  - `class1` — 仅 1 个性状显著的 SNP ID 字符向量
  - `layer2` — ≥2 个性状显著的 SNP ID（多效性候选）
  - `snp_summary` — 宽格式 data.frame，每个显著 SNP 一行；列 `snp_id, n_marg` 及每性状 `beta_*/se_*/p_*/sig_*`，按 `n_marg` 降序
  - `scan_result` — 长格式 data.frame，每个 (SNP, 性状) 对一行（含不显著），列 `snp_id, trait, beta, se, p_value, sig`，供 Layer 3 MR 筛工具变量
  - `threshold` — 实际 Bonferroni 阈值
  - `n_tested` — 通过 MAF 过滤后扫描的 SNP 数
  - `alpha` — 传入的显著性水平
- **功能**：对每个 SNP×性状做 OLS 扫描，Bonferroni 校正后按显著性状数把位点分为 class1（n_marg=1）与 layer2 候选（n_marg≥2）。

### `validate_layer1()` — 模拟验证

```r
validate_layer1(l1, truth,
                true_classes = c("class1", "class2", "class3", "class4", "class5"))
```

| 参数 | 含义 |
|---|---|
| `l1` | `layer1_marginal_scan()` 的输出 list |
| `truth` | 模拟真值表（`generate_dataset_*()$truth`），需含 `SNP`、`class` 列 |
| `true_classes` | 视为真实功能位点的 class 标签 |

- **返回**：`invisible` list（`power`、`fdr`、`per_locus` 逐位点 data.frame），并打印格式化验证报告。

### `.ols_scan_one_trait()` — 内部函数（未导出）

```r
.ols_scan_one_trait(y, X, Z = NULL)
```

| 参数 | 含义 |
|---|---|
| `y` | 长度 n 的单性状表型向量 |
| `X` | n×p 基因型矩阵（加性编码） |
| `Z` | n×q 协变量矩阵（不含截距）；`NULL` 表示仅截距 |

- **返回**：p 行 data.frame，列 `beta, se, t_stat, p_value`（双侧 p 值）。
- **功能**：单性状向量化单 SNP-OLS 扫描；用 QR 偏出协变量（FWL 定理）+ 闭合公式，避免逐 SNP 调 `lm()`。

---

## 3. Layer 2 条件投影（`R/conditional.R`）

### `compute_conditional_phenotype()` — 全条件表型构造

```r
compute_conditional_phenotype(Y_residual, V = NULL, tol_orthogonality = 0.01)
```

| 参数 | 默认值 | 类型 / 含义 |
|---|---|---|
| `Y_residual` | — | n×m 数值残差表型矩阵（须已去均值/协变量效应，保留遗传效应） |
| `V` | `NULL` | 可选 m×m 协方差矩阵；`NULL` 时内部用 `estimate_phenotypic_covariance()` 估计 |
| `tol_orthogonality` | `0.01` | 正交性校验容差，\|Cov\| 超阈值时警告 |

- **返回**：list：
  - `Y_cond` — n×m 条件表型矩阵（核心输出，列名同输入）
  - `gamma` — 长度 m 的 list，`gamma[[i]]` 为性状 i 的投影系数向量
  - `V` — 实际使用的 m×m 协方差矩阵
  - `schur` — 长度 m 向量，各性状 Schur 补（理论条件方差）
  - `diagnostics` — data.frame（列：`trait, max_cov, var_obs, schur, var_diff`）
- **功能**：实现论文 Eq.4，对每个性状构造全条件表型 ỹ*_i = ỹ_i − ỹ_{−i}ᵀγ，保证与其余性状正交。m=1 时直接返回原始残差并提示。

### `fit_conditional_model()` — 条件效应估计

```r
fit_conditional_model(Y_cond, X_loci, alpha2 = NULL,
                      method = c("OLS", "LMM"), env = NULL, epi_pairs = NULL)
```

| 参数 | 默认值 | 类型 / 含义 |
|---|---|---|
| `Y_cond` | — | n×m（多环境时 nK×m）条件表型矩阵，来自 `compute_conditional_phenotype()$Y_cond` |
| `X_loci` | — | n×k（或 nK×k）候选位点基因型矩阵，列名为 SNP ID |
| `alpha2` | `NULL` | Bonferroni 显著性阈值；`NULL` 时自动 = `0.05 / (k × m)` |
| `method` | `"OLS"` | `"OLS"`（单环境）或 `"LMM"`（多环境，需 lme4） |
| `env` | `NULL` | 因子/字符向量，环境标签；`method = "LMM"` 时必须提供 |
| `epi_pairs` | `NULL` | list，元素为 `c(l1, l2)` 列索引对，指定加性×加性上位性位点对（作为额外固定效应） |

- **返回**：data.frame，每行一个 (性状, 位点) 对，列：`trait, locus, theta_cond`（条件效应点估计）`, se_cond, pval_cond`（OLS 用 t 检验 / LMM 用 Satterthwaite）`, sig_cond`（pval < alpha2）`, method`（`"OLS"` / `"LMM"` / `"OLS_fallback"`）。
- **细节**：LMM 公式为 `y ~ x_snp (+ 上位性项) + (1 + x_snp | env)`，REML + bobyqa；LMM 失败退回 OLS 并标记 `OLS_fallback`；无 lmerTest 时用正态近似 p 值并警告；无变异位点（var < 1e-10）跳过。

### `compute_pairwise_conditional()` — 成对条件表型

```r
compute_pairwise_conditional(Y_residual, exposure, outcome,
                             V = NULL, tol_orthogonality = 0.01)
```

| 参数 | 含义 |
|---|---|
| `Y_residual` | n×m 残差表型矩阵 |
| `exposure` | 暴露性状，列名（字符）或列索引（整数） |
| `outcome` | 结局性状，列名（字符）或列索引（整数） |
| `V` | 可选 m×m 协方差矩阵；`NULL` 时用 `cov(Y_residual)` 估计 |
| `tol_orthogonality` | 正交性校验容差 |

- **返回**：list（`y_cond` 成对条件表型向量、`gamma` 标量投影系数 Cov(out,exp)/Var(exp)、`schur` 理论条件方差、`max_cov` 正交性校验指标）。
- **功能**：只把结局对单一暴露做投影（2×2 子结构），供 Layer 3 MR 的 IV 排他性筛选；m≥3 时保留混杂经后门路径的残余信号以挡住污染工具。exposure = outcome 或暴露方差 < 1e-12 时报错。

### `print_conditional_diagnostics()` — 诊断打印

```r
print_conditional_diagnostics(result)
```

| 参数 | 含义 |
|---|---|
| `result` | `compute_conditional_phenotype()` 的返回值 |

- **返回**：无返回值（仅 `cat` 打印到控制台）：逐性状正交性（max\|Cov\|）、方差等价性（观测 vs Schur 补）及各性状投影系数 γ。

### `estimate_phenotypic_covariance()` — 内部函数（未导出）

```r
estimate_phenotypic_covariance(Y_residual, method = "pearson")
```

- **参数**：`Y_residual` n×m 残差表型；`method` 目前仅支持 `"pearson"`。
- **返回**：`list(V, Cor, n, m)` — 样本协方差矩阵、相关矩阵、样本量、性状数。
- **细节**：n ≤ m 报错；n/m < 10 警告估计不稳定。

### `compute_projection_coefficients()` — 内部函数（未导出）

```r
compute_projection_coefficients(V, i)
```

- **参数**：`V` m×m 协方差矩阵；`i` 目标性状索引（1-based）。
- **返回**：长度 m−1 的投影系数向量 γ_{i,−i} = V_{−i}⁻¹ C_{−i,i}。
- **细节**：用 `solve(A, b)` 保证数值稳定；条件数 > 1e10 时警告近奇异。

### `.project_outcome_on()` — 内部函数（未导出）

```r
.project_outcome_on(Y_residual, outcome_idx, cond_idx, V)
```

- **参数**：`Y_residual` n×m 残差表型；`outcome_idx` 结局列索引；`cond_idx` 条件性状集列索引（pairwise = {exposure}，full = 除 outcome 外全部）；`V` m×m 协方差矩阵。
- **返回**：长度 n 的数值向量 —— 结局对条件性状集 GLS 投影后的残差表型；`cond_idx` 为空时直接返回原始 outcome 列。

---

## 4. Layer 3 MR（`R/mr_v2.R`）

### 4.1 工具变量选择

#### `select_instruments_pairwise_filter()` — 按方向选择 IV

```r
select_instruments_pairwise_filter(
    marginal_effects, Y_residual, X_loci, exposure, outcome, alpha1,
    alpha_filter = 0.05, F_threshold = 10, r2_threshold = 0.1,
    ld_matrix = NULL, V = NULL,
    screen_mode = c("pairwise", "full", "none"), min_iv = 3L)
```

| 参数 | 默认值 | 类型 / 含义 |
|---|---|---|
| `marginal_effects` | — | data.frame，须含列 `TRAIT, SNPID, A, SE, P_Value` |
| `Y_residual` | — | n×m 残差表型矩阵（须列名） |
| `X_loci` | — | n×p 候选位点基因型矩阵（须 SNP 列名） |
| `exposure` / `outcome` | — | 暴露 / 结局性状名 |
| `alpha1` | — | 边际相关性阈值，无默认，须 ∈ (0,1) |
| `alpha_filter` | `0.05` | 条件结局效应过滤阈值 |
| `F_threshold` | `10` | IV 强度 F 统计量阈值 |
| `r2_threshold` | `0.1` | LD 剪枝 r² 阈值 |
| `ld_matrix` | `NULL` | 可选 SNP 相关矩阵 R |
| `V` | `NULL` | 可选表型协方差矩阵 |
| `screen_mode` | `"pairwise"` | 条件过滤口径：`"pairwise"`（推荐）/ `"full"`（仅方法比较）/ `"none"` |
| `min_iv` | `3L` | 最小 IV 数 |

- **返回**：结构化 list：`ivs`、`n_iv`、`status`（`"ok"` / `"no_relevant_instruments"` / `"no_instruments_after_filtering"` / `"insufficient_iv"`）、`exposure`/`outcome`、`screen_mode`、`audit`（逐候选 SNP 审计表：`beta/se/p_exposure, F_stat, strength_pass, p_cond_outcome, conditional_filter_pass, ld_pass, selected`）、`projection`、`thresholds`。
- **功能**：按方向选择 IV：相关性 → F 统计量 → 成对/全条件结局效应过滤 → LD 剪枝。

#### `select_instruments_pairwise()` — 向后兼容别名

与 `select_instruments_pairwise_filter()` 完全同参同返回。

### 4.2 因果估计与诊断

#### `estimate_causal_ld_gls()` — LD 感知 GLS/IVW 估计

```r
estimate_causal_ld_gls(theta_exp, theta_out, se_out, ld_matrix = NULL,
                       Y = NULL, X_iv = NULL, exp_col = 1L, out_col = 2L,
                       n_boot = 200L, seed = NULL)
```

| 参数 | 默认值 | 类型 / 含义 |
|---|---|---|
| `theta_exp` / `theta_out` / `se_out` | — | 命名数值向量：IV 的暴露效应 / 结局效应 / 结局标准误（按名字取共有 IV） |
| `ld_matrix` | `NULL` | 可选 SNP 相关矩阵 R（IV 间 LD） |
| `Y` / `X_iv` | `NULL` | 个体水平表型 / IV 基因型数据，供 bootstrap |
| `exp_col` / `out_col` | `1L` / `2L` | Y 中暴露 / 结局列索引 |
| `n_boot` | `200L` | bootstrap 次数；≤1 时仅模型 SE |
| `seed` | `NULL` | 随机种子 |

- **返回**：list（`gamma, se, se_model, z, pval, n_iv, method`（`"model_based_gls"` / `"individual_bootstrap_gls"`）`, ivs, q`（Q 诊断 list）`, weighted_median, bootstrap_values, status = "ok"`）。共有 IV < 3 时报错。

#### `estimate_causal_gls()` — 向后兼容包装

```r
estimate_causal_gls(theta_exp, theta_out, ld_matrix = NULL,
                    Y = NULL, X_iv = NULL, exp_col = 1L, out_col = 2L,
                    n_boot = 200L, se_out = NULL, seed = NULL)
```

- 返回同 `estimate_causal_ld_gls()`；未给 `se_out` 时警告并用等权近似，转调 `estimate_causal_ld_gls()`。

#### `cochran_q_diagnostic()` — 异质性诊断

```r
cochran_q_diagnostic(beta_exp, beta_out, se_out, gamma, R = NULL)
```

- **参数**：各 IV 暴露/结局效应、结局 SE、因果估计 `gamma`、可选相关矩阵 `R`。
- **返回**：list（`Q, df = k−1, pval, heterogeneous` 逻辑值）；k < 3 或有缺失时全 NA。

#### `weighted_median_mr()` — 加权中位数稳健估计

```r
weighted_median_mr(beta_exp, beta_out, se_out)
```

- **返回**：list（`gamma, n_iv, status` = `"ok"` / `"insufficient_iv"`）；比率 β_out/β_exp，权重 β_exp²/se_out²，有效 IV < 3 则不足。

#### `steiger_direction_diagnostic()` — Steiger 方向诊断

```r
steiger_direction_diagnostic(beta_exp, beta_out, var_g, var_exp, var_out,
                             consistency_threshold = 0.5)
```

- **参数**：各 IV 的暴露/结局效应、SNP 基因型方差、暴露/结局表型方差、一致性比例阈值。
- **返回**：list（`consistent` 逻辑值、`proportion_consistent`、`n_iv`、`per_iv` data.frame：`r2_exposure, r2_outcome, consistent`）。

### 4.3 双向 MR

#### `bidirectional_mr()` — 单个性状对双向 MR（实际未导出，见文末备注）

```r
bidirectional_mr(marginal_effects, traits, Y_residual, Y = NULL, X_all,
                 ld_matrix = NULL, alpha1, alpha_filter = 0.05,
                 F_threshold = 10, r2_threshold = 0.1, alpha_mr = 0.05,
                 n_boot = 200L, V = NULL,
                 screen_mode = c("pairwise", "full", "none"), min_iv = 3L,
                 steiger_threshold = 0.5, robust_tolerance = 1.0, seed = NULL)
```

- **参数**：`traits` 两个不同性状名；`X_all` 全部候选位点基因型；`steiger_threshold` Steiger 一致性阈值；`robust_tolerance` GLS 与加权中位数方向一致的相对容差；其余同 4.1/4.2。
- **返回**：list（`traits`、`AB`（A→B 方向结果）、`BA`（B→A）、`iv_AB`、`iv_BA`、`settings`）。每个方向结果含 `estimate_causal_ld_gls()` 全部字段 + `sig`（pval < alpha_mr）、`selection`、`steiger`、`robust_direction_consistent`；失败方向为标准空结构（status 如 `"insufficient_complete_effects"` / `"estimation_failed"`，附 `error_message`）。

#### `run_all_trait_pairs_mr()` — 全部性状对 MR

```r
run_all_trait_pairs_mr(marginal_effects, traits, Y_residual, Y = NULL, X_all,
                       ld_matrix = NULL, alpha1, alpha_filter = 0.05,
                       F_threshold = 10, r2_threshold = 0.1,
                       alpha_mr = NULL, adjust_mr = TRUE, n_boot = 200L,
                       V = NULL, screen_mode = c("pairwise", "full", "none"),
                       min_iv = 3L, steiger_threshold = 0.5,
                       robust_tolerance = 1.0, seed = NULL)
```

- **与 `bidirectional_mr()` 的差异**：`traits` 为 ≥2 个性状向量；`alpha_mr = NULL` 且 `adjust_mr = TRUE` 时按 Bonferroni 自动取 `0.05 / (m(m−1))`。
- **返回**：类为 `condped_pairwise_mr_list` 的命名 list，每元素为一对性状的 `bidirectional_mr()` 结果（名如 `"T1__T2"`），附属性 `alpha_mr`、`traits`。共 C(m,2) 对、m(m−1) 个方向。

### 4.4 MVMR 扩展

#### `mvmr_estimate()` — 多变量 MR（可选扩展，不在默认流程）

```r
mvmr_estimate(marginal_effects, exposures, outcome, alpha1,
              F_threshold = 10, r2_threshold = 0.1, ld_matrix = NULL,
              focal = NULL)
```

| 参数 | 含义 |
|---|---|
| `exposures` | ≥2 个共暴露性状名 |
| `focal` | 关注暴露（缺省第一个） |
| 其余 | 同 `select_instruments_pairwise_filter()` |

- **返回**：成功时 list（`status = "ok"`，focal 的 `gamma/se/z/pval`，`n_iv, ivs, gamma_all, se_all, exposures, focal, outcome`）；失败时 `status = "insufficient_iv"` / `"insufficient_complete_effects"`，`gamma/se/pval` 为 NA。
- **功能**：多变量 MR（GLS），IV 取各暴露相关性并集 + 任一暴露 F 过滤 + LD 剪枝。

### 4.5 内部辅助函数（均未导出）

| 函数 | 签名 | 功能 / 返回 |
|---|---|---|
| `.validate_marginal_effects()` | `(marginal_effects)` | 校验并标准化边际效应表（列 `TRAIT, SNPID, A, SE, P_Value`；重复记录报错） |
| `.safe_inverse()` | `(M, ridge = 1e-8)` | 对称化求逆，失败加 ridge 再试；返回逆矩阵 |
| `compute_f_stat()` | `(A, SE)` | F = (A/SE)² 向量，非法位置 NA |
| `.fast_lm_effect()` | `(x, y)` | 手写一元回归，返回 `c(beta, se, pval)`（n<4 或 x 方差过小时全 NA） |
| `.project_outcome_pairwise()` | `(Y_residual, exposure, outcome, V = NULL)` | 成对条件结局投影，返回 `list(y_cond, gamma, exposure, outcome)` |
| `.project_outcome_full()` | `(Y_residual, outcome, V = NULL)` | 全条件投影（除 outcome 外全部性状），返回同结构 list |
| `ld_prune_ivs()` | `(loci, ld_matrix, r2_threshold = 0.1, priority = NULL)` | 按优先级贪心剔除高 LD 位点，返回保留的 SNP ID |
| `.build_outcome_covariance()` | `(se_out, R = NULL, ridge = 1e-8)` | 构造 GLS 结局协方差矩阵 Ω = D R D（加 ridge） |
| `.estimate_gls_core()` | `(beta_exp, beta_out, se_out, R = NULL)` | 单参数 GLS/IVW 核心，返回 `list(gamma, se, W, Omega)` |
| `.weighted_median()` | `(values, weights)` | 加权中位数标量（无有效值时 NA） |
| `.empty_direction_result()` | `(exposure, outcome, status, selection = NULL)` | 构造失败方向结果的标准空模板 |
| `.run_mr_direction()` | `(..., seed = NULL)`（19 个参数，透传所有阈值设置） | 跑一个方向：选 IV → 取效应 → GLS → Steiger → 稳健性一致性 |

---

## 5. 模式分类（`R/classify_v2.R`）

### `classify_pair_pattern()` — 单个「位点 × 性状对」模式判定

```r
classify_pair_pattern(marg_sig_trait1, marg_sig_trait2,
                      cond_sig_trait1, cond_sig_trait2,
                      mr_results = NULL, traits = c("Trait1", "Trait2"),
                      min_iv = 3L, require_steiger = TRUE,
                      reject_heterogeneity = TRUE,
                      require_robust_direction = FALSE,
                      return_details = TRUE)
```

| 参数 | 默认值 | 类型 / 含义 |
|---|---|---|
| `marg_sig_trait1` / `marg_sig_trait2` | — | 逻辑标量，两性状的边际显著性 |
| `cond_sig_trait1` / `cond_sig_trait2` | — | 逻辑标量，两性状的条件显著性 |
| `mr_results` | `NULL` | `bidirectional_mr()` 输出，含 `$AB`（trait1→trait2）和 `$BA` 两方向 |
| `traits` | `c("Trait1", "Trait2")` | 长度 2 的字符向量，必须不同名 |
| `min_iv` | `3L` | 方向可估计所需的最小 IV 数 |
| `require_steiger` | `TRUE` | 判定 pattern3/4 前要求 Steiger 一致 |
| `reject_heterogeneity` | `TRUE` | 支持方向有显著 Cochran-Q 异质性时返回 unresolved |
| `require_robust_direction` | `FALSE` | 判定 pattern3/4 前要求 GLS 与 weighted-median 方向一致 |
| `return_details` | `TRUE` | `FALSE` 时只返回模式字符串 |

- **返回**：`return_details = TRUE` 时为单行 data.frame，列 `trait1, trait2, marginal_sig_trait1/2, conditional_sig_trait1/2, pattern, direction, confidence, unresolved_reason, mr_forward_status, mr_reverse_status`；否则为单个模式字符串。
- **模式取值**：`not_detected` / `pattern1`–`pattern5` / `unresolved`。

### `classify_multitrait_patterns()` — 批量模式标注

```r
classify_multitrait_patterns(marginal_results, conditional_results,
                             mr_by_pair, traits, alpha1, alpha2,
                             loci = NULL, include_not_detected = FALSE,
                             min_iv = 3L, require_steiger = TRUE,
                             reject_heterogeneity = TRUE,
                             require_robust_direction = FALSE)
```

| 参数 | 默认值 | 类型 / 含义 |
|---|---|---|
| `marginal_results` | — | Layer-1 结果表；支持列名 `TRAIT/trait`、`SNPID/snp_id/locus`、`P_Value/p_value/pval/p` |
| `conditional_results` | — | Layer-2 结果表；支持列名 `trait/locus/pval_cond/sig_cond`（pval_cond 与 sig_cond 至少其一） |
| `mr_by_pair` | — | `run_all_trait_pairs_mr()` 返回的命名 list，key 为 `"trait1__trait2"`（须与 traits 顺序一致，反向 key 报错） |
| `traits` | — | 性状名（决定 pair 组合与排序），至少 2 个且须同时出现在两表中 |
| `alpha1` / `alpha2` | — | Layer-1 / Layer-2 显著性阈值（无默认）；alpha2 仅在缺 `sig_cond` 列时使用 |
| `loci` | `NULL` | 可选 SNP 子集；默认取 marginal_results 中 SNP 并集 |
| `include_not_detected` | `FALSE` | 是否保留两性状均不边际显著的 pair 行 |
| `min_iv` 等 4 个 | 同上 | 透传给 `classify_pair_pattern()` |

- **返回**：长格式 data.frame，每行一个 locus × 无序性状对（在 `classify_pair_pattern()` 列基础上前置 `SNPID` 列）；无行时返回仅含 `SNPID/trait1/trait2/pattern/direction` 的空 data.frame。

### `classify_locus()` — 已弃用的两性状兼容包装

```r
classify_locus(n_marg_sig, cond_sig_by_trait, mr_results = NULL, traits = NULL)
```

| 参数 | 含义 |
|---|---|
| `n_marg_sig` | 边际显著性状计数（NA 或 ≤0 → `"not_detected"`，=1 → `"pattern1"`） |
| `cond_sig_by_trait` | 长度 2 的命名向量/list，各性状条件显著性 |
| `mr_results` | 同 `classify_pair_pattern()` |
| `traits` | 默认取 `names(cond_sig_by_trait)`；仅支持两个性状 |

- **返回**：单个模式字符串（内部以 `return_details = FALSE` 转发 `classify_pair_pattern()`，调用时发出弃用 warning）。

### 内部辅助函数（均未导出）

| 函数 | 签名 | 功能 / 返回 |
|---|---|---|
| `.normalize_marginal_results()` | `(marginal_results)` | 按别名映射标准化 Layer-1 表为 `TRAIT/SNPID/P_Value`；列名无法匹配或重复记录时 `stop` |
| `.normalize_conditional_results()` | `(conditional_results)` | 标准化 Layer-2 表为 `trait/locus/pval_cond/sig_cond_input`；缺必需列或重复记录时 `stop` |
| `.mr_direction_valid()` | `(mr_result, min_iv = 3L)` | 方向可估计性检查：`status == "ok"` 且 `pval/gamma/n_iv` 为有限标量、`n_iv >= min_iv`；返回逻辑标量 |
| `.mr_direction_significant()` | `(mr_result, min_iv = 3L)` | 在 valid 基础上再要求 `isTRUE(mr_result$sig)` |
| `.direction_diagnostic_status()` | `(mr_result, require_steiger = TRUE, reject_heterogeneity = TRUE, require_robust_direction = FALSE)` | 汇总诊断项（Steiger 不一致/不可用、IV 异质性、GLS 与 weighted-median 冲突），返回 `list(ok, reasons)` |

---

## 6. 数据模拟（`R/simulate_data.R`，13 个函数全部导出）

### 6.1 基础工具

#### `sim_genotype()` — 基因型生成

```r
sim_genotype(n, p, maf = 0.3, population = c("RIL", "F2"), seed = NULL)
```

- **参数**：`n` 样本量；`p` SNP 数；`maf` 次等位基因频率；`population` `"RIL"`（-1/+1 编码）或 `"F2"`（-1/0/+1 编码）；`seed` 随机种子。
- **返回**：n×p 数值矩阵，行名 `Ind1...`，列名 `SNP1...`。

#### `beta_from_pve()` — PVE 反算效应量

```r
beta_from_pve(pve, maf = 0.3, population = c("RIL", "F2"))
```

- **参数**：`pve` 位点解释的表型方差比例（0–1）；`population` 决定 Var(X)：RIL 为 4·maf·(1−maf)，F2 为 2·maf·(1−maf)。
- **返回**：数值标量 β = sqrt(PVE / Var(X))（Var(Y)=1 标准化下）。

### 6.2 方差标准化

#### `compute_residual_variance()`

```r
compute_residual_variance(B, Tau, var_x, rho = 0, max_iter = 200L, tol = 1e-8)
```

| 参数 | 含义 |
|---|---|
| `B` | p×m SNP 效应矩阵 |
| `Tau` | m×m 因果效应矩阵（`Tau[j,i]` = 性状 i→j 的效应） |
| `var_x` | 基因型方差 |
| `rho` | 残差相关系数 |
| `max_iter` / `tol` | 迭代上限 / 收敛阈值 |

- **返回**：m×m 残差协方差矩阵 Sigma_E。Tau=0 时解析求解、Tau≠0 时迭代调整，使 diag(Var(Y)) = 1；不收敛时 warning。

### 6.3 表型生成模型

#### `sim_phenotype_M1()` — 无因果链模型

```r
sim_phenotype_M1(X, B, Sigma_E, epi_pairs = NULL, epi_effects = NULL, seed = NULL)
```

- **参数**：`X` n×p 基因型；`B` p×m 加性效应矩阵；`Sigma_E` m×m 残差协方差；`epi_pairs` 上位性位点对 list（元素 `c(l1, l2)`）；`epi_effects` 行数=上位对数、列数=m 的效应矩阵。
- **返回**：`list(Y, G, eps)` — n×m 表型 / 遗传值 / 残差。
- **功能**：Y = X·B + eps，可选加性×加性上位性。

#### `sim_phenotype_M2()` — 有因果链模型

```r
sim_phenotype_M2(X, B, Tau, Sigma_E, epi_pairs = NULL, epi_effects = NULL, seed = NULL)
```

- **参数**：同 M1，另加 `Tau` m×m 因果效应矩阵。
- **返回**：`list(Y, G, eps, Tau, I_T_inv)`（`I_T_inv` 为因果传播算子 (I − Tau′)⁻¹）。
- **功能**：内部调用 M1 生成基础表型后乘传播算子。

#### `sim_phenotype_multienv()` — 多环境扩展（补充模拟 S5）

```r
sim_phenotype_multienv(X, B, Sigma_E, K = 3L, sigma_env = 0.3,
                       sigma_gxe = 0.1, gxe_loci = NULL, Tau = NULL, seed = NULL)
```

| 参数 | 默认值 | 含义 |
|---|---|---|
| `K` | `3L` | 环境数 |
| `sigma_env` | `0.3` | 环境主效应 SD |
| `sigma_gxe` | `0.1` | G×E 互作 SD（NULL/0 = 无互作） |
| `gxe_loci` | `NULL` | 参与 G×E 的位点索引（NULL = 所有非零效应位点） |
| `Tau` | `NULL` | 因果矩阵（NULL = M1 模型） |

- **返回**：`list(Y_long, Y_by_env, Y_mean, X_long, env_id, ind_id, G_main, h2_obs, K, n, sigma_env, sigma_gxe)` — 长格式表型 (nK)×m、n×m×K 数组、均值表型、长格式基因型等。

### 6.4 主模拟数据集

#### `generate_dataset_I()` — 无因果（tau=0），验证 Class 1/2

```r
generate_dataset_I(n, p = 1000L, class1_pve = 0.05, class2_pve = 0.05,
                   rho = 0, epi_pairs = NULL, epi_effects = NULL,
                   maf = 0.3, population = "RIL", seed = NULL)
```

- **位点布局**：SNP1–2 = Class 1（仅影响 TraitA）；SNP3–4 = Class 2（水平多效）；其余 null。
- **返回**：`list(Y, G, eps, X, truth, dataset = "I", seed, h2_obs)`；`truth` 为 data.frame（`SNP, class, beta_A, beta_B, is_IV`）。

#### `generate_dataset_II()` — 单向因果 A→B，验证 Class 1/3/4

```r
generate_dataset_II(n, p = 1000L, class1_pve = 0.05, class3_pve = 0.05,
                    class4_pve_A = 0.05, class4_pve_B = 0.05, tau = 0.3,
                    iv_pve = NULL, n_iv = 20L, rho = 0, epi_pairs = NULL,
                    epi_effects = NULL, maf = 0.3, population = "RIL", seed = NULL)
```

- **参数**：`class1_pve` Class1 对 TraitB 的 PVE；`class3_pve` Class3 对 TraitA；`class4_pve_A/B` Class4 对 A/B 的 PVE；`tau` A→B 因果效应；`iv_pve` IV 位点 PVE（NULL = `max(0.02, 12/n)`）；`n_iv` IV 位点数。
- **位点布局**：SNP1–2 Class1（仅 B）；SNP3–4 Class3（仅 A，完全中介）；SNP5–6 Class4（部分中介）；SNP7–(6+n_iv) IV 位点（仅 A，truth 标记 class3 + `is_IV = TRUE`）；其余 null。
- **返回**：`list(Y, G, eps, X, truth, dataset = "II", Tau, seed, h2_obs, iv_pve_used)`；truth 含额外列 `tau, is_IV`。

#### `generate_dataset_III()` — 双向因果 A↔B，验证 Class 5

```r
generate_dataset_III(n, p = 1000L, class5_pve = 0.05, tau_AB = 0.2,
                     tau_BA = 0.2, iv_pve = NULL, n_iv = 20L, rho = 0,
                     epi_pairs = NULL, epi_effects = NULL, maf = 0.3,
                     population = "RIL", seed = NULL)
```

- **位点布局**：SNP1–2 Class5；SNP3–(2+n_iv) IV_A（仅 A）；之后 n_iv 个 IV_B（仅 B）；其余 null。
- **返回**：`list(Y, G, eps, X, truth, dataset = "III", Tau, seed, h2_obs, iv_pve_used)`；truth 含 `tau_AB, tau_BA, is_IV` 列。

### 6.5 补充模拟数据集

#### `generate_dataset_IV()` — 三性状因果拓扑，对比 pairwise MR 与 MVMR

```r
generate_dataset_IV(n, p = 1000L,
                    topology = c("chain", "confound_obs", "confound_latent"),
                    tau_AB = 0.4, tau_BC = 0.4, tau_CA = 0.5, tau_CB = 0.5,
                    tau_UA = 0.5, tau_UB = 0.5, conf_pve = 0.03,
                    class_pve = 0.05, iv_pve = NULL, n_iv = 20L,
                    n_conf = 10L, rho = 0, maf = 0.3,
                    population = "RIL", seed = NULL)
```

- **参数**：`topology` 三种因果拓扑（`chain` 因果链 A→B→C / `confound_obs` 可观测混杂 C→A,C→B / `confound_latent` 潜在混杂 U→A,U→B）；`tau_*` 各路径因果效应；`conf_pve` 混杂工具单点 PVE（要求 `n_conf × conf_pve < 0.9`）；`class_pve` 性状特异位点 PVE；`n_conf` 混杂工具数。
- **返回**：`list(Y, G, eps, X, truth, dataset = "IV", topology, Tau_full, Tau_obs, tau_AB_true, latent, seed, h2_obs, iv_pve_used)`；Y 为 n×3（TraitA/B/C），truth 含 `role/iv_for/beta_A/B/C(/beta_U)/class/is_IV`。

#### `generate_dataset_IV_misspec()` — 四性状，验证 MVMR 错误指定混杂时的失效

```r
generate_dataset_IV_misspec(n, p = 1000L, tau_CA = 0.5, tau_CB = 0.5,
                            tau_CD = 0.5, conf_pve = 0.04, iv_pve = NULL,
                            n_iv = 20L, n_conf = 15L, rho = 0, maf = 0.3,
                            population = "RIL", seed = NULL)
```

- **参数**：`tau_CA`/`tau_CB` 真混杂路径 C→A、C→B；`tau_CD` 诱饵路径 C→D（D 是 C 下游）。
- **返回**：`list(Y, G, eps, X, truth, dataset = "IV_misspec", Tau_full, Tau_obs, tau_AB_true = 0, seed, h2_obs, iv_pve_used)`；Y 为 n×4（TraitA–D），真值 tau_AB = 0。

### 6.6 统一入口与诊断

#### `generate_condped()` — 统一数据集生成入口

```r
generate_condped(dataset = c("I", "II", "III", "IV"), ...)
```

- **参数**：`dataset` 数据集类型；`...` 透传给对应生成函数（Dataset IV 需传 `topology`）。**注意：不支持 `"IV_misspec"`**。
- **返回**：对应 `generate_dataset_*()` 的输出 list。

#### `check_dataset()` — 数据诊断

```r
check_dataset(dat)
```

- **参数**：`dat` 任意 `generate_dataset_*()` 的返回值。
- **返回**：`invisible(dat)`；副作用是打印诊断报告（性状方差、h²、功能/IV/Null 位点计数、IV 期望 F 统计量、truth 表前几行），目标 Var(Y)≈1、E[F]>12。

---

## 附：已知问题 / 备注

1. **`bidirectional_mr()` 未实际导出**：`R/mr_v2.R:835` 的 roxygen 注释有笔误（`#' #' @export`），导致该函数不在 NAMESPACE 导出列表中；它仍被导出的 `run_all_trait_pairs_mr()` 内部调用，功能正常，但用户无法直接调用 `bidirectional_mr()`。
2. **`fit_conditional_model()` 注释与实现不一致**：函数体内注释声称 OLS 分支"对所有候选位点做联合 OLS"，但实际代码是逐位点 `lm(y ~ x_snp)`（除上位性项外只含目标位点一个自变量）。
3. **`generate_condped()` 入口不支持 `"IV_misspec"`**：`dataset` 参数只允许 `"I"/"II"/"III"/"IV"`，四性状 misspec 数据集须直接调用 `generate_dataset_IV_misspec()`。
