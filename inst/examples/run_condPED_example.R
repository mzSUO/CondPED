# =============================================================================
# CondPED Example Analysis Pipeline
# =============================================================================
# 本脚本使用 CondPED 自带的模拟数据生成器，演示完整分析流程。
# 运行前请确保：
#   1. 处于 CondPED 包根目录
#   2. 已安装/加载 CondPED（开发阶段用 devtools::load_all）
#
# Usage:
#   Rscript inst/examples/run_condPED_example.R
#   # 或在 RStudio 里:
#   # setwd("/Users/mingzhe/Desktop/CondPED")
#   # devtools::load_all(".")
#   # source("inst/examples/run_condPED_example.R")

# ------------------------------------------------------------------------------
# 0. 加载环境
# ------------------------------------------------------------------------------
if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".")
} else {
  stop("Please install devtools: install.packages('devtools')")
}

cat("
================================================================================
              CondPED (Conditional Projection-based Effect Decomposition)
                          Real Data Example
================================================================================

This example demonstrates the complete CondPED analysis pipeline using
internally simulated multi-trait data with known ground truth.

Analysis Steps:
  Step 1: Simulate Dataset II (unidirectional A -> B, tau = 0.3)
  Step 2: Layer 1 - Marginal Multi-Trait Scan
  Step 3: Layer 2 - Conditional Projection
  Step 4: Layer 3 - Bidirectional Mendelian Randomization
  Step 5: QTL Classification (5-class mechanism assignment)
  Step 6: Validation against Ground Truth

================================================================================
")

# =============================================================================
# Step 1: Simulate Multi-Trait Data (Dataset II)
# =============================================================================
cat("=== Step 1: Generating Dataset II (A -> B, tau = 0.3) ===\n\n")

set.seed(42)

dat <- generate_dataset_II(
  n            = 1000,
  p            = 1000,           
  class1_pve   = 0.01,
  class3_pve_A = 0.02,
  class4_pve_A = 0.02,
  class4_pve_B = 0.003,
  tau          = 0.3,
  iv_pve       = 0.02,
  maf          = 0.3,
  population   = "RIL",
  seed         = 42
)


cat("Simulated Data Summarys:\n")
cat("  - Individuals (n):", nrow(dat$Y), "\n")
cat("  - SNPs (p):", ncol(dat$X), "\n")
cat("  - Traits (m):", ncol(dat$Y), "\n")
cat("  - True causal effect (tau):", dat$Tau[2, 1], "\n")
cat("  - Phenotypic variance (Trait 1):", round(var(dat$Y[, 1]), 4), "\n")
cat("  - Phenotypic variance (Trait 2):", round(var(dat$Y[, 2]), 4), "\n")
cat("  - Class 1 loci (B-specific):", sum(dat$truth$class == "class1"), "\n")
cat("  - Class 3 loci (complete mediation):", sum(dat$truth$class == "class3"), "\n")
cat("  - Class 4 loci (partial mediation):", sum(dat$truth$class == "class4"), "\n")
cat("  - IV auxiliary loci:", sum(dat$truth$is_IV), "\n")
cat("  - Null loci:", sum(dat$truth$class == "null"), "\n\n")

# 提取真值供后续验证
ground_truth <- dat$truth
iv_idx <- which(dat$truth$is_IV)

# =============================================================================
# Step 2: Layer 1 - Marginal Multi-Trait Scan
# =============================================================================
cat("=== Step 2: Layer 1 - Marginal Multi-Trait Scan ===\n\n")

# --- 2.1 确定输出路径（当前工作目录下创建临时目录）---
output_dir <- file.path(getwd(), "inst", "simulation", "results", "pilot")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

prefix <- file.path(output_dir, "dataset_II")

# --- 2.2 转换数据为 QTLNetwork 格式 ---
condped_to_qtlnetwork(
  sim_result          = dat,
  geno_output_prefix  = prefix,
  pheno_output_prefix = prefix,
  population          = "RIL"
)

# 验证文件是否生成
list.files(output_dir, pattern = "\\.(gen|phe)$")

# --- 2.3 运行 QTLNetwork ---
# ⚠️ 注意：这里需要服务器上的 QTLNetwork 可执行文件路径
# 如果你在 Mac 上开发，但 QTLNetwork 在 Linux 服务器上，需要把 .gen/.phe 传到服务器跑

qtxnetwork_path <- "/data2/smz/QTXNetwork_4.0_MT/build/QTXNetwork"  # 改成你的实际路径

qtxnetwork.perform(
  qtxnetwork_path = qtxnetwork_path,
  gen_file        = paste0(prefix, ".gen"),
  phe_file        = paste0(prefix, ".phe"),
  pre_file        = paste0(prefix, ".pre"),
  scan_2d         = FALSE,   # 先只跑 1D，快
  qtx_mode        = 2        # 多性状模式
)

# --- 2.4 解析 .pre 文件 ---
# pheno_data 需要与 .phe 文件中的 trait 名称一致
pheno_for_parse <- data.frame(
  id = rownames(dat$Y),
  dat$Y,
  check.names = FALSE
)
colnames(pheno_for_parse) <- c("id", colnames(dat$Y))

qtl_out <- qtxnetwork.output.trans(
  pheno_data = pheno_for_parse,
  pre_file   = paste0(prefix, ".pre"),
  scan_2d    = FALSE
)

# 查看检出结果
cat("QTLNetwork 检出位点:\n")
print(head(qtl_out$qtl_data, 20))

# --- 2.5 Layer 1 筛选 ---
# alpha 建议用 Bonferroni 校正：0.05 / p
alpha1 <- 0.05 / ncol(dat$X)

l1 <- qtxnetwork.layer1.screen(qtl_out, alpha = alpha1)

cat("\nLayer 1 筛选结果:\n")
cat("  - Class 1 候选 (仅1个性状显著):", length(l1$class1), "\n")
cat("  - Layer 2 候选 (≥2个性状显著):", length(l1$layer2), "\n")
cat("  - 各 SNP 显著性状数:\n")
print(head(l1$snp_count, 20))

# =============================================================================
# Step 3: Layer 2 - Conditional Projection
# =============================================================================
cat("=== Step 3: Layer 2 - Conditional Projection ===\n\n")

# ---------------------------------------------------------------------------
# NOTE: 以下调用假设你的条件投影函数已实现。
# 核心逻辑：对每个候选位点，控制其他性状后检验残差关联。
# ---------------------------------------------------------------------------

# 示例：假设 layer2_project() 已存在
# l2 <- layer2_project(X = dat$X, Y = dat$Y, candidates = l1$candidates)

# 占位：模拟条件投影结果
n_cond <- sapply(l1$candidates, function(j) {
  true_class <- dat$truth$class[j]
  if (true_class == "class1") return(1)   # 仅 B 显著
  if (true_class == "class3") return(1)   # 仅 A 显著（条件后）
  if (true_class == "class4") return(2)   # A,B 均显著
  return(1)
})

l2 <- data.frame(
  SNP      = dat$truth$SNP[l1$candidates],
  n_cond   = n_cond,
  theta_A  = rnorm(length(l1$candidates), 0, 0.1),
  theta_B  = rnorm(length(l1$candidates), 0, 0.1),
  p_cond_A = runif(length(l1$candidates)),
  p_cond_B = runif(length(l1$candidates))
)

cat("Layer 2 Results (Placeholder):\n")
print(head(l2, 10))
cat("\n")

# =============================================================================
# Step 4: Layer 3 - Bidirectional Mendelian Randomization
# =============================================================================
cat("=== Step 4: Layer 3 - Bidirectional MR ===\n\n")

# ---------------------------------------------------------------------------
# NOTE: 以下调用假设你的 MR 函数已实现。
# 需要 IV 位点索引 (iv_idx) 从 truth 中提取。
# ---------------------------------------------------------------------------

# 示例：假设 layer3_mr() 已存在
# l3 <- layer3_mr(X = dat$X, Y = dat$Y, iv_idx = iv_idx, l2 = l2)

# 占位：模拟 MR 结果
l3 <- list(
  AB = list(gamma = 0.28, se = 0.05, pval = 0.001, sig = TRUE),
  BA = list(gamma = 0.02, se = 0.06, pval = 0.75,  sig = FALSE)
)

cat("MR Results (Placeholder):\n")
cat("  - A -> B: gamma =", round(l3$AB$gamma, 3),
    ", SE =", round(l3$AB$se, 3),
    ", p =", format(l3$AB$pval, digits = 3),
    ", significant =", l3$AB$sig, "\n")
cat("  - B -> A: gamma =", round(l3$BA$gamma, 3),
    ", SE =", round(l3$BA$se, 3),
    ", p =", format(l3$BA$pval, digits = 3),
    ", significant =", l3$BA$sig, "\n")
cat("  - True tau (A->B):", dat$Tau[2, 1], "\n\n")

# =============================================================================
# Step 5: QTL Classification
# =============================================================================
cat("=== Step 5: QTL Classification ===\n\n")

# ---------------------------------------------------------------------------
# NOTE: 以下调用假设你的 classify() 函数已实现。
# 分类规则基于 n_cond 和 MR 方向性。
# ---------------------------------------------------------------------------

# 示例：假设 classify() 已存在
# cls <- classify(l1 = l1, l2 = l2, l3 = l3)

# 占位：根据真值直接分类（用于演示输出格式）
classification <- data.frame(
  SNP            = dat$truth$SNP[l1$candidates],
  true_class     = dat$truth$class[l1$candidates],
  predicted_class = dat$truth$class[l1$candidates],  # 占位：假设预测完美
  n_cond         = l2$n_cond,
  stringsAsFactors = FALSE
)

cat("Classification Results:\n")
print(classification)
cat("\n")

# 混淆矩阵（真值 vs 预测）
cat("Confusion Matrix (True vs Predicted):\n")
if (requireNamespace("stats", quietly = TRUE)) {
  print(table(True = classification$true_class,
              Pred = classification$predicted_class))
}
cat("\n")

# =============================================================================
# Step 6: Validation Against Ground Truth
# =============================================================================
cat("=== Step 6: Validation Against Ground Truth ===\n\n")

# 逐位点对比
cat("Locus-by-Locus Validation:\n")
for (i in seq_len(nrow(classification))) {
  snp   <- classification$SNP[i]
  true  <- classification$true_class[i]
  pred  <- classification$predicted_class[i]
  match <- ifelse(true == pred, "✅ CORRECT", "❌ MISMATCH")
  cat(sprintf("  %s: true = %-20s pred = %-20s %s\n", snp, true, pred, match))
}
cat("\n")

# 计算准确率
accuracy <- mean(classification$true_class == classification$predicted_class)
cat(sprintf("Overall Accuracy: %.1f%%\n\n", accuracy * 100))

# =============================================================================
# Step 7: Summary & Output Structure
# =============================================================================
cat("=== Step 7: Summary ===\n\n")

cat("Expected Output Structure from condped():\n\n")
cat("results/\n")
cat("├── layer1/\n")
cat("│   ├── marginal_pvals.rds      # Layer 1 p-values (p x m)\n")
cat("│   └── candidates.rds          # Candidate locus indices\n")
cat("├── layer2/\n")
cat("│   ├── conditional_effects.csv # theta_cond per locus-trait\n")
cat("│   └── n_cond.csv              # Number of significant conditional traits\n")
cat("├── layer3/\n")
cat("│   ├── mr_AB.rds               # MR A->B results (gamma, SE, p)\n")
cat("│   └── mr_BA.rds               # MR B->A results\n")
cat("├── classification/\n")
cat("│   └── qtl_classes.csv         # Final 5-class assignment\n")
cat("└── figures/\n")
cat("    ├── figure1_power_fdr.pdf\n")
cat("    ├── figure2_confusion.pdf\n")
cat("    └── figure3_tau_distribution.pdf\n\n")

# =============================================================================
# Full Pipeline Call (for reference)
# =============================================================================
cat("================================================================================
              Full condped() Call (when all functions are ready)
================================================================================

# One-line analysis after functions are implemented:
res <- condped(
  X         = dat$X,
  Y         = dat$Y,
  iv_idx    = which(dat$truth$is_IV),
  alpha1    = 0.05,
  alpha2    = 0.05 / (ncol(dat$X) * ncol(dat$Y)),
  alpha_mr  = 0.05
)

# Access results:
# res$layer1$candidates      # Layer 1 significant loci
# res$layer2$theta_cond      # Conditional effect estimates
# res$layer3$AB$gamma        # Causal effect A -> B
# res$classification         # Data frame with predicted classes

================================================================================
")

cat("Example completed successfully!\n")
cat("Next step: Replace placeholder sections (Steps 2-5) with actual function calls.\n")