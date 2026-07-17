# =============================================================================
# CondPED Example Analysis Pipeline (Server Version)
# =============================================================================
# Environment : Linux server
# R package   : /data2/smz/CondPED
#
# Usage:
#   cd /data2/smz/CondPED
#   Rscript inst/examples/run_condPED_example.R
#
# Analysis Steps:
#   Step 1: Simulate Dataset II (unidirectional A -> B, tau = 0.3)
#   Step 2: Layer 1 - Marginal Multi-Trait Scan (OLS + Bonferroni)
#   Step 3: Layer 2 - Conditional Projection
#   Step 4: Layer 3 - Bidirectional Mendelian Randomization
#   Step 5: QTL Classification (5-class mechanism assignment)
#   Step 6: Validation against Ground Truth
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Load environment
# -----------------------------------------------------------------------------
pkg_root <- "/data2/smz/CondPED"
if (getwd() != pkg_root) {
  setwd(pkg_root)
  message("Working directory set to: ", pkg_root)
}

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".")
} else {
  source("R/simulate_data.R")
  source("R/gwas_marginal.R")
  source("R/conditional.R")
  source("R/mr.R")
  source("R/classify.R")
}

cat("
================================================================================
         CondPED (Conditional Projection-based Effect Decomposition)
                         Server Example
================================================================================

Environment:
  R package root : /data2/smz/CondPED

Analysis Steps:
  Step 1 : Simulate Dataset II  (A -> B, tau = 0.3)
  Step 2 : Layer 1 - Marginal OLS scan with Bonferroni correction
  Step 3 : Layer 2 - Conditional projection
  Step 4 : Layer 3 - Bidirectional Mendelian Randomization
  Step 5 : QTL classification (Class 1 - 5)
  Step 6 : Validation against ground truth

================================================================================
")


# =============================================================================
# Step 1: Simulate Dataset II
# =============================================================================
cat("=== Step 1: Generating Dataset II (A -> B, tau = 0.3) ===\n\n")

dat <- generate_dataset_II(
  n = 1000,
  p = 1000,
  class1_pve = 0.05,
  class3_pve = 0.05,
  class4_pve_A = 0.05,
  class4_pve_B = 0.05,
  tau = 0.3,
  seed = 42
)

cat("Simulated Data Summary:\n")
cat("  - Individuals (n)       :", nrow(dat$Y), "\n")
cat("  - SNPs (p)              :", ncol(dat$X), "\n")
cat("  - Traits (m)            :", ncol(dat$Y), "\n")
cat("  - True causal effect tau:", dat$Tau[2, 1], "\n")
cat("  - Var(Trait1)           :", round(var(dat$Y[, 1]), 4), "\n")
cat("  - Var(Trait2)           :", round(var(dat$Y[, 2]), 4), "\n")
cat("  - Observed h2           :", round(dat$h2_obs, 4), "\n")
cat("  - Class 1 loci          :", sum(dat$truth$class == "class1"), "\n")
cat("  - Class 3 loci          :", sum(dat$truth$class == "class3"), "\n")
cat("  - Class 4 loci          :", sum(dat$truth$class == "class4"), "\n")
cat("  - IV loci               :", sum(dat$truth$is_IV), "\n")
cat("  - Null loci             :", sum(dat$truth$class == "null"), "\n\n")

cat("True effect sizes (loci 1-16):\n")
print(dat$truth[1:16, c("SNP", "class", "beta_A", "beta_B")])
cat("\n")


# =============================================================================
# Step 2: Layer 1 - Marginal Multi-Trait Scan (OLS + Bonferroni)
# =============================================================================
cat("=== Step 2: Layer 1 - Marginal OLS Scan ===\n\n")

# Statistical model (Eq. 1, additive fixed effects, single environment):
#   y_ij = mu_i + a_il * x_jl + eps_ij,   eps_ij ~ N(0, sigma^2_i)
#
# Bonferroni threshold: alpha / (p * m)
# Loci with n_marg = 1  -> Class 1 (trait-specific), no further testing
# Loci with n_marg >= 2 -> pleiotropy candidates, forwarded to Layer 2

alpha1 <- 0.05

l1 <- layer1_marginal_scan(
  X       = dat$X,
  Y       = dat$Y,
  alpha   = alpha1,
  verbose = TRUE
)

cat("\nLayer 1 Screening Results:\n")
cat("  - Bonferroni threshold  :", formatC(l1$threshold, format = "e", digits = 2), "\n")
cat("  - SNPs tested           :", l1$n_tested, "\n")
cat("  - Total candidates      :", nrow(l1$snp_summary), "\n")
cat("  - Class 1 candidates    :", length(l1$class1), "\n")
cat("  - Layer 2 candidates    :", length(l1$layer2), "\n")

cat("\nSignificant loci summary:\n")
print(l1$snp_summary)

# Validate against ground truth
cat("\n")
val_l1 <- validate_layer1(
  l1,
  dat$truth,
  true_classes = c("class1", "class3", "class4")
)


# =============================================================================
# Step 3: Layer 2 - Conditional Projection
# =============================================================================
cat("=== Step 3: Layer 2 - Conditional Projection ===\n\n")

# Construct conditional phenotypes:
#   y*_ijk = y_ijk - y_{-i,jk}^T * gamma_{i,-i}
# where gamma_{i,-i} = V_{-i}^{-1} C_{-i,i} (projection coefficients)

Y_residual <- scale(dat$Y, center = TRUE, scale = FALSE)

cond_pheno <- compute_conditional_phenotype(Y_residual)

cat("Conditional phenotype diagnostics:\n")
print(cond_pheno$diagnostics)
cat("\n")

# Fit conditional model for all significant loci (Class 1 + Layer 2 candidates)
all_sig_loci <- unique(c(l1$class1, l1$layer2))
X_layer2 <- dat$X[, all_sig_loci, drop = FALSE]

# Bonferroni correction for conditional tests: alpha / (|L| * m)
alpha2 <- 0.05 / (length(all_sig_loci) * ncol(dat$Y))

cond_eff <- fit_conditional_model(
  Y_cond = cond_pheno$Y_cond,
  X_loci = X_layer2,
  alpha2 = alpha2,
  method = "OLS" # single environment simulation; switch to "LMM" for
  # multi-environment real data (requires env argument)
)

n_cond_by_snp <- tapply(cond_eff$sig_cond, cond_eff$locus, sum)
cat("Conditional effects:\n")
print(cond_eff[, c(
  "trait", "locus", "theta_cond", "se_cond",
  "pval_cond", "sig_cond", "method"
)])


# =============================================================================
# Step 4: Layer 3 - Bidirectional Mendelian Randomization
# =============================================================================
cat("=== Step 4: Layer 3 - Bidirectional MR ===\n\n")

# IV selection criteria (Section 2.3.1):
#   (1) Relevance    : P_marg < l1$threshold  (Bonferroni-corrected)
#   (2) Exclusion    : P_cond > 0.05  (no significant direct effect on outcome)
#   (3) LD pruning   : |r_ll'| < 0.1
#   (4) Instrument   : F-statistic > 10
# 注意：IV 排他性筛选现在由 bidirectional_mr 内部按 MR 方向现算成对条件投影
# （outcome | exposure），不再消费 Step 3 的全条件 cond_eff。
# cond_eff 仅用于 Step 5 的 5 类分类（n_cond 指纹）。

marginal_df <- data.frame(
  TRAIT = l1$scan_result$trait,
  SNPID = l1$scan_result$snp_id,
  A = l1$scan_result$beta,
  SE = l1$scan_result$se,
  P_Value = l1$scan_result$p_value,
  stringsAsFactors = FALSE
)

traits <- colnames(dat$Y)

mr_res <- bidirectional_mr(
  marginal_effects = marginal_df,
  traits           = traits,
  Y_residual       = Y_residual,
  Y                = dat$Y,
  X_all            = dat$X,
  ld_matrix        = NULL,
  alpha1           = l1$threshold,
  alpha2           = 0.05,
  F_threshold      = 10,
  r2_threshold     = 0.1,
  alpha_mr         = 0.05,
  n_boot           = 200
)
cat("MR Results:\n")
cat(sprintf(
  "  A -> B : gamma = %6.4f  SE = %6.4f  p = %s  sig = %s  n_IV = %d\n",
  mr_res$AB$gamma, mr_res$AB$se,
  format(mr_res$AB$pval, digits = 3),
  mr_res$AB$sig, mr_res$AB$n_iv
))
cat(sprintf(
  "  B -> A : gamma = %6.4f  SE = %6.4f  p = %s  sig = %s  n_IV = %d\n",
  mr_res$BA$gamma, mr_res$BA$se,
  format(mr_res$BA$pval, digits = 3),
  mr_res$BA$sig, mr_res$BA$n_iv
))
cat("  True tau (A->B):", dat$Tau[2, 1], "\n\n")


# =============================================================================
# Step 5: QTL Classification
# =============================================================================
cat("=== Step 5: QTL Classification ===\n\n")

classification_l2 <- vector("list", length(l1$layer2))

for (i in seq_along(l1$layer2)) {
  snp <- l1$layer2[i]
  n_marg <- l1$snp_summary$n_marg[l1$snp_summary$snp_id == snp]
  cond_sub <- cond_eff[cond_eff$locus == snp, ]
  cond_sig <- setNames(cond_sub$sig_cond, cond_sub$trait)

  cls <- classify_locus(n_marg, cond_sig, mr_res, traits)
  true_cls <- dat$truth$class[dat$truth$SNP == snp]

  classification_l2[[i]] <- data.frame(
    SNP = snp,
    true_class = true_cls,
    predicted = cls,
    n_marg = n_marg,
    n_cond = sum(cond_sig),
    route = "layer2",
    stringsAsFactors = FALSE
  )
}

classification_l1 <- vector("list", length(l1$class1))

for (i in seq_along(l1$class1)) {
  snp <- l1$class1[i]
  cond_sub <- cond_eff[cond_eff$locus == snp, ]
  cond_sig <- setNames(cond_sub$sig_cond, cond_sub$trait)

  cls <- classify_locus(1L, cond_sig, mr_res, traits)
  true_cls <- dat$truth$class[dat$truth$SNP == snp]

  classification_l1[[i]] <- data.frame(
    SNP = snp,
    true_class = true_cls,
    predicted = cls,
    n_marg = 1L,
    n_cond = sum(cond_sig, na.rm = TRUE),
    route = "layer1_to_mr",
    stringsAsFactors = FALSE
  )
}

classification <- rbind(
  do.call(rbind, classification_l2),
  do.call(rbind, classification_l1)
)

snp_order <- order(as.integer(gsub("SNP", "", classification$SNP)))
classification <- classification[snp_order, ]
rownames(classification) <- NULL

cat("Classification Results:\n")
print(classification)

cat("\nConfusion Matrix (True vs Predicted):\n")
print(table(True = classification$true_class, Pred = classification$predicted))
cat("  Note: IV_A rows = Layer 1 false positives, not framework misclassifications.\n")

# Accuracy at three levels
func_cls <- classification[
  classification$true_class %in% c("class1", "class2", "class3", "class4", "class5"),
]
acc_func <- mean(func_cls$true_class == func_cls$predicted, na.rm = TRUE)

cls_l2 <- do.call(rbind, classification_l2)
acc_l2 <- mean(cls_l2$true_class == cls_l2$predicted, na.rm = TRUE)

acc_all <- mean(classification$true_class == classification$predicted, na.rm = TRUE)

cat(sprintf("\nFunctional loci accuracy   : %.1f%%\n", acc_func * 100))
cat(sprintf("Layer 2 accuracy           : %.1f%%\n", acc_l2 * 100))
cat(sprintf("Overall accuracy           : %.1f%%  (includes Layer 1 FP)\n\n", acc_all * 100))


# =============================================================================
# Step 6: Validation Summary
# =============================================================================
cat("=== Step 6: Validation Summary ===\n\n")

cat("Layer 1 Performance:\n")
cat(sprintf("  Power : %.1f%%\n", val_l1$power * 100))
cat(sprintf("  FDR   : %.1f%%\n", val_l1$fdr * 100))

cat("\nLocus-by-Locus Classification:\n")
for (i in seq_len(nrow(classification))) {
  snp <- classification$SNP[i]
  true <- classification$true_class[i]
  pred <- classification$predicted[i]
  route <- classification$route[i]
  status <- ifelse(true == pred, "CORRECT", "MISMATCH")
  cat(sprintf(
    "  %-6s  true = %-12s  pred = %-12s  route = %-14s  [%s]\n",
    snp, true, pred, route, status
  ))
}
cat("\n")


# =============================================================================
# Save Results
# =============================================================================
cat("=== Saving Results ===\n\n")

output_dir <- "/data2/smz/CondPED/inst/results/files"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

save(
  dat, l1, cond_pheno, cond_eff, mr_res, classification, val_l1,
  file = file.path(output_dir, "condped_example_results.RData")
)

write.csv(classification,
  file = file.path(output_dir, "classification.csv"),
  row.names = FALSE
)

write.csv(l1$snp_summary,
  file = file.path(output_dir, "layer1_summary.csv"),
  row.names = FALSE
)

cat("Results saved to:", output_dir, "\n")
cat("  - condped_example_results.RData\n")
cat("  - classification.csv\n")
cat("  - layer1_summary.csv\n")
cat("\nExample pipeline completed.\n")
