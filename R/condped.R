# condped.R
# CondPED 主控函数：Layer 1 → Layer 2 → Layer 3 → 分类
# ============================================================================

#' condped
#' CondPED 完整分析流程
#' 
#' 执行三层分析框架：
#'   Layer 1: 边际多性状 GWAS 筛选（调用 QTXNetwork）
#'   Layer 2: 条件投影分析（条件表型 + 条件效应估计）
#'   Layer 3: 双向孟德尔随机化因果推断
#'   Classify: 六类多效性分类决策
#' 
#' @param geno 基因型矩阵或数据框，n × p。
#'   RIL: 编码为 -1/1；F2: 编码为 -1/0/1。
#'   行名 = 个体ID，列名 = SNP ID。
#' @param pheno 表型数据框，n × m。纯数值，列名 = 性状名。
#'   若第一列为字符型ID，会自动识别并排除。
#' @param qtxnetwork_path QTXNetwork 可执行文件绝对路径。Layer 1 必需。
#' @param covariates 协变量数据框，n × c。可选。若提供，会先对表型回归去残差。
#' @param alpha1 Layer 1 显著性阈值。默认 0.05。
#'   也用作 MR 工具变量筛选的边际显著性截断。
#' @param alpha2 Layer 2 Bonferroni 校正阈值。默认 NULL（自动计算 = 0.05/(k×m)）。
#' @param scan_2d 是否进行 2D 上位性扫描。默认 FALSE。
#' @param ld_matrix 位点间 LD 相关矩阵（|r| 值），p × p。可选，用于 MR 的 LD 修剪。
#'   行名/列名须与 geno 的列名一致。
#' @param traits 性状名称向量，长度 = ncol(pheno)。默认取 pheno 列名。
#' @param population 群体类型："RIL"（默认，纯合系，无显性）或 "F2"。
#' @param output_dir 结果输出目录。默认 "./condped_output"。
#' @param n_boot Bootstrap 次数，用于 MR 标准误估计。默认 1000。
#' @param verbose 是否打印进度。默认 TRUE。
#' 
#' @return 列表（invisible），包含：
#'   \describe{
#'     \item{layer1}{Layer 1 筛选结果：class1, layer2, snp_count, full_table}
#'     \item{layer2}{条件投影结果：Y_cond, gamma, V, schur, diagnostics}
#'     \item{cond_effects}{条件效应估计数据框：trait, locus, theta_cond, se_cond, pval_cond, sig_cond}
#'     \item{layer3}{双向 MR 结果：AB, BA, iv_AB, iv_BA}
#'     \item{classification}{最终分类数据框：snp_id, class, n_marg_sig, 条件效应, MR 结果}
#'     \item{params}{分析参数记录}
#'   }
#'   同时写入 output_dir：
#'   - condped_result.rds（完整结果）
#'   - classification.csv（分类表格）
#' 
#' @export
#' @examples
#' \dontrun{
#' # 1. 模拟数据（已在 Linux 生成或 R 内模拟）
#' sim <- generate_condped("class4", n = 200, p = 100, pve_A = 0.02, tau = 0.3)
#' 
#' # 2. 直接调用（最简单）
#' res <- condped(
#'   geno = sim$X,
#'   pheno = as.data.frame(sim$Y),
#'   qtxnetwork_path = "/data2/smz/QTXNetwork_4.0_MT/build/QTXNetwork",
#'   traits = c("TraitA", "TraitB"),
#'   output_dir = "./output"
#' )
#' 
#' # 3. 查看分类结果
#' print(res$classification)
#' }
condped <- function(geno,
                    pheno,
                    qtxnetwork_path = NULL,
                    covariates = NULL,
                    alpha1 = 0.05,
                    alpha2 = NULL,
                    scan_2d = FALSE,
                    ld_matrix = NULL,
                    traits = NULL,
                    population = c("RIL", "F2"),
                    output_dir = "./condped_output",
                    n_boot = 1000,
                    verbose = TRUE) {

  population <- match.arg(population)

  # --------------------------------------------------------------------------
  # 0. 输入校验与预处理
  # --------------------------------------------------------------------------
  geno <- as.matrix(geno)
  pheno <- as.data.frame(pheno)

  stopifnot(nrow(geno) == nrow(pheno))
  n <- nrow(geno)
  p <- ncol(geno)
  m <- ncol(pheno)

  # 自动识别 pheno 第一列是否为 ID（字符型或非数值）
  if (ncol(pheno) > 1 && !is.numeric(pheno[, 1])) {
    pheno <- pheno[, -1, drop = FALSE]
    m <- ncol(pheno)
  }

  if (is.null(traits)) {
    traits <- colnames(pheno)
    if (is.null(traits)) traits <- paste0("Trait", seq_len(m))
  }
  colnames(pheno) <- traits

  if (is.null(colnames(geno))) colnames(geno) <- paste0("SNP", seq_len(p))
  if (is.null(rownames(geno))) rownames(geno) <- paste0("Ind", seq_len(n))

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # 去除协变量（如果有）
  if (!is.null(covariates)) {
    Y_residual <- .remove_covariates(pheno, covariates)
  } else {
    Y_residual <- as.matrix(pheno)
  }

  if (verbose) {
    cat("=== CondPED Analysis ===\n")
    cat(sprintf("Population: %s | Samples: %d | SNPs: %d | Traits: %d\n",
                population, n, p, m))
    if (!is.null(covariates)) cat("Covariates: removed via OLS regression\n")
  }

  # --------------------------------------------------------------------------
  # 1. Layer 1: 边际多性状 GWAS（QTXNetwork）
  # --------------------------------------------------------------------------
  if (verbose) cat("\n--- Layer 1: Marginal Multi-trait GWAS ---\n")

  prefix <- file.path(output_dir, "layer1")

  # 构造标准输入格式
  geno_df <- data.frame(
    chr = rep("1", p),
    snp_id = colnames(geno),
    pos = seq_len(p),
    t(geno),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  colnames(geno_df)[4:ncol(geno_df)] <- paste0("Ind", seq_len(n))

  pheno_df <- data.frame(
    id = paste0("Ind", seq_len(n)),
    pheno,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  qtxnetwork.input.trans(geno_df, pheno_df, prefix, prefix, population = population)

  if (is.null(qtxnetwork_path)) {
    stop("qtxnetwork_path is required for Layer 1 (QTXNetwork executable).")
  }

  qtxnetwork.perform(
    qtxnetwork_path = qtxnetwork_path,
    gen_file = paste0(prefix, ".gen"),
    phe_file = paste0(prefix, ".phe"),
    pre_file = paste0(prefix, ".pre"),
    scan_2d = scan_2d
  )

  qtl <- qtxnetwork.output.trans(pheno_df, paste0(prefix, ".pre"), scan_2d = scan_2d)
  layer1 <- qtxnetwork.layer1.screen(qtl, alpha = alpha1)

  if (verbose) {
    cat(sprintf("  Class 1 (trait-specific): %d SNPs\n", length(layer1$class1)))
    cat(sprintf("  Layer 2 candidates (multi-trait): %d SNPs\n", length(layer1$layer2)))
  }

  # 无多效性候选时提前返回
  if (length(layer1$layer2) == 0) {
    if (verbose) cat("\nNo multi-trait significant SNPs. Analysis complete.\n")
    result <- data.frame(
      snp_id = layer1$class1,
      class = rep("class1", length(layer1$class1)),
      stringsAsFactors = FALSE
    )
    out <- list(layer1 = layer1, classification = result)
    saveRDS(out, file.path(output_dir, "condped_result.rds"))
    return(invisible(out))
  }

  # --------------------------------------------------------------------------
  # 2. Layer 2: 条件投影分析
  # --------------------------------------------------------------------------
  if (verbose) cat("\n--- Layer 2: Conditional Projection ---\n")

  # 2a. 估计表型协方差并构造条件表型
  cond_pheno <- compute_conditional_phenotype(Y_residual, tol_orthogonality = 0.01)

  # 2b. 提取候选位点基因型
  loci_idx <- match(layer1$layer2, colnames(geno))
  if (any(is.na(loci_idx))) {
    missing <- layer1$layer2[is.na(loci_idx)]
    warning(sprintf("Layer 2 SNPs missing in genotype: %s", paste(missing, collapse = ", ")))
    layer1$layer2 <- layer1$layer2[!is.na(loci_idx)]
    loci_idx <- loci_idx[!is.na(loci_idx)]
  }
  X_loci <- geno[, loci_idx, drop = FALSE]
  colnames(X_loci) <- layer1$layer2

  # 2c. 拟合条件模型
  if (is.null(alpha2)) alpha2 <- 0.05 / (length(layer1$layer2) * m)
  cond_effects <- fit_conditional_model(cond_pheno$Y_cond, X_loci, alpha2 = alpha2)

  if (verbose) {
    cat(sprintf("  Conditional effects tested: %d (SNP × trait pairs)\n", nrow(cond_effects)))
    cat(sprintf("  Significant after Bonferroni correction: %d\n", sum(cond_effects$sig_cond)))
  }

  # --------------------------------------------------------------------------
  # 3. Layer 3: 双向孟德尔随机化
  # --------------------------------------------------------------------------
  if (verbose) cat("\n--- Layer 3: Bidirectional MR ---\n")

  mr <- bidirectional_mr(
    marginal_effects = qtl$qtl_data,
    cond_effects = cond_effects,
    traits = traits,
    Y = as.matrix(pheno),
    X_all = geno,
    ld_matrix = ld_matrix,
    alpha1 = alpha1,
    alpha2 = alpha2,
    n_boot = n_boot
  )

  if (verbose) {
    cat(sprintf("  A→B: gamma = %.3f, SE = %.3f, p = %.3g, n_iv = %d\n",
                mr$AB$gamma, mr$AB$se, mr$AB$pval, mr$AB$n_iv))
    cat(sprintf("  B→A: gamma = %.3f, SE = %.3f, p = %.3g, n_iv = %d\n",
                mr$BA$gamma, mr$BA$se, mr$BA$pval, mr$BA$n_iv))
  }

  # --------------------------------------------------------------------------
  # 4. 六类分类决策
  # --------------------------------------------------------------------------
  if (verbose) cat("\n--- Classification ---\n")

  # 按位点整理条件显著性（命名逻辑向量）
  cond_by_locus <- lapply(layer1$layer2, function(snp) {
    sub <- cond_effects[cond_effects$locus == snp, ]
    out <- setNames(rep(FALSE, m), traits)
    if (nrow(sub) > 0) {
      sig <- setNames(sub$sig_cond, sub$trait)
      out[names(sig)] <- sig
    }
    out
  })
  names(cond_by_locus) <- layer1$layer2

  # 逐位点分类
  classes <- sapply(layer1$layer2, function(snp) {
    n_marg <- layer1$snp_count$n_sig[layer1$snp_count$snp_id == snp]
    if (length(n_marg) == 0) n_marg <- 0
    classify_locus(
      n_marg_sig = n_marg,
      cond_sig_by_trait = cond_by_locus[[snp]],
      mr_results = mr,
      traits = traits
    )
  })

  # 构造主结果表
  result <- data.frame(
    snp_id = layer1$layer2,
    class = classes,
    n_marg_sig = layer1$snp_count$n_sig[match(layer1$layer2, layer1$snp_count$snp_id)],
    stringsAsFactors = FALSE
  )

  # 合并条件效应（宽格式）
  if (nrow(cond_effects) > 0) {
    cond_wide <- stats::reshape(
      cond_effects[, c("trait", "locus", "theta_cond", "pval_cond", "sig_cond")],
      idvar = "locus", timevar = "trait", direction = "wide"
    )
    result <- merge(result, cond_wide, by.x = "snp_id", by.y = "locus", all.x = TRUE)
  }

  # 合并 MR 结果（所有位点共享同一性状对层面的因果估计）
  result$mr_AB_gamma <- mr$AB$gamma
  result$mr_AB_pval  <- mr$AB$pval
  result$mr_BA_gamma <- mr$BA$gamma
  result$mr_BA_pval  <- mr$BA$pval

  # --------------------------------------------------------------------------
  # 5. 保存与返回
  # --------------------------------------------------------------------------
  full_result <- list(
    layer1 = layer1,
    layer2 = cond_pheno,
    cond_effects = cond_effects,
    layer3 = mr,
    classification = result,
    params = list(
      alpha1 = alpha1,
      alpha2 = alpha2,
      traits = traits,
      population = population,
      n_boot = n_boot,
      scan_2d = scan_2d
    )
  )

  saveRDS(full_result, file.path(output_dir, "condped_result.rds"))
  write.csv(result, file.path(output_dir, "classification.csv"), row.names = FALSE)

  if (verbose) {
    cat(sprintf("\nResults saved to:\n  %s\n  %s\n",
                file.path(output_dir, "condped_result.rds"),
                file.path(output_dir, "classification.csv")))
    summarise_classification(result)
  }

  invisible(full_result)
}


# ============================================================================
# 内部辅助函数
# ============================================================================

#' 去除协变量（OLS 回归残差）
#' @keywords internal
.remove_covariates <- function(Y, cov) {
  Y <- as.matrix(Y)
  cov <- as.matrix(cov)
  resid <- matrix(NA, nrow(Y), ncol(Y))
  for (j in seq_len(ncol(Y))) {
    fit <- lm(Y[, j] ~ cov)
    resid[, j] <- residuals(fit)
  }
  colnames(resid) <- colnames(Y)
  resid
}


# ============================================================================
# 兼容运算符（避免依赖 rlang）
# ============================================================================
`%||%` <- function(x, y) if (is.null(x)) y else x