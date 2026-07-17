# ==============================================================================
# condped.R — CondPED minimal main workflow
# Layer 1: marginal OLS scan
# Layer 2: conditional projection
# Layer 3: pairwise bidirectional MR
# Output : locus-trait-pair association patterns
# ==============================================================================

#' Run the CondPED workflow
#'
#' @param geno Numeric genotype matrix, individuals by SNPs.
#' @param pheno Numeric phenotype matrix or data frame, individuals by traits.
#' @param covariates Optional covariate data frame. Phenotypes are residualized
#'   before all three layers.
#' @param alpha1 Family-wise alpha level for Layer 1. The actual SNP-trait
#'   threshold is calculated inside `layer1_marginal_scan()`.
#' @param alpha2 Layer-2 threshold. If NULL, uses
#'   0.05 / (number of Layer-2 loci * number of traits).
#' @param min_maf Minimum MAF used by Layer 1.
#' @param ld_matrix Optional SNP correlation matrix R for MR LD pruning.
#' @param traits Optional trait names. Defaults to phenotype column names.
#' @param output_dir Output directory.
#' @param n_boot Number of MR bootstrap replicates.
#' @param verbose Whether to print progress.
#'
#' @return Invisibly returns a list containing Layer 1 results, conditional
#' results, pairwise MR results, classifications, and analysis parameters.
#'
#' @export
condped <- function(geno,
                    pheno,
                    covariates = NULL,
                    alpha1 = 0.05,
                    alpha2 = NULL,
                    min_maf = 0.01,
                    ld_matrix = NULL,
                    traits = NULL,
                    output_dir = "./condped_output",
                    n_boot = 200L,
                    verbose = TRUE) {
  # --------------------------------------------------------------------------
  # 0. Input preparation
  # --------------------------------------------------------------------------
  geno <- as.matrix(geno)
  pheno <- as.data.frame(pheno)

  if (ncol(pheno) > 1L && !is.numeric(pheno[[1L]])) {
    pheno <- pheno[, -1L, drop = FALSE]
  }

  if (nrow(geno) != nrow(pheno)) {
    stop("geno 与 pheno 的个体数不一致。", call. = FALSE)
  }
  if (!is.numeric(geno) ||
    !all(vapply(pheno, is.numeric, logical(1L)))) {
    stop("geno 和 pheno 必须为数值型。", call. = FALSE)
  }
  if (anyNA(geno) || anyNA(pheno)) {
    stop("geno 和 pheno 不能包含缺失值。", call. = FALSE)
  }

  n <- nrow(geno)
  p <- ncol(geno)
  m <- ncol(pheno)

  if (m < 2L) {
    stop("CondPED 至少需要两个性状。", call. = FALSE)
  }

  if (is.null(colnames(geno))) {
    colnames(geno) <- paste0("SNP", seq_len(p))
  }
  if (anyDuplicated(colnames(geno))) {
    stop("geno 的 SNP 列名不能重复。", call. = FALSE)
  }

  if (is.null(traits)) {
    traits <- colnames(pheno)
    if (is.null(traits)) {
      traits <- paste0("Trait", seq_len(m))
    }
  }
  traits <- as.character(traits)

  if (length(traits) != m || anyDuplicated(traits)) {
    stop("traits 必须与 pheno 列数一致且不能重复。", call. = FALSE)
  }

  Y <- as.matrix(pheno)
  colnames(Y) <- traits

  if (!is.null(covariates)) {
    if (nrow(covariates) != n) {
      stop("covariates 与 geno/pheno 的个体数不一致。", call. = FALSE)
    }
    if (anyNA(covariates)) {
      stop("covariates 不能包含缺失值。", call. = FALSE)
    }
    Y <- .remove_covariates(Y, covariates)
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  if (verbose) {
    cat("=== CondPED ===\n")
    cat(sprintf("Samples: %d | SNPs: %d | Traits: %d\n", n, p, m))
  }

  # --------------------------------------------------------------------------
  # 1. Layer 1: marginal OLS scan
  # --------------------------------------------------------------------------
  if (verbose) cat("\n--- Layer 1: Marginal scan ---\n")

  layer1 <- layer1_marginal_scan(
    X = geno,
    Y = Y,
    alpha = alpha1,
    min_maf = min_maf,
    verbose = verbose
  )

  marginal_effects <- data.frame(
    TRAIT = layer1$scan_result$trait,
    SNPID = layer1$scan_result$snp_id,
    A = layer1$scan_result$beta,
    SE = layer1$scan_result$se,
    P_Value = layer1$scan_result$p_value,
    stringsAsFactors = FALSE
  )

  detected_loci <- unique(c(layer1$class1, layer1$layer2))
  alpha1_used <- layer1$threshold

  if (verbose) {
    cat(sprintf("  Trait-specific loci: %d\n", length(layer1$class1)))
    cat(sprintf("  Layer-2 candidates: %d\n", length(layer1$layer2)))
  }

  # No detected locus
  if (length(detected_loci) == 0L) {
    classification <- data.frame(
      SNPID = character(),
      trait1 = character(),
      trait2 = character(),
      pattern = character(),
      direction = character(),
      stringsAsFactors = FALSE
    )

    result <- list(
      layer1 = layer1,
      layer2 = NULL,
      cond_effects = NULL,
      layer3 = NULL,
      classification = classification,
      params = list(
        alpha1 = alpha1,
        alpha1_used = alpha1_used,
        alpha2 = alpha2,
        min_maf = min_maf,
        traits = traits,
        n_boot = n_boot
      )
    )

    saveRDS(result, file.path(output_dir, "condped_result.rds"))
    utils::write.csv(
      classification,
      file.path(output_dir, "classification.csv"),
      row.names = FALSE
    )

    if (verbose) cat("\nNo significant loci.\n")
    return(invisible(result))
  }

  # --------------------------------------------------------------------------
  # 2. Layer 2: conditional projection
  # --------------------------------------------------------------------------
  cond_pheno <- NULL
  cond_effects <- data.frame(
    trait = character(),
    locus = character(),
    pval_cond = numeric(),
    sig_cond = logical(),
    stringsAsFactors = FALSE
  )

  if (length(layer1$layer2) > 0L) {
    if (verbose) cat("\n--- Layer 2: Conditional projection ---\n")

    cond_pheno <- compute_conditional_phenotype(Y)

    X_loci <- geno[, layer1$layer2, drop = FALSE]

    if (is.null(alpha2)) {
      alpha2 <- 0.05 / (length(layer1$layer2) * m)
    }

    cond_effects <- fit_conditional_model(
      Y_cond = cond_pheno$Y_cond,
      X_loci = X_loci,
      alpha2 = alpha2
    )
  } else if (is.null(alpha2)) {
    alpha2 <- 0.05
  }

  # Add non-tested loci as conditional non-significant rows.
  # This is only needed so Pattern 1 loci can enter the common classifier.
  cond_all <- expand.grid(
    trait = traits,
    locus = detected_loci,
    stringsAsFactors = FALSE
  )
  cond_all$pval_cond <- NA_real_
  cond_all$sig_cond <- FALSE

  if (nrow(cond_effects) > 0L) {
    key_all <- paste(cond_all$trait, cond_all$locus, sep = "\r")
    key_eff <- paste(cond_effects$trait, cond_effects$locus, sep = "\r")
    pos <- match(key_eff, key_all)

    cond_all$pval_cond[pos] <- cond_effects$pval_cond
    cond_all$sig_cond[pos] <- cond_effects$sig_cond
  }

  # --------------------------------------------------------------------------
  # 3. Layer 3: pairwise bidirectional MR
  # --------------------------------------------------------------------------
  mr <- list()

  if (length(layer1$layer2) > 0L) {
    if (verbose) cat("\n--- Layer 3: Pairwise bidirectional MR ---\n")

    mr <- run_all_trait_pairs_mr(
      marginal_effects = marginal_effects,
      traits = traits,
      Y_residual = Y,
      Y = Y,
      X_all = geno,
      ld_matrix = ld_matrix,
      alpha1 = alpha1_used,
      n_boot = n_boot
    )
  }

  # --------------------------------------------------------------------------
  # 4. Pair-level classification
  # --------------------------------------------------------------------------
  if (verbose) cat("\n--- Classification ---\n")

  classification <- classify_multitrait_patterns(
    marginal_results = marginal_effects,
    conditional_results = cond_all,
    mr_by_pair = mr,
    traits = traits,
    alpha1 = alpha1_used,
    alpha2 = alpha2,
    loci = detected_loci,
    include_not_detected = FALSE,
    require_steiger = FALSE,
    reject_heterogeneity = FALSE,
    require_robust_direction = FALSE
  )

  # --------------------------------------------------------------------------
  # 5. Save and return
  # --------------------------------------------------------------------------
  result <- list(
    layer1 = layer1,
    layer2 = cond_pheno,
    cond_effects = cond_effects,
    layer3 = mr,
    classification = classification,
    params = list(
      alpha1 = alpha1,
      alpha1_used = alpha1_used,
      alpha2 = alpha2,
      alpha_mr = attr(mr, "alpha_mr"),
      min_maf = min_maf,
      traits = traits,
      n_boot = n_boot
    )
  )

  saveRDS(result, file.path(output_dir, "condped_result.rds"))
  utils::write.csv(
    classification,
    file.path(output_dir, "classification.csv"),
    row.names = FALSE
  )

  if (verbose) {
    print(table(classification$pattern, useNA = "ifany"))
    cat("\nResults saved in: ", output_dir, "\n", sep = "")
  }

  invisible(result)
}


# Remove covariate effects from all traits
.remove_covariates <- function(Y, covariates) {
  design <- stats::model.matrix(~., data = as.data.frame(covariates))

  residuals <- apply(
    Y,
    2L,
    function(y) stats::lm.fit(design, y)$residuals
  )

  residuals <- as.matrix(residuals)
  colnames(residuals) <- colnames(Y)
  residuals
}
