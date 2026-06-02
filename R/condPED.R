#' condPED: Conditional Projection for Dissecting Pleiotropic Effects
#'
#' Main user interface function that runs the complete CondPED analysis pipeline.
#'
#' @param pheno_data Phenotype data frame (with GID and trait columns)
#' @param geno_data Genotype matrix (n_ind x n_snp) or list with G and snp_info
#' @param qtxnetwork_path Path to QTLNetwork executable (if method = "qtlnetwork")
#' @param method GWAS method: "qtlnetwork" (default)
#' @param output_dir Output directory for results
#' @param scan_2d Whether to perform 2D scanning for epistasis
#' @param run_classification Whether to run QTL classification
#' @param run_causal Whether to run MR causal inference
#' @param alpha_marginal Marginal significance threshold (used in classification, default: 0.05/n_snp)
#' @param alpha_forward Forward conditional significance threshold (default: 0.05/(n_snp*n_trait))
#' @param alpha_reverse Reverse conditional significance threshold (default: 0.05)
#' @param ... Additional arguments
#'
#' @return List with:
#'   \item{gwas_results}{GWAS results from all steps}
#'   \item{classification}{QTL classification results}
#'   \item{causal_results}{Causal inference results (if run_causal = TRUE)}
#'
#' @export
#' @examples
#' \dontrun{
#' # Run complete CondPED analysis
#' results <- condPED(pheno_data, geno_data, qtxnetwork_path = "/path/to/QTLNetwork")
#' }
#'
condPED <- function(
  pheno_data,
  geno_data,
  qtxnetwork_path = NULL,
  method = "qtlnetwork",
  output_dir = "./condPED_output",
  scan_2d = FALSE,
  run_classification = TRUE,
  run_causal = TRUE,
  alpha_marginal = NULL,
  alpha_forward = NULL,
  alpha_reverse = 0.05,
  ...
) {

  # ==========================================================================
  # Step 0: Prepare data
  # ==========================================================================
  cat("=== CondPED Analysis ===\n\n")

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Process genotype data
  if (is.list(geno_data) && !is.null(geno_data$G)) {
    sim_data <- geno_data
  } else if (is.matrix(geno_data) || is.data.frame(geno_data)) {
    # Convert to standard format
    G <- as.matrix(geno_data)
    n_ind <- nrow(G)
    n_snp <- ncol(G)

    snp_info <- data.frame(
      CHR = rep(1, n_snp),
      SNP = colnames(G),
      BP = 1:n_snp
    )

    sim_data <- list(
      G = scale(G),
      G_raw = round(G),
      Y = as.matrix(pheno_data[, -1]),
      snp_info = snp_info
    )
    rownames(sim_data$Y) <- pheno_data$GID
    colnames(sim_data$Y) <- colnames(pheno_data)[-1]
  } else {
    stop("Invalid geno_data format")
  }

  # ==========================================================================
  # Step 1: Marginal GWAS
  # ==========================================================================
  cat("Step 1: Marginal GWAS\n")

  if (method == "qtlnetwork") {
    if (is.null(qtxnetwork_path)) {
      stop("qtxnetwork_path is required for method = 'qtlnetwork'")
    }

    gwas_results <- marginal_gwas(
      sim_data = sim_data,
      qtxnetwork_path = qtxnetwork_path,
      output_dir = output_dir,
      scan_2d = scan_2d
    )
  } else {
    stop("Only 'qtlnetwork' method is currently supported")
  }

  # ==========================================================================
  # Step 2: Forward Conditional GWAS
  # ==========================================================================
  cat("Step 2: Forward Conditional GWAS\n")

  candidate_snps <- unique(gwas_results$gwas_results$SNP)

  forward_results <- list()
  for (trait in colnames(sim_data$Y)) {
    forward_results[[trait]] <- forward_conditional_gwas(
      sim_data = sim_data,
      target_trait = trait,
      candidate_snps = candidate_snps
    )
  }

  # ==========================================================================
  # Step 3: Reverse Conditional GWAS
  # ==========================================================================
  cat("Step 3: Reverse Conditional GWAS\n")

  n_trait <- ncol(sim_data$Y)
  if (is.null(alpha_forward)) alpha_forward <- 0.05 / (length(candidate_snps) * n_trait)

  reverse_results <- list()

  for (target_trait in colnames(sim_data$Y)) {
    sig_snps <- forward_results[[target_trait]]$SNP[
      forward_results[[target_trait]]$P_Value < alpha_forward
    ]

    for (test_trait in colnames(sim_data$Y)) {
      if (test_trait != target_trait && length(sig_snps) > 0) {
        key <- paste0(test_trait, "_excl_", target_trait)
        reverse_results[[key]] <- reverse_conditional_gwas(
          sim_data = sim_data,
          target_trait = target_trait,
          test_trait = test_trait,
          candidate_snps = sig_snps
        )
      }
    }
  }

  # Combine results
  full_gwas_results <- list(
    step1_marginal = gwas_results,
    step2_forward = forward_results,
    step3_reverse = reverse_results,
    candidate_snps = candidate_snps
  )

  # ==========================================================================
  # Step 4: QTL Classification
  # ==========================================================================
  classification_results <- NULL
  if (run_classification) {
    cat("Step 4: QTL Classification\n")
    classification_results <- classify_qtl(
      full_gwas_results,
      alpha1 = alpha_marginal,
      alpha2 = alpha_forward,
      alpha3 = alpha_reverse
    )
  }

  # ==========================================================================
  # Step 5: Causal Inference (MR)
  # ==========================================================================
  causal_results <- NULL
  if (run_causal && n_trait >= 2) {
    cat("Step 5: Causal Inference\n")

    # Run bidirectional causal tests for each trait pair
    causal_results <- list()

    traits <- colnames(sim_data$Y)
    for (i in 1:(length(traits) - 1)) {
      for (j in (i + 1):length(traits)) {
        key <- paste0(traits[i], "_vs_", traits[j])
        causal_results[[key]] <- bidirectional_causal_test(
          full_gwas_results,
          traits[i],
          traits[j],
          G = sim_data$G
        )
      }
    }
  }

  # ==========================================================================
  # Return results
  # ==========================================================================
  cat("\n=== Analysis Complete ===\n")

  return(invisible(list(
    gwas_results = full_gwas_results,
    classification = classification_results,
    causal_results = causal_results,
    sim_data = sim_data,
    parameters = list(
      method = method,
      scan_2d = scan_2d,
      output_dir = output_dir,
      alpha_marginal = alpha_marginal,
      alpha_forward = alpha_forward,
      alpha_reverse = alpha_reverse
    )
  )))
}


#' Print CondPED results
#'
#' @param x CondPED results object
#' @param ... Additional arguments
#'
#' @export
#'
print.condped <- function(x, ...) {

  cat("=== CondPED Results ===\n\n")

  if (!is.null(x$classification)) {
    cat("QTL Classification:\n")
    print(table(x$classification$classification))
    cat("\n")
  }

  if (!is.null(x$causal_results)) {
    cat("Causal Inference Results:\n")
    for (key in names(x$causal_results)) {
      cat("\n", key, ":\n")
      res <- x$causal_results[[key]]
      cat("  ", x$traits[1], "->", x$traits[2], ": ",
          res$A_to_B$causal_effect, " (p=", res$A_to_B$p_value, ")\n", sep = "")
      cat("  ", x$traits[2], "->", x$traits[1], ": ",
          res$B_to_A$causal_effect, " (p=", res$B_to_A$p_value, ")\n", sep = "")
    }
  }
}
