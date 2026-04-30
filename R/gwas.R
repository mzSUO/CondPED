#' GWAS Interface for CondPED using QTLNetwork
#'
#' This module provides integration between the conditional projection framework
#' and QTLNetwork for multi-trait GWAS analysis.
#'
#' @name gwas-functions
#' @docType package
NULL


#' Convert genotype data to QTLNetwork format
#'
#' Follows qtxnetwork.input.trans format:
#' - Uses original SNP names from snp_info$SNP
#' - Uses original individual IDs from rownames of G
#'
#' @param snp_info SNP information data frame (with CHR, SNP, BP columns)
#' @param G Genotype matrix (n_ind x n_snp), values -1, 0, 1
#' @param output_prefix Output file prefix
#'
#' @return NULL (writes files)
#' @export
#' @importFrom utils write.table
#'
convert_geno_to_qtlnetwork <- function(snp_info, G, output_prefix) {

  # QTLNetwork .gen format (following qtxnetwork.input.trans):
  # #Ind  SNP1_name  SNP2_name  ...  SNPn_name
  # where each row is an individual

  n_ind <- nrow(G)
  n_snp <- ncol(G)

  # Get individual IDs from rownames of G
  gid <- rownames(G)

  # Convert standardized G back to -1, 0, 1 coding
  G_raw <- round(G)

  # Use original SNP names from snp_info
  # If SNP column exists, use it; otherwise use column names from G
  if (!is.null(snp_info$SNP) && length(snp_info$SNP) == n_snp) {
    snp_names <- snp_info$SNP
  } else {
    snp_names <- colnames(G)
  }

  # Create genotype matrix with original SNP names
  colnames(G_raw) <- snp_names

  # Create data frame: each row is an individual
  # First column is ID (#Ind), then genotypes for all SNPs
  result <- data.frame(
    GID = gid,
    stringsAsFactors = FALSE
  )

  # Add genotype columns
  result <- cbind(result, G_raw)

  # Rename first column to #Ind (QTXNetwork format)
  colnames(result)[1] <- "#Ind"

  # Write .gen file
  gen_file <- paste0(output_prefix, ".gen")
  printer <- file(gen_file, "w")
  write.table(result, file = printer, append = TRUE, na = ".",
              row.names = FALSE, col.names = TRUE, quote = FALSE,
              sep = "\t", eol = "\n")
  close(printer)

  # Return chromosome summary for phenotype conversion
  chr_summary <- table(snp_info$CHR)
  chr_names <- names(chr_summary)

  return(invisible(list(
    gen_file = gen_file,
    chr = chr_summary,
    chr_names = chr_names,
    n_snp = n_snp,
    gid = gid,
    snp_names = snp_names
  )))
}


#' Convert phenotype data to QTLNetwork format
#'
#' @param pheno_data Phenotype data frame with GID and trait columns
#' @param trait_name Trait name to process
#' @param output_prefix Output file prefix
#' @param chr_summary Chromosome summary from genotype conversion (optional)
#'
#' @return NULL (writes files)
#' @export
#'
convert_pheno_to_qtlnetwork <- function(pheno_data, trait_name, output_prefix, chr_summary = NULL) {

  # Filter and format phenotype - keep original GID format
  phe <- pheno_data[, c("GID", trait_name)]
  phe <- phe[complete.cases(phe), ]

  n_ind <- nrow(phe)
  n_trait <- 1

  # Use provided chromosome info or default
  if (is.null(chr_summary)) {
    n_chr <- 1
    chr_names <- "1"
    n_markers <- 100
    chr_markers <- "100"
  } else {
    n_chr <- length(chr_summary)
    chr_names <- paste(chr_names <- names(chr_summary), collapse = " ")
    n_markers <- sum(chr_summary)
    chr_markers <- paste(chr_summary, collapse = " ")
  }

  output_file <- paste0(output_prefix, "_", trait_name, ".phe")
  printer <- file(output_file, "w")

  # Header (following reference code format exactly)
  write(paste0("_Population\tRIL"), printer)
  write(paste0("_Genotypes\t", n_ind), printer, append = TRUE)
  write(paste0("_Observations\t", n_ind), printer, append = TRUE)
  write("_Environments\tno", printer, append = TRUE)
  write("_Replications\tno", printer, append = TRUE)
  write(paste0("_TraitNumber\t", n_trait), printer, append = TRUE)
  write(paste0("_Chromosomes\t", n_chr, "\t", chr_names), printer, append = TRUE)
  write(paste0("_TotalMarker\t", n_markers, "\t", chr_markers), printer, append = TRUE)
  write("_MarkerCode\tP1=1\tP2=-1\tF1=0\n", printer, append = TRUE)
  write("*TraitBegin*", printer, append = TRUE)

  # Data rows (with ";" as end-of-line, matching reference)
  write.table(phe[, c("GID", trait_name)], file = printer, append = TRUE, na = ".",
              row.names = FALSE, col.names = TRUE, quote = FALSE,
              sep = "\t", eol = ";\n")

  write("*TraitEnd*", printer, append = TRUE)
  close(printer)

  return(invisible(output_file))
}


#' Run QTLNetwork GWAS
#'
#' @param qtxnetwork_path Path to QTLNetwork executable
#' @param gen_file Path to .gen file
#' @param phe_file Path to .phe file
#' @param output_prefix Output file prefix
#' @param scan_2d Whether to perform 2D scanning for epistasis (default: FALSE)
#' @param verbose Print command (default: TRUE)
#'
#' @return NULL (writes output file)
#' @export
#'
run_qtlnetwork <- function(qtxnetwork_path, gen_file, phe_file, output_prefix,
                          scan_2d = FALSE, verbose = TRUE) {

  # Build command
  if (scan_2d) {
    # 2D scan: include epistasis
    cmd <- paste0(qtxnetwork_path,
                  " --map ", gen_file,
                  " --txt ", phe_file,
                  " --out ", output_prefix,
                  " --QTX 1")
  } else {
    # 1D scan only: no epistasis
    cmd <- paste0(qtxnetwork_path,
                  " --map ", gen_file,
                  " --txt ", phe_file,
                  " --out ", output_prefix,
                  " --QTX 1 --only-1D 1")
  }

  if (verbose) {
    cat("*** QTXNetwork command:\n")
    cat(cmd, "\n")
  }

  # Run QTLNetwork
  system(cmd)

  return(invisible(paste0(output_prefix, ".pre")))
}


#' Parse QTLNetwork output (following qtxnetwork.output.trans format)
#'
#' Parses 1D (additive, type A) and 2D (epistatic, type AA) effects.
#'
#' @param pre_file Path to .pre output file
#' @param trait_name Trait name
#'
#' @return List with:
#'   \item{qtl_data}{Data frame with 1D additive QTL (A type)}
#'   \item{qtl_2d_data}{Data frame with 2D epistatic QTL (AA type, if present)}
#' @export
#'
parse_qtlnetwork_output <- function(pre_file, trait_name) {

  if (!file.exists(pre_file)) {
    stop("Output file not found: ", pre_file)
  }

  lines <- readLines(pre_file)

  # ==========================================================================
  # Parse 1D effects: columns are "QTL", "SNPID", "A", "SE", "P-Value"
  # ==========================================================================
  start_idx <- which(lines == "_1D_effect")
  end_idx <- which(lines == "_1D_heritability")

  qtl_data <- data.frame(
    TRAIT = character(),
    QTL = character(),
    SNPID = character(),
    A = numeric(),
    SE = numeric(),
    P_Value = numeric(),
    stringsAsFactors = FALSE
  )

  if (length(start_idx) > 0 && length(end_idx) > 0) {
    effect_lines <- lines[(start_idx + 2):(end_idx - 2)]

    for (line in effect_lines) {
      if (nchar(trimws(line)) == 0) next

      parts <- strsplit(trimws(line), "\\s+")[[1]]
      # Format: "QTL", "SNPID", "A", "SE", "P-Value"
      if (length(parts) >= 5) {
        qtl_data <- rbind(qtl_data, data.frame(
          TRAIT = trait_name,
          QTL = parts[1],
          SNPID = parts[1],
          A = as.numeric(parts[2]),
          SE = as.numeric(parts[3]),
          P_Value = as.numeric(parts[4]),
          stringsAsFactors = FALSE
        ))
      }
    }
  } else {
    warning("No 1D QTL effects found in output")
  }

  # ==========================================================================
  # Parse 2D effects: columns are "QTL", "SNPID", "AA", "SE", "P-Value"
  # ==========================================================================
  qtl_2d_data <- data.frame(
    TRAIT = character(),
    QTL = character(),
    SNPID = character(),
    AA = numeric(),
    SE = numeric(),
    P_Value = numeric(),
    stringsAsFactors = FALSE
  )

  start_idx_2d <- which(lines == "_2D_effect")
  end_idx_2d <- which(lines == "_2D_heritability")

  if (length(start_idx_2d) > 0 && length(end_idx_2d) > 0) {
    effect_lines_2d <- lines[(start_idx_2d + 2):(end_idx_2d - 2)]

    for (line in effect_lines_2d) {
      if (nchar(trimws(line)) == 0) next

      parts <- strsplit(trimws(line), "\\s+")[[1]]
      # Format: "QTL", "SNPID", "AA", "SE", "P-Value"
      if (length(parts) >= 5) {
        qtl_2d_data <- rbind(qtl_2d_data, data.frame(
          TRAIT = trait_name,
          QTL = paste(parts[1], parts[2], sep = "_x_"),
          SNPID = paste(parts[1], parts[2], sep = "_x_"),
          AA = as.numeric(parts[2]),
          SE = as.numeric(parts[3]),
          P_Value = as.numeric(parts[4]),
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  return(list(
    qtl_data = qtl_data,
    qtl_2d_data = qtl_2d_data
  ))
}


#' Marginal GWAS using QTLNetwork
#'
#' Performs genome-wide scan on original phenotypes to detect significant loci.
#'
#' @param sim_data Output from simulate_multitrait_data()
#' @param qtxnetwork_path Path to QTLNetwork executable
#' @param output_dir Output directory (default: tempdir())
#' @param alpha1 Significance threshold (default: 0.05/n_snp)
#' @param scan_2d Whether to perform 2D scanning for epistasis (default: FALSE)
#'
#' @return List with:
#'   \item{gwas_results}{Data frame with all SNP 1D effects and p-values}
#'   \item{gwas_results_2d}{Data frame with 2D epistatic effects (if scan_2d=TRUE)}
#'   \item{significant_snps}{Vector of significant SNP indices}
#'
#' @export
#'
marginal_gwas <- function(sim_data,
                          qtxnetwork_path,
                          output_dir = tempdir(),
                          alpha1 = NULL,
                          scan_2d = FALSE) {

  n_snp <- ncol(sim_data$G)

  if (is.null(alpha1)) {
    alpha1 <- 0.05 / n_snp
  }

  # Create output prefix
  output_prefix <- file.path(output_dir, "marginal_gwas")

  # Convert genotype (use G_raw for QTLNetwork, not standardized G)
  geno_result <- convert_geno_to_qtlnetwork(
    sim_data$snp_info,
    sim_data$G_raw,
    output_prefix
  )
  gen_file <- geno_result$gen_file
  chr_summary <- geno_result$chr

  # Run GWAS for each trait
  all_results_1d <- data.frame()
  all_results_2d <- data.frame()

  for (trait in colnames(sim_data$Y)) {
    # Convert phenotype
    pheno_df <- data.frame(GID = rownames(sim_data$Y), sim_data$Y)
    phe_file <- convert_pheno_to_qtlnetwork(pheno_df, trait, output_prefix, chr_summary)

    # Run QTLNetwork with optional 2D scanning
    pre_file <- run_qtlnetwork(qtxnetwork_path, gen_file, phe_file, output_prefix, scan_2d = scan_2d)

    # Parse results (1D always, 2D if scan_2d=TRUE)
    parsed <- parse_qtlnetwork_output(pre_file, trait)
    all_results_1d <- rbind(all_results_1d, parsed$qtl_data)

    # Add 2D results if present
    if (nrow(parsed$qtl_2d_data) > 0) {
      all_results_2d <- rbind(all_results_2d, parsed$qtl_2d_data)
    }
  }

  # Identify significant SNPs (use P_Value to match new format)
  sig_snps <- which(all_results_1d$P_Value < alpha1)

  result <- list(
    gwas_results = all_results_1d,
    significant_snps = sig_snps,
    threshold = alpha1
  )

  if (nrow(all_results_2d) > 0) {
    result$gwas_results_2d <- all_results_2d
  }

  return(result)
}


#' Forward Conditional GWAS (Step 2)
#'
#' Performs conditional GWAS on candidate SNP set using conditional phenotypes.
#' This implements Step 2 of the bidirectional conditional GWAS framework.
#' The conditional phenotype is constructed as:
#'   y_i* = y_i - y_{-i} * gamma
#' where gamma is estimated from the phenotypic covariance matrix.
#'
#' @param sim_data Output from simulate_multitrait_data()
#' @param target_trait Target trait name
#' @param candidate_snps Character vector of candidate SNP names (from Step 1)
#' @param alpha2 Significance threshold (default: 0.05 / (|L| * m))
#'
#' @return Data frame with conditional GWAS results
#' @export
#'
forward_conditional_gwas <- function(sim_data,
                                     target_trait,
                                     candidate_snps,
                                     alpha2 = NULL) {

  n_snp <- ncol(sim_data$G)
  n_trait <- ncol(sim_data$Y)
  n_ind <- nrow(sim_data$Y)

  # Get candidate SNP indices
  snp_names <- colnames(sim_data$G)
  snp_idx <- which(snp_names %in% candidate_snps)

  if (length(snp_idx) == 0) {
    warning("No candidate SNPs found in genotype matrix")
    return(data.frame(
      TRAIT = character(),
      SNP = character(),
      Effect = numeric(),
      SE = numeric(),
      P_Value = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  if (is.null(alpha2)) {
    alpha2 <- 0.05 / (length(snp_idx) * n_trait)
  }

  # Step 1: Construct conditional phenotype
  # y_i* = y_i - y_{-i} * gamma
  # where gamma = V_{-i}^{-1} * Cov(y_i, y_{-i})
  pheno_df <- data.frame(GID = rownames(sim_data$Y), sim_data$Y)
  cond_result <- conditional_phenotype(pheno_df, target_trait)

  # Get conditional phenotype values
  y_cond <- cond_result$pheno_cond[, 2]
  G <- sim_data$G

  # Step 2: Run GWAS on conditional phenotype for candidate SNPs only
  results <- data.frame()

  for (snp in snp_names[snp_idx]) {
    x <- G[, snp]

    # Simple linear regression: y_cond ~ x
    fit <- lm(y_cond ~ x)

    # Extract results
    coef_summary <- summary(fit)$coefficients
    if (nrow(coef_summary) >= 2) {
      beta <- coef_summary["x", "Estimate"]
      se <- coef_summary["x", "Std. Error"]
      pval <- coef_summary["x", "Pr(>|t|)"]

      results <- rbind(results, data.frame(
        TRAIT = target_trait,
        SNP = snp,
        Effect = beta,
        SE = se,
        P_Value = pval,
        stringsAsFactors = FALSE
      ))
    }
  }

  # Add threshold to result
  attr(results, "threshold") <- alpha2

  return(results)
}


#' Reverse Conditional GWAS (Step 3)
#'
#' Tests whether SNP effects on other traits are mediated through target trait.
#' This implements Step 3 of the bidirectional conditional GWAS framework.
#' For each (l, i) in S1 and each t != i, we construct reverse conditional phenotype
#' by conditioning on all traits except i.
#'
#' @param sim_data Output from simulate_multitrait_data()
#' @param target_trait Target trait (the trait we want to test specificity for)
#' @param test_trait Trait to test for mediation
#' @param candidate_snps Character vector of candidate SNP names (from forward step)
#' @param alpha3 Significance threshold
#'
#' @return Data frame with reverse conditional GWAS results
#' @export
#'
reverse_conditional_gwas <- function(sim_data,
                                     target_trait,
                                     test_trait,
                                     candidate_snps,
                                     alpha3 = NULL) {

  n_snp <- ncol(sim_data$G)
  n_trait <- ncol(sim_data$Y)
  n_ind <- nrow(sim_data$Y)

  # Get candidate SNP indices
  snp_names <- colnames(sim_data$G)
  snp_idx <- which(snp_names %in% candidate_snps)

  if (length(snp_idx) == 0) {
    warning("No candidate SNPs found in genotype matrix")
    return(data.frame(
      TRAIT = character(),
      SNP = character(),
      Effect = numeric(),
      SE = numeric(),
      P_Value = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  if (is.null(alpha3)) {
    # Default: FDR controlled at 0.05
    alpha3 <- 0.05
  }

  # Construct reverse conditional phenotype:
  # y_{t|(-i)} = y_t - y_{-i} * gamma_{t,-i}
  # where -i means all traits except target_trait
  pheno_df <- data.frame(GID = rownames(sim_data$Y), sim_data$Y)

  # Use reverse_conditional to get phenotype conditioned on all but target_trait
  reverse_result <- reverse_conditional(pheno_df,
                                         target_trait = test_trait,
                                         exclude_trait = target_trait)

  # Get reverse conditional phenotype
  y_rev <- reverse_result$pheno_cond[, 2]
  G <- sim_data$G

  # Run GWAS on reverse conditional phenotype
  results <- data.frame()

  for (snp in snp_names[snp_idx]) {
    x <- G[, snp]

    # Simple linear regression: y_rev ~ x
    fit <- lm(y_rev ~ x)

    # Extract results
    coef_summary <- summary(fit)$coefficients
    if (nrow(coef_summary) >= 2) {
      beta <- coef_summary["x", "Estimate"]
      se <- coef_summary["x", "Std. Error"]
      pval <- coef_summary["x", "Pr(>|t|)"]

      results <- rbind(results, data.frame(
        TRAIT = test_trait,
        SNP = snp,
        Effect = beta,
        SE = se,
        P_Value = pval,
        stringsAsFactors = FALSE
      ))
    }
  }

  # Add threshold to result
  attr(results, "threshold") <- alpha3

  return(results)
}


#' Complete Bidirectional Conditional GWAS Pipeline
#'
#' Implements the full workflow from method.md Section 0.2:
#' 1. Step 1: Marginal multivariate GWAS (using QTLNetwork) -> significant SNP set L
#' 2. Step 2: Forward conditional GWAS (R implementation) -> independent effects S1
#' 3. Step 3: Reverse conditional GWAS (R implementation) -> trait specificity test
#'
#' @param sim_data Output from simulate_multitrait_data()
#' @param qtxnetwork_path Path to QTLNetwork executable (for Step 1 only)
#' @param output_dir Output directory
#' @param scan_2d Whether to perform 2D scanning in Step 1 (default: FALSE)
#'
#' @return List with all GWAS results and classifications
#'
#' @export
#'
run_bidirectional_gwas <- function(sim_data,
                                    qtxnetwork_path,
                                    output_dir = tempdir(),
                                    scan_2d = FALSE) {
  n_snp <- ncol(sim_data$G)
  n_trait <- ncol(sim_data$Y)

  # ==========================================================================
  # Step 1: Marginal multivariate GWAS (using QTLNetwork)
  # ==========================================================================
  cat("=== Step 1: Marginal GWAS (QTLNetwork) ===\n")
  marginal <- marginal_gwas(sim_data, qtxnetwork_path, output_dir, scan_2d = scan_2d)

  # Get candidate SNP set L
  candidate_snps <- unique(marginal$gwas_results$SNP)
  cat("Candidate SNP set L:", length(candidate_snps), "SNPs\n")

  # ==========================================================================
  # Step 2: Forward conditional GWAS (R implementation)
  # ==========================================================================
  cat("\n=== Step 2: Forward Conditional GWAS (R) ===\n")
  forward_results <- list()

  for (trait in colnames(sim_data$Y)) {
    cat("Processing trait:", trait, "\n")
    forward_results[[trait]] <- forward_conditional_gwas(
      sim_data = sim_data,
      target_trait = trait,
      candidate_snps = candidate_snps
    )
  }

  # ==========================================================================
  # Step 3: Reverse conditional GWAS (R implementation)
  # ==========================================================================
  cat("\n=== Step 3: Reverse Conditional GWAS (R) ===\n")
  reverse_results <- list()

  # Get significant pairs from forward step
  alpha2 <- 0.05 / (length(candidate_snps) * n_trait)

  for (target_trait in colnames(sim_data$Y)) {
    # Get significant SNPs for this trait from forward step
    sig_snps <- forward_results[[target_trait]]$SNP[
      forward_results[[target_trait]]$P_Value < alpha2
    ]

    for (test_trait in colnames(sim_data$Y)) {
      if (test_trait != target_trait) {
        cat("Testing:", test_trait, "| exclude:", target_trait, "\n")
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

  return(list(
    step1_marginal = marginal,
    step2_forward = forward_results,
    step3_reverse = reverse_results,
    candidate_snps = candidate_snps
  ))
}
