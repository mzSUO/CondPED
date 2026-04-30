#' CondPED: Visualization Functions
#'
#' This module provides visualization functions for GWAS results.
#'
#' @name visualize-functions
#' @docType package
NULL


#' Manhattan Plot
#'
#' @param gwas_results GWAS results data frame with SNP, P_Value columns
#' @param snp_info SNP information data frame (with CHR, BP columns)
#' @param sig_threshold Significance threshold (default: 0.05/n_snp)
#' @param title Plot title
#'
#' @export
#'
manhattan_plot <- function(gwas_results, snp_info, sig_threshold = NULL, title = "Manhattan Plot") {

  # Merge with SNP info
  df <- merge(gwas_results, snp_info, by.x = "SNP", by.y = "SNP", all.x = TRUE)

  # Set default threshold
  if (is.null(sig_threshold)) {
    sig_threshold <- 0.05 / nrow(snp_info)
  }

  # Calculate -log10(p)
  df$neg_log10_p <- -log10(df$P_Value)

  # Prepare chromosome positions
  df <- df[order(df$CHR, df$BP), ]

  # Create cumulative position
  df$pos_cum <- 0
  chr_pos <- 0
  for (chr in unique(df$CHR)) {
    chr_max <- max(df$BP[df$CHR == chr])
    df$pos_cum[df$CHR == chr] <- df$BP[df$CHR == chr] + chr_pos
    chr_pos <- chr_pos + chr_max
  }

  # Colors
  df$color <- ifelse(df$CHR %% 2 == 0, "gray1", "gray2")

  # Significance line
  sig_line <- -log10(sig_threshold)

  # Plot
  par(mar = c(5, 5, 4, 2))
  plot(df$pos_cum, df$neg_log10_p,
       col = df$color, pch = 16, cex = 0.6,
       xlab = "Chromosome", ylab = expression(-log[10](P)),
       main = title, cex.main = 1.2, cex.lab = 1.1)

  # Add significance line
  abline(h = sig_line, col = "red", lty = 2, lwd = 1.5)

  # Add chromosome labels
  chr_centers <- tapply(df$pos_cum, df$CHR, mean)
  axis(side = 1, at = chr_centers, labels = names(chr_centers), cex.axis = 0.8)

  # Highlight significant SNPs
  sig_snps <- df[df$P_Value < sig_threshold, ]
  if (nrow(sig_snps) > 0) {
    points(sig_snps$pos_cum, sig_snps$neg_log10_p, col = "red", pch = 16, cex = 0.8)
  }
}


#' QQ Plot
#'
#' @param pvalues Vector of p-values
#' @param title Plot title
#'
#' @export
#'
qq_plot <- function(pvalues, title = "Q-Q Plot") {

  # Observed -log10(p)
  obs <- -log10(sort(pvalues))

  # Expected under null
  n <- length(pvalues)
  exp <- -log10(ppoints(n))

  # Plot
  par(mar = c(5, 5, 4, 2))
  plot(exp, obs, pch = 16, cex = 0.6,
       xlab = expression("Expected -log"[10](P)),
       ylab = expression("Observed -log"[10](P)),
       main = title, cex.main = 1.2, cex.lab = 1.1)

  # Add confidence interval (95%)
  n <- length(pvalues)
  ci <- 0.05
  lines(exp, exp - qnorm(1 - ci/2) / sqrt(n), col = "gray", lty = 2)
  lines(exp, exp + qnorm(1 - ci/2) / sqrt(n), col = "gray", lty = 2)

  # Add diagonal
  abline(0, 1, col = "red", lwd = 2)
}


#' Classification Summary Plot
#'
#' @param classification_results Output from classify_qtl()
#' @param title Plot title
#'
#' @export
#'
plot_classification_summary <- function(classification_results, title = "QTL Classification Summary") {

  # Summary table
  summary_df <- table(classification_results$classification, classification_results$Trait)

  # Colors for categories
  colors <- c("not_detected" = "gray",
              "shared" = "blue",
              "independent" = "green",
              "trait_specific" = "red")

  # Bar plot
  par(mar = c(8, 5, 4, 2))
  barplot(summary_df, beside = TRUE, col = rownames(summary_df),
          xlab = "Trait", ylab = "Count",
          main = title, cex.main = 1.2, cex.lab = 1.1,
          legend = TRUE, args.legend = list(x = "topright", cex = 0.8))
}


#' Effect Size Comparison Plot
#'
#' @param marginal_effects Marginal effect sizes
#' @param conditional_effects Conditional effect sizes
#' @param snp_names SNP names
#' @param title Plot title
#'
#' @export
#'
plot_effect_comparison <- function(marginal_effects, conditional_effects,
                                   snp_names = NULL, title = "Effect Size Comparison") {

  # Create data frame
  df <- data.frame(
    Marginal = marginal_effects,
    Conditional = conditional_effects
  )

  if (!is.null(snp_names)) {
    df$SNP <- snp_names
  }

  # Plot
  par(mar = c(5, 5, 4, 2))
  plot(df$Marginal, df$Conditional, pch = 16, cex = 0.7,
       xlab = "Marginal Effect", ylab = "Conditional Effect",
       main = title, cex.main = 1.2, cex.lab = 1.1)

  # Add diagonal line
  lim <- c(min(c(df$Marginal, df$Conditional), na.rm = TRUE),
           max(c(df$Marginal, df$Conditional), na.rm = TRUE))
  lines(lim, lim, col = "red", lty = 2, lwd = 2)

  # Label notable SNPs
  if (!is.null(snp_names)) {
    diff <- abs(df$Marginal - df$Conditional)
    notable <- diff > quantile(diff, 0.9, na.rm = TRUE)
    if (any(notable)) {
      text(df$Marginal[notable], df$Conditional[notable],
           df$SNP[notable], cex = 0.7, pos = 4)
    }
  }
}


#' Detection Gain Visualization
#'
#' @param marginal_count Number of significant SNPs from marginal GWAS
#' @param conditional_count Number of significant SNPs from conditional GWAS
#' @param title Plot title
#'
#' @export
#'
plot_detection_gain <- function(marginal_count, conditional_count, title = "Detection Gain") {

  df <- data.frame(
    Method = c("Marginal GWAS", "Conditional GWAS"),
    Count = c(marginal_count, conditional_count)
  )

  par(mar = c(5, 5, 4, 2))
  barplot(df$Count, names.arg = df$Method, col = c("steelblue", "darkgreen"),
          xlab = "Method", ylab = "Significant SNP Count",
          main = title, cex.main = 1.2, cex.lab = 1.1)

  # Add count labels
  text(c(1, 2), df$Count + 0.5, labels = df$Count, cex = 1.2)

  # Calculate gain
  if (marginal_count > 0) {
    gain <- (conditional_count - marginal_count) / marginal_count * 100
    mtext(sprintf("Gain: %.1f%%", gain), side = 3, cex = 0.9)
  }
}


#' Multi-trait Effect Heatmap
#'
#' @param effect_matrix Effect size matrix (n_snp x n_trait)
#' @param snp_names SNP names
#' @param trait_names Trait names
#' @param title Plot title
#'
#' @export
#'
plot_effect_heatmap <- function(effect_matrix, snp_names = NULL, trait_names = NULL,
                               title = "Effect Size Heatmap") {

  # Set names
  if (is.null(rownames(effect_matrix)) && !is.null(snp_names)) {
    rownames(effect_matrix) <- snp_names
  }
  if (is.null(colnames(effect_matrix)) && !is.null(trait_names)) {
    colnames(effect_matrix) <- trait_names
  }

  # Load or install heatmapr if needed
  if (!requireNamespace("gplots", quietly = TRUE)) {
    message("Installing gplots...")
    install.packages("gplots")
  }

  library(gplots)

  # Create heatmap
  par(mar = c(10, 10, 4, 2))
  heatmap.2(effect_matrix, trace = "none", dendrogram = "none",
            col = colorRampPalette(c("blue", "white", "red"))(100),
            main = title, cex.main = 1.2,
            xlab = "Trait", ylab = "SNP",
            cex.lab = 0.8, cex.axis = 0.7)
}
