#' Compute contribution ratio
#'
#' The relative contribution of independent genetic effects:
#' rho = (theta_ind^2) / (theta_ind^2 + theta_shared^2)
#'
#' @param theta_ind Independent effect estimates
#' @param theta_shared Shared effect estimates (optional)
#'
#' @return Data frame with contribution ratios
#' @export
#'
compute_contribution_ratio <- function(theta_ind, theta_shared = NULL) {

  if (is.null(theta_shared)) {
    # If no shared effects, all variance is from independent
    return(rep(1, length(theta_ind)))
  }

  ind_sq <- theta_ind^2
  shared_sq <- theta_shared^2

  rho <- ind_sq / (ind_sq + shared_sq)
  rho[is.na(rho)] <- 0
  rho[is.nan(rho)] <- 0

  return(rho)
}


#' Compute conditional variance ratio (CVR)
#'
#' CVR = Var(y_i| -i) / Var(y_i)
#' Quantifies proportion of genetic variation independent of other traits.
#'
#' @param Y Phenotype matrix (n x m)
#'
#' @return Data frame with CVR for each trait
#' @export
#'
compute_cvr <- function(Y) {

  n <- nrow(Y)
  m <- ncol(Y)
  trait_names <- colnames(Y)

  cvr_results <- data.frame(
    Trait = trait_names,
    CVR = NA,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(m)) {
    # Marginal variance
    var_marginal <- var(Y[, i])

    # Conditional variance using conditional_phenotype
    pheno_df <- data.frame(GID = paste0("ind_", 1:n), Y)
    colnames(pheno_df)[-1] <- trait_names

    cond_result <- tryCatch({
      conditional_phenotype(pheno_df, trait_names[i])
    }, error = function(e) NULL)

    if (!is.null(cond_result)) {
      y_cond <- cond_result$pheno_cond[, 2]
      var_cond <- var(y_cond, na.rm = TRUE)

      if (var_marginal > 0) {
        cvr_results$CVR[i] <- var_cond / var_marginal
      }
    }
  }

  return(cvr_results)
}


#' Compute detection gain
#'
#' Gain = |{l: P_cond < alpha2, P_marg >= alpha1}| / |{l: P_marg < alpha1}|
#' Measures proportion of loci detectable only after removing cross-trait confounding.
#'
#' @param pval_marginal Marginal p-values
#' @param pval_cond Conditional p-values
#' @param alpha1 Marginal significance threshold
#' @param alpha2 Conditional significance threshold
#'
#' @return Detection gain value
#' @export
#'
compute_detection_gain <- function(pval_marginal, pval_cond,
                                  alpha1 = 0.05, alpha2 = 0.05) {

  n <- length(pval_marginal)

  # Significant in marginal
  sig_marg <- pval_marginal < alpha1
  n_marg <- sum(sig_marg)

  # Significant in conditional but not marginal
  sig_cond_not_marg <- (pval_cond < alpha2) & (pval_marginal >= alpha1)
  n_cond_not_marg <- sum(sig_cond_not_marg)

  if (n_marg == 0) {
    return(0)
  }

  gain <- n_cond_not_marg / n_marg
  return(gain)
}


#' Compute power and FDR for QTL detection
#'
#' @param pval_obs Observed p-values
#' @param true_labels True classification labels
#' @param alpha Significance threshold
#'
#' @return List with power, FDR, and other metrics
#' @export
#'
compute_power_fdr <- function(pval_obs, true_labels, alpha = 0.05) {

  # True positives: significant and truly associated
  tp <- sum(pval_obs < alpha & true_labels %in% c("trait_specific", "shared", "independent"))

  # False positives: significant but not truly associated
  fp <- sum(pval_obs < alpha & true_labels == "null")

  # True negatives
  tn <- sum(pval_obs >= alpha & true_labels == "null")

  # False negatives
  fn <- sum(pval_obs >= alpha & true_labels %in% c("trait_specific", "shared", "independent"))

  # Power = TP / (TP + FN)
  power <- if ((tp + fn) > 0) tp / (tp + fn) else 0

  # FDR = FP / (TP + FP)
  fdr <- if ((tp + fp) > 0) fp / (tp + fp) else 0

  # Type I error = FP / (FP + TN)
  typeI <- if ((fp + tn) > 0) fp / (fp + tn) else 0

  return(list(
    power = power,
    fdr = fdr,
    typeI_error = typeI,
    tp = tp,
    fp = fp,
    tn = tn,
    fn = fn
  ))
}


#' Compute classification accuracy with F1 scores
#'
#' Computes overall accuracy and per-class F1 scores (method.md Section 0.3).
#'
#' @param predicted_labels Predicted QTL categories (5 classes: not_detected, covariance_induced, trait_specific, horizontal_pleiotropy, vertical_pleiotropy)
#' @param true_labels True QTL categories
#'
#' @return List with accuracy, confusion matrix, per-class precision/recall/F1
#' @export
#'
compute_classification_accuracy <- function(predicted_labels, true_labels) {

  # Overall accuracy: Acc = (1/N) * sum(I(C_hat == C_true))
  accuracy <- sum(predicted_labels == true_labels) / length(true_labels)

  # Confusion matrix
  confusion <- table(True = true_labels, Predicted = predicted_labels)

  # Per-class metrics
  classes <- unique(c(true_labels, predicted_labels))

  per_class_metrics <- data.frame(
    class = classes,
    precision = NA,
    recall = NA,
    f1_score = NA,
    stringsAsFactors = FALSE
  )

  for (c in classes) {
    # TP: predicted c, true c
    tp <- sum(predicted_labels == c & true_labels == c)

    # FP: predicted c, true not c
    fp <- sum(predicted_labels == c & true_labels != c)

    # FN: predicted not c, true c
    fn <- sum(predicted_labels != c & true_labels == c)

    # Precision = TP / (TP + FP)
    precision <- if ((tp + fp) > 0) tp / (tp + fp) else 0

    # Recall = TP / (TP + FN)
    recall <- if ((tp + fn) > 0) tp / (tp + fn) else 0

    # F1 = 2 * Precision * Recall / (Precision + Recall)
    f1 <- if ((precision + recall) > 0) 2 * precision * recall / (precision + recall) else 0

    per_class_metrics$precision[per_class_metrics$class == c] <- precision
    per_class_metrics$recall[per_class_metrics$class == c] <- recall
    per_class_metrics$f1_score[per_class_metrics$class == c] <- f1
  }

  return(list(
    accuracy = accuracy,
    confusion_matrix = confusion,
    per_class_metrics = per_class_metrics
  ))
}


#' Compare estimated effects with oracle (ground truth)
#'
#' Computes bias, RMSE, and correlation for effect decomposition.
#'
#' @param beta_ind_est Estimated independent effects (n_snp x n_trait)
#' @param beta_ind_oracle Oracle independent effects (n_snp x n_trait)
#' @param beta_shared_est Estimated shared effects (optional)
#' @param beta_shared_oracle Oracle shared effects (optional)
#'
#' @return List with evaluation metrics
#' @export
#'
compare_estimation_accuracy <- function(
    beta_ind_est,
    beta_ind_oracle,
    beta_shared_est = NULL,
    beta_shared_oracle = NULL
) {

  # Compute bias
  bias_ind <- beta_ind_est - beta_ind_oracle

  # Metrics for beta_ind
  metrics_ind <- data.frame(
    trait = colnames(beta_ind_est),
    mean_abs_bias = colMeans(abs(bias_ind)),
    rmse = sqrt(colMeans(bias_ind^2)),
    correlation = sapply(seq_len(ncol(beta_ind_est)), function(i) {
      cor(beta_ind_est[, i], beta_ind_oracle[, i])
    })
  )

  # If beta_shared provided, compute metrics
  metrics_shared <- NULL
  if (!is.null(beta_shared_est) && !is.null(beta_shared_oracle)) {
    bias_shared <- beta_shared_est - beta_shared_oracle

    metrics_shared <- data.frame(
      trait = colnames(beta_shared_est),
      mean_abs_bias = colMeans(abs(bias_shared)),
      rmse = sqrt(colMeans(bias_shared^2)),
      correlation = sapply(seq_len(ncol(beta_shared_est)), function(i) {
        cor(beta_shared_est[, i], beta_shared_oracle[, i])
      })
    )
  }

  # Return
  list(
    beta_ind = metrics_ind,
    beta_shared = metrics_shared,
    bias_ind_matrix = bias_ind
  )
}


#' Plot estimated vs oracle effects
#'
#' @param beta_est Estimated effects (n_snp x n_trait)
#' @param beta_oracle Oracle effects (n_snp x n_trait)
#' @param true_labels Ground truth labels (optional, for coloring)
#' @param title Plot title
#'
#' @export
plot_estimation_accuracy <- function(
    beta_est,
    beta_oracle,
    true_labels = NULL,
    title = "Estimated vs Oracle Effects"
) {
  
  m <- ncol(beta_est)
  
  par(mfrow = c(1, min(m, 3)))
  
  for (i in seq_len(min(m, 3))) {
    
    # Colors
    if (!is.null(true_labels)) {
      colors <- ifelse(true_labels == "trait_specific", "red",
                ifelse(true_labels == "shared", "blue",
                ifelse(true_labels == "independent", "green", "gray")))
    } else {
      colors <- "black"
    }
    
    # Plot
    plot(
      beta_oracle[, i],
      beta_est[, i],
      xlab = "Oracle Effect",
      ylab = "Estimated Effect",
      main = sprintf("%s - Trait %d", title, i),
      pch = 16,
      col = colors,
      cex = 0.7
    )
    
    # Identity line
    abline(0, 1, col = "black", lty = 2, lwd = 2)
    
    # Correlation
    r <- cor(beta_oracle[, i], beta_est[, i])
    text(
      x = min(beta_oracle[, i]), 
      y = max(beta_est[, i]),
      labels = sprintf("r = %.3f", r),
      adj = c(0, 1),
      cex = 0.9
    )
    
    # Legend (first plot only)
    if (i == 1 && !is.null(true_labels)) {
      legend("topleft",
             legend = c("Trait-specific", "Shared", "Independent", "Null"),
             col = c("red", "blue", "green", "gray"),
             pch = 16,
             cex = 0.6)
    }
  }
}
