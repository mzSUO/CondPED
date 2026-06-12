# qtlnetwork.R
# QTLNetwork wrapper for RIL GWAS (additive + AA epistasis, no dominance)
# Compatible with simulation.R output (RIL genotype: -1/1)

# ============================================================================
# 1. Execution
# ============================================================================

#' qtxnetwork.perform
#' Execute QTXNetwork for 1D and/or 2D scanning.
#' 
#' @param qtxnetwork_path Path to QTXNetwork executable.
#' @param gen_file Input .gen file.
#' @param phe_file Input .phe file.
#' @param pre_file Output .pre file.
#' @param scan_2d Logical. TRUE enables 2D AA epistasis scanning. Default FALSE.
#' @param qtx_mode QTX mode. Default 2 (multi-trait).
#' 
#' @return NULL (invisible)
#' @export
qtxnetwork.perform <- function(qtxnetwork_path, gen_file, phe_file, pre_file,
                               scan_2d = FALSE, qtx_mode = 2) {

  if (scan_2d) {
    cmd <- paste0(qtxnetwork_path, " --map ", gen_file,
                  " --txt ", phe_file, " --out ", pre_file, " --QTX ", qtx_mode)
  } else {
    cmd <- paste0(qtxnetwork_path, " --map ", gen_file,
                  " --txt ", phe_file, " --out ", pre_file,
                  " --QTX ", qtx_mode, " --only-1D 1")
  }

  cat("*** QTXNetwork command:\n", cmd, "\n", sep = "")
  system(cmd)
  invisible(NULL)
}


# ============================================================================
# 2. Input Conversion
# ============================================================================

#' qtxnetwork.input.trans
#' Convert data frames to QTXNetwork .gen + .phe format.
#' 
#' @param geno_data Data frame: chr, snp_id, pos, Ind1, Ind2, ...
#'   Genotype coding must match population type.
#' @param pheno_data Data frame: id, Trait1, Trait2, ...
#' @param geno_output_prefix Prefix for .gen file.
#' @param pheno_output_prefix Prefix for .phe file.
#' @param population "RIL" (default, no dominance) or "F2".
#' 
#' @return NULL (invisible)
#' @export
qtxnetwork.input.trans <- function(geno_data, pheno_data,
                                   geno_output_prefix, pheno_output_prefix,
                                   population = c("RIL", "F2")) {

  population <- match.arg(population)

  chr <- table(geno_data[, 1])
  chrName <- names(chr)
  chrCounts <- as.vector(chr)

  # 提取基因型矩阵（去掉前3列：chr, snp_id, pos），行=SNP，列=个体
  geno_mat <- as.matrix(geno_data[, -c(1:3)])
  rownames(geno_mat) <- geno_data[, 2]  # SNP名作为行名
  
  # 转置为：行=个体，列=SNP
  geno_t <- t(geno_mat)
  
  # 构造 .gen 数据框：第一列是个体名，列名行第一列为 _Genotype
  geno_df <- as.data.frame(geno_t, stringsAsFactors = FALSE)
  geno_df <- cbind(rownames(geno_df), geno_df)
  colnames(geno_df) <- c("_Genotype", colnames(geno_t))

  traitNum <- ncol(pheno_data) - 1
  n_ind <- nrow(pheno_data)

  # --- .gen file ---
  gen_file <- paste0(geno_output_prefix, ".gen")
  cat("Writing .gen:", gen_file, "\n")

  con <- file(gen_file, "w")
  on.exit(close(con))

  write.table(geno_df, file = con, append = TRUE, na = ".",
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

  close(con)

  # --- .phe file ---
  phe_file <- paste0(pheno_output_prefix, ".phe")
  cat("Writing", traitNum, "traits to .phe:", phe_file, "\n")

  con <- file(phe_file, "w")
  on.exit(close(con))

  writeLines(paste("_Population", population, sep = "\t"), con)
  writeLines(paste("_Genotypes", n_ind, sep = "\t"), con)
  writeLines(paste("_Observations", n_ind, sep = "\t"), con)
  writeLines("_Environments\tno", con)
  writeLines("_Replications\tno", con)
  writeLines(paste("_TraitNumber", traitNum, paste(seq_len(traitNum), collapse = "+"), sep = "\t"), con)
  writeLines(paste("_Chromosomes", length(chrName), paste(chrName, collapse = " "), sep = "\t"), con)
  writeLines(paste("_TotalMarker", sum(chrCounts), paste(chrCounts, collapse = " "), sep = "\t"), con)
  writeLines(ifelse(population == "RIL",
                    "_MarkerCode\tP1=1\tP2=-1",
                    "_MarkerCode\tP1=1\tP2=-1\tF1=0"), con)
  writeLines("*TraitBegin*", con)

  phe_out <- pheno_data
  colnames(phe_out)[1] <- "Geno#"
  write.table(phe_out, file = con, append = TRUE, na = ".",
              row.names = FALSE, col.names = TRUE, quote = FALSE,
              sep = "\t", eol = ";\n")
  writeLines("*TraitEnd*", con)
  invisible(NULL)
}


# ============================================================================
# 3. Adapter: simulation.R -> QTXNetwork
# ============================================================================

#' condped_to_qtlnetwork
#' Bridge simulation.R output to QTXNetwork input.
#' 
#' simulation.R RIL output is already -1/1 coded, so no conversion needed.
#' All markers placed on a single virtual chromosome for GWAS.
#' 
#' @param sim_result Output from generate_condped() or generate_class*().
#' @param geno_output_prefix Prefix for .gen.
#' @param pheno_output_prefix Prefix for .phe.
#' @param population "RIL" (default) or "F2".
#' 
#' @return NULL (invisible)
#' @export
condped_to_qtlnetwork <- function(sim_result,
                                   geno_output_prefix,
                                   pheno_output_prefix,
                                   population = c("RIL", "F2")) {

  population <- match.arg(population)
  X <- sim_result$X    # n x p, already -1/1 for RIL
  Y <- sim_result$Y    # n x m

  n_ind <- nrow(X)
  p <- ncol(X)
  m <- ncol(Y)

  # Genotype: single virtual chromosome (GWAS)
  # 注意：chr 列的值会原样写入 _Chromosomes 行。示例用 "chr1"，你也可传 "1"
  geno_data <- data.frame(
    chr = rep("chr1", p),
    snp_id = colnames(X) %||% paste0("SNP", seq_len(p)),
    pos = seq_len(p),
    t(X),  # transpose to p x n
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  colnames(geno_data)[4:ncol(geno_data)] <- paste0("Ind", seq_len(n_ind))

  # Phenotype
  trait_names <- colnames(Y) %||% paste0("Trait", seq_len(m))
  pheno_data <- data.frame(
    id = paste0("Ind", seq_len(n_ind)),
    Y,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  colnames(pheno_data) <- c("id", trait_names)

  qtxnetwork.input.trans(geno_data, pheno_data,
                        geno_output_prefix, pheno_output_prefix,
                        population = population)
}

# ============================================================================
# 3. Adapter: simulation.R -> QTXNetwork
# ============================================================================

#' condped_to_qtlnetwork
#' Bridge simulation.R output to QTXNetwork input.
#' 
#' simulation.R RIL output is already -1/1 coded, so no conversion needed.
#' All markers placed on a single virtual chromosome for GWAS.
#' 
#' @param sim_result Output from generate_condped() or generate_class*().
#' @param geno_output_prefix Prefix for .gen.
#' @param pheno_output_prefix Prefix for .phe.
#' @param population "RIL" (default) or "F2".
#' 
#' @return NULL (invisible)
#' @export
condped_to_qtlnetwork <- function(sim_result,
                                   geno_output_prefix,
                                   pheno_output_prefix,
                                   population = c("RIL", "F2")) {

  population <- match.arg(population)
  X <- sim_result$X    # n x p, already -1/1 for RIL
  Y <- sim_result$Y    # n x m

  n_ind <- nrow(X)
  p <- ncol(X)
  m <- ncol(Y)

  # Genotype: single virtual chromosome (GWAS)
  geno_data <- data.frame(
    chr = rep("1", p),
    snp_id = colnames(X) %||% paste0("SNP", seq_len(p)),
    pos = seq_len(p),
    t(X),  # transpose to p x n
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  colnames(geno_data)[4:ncol(geno_data)] <- paste0("Ind", seq_len(n_ind))

  # Phenotype
  trait_names <- colnames(Y) %||% paste0("Trait", seq_len(m))
  pheno_data <- data.frame(
    id = paste0("Ind", seq_len(n_ind)),
    Y,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  colnames(pheno_data) <- c("id", trait_names)

  qtxnetwork.input.trans(geno_data, pheno_data,
                        geno_output_prefix, pheno_output_prefix,
                        population = population)
}


# ============================================================================
# 4. Output Parsing
# ============================================================================

#' qtxnetwork.output.trans
#' Parse .pre file into standardized data frames.
#' 
#' @param pheno_data Phenotype data frame (fallback for trait names).
#' @param pre_file Path to .pre file.
#' @param scan_2d Whether 2D scanning was performed. Default FALSE.
#' 
#' @return List: qtl_data, qtl_dom_data, qtl_aa_data (if scan_2d).
#' @export
qtxnetwork.output.trans <- function(pheno_data, pre_file, scan_2d = FALSE) {

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.")
  }

  lines <- readLines(pre_file)

  # --- Parse _trait blocks ---
  trait_blocks <- list()
  tl <- grep("^_trait\\s+\\d+", lines)

  for (i in seq_along(tl)) {
    parts <- strsplit(trimws(lines[tl[i]]), "\\s+")[[1]]
    name <- paste(parts[-(1:2)], collapse = " ")
    end <- if (i < length(tl)) tl[i + 1] - 1 else length(lines)
    trait_blocks[[i]] <- list(name = name, start = tl[i], end = end)
  }

  if (length(trait_blocks) == 0) {
    tnames <- colnames(pheno_data)[-1]
  } else {
    tnames <- sapply(trait_blocks, function(x) x$name)
  }

  # --- Helpers ---
  extract <- function(lines, sm, em) {
    secs <- list(); si <- 0
    for (i in seq_along(lines)) {
      if (trimws(lines[i]) == sm) { si <- i
      } else if (trimws(lines[i]) == em && si > 0) {
        if (i - si > 2) secs[[length(secs) + 1]] <- lines[(si + 2):(i - 1)]
        si <- 0
      }
    }
    secs
  }

  parse_sec <- function(sec) {
    if (length(sec) < 2) return(data.frame())
    hdr <- make.names(strsplit(trimws(sec[1]), "\\s+")[[1]], unique = TRUE)
    rows <- list()
    for (ln in sec[-1]) {
      ln <- trimws(ln); if (nchar(ln) == 0) next
      pt <- strsplit(ln, "\\s+")[[1]]
      if (length(pt) < length(hdr)) pt <- c(pt, rep(NA, length(hdr) - length(pt)))
      if (length(pt) > length(hdr)) pt <- pt[1:length(hdr)]
      rows[[length(rows) + 1]] <- pt
    }
    if (length(rows) == 0) return(data.frame())
    df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
    colnames(df) <- hdr
    df
  }

  # --- 1D parsing ---
  df1 <- data.frame()
  for (tb in trait_blocks) {
    sub <- lines[tb$start:tb$end]
    for (sec in extract(sub, "_1D_effect", "_1D_heritability")) {
      d <- parse_sec(sec)
      if (nrow(d) > 0) { d$TRAIT <- tb$name; df1 <- dplyr::bind_rows(df1, d) }
    }
  }

  if (nrow(df1) == 0) {
    qtl_data <- data.frame(TRAIT = character(), QTL = character(), A = character(),
                           SE = character(), P_Value = character(), stringsAsFactors = FALSE)
    qtl_dom <- data.frame(TRAIT = character(), QTL = character(), stringsAsFactors = FALSE)
  } else {
    qc <- grep("^(QTL|SNPID)$", colnames(df1), value = TRUE)[1]
    ac <- grep("^A$", colnames(df1), value = TRUE)[1]
    pc <- grep("^P\\.Value", colnames(df1), value = TRUE)[1]
    sc <- grep("^SE$", colnames(df1), value = TRUE)[1]

    vr <- !is.na(df1[[ac]]) & df1[[ac]] != "---"
    qtl_data <- df1[vr, c("TRAIT", qc, ac, sc, pc), drop = FALSE]
    colnames(qtl_data) <- c("TRAIT", "QTL", "A", "SE", "P_Value")

    dc <- grep("^D\\.", colnames(df1), value = TRUE)[1]
    if (is.na(dc)) dc <- grep("^D$", colnames(df1), value = TRUE)[1]
    if (!is.na(dc)) {
      vr <- !is.na(df1[[dc]]) & df1[[dc]] != "---"
      qtl_dom <- if (any(vr)) {
        d <- df1[vr, c("TRAIT", qc), drop = FALSE]
        colnames(d) <- c("TRAIT", "QTL"); d
      } else data.frame(TRAIT = character(), QTL = character(), stringsAsFactors = FALSE)
    } else {
      qtl_dom <- data.frame(TRAIT = character(), QTL = character(), stringsAsFactors = FALSE)
    }
  }

  res <- list(qtl_data = qtl_data, qtl_dom_data = qtl_dom)

  # --- 2D parsing ---
  if (scan_2d) {
    df2 <- data.frame()
    for (tb in trait_blocks) {
      sub <- lines[tb$start:tb$end]
      for (sec in extract(sub, "_2D_effect", "_2D_heritability")) {
        d <- parse_sec(sec)
        if (nrow(d) > 0) { d$TRAIT <- tb$name; df2 <- dplyr::bind_rows(df2, d) }
      }
    }

    if (nrow(df2) > 0) {
      q1 <- grep("^(QTL1|SNPID1)$", colnames(df2), value = TRUE)[1]
      q2 <- grep("^(QTL2|SNPID2)$", colnames(df2), value = TRUE)[1]
      aa <- grep("^AA", colnames(df2), value = TRUE)[1]
      if (!is.na(aa) && !is.na(q1) && !is.na(q2)) {
        vr <- !is.na(df2[[aa]]) & df2[[aa]] != "---"
        if (any(vr)) {
          d <- df2[vr, c("TRAIT", q1, q2), drop = FALSE]
          colnames(d) <- c("TRAIT", "QTL1", "QTL2")
          res$qtl_aa_data <- d
        } else {
          res$qtl_aa_data <- data.frame(TRAIT = character(), QTL1 = character(), QTL2 = character(), stringsAsFactors = FALSE)
        }
      } else {
        res$qtl_aa_data <- data.frame(TRAIT = character(), QTL1 = character(), QTL2 = character(), stringsAsFactors = FALSE)
      }
    } else {
      res$qtl_aa_data <- data.frame(TRAIT = character(), QTL1 = character(), QTL2 = character(), stringsAsFactors = FALSE)
    }
  }

  res
}


# ============================================================================
# 5. Layer 1 Screening
# ============================================================================

#' qtxnetwork.layer1.screen
#' Layer 1: classify SNPs by number of significant traits.
#' 
#' @param qtl_output Output from qtxnetwork.output.trans().
#' @param alpha Significance threshold. Default 0.05.
#' 
#' @return List: class1, layer2, snp_count, full_table.
#' @export
qtxnetwork.layer1.screen <- function(qtl_output, alpha = 0.05) {

  df <- qtl_output$qtl_data
  if (nrow(df) == 0) {
    return(list(class1 = character(), layer2 = character(),
                snp_count = data.frame(snp_id = character(), n_sig = integer()),
                full_table = df))
  }

  df$p_numeric <- suppressWarnings(as.numeric(df$P_Value))
  sig <- df[!is.na(df$p_numeric) & df$p_numeric < alpha, ]

  if (nrow(sig) == 0) {
    return(list(class1 = character(), layer2 = character(),
                snp_count = data.frame(snp_id = character(), n_sig = integer()),
                full_table = df))
  }

  cnt <- as.data.frame(table(sig$QTL), stringsAsFactors = FALSE)
  colnames(cnt) <- c("snp_id", "n_sig")
  cnt$n_sig <- as.integer(cnt$n_sig)

  list(
    class1 = cnt$snp_id[cnt$n_sig == 1],
    layer2 = cnt$snp_id[cnt$n_sig >= 2],
    snp_count = cnt,
    full_table = df
  )
}


# ============================================================================
# 6. Utility: %||% (from rlang, redefined to avoid dependency)
# ============================================================================
`%||%` <- function(x, y) if (is.null(x)) y else x