# qtlnetwork.R
# QTLNetwork wrapper for RIL GWAS (additive + AA epistasis, no dominance)

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
  
  geno_mat <- as.matrix(geno_data[, -c(1:3)])
  rownames(geno_mat) <- geno_data[, 2]
  geno_t <- t(geno_mat)
  
  geno_df <- as.data.frame(geno_t, stringsAsFactors = FALSE)
  geno_df <- cbind(rownames(geno_df), geno_df)
  colnames(geno_df) <- c("_Genotype", colnames(geno_t))
  
  traitNum <- ncol(pheno_data) - 1
  n_ind <- nrow(pheno_data)
  
  # --- .gen file ---
  gen_file <- paste0(geno_output_prefix, ".gen")
  cat("Writing .gen:", gen_file, "\n")
  
  con_gen <- file(gen_file, "w")  # ← 独立变量名
  write.table(geno_df, file = con_gen, append = FALSE, na = ".",
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  close(con_gen)  # ← 显式关闭，无 on.exit
  
  # --- .phe file ---
  phe_file <- paste0(pheno_output_prefix, ".phe")
  cat("Writing", traitNum, "traits to .phe:", phe_file, "\n")
  
  con_phe <- file(phe_file, "w")  # ← 独立变量名，不会覆盖 con_gen
  
  writeLines(paste("_Population", population, sep = "\t"), con_phe)
  writeLines(paste("_Genotypes", n_ind, sep = "\t"), con_phe)
  writeLines(paste("_Observations", n_ind, sep = "\t"), con_phe)
  writeLines("_Environments\tno", con_phe)
  writeLines("_Replications\tno", con_phe)
  writeLines(paste("_TraitNumber", traitNum, paste(seq_len(traitNum), collapse = "+"), sep = "\t"), con_phe)
  writeLines(paste("_Chromosomes", length(chrName), paste(chrName, collapse = " "), sep = "\t"), con_phe)
  writeLines(paste("_TotalMarker", sum(chrCounts), paste(chrCounts, collapse = " "), sep = "\t"), con_phe)
  writeLines(ifelse(population == "RIL",
                    "_MarkerCode\tP1=1\tP2=-1",
                    "_MarkerCode\tP1=1\tP2=-1\tF1=0"), con_phe)
  writeLines("*TraitBegin*", con_phe)
  
  phe_out <- pheno_data
  colnames(phe_out)[1] <- "Geno#"
  write.table(phe_out, file = con_phe, append = TRUE, na = ".",
              row.names = FALSE, col.names = TRUE, quote = FALSE,
              sep = "\t", eol = ";\n")
  writeLines("*TraitEnd*", con_phe)
  
  close(con_phe)  # ← 显式关闭
  invisible(NULL)
}


# ============================================================================
# 3. Adapter: simulation.R -> QTXNetwork
# ============================================================================

#' condped_to_qtlnetwork
#' Bridge simulation.R output to QTXNetwork input.
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
  X <- sim_result$X
  Y <- sim_result$Y
  
  # 确保行列名完整
  if (is.null(colnames(Y))) colnames(Y) <- paste0("Trait", seq_len(ncol(Y)))
  if (is.null(rownames(Y))) rownames(Y) <- paste0("Ind", seq_len(nrow(Y)))
  if (is.null(rownames(X))) rownames(X) <- rownames(Y)
  stopifnot(identical(rownames(X), rownames(Y)))
  
  n_ind <- nrow(X)
  p <- ncol(X)
  m <- ncol(Y)
  
  ind_names <- rownames(Y)
  geno_data <- data.frame(
    chr = rep("1", p),
    snp_id = colnames(X) %||% paste0("SNP", seq_len(p)),
    pos = seq_len(p),
    t(X),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  colnames(geno_data)[4:ncol(geno_data)] <- ind_names
  
  trait_names <- colnames(Y)
  pheno_data <- data.frame(
    id = ind_names,
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
#' @param pheno_data Phenotype data frame or matrix.
#' @param pre_file Path to .pre file.
#' @param scan_2d Whether 2D scanning was performed. Default FALSE.
#' 
#' @return List: qtl_data, qtl_dom_data, qtl_aa_data (if scan_2d).
#' @export
qtxnetwork.output.trans <- function(pheno_data, pre_file, scan_2d = FALSE)
{
  # --- 兼容矩阵输入 ---
  if (is.matrix(pheno_data)) {
    if (is.null(rownames(pheno_data))) stop("pheno_data matrix must have row names")
    if (is.null(colnames(pheno_data))) colnames(pheno_data) <- paste0("Trait", seq_len(ncol(pheno_data)))
    pheno_data <- as.data.frame(pheno_data)
    pheno_data <- cbind(id = rownames(pheno_data), pheno_data)
    rownames(pheno_data) <- NULL
  }
  
  traitName <- colnames(pheno_data)[-c(1)]
  traitNum <- length(traitName)
  
  file_lines <- readLines(pre_file)
  
  # --- 先定位所有 _trait 标记，确定每个 trait 的文本范围 ---
  trait_idx <- grep("^_trait\\s+\\d+", file_lines)
  trait_ranges <- list()
  for (i in seq_along(trait_idx)) {
    parts <- strsplit(trimws(file_lines[trait_idx[i]]), "\\s+")[[1]]
    name <- paste(parts[-(1:2)], collapse = " ")
    start_line <- trait_idx[i]
    end_line <- if (i < length(trait_idx)) trait_idx[i + 1] - 1 else length(file_lines)
    trait_ranges[[i]] <- list(name = name, start = start_line, end = end_line)
  }
  
  if (length(trait_ranges) == 0) {
    trait_ranges <- list(list(name = traitName[1], start = 1, end = length(file_lines)))
  }
  
  # ---- 1D parsing ----
  start_title <- "_1D_effect"
  end_title <- "_1D_heritability"
  df.qtl <- data.frame()
  
  for (c in seq_along(trait_ranges)) {
    trait <- trait_ranges[[c]]$name
    block_lines <- file_lines[trait_ranges[[c]]$start:trait_ranges[[c]]$end]
    
    start_line <- 0
    end_line <- 0
    for (i in seq_along(block_lines)) {
      line <- block_lines[i]
      if (trimws(line) == start_title) {
        start_line <- i
        next
      }
      if (trimws(line) == end_title) {
        end_line <- i
        break
      }
    }
    
    if (start_line > 0 && end_line > start_line) {
      raw <- block_lines[(start_line + 1):(end_line - 1)]
      raw <- raw[nchar(trimws(raw)) > 0]
      
      if (length(raw) >= 2) {
        hdr <- make.names(strsplit(trimws(raw[1]), "\\s+")[[1]], unique = TRUE)
        rows <- list()
        for (ln in raw[-1]) {
          ln <- trimws(ln)
          if (nchar(ln) == 0) next
          pt <- strsplit(ln, "\\s+")[[1]]
          if (length(pt) < length(hdr)) pt <- c(pt, rep(NA, length(hdr) - length(pt)))
          if (length(pt) > length(hdr)) pt <- pt[1:length(hdr)]
          rows[[length(rows) + 1]] <- pt
        }
        if (length(rows) > 0) {
          tmp <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
          colnames(tmp) <- hdr
          tmp$TRAIT <- trait
          df.qtl <- dplyr::bind_rows(df.qtl, tmp)
        }
      }
    }
  }
  
  if (nrow(df.qtl) == 0) {
    return(list(
      qtl_data = data.frame(TRAIT = character(), QTL = character(), A = character(),
                            SE = character(), P_Value = character(), stringsAsFactors = FALSE),
      qtl_dom_data = data.frame(TRAIT = character(), QTL = character(), stringsAsFactors = FALSE)
    ))
  }
  
  # --- 容错列名匹配 ---
  qc <- grep("SNPID", colnames(df.qtl), value = TRUE, ignore.case = TRUE)[1]
  ac <- grep("^A$|Add|ADD|Effect", colnames(df.qtl), value = TRUE, ignore.case = TRUE)[1]
  pc <- grep("P.*Value|Pval|p.*value", colnames(df.qtl), value = TRUE, ignore.case = TRUE)[1]
  sc <- grep("^SE$|Std.*Err|StdErr", colnames(df.qtl), value = TRUE, ignore.case = TRUE)[1]
  
  if (is.na(qc)) stop("Cannot find QTL ID column. Available: ", paste(colnames(df.qtl), collapse = ", "))
  if (is.na(ac)) stop("Cannot find additive effect (A) column. Available: ", paste(colnames(df.qtl), collapse = ", "))
  if (is.na(pc)) stop("Cannot find P-value column. Available: ", paste(colnames(df.qtl), collapse = ", "))
  
  if (is.na(sc)) {
    warning("Cannot find SE column; filling with NA")
    df.qtl[["SE"]] <- NA_character_
    sc <- "SE"
  }
  
  vr <- !is.na(df.qtl[[ac]]) & df.qtl[[ac]] != "---"
  qtl_data <- df.qtl[vr, c("TRAIT", qc, ac, sc, pc), drop = FALSE]
  colnames(qtl_data) <- c("TRAIT", "SNPID", "A", "SE", "P_Value")
  
  dc <- grep("^D$|Dom|DOM|Dominance", colnames(df.qtl), value = TRUE, ignore.case = TRUE)[1]
  if (!is.na(dc)) {
    vr <- !is.na(df.qtl[[dc]]) & df.qtl[[dc]] != "---"
    qtl_dom_data <- if (any(vr)) {
      d <- df.qtl[vr, c("TRAIT", qc), drop = FALSE]
      colnames(d) <- c("TRAIT", "QTL")
      d
    } else {
      data.frame(TRAIT = character(), QTL = character(), stringsAsFactors = FALSE)
    }
  } else {
    qtl_dom_data <- data.frame(TRAIT = character(), QTL = character(), stringsAsFactors = FALSE)
  }
  
  if (nrow(qtl_data) > 0) rownames(qtl_data) <- seq_len(nrow(qtl_data))
  if (nrow(qtl_dom_data) > 0) rownames(qtl_dom_data) <- seq_len(nrow(qtl_dom_data))
  
  res <- list(qtl_data = qtl_data, qtl_dom_data = qtl_dom_data)
  
  # ---- 2D parsing ----
  if (scan_2d) {
    start_2d <- "_2D_effect"
    end_2d <- "_2D_heritability"
    df.2d <- data.frame()
    
    for (c in seq_along(trait_ranges)) {
      trait <- trait_ranges[[c]]$name
      block_lines <- file_lines[trait_ranges[[c]]$start:trait_ranges[[c]]$end]
      
      start_line <- 0
      end_line <- 0
      for (i in seq_along(block_lines)) {
        line <- block_lines[i]
        if (trimws(line) == start_2d) {
          start_line <- i
          next
        }
        if (trimws(line) == end_2d) {
          end_line <- i
          break
        }
      }
      
      if (start_line > 0 && end_line > start_line) {
        raw <- block_lines[(start_line + 1):(end_line - 1)]
        raw <- raw[nchar(trimws(raw)) > 0]
        if (length(raw) >= 2) {
          hdr <- make.names(strsplit(trimws(raw[1]), "\\s+")[[1]], unique = TRUE)
          rows <- list()
          for (ln in raw[-1]) {
            ln <- trimws(ln)
            if (nchar(ln) == 0) next
            pt <- strsplit(ln, "\\s+")[[1]]
            if (length(pt) < length(hdr)) pt <- c(pt, rep(NA, length(hdr) - length(pt)))
            if (length(pt) > length(hdr)) pt <- pt[1:length(hdr)]
            rows[[length(rows) + 1]] <- pt
          }
          if (length(rows) > 0) {
            tmp <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
            colnames(tmp) <- hdr
            tmp$TRAIT <- trait
            df.2d <- dplyr::bind_rows(df.2d, tmp)
          }
        }
      }
    }
    
    if (nrow(df.2d) > 0) {
      q1 <- grep("QTL1|SNPID1|SNP1", colnames(df.2d), value = TRUE, ignore.case = TRUE)[1]
      q2 <- grep("QTL2|SNPID2|SNP2", colnames(df.2d), value = TRUE, ignore.case = TRUE)[1]
      aa <- grep("AA|Epistasis|AA_effect", colnames(df.2d), value = TRUE, ignore.case = TRUE)[1]
      
      if (!is.na(aa) && !is.na(q1) && !is.na(q2)) {
        vr <- !is.na(df.2d[[aa]]) & df.2d[[aa]] != "---"
        if (any(vr)) {
          d <- df.2d[vr, c("TRAIT", q1, q2), drop = FALSE]
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
  
  cnt <- as.data.frame(table(sig$SNPID), stringsAsFactors = FALSE)
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
# 6. Utility
# ============================================================================
`%||%` <- function(x, y) if (is.null(x)) y else x