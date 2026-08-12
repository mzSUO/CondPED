#' Define associated loci from selected markers
#'
#' Organises an explicitly supplied set of selected markers into
#' associated loci. Selected markers only CONSTRUCT loci and elect
#' lead SNPs; once a locus boundary is defined, membership includes
#' every QC-passed marker inside the region, including markers that
#' were not significant in the marginal scan. This function never
#' performs genome-wide selection itself.
#'
#' @details
#' Methods:
#' \describe{
#'   \item{`physical`}{Per chromosome, each selected marker opens the
#'     interval `position +/- window_bp`; overlapping intervals are
#'     merged when `merge_overlaps = TRUE`.}
#'   \item{`ld_clump`}{Selected markers are processed in ascending
#'     marginal p-value order: the best remaining marker becomes the
#'     lead of a new clump; every remaining selected marker within
#'     `window_bp` of the lead AND with squared genotype correlation
#'     `r^2 >= r2_threshold` to the lead joins that clump. The locus
#'     region is `lead position +/- window_bp`.}
#'   \item{`custom`}{`locus_map` (a data.frame with `marker_id` and
#'     `locus_id`) defines membership directly; it must cover every
#'     selected marker.}
#' }
#'
#' For `physical` and `ld_clump`, membership is then every non-filtered
#' (QC-passed) marker in the omnibus table inside the region. The lead
#' SNP is the selected marker with the smallest marginal omnibus
#' p-value in the locus. Results are invariant to the input order of
#' markers.
#'
#' @param omnibus Omnibus scan result returned by [scan_mt_omnibus()],
#'   or its `omnibus` data.frame (columns `marker_id`, `p_value`,
#'   `status`).
#' @param G Numeric `n x p` genotype dosage matrix, used for
#'   `r^2` computations.
#' @param chromosome Length-`p` chromosome annotation per marker.
#' @param position Length-`p` physical position per marker (bp).
#' @param selected_markers Character vector of marker ids selected by
#'   the genome-wide screen. Required; no default.
#' @param marker_ids Marker order of `chromosome`/`position`/columns
#'   of `G`; defaults to `colnames(G)`.
#' @param method Locus construction method.
#' @param window_bp Window size in bp (required for `physical` and
#'   `ld_clump`).
#' @param r2_threshold Squared-correlation threshold for `ld_clump`.
#' @param locus_map data.frame with columns `marker_id`, `locus_id`
#'   (required for `custom`).
#' @param merge_overlaps Logical; merge overlapping intervals
#'   (`physical` method).
#'
#' @return A list with components
#' \describe{
#'   \item{loci}{data.frame with columns `locus_id` (deterministic
#'     `chromosome:start-end`, built from member positions),
#'     `chromosome`, `start`, `end`, `lead_snp`, `lead_p`,
#'     `n_significant_markers`, `n_region_markers`, `status`.}
#'   \item{membership}{data.frame with columns `locus_id`,
#'     `marker_id`, `significant_in_marginal_scan`, `r2_to_lead`,
#'     `position`.}
#'   \item{status}{Standard status list (`ok` or `empty_selection`).}
#'   \item{diagnostics}{List with `method`, `n_selected`,
#'     `n_loci`, `n_members` and parameter echoes.}
#' }
#' @export
define_associated_loci <- function(
  omnibus,
  G,
  chromosome,
  position,
  selected_markers,
  marker_ids = colnames(G),
  method = c("ld_clump", "physical", "custom"),
  window_bp = NULL,
  r2_threshold = NULL,
  locus_map = NULL,
  merge_overlaps = TRUE
) {
  method <- match.arg(method)
  if (missing(selected_markers) || is.null(selected_markers)) {
    .stop_invalid_input(
      "selected_markers must be supplied explicitly; define_associated_loci() never performs genome-wide selection."
    )
  }
  if (!is.character(selected_markers) || anyNA(selected_markers)) {
    .stop_invalid_input("selected_markers must be a character vector of marker ids.")
  }
  selected_markers <- unique(selected_markers)

  # ---- marker annotation ---------------------------------------------------
  if (is.null(marker_ids)) marker_ids <- paste0("marker", seq_len(ncol(G)))
  p <- length(marker_ids)
  if (anyDuplicated(marker_ids)) {
    .stop_invalid_input("marker_ids must be unique.")
  }
  if (length(chromosome) != p || anyNA(chromosome)) {
    .stop_invalid_input("chromosome must have one non-missing value per marker.")
  }
  if (!is.numeric(position) || length(position) != p || anyNA(position)) {
    .stop_invalid_input("position must be a numeric value per marker.")
  }
  names(chromosome) <- names(position) <- marker_ids
  unknown <- setdiff(selected_markers, marker_ids)
  if (length(unknown) > 0L) {
    .stop_invalid_input(
      "selected_markers not found among marker_ids: %s.",
      paste(utils::head(unknown, 3L), collapse = ", ")
    )
  }

  # ---- method requirements ----------------------------------------------------
  if (method %in% c("physical", "ld_clump")) {
    if (is.null(window_bp) || !is.numeric(window_bp) ||
        length(window_bp) != 1L || !is.finite(window_bp) ||
        window_bp <= 0) {
      .stop_invalid_input(
        "window_bp must be a single positive number for method = \"%s\".",
        method
      )
    }
  }
  if (method == "ld_clump") {
    if (is.null(r2_threshold) || !is.numeric(r2_threshold) ||
        length(r2_threshold) != 1L || !is.finite(r2_threshold) ||
        r2_threshold < 0 || r2_threshold > 1) {
      .stop_invalid_input(
        "r2_threshold must be a single number in [0, 1] for method = \"ld_clump\"."
      )
    }
  }
  if (method == "custom") {
    if (is.null(locus_map) || !is.data.frame(locus_map) ||
        !all(c("marker_id", "locus_id") %in% names(locus_map))) {
      .stop_invalid_input(
        "locus_map (data.frame with marker_id and locus_id) is required for method = \"custom\"."
      )
    }
    locus_map$marker_id <- as.character(locus_map$marker_id)
    uncovered <- setdiff(selected_markers, locus_map$marker_id)
    if (length(uncovered) > 0L) {
      .stop_invalid_input(
        "locus_map does not cover every selected marker (e.g. \"%s\").",
        uncovered[1L]
      )
    }
  }
  if (!is.logical(merge_overlaps) || length(merge_overlaps) != 1L ||
      is.na(merge_overlaps)) {
    .stop_invalid_input("merge_overlaps must be TRUE or FALSE.")
  }

  om <- .as_omnibus_table(omnibus)
  if (!"status" %in% names(om)) om$status <- "ok"
  p_of <- stats::setNames(om$p_value, om$marker_id)
  qc_markers <- om$marker_id[om$status != "filtered"]
  p_of_selected <- p_of[selected_markers]
  if (anyNA(p_of_selected)) {
    .stop_invalid_input(
      "selected_markers must have non-NA marginal omnibus p-values (e.g. \"%s\").",
      selected_markers[which(is.na(p_of_selected))[1L]]
    )
  }

  empty_loci <- data.frame(
    locus_id = character(), chromosome = character(),
    start = numeric(), end = numeric(), lead_snp = character(),
    lead_p = numeric(), n_significant_markers = integer(),
    n_region_markers = integer(), status = character(),
    stringsAsFactors = FALSE
  )
  empty_membership <- data.frame(
    locus_id = character(), marker_id = character(),
    significant_in_marginal_scan = logical(), r2_to_lead = numeric(),
    position = numeric(), stringsAsFactors = FALSE
  )

  if (length(selected_markers) == 0L) {
    return(list(
      loci = empty_loci,
      membership = empty_membership,
      status = .new_status(
        ok = TRUE, code = "empty_selection",
        message = "No selected markers; returning empty locus structures."
      ),
      diagnostics = list(
        method = method, n_selected = 0L, n_loci = 0L, n_members = 0L,
        window_bp = window_bp, r2_threshold = r2_threshold
      )
    ))
  }

  # ---- region construction ------------------------------------------------------
  # intervals: data.frame(chromosome, lo, hi, selected (list of marker ids))
  if (method == "custom") {
    locus_map <- locus_map[locus_map$marker_id %in% marker_ids, ,
                           drop = FALSE]
    lids <- unique(as.character(locus_map$locus_id))
    intervals <- lapply(lids, function(lid) {
      mk <- locus_map$marker_id[locus_map$locus_id == lid]
      list(locus_id = lid,
           chromosome = unique(as.character(chromosome[mk])),
           lo = min(position[mk]), hi = max(position[mk]),
           selected = intersect(selected_markers, mk),
           members = mk)
    })
    intervals <- intervals[lengths(lapply(intervals, `[[`, "selected")) > 0L]
    intervals <- lapply(intervals, function(iv) {
      if (length(iv$chromosome) != 1L) {
        .stop_invalid_input(
          "custom locus \"%s\" spans multiple chromosomes.", iv$locus_id
        )
      }
      iv
    })
  } else if (method == "physical") {
    sel <- data.frame(
      marker_id = selected_markers,
      chromosome = as.character(chromosome[selected_markers]),
      position = position[selected_markers],
      lo = position[selected_markers] - window_bp,
      hi = position[selected_markers] + window_bp,
      stringsAsFactors = FALSE
    )
    sel <- sel[order(sel$chromosome, sel$position), , drop = FALSE]
    intervals <- list()
    for (chr in unique(sel$chromosome)) {
      s <- sel[sel$chromosome == chr, , drop = FALSE]
      cur <- NULL
      for (i in seq_len(nrow(s))) {
        if (is.null(cur)) {
          cur <- list(locus_id = NULL, chromosome = chr,
                      lo = s$lo[i], hi = s$hi[i],
                      selected = s$marker_id[i], members = NULL)
          next
        }
        if (isTRUE(merge_overlaps) && s$lo[i] <= cur$hi) {
          cur$hi <- max(cur$hi, s$hi[i])
          cur$selected <- c(cur$selected, s$marker_id[i])
        } else {
          intervals <- c(intervals, list(cur))
          cur <- list(locus_id = NULL, chromosome = chr,
                      lo = s$lo[i], hi = s$hi[i],
                      selected = s$marker_id[i], members = NULL)
        }
      }
      if (!is.null(cur)) intervals <- c(intervals, list(cur))
    }
  } else { # ld_clump
    sel <- data.frame(
      marker_id = selected_markers,
      chromosome = as.character(chromosome[selected_markers]),
      position = position[selected_markers],
      p = unname(p_of_selected),
      stringsAsFactors = FALSE
    )
    sel <- sel[order(sel$p, sel$marker_id), , drop = FALSE]
    intervals <- list()
    pool <- seq_len(nrow(sel))
    k <- 0L
    while (length(pool) > 0L) {
      k <- k + 1L
      lead_i <- pool[1L]
      lead_chr <- sel$chromosome[lead_i]
      lead_pos <- sel$position[lead_i]
      in_clump <- lead_i
      rest <- pool[-1L]
      for (i in rest) {
        if (sel$chromosome[i] != lead_chr) next
        if (abs(sel$position[i] - lead_pos) > window_bp) next
        r2 <- .r2_between(G, marker_ids, sel$marker_id[lead_i],
                          sel$marker_id[i])
        if (is.finite(r2) && r2 >= r2_threshold) {
          in_clump <- c(in_clump, i)
        }
      }
      intervals <- c(intervals, list(list(
        locus_id = NULL, chromosome = lead_chr,
        lo = lead_pos - window_bp, hi = lead_pos + window_bp,
        selected = sel$marker_id[in_clump], members = NULL
      )))
      pool <- setdiff(pool, in_clump)
    }
    if (isTRUE(merge_overlaps)) {
      # clumping partitions selected markers by LD, but the resulting
      # search windows may overlap; in the formal analysis they are
      # merged into one locus
      intervals <- .merge_interval_list(intervals)
    }
  }

  # ---- membership, lead SNP, assembly -----------------------------------------
  r2_to_lead <- function(lead, member) {
    if (lead == member) return(1)
    .r2_between(G, marker_ids, lead, member)
  }
  loci_rows <- list()
  member_rows <- list()
  for (li in seq_along(intervals)) {
    iv <- intervals[[li]]
    chr <- iv$chromosome
    if (method == "custom") {
      members <- iv$members
      members <- intersect(members, qc_markers)
    } else {
      in_region <- qc_markers[
        as.character(chromosome[qc_markers]) == chr &
          position[qc_markers] >= iv$lo & position[qc_markers] <= iv$hi
      ]
      members <- union(iv$selected, in_region)
    }
    members <- members[order(position[members])]
    lead <- iv$selected[which.min(p_of[iv$selected])]
    start <- min(position[members])
    end <- max(position[members])
    # self-describing, deterministic locus id carrying the genomic
    # coordinates (truth/estimated locus MATCHING still uses
    # chromosome + overlap, never string equality)
    locus_id <- paste0(
      chr, ":",
      format(start, scientific = FALSE, trim = TRUE), "-",
      format(end, scientific = FALSE, trim = TRUE)
    )
    loci_rows[[li]] <- data.frame(
      locus_id = locus_id,
      chromosome = chr,
      start = start, end = end,
      lead_snp = lead,
      lead_p = unname(p_of[lead]),
      n_significant_markers = length(iv$selected),
      n_region_markers = length(members),
      status = "ok",
      stringsAsFactors = FALSE
    )
    member_rows[[li]] <- data.frame(
      locus_id = locus_id,
      marker_id = members,
      significant_in_marginal_scan = members %in% selected_markers,
      r2_to_lead = vapply(members, function(mk) r2_to_lead(lead, mk),
                          numeric(1)),
      position = unname(position[members]),
      stringsAsFactors = FALSE
    )
  }

  loci <- do.call(rbind, loci_rows)
  # identical spans (possible when overlapping search regions are
  # intentionally kept) get a deterministic suffix
  loci$locus_id <- make.unique(loci$locus_id, sep = "#")
  member_rows <- Map(function(df, id) {
    df$locus_id <- id
    df
  }, member_rows, loci$locus_id)

  list(
    loci = loci,
    membership = do.call(rbind, member_rows),
    status = .new_status(ok = TRUE, code = "ok"),
    diagnostics = list(
      method = method,
      n_selected = length(selected_markers),
      n_loci = length(intervals),
      n_members = sum(vapply(member_rows, nrow, integer(1))),
      window_bp = window_bp,
      r2_threshold = r2_threshold,
      merge_overlaps = merge_overlaps
    )
  )
}

#' Squared genotype correlation between two markers
#'
#' @param G n x p genotype matrix.
#' @param marker_ids Column order of `G`.
#' @param a,b Marker ids.
#' @return Numeric r^2, or `NA_real_` when undefined (zero variance or
#'   missing marker).
#' @keywords internal
.r2_between <- function(G, marker_ids, a, b) {
  ia <- match(a, marker_ids)
  ib <- match(b, marker_ids)
  if (is.na(ia) || is.na(ib)) return(NA_real_)
  xa <- G[, ia]
  xb <- G[, ib]
  keep <- !is.na(xa) & !is.na(xb)
  if (sum(keep) < 3L) return(NA_real_)
  xa <- xa[keep]
  xb <- xb[keep]
  if (stats::var(xa) <= 0 || stats::var(xb) <= 0) return(NA_real_)
  stats::cor(xa, xb)^2
}

#' Merge overlapping same-chromosome intervals
#'
#' @param intervals List of interval lists with `chromosome`, `lo`,
#'   `hi`, `selected`.
#' @return Merged list (selected markers unioned).
#' @keywords internal
.merge_interval_list <- function(intervals) {
  if (length(intervals) <= 1L) return(intervals)
  ord <- order(vapply(intervals, `[[`, character(1), "chromosome"),
               vapply(intervals, `[[`, numeric(1), "lo"))
  intervals <- intervals[ord]
  out <- list()
  cur <- intervals[[1L]]
  for (iv in intervals[-1L]) {
    if (iv$chromosome == cur$chromosome && iv$lo <= cur$hi) {
      cur$hi <- max(cur$hi, iv$hi)
      cur$selected <- union(cur$selected, iv$selected)
    } else {
      out <- c(out, list(cur))
      cur <- iv
    }
  }
  c(out, list(cur))
}
