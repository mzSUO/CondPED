# Tests for define_associated_loci() (Stage 6B-1).

dal <- define_associated_loci

mk_inputs <- function() {
  set.seed(1)
  n <- 60
  p <- 8
  G <- matrix(stats::rbinom(n * p, 2, 0.3), n, p,
              dimnames = list(NULL, paste0("M", 1:p)))
  chromosome <- rep("chr1", p)
  position <- c(1e5, 1.1e5, 1.2e5, 5e6, 5.05e6, 9e6, 9.02e6, 1e7)
  om <- data.frame(
    marker_id = paste0("M", 1:8),
    p_value = c(0.40, 1e-8, 0.30, 2e-6, 0.50, 0.70, 0.80, 0.001),
    status = c("ok", "ok", "ok", "ok", "ok", "filtered", "ok", "ok"),
    stringsAsFactors = FALSE
  )
  list(G = G, chromosome = chromosome, position = position, om = om)
}

test_that("missing selected_markers is invalid input", {
  x <- mk_inputs()
  expect_error(
    dal(x$om, x$G, x$chromosome, x$position,
        method = "physical", window_bp = 1e6),
    class = "condped_invalid_input"
  )
})

test_that("one peak gives one locus with a complete region membership", {
  x <- mk_inputs()
  res <- dal(x$om, x$G, x$chromosome, x$position,
             selected_markers = "M2",
             method = "physical", window_bp = 1e6)
  expect_identical(nrow(res$loci), 1L)
  expect_identical(res$loci$lead_snp, "M2")
  expect_identical(res$loci$locus_id, "chr1:100000-120000")
  expect_equal(res$loci$lead_p, 1e-8)
  # region = M2 +/- 1e6 -> covers M1, M2, M3 (not M6: filtered)
  expect_identical(res$membership$marker_id, c("M1", "M2", "M3"))
  expect_identical(res$membership$significant_in_marginal_scan,
                   c(FALSE, TRUE, FALSE))
  expect_identical(res$loci$n_significant_markers, 1L)
  expect_identical(res$loci$n_region_markers, 3L)
  # marginal non-significant QC markers are kept; filtered ones are not
  expect_true("M1" %in% res$membership$marker_id)
  expect_false("M6" %in% res$membership$marker_id)
})

test_that("two distant peaks give two loci", {
  x <- mk_inputs()
  res <- dal(x$om, x$G, x$chromosome, x$position,
             selected_markers = c("M2", "M8"),
             method = "physical", window_bp = 1e6)
  expect_identical(nrow(res$loci), 2L)
  expect_identical(res$loci$lead_snp, c("M2", "M8"))
})

test_that("merge_overlaps controls interval merging", {
  x <- mk_inputs()
  merged <- dal(x$om, x$G, x$chromosome, x$position,
                selected_markers = c("M2", "M3"),
                method = "physical", window_bp = 1e6,
                merge_overlaps = TRUE)
  expect_identical(nrow(merged$loci), 1L)
  unmerged <- dal(x$om, x$G, x$chromosome, x$position,
                  selected_markers = c("M2", "M3"),
                  method = "physical", window_bp = 1e6,
                  merge_overlaps = FALSE)
  expect_identical(nrow(unmerged$loci), 2L)
})

test_that("method-specific required arguments are enforced", {
  x <- mk_inputs()
  expect_error(dal(x$om, x$G, x$chromosome, x$position,
                   selected_markers = "M2", method = "physical"),
               class = "condped_invalid_input")
  expect_error(dal(x$om, x$G, x$chromosome, x$position,
                   selected_markers = "M2", method = "ld_clump"),
               class = "condped_invalid_input")
  expect_error(dal(x$om, x$G, x$chromosome, x$position,
                   selected_markers = "M2", method = "ld_clump",
                   window_bp = 1e6),
               class = "condped_invalid_input")
  expect_error(dal(x$om, x$G, x$chromosome, x$position,
                   selected_markers = "M2", method = "custom"),
               class = "condped_invalid_input")
})

test_that("lead SNP is the selected marker with the smallest marginal p", {
  x <- mk_inputs()
  # M2 (p = 1e-8) and M3 (p = 0.30 selected manually); M3 comes first
  # by position but M2 must lead
  x$om$p_value[x$om$marker_id == "M3"] <- 1e-4
  res <- dal(x$om, x$G, x$chromosome, x$position,
             selected_markers = c("M3", "M2"),
             method = "physical", window_bp = 1e6)
  expect_identical(nrow(res$loci), 1L)
  expect_identical(res$loci$lead_snp, "M2")
  expect_equal(res$loci$lead_p, 1e-8)
})

test_that("input marker order does not change locus definitions", {
  x <- mk_inputs()
  sel <- c("M2", "M8")
  base <- dal(x$om, x$G, x$chromosome, x$position,
              selected_markers = sel, method = "physical",
              window_bp = 1e6)
  perm <- sample.int(8)
  shuffled <- dal(
    x$om[perm, ], x$G[, perm], x$chromosome[perm], x$position[perm],
    selected_markers = sel, marker_ids = colnames(x$G)[perm],
    method = "physical", window_bp = 1e6
  )
  expect_identical(base$loci[, c("chromosome", "start", "end",
                                 "lead_snp", "n_significant_markers",
                                 "n_region_markers")],
                   shuffled$loci[, c("chromosome", "start", "end",
                                     "lead_snp", "n_significant_markers",
                                     "n_region_markers")])
  expect_identical(base$membership[order(base$membership$locus_id,
                                         base$membership$position), ],
                   shuffled$membership[order(shuffled$membership$locus_id,
                                             shuffled$membership$position), ],
                   ignore_attr = TRUE)
})

test_that("ld_clump groups by p order, window and r2", {
  n <- 200
  x1 <- stats::rbinom(n, 2, 0.3)
  x2 <- x1 + stats::rbinom(n, 1, 0.05)   # high r2 with x1
  x3 <- stats::rbinom(n, 2, 0.3)         # independent
  G <- cbind(A1 = x1, A2 = x2, B1 = x3)
  om <- data.frame(marker_id = c("A1", "A2", "B1"),
                   p_value = c(1e-6, 1e-4, 1e-5),
                   status = "ok", stringsAsFactors = FALSE)
  chr <- rep("chr1", 3)
  pos <- c(1e5, 1.2e5, 1.1e5)
  # merge_overlaps = FALSE: A1 leads; A2 joins (close + high r2);
  # B1 forms its own locus
  res <- dal(om, G, chr, pos, selected_markers = c("A1", "A2", "B1"),
             method = "ld_clump", window_bp = 5e5, r2_threshold = 0.5,
             merge_overlaps = FALSE)
  expect_identical(nrow(res$loci), 2L)
  expect_identical(res$loci$lead_snp[1], "A1")
  expect_identical(res$loci$n_significant_markers, c(2L, 1L))
  # the two overlapping search regions share regional QC members
  expect_true("B1" %in% res$membership$marker_id[
    res$membership$locus_id == res$loci$locus_id[1]])
  # identical spans get deterministic suffixes, never a collision
  expect_length(unique(res$loci$locus_id), 2L)
  # merge_overlaps = TRUE (formal default): overlapping regions merge
  res_m <- dal(om, G, chr, pos, selected_markers = c("A1", "A2", "B1"),
               method = "ld_clump", window_bp = 5e5, r2_threshold = 0.5)
  expect_identical(nrow(res_m$loci), 1L)
  expect_identical(res_m$loci$lead_snp, "A1")
  expect_identical(res_m$loci$n_significant_markers, 3L)
  # with a strict r2 threshold A2 also splits off
  res2 <- dal(om, G, chr, pos, selected_markers = c("A1", "A2", "B1"),
              method = "ld_clump", window_bp = 5e5,
              r2_threshold = 0.999, merge_overlaps = FALSE)
  expect_identical(nrow(res2$loci), 3L)
  # r2_to_lead is recorded and equals 1 for the lead itself
  expect_equal(res$membership$r2_to_lead[
    res$membership$locus_id == res$loci$locus_id[1] &
      res$membership$marker_id == "A1"], 1)
})

test_that("custom mode uses locus_map and validates coverage", {
  x <- mk_inputs()
  map <- data.frame(
    marker_id = c("M1", "M2", "M3", "M8"),
    locus_id = c("locA", "locA", "locA", "locB"),
    stringsAsFactors = FALSE
  )
  res <- dal(x$om, x$G, x$chromosome, x$position,
             selected_markers = c("M2", "M8"),
             method = "custom", locus_map = map)
  expect_identical(nrow(res$loci), 2L)
  expect_identical(res$membership$marker_id[
    res$membership$locus_id == res$loci$locus_id[1]], c("M1", "M2", "M3"))
  bad_map <- map[map$marker_id != "M8", ]
  expect_error(
    dal(x$om, x$G, x$chromosome, x$position,
        selected_markers = c("M2", "M8"),
        method = "custom", locus_map = bad_map),
    class = "condped_invalid_input"
  )
})

test_that("empty selected_markers returns a complete empty result", {
  x <- mk_inputs()
  res <- dal(x$om, x$G, x$chromosome, x$position,
             selected_markers = character(),
             method = "physical", window_bp = 1e6)
  expect_identical(res$status$code, "empty_selection")
  expect_true(res$status$ok)
  expect_identical(nrow(res$loci), 0L)
  expect_identical(nrow(res$membership), 0L)
  expect_identical(names(res$loci),
                   c("locus_id", "chromosome", "start", "end",
                     "lead_snp", "lead_p", "n_significant_markers",
                     "n_region_markers", "status"))
})

test_that("selected markers must exist and have marginal p-values", {
  x <- mk_inputs()
  expect_error(dal(x$om, x$G, x$chromosome, x$position,
                   selected_markers = "nope",
                   method = "physical", window_bp = 1e6),
               class = "condped_invalid_input")
  x$om$p_value[x$om$marker_id == "M2"] <- NA
  expect_error(dal(x$om, x$G, x$chromosome, x$position,
                   selected_markers = "M2",
                   method = "physical", window_bp = 1e6),
               class = "condped_invalid_input")
})
