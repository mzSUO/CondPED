# tests/testthat/test-qtlnetwork.R
# Test suite for QTXNetwork wrapper functions

# -----------------------------------------------------------------------------
# qtxnetwork.input.trans
# -----------------------------------------------------------------------------

test_that("qtxnetwork.input.trans produces correct .gen header", {
  geno <- data.frame(
    chr = c("1", "1", "2"),
    snp_id = c("SNP1", "SNP2", "SNP3"),
    pos = c(1, 2, 1),
    Ind1 = c(1, -1, 1),
    Ind2 = c(-1, 1, -1),
    stringsAsFactors = FALSE
  )
  pheno <- data.frame(id = c("Ind1", "Ind2"), TraitA = c(1.5, 2.3), TraitB = c(0.5, 1.2))

  prefix <- tempfile()
  qtxnetwork.input.trans(geno, pheno, prefix, prefix, population = "RIL")

  # Check .gen file exists and has header
  gen_file <- paste0(prefix, ".gen")
  expect_true(file.exists(gen_file))
  lines <- readLines(gen_file)
  expect_true(any(grepl("_Chromosomes", lines)))
  expect_true(any(grepl("_TotalMarker", lines)))
  expect_true(any(grepl("_MarkerCode", lines)))
  expect_true(any(grepl("P1=1.*P2=-1", lines)))  # RIL: no F1=0
  expect_false(any(grepl("F1=0", lines[1:3])))    # Header should not have F1=0 for RIL

  # Check .phe file
  phe_file <- paste0(prefix, ".phe")
  expect_true(file.exists(phe_file))
  phe_lines <- readLines(phe_file)
  expect_true(any(grepl("_TraitNumber", phe_lines)))
  expect_true(any(grepl("TraitA.*TraitB", phe_lines)))
  expect_true(any(grepl("\\*TraitBegin\\*", phe_lines)))
  expect_true(any(grepl("\\*TraitEnd\\*", phe_lines)))

  # Cleanup
  unlink(c(gen_file, phe_file))
})

test_that("qtxnetwork.input.trans produces F2 header with F1=0", {
  geno <- data.frame(chr = "1", snp_id = "SNP1", pos = 1, Ind1 = 1, stringsAsFactors = FALSE)
  pheno <- data.frame(id = "Ind1", TraitA = 1.0, stringsAsFactors = FALSE)

  prefix <- tempfile()
  qtxnetwork.input.trans(geno, pheno, prefix, prefix, population = "F2")

  lines <- readLines(paste0(prefix, ".gen"))
  expect_true(any(grepl("F1=0", lines[1:3])))

  unlink(paste0(prefix, c(".gen", ".phe")))
})


# -----------------------------------------------------------------------------
# condped_to_qtlnetwork
# -----------------------------------------------------------------------------

test_that("condped_to_qtlnetwork handles RIL simulation output", {
  sim <- generate_class1(n = 10, p = 5, seed = 42)
  prefix <- tempfile()

  condped_to_qtlnetwork(sim, prefix, prefix, population = "RIL")

  gen_file <- paste0(prefix, ".gen")
  phe_file <- paste0(prefix, ".phe")
  expect_true(file.exists(gen_file))
  expect_true(file.exists(phe_file))

  # Verify geno data in .gen file
  gen_lines <- readLines(gen_file)
  expect_true(any(grepl("SNP1", gen_lines)))

  unlink(c(gen_file, phe_file))
})

test_that("condped_to_qtlnetwork handles missing colnames", {
  sim <- generate_class1(n = 10, p = 5, seed = 42)
  colnames(sim$X) <- NULL  # Remove colnames
  prefix <- tempfile()

  condped_to_qtlnetwork(sim, prefix, prefix, population = "RIL")

  gen_lines <- readLines(paste0(prefix, ".gen"))
  expect_true(any(grepl("SNP1", gen_lines)))  # Should auto-generate

  unlink(paste0(prefix, c(".gen", ".phe")))
})


# -----------------------------------------------------------------------------
# qtxnetwork.output.trans
# -----------------------------------------------------------------------------

# Create a mock .pre file for testing
make_mock_pre <- function(path, with_2d = FALSE) {
  lines <- c(
    "_trait 1 TraitA",
    "_variance_components ( population mean: 0.0 phenotypic variance: 1.0 Total Heritability: 0.5 )",
    "h^2(A)         h^2(I)         V(e)/V(P)      ",
    "0.5000         0.0000         0.5000         ",
    "",
    "_1D_effect",
    "QTL       SNPID            A           SE          P-Value     ",
    "1-1       SNP1             1.2345      0.1234      1.234e-10   ",
    "1-2       SNP2             0.9876      0.0987      2.345e-08   ",
    "",
    "_1D_heritability",
    "QTL              h^2(a)         ",
    "1-1              0.3000         ",
    "1-2              0.2000         ",
    "",
    "_trait 2 TraitB",
    "_variance_components ( population mean: 0.0 phenotypic variance: 1.0 Total Heritability: 0.5 )",
    "h^2(A)         h^2(I)         V(e)/V(P)      ",
    "0.5000         0.0000         0.5000         ",
    "",
    "_1D_effect",
    "QTL       SNPID            A           SE          P-Value     ",
    "1-1       SNP1             0.8765      0.0876      3.456e-09   ",
    "1-3       SNP3             0.7654      0.0765      4.567e-07   ",
    "",
    "_1D_heritability",
    "QTL              h^2(a)         ",
    "1-1              0.2500         ",
    "1-3              0.1500         "
  )

  if (with_2d) {
    lines <- c(lines, "",
      "_2D_effect",
      "QTL1      QTL2      SNPID1    SNPID2    AA          SE          P-Value     ",
      "1-1       1-2       SNP1      SNP2      0.5432      0.0543      5.678e-06   ",
      "",
      "_2D_heritability",
      "QTL1             QTL2             h^2(aa)        ",
      "1-1              1-2              0.1000         "
    )
  }

  writeLines(lines, path)
}

test_that("qtxnetwork.output.trans parses 1D results correctly", {
  pre_file <- tempfile()
  make_mock_pre(pre_file, with_2d = FALSE)

  pheno <- data.frame(id = 1:2, TraitA = c(1, 2), TraitB = c(3, 4))
  qtl <- qtxnetwork.output.trans(pheno, pre_file, scan_2d = FALSE)

  expect_equal(nrow(qtl$qtl_data), 4)  # 2 SNPs × 2 traits
  expect_equal(ncol(qtl$qtl_data), 5)   # TRAIT, QTL, A, SE, P_Value
  expect_equal(sort(unique(qtl$qtl_data$TRAIT)), c("TraitA", "TraitB"))
  expect_true(all(c("SNP1", "SNP2", "SNP3") %in% qtl$qtl_data$QTL))

  # SNP1 appears in both traits
  expect_equal(sum(qtl$qtl_data$QTL == "SNP1"), 2)

  # Dominance should be empty for RIL output
  expect_equal(nrow(qtl$qtl_dom_data), 0)

  unlink(pre_file)
})

test_that("qtxnetwork.output.trans parses 2D results when scan_2d=TRUE", {
  pre_file <- tempfile()
  make_mock_pre(pre_file, with_2d = TRUE)

  pheno <- data.frame(id = 1:2, TraitA = c(1, 2), TraitB = c(3, 4))
  qtl <- qtxnetwork.output.trans(pheno, pre_file, scan_2d = TRUE)

  expect_true("qtl_aa_data" %in% names(qtl))
  expect_equal(nrow(qtl$qtl_aa_data), 1)
  expect_equal(qtl$qtl_aa_data$QTL1[1], "SNP1")
  expect_equal(qtl$qtl_aa_data$QTL2[1], "SNP2")

  unlink(pre_file)
})

test_that("qtxnetwork.output.trans handles empty .pre file", {
  pre_file <- tempfile()
  writeLines("_trait 1 TraitA\n_1D_effect\n_1D_heritability", pre_file)

  pheno <- data.frame(id = 1, TraitA = 1)
  qtl <- qtxnetwork.output.trans(pheno, pre_file)

  expect_equal(nrow(qtl$qtl_data), 0)
  expect_equal(nrow(qtl$qtl_dom_data), 0)

  unlink(pre_file)
})

test_that("qtxnetwork.output.trans handles missing _trait lines", {
  pre_file <- tempfile()
  writeLines("_1D_effect\nQTL SNPID A\n1-1 SNP1 1.0\n_1D_heritability", pre_file)

  pheno <- data.frame(id = 1, TraitA = 1)
  qtl <- qtxnetwork.output.trans(pheno, pre_file)

  expect_equal(nrow(qtl$qtl_data), 1)
  expect_equal(qtl$qtl_data$TRAIT[1], "TraitA")  # Fallback to pheno column name

  unlink(pre_file)
})


# -----------------------------------------------------------------------------
# qtxnetwork.layer1.screen
# -----------------------------------------------------------------------------

test_that("layer1.screen counts significant traits per SNP", {
  qtl_data <- data.frame(
    TRAIT = c("TraitA", "TraitA", "TraitB", "TraitB", "TraitC"),
    QTL = c("SNP1", "SNP2", "SNP1", "SNP3", "SNP4"),
    A = c("1.0", "0.5", "0.8", "0.6", "0.3"),
    SE = c("0.1", "0.1", "0.1", "0.1", "0.1"),
    P_Value = c("1e-10", "1e-3", "1e-8", "1e-6", "0.1"),  # SNP4 not significant
    stringsAsFactors = FALSE
  )
  qtl <- list(qtl_data = qtl_data, qtl_dom_data = data.frame())

  screen <- qtxnetwork.layer1.screen(qtl, alpha = 0.05)

  expect_equal(sort(screen$class1), c("SNP2", "SNP3", "SNP4"))
  expect_equal(screen$layer2, "SNP1")  # SNP1 significant in 2 traits
  expect_equal(screen$snp_count$n_sig[screen$snp_count$snp_id == "SNP1"], 2)
})

test_that("layer1.screen handles empty qtl_data", {
  qtl <- list(qtl_data = data.frame(TRAIT = character(), QTL = character(), 
                                     A = character(), SE = character(), 
                                     P_Value = character(), stringsAsFactors = FALSE),
              qtl_dom_data = data.frame())
  screen <- qtxnetwork.layer1.screen(qtl)
  expect_equal(length(screen$class1), 0)
  expect_equal(length(screen$layer2), 0)
})

test_that("layer1.screen respects alpha threshold", {
  qtl_data <- data.frame(
    TRAIT = c("TraitA", "TraitB"),
    QTL = c("SNP1", "SNP1"),
    A = c("1.0", "0.5"),
    SE = c("0.1", "0.1"),
    P_Value = c("0.01", "0.1"),  # Only first is significant at alpha=0.05
    stringsAsFactors = FALSE
  )
  qtl <- list(qtl_data = qtl_data, qtl_dom_data = data.frame())

  screen <- qtxnetwork.layer1.screen(qtl, alpha = 0.05)
  expect_equal(screen$class1, "SNP1")  # Only 1 significant trait
  expect_equal(length(screen$layer2), 0)
})