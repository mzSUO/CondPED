# Tests for Stage-2 subset machinery (v1.0):
# .enumerate_subsets(), .extract_minimum_representative_sets(),
# .extract_irreducible_modules(), .build_tolerance_path(),
# .check_loss_monotonicity().
#
# Structure anchors A-D (Runbook 4.7) use hand-built subset tables so
# the extraction logic is tested against exact expectations, not
# against losses computed by the code under test.

enum <- CondPED:::.enumerate_subsets
ext_rep <- CondPED:::.extract_minimum_representative_sets
ext_mod <- CondPED:::.extract_irreducible_modules
mono <- CondPED:::.check_loss_monotonicity

# Hand-build a per-locus subset table from key -> loss assignments.
# Unlisted subsets get loss_fill.
make_tab <- function(A, losses, loss_fill = 0) {
  sets <- unlist(lapply(0:length(A),
                        function(k) utils::combn(A, k, simplify = FALSE)),
                 recursive = FALSE)
  keys <- vapply(sets, CondPED:::.trait_set_key, character(1))
  data.frame(
    representing_set = I(sets),
    representing_key = keys,
    set_size = lengths(sets),
    representation_loss = ifelse(keys %in% names(losses),
                                 losses[keys], loss_fill),
    status = "ok",
    stringsAsFactors = FALSE
  )
}

# ---- enumeration ------------------------------------------------------------

test_that("all mode enumerates exactly 2^k subsets", {
  A <- c("A", "B", "C", "D")
  for (k in 1:4) {
    e <- enum(A[seq_len(k)], "all")
    expect_length(e$sets, 2^k)
    expect_length(unique(e$keys), 2^k)
  }
  e3 <- enum(c("A", "B", "C"), "all")
  expect_identical(e3$keys[1], "<empty>")
  expect_identical(e3$keys[length(e3$keys)], "A|B|C")
})

test_that("singleton mode: empty, full, singletons and complements, deduplicated", {
  e3 <- enum(c("A", "B", "C"), "singleton")
  expect_identical(sort(e3$keys),
                   sort(c("<empty>", "A|B|C", "A", "B", "C",
                          "B|C", "A|C", "A|B")))
  # k = 2: singleton complements coincide with singletons -> dedup to 4
  e2 <- enum(c("A", "B"), "singleton")
  expect_identical(sort(e2$keys), sort(c("<empty>", "A|B", "A", "B")))
  # k = 1: empty, full, singleton, singleton-complement(=empty) -> dedup
  # leaves exactly the empty set and {A}
  e1 <- enum("A", "singleton")
  expect_identical(sort(e1$keys), c("<empty>", "A"))
})

test_that("custom mode adds empty and full sets and validates membership", {
  e <- enum(c("A", "B", "C"), "custom",
            custom_sets = list(c("A", "B"), "C"))
  expect_identical(sort(e$keys), sort(c("<empty>", "A|B", "C", "A|B|C")))
  # custom set equal to the full set deduplicates
  e2 <- enum(c("A", "B"), "custom", custom_sets = list(c("A", "B")))
  expect_length(e2$sets, 2L)
  expect_error(enum(c("A", "B"), "custom", custom_sets = list("Z")),
               class = "condped_invalid_input")
  expect_error(enum(c("A", "B"), "custom", custom_sets = NULL),
               class = "condped_invalid_input")
})

test_that("enumeration covers the same family regardless of trait order", {
  e1 <- enum(c("A", "B", "C"), "all")
  e2 <- enum(c("C", "A", "B"), "all")
  # same family of sets (order-insensitive comparison)
  fam <- function(e) sort(vapply(e$sets, function(s) {
    paste(sort(s), collapse = "|")
  }, character(1)))
  expect_identical(fam(e1), fam(e2))
  # within one call, keys follow that call's reference order
  expect_true(all(e2$keys == vapply(e2$sets, CondPED:::.trait_set_key,
                                    character(1))))
})

# ---- structure anchors (Runbook 4.7) ------------------------------------------

test_that("anchor A: the unique minimum feasible set is returned alone", {
  tab <- make_tab(c("A", "B", "C"),
                  c("<empty>" = 1, A = 0.3, B = 0.4, C = 0.5,
                    `A|B` = 0.08, `A|C` = 0.3, `B|C` = 0.2, `A|B|C` = 0))
  rep <- ext_rep(tab, 0.10)
  expect_identical(nrow(rep), 1L)
  expect_identical(rep$trait_key, "A|B")
  expect_identical(rep$n_tied_solutions, 1L)
  expect_identical(rep$set_size, 2L)
})

test_that("anchor B: tied minimum sets are all returned", {
  tab <- make_tab(c("A", "B", "C"),
                  c(`<empty>` = 1, A = 0.3, B = 0.3, C = 0.3,
                    `A|B` = 0.05, `A|C` = 0.05, `B|C` = 0.3, `A|B|C` = 0))
  rep <- ext_rep(tab, 0.10)
  expect_identical(nrow(rep), 2L)
  expect_identical(sort(rep$trait_key), c("A|B", "A|C"))
  expect_identical(rep$n_tied_solutions, c(2L, 2L))
})

test_that("anchor C: a joint irreducible module is not split", {
  # Infeasible sets: empty, {C}, {D}, {C,D}; anything containing A or B
  # is feasible. The unique MAXIMAL infeasible set is {C,D}, so the
  # unique module is its complement {A,B}: removing A and B jointly
  # loses too much, but deleting either alone (complements {B,C,D},
  # {A,C,D}, loss 0.05) is fine, so the module is not split.
  tab <- make_tab(c("A", "B", "C", "D"),
                  c(`<empty>` = 1, C = 0.5, D = 0.5, `C|D` = 0.5),
                  loss_fill = 0.05)
  mod <- ext_mod(tab, 0.10)
  expect_identical(nrow(mod), 1L)
  expect_identical(mod$trait_key, "A|B")
  expect_identical(mod$complement_key, "C|D")
  expect_equal(mod$complement_loss, 0.5)
})

test_that("anchor D: multiple non-containing modules are all returned", {
  tab <- make_tab(c("A", "B", "C"),
                  c(`<empty>` = 1, A = 0.5, B = 0.05, C = 0.05,
                    `A|B` = 0.02, `A|C` = 0.02, `B|C` = 0.5, `A|B|C` = 0))
  mod <- ext_mod(tab, 0.10)
  expect_identical(sort(mod$trait_key), c("A", "B|C"))
  expect_equal(mod$complement_loss[mod$trait_key == "A"], 0.5)      # loss(B|C)
  expect_equal(mod$complement_loss[mod$trait_key == "B|C"], 0.5)    # loss(A)
})

test_that("full set as the unique representative; singleton module", {
  tab <- make_tab(c("A", "B", "C"), c(`A|B|C` = 0), loss_fill = 0.5)
  rep <- ext_rep(tab, 0.10)
  expect_identical(rep$trait_key, "A|B|C")
  expect_identical(rep$set_size, 3L)
  mod <- ext_mod(tab, 0.10)
  # every proper subset is infeasible -> every singleton's complement has
  # loss 0.5 > tol, so each singleton is a module
  expect_identical(sort(mod$trait_key), c("A", "B", "C"))

  tab2 <- make_tab(c("A", "B", "C"),
                   c(`<empty>` = 1, `B|C` = 0.5, A = 0.02,
                     `A|B` = 0.01, `A|C` = 0.01, `A|B|C` = 0),
                   loss_fill = 0.4)
  mod2 <- ext_mod(tab2, 0.10)
  expect_identical(mod2$trait_key, "A")   # only {A} is irreducible
})

test_that("feasible sets are upward closed on monotone tables", {
  tab <- make_tab(c("A", "B", "C"),
                  c(`<empty>` = 1, A = 0.3, B = 0.5, C = 0.6,
                    `A|B` = 0.08, `A|C` = 0.2, `B|C` = 0.3, `A|B|C` = 0))
  rep <- ext_rep(tab, 0.10)
  expect_identical(rep$trait_key, "A|B")
  # every superset of a feasible set is feasible on this table
  feas <- tab[tab$representation_loss <= 0.10, ]
  for (i in seq_len(nrow(feas))) {
    S <- feas$representing_set[[i]]
    supersets <- tab[vapply(tab$representing_set,
                            function(T) all(S %in% T), logical(1)), ]
    expect_true(all(supersets$representation_loss <= 0.10))
  }
})

test_that("modules equal complements of maximal infeasible sets (cross-check)", {
  tab <- make_tab(c("A", "B", "C", "D"),
                  c(`<empty>` = 1, A = 0.6, B = 0.05, C = 0.7, D = 0.05,
                    `A|B` = 0.02, `A|C` = 0.55, `A|D` = 0.02, `B|C` = 0.03,
                    `B|D` = 0.01, `C|D` = 0.5, `A|B|C` = 0.01,
                    `A|B|D` = 0.01, `A|C|D` = 0.5, `B|C|D` = 0.02,
                    `A|B|C|D` = 0))
  mod <- ext_mod(tab, 0.10)
  # independent algorithm: complements of inclusion-maximal infeasible sets
  infeas <- tab[tab$representation_loss > 0.10, ]
  maximal <- infeas[vapply(seq_len(nrow(infeas)), function(i) {
    Si <- infeas$representing_set[[i]]
    !any(vapply(seq_len(nrow(infeas)), function(j) {
      if (i == j) return(FALSE)
      all(Si %in% infeas$representing_set[[j]])
    }, logical(1)))
  }, logical(1)), ]
  expected <- sort(vapply(
    seq_len(nrow(maximal)),
    function(i) {
      CondPED:::.trait_set_key(setdiff(c("A", "B", "C", "D"),
                                       maximal$representing_set[[i]]))
    }, character(1)))
  expect_identical(sort(mod$trait_key), expected)
})

# ---- monotonicity checker ------------------------------------------------------

test_that("monotonicity checker detects violations and passes clean tables", {
  clean <- make_tab(c("A", "B"),
                    c(`<empty>` = 1, A = 0.4, B = 0.5, `A|B` = 0.1))
  m1 <- mono(clean)
  expect_identical(m1$n_violations, 0L)
  expect_equal(m1$max_violation, 0)

  dirty <- make_tab(c("A", "B"),
                    c(`<empty>` = 1, A = 0.05, B = 0.5, `A|B` = 0.4))
  m2 <- mono(dirty)
  expect_gt(m2$n_violations, 0L)
  expect_equal(m2$max_violation, 0.4 - 0.05, tolerance = 1e-12)

  # NA rows are skipped
  na_tab <- clean
  na_tab$representation_loss[na_tab$representing_key == "A"] <- NA
  expect_identical(mono(na_tab)$n_violations, 0L)
})

# ---- tolerance path --------------------------------------------------------------

test_that("tolerance path summarises counts per tolerance", {
  tab <- make_tab(c("A", "B", "C"),
                  c(`<empty>` = 1, A = 0.3, B = 0.3, C = 0.3,
                    `A|B` = 0.05, `A|C` = 0.05, `B|C` = 0.3, `A|B|C` = 0))
  tols <- c(0.05, 0.10)
  mr <- mo <- setNames(vector("list", 2), format(tols))
  for (t in tols) {
    mr[[format(t)]] <- ext_rep(tab, t)
    mo[[format(t)]] <- ext_mod(tab, t)
  }
  path <- CondPED:::.build_tolerance_path(tols, mr, mo)
  expect_identical(nrow(path), 2L)
  expect_identical(path$minimum_set_size, c(2L, 2L))
  expect_identical(path$n_minimum_sets, c(2L, 2L))
})

test_that("no feasible set yields an empty result, not an error", {
  tab <- make_tab(c("A", "B"), c(`<empty>` = 1, A = 0.5, B = 0.5,
                                 `A|B` = 0.5))
  rep <- ext_rep(tab, 0.10)
  expect_identical(nrow(rep), 0L)
  # loss(B) = 0.5 > tol -> {A} irreducible; loss(A) likewise -> {B}.
  # Both singletons are modules, so {A,B} is not minimal.
  mod <- ext_mod(tab, 0.10)
  expect_identical(sort(mod$trait_key), c("A", "B"))
})
