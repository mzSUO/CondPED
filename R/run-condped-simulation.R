# Reproducible simulation runner (Stage 7).
#
# The runner only orchestrates: scenario grid -> canonical scenario id
# -> deterministic per-replicate seed -> simulate_condped_data() ->
# the frozen Stage 6B/6C analysis and comparison functions -> the
# frozen evaluator -> atomic per-replicate saves -> resume/parallel.
# It never recomputes truth, eta, rho, Rep/Irr, signal matching,
# ASSET statistics or any CondPED statistic.

#' Canonical scenario id from formal settings
#'
#' Only formal data-generation/analysis parameters participate.
#' Serialisation is order-invariant (sorted by name), so grid row
#' order and parameter order never change the id. Runtime, worker id
#' and output paths are excluded by construction.
#'
#' @param settings Named list of formal settings.
#' @return A single canonical string.
#' @keywords internal
.canonical_scenario_id <- function(settings) {
  drop <- c("master_seed", "rep_id", "workers", "out_dir", "backend",
            "resume", "overwrite", "save_data", "save_fit",
            "fail_policy", "progress")
  settings <- settings[setdiff(names(settings), drop)]
  settings <- settings[sort(names(settings))]
  fmt <- function(v) {
    if (is.null(v)) return("<null>")
    if (is.list(v)) {
      return(paste0("(", paste(vapply(v, fmt, character(1)),
                               collapse = ","), ")"))
    }
    if (length(v) > 1L) return(paste(v, collapse = ","))
    if (is.na(v)) return("<NA>")
    format(v, scientific = FALSE, trim = TRUE)
  }
  paste(sprintf("%s=%s", names(settings),
                vapply(settings, fmt, character(1))),
        collapse = "|")
}

#' 32-bit FNV-1a hash as a non-negative integer
#' @keywords internal
.fnv1a_int <- function(x) {
  bytes <- utf8ToInt(x)
  h <- 2166136261                    # offset basis (double, exact < 2^53)
  for (b in bytes) {
    # 32-bit XOR via signed-integer round trip
    h_signed <- if (h >= 2147483648) h - 4294967296 else h
    x <- bitwXor(as.integer(h_signed), as.integer(b))
    h <- if (x < 0) x + 4294967296 else as.numeric(x)
    h <- (h * 16777619) %% 4294967296
  }
  h
}

#' Deterministic per-replicate seed
#'
#' The seed depends only on (master_seed, scenario_id, rep_id) —
#' never on worker number, grid row or execution order.
#'
#' @param master_seed Master seed.
#' @param scenario_id Canonical scenario id.
#' @param rep_id Replicate number.
#' @return Integer seed between 1 and 2^31 - 1.
#' @keywords internal
.seed_for_rep <- function(master_seed, scenario_id, rep_id) {
  h <- .fnv1a_int(paste(master_seed, scenario_id, rep_id, sep = "\r"))
  as.integer(h %% 2147483647) + 1L
}

#' Short filesystem-safe scenario key
#' @keywords internal
.scenario_dir_key <- function(scenario_id) {
  sprintf("%010.0f", .fnv1a_int(scenario_id))
}

#' The formal pilot grid (15 scenarios, fixed mode/pipeline mapping)
#'
#' @return A data.frame with one row per formal scenario.
#' @keywords internal
.pilot_grid <- function() {
  rows <- list(
    # Simulation I: signal_resolution, full mode, condped_full
    c("signal_resolution", "null", "full", "condped_full"),
    c("signal_resolution", "single_multi_trait", "full", "condped_full"),
    c("signal_resolution", "two_linked_trait_specific", "full",
      "condped_full"),
    c("signal_resolution", "two_heterogeneous", "full", "condped_full"),
    # Simulation II Part A: signal_oracle
    c("trait_representation", "trait_specific", "signal_oracle",
      "condped_full"),
    c("trait_representation", "two_trait_concordant", "signal_oracle",
      "condped_full"),
    c("trait_representation", "two_trait_antagonistic", "signal_oracle",
      "condped_full"),
    c("trait_representation", "broad_concordant", "signal_oracle",
      "condped_full"),
    # Simulation II Part B: signal_trait_oracle
    c("trait_representation", "highly_representable",
      "signal_trait_oracle", "condped_full"),
    c("trait_representation", "partially_representable",
      "signal_trait_oracle", "condped_full"),
    c("trait_representation", "strongly_nonredundant",
      "signal_trait_oracle", "condped_full"),
    # Simulation III: end_to_end, full mode, paired three pipelines
    c("end_to_end", "single_highly_representable", "full", "paired"),
    c("end_to_end", "single_nonredundant", "full", "paired"),
    c("end_to_end", "linked_pseudo_multitrait", "full", "paired"),
    c("end_to_end", "mixed_multisignal", "full", "paired")
  )
  grid <- do.call(rbind, lapply(rows, function(r) {
    data.frame(
      experiment = r[1], scenario = r[2], analysis_mode = r[3],
      comparison_pipeline = r[4],
      n = 1000L, m = 4L, p = 1000L,
      locus_pve = 0.01, secondary_signal_pve = 0.005,
      local_ld = NA_character_, target_r2 = NA_real_,
      correlation = "block", target_loss = NA_real_,
      tolerance = 0.10, alpha_omnibus = 0.05,
      stringsAsFactors = FALSE
    )
  }))
  rownames(grid) <- NULL
  grid
}

#' Run CondPED simulation scenarios
#'
#' Executes a grid of simulation scenarios with reproducible random
#' seeds. Each replicate is uniquely determined by `(master_seed,
#' scenario_id, rep_id)`; parallel scheduling never changes the random
#' stream. Failed replicates are saved with their status instead of
#' being silently dropped.
#'
#' @param grid Data.frame describing the scenario grid (one row per
#'   scenario); see [`.pilot_grid()`] for the formal columns.
#' @param reps Integer; number of replicates per scenario.
#' @param out_dir Output directory for per-replicate RDS files.
#' @param master_seed Master random seed.
#' @param workers Integer; number of parallel workers.
#' @param backend Execution backend: `"sequential"` or `"parallel"`
#'   (fork-based; on Windows it falls back to sequential with a
#'   warning).
#' @param resume Logical; skip replicates whose output files already
#'   exist with `status = "ok"`.
#' @param overwrite Logical; rerun everything, including failed
#'   replicates.
#' @param save_data Logical; also save the simulated `G`/`Y`.
#' @param save_fit Logical; save the fitted null object.
#' @param fail_policy `"save_unstable"` saves failed replicates with
#'   status information; `"stop"` aborts on the first failure.
#' @param progress Logical; show a progress indicator.
#'
#' @return A list with components `manifest` (data.frame), `completed`,
#'   `failed`, `skipped`, `elapsed` and `status`.
#' @export
run_condped_simulation <- function(
  grid,
  reps,
  out_dir,
  master_seed = 20260804L,
  workers = 1L,
  backend = c("sequential", "parallel"),
  resume = TRUE,
  overwrite = FALSE,
  save_data = FALSE,
  save_fit = FALSE,
  fail_policy = c("save_unstable", "stop"),
  progress = interactive()
) {
  backend <- match.arg(backend)
  fail_policy <- match.arg(fail_policy)
  if (!is.data.frame(grid) || nrow(grid) == 0L) {
    .stop_invalid_input("grid must be a non-empty data.frame.")
  }
  need_cols <- c("experiment", "scenario", "analysis_mode",
                 "comparison_pipeline", "n", "m", "p", "locus_pve",
                 "secondary_signal_pve", "correlation", "tolerance")
  missing_cols <- setdiff(need_cols, names(grid))
  if (length(missing_cols) > 0L) {
    .stop_invalid_input("grid is missing column(s): %s.",
                        paste(missing_cols, collapse = ", "))
  }
  .check_count(reps, "reps", min = 1L)
  .check_count(workers, "workers", min = 1L)
  .check_count(master_seed, "master_seed", min = 1L)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # ---- task list ------------------------------------------------------------
  tasks <- list()
  for (gi in seq_len(nrow(grid))) {
    row <- grid[gi, , drop = FALSE]
    settings <- as.list(row)
    settings$master_seed <- master_seed
    scenario_id <- .canonical_scenario_id(settings)
    for (rep_id in seq_len(reps)) {
      tasks[[length(tasks) + 1L]] <- list(
        row = row, scenario_id = scenario_id, rep_id = rep_id
      )
    }
  }

  run_task <- function(task) {
    file <- file.path(out_dir, .scenario_dir_key(task$scenario_id),
                      sprintf("rep_%04d.rds", task$rep_id))
    if (!overwrite && file.exists(file)) {
      prev <- tryCatch(readRDS(file), error = function(e) NULL)
      if (!is.null(prev)) {
        # existing ok replicates are skipped on resume; failed ones are
        # kept and reported (rerun only when overwrite = TRUE)
        kept <- if (isTRUE(prev$status$ok)) {
          "skipped_ok"
        } else {
          "skipped_failed"
        }
        return(list(scenario_id = task$scenario_id, rep_id = task$rep_id,
                    status = kept, runtime = 0, file = file))
      }
    }
    .run_one_replicate(task$row, task$scenario_id, task$rep_id,
                       master_seed, file, save_data, save_fit)
  }

  t0 <- proc.time()[["elapsed"]]
  results <- if (workers > 1L && backend == "parallel" &&
                 .Platform$OS.type != "windows") {
    parallel::mclapply(tasks, run_task, mc.cores = workers)
  } else {
    if (workers > 1L && backend == "parallel") {
      warning("parallel backend is fork-based; falling back to sequential on this platform.")
    }
    lapply(tasks, run_task)
  }

  manifest <- do.call(rbind, lapply(results, function(r) {
    data.frame(scenario_id = r$scenario_id, rep_id = r$rep_id,
               status = r$status, runtime = r$runtime, file = r$file,
               stringsAsFactors = FALSE)
  }))
  failed <- manifest[!manifest$status %in%
                       c("ok", "skipped_ok", "skipped_failed"),
                     , drop = FALSE]
  if (fail_policy == "stop" && nrow(failed) > 0L) {
    stop("run_condped_simulation stopped on failed replicate: ",
         failed$scenario_id[1L], " rep ", failed$rep_id[1L],
         call. = FALSE)
  }
  list(
    manifest = manifest,
    completed = sum(manifest$status == "ok"),
    failed = failed,
    skipped = sum(manifest$status %in% c("skipped_ok", "skipped_failed")),
    elapsed = proc.time()[["elapsed"]] - t0,
    status = .new_status(
      ok = nrow(failed) == 0L,
      code = if (nrow(failed) == 0L) "ok" else "unstable",
      message = if (nrow(failed) == 0L) "" else
        sprintf("%d replicate(s) failed or were unstable.", nrow(failed))
    )
  )
}

#' Execute one replicate and save it atomically
#'
#' @param row One grid row (data.frame).
#' @param scenario_id Canonical id.
#' @param rep_id Replicate number.
#' @param master_seed Master seed.
#' @param file Target RDS path.
#' @param save_data,save_fit Saving switches.
#' @return A compact per-task list (status/runtime/file).
#' @keywords internal
.run_one_replicate <- function(row, scenario_id, rep_id, master_seed,
                               file, save_data, save_fit) {
  t0 <- proc.time()[["elapsed"]]
  seed <- .seed_for_rep(master_seed, scenario_id, rep_id)
  set.seed(seed)
  settings <- as.list(row)
  settings$master_seed <- master_seed

  result <- tryCatch(
    .execute_replicate(row, settings, seed, rep_id, save_data, save_fit),
    error = function(e) e
  )
  runtime <- proc.time()[["elapsed"]] - t0

  if (inherits(result, "error")) {
    obj <- list(
      scenario_id = scenario_id, rep_id = rep_id, seed = seed,
      settings = settings,
      status = .new_status(ok = FALSE, code = "unstable",
                           message = conditionMessage(result)),
      run_status = "failed",
      failure_stage = "runtime",
      error_class = class(result)[1L],
      error_message = conditionMessage(result),
      runtime = runtime
    )
  } else {
    obj <- result
    obj$runtime <- runtime
  }
  .save_replicate_atomic(obj, file)
  run_status <- if (!is.null(obj$run_status)) {
    obj$run_status
  } else if (isTRUE(obj$status$ok)) {
    "ok"
  } else {
    obj$status$code
  }
  list(scenario_id = scenario_id, rep_id = rep_id,
       status = run_status, runtime = runtime, file = file)
}

#' The frozen replicate pipeline (calls only frozen functions)
#' @keywords internal
.execute_replicate <- function(row, settings, seed, rep_id,
                               save_data, save_fit) {
  mode <- row$analysis_mode
  pipeline <- row$comparison_pipeline
  scenario_id <- .canonical_scenario_id(settings)

  # ---- simulate ---------------------------------------------------------------
  sim_args <- list(
    n = row$n, m = row$m, p = row$p,
    experiment = row$experiment, scenario = row$scenario,
    locus_pve = row$locus_pve,
    secondary_signal_pve = row$secondary_signal_pve,
    correlation = row$correlation, tolerance = row$tolerance
  )
  if (!is.null(row$local_ld) && !is.na(row$local_ld)) {
    sim_args$local_ld <- row$local_ld
  }
  if (!is.null(row$target_r2) && !is.na(row$target_r2)) {
    sim_args$target_r2 <- row$target_r2
  }
  if (!is.null(row$target_loss) && !is.na(row$target_loss)) {
    sim_args$target_loss <- row$target_loss
  }
  sim <- do.call(simulate_condped_data, sim_args)
  if (!isTRUE(sim$status$ok)) {
    return(list(
      scenario_id = scenario_id,
      rep_id = rep_id, seed = seed, settings = settings,
      truth = sim$truth,
      status = .new_status(ok = FALSE, code = "unstable",
                           message = sim$status$message),
      run_status = "unstable",
      failure_stage = "simulation",
      error_class = sim$status$code,
      error_message = sim$status$message,
      runtime = NA_real_
    ))
  }
  G <- sim$G
  marker_ids <- colnames(G)
  position <- stats::setNames(seq_along(marker_ids) * 1000, marker_ids)

  # ---- one null fit --------------------------------------------------------------
  fit <- fit_mt_null(sim$Y, K = sim$K_bg,
                     control = list(maxit = 300L))
  if (!isTRUE(fit$status$ok)) {
    return(list(
      scenario_id = scenario_id,
      rep_id = rep_id, seed = seed, settings = settings,
      truth = sim$truth,
      status = .new_status(ok = FALSE, code = "non_convergence",
                           message = "fit_mt_null() did not converge."),
      run_status = "failed",
      failure_stage = "null_fit",
      error_class = "non_convergence",
      error_message = "fit_mt_null() did not converge.",
      runtime = NA_real_
    ))
  }
  basis <- derive_conditional_contrasts(fit$Sigma_P_ref)

  empty_long <- data.frame(
    locus_id = character(), signal_id = character(),
    representative_snp = character(), trait = character(),
    beta = numeric(), se = numeric(), p_value = numeric(),
    stringsAsFactors = FALSE
  )

  # ---- analysis ----------------------------------------------------------------
  scan <- NULL
  resolution <- NULL
  loci_obj <- NULL
  asset_results <- NULL
  if (mode == "full") {
    scan <- scan_mt_omnibus(fit, G)
    om <- scan$omnibus
    valid <- !is.na(om$p_value)
    p_adj <- rep(NA_real_, nrow(om))
    p_adj[valid] <- stats::p.adjust(om$p_value[valid], method = "BH")
    selected_markers <- om$marker_id[valid][
      p_adj[valid] <= (row$alpha_omnibus %||% 0.05)]
    if (length(selected_markers) > 0L) {
      loci_obj <- define_associated_loci(
        scan, G, chromosome = rep("chr1", ncol(G)),
        position = unname(position),
        selected_markers = selected_markers,
        marker_ids = marker_ids,
        method = "physical", window_bp = 5000
      )
      resolution <- resolve_locus_signals(fit, G, loci_obj,
                                          marker_ids = marker_ids)
      if (pipeline == "paired") {
        pipes <- .run_comparison_pipelines(
          fit, G, loci_obj, scan,
          trait_names = fit$trait_names,
          tolerance = row$tolerance
        )
        asset_results <- list(
          lead_asset = pipes$lead_asset,
          resolved_asset = pipes$resolved_asset
        )
        condped_out <- pipes$condped_full
      } else {
        condped_out <- .pipeline_condped_full(
          fit, resolution, alpha_trait = 0.05,
          tolerance = row$tolerance
        )
      }
      attr <- list(
        candidate_sets = condped_out$candidate_sets,
        signal_table = condped_out$signal_table,
        trait_table = condped_out$trait_table
      )
      dec <- condped_out$subset_analysis
      effects_obj <- resolution
    } else {
      attr <- attribute_traits(empty_long, candidate_mode = "all_traits")
      dec <- decompose_conditional_effects(
        empty_long, attr, basis, tolerance = row$tolerance
      )
      effects_obj <- NULL
      condped_out <- NULL
    }
  } else {
    # oracle modes: truth-defined signals, no selection
    ps <- .oracle_signals_from_truth(sim$truth)
    effects_obj <- .effects_from_predefined(
      fit, G, ps, marker_ids, sqrt(.Machine$double.eps)
    )
    attr <- if (mode == "signal_trait_oracle") {
      attribute_traits(effects_obj, candidate_mode = "predefined",
                       predefined_sets = sim$truth$candidate_traits)
    } else {
      attribute_traits(effects_obj, candidate_mode = "holm_fwer")
    }
    dec <- decompose_conditional_effects(
      effects_obj, attr, basis, tolerance = row$tolerance
    )
    condped_out <- NULL
  }

  # ---- evaluator-facing estimates --------------------------------------------------
  signals_df <- if (!is.null(effects_obj) &&
                    !is.null(effects_obj$signals)) {
    s <- effects_obj$signals
    if (nrow(s) > 0L) {
      s$position <- unname(position[s$representative_snp])
    }
    s
  } else {
    data.frame()
  }
  estimates <- list(
    signals = signals_df,
    beta = if (!is.null(effects_obj)) effects_obj$beta else NULL,
    candidate_sets = attr$candidate_sets,
    minimum_representative_sets = dec$minimum_representative_sets,
    irreducible_modules = dec$irreducible_modules,
    subset_table = dec$subset_table,
    omnibus_summary = if (!is.null(scan)) scan$omnibus else NULL,
    positions = position,
    asset = .asset_estimates_for_eval(asset_results, condped_out)
  )
  eval_settings <- list(
    experiment = row$experiment,
    architecture = row$scenario,
    n = row$n,
    locus_pve = row$locus_pve,
    correlation = row$correlation,
    target_loss = if (!is.null(row$target_loss)) row$target_loss else NA,
    tolerance = row$tolerance,
    alpha_omnibus = row$alpha_omnibus %||% 0.05
  )
  settings$architecture <- row$scenario   # evaluator grouping key
  evaluation <- tryCatch(
    .evaluate_one_replicate(sim$truth, estimates, eval_settings,
                            runtime = NA_real_, G = G),
    error = function(e) e
  )
  if (inherits(evaluation, "error")) {
    evaluation <- list(error = conditionMessage(evaluation))
  }

  out <- list(
    scenario_id = scenario_id,
    rep_id = rep_id,
    seed = seed,
    settings = settings,
    truth = sim$truth,
    omnibus_summary = if (!is.null(scan)) scan$omnibus else NULL,
    loci = if (!is.null(loci_obj)) loci_obj$loci else NULL,
    signals = signals_df,
    effect_estimates = if (!is.null(effects_obj)) {
      effects_obj$beta
    } else {
      NULL
    },
    candidate_trait_estimates = attr,
    subset_table = dec$subset_table,
    minimum_representative_sets = dec$minimum_representative_sets,
    irreducible_modules = dec$irreducible_modules,
    tolerance_path = dec$tolerance_path,
    asset_results = asset_results,
    comparison_pipeline = pipeline,
    Sigma_P_ref = fit$Sigma_P_ref,
    evaluation = evaluation,
    estimates = estimates,
    status = .new_status(ok = TRUE, code = "ok"),
    run_status = "ok",
    diagnostics = list(
      n_signals = nrow(signals_df),
      mode = mode,
      pipeline = pipeline
    ),
    runtime = NA_real_
  )
  if (!save_data) {
    out$G <- NULL
  } else {
    out$G <- G
  }
  if (save_fit) out$fit <- fit
  out
}

#' Predefined-signals table from simulation truth (oracle modes)
#' @keywords internal
.oracle_signals_from_truth <- function(truth) {
  sig <- truth$signals
  if (nrow(sig) == 0L) {
    .stop_invalid_input("truth contains no signals for the oracle mode.")
  }
  conditioning <- vector("list", nrow(sig))
  for (i in seq_len(nrow(sig))) {
    others <- sig$representative_snp[
      sig$locus_id == sig$locus_id[i] &
        sig$representative_snp != sig$representative_snp[i]]
    conditioning[[i]] <- others
  }
  data.frame(
    locus_id = sig$locus_id,
    signal_id = sig$signal_id,
    representative_snp = sig$representative_snp,
    conditioning_snps = I(conditioning),
    stringsAsFactors = FALSE
  )
}

#' Assemble the evaluator-facing asset structure
#' @keywords internal
.asset_estimates_for_eval <- function(asset_results, condped_out) {
  if (is.null(asset_results) && is.null(condped_out)) return(NULL)
  pack <- function(entries, id_field) {
    if (is.null(entries)) return(NULL)
    items <- if (!is.null(entries$loci)) {
      entries$loci
    } else {
      entries$signals
    }
    if (is.null(items) || length(items) == 0L) {
      return(list(sets = stats::setNames(list(), character()),
                  p = stats::setNames(numeric(), character()),
                  direction = stats::setNames(character(), character())))
    }
    sets <- lapply(items, function(x) x$asset_best_subset)
    p <- vapply(items, function(x) {
      if (is.null(x$asset_p)) NA_real_ else x$asset_p
    }, numeric(1))
    # direction representation derived only from the real ASSET
    # positive/negative subsets (never recomputing the ASSET statistic)
    dir <- vapply(items, function(x) {
      .direction_from_subsets(x$asset_positive_subset,
                              x$asset_negative_subset)
    }, character(1))
    names(sets) <- names(p) <- names(dir) <- names(items)
    list(sets = sets, p = p, direction = dir)
  }
  out <- list()
  if (!is.null(asset_results)) {
    out$lead <- pack(asset_results$lead_asset, "loci")
    out$resolved <- pack(asset_results$resolved_asset, "signals")
  }
  if (!is.null(condped_out)) {
    cs <- condped_out$candidate_sets
    dir <- NULL
    st <- condped_out$signal_table
    if (!is.null(st) &&
        all(c("signal_id", "direction_pattern") %in% names(st))) {
      dir <- stats::setNames(as.character(st$direction_pattern),
                             as.character(st$signal_id))
    }
    out$condped <- list(
      sets = cs,
      p = stats::setNames(rep(NA_real_, length(cs)), names(cs)),
      direction = dir
    )
  }
  if (length(out) == 0L) NULL else out
}

#' Direction label from ASSET positive/negative trait subsets
#'
#' Maps the real ASSET subset output onto the frozen direction
#' vocabulary used by signal summaries and simulation truth.
#' @keywords internal
.direction_from_subsets <- function(positive, negative) {
  pos <- if (is.null(positive)) character() else as.character(positive)
  neg <- if (is.null(negative)) character() else as.character(negative)
  u <- union(pos, neg)
  if (length(u) == 0L) return("not_applicable")
  if (length(u) == 1L) return("single_trait")
  if (length(pos) > 0L && length(neg) > 0L) {
    if (length(u) == 2L) "antagonistic" else "mixed"
  } else {
    "concordant"
  }
}

#' Atomic save: write to a temporary file, then rename
#' @keywords internal
.save_replicate_atomic <- function(obj, file) {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  tmp <- tempfile(pattern = ".rep_", tmpdir = dirname(file),
                  fileext = ".tmp")
  saveRDS(obj, tmp)
  if (!file.rename(tmp, file)) {
    unlink(tmp)
    .stop_invalid_input("failed to rename temporary replicate file.")
  }
  invisible(file)
}
