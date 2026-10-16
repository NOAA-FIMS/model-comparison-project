#' Summarize optimization and uncertainty-reporting runtimes
#'
#' Reads timing tables (component, time, unit) or named numeric vectors in
#' seconds. ADMB uncertainty reporting uses fit_sdreport: Hessian evaluation
#' and inversion, delta-method gradient sweeps, and .std/.cor writing, but not
#' .rep reporting. TMB uses tmb_sdreport: Hessian calculation without the
#' delta method. These stages therefore do not represent equivalent workloads.
#' Missing files, components, and NA timings remain NA in the returned details
#' and are excluded from summary statistics. Entirely missing groups are omitted.
#' Unreadable files and negative or infinite timings warn and are treated as NA.
#'
#' @param em_names Character vector drawn from ASAP, BAM, SS, WHAM, FIMS,
#'   and Rceattle. SS is labeled SS3 in the output.
#' @param n_sim Positive integer number of operating-model simulations.
#' @param case_dir Character path to the case directory.
#' @return A long-format tibble with sim, model, stage, component, and
#'   runtime_sec. Also writes figure/runtime_optimization_summary.csv and
#'   figure/runtime_uncertainty_reporting_summary.csv with counts, mean, median,
#'   minimum, and maximum seconds for each model/scenario. The CSV files omit
#'   stage and component columns; these remain in the returned details.
#' @examples
#' source("scripts/check_runtime.R")
#' runtimes <- check_runtime(c("ASAP", "BAM", "SS", "FIMS", "Rceattle", "WHAM"),
#'                          n_sim = 10, case_dir = "C0")
check_runtime <- function(em_names, n_sim, case_dir) {
  supported_models <- c("ASAP", "BAM", "SS", "FIMS", "Rceattle", "WHAM")
  if (!is.character(em_names) || length(em_names) == 0L ||
      anyNA(em_names) || any(!em_names %in% supported_models) ||
      anyDuplicated(em_names)) {
    stop("em_names must contain unique names from: ",
         paste(supported_models, collapse = ", "), ".")
  }
  if (!is.numeric(n_sim) || length(n_sim) != 1L || is.na(n_sim) ||
      !is.finite(n_sim) || n_sim < 1 || n_sim != floor(n_sim)) {
    stop("n_sim must be a positive integer.")
  }

  scenarios <- c("fixed_effects", "random_effects", "random_effects_sigmaR_constant")
  labels <- c("penalized_MLE", "state_space_unknown_sigmaR", "state_space_known_sigmaR")
  model_order <- c(
    "ASAP", "BAM", "SS3",
    paste0("FIMS_", labels),
    paste0("Rceattle_", labels),
    paste0("WHAM_", labels)
  )
  stages <- c("optimization", "uncertainty_reporting")

  read_timing <- function(path, components) {
    missing <- rep(NA_real_, length(components))
    if (!file.exists(path)) return(missing)
    timing <- tryCatch(readRDS(path), error = function(e) {
      warning("check_runtime(): cannot read ", path, ": ",
              conditionMessage(e), call. = FALSE)
      NULL
    })
    if (is.null(timing)) return(missing)

    if (is.data.frame(timing) &&
        all(c("component", "time", "unit") %in% names(timing)) &&
        is.character(timing[["component"]]) && is.numeric(timing[["time"]])) {
      selected <- timing[["component"]] %in% components
      if (anyDuplicated(timing[["component"]][selected]) ||
          anyNA(timing[["unit"]][selected]) ||
          any(!timing[["unit"]][selected] %in% c("seconds", "secs", "s"))) {
        stop("check_runtime(): duplicate components or invalid time units in ", path)
      }
      values <- timing[["time"]][match(components, timing[["component"]])]
    } else if (is.numeric(timing) && is.null(dim(timing)) &&
               !is.null(names(timing))) {
      if (anyDuplicated(names(timing)[names(timing) %in% components])) {
        stop("check_runtime(): duplicate timing components in ", path)
      }
      values <- unname(timing[components])
    } else {
      stop("check_runtime(): expected a component/time/unit table or a named ",
           "numeric vector in ", path)
    }
    invalid <- !is.na(values) & (!is.finite(values) | values < 0)
    if (any(invalid)) {
      warning("check_runtime(): excluding negative or infinite timing(s) for ",
              paste(components[invalid], collapse = ", "), " in ", path,
              call. = FALSE)
      values[invalid] <- NA_real_
    }
    values
  }

  results <- list()
  for (model in supported_models[supported_models %in% em_names]) {
    is_admb <- model %in% c("ASAP", "BAM", "SS")
    components <- c("fit_optimization", if (is_admb) "fit_sdreport" else "tmb_sdreport")
    model_scenarios <- if (is_admb) NA_character_ else scenarios
    for (i in seq_along(model_scenarios)) {
      scenario <- model_scenarios[[i]]
      if (is_admb) {
        filename <- switch(model,
          ASAP = "runtime_asap.RDS",
          BAM = "runtime_bam.RDS",
          SS = "runtime_ss3.RDS"
        )
        output_model <- if (model == "SS") "SS3" else model
      } else {
        filename <- if (model == "WHAM") {
          wham_scenario <- if (scenario == "fixed_effects") "fixed_effects_logNAA" else scenario
          paste0("runtime_", wham_scenario, ".RDS")
        } else {
          paste0("run_time_", tolower(model), "_", scenario, ".RDS")
        }
        output_model <- paste0(model, "_", labels[[i]])
      }
      for (sim in seq_len(n_sim)) {
        path <- file.path(case_dir, "output", model, paste0("s", sim), filename)
        results[[length(results) + 1L]] <- tibble::tibble(
          sim = sim,
          model = output_model,
          stage = stages,
          component = components,
          runtime_sec = read_timing(path, components)
        )
      }
    }
  }
  results <- dplyr::bind_rows(results) |>
    dplyr::arrange(match(model, model_order), sim, match(stage, stages))

  runtime_summary <- results |>
    dplyr::filter(!is.na(runtime_sec)) |>
    dplyr::mutate(
      model = factor(model, levels = model_order),
      stage = factor(stage, levels = stages)
    ) |>
    dplyr::group_by(model, stage, component) |>
    dplyr::summarize(
      total_simulations = n_sim,
      recorded_runtimes = dplyr::n(),
      mean_runtime_sec = mean(runtime_sec, na.rm = TRUE),
      median_runtime_sec = median(runtime_sec, na.rm = TRUE),
      min_runtime_sec = if (dplyr::n() > 0L) min(runtime_sec, na.rm = TRUE) else NA_real_,
      max_runtime_sec = if (dplyr::n() > 0L) max(runtime_sec, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) |>
    dplyr::mutate(dplyr::across(dplyr::where(is.double), ~ round(.x, 2))) |>
    dplyr::arrange(model, stage)

  dir.create(file.path(case_dir, "figure"), showWarnings = FALSE, recursive = TRUE)
  for (runtime_stage in stages) {
    stage_summary <- runtime_summary |>
      dplyr::filter(stage == runtime_stage) |>
      dplyr::select(-stage, -component)
    write.csv(
      stage_summary,
      file.path(case_dir, "figure", paste0("runtime_", runtime_stage, "_summary.csv")),
      row.names = FALSE
    )
  }
  results
}
