#' Check execution run time for ASAP, BAM, SS3, WHAM, and FIMS estimation models
#' @param em_names Character vector of estimation model names (e.g., c("ASAP", "BAM", "SS", "WHAM", "FIMS"))
#' @param n_sim Integer, number of operating model simulations
#' @param case_dir Character, path to case working directory
#' @return A long-format tibble containing run time details for all simulations and models
check_runtime <- function(em_names, n_sim, case_dir) {
  
  # Separate models by type for different file path approaches
  admb_models <- em_names[em_names %in% c("ASAP", "BAM", "SS")]
  wham_models <- em_names[em_names == "WHAM"]
  fims_models <- em_names[em_names == "FIMS"]
  
  # Initialize empty results tibble to store run time info for all models
  results <- tibble::tibble()
  
  # Process ADMB models (ASAP, BAM, SS3)
  if (length(admb_models) > 0) {
    admb_results <- expand.grid(
      sim = seq_len(n_sim),
      model = admb_models,
      stringsAsFactors = FALSE
    ) |> 
      tibble::as_tibble() |> 
      dplyr::mutate(
        sim_dir = file.path(case_dir, "output", model, paste0("s", sim)),
        # Map to specific run time RDS files
        runtime_file = dplyr::case_when(
          model == "ASAP" ~ file.path(sim_dir, "runtime_asap.RDS"),
          model == "BAM"  ~ file.path(sim_dir, "runtime_bam.RDS"),
          model == "SS"   ~ file.path(sim_dir, "runtime_ss3.RDS"),
          TRUE ~ NA_character_
        ),
        # Safely read RDS run time value
        runtime = purrr::map_dbl(runtime_file, ~ {
          if (!is.na(.x) && file.exists(.x)) {
            tryCatch({
              readRDS(.x)
            }, error = function(e) NA_real_)
          } else {
            NA_real_
          }
        }),
        # Convert SS to SS3 for unified reporting naming convention
        model = dplyr::if_else(model == "SS", "SS3", model)
      ) |>
      dplyr::select(sim, model, runtime)
    
    results <- dplyr::bind_rows(results, admb_results)
  }
  
  # Process WHAM models
  if (length(wham_models) > 0) {
    wham_results <- expand.grid(
      sim = seq_len(n_sim),
      model = wham_models,
      stringsAsFactors = FALSE
    ) |> 
      tibble::as_tibble() |> 
      dplyr::mutate(
        sim_dir = file.path(case_dir, "output", model, paste0("s", sim)),
        fit_file_fixed = file.path(sim_dir, "runtime_fixed_effects_logNAA.RDS"),
        fit_file_random = file.path(sim_dir, "runtime_random_effects.RDS"),
        fit_file_random_sigmaR_constant = file.path(sim_dir, "runtime_random_effects_sigmaR_constant.RDS"),
        
        # Read run times for all three variations safely
        runtime_fixed = purrr::map_dbl(fit_file_fixed, ~ {
          if (file.exists(.x)) tryCatch(readRDS(.x), error = function(e) NA_real_) else NA_real_
        }),
        runtime_random = purrr::map_dbl(fit_file_random, ~ {
          if (file.exists(.x)) tryCatch(readRDS(.x), error = function(e) NA_real_) else NA_real_
        }),
        runtime_random_sigmaR_constant = purrr::map_dbl(fit_file_random_sigmaR_constant, ~ {
          if (file.exists(.x)) tryCatch(readRDS(.x), error = function(e) NA_real_) else NA_real_
        })
      ) |>
      dplyr::select(sim, runtime_fixed, runtime_random, runtime_random_sigmaR_constant) |>
      # Pivot to long format using the tight anchored regex fix
      tidyr::pivot_longer(
        cols = c(runtime_fixed, runtime_random, runtime_random_sigmaR_constant),
        names_to = c(".value", "effect_type"),
        names_pattern = "^(.+)_(random_sigmaR_constant|random|fixed)$"
      ) |>
      dplyr::mutate(
        model = dplyr::case_when(
          effect_type == "fixed" ~ "WHAM_fixed_effects",
          effect_type == "random" ~ "WHAM_random_effects",
          effect_type == "random_sigmaR_constant" ~ "WHAM_random_effects_sigmaR_constant",
          TRUE ~ "WHAM"
        )
      ) |>
      dplyr::select(sim, model, runtime)
    
    results <- dplyr::bind_rows(results, wham_results)
  }
  
  # Process FIMS models
  if (length(fims_models) > 0) {
    fims_results <- expand.grid(
      sim = seq_len(n_sim),
      model = fims_models,
      stringsAsFactors = FALSE
    ) |> 
      tibble::as_tibble() |> 
      dplyr::mutate(
        sim_dir = file.path(case_dir, "output", model, paste0("s", sim)),
        fit_file_fixed = file.path(sim_dir, "run_time_fims_fixed_effects.RDS"),
        fit_file_random = file.path(sim_dir, "run_time_fims_random_effects.RDS"),
        fit_file_random_sigmaR_constant = file.path(sim_dir, "run_time_fims_random_effects_sigmaR_constant.RDS"),
        
        # Read run times safely
        runtime_fixed = purrr::map_dbl(fit_file_fixed, ~ {
          if (file.exists(.x)) tryCatch(readRDS(.x)[["fit_total"]], error = function(e) NA_real_) else NA_real_
        }),
        runtime_random = purrr::map_dbl(fit_file_random, ~ {
          if (file.exists(.x)) tryCatch(readRDS(.x)[["fit_total"]], error = function(e) NA_real_) else NA_real_
        }),
        runtime_random_sigmaR_constant = purrr::map_dbl(fit_file_random_sigmaR_constant, ~ {
          if (file.exists(.x)) tryCatch(readRDS(.x)[["fit_total"]], error = function(e) NA_real_) else NA_real_
        })
      ) |>
      dplyr::select(sim, runtime_fixed, runtime_random, runtime_random_sigmaR_constant) |>
      tidyr::pivot_longer(
        cols = c(runtime_fixed, runtime_random, runtime_random_sigmaR_constant),
        names_to = c(".value", "effect_type"),
        names_pattern = "^(.+)_(random_sigmaR_constant|random|fixed)$"
      ) |>
      dplyr::mutate(
        model = dplyr::case_when(
          effect_type == "fixed" ~ "FIMS_fixed_effects",
          effect_type == "random" ~ "FIMS_random_effects",
          effect_type == "random_sigmaR_constant" ~ "FIMS_random_effects_sigmaR_constant",
          TRUE ~ "FIMS"
        )
      ) |>
      dplyr::select(sim, model, runtime)
    
    results <- dplyr::bind_rows(results, fims_results)
  }
  
  # Calculate summary metrics for model runtimes
  runtime_summary <- results |>
    dplyr::mutate(
      model = factor(
        model, 
        levels = c(
          "ASAP", "BAM", "SS3", 
          "WHAM_fixed_effects", "WHAM_random_effects", "WHAM_random_effects_sigmaR_constant",
          "FIMS_fixed_effects", "FIMS_random_effects", "FIMS_random_effects_sigmaR_constant"
        )
      )
    ) |>
    dplyr::group_by(model) |>
    dplyr::summarize(
      total_simulations  = dplyr::n(),
      recorded_runtimes  = sum(!is.na(runtime)),
      mean_runtime_sec   = mean(runtime, na.rm = TRUE),
      median_runtime_sec = median(runtime, na.rm = TRUE),
      min_runtime_sec    = min(runtime, na.rm = TRUE),
      max_runtime_sec    = max(runtime, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(dplyr::across(dplyr::where(is.double), ~ round(., 2))) |>
    dplyr::arrange(model)
  
  # Save the run time summaries to a CSV file matching your pipeline pattern
  write.csv(runtime_summary, file.path(case_dir, "figure", "runtime_summary.csv"), row.names = FALSE)
  
  # Return the detailed raw results dataset for plotting functions (like boxplots/violins)
  return(results)
}