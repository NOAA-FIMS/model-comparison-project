#' Check convergence for ASAP, BAM, SS, WHAM, and FIMS estimation models
#' @param em_names Character vector of estimation model names (e.g., c("ASAP", "BAM", "SS", "WHAM", "FIMS"))
#' @param n_sim Integer, number of operating model simulations
#' @param case_dir Character, path to case working directory
#' @param gradient_threshold Numeric, threshold for gradient convergence (default = 0.01)
#' @return A vector of simulation IDs that converged for ALL models
check_convergence <- function(em_names, n_sim, case_dir, gradient_threshold = 0.01) {
  
  # Separate models by type for different convergence checking approaches
  # ADMB models (ASAP, BAM, SS) use .par files and admodel.cov
  # WHAM uses convergence status from fit
  # FIMS uses fit objects with FIMS::get_max_gradient() function
  admb_models <- em_names[em_names %in% c("ASAP", "BAM", "SS")]
  wham_models <- em_names[em_names == "WHAM"]
  fims_models <- em_names[em_names == "FIMS"]
  
  # Initialize empty results tibble to store convergence info for all models
  results <- tibble::tibble()
  
  # Process ADMB models (ASAP, BAM, SS3)
  if (length(admb_models) > 0) {
    # Helper function to find the correct .par file for each model type
    get_par_file <- function(sim_dir, model) {
      par_files <- c("asap3.par", "BAM-Sim.par", "ss.par")
      found <- par_files[file.exists(file.path(sim_dir, par_files))]
      if (length(found) > 0) file.path(sim_dir, found[1]) else NA_character_
    }
    
    # Create a grid of all sim-model combinations
    admb_results <- expand.grid(
      sim = seq_len(n_sim),
      model = admb_models,
      stringsAsFactors = FALSE
    ) |> 
      tibble::as_tibble() |> 
      dplyr::mutate(
        # Build paths to output directories for each simulation
        sim_dir = file.path(case_dir, "output", model, paste0("s", sim)),
        cov_path = file.path(sim_dir, "admodel.cov"),
        par_path = purrr::map2_chr(sim_dir, model, get_par_file),
        # positive_hessian = TRUE if admodel.cov exists
        positive_hessian = file.exists(cov_path),
        # Extract gradient from .par file (16th value after extracting positions 6, 11, 16)
        # This is the maximum gradient at convergence
        gradient = purrr::map_dbl(par_path, ~{
          if (file.exists(.x)) {
            as.numeric(scan(.x, what = '', n = 16, quiet = TRUE)[c(6, 11, 16)])[3]
          } else {
            NA_real_
          }
        })
      ) |>
      dplyr::select(sim, model, positive_hessian, gradient) |>
      dplyr::mutate(
        model = dplyr::if_else(model == "SS", "SS3", model)
      )
    
    results <- dplyr::bind_rows(results, admb_results)
  }
  
  # Process WHAM models
  # positive_hessian: based on na_sdrep (FALSE if na_sdrep=TRUE, TRUE if na_sdrep=FALSE)
  # Create separate rows for fixed and random effects to track them independently
  if (length(wham_models) > 0) {
    wham_results <- expand.grid(
      sim = seq_len(n_sim),
      model = wham_models,
      stringsAsFactors = FALSE
    ) |> 
      tibble::as_tibble() |> 
      dplyr::mutate(
        sim_dir = file.path(case_dir, "output", model, paste0("s", sim)),
        fit_file_fixed = file.path(sim_dir, "fit_wham_fixed_effects.RDS"),
        fit_file_random = file.path(sim_dir, "fit_wham_random_effects.RDS"),
        # Load WHAM fixed effects gradient and positive_hessian
        convergence_fixed = purrr::map_dbl(fit_file_fixed, ~{
          if (file.exists(.x)) {
            tryCatch({
              fit <- readRDS(.x)
              max(abs(fit[["final_gradient"]]))
            }, error = function(e) NA_real_)
          } else {
            NA_real_
          }
        }),
        positive_hessian_fixed = purrr::map_lgl(fit_file_fixed, ~{
          if (file.exists(.x)) {
            tryCatch({
              fit <- readRDS(.x)
              # positive_hessian is TRUE if na_sdrep is FALSE
              !fit[["na_sdrep"]]
            }, error = function(e) NA)
          } else {
            NA
          }
        }),
        # Load WHAM random effects gradient and positive_hessian
        convergence_random = purrr::map_dbl(fit_file_random, ~{
          if (file.exists(.x)) {
            tryCatch({
              fit <- readRDS(.x)
              max(abs(fit[["final_gradient"]]))
            }, error = function(e) NA_real_)
          } else {
            NA_real_
          }
        }),
        positive_hessian_random = purrr::map_lgl(fit_file_random, ~{
          if (file.exists(.x)) {
            tryCatch({
              fit <- readRDS(.x)
              # positive_hessian is TRUE if na_sdrep is FALSE
              !fit[["na_sdrep"]]
            }, error = function(e) NA)
          } else {
            NA
          }
        })
      ) |>
      dplyr::select(sim, convergence_fixed, convergence_random, positive_hessian_fixed, positive_hessian_random) |>
      # Reshape to long format: one row per effect type per simulation
      # This allows tracking fixed and random effects convergence separately
      # Uses names_pattern to reshape both convergence and positive_hessian columns together
      tidyr::pivot_longer(
        cols = c(convergence_fixed, convergence_random, positive_hessian_fixed, positive_hessian_random),
        names_to = c(".value", "effect_type"),
        names_pattern = "(.+)_(fixed|random)"
      ) |>
      dplyr::mutate(
        # Create descriptive model names for each effect type
        model = dplyr::case_when(
          effect_type == "fixed" ~ "WHAM_fixed_effects",
          effect_type == "random" ~ "WHAM_random_effects",
          TRUE ~ "WHAM"
        ),
        # Rename convergence column to gradient for consistency with other models
        gradient = convergence
      ) |>
      dplyr::select(sim, model, positive_hessian, gradient)
    
    results <- dplyr::bind_rows(results, wham_results)
  }
  
  # Process FIMS models
  # FIMS models store fit objects in RDS files for both fixed and random effects
  # Use FIMS::get_max_gradient() to extract the maximum gradient from each fit
  # Like WHAM, we track fixed and random effects separately
  if (length(fims_models) > 0) {
    fims_results <- expand.grid(
      sim = seq_len(n_sim),
      model = fims_models,
      stringsAsFactors = FALSE
    ) |> 
      tibble::as_tibble() |> 
      dplyr::mutate(
        sim_dir = file.path(case_dir, "output", model, paste0("s", sim)),
        fit_file_fixed = file.path(sim_dir, "fit_fims_fixed_effects.RDS"),
        fit_file_random = file.path(sim_dir, "fit_fims_random_effects.RDS"),
        gradient_file_fixed = file.path(sim_dir, "max_gradient_fims_fixed_effects.RDS"),
        gradient_file_random = file.path(sim_dir, "max_gradient_fims_random_effects.RDS"),
        hessian_file_fixed = file.path(sim_dir, "hessian_fims_fixed_effects.RDS"),
        hessian_file_random = file.path(sim_dir, "hessian_fims_random_effects.RDS"),
        # Load FIMS fixed effects maximum gradient
        gradient_fixed = purrr::map_dbl(gradient_file_fixed, ~{
          if (file.exists(.x)) {
            tryCatch({
              readRDS(.x)
            }, error = function(e) NA_real_)
          } else {
            NA_real_
          }
        }),
        # Load Hession for fix effects models
        positive_hessian_fixed = purrr::map_lgl(hessian_file_fixed, ~{
          if (file.exists(.x)) {
            tryCatch({
              readRDS(.x)
            }, error = function(e) NA)
          } else {
            NA
          }
        }),
        # Load FIMS random effects maximum gradient
        gradient_random = purrr::map_dbl(gradient_file_random, ~{
          if (file.exists(.x)) {
            tryCatch({
              readRDS(.x)
            }, error = function(e) NA_real_)
          } else {
            NA_real_
          }
        }),
        # Load Hession for random effects models
        positive_hessian_random = purrr::map_lgl(hessian_file_random, ~{
          if (file.exists(.x)) {
            tryCatch({
              readRDS(.x)
            }, error = function(e) NA)
          } else {
            NA
          }
        })
      ) |>
      dplyr::select(
        sim, 
        gradient_fixed, gradient_random, 
        positive_hessian_fixed, positive_hessian_random
      ) |>
      # Reshape to long format: one row per effect type per simulation
      tidyr::pivot_longer(
        cols = c(
          gradient_fixed, gradient_random,
          positive_hessian_fixed, positive_hessian_random
        ),
        names_to = c(".value", "effect_type"),
        names_pattern = "(.+)_(fixed|random)"
      ) |>
      dplyr::mutate(
        # Create descriptive model names for each effect type
        model = dplyr::case_when(
          effect_type == "fixed" ~ "FIMS_fixed_effects",
          effect_type == "random" ~ "FIMS_random_effects",
          TRUE ~ "FIMS"
        )
      ) |>
      dplyr::select(sim, model, positive_hessian, gradient)
    
    results <- dplyr::bind_rows(results, fims_results)
  }
  
  # Make a table of convergence rates for each estimation model
  convergence_rates <- results |>
    dplyr::mutate(
      model = factor(
        model, 
        levels = c(
          "ASAP", "BAM", "SS3", "WHAM_fixed_effects", "WHAM_random_effects", 
          "FIMS_fixed_effects", "FIMS_random_effects"
        )
      )
    ) |>
    dplyr::group_by(model) |>
    dplyr::summarize(
      total_simulations = dplyr::n(),
      converged_simulations = sum(positive_hessian == TRUE & gradient <= ifelse(is.null(names(gradient_threshold)), gradient_threshold, gradient_threshold[model]), na.rm = TRUE),
      convergence_rate = converged_simulations / total_simulations * 100,
      .groups = "drop"
    ) |>
    dplyr::arrange(model)
  
  # Save convergence rates to a CSV file for reporting
  write.csv(convergence_rates, file.path(case_dir, "figure", "convergence_rates.csv"))
  
  # Handle model-specific gradient thresholds
  if (length(gradient_threshold) == 1 && is.null(names(gradient_threshold))) {
    threshold_df <- tibble::tibble(model = unique(results[["model"]]), threshold = gradient_threshold)
  } else if (!is.null(names(gradient_threshold))) {
    threshold_df <- tibble::tibble(model = names(gradient_threshold), threshold = unname(gradient_threshold))
    
    # Add default 0.01 for any models present in results but missing from gradient_threshold
    missing_models <- setdiff(unique(results[["model"]]), threshold_df[["model"]])
    if (length(missing_models) > 0) {
      threshold_df <- dplyr::bind_rows(
        threshold_df,
        tibble::tibble(model = missing_models, threshold = 0.01)
      )
    }
  } else {
    stop("gradient_threshold must be a single numeric value or a named numeric vector.")
  }
  
  # Identify converged simulations for each model
  # Convergence criteria: positive_hessian = TRUE AND gradient <= threshold
  # Group by model to get a list of converged simulation IDs for each
  converged_by_model <- results |> 
    dplyr::left_join(threshold_df, by = "model") |>
    dplyr::filter(positive_hessian == TRUE & gradient <= threshold) |>
    dplyr::group_by(model) |> 
    dplyr::summarise(converged_sims = list(sim), .groups = "drop")
  
  # Find the intersection of converged simulations across ALL models
  # Only simulations that converged for every model will be in common_sim_ids
  common_sim_ids <- Reduce(intersect, converged_by_model[["converged_sims"]])
}
