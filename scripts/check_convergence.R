#' Check convergence for ASAP, BAM, SS3, WHAM, and FIMS estimation models
#' 
#' Evaluates convergence based on two criteria: a positive definite Hessian and
#' a maximum gradient component below a specified threshold.
#' 
#' @param em_names Character vector of model names (e.g., \code{c("ASAP", "BAM", "SS", "WHAM", "FIMS")}).
#' @param n_sim Integer, number of operating model simulations.
#' @param case_dir Character, path to case working directory.
#' @param gradient_threshold Numeric. A single value or a named vector 
#'   (e.g., \code{c(ASAP=0.001, WHAM_fixed_effects=1e-4)}). Defaults to 0.01.
#' 
#' @return A vector of simulation IDs that converged for ALL estimation models.
#' 
#' @examples
#' \dontrun{
#' converged_ids <- check_convergence(
#'   em_names = c("ASAP", "WHAM"),
#'   n_sim = 100,
#'   case_dir = "C0",
#'   gradient_threshold = 0.01
#' )
#' 
#' converged_ids <- check_convergence(
#'   em_names = c("ASAP", "WHAM"),
#'   n_sim = 100,
#'   case_dir = "C0",
#'   gradient_threshold = c("ASAP" = 0.005, "WHAM_fixed_effects" = 0.0001)
#' )
#' }
check_convergence <- function(em_names, n_sim, case_dir, gradient_threshold = 0.01) {
  
  # Internal Helper: Safe RDS Reader
  # Centralizes the logic for reading RDS files and extracting model-specific metrics
  get_rds_value <- function(file_path, target_type = "direct", is_logical_type = FALSE) {
    if (!file.exists(file_path)) {
      return(if (is_logical_type) NA else NA_real_)
    }
    
    rds_content <- tryCatch(readRDS(file_path), error = function(e) NULL)
    if (is.null(rds_content)) return(if (is_logical_type) NA else NA_real_)
    
    # Model-specific extraction logic
    if (target_type == "wham_gradient") {
      # Use the absolute maximum of the final gradient vector
      return(max(abs(rds_content[["final_gradient"]]), na.rm = TRUE))
    }
    if (target_type == "wham_hessian") {
      # In WHAM, a successful Hessian is indicated by na_sdrep being FALSE
      return(!rds_content[["na_sdrep"]])
    }
    
    # Default for FIMS results or direct logical/numeric flags
    return(rds_content)
  }
  
  # Process ADMB Models (ASAP, BAM, SS3)
  admb_models <- intersect(em_names, c("ASAP", "BAM", "SS"))
  admb_results <- tibble::tibble()
  
  if (length(admb_models) > 0) {
    admb_results <- expand.grid(
      simulation = seq_len(n_sim), 
      model = admb_models, 
      stringsAsFactors = FALSE
    ) |>
      dplyr::mutate(
        simulation_directory = file.path(case_dir, "output", model, paste0("s", simulation)),
        hessian_path = file.path(simulation_directory, "admodel.cov"),
        parameter_path = dplyr::case_when(
          model == "ASAP" ~ file.path(simulation_directory, "asap3.par"),
          model == "BAM"  ~ file.path(simulation_directory, "BAM-Sim.par"),
          model == "SS"   ~ file.path(simulation_directory, "ss.par"),
          TRUE ~ NA_character_
        ),
        positive_hessian = file.exists(hessian_path),
        gradient = purrr::map_dbl(parameter_path, ~ {
          if (file.exists(.x)) {
            # Extract the maximum gradient (16th value) from the .par file header
            as.numeric(scan(.x, what = '', n = 16, quiet = TRUE)[c(6, 11, 16)])[3]
          } else NA_real_
        })
      ) |>
      dplyr::select(simulation, model, positive_hessian, gradient)
  }
  
  # Process TMB Models (WHAM & FIMS)
  # These models use a dual Fixed/Random effects structure
  tmb_models <- intersect(em_names, c("WHAM", "FIMS"))
  tmb_results <- tibble::tibble()
  
  if (length(tmb_models) > 0) {
    tmb_results <- expand.grid(
      simulation = seq_len(n_sim), 
      model = tmb_models, 
      stringsAsFactors = FALSE
    ) |>
      dplyr::mutate(
        simulation_directory = file.path(case_dir, "output", model, paste0("s", simulation)),
        # Map paths based on model naming conventions
        fixed_path  = ifelse(
          model == "WHAM", 
          file.path(simulation_directory, "fit_wham_fixed_effects.RDS"), 
          file.path(simulation_directory, "max_gradient_fims_fixed_effects.RDS")
        ),
        random_path = ifelse(
          model == "WHAM", 
          file.path(simulation_directory, "fit_wham_random_effects.RDS"), 
          file.path(simulation_directory, "max_gradient_fims_random_effects.RDS")
        ),
        # FIMS has separate Hessian files; WHAM checks the fit object itself
        fixed_hessian_path  = ifelse(
          model == "WHAM", 
          fixed_path, 
          file.path(simulation_directory, "hessian_fims_fixed_effects.RDS")
        ),
        random_hessian_path = ifelse(
          model == "WHAM", 
          random_path, 
          file.path(simulation_directory, "hessian_fims_random_effects.RDS")
        )
      ) |>
      dplyr::mutate(
        gradient_fixed = purrr::map_dbl(fixed_path, ~ {
          if(grepl("wham", .x)) {
            get_rds_value(.x, "wham_gradient")
          } else {
            get_rds_value(.x)
          }
        }
        ),
        gradient_random = purrr::map_dbl(random_path, ~ {
          if(grepl("wham", .x)) {
            get_rds_value(.x, "wham_gradient")
          } else {
            get_rds_value(.x)
          }
        }
        ),
        hessian_fixed   = purrr::map_lgl(fixed_hessian_path, ~ {
          if(grepl("wham", .x)) {
            get_rds_value(.x, "wham_hessian", TRUE)
          } else {
            get_rds_value(.x, is_logical_type = TRUE)
          }
        }
        ),
        hessian_random  = purrr::map_lgl(random_hessian_path, ~ {
          if(grepl("wham", .x)) {
            get_rds_value(.x, "wham_hessian", TRUE)
          } else {
            get_rds_value(.x, is_logical_type = TRUE)
          }
        }
        )
      ) |>
      # Pivot using a pattern
      tidyr::pivot_longer(
        cols = matches("gradient|hessian"), 
        names_to = c(".value", "effect_type"), 
        names_pattern = "(gradient|hessian)_(fixed|random)"
      ) |>
      dplyr::mutate(model = paste0(model, "_", effect_type, "_effects")) |>
      dplyr::select(
        simulation, model, positive_hessian = hessian, gradient = gradient
      )
  }
  
  # Threshold Evaluation & Intersection
  all_results <- dplyr::bind_rows(admb_results, tmb_results)
  unique_models_found <- unique(all_results$model)
  
  # Standardize thresholds into a lookup table
  if (is.null(names(gradient_threshold))) {
    threshold_lookup <- tibble::tibble(
      model = unique_models_found, threshold_value = gradient_threshold
    )
  } else {
    threshold_lookup <- tibble::tibble(
      model = names(gradient_threshold), 
      threshold_value = unname(gradient_threshold)
    ) |>
      dplyr::full_join(tibble::tibble(model = unique_models_found), by = "model") |>
      dplyr::mutate(threshold_value = tidyr::replace_na(threshold_value, 0.001))
  }
  
  # Return IDs that converged across every model type requested
  converged_simulation_ids <- all_results |>
    dplyr::left_join(threshold_lookup, by = "model") |>
    dplyr::filter(positive_hessian == TRUE, gradient <= threshold_value) |>
    dplyr::group_by(model) |>
    dplyr::summarize(valid_simulations = list(simulation), .groups = "drop") |>
    purrr::pluck("valid_simulations") |>
    Reduce(f = intersect)
}
