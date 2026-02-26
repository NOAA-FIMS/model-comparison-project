#' Tidyverse plotting functions for SSB, Recruitment, and F
#' Supports: ASAP, BAM, SS, WHAM (fixed/random), FIMS (fixed/random)
#' 
#' @author Model Comparison Project

# Model color scheme using Okabe-Ito colorblind-friendly palette
# Reference: https://jfly.uni-koeln.de/color/
model_colors <- c(
  "OM" = "#000000",           # Black (reference/truth)
  "ASAP" = "#009E73",         # Bluish green
  "BAM" = "#D55E00",          # Vermillion
  "SS" = "#0072B2",           # Blue
  "WHAM_fixed_effects" = "#CC79A7",   # Reddish purple
  "WHAM_random_effects" = "#882255",  # Dark reddish purple
  "FIMS_fixed_effects" = "#E69F00",   # Orange
  "FIMS_random_effects" = "#CC7700"   # Dark orange
)

#' Read SB, recruitment, and F data in tidy format
#'
#' @param main_dir Case directory path
#' @param em_names Estimation model names
#' @param sim_ids Simulation IDs to load (e.g., converged simulations only). If NULL, loads all available.
#' @return Tibble with columns: case, model, metric, year, simulation, value
read_ssb_r_f_data <- function(
    main_dir, 
    em_names = c("ASAP", "BAM", "SS", "WHAM", "FIMS"),
    sim_ids = NULL
) {
  
  case_name <- basename(main_dir)
  output_dir <- file.path(main_dir, "output")
  
  # Get simulation IDs if not provided
  if (is.null(sim_ids)) {
    om_files <- list.files(file.path(output_dir, "OM"), pattern = "^OM[0-9]+\\.RData$", full.names = FALSE)
    sim_ids <- as.integer(stringr::str_extract(om_files, "[0-9]+"))
    sim_ids <- sim_ids[!is.na(sim_ids)]
  }
  
  if (length(sim_ids) == 0) {
    warning(paste0("No simulation IDs provided or found in: ", file.path(output_dir, "OM")))
    return(NULL)
  }
  
  # Read OM data
  om_data <- purrr::map_dfr(sim_ids, function(sim_id) {
    env <- new.env()
    load(file.path(output_dir, "OM", paste0("OM", sim_id, ".RData")), envir = env)
    
    tibble::tibble(
      case = case_name,
      model = "OM",
      simulation = sim_id,
      year = seq_len(env[["om_input"]][["nyr"]]),
      spawning_biomass = env[["om_output"]][["SSB"]],
      recruitment = env[["om_output"]][["N.age"]][, 1] / 1000,
      fishing_mortality = apply(env[["om_output"]][["FAA"]], 1, max)
    ) |>
      tidyr::pivot_longer(
        cols = c(spawning_biomass, recruitment, fishing_mortality),
        names_to = "metric",
        values_to = "value"
      )
  })
  
  # Initialize all model data to NULL
  # asap_data <- NULL
  # bam_data <- NULL
  # ss_data <- NULL
  # wham_fixed_data <- NULL
  # wham_random_data <- NULL
  # fims_fixed_data <- NULL
  # fims_random_data <- NULL
  
  # Read ASAP data
  if ("ASAP" %in% em_names) {
    asap_data <- purrr::map_dfr(sim_ids, function(sim_id) {
      rdat_file <- file.path(output_dir, "ASAP", paste0("s", sim_id), "asap3.rdat")
      if (!file.exists(rdat_file)) return(NULL)
      
      asap_output <- dget(rdat_file)
      
      tibble::tibble(
        case = case_name,
        model = "ASAP",
        simulation = sim_id,
        year = seq_along(asap_output[["SSB"]]),
        spawning_biomass = asap_output[["SSB"]],
        recruitment = asap_output[["N.age"]][, 1],
        fishing_mortality = apply(asap_output[["fleet.FAA"]][["FAA.directed.fleet1"]], 1, max)
      ) |>
        tidyr::pivot_longer(
          cols = c(spawning_biomass, recruitment, fishing_mortality),
          names_to = "metric",
          values_to = "value"
        )
    })
  }
  
  # Read BAM data
  if ("BAM" %in% em_names) {
    bam_data <- purrr::map_dfr(sim_ids, function(sim_id) {
      rdat_file <- file.path(output_dir, "BAM", paste0("s", sim_id), "BAM-Sim.rdat")
      if (!file.exists(rdat_file)) return(NULL)
      
      bam_output <- dget(rdat_file)
      env <- new.env()
      load(file.path(output_dir, "OM", paste0("OM", sim_ids[1], ".RData")), envir = env)
      nyr <- env[["om_input"]][["nyr"]]
      
      tibble::tibble(
        case = case_name,
        model = "BAM",
        simulation = sim_id,
        year = seq_len(nyr),
        spawning_biomass = bam_output[["t.series"]][["SSB"]][1:nyr],
        recruitment = bam_output[["t.series"]][["recruits"]][1:nyr] / 1000,
        fishing_mortality = bam_output[["t.series"]][["F.full"]][1:nyr]
      ) |>
        tidyr::pivot_longer(
          cols = c(spawning_biomass, recruitment, fishing_mortality),
          names_to = "metric",
          values_to = "value"
        )
    })
  }
  
  # Read SS data
  if ("SS" %in% em_names && requireNamespace("r4ss", quietly = TRUE)) {
    ss_data <- purrr::map_dfr(sim_ids, function(sim_id) {
      ss_dir <- file.path(output_dir, "SS", paste0("s", sim_id))
      if (!dir.exists(ss_dir)) return(NULL)
      
      ss_output <- r4ss::SS_output(dir = ss_dir, ncols = 300, verbose = FALSE, printstats = FALSE)
      
      # Get year range from first OM file
      env <- new.env()
      load(file.path(output_dir, "OM", paste0("OM", sim_ids[1], ".RData")), envir = env)
      year_range <- env[["om_input"]][["year"]]
      
      ss_ts <- ss_output[["timeseries"]] |>
        dplyr::filter(Yr >= min(year_range), Yr <= max(year_range))
      
      ss_natage <- ss_output[["natage_annual_2_with_fishery"]] |>
        dplyr::filter(Yr >= min(year_range), Yr <= max(year_range))
      
      tibble::tibble(
        case = case_name,
        model = "SS",
        simulation = sim_id,
        year = seq_along(ss_ts[["SpawnBio"]]),
        spawning_biomass = ss_ts[["SpawnBio"]],
        recruitment = ss_natage[[5]],
        fishing_mortality = ss_ts[["F:_1"]]
      ) |>
        tidyr::pivot_longer(
          cols = c(spawning_biomass, recruitment, fishing_mortality),
          names_to = "metric",
          values_to = "value"
        )
    })
  }
  
  # Read WHAM data (fixed effects)
  if ("WHAM" %in% em_names) {
    wham_fixed_data <- purrr::map_dfr(sim_ids, function(sim_id) {
      rds_file <- file.path(output_dir, "WHAM", paste0("s", sim_id), "fit_wham_fixed_effects.RDS")
      if (!file.exists(rds_file)) return(NULL)
      
      fit_wham <- readRDS(rds_file)
      
      # Extract F from FAA (fishing mortality at age)
      f_mort <- exp(fit_wham[["rep"]][["log_F_tot"]])
      
      tibble::tibble(
        case = case_name,
        model = "WHAM_fixed_effects",
        simulation = sim_id,
        year = seq_along(fit_wham[["rep"]][["SSB"]]),
        spawning_biomass = fit_wham[["rep"]][["SSB"]],
        recruitment = fit_wham[["rep"]][["NAA"]][1, 1, , 1],
        fishing_mortality = f_mort
      ) |>
        tidyr::pivot_longer(
          cols = c(spawning_biomass, recruitment, fishing_mortality),
          names_to = "metric",
          values_to = "value"
        )
    })
    
    # Read WHAM data (random effects)
    wham_random_data <- purrr::map_dfr(sim_ids, function(sim_id) {
      rds_file <- file.path(output_dir, "WHAM", paste0("s", sim_id), "fit_wham_random_effects.RDS")
      if (!file.exists(rds_file)) return(NULL)
      
      fit_wham <- readRDS(rds_file)
      
      # Extract F from FAA (fishing mortality at age)
      f_mort <- exp(fit_wham[["rep"]][["log_F_tot"]])
      
      tibble::tibble(
        case = case_name,
        model = "WHAM_random_effects",
        simulation = sim_id,
        year = seq_along(fit_wham[["rep"]][["SSB"]]),
        spawning_biomass = fit_wham[["rep"]][["SSB"]],
        recruitment = fit_wham[["rep"]][["NAA"]][1, 1, , 1],
        fishing_mortality = f_mort
      ) |>
        tidyr::pivot_longer(
          cols = c(spawning_biomass, recruitment, fishing_mortality),
          names_to = "metric",
          values_to = "value"
        )
    })
  }
  
  # Read FIMS data (fixed effects)
  if ("FIMS" %in% em_names) {
    fims_fixed_data <- purrr::map_dfr(sim_ids, function(sim_id) {
      rds_file <- file.path(output_dir, "FIMS", paste0("s", sim_id), "fit_fims_fixed_effects.RDS")
      if (!file.exists(rds_file)) return(NULL)
      
      fims_output <- readRDS(rds_file)
        # TODO: remove the code below because run_fims() will return estimates instead of a full fit
        # FIMS::get_estimates()
      
      # Get year range from first OM file
      env <- new.env()
      load(file.path(output_dir, "OM", paste0("OM", sim_ids[1], ".RData")), envir = env)
      nyr <- env[["om_input"]][["nyr"]]
      
      # Extract spawning biomass
      sb_data <- fims_output |>
        dplyr::filter(label == "spawning_biomass", year_i %in% 1:nyr) |>
        dplyr::arrange(year_i) |>
        dplyr::pull(estimated)
      
      # Extract recruitment (first age)
      recruit_data <- fims_output |>
        dplyr::filter(label == "expected_recruitment", year_i %in% 1:nyr) |>
        dplyr::arrange(year_i) |>
        dplyr::pull(estimated) / 1000
      
      # Extract fishing mortality
      f_data <- fims_output |>
        dplyr::filter(label == "log_Fmort", year_i %in% 1:nyr, module_id == 1) |>
        dplyr::arrange(year_i) |>
        dplyr::pull(estimated) |>
        exp()
      
      tibble::tibble(
        case = case_name,
        model = "FIMS_fixed_effects",
        simulation = sim_id,
        year = 1:nyr,
        spawning_biomass = sb_data,
        recruitment = recruit_data,
        fishing_mortality = f_data
      ) |>
        tidyr::pivot_longer(
          cols = c(spawning_biomass, recruitment, fishing_mortality),
          names_to = "metric",
          values_to = "value"
        )
    })
    
    # Read FIMS data (random effects)
    fims_random_data <- purrr::map_dfr(sim_ids, function(sim_id) {
      rds_file <- file.path(output_dir, "FIMS", paste0("s", sim_id), "fit_fims_random_effects.RDS")
      if (!file.exists(rds_file)) return(NULL)
      
      fims_output <- readRDS(rds_file) |>
        purrr::pluck("sdr_report") |>
        tibble::as_tibble(rownames = "label")
      
      # Get year range from first OM file
      env <- new.env()
      load(file.path(output_dir, "OM", paste0("OM", sim_ids[1], ".RData")), envir = env)
      nyr <- env[["om_input"]][["nyr"]]
      
      # Extract spawning biomass
      sb_data <- fims_output |>
        dplyr::filter(label == "spawning_biomass") |>
        dplyr::slice(1:nyr) |>
        dplyr::pull(Estimate)
      
      # Extract recruitment (first age)
      recruit_data <- fims_output |>
        dplyr::filter(label == "expected_recruitment") |>
        dplyr::slice(1:nyr) |>
        dplyr::pull(Estimate) / 1000
      
      # Extract fishing mortality
      f_data <- fims_output |>
        dplyr::filter(label == "log_Fmort") |>
        dplyr::slice(1:nyr) |>
        dplyr::pull(Estimate) |>
        exp()
      
      tibble::tibble(
        case = case_name,
        model = "FIMS_random_effects",
        simulation = sim_id,
        year = 1:nyr,
        spawning_biomass = sb_data,
        recruitment = recruit_data,
        fishing_mortality = f_data
      ) |>
        tidyr::pivot_longer(
          cols = c(spawning_biomass, recruitment, fishing_mortality),
          names_to = "metric",
          values_to = "value"
        )
    })
  }
  
  # Combine all data
  dplyr::bind_rows(
    om_data,
    asap_data,
    bam_data,
    ss_data,
    wham_fixed_data,
    wham_random_data,
    fims_fixed_data,
    fims_random_data
  )
}

#' Calculate 95% confidence interval for median using binomial distribution
#'
#' @param x Numeric vector
#' @return Named vector with lower and upper bounds
#' @references Thompson (1936). Ann. Math. Stat. 7, 122-128. doi:10.1214/aoms/1177732502
median_ci_binomial <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  
  if (n < 2) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  
  # Calculate order statistics for 95% CI using binomial distribution
  k_lower <- stats::qbinom(0.025, n, 0.5)
  k_upper <- stats::qbinom(0.975, n, 0.5)
  
  sorted_x <- sort(x)
  
  # Ensure valid indices (1-based indexing in R)
  k_lower <- max(1, k_lower + 1)
  k_upper <- min(n, k_upper + 1)
  
  return(c(
    lower = sorted_x[k_lower],
    upper = sorted_x[k_upper]
  ))
}

#' Calculate summary statistics with uncertainty
#'
#' @param data Tidy data frame
#' @return Tibble with median, 95% CI for median, and percentiles
calc_summary_stats <- function(data) {
  data |>
    dplyr::group_by(case, model, metric, year) |>
    dplyr::summarize(
      median = stats::median(value, na.rm = TRUE),
      # 95% CI for median using binomial distribution (Thompson, 1936)
      median_lower = median_ci_binomial(value)["lower"],
      median_upper = median_ci_binomial(value)["upper"],
      # Use 5th and 95th percentiles for uncertainty intervals
      lower = stats::quantile(value, 0.05, na.rm = TRUE),
      upper = stats::quantile(value, 0.95, na.rm = TRUE),
      # Also get 25th and 75th percentiles for potential use in boxplots or ribbons
      q25 = stats::quantile(value, 0.25, na.rm = TRUE),
      q75 = stats::quantile(value, 0.75, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    )
}

#' Plot SSB, recruitment, or F with uncertainty
#'
#' @param summary_data Summary statistics data
#' @param metric_name Metric to plot ("spawning_biomass", "recruitment", or "fishing_mortality")
#' @param facet_by Facet by "model", "case", or NULL
#' @param include_om Include OM as background
#' @return ggplot object
plot_metric <- function(summary_data, 
                        metric_name = "spawning_biomass",
                        facet_by = "model",
                        include_om = TRUE) {
  
  # Filter to metric
  plot_data <- summary_data |> dplyr::filter(metric == metric_name)
  
  # Define model ordering
  model_order <- c("ASAP", "BAM", "SS", "WHAM_fixed_effects", "WHAM_random_effects", "FIMS_fixed_effects", "FIMS_random_effects")
  
  # Separate OM
  if (include_om && !is.null(facet_by)) {
    # Get OM data without model column for replication
    om_data_template <- plot_data |> 
      dplyr::filter(model == "OM") |>
      dplyr::select(-model)
    
    # Get EM data with proper factor ordering
    em_data <- plot_data |> 
      dplyr::filter(model != "OM") |>
      dplyr::mutate(model = factor(model, levels = model_order))
    
    # Replicate OM data for each EM model to show in all facets as background
    em_models <- levels(em_data$model)
    om_data <- purrr::map_dfr(em_models, ~{
      om_data_template |>
        dplyr::mutate(model = factor(.x, levels = model_order))
    })
  } else {
    em_data <- plot_data |>
      dplyr::mutate(model = factor(model, levels = c("OM", model_order)))
    om_data <- NULL
  }
  
  # Axis labels
  y_labels <- c(
    spawning_biomass = "Spawning Stock Biomass (mt)",
    recruitment = "Recruitment (×1000 fish)",
    fishing_mortality = "Fishing Mortality (F)"
  )
  
  y_label <- y_labels[[metric_name]]
  
  # Create plot
  p <- ggplot2::ggplot(em_data, ggplot2::aes(x = year, y = median, color = model, fill = model))
  
  # Add OM background
  if (!is.null(om_data) && nrow(om_data) > 0) {
    p <- p +
      ggplot2::geom_ribbon(
        data = om_data,
        ggplot2::aes(ymin = lower, ymax = upper),
        alpha = 0.2,
        color = NA,
        fill = model_colors[["OM"]]
      ) +
      ggplot2::geom_line(
        data = om_data,
        ggplot2::aes(y = median),
        color = model_colors[["OM"]],
        linewidth = 0.7
      )
  }
  
  # Add EM data
  p <- p +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_color_manual(values = model_colors, drop = FALSE) +
    ggplot2::scale_fill_manual(values = model_colors, drop = FALSE) +
    ggplot2::labs(
      x = "Year",
      y = y_label,
      color = "Model",
      fill = "Model"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "white"),
      strip.text = ggplot2::element_text(face = "bold", size = 10)
    )
  
  # Add faceting
  if (!is.null(facet_by)) {
    if (facet_by == "model") {
      p <- p + ggplot2::facet_wrap(~model, ncol = 2, scales = "free_y")
    } else if (facet_by == "case") {
      p <- p + ggplot2::facet_wrap(~case, ncol = 3, scales = "free_y")
    }
  }
  
  return(p)
}

#' Calculate relative errors from model output data
#'
#' @param all_data Tidy data frame with OM and EM outputs
#' @return Tibble with relative errors (RE) and absolute relative errors (ARE)
#' @details Calculates RE = (EM - OM) / OM * 100 for each metric, year, simulation
calculate_relative_errors <- function(all_data) {
  
  if (is.null(all_data) || nrow(all_data) == 0) {
    warning("No data provided for relative error calculation")
    return(NULL)
  }
  
  # Separate OM and EM data
  om_data <- all_data |>
    dplyr::filter(model == "OM") |>
    dplyr::select(case, metric, year, simulation, om_value = value)
  
  em_data <- all_data |>
    dplyr::filter(model != "OM")
  
  # Join EM with OM and calculate relative errors
  re_data <- em_data |>
    dplyr::left_join(
      om_data,
      by = c("case", "metric", "year", "simulation")
    ) |>
    dplyr::mutate(
      # Relative Error: (EM - OM) / OM * 100
      re = (value - om_value) / om_value * 100,
      # Absolute Relative Error
      are = abs(re)
    ) |>
    dplyr::select(case, model, metric, year, simulation, re, are)
  
  return(re_data)
}

#' Calculate median relative error summary table
#'
#' @param re_data Relative error data from calculate_relative_errors()
#' @return Tibble with mean of median RE by case, model, and metric
#' @details For each simulation, calculates median RE across years, then takes
#'   mean of medians across simulations with 95% CI (Thompson, 1936)
calculate_mre_summary <- function(re_data) {
  
  if (is.null(re_data) || nrow(re_data) == 0) {
    return(NULL)
  }
  
  # Step 1: Calculate median RE across years for each simulation
  sim_medians <- re_data |>
    dplyr::group_by(case, model, metric, simulation) |>
    dplyr::summarize(
      median_re = stats::median(re, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Step 2: Calculate mean of medians across simulations with 95% CI
  mre_summary <- sim_medians |>
    dplyr::group_by(case, model, metric) |>
    dplyr::summarize(
      mre = mean(median_re, na.rm = TRUE),
      mre_lower = median_ci_binomial(median_re)["lower"],
      mre_upper = median_ci_binomial(median_re)["upper"],
      sd_mre = stats::sd(median_re, na.rm = TRUE),
      n_sims = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      model = factor(model, levels = c("ASAP", "BAM", "SS", "WHAM_fixed_effects", "WHAM_random_effects", "FIMS_fixed_effects", "FIMS_random_effects")),
      mre = round(mre, 2),
      mre_lower = round(mre_lower, 2),
      mre_upper = round(mre_upper, 2),
      sd_mre = round(sd_mre, 2)
    ) |>
    dplyr::arrange(case, metric, model)
  
  return(mre_summary)
}

#' Plot relative error boxplots over time by model
#'
#' @param re_data Relative error data (with RE in percentage)
#' @param metric_filter Metric to plot
#' @return ggplot object
#' @details Converts RE from percentage to ratio: (EM - OM) / OM. Facets by year with models on x-axis.
plot_re_boxplot_by_year <- function(re_data, metric_filter = "spawning_biomass") {
  
  # Define model ordering
  model_order <- c("ASAP", "BAM", "SS", "WHAM_fixed_effects", "WHAM_random_effects", "FIMS_fixed_effects", "FIMS_random_effects")
  
  plot_data <- re_data |> 
    dplyr::filter(metric %in% metric_filter) |>
    dplyr::mutate(
      model = factor(model, levels = model_order),
      # Convert from percentage to ratio
      re_ratio = re / 100,
      metric_label = dplyr::case_when(
        metric == "spawning_biomass" ~ "Spawning Biomass",
        metric == "recruitment" ~ "Recruitment",
        metric == "fishing_mortality" ~ "Fishing Mortality",
        TRUE ~ metric
      )
    )
  
  ggplot2::ggplot(plot_data, ggplot2::aes(x = model, y = re_ratio, fill = model)) +
    ggplot2::geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
    ggplot2::facet_wrap(~year, ncol = 5) +
    ggplot2::scale_fill_manual(values = model_colors) +
    ggplot2::coord_cartesian(ylim = c(-0.4, 0.4)) +
    ggplot2::labs(
      x = "Model",
      y = "Relative Error",
      title = paste0("Relative Error Over Time: ", unique(plot_data$metric_label))
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "white"),
      strip.text = ggplot2::element_text(face = "bold", size = 9),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7)
    )
}

#' Plot relative error distributions
#'
#' @param re_data Relative error data
#' @param metric_filter Metric to plot
#' @return ggplot object
plot_re_violin <- function(re_data, metric_filter = "spawning_biomass") {
  
  # Define model ordering
  model_order <- c("ASAP", "BAM", "SS", "WHAM_fixed_effects", "WHAM_random_effects", "FIMS_fixed_effects", "FIMS_random_effects")
  
  plot_data <- re_data |> dplyr::filter(metric %in% metric_filter)
  
  plot_data <- plot_data |>
    dplyr::mutate(
      model = factor(model, levels = model_order),
      # Convert from percentage to ratio
      re_ratio = re / 100,
      case_label = case |>
        stringr::str_remove("^case") |>
        stringr::str_replace("^(\\d+)$", paste0("C", "\\1")),
      metric_label = dplyr::case_when(
        metric == "spawning_biomass" ~ "Spawning Biomass",
        metric == "recruitment" ~ "Recruitment",
        metric == "fishing_mortality" ~ "Fishing Mortality",
        TRUE ~ metric
      )
    )
  
  ggplot2::ggplot(plot_data, ggplot2::aes(x = model, y = re_ratio, fill = model)) +
    ggplot2::geom_violin(alpha = 0.7, quantiles = c(0.25, 0.5, 0.75)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
    ggplot2::facet_grid(metric_label ~ case_label, scales = "free_y") +
    ggplot2::scale_fill_manual(values = model_colors) +
    ggplot2::labs(
      x = "Model",
      y = "Relative Error",
      title = "Relative Error Distributions"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "gray95"),
      strip.text = ggplot2::element_text(face = "bold", size = 9),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8)
    )
}

#' Create all SSB, R, F plots for multiple cases
#'
#' @param main_dir main directory for a case
#' @param output_dir Output directory for figures and tables
#' @param em_names Estimation model names
#' @param sim_ids Simulation IDs to load (e.g., converged simulations). If NULL, loads all available.
#' @return List with data, summary, relative_errors, mre_summary, and plots
#' @details Calculates relative errors as (EM - OM) / OM * 100. Saves plots as PNG and 
#'   median relative error summary as CSV if output_dir is specified.
create_ssb_r_f_plots <- function(main_dir,
                                 output_dir = NULL,
                                 em_names = c("ASAP", "BAM", "SS", "WHAM", "FIMS"),
                                 sim_ids = NULL) {
  
  all_data <- read_ssb_r_f_data(
    main_dir = main_dir,
    em_names = em_names,
    sim_ids = sim_ids
  )
  
  if (is.null(all_data) || nrow(all_data) == 0) {
    stop("No data loaded")
  }
  
  summary_stats <- calc_summary_stats(all_data)
  
  # Calculate relative errors
  all_re <- calculate_relative_errors(all_data)
  mre_summary <- calculate_mre_summary(all_re)
  
  plots <- list()
  
  # Combined plots: All EMs vs OM with subplots by model
  # Filter data to include OM for background but exclude from facets
  plots[["sb_all_em_vs_om"]] <- summary_stats |>
    dplyr::filter(metric == "spawning_biomass") |>
    plot_metric(
      metric_name = "spawning_biomass",
      facet_by = "model",
      include_om = TRUE
    ) +
    ggplot2::labs(title = "Spawning Biomass: EMs vs OM")
  
  plots[["recruit_all_em_vs_om"]] <- summary_stats |>
    dplyr::filter(metric == "recruitment") |>
    plot_metric(
      metric_name = "recruitment",
      facet_by = "model",
      include_om = TRUE
    ) +
    ggplot2::labs(title = "Recruitment: EMs vs OM")
  
  plots[["f_all_em_vs_om"]] <- summary_stats |>
    dplyr::filter(metric == "fishing_mortality") |>
    plot_metric(
      metric_name = "fishing_mortality",
      facet_by = "model",
      include_om = TRUE
    ) +
    ggplot2::labs(title = "Fishing Mortality: EMs vs OM")
  
  # Relative error boxplot plots over time
  if (!is.null(all_re) && nrow(all_re) > 0) {
    plots[["re_boxplot_sb"]] <- plot_re_boxplot_by_year(all_re, "spawning_biomass")
    plots[["re_boxplot_recruit"]] <- plot_re_boxplot_by_year(all_re, "recruitment")
    plots[["re_boxplot_f"]] <- plot_re_boxplot_by_year(all_re, "fishing_mortality")
  }
  
  # Relative error plots
  if (!is.null(all_re) && nrow(all_re) > 0) {
    plots[["re_sb"]] <- plot_re_violin(all_re, "spawning_biomass")
    plots[["re_recruit"]] <- plot_re_violin(all_re, "recruitment")
    plots[["re_f"]] <- plot_re_violin(all_re, "fishing_mortality")
  }
  
  # Save plots and tables
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Save plots
    purrr::iwalk(plots, function(plot, name) {
      filename <- file.path(output_dir, paste0(name, ".png"))
      
      # Set dimensions based on plot type
      if (stringr::str_detect(name, "^re_boxplot")) {
        width <- 14
        height <- 12
      } else if (stringr::str_detect(name, "^re_")) {
        width <- 14
        height <- 12
      } else {
        width <- 14
        height <- 10
      }
      
      ggplot2::ggsave(
        filename = filename,
        plot = plot,
        width = width,
        height = height,
        dpi = 1200
      )
    })
    
    # Save MRE summary table
    if (!is.null(mre_summary)) {
      mre_file <- file.path(output_dir, "median_relative_error_summary.csv")
      readr::write_csv(mre_summary, mre_file)
    }
  }
  
  return(list(
    data = all_data,
    summary = summary_stats,
    relative_errors = all_re,
    mre_summary = mre_summary,
    plots = plots
  ))
}