# Model Comparison Project 
# The code runs operating model (OM), and 7 estimation models (EM), and conduct
# convergence analysis as well as compare estimates.


# System and directory setup ----------------------------------------------

# Load custom helper functions and plotting utilities
source("scripts/utils.R")
source("scripts/check_convergence.R")
source("scripts/plot_ssb_r_f.R")

# Use the sequence of shell commands to install Wine to run .exe on Ubuntu
install_wine_on_linux()
# Install required R packages
install_required_packages()

library(parallel)
library(doParallel)
library(FIMS)

# Define directories
maindir <- getwd()
project_dir <- getwd()
case_name   <- "C0"
output_dir  <- file.path(project_dir, case_name)
figure_dir  <- file.path(output_dir, "figure")

# Create figure directory if it doesn't exist
if (!dir.exists(figure_dir)) dir.create(figure_dir, recursive = TRUE)

# Simulation configuration ------------------------------------------------

# Save the initial OM input using {ASSAMC} (sigmaR = 0.4)
model_input <- ASSAMC::save_initial_input()

sim_num <- 500
sim_input <- ASSAMC::save_initial_input(
  base_case = TRUE,
  input_list = model_input,
  maindir = maindir,
  om_sim_num = sim_num,
  keep_sim_num = sim_num,
  figure_number = 100,
  seed_num = 9924,
  case_name = case_name
)

# Run OM and generate om_input, om_output, and em_input
ASSAMC::run_om(input_list = sim_input)

# Run EMs
# Mapping specific config files to the respective ADMB models
ASSAMC::run_em(
  em_names = c("ASAP", "BAM", "SS", "WHAM", "FIMS"),
  input_list = sim_input,
  em_input_filenames = data.frame(
    ASAP = "C0",
    SS = "C1", 
    BAM = "C0"
  )
)

# Reset the directory
setwd(project_dir)

# Load outputs from all simulations
all_data <- read_output_data(
  main_dir = file.path(maindir, "C0"),
  em_names = c("ASAP", "BAM", "SS", "WHAM", "FIMS"),
  sim_ids = 1:sim_num
)
saveRDS(all_data, file.path(project_dir, "C0", "figure", "all_data.RDS"))

# Calculate relative errors (RE) for all key quantities:
# spawning biomass, recruitment, abundance, fishing_mortality, catchability
all_re <- calculate_relative_errors(all_data)

# Convergence analysis ----------------------------------------------------

# Convergence of Median Absolute Relative Error (MARE) with respect to number
# of iterations (i.e. 500)
# Convergence analysis for catchability
q_sim_are <- all_re |>
  dplyr::filter(metric == "catchability") |>
  dplyr::mutate(are = are/100) |>
  dplyr::group_by(case, model, simulation) |>
  dplyr::summarize(are_sim = median(are, na.rm = TRUE), .groups = "drop")

# Calculate cumulative MARE for varying number of iterations
q_convergence_results <- calc_centered_mare(q_sim_are)
q_convergence_plot <- plot_convergence_analysis(
  q_convergence_results, 
  metric_label = "Catchability"
)
ggplot2::ggsave(
  filename = file.path(project_dir, "C0", "figure", "q_convergence.png"),
  plot = q_convergence_plot,
  width = 14,
  height = 12,
  dpi = 1200
)

# convergence analysis for spawning_biomass
sb_sum_data <- all_data |>
  dplyr::filter(metric == "spawning_biomass") |>
  dplyr::group_by(case, model, simulation) |>
  dplyr::summarize(total_sb = sum(value, na.rm = TRUE), .groups = "drop")

# Get OM summed baseline
sb_om_sum <- sb_sum_data |>
  dplyr::filter(model == "OM") |>
  dplyr::select(case, simulation, om_sb = total_sb)

# Calculate absolute relative error across the summed measurements
sb_sim_are <- sb_sum_data |>
  dplyr::filter(model != "OM") |>
  dplyr::left_join(sb_om_sum, by = c("case", "simulation")) |>
  dplyr::mutate(
    re = (total_sb - om_sb) / om_sb,
    are_sim = abs(re)
  ) |>
  dplyr::select(case, model, simulation, are_sim)

# Calculate cumulative MARE for varying number of iterations
sb_convergence_results <- calc_centered_mare(sb_sim_are)
sb_convergence_plot <- plot_convergence_analysis(sb_convergence_results, metric_label = "Total Spawning Biomass")
ggplot2::ggsave(
  filename =file.path(project_dir, "C0", "figure", "sb_convergence.png"),
  plot = sb_convergence_plot,
  width = 14,
  height = 12,
  dpi = 1200
)

# Convergence check
# Filter simulations based on model maximum gradient and Hessian
converged_sims <- check_convergence(
  em_names = c("ASAP", "BAM", "SS", "WHAM", "FIMS"),
  n_sim = sim_input[["om_sim_num"]],
  case_dir = "C0",
  gradient_threshold = 0.004
)
saveRDS(converged_sims, file.path(project_dir, "C0", "figure", "converged_sims.RDS"))

# Identify bad simulations where WHAM fixed effects models exploded 
# These are outliers that likely didn't converge realistically.
bad_sim_ids <- all_re |>
  dplyr::filter(
    (metric == "spawning_biomass") & 
    (as.vector(are) > 15) & 
    (model == "WHAM_fixed_effects")
  ) |>
  dplyr::pull(simulation) |>
  unique() |>
  sort()

# Final set: Successfully converged simulations minus the outliers 
# Limit to first 100 iterations
converged_sim_ids <- converged_sims |>
  setdiff(bad_sim_ids) |>
  head(100)

# Generate final figures
results <- create_ssb_r_f_plots(
  main_dir = file.path(maindir, "C0"),
  output_dir = file.path(maindir, "C0", "figure"),
  data = all_data,
  em_names = c("ASAP", "BAM", "SS", "WHAM", "FIMS"),
  sim_ids = converged_sim_ids
)
