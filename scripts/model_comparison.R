# Use the sequence of shell commands to install Wine to run .exe on Ubuntu
source("scripts/utils.R")
source("scripts/check_convergence.R")
source("scripts/plot_ssb_r_f.R")

install_wine_on_linux()
install_required_packages()

library(parallel)
library(doParallel)
library(FIMS)

# Save the initial OM input using ASSAMC package (sigmaR = 0.4)
model_input <- ASSAMC::save_initial_input()
maindir <- getwd()
project_dir <- getwd()
# Configure the input parameters for the simulation
sim_num <- 250
sim_input <- ASSAMC::save_initial_input(
  base_case = TRUE,
  input_list = model_input,
  maindir = maindir,
  om_sim_num = sim_num,
  keep_sim_num = sim_num,
  figure_number = 100,
  seed_num = 9924,
  case_name = "C0"
)
# sim_input$r_dev_sum2zero = TRUE
# Run OM and generate om_input, om_output, and em_input
# using function from the model comparison project
ASSAMC::run_om(input_list = sim_input)
# sim_input[["om_sim_num"]] <- 1
# sim_input[["keep_sim_num"]] <- 1
# Run EMs
ASSAMC::run_em(
  em_names = c("ASAP", "BAM", "SS", "WHAM", "FIMS"),
  input_list = sim_input,
  em_input_filenames = data.frame(
    ASAP = "C0",
    SS = "C1", 
    BAM = "C0"
  )
)

setwd(project_dir)
converged_sim_ids <- check_convergence(
  em_names = c("ASAP", "BAM", "SS", "WHAM", "FIMS"),
  n_sim = sim_input[["om_sim_num"]],
  case_dir = "C0",
  gradient_threshold = 0.004
) |>
  head(100)

results <- create_ssb_r_f_plots(
  main_dir = file.path(maindir, "C0"),
  output_dir = file.path(maindir, "C0", "figure"),
  em_names = c("ASAP", "BAM", "SS", "WHAM", "FIMS"),
  sim_ids = converged_sim_ids
)
