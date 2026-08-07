# Use the sequence of shell commands to install Wine to run .exe on Ubuntu
install_wine_on_linux <- function() {
  # Check if the operating system is Linux
  if (Sys.info()["sysname"] != "Linux") {
    message("This script is for Linux only. Skipping.")
    return()
  }
  
  # These are the commands to run, in order.
  # They are specific to Ubuntu 22.04 (Jammy).
  commands <- c(
    "sudo wget -O /etc/apt/keyrings/winehq-archive.key https://dl.winehq.org/wine-builds/winehq.key",
    "sudo wget -NP /etc/apt/sources.list.d/ https://dl.winehq.org/wine-builds/ubuntu/dists/jammy/winehq-jammy.sources",
    "sudo dpkg --add-architecture i386",
    "sudo apt update",
    "sudo apt-get install -f -y",
    "sudo apt install --install-recommends winehq-stable -y"
  )
  
  # Loop through each command and run it
  for (cmd in commands) {
    cli::cli_text(paste("Running:", cmd))
    status <- system(cmd) # This runs the command

    # Check if the command failed (status != 0)
    if (status != 0) {
      cli::cli_abort("A command failed. Stopping the script.")
      return() # Stop the function
    }
  }

  cli::cli_text("All commands finished successfully.")
}

# Install R packages that are required for the project
install_required_packages <- function() {
  # Install the required R packages
  install.packages("pak")
  pak::pkg_install(
    c(
      # Requires run_rceattle() on this branch. Until that PR merges, run_em()
      # will skip "Rceattle" silently and its panels will be absent from every
      # figure.
      "NOAA-FIMS/Age_Structured_Stock_Assessment_Model_Comparison@update-ems",
      "NOAA-FIMS/FIMS",
      "timjmiller/wham",
      # Rceattle is a TMB model, so installing it compiles C++ and needs a
      # working toolchain. Pinned to a release tag rather than a branch: a
      # branch moves, and the manuscript's numbers must not move with it.
      # 4.9.1 is commit 7d873ca7.
      "grantdadams/Rceattle@4.9.1",
      # TMBhelper is only in Rceattle's Suggests, so pak will NOT pull it in --
      # and Rceattle takes a different optimizer path when it is absent
      # (R/0-tmb_helpers.R), which would make fits depend on what happens to be
      # installed. Request it explicitly.
      "kaskr/TMB_contrib_R/TMBhelper",
      # NOTE: r4ss is not listed here, so read_output_data() cannot read the SS3
      # output and SS3 is absent from the figures. It now warns rather than
      # skipping quietly. Adding r4ss changes which simulations survive the
      # shared outlier filter, so it is left as a separate decision.
      "httr",
      "parallel",
      "doParallel",
      "googledrive",
      "gargle"
    ), 
    ask = FALSE
  )
  
  # Download the EM input files required for the model comparison project, including
  # source code .tpl files for model compilation, input data files, and configuration
  # files.
  # Create the API URL
  api_url <- "https://api.github.com/repos/NOAA-FIMS/Age_Structured_Stock_Assessment_Model_Comparison/contents/example/em_input"
  # Make the API request
  response <- httr::GET(api_url)
  # Parse the JSON response into a data frame
  file_list <- jsonlite::fromJSON(httr::content(response, "text"), flatten = TRUE)
  # Create a local directory to save the files
  main_dir <- getwd()
  if (!dir.exists(main_dir)) {
    dir.create(main_dir)
  }
  local_dir <- file.path(main_dir, "em_input")
  if (!dir.exists(local_dir)) {
    dir.create(local_dir)
  }
  # Filter for files only (type == "file") and get their download URLs
  file_urls <- subset(file_list, type == "file")[["download_url"]]
  # Loop through the URLs and download each file
  for (url in file_urls) {
    # Define the local path to save the file
    file_name <- basename(url)
    destination_path <- file.path(local_dir, file_name)
    # Download the file
    download.file(url, destfile = destination_path, mode = "wb")
  }
  # Download Stock Synthesis 3 executable file
  exe_urls <- c(
    "https://github.com/Bai-Li-NOAA/Model_Comparison_Paper/raw/refs/heads/master/em/em_input_raw/amak.exe",
    "https://github.com/Bai-Li-NOAA/Model_Comparison_Paper/raw/refs/heads/master/em/em_input_raw/ASAP3.exe",
    "https://github.com/Bai-Li-NOAA/Model_Comparison_Paper/raw/refs/heads/master/em/em_input_raw/BAM-Sim_C0.exe",
    "https://github.com/Bai-Li-NOAA/Model_Comparison_Paper/raw/refs/heads/master/em/em_input_raw/BAM-Sim_C8.exe",
    "https://github.com/Bai-Li-NOAA/Model_Comparison_Paper/raw/refs/heads/master/em/em_input_raw/BAM-Sim_C9.exe",
    "https://github.com/Bai-Li-NOAA/Model_Comparison_Paper/raw/refs/heads/master/em/em_input_raw/BAM-Sim_C13.exe",
    "https://github.com/Bai-Li-NOAA/Model_Comparison_Paper/raw/refs/heads/master/em/em_input_raw/ss.exe"
  )
  for (url in exe_urls) {
    # Define the local path to save the file
    file_name <- basename(url)
    destination_path <- file.path(local_dir, file_name)
    # Download the file
    download.file(url, destfile = destination_path, mode = "wb")
  }
}
