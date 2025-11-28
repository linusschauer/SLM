# =============================================================================
# EXECUTE SLM
# =============================================================================
rm(list = ls())

config_file <- "path\\to\\project\\simulations\\config_example_small.json"

library(OCNet)
library(reticulate)
library(tidyverse)
library(glue)
library(matrixStats)
library(doFuture)
library(progressr)
library(tictoc)
library(arrow)
library(jsonlite)

# Configure progress bar appearance
handlers(handler_txtprogressbar(char = cli::col_red(cli::symbol$heart)))

config_data <- fromJSON(config_file)

# Extract parameter sets array from configuration
param_sets <- config_data$parameter_sets

# Convert parameter sets to dataframe
param_df <- as.data.frame(param_sets)

# =============================================================================
# GENERAL SIMULATION SETTINGS
# =============================================================================

network <- config_data$general$network

set_build <- config_data$general$set_build

# =============================================================================
# SYSTEM PATH CONFIGURATION
# =============================================================================

# Python interpreter path
path_to_python <- config_data$system$path_to_python

# Python functions file
path_to_python_functions <- config_data$system$path_to_python_functions

# R functions file
path_to_R_functions <- config_data$system$path_to_R_functions

# Optimal Channel Network (OCN) data file location
path_to_OCN <- config_data$system$path_to_OCN

# Output directory for simulation results
path_to_output <- config_data$system$path_to_output

# Load SLM model implementation from functions file
source(path_to_R_functions)

options(warn = 1)
options(progress_enabled = TRUE)


# Custom logging function for parallel processing environments
log_warning_parallel <- function(expr, log_file) {
  tryCatch(
    expr,
    warning = function(w) {
      # Open log file in append mode for concurrent writing
      con <- file(log_file, open = "a")
      # Ensure file is closed even if error occurs
      on.exit(close(con))
      # Write timestamped warning message to log file
      writeLines(paste("Warning: ", conditionMessage(w)), con)
    }
  )
}

print(glue("number of simulations: {nrow(param_sets)}"))

# Configure parallel execution plan
plan(multisession, workers = 2)

# Enable strict global variable checking to prevent parallel execution errors
options(future.globals.onReference = "error")

# Create parquet subdirectory
dir.create(file.path(path_to_output, "data"), recursive = TRUE)

# Initialize log file for tracking simulation progress and warnings
log_file <- file.path(path_to_output, glue("log_SLM_{network}.log"))
file.create(log_file)

# Start timing for performance monitoring
time_elapsed <- tic(msg = "parallel")

with_progress({
  p <- progressor(
    steps = nrow(param_sets)
  )

  foreach(
    gamma_sc = param_sets$gamma_sc,
    sigma_w_sc = param_sets$sigma_w_sc,
    damkohler_transport = param_sets$damkohler_transport,
    gamma_ls = param_sets$gamma_ls,
    sigma_w_ls = param_sets$sigma_w_ls,
    tmax_yrs = param_sets$tmax_yrs,
    warm_up = param_sets$warm_up,
    overall_mean_concentration = param_sets$overall_mean_concentration,
    oro_scaling = param_sets$oro_scaling,
    rho = param_sets$rho,
    damkohler_interarrival = param_sets$damkohler_interarrival,
    damkohler_longterm = param_sets$damkohler_longterm,
    rain_per_year = param_sets$rain_per_year,
    interarrival_time_mean = param_sets$interarrival_time_mean,
    ET_max = param_sets$ET_max,
    z_r = param_sets$z_r,
    z_vz = param_sets$z_vz,
    theta_fc = param_sets$theta_fc,
    theta_wp = param_sets$theta_wp,
    theta_res = param_sets$theta_res,
    R = param_sets$R,
    mean_tr = param_sets$mean_tr,
    theta_sat = param_sets$theta_sat,
    aquifer_z = param_sets$aquifer_z,
    vf = param_sets$vf,
    mannings_n = param_sets$mannings_n,
    side_slope = param_sets$side_slope,
    simulation_identifier = param_sets$simulation_identifier,
    .options.future = list(
      seed = TRUE, # Enable reproducible random number generation across workers
      packages = c("OCNet", "progress", "reticulate", "tidyverse", "arrow", "matrixStats")
    ),
    # .combine = rbind,  # Disabled: simulations save results directly to files
    .errorhandling = "pass" # Continue execution if individual simulations fail
  ) %dofuture% {
    # Report progress
    p()
    gc()
    set.seed(42)
    run_SLM(
      network = network,
      gamma_sc = gamma_sc,
      sigma_w_sc = sigma_w_sc,
      damkohler_transport = damkohler_transport,
      gamma_ls = gamma_ls,
      sigma_w_ls = sigma_w_ls,
      tmax_yrs = tmax_yrs,
      warm_up = warm_up,
      overall_mean_concentration = overall_mean_concentration,
      oro_scaling = oro_scaling,
      rho = rho,
      damkohler_interarrival = damkohler_interarrival,
      damkohler_longterm = damkohler_longterm,
      rain_per_year = rain_per_year,
      interarrival_time_mean = interarrival_time_mean,
      ET_max = ET_max,
      z_r = z_r,
      z_vz = z_vz,
      theta_fc = theta_fc,
      theta_wp = theta_wp,
      theta_res = theta_res,
      R = R,
      mean_tr = mean_tr,
      theta_sat = theta_sat,
      aquifer_z = aquifer_z,
      vf = vf,
      mannings_n = mannings_n,
      side_slope = param_sets$side_slope,
      simulation_identifier = simulation_identifier,
      print_progress = FALSE, # Disable individual progress bars in parallel mode
      log_file = log_file,
      path_to_python = path_to_python,
      path_to_python_functions = path_to_python_functions,
      path_to_OCN = path_to_OCN,
      path_to_output = path_to_output,
      log_warning_parallel = log_warning_parallel
    )
  }
})

# Report total execution time
toc()
time_elapsed <- as.integer(time_elapsed)
time_elapsed <- round(time_elapsed / 3600, 4)
log_warning_parallel(warning(glue("Ellapsed time: {time_elapsed} hours")), log_file)

# Capture all accumulated warnings from R session
warning_message <- capture.output(warnings())

# Write warnings to log file for post-execution debugging and analysis
log_warning_parallel(warning(warning_message), log_file)
