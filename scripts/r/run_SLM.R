# =============================================================================
# SLM SIMULATION RUNNER SCRIPT
# =============================================================================
# This script orchestrates parallel execution of SLM river network model simulations
# Handles configuration loading, parameter sweeps, parallel processing, and results aggregation

# Clear workspace to ensure clean simulation environment
rm(list = ls())

# =============================================================================
# CONFIGURATION FILE SETUP
# =============================================================================
# Configuration file contains all simulation parameters and system paths
# Can be provided via command line argument or set manually for development

config_file = "path\\to\\project\\simulations\\config_monte-carlo_local_test_just1.json"

# Verify working directory for relative path resolution
getwd()

# =============================================================================
# LIBRARY LOADING AND DEPENDENCIES
# =============================================================================
# Load all required packages for SLM simulation system

# Core river network modeling and spatial analysis
library(OCNet)          # Optimal Channel Networks for realistic river network generation
library(reticulate)     # Python integration for computationally intensive subcatchment modeling
library(fields)         # Spatial statistics and interpolation functions

# Data manipulation and visualization ecosystem  
library(tidyverse)      # Data wrangling, transformation, and piping operations
library(ggplot2)        # Advanced plotting and visualization (part of tidyverse but explicit)
library(glue)           # String interpolation and formatting

# Statistical analysis and time series processing
library(synchrony)      # Spatial-temporal synchrony analysis for ecological time series
library(matrixStats)    # Fast matrix statistical operations

# High-performance computing and parallel processing
library(doFuture)       # Parallel foreach backend using future framework
library(progressr)      # Progress reporting for long-running parallel computations
library(tictoc)         # Timing and performance monitoring

# Data storage and configuration management
library(arrow)          # High-performance columnar data storage (parquet format)
library(jsonlite)       # JSON configuration file parsing

# Configure progress bar appearance with custom styling
handlers(handler_txtprogressbar(char = cli::col_red(cli::symbol$heart)))

# =============================================================================
# CONFIGURATION FILE PROCESSING
# =============================================================================
# Parse JSON configuration file containing simulation parameters and system settings

# Load and parse JSON configuration file
# Configuration contains parameter sweeps, network specifications, and file paths
config_data <- fromJSON(config_file)

# =============================================================================
# PARAMETER SWEEP EXTRACTION
# =============================================================================
# Extract parameter combinations for Monte Carlo or sensitivity analysis
# Each row represents one complete simulation with unique parameter combination

# Extract parameter sets array from configuration
# This defines the parameter space to be explored (e.g., Monte Carlo sampling)
param_sets <- config_data$parameter_sets

# Convert parameter sets to dataframe for easy iteration and indexing
# Each column represents a model parameter, each row a simulation instance
param_df <- as.data.frame(param_sets)

# =============================================================================
# GENERAL SIMULATION SETTINGS
# =============================================================================
# Extract general configuration parameters that apply to all simulations

# River network specification (e.g., network topology, subcatchment threshold)
network <- config_data$general$network

# Spatial analysis grouping method (e.g., "order-based", "upstream-area-based")  
set_build <- config_data$general$set_build

# =============================================================================
# SYSTEM PATH CONFIGURATION
# =============================================================================
# Extract file paths for executables, data, and output directories
# These paths enable cross-platform compatibility and flexible deployment

# Python interpreter path for subcatchment-scale modeling
path_to_python <- config_data$system$path_to_python

# Python functions file containing computational modules
path_to_python_functions <- config_data$system$path_to_python_functions

# R functions file containing SLM model implementation
path_to_R_functions <- config_data$system$path_to_R_functions

# Optimal Channel Network (OCN) data file location
path_to_OCN <- config_data$system$path_to_OCN

# Output directory for simulation results and analysis products
path_to_output <- config_data$system$path_to_output

# =============================================================================
# FUNCTION LOADING AND ENVIRONMENT SETUP  
# =============================================================================
# Load SLM model functions and configure execution environment

# Load SLM model implementation from functions file
# This includes all core modeling functions (routing, mixing, analysis, etc.)
source(path_to_R_functions)

# =============================================================================
# EXECUTION ENVIRONMENT CONFIGURATION
# =============================================================================
# Configure R environment settings for simulation execution

# Set warning reporting to immediate display (warn = 1)
# This ensures warnings are visible during parallel execution
options(warn = 1)

# Enable progress bar reporting for long-running computations
# Progress feedback is essential for monitoring large parameter sweeps
options(progress_enabled = TRUE)

# =============================================================================
# PARALLEL-SAFE LOGGING SYSTEM
# =============================================================================
# Define logging function that handles concurrent writes from multiple workers

# Custom logging function for parallel processing environments
# Safely captures warnings from multiple parallel workers to single log file
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

# =============================================================================
# SIMULATION BATCH OVERVIEW
# =============================================================================
# Display simulation batch information for user confirmation

# Report total number of simulations to be executed
# This helps users verify parameter sweep size before execution
print(glue("number of simulations: {nrow(param_sets)}"))

# =============================================================================
# PARALLEL PROCESSING INFRASTRUCTURE SETUP
# =============================================================================
# Configure parallel computing environment for efficient parameter sweep execution

# =============================================================================
# PARALLEL EXECUTION PLAN CONFIGURATION
# =============================================================================
# Set up parallel processing backend using future framework
# Multisession provides true parallel execution across CPU cores

# Configure parallel execution plan (currently set to 1 worker for debugging)
# Increase workers count for production runs based on available CPU cores
plan(multisession, workers = 1)

# Enable strict global variable checking to prevent parallel execution errors
# This catches variable scoping issues that could cause silent failures
options(future.globals.onReference = "error")

# =============================================================================
# OUTPUT DIRECTORY AND FILE MANAGEMENT
# =============================================================================
# Prepare output infrastructure for simulation results storage

# Create parquet subdirectory for efficient columnar data storage
# Parquet format enables fast I/O for large time series datasets
dir.create(file.path(path_to_output, "data"), recursive = TRUE)

# Initialize log file for tracking simulation progress and warnings
# Network identifier in filename enables parallel execution of different networks
log_file <- file.path(path_to_output, glue("log_SLM_{network}.log"))
file.create(log_file)

# =============================================================================
# MAIN SIMULATION BATCH EXECUTION
# =============================================================================
# Execute parameter sweep using parallel foreach loop with progress monitoring

# Start timing for performance monitoring
time_elapsed <- tic(msg = "parallel")

# Execute simulation batch with progress reporting
with_progress({
  # Initialize progress reporter for user feedback
  # Steps = total number of parameter combinations to simulate
  p <- progressor(
    steps = nrow(param_sets)
  )
  
  # =============================================================================
  # PARALLEL FOREACH PARAMETER SWEEP
  # =============================================================================
  # Iterate through all parameter combinations using parallel workers
  # Each parameter vector becomes a separate simulation instance
  
  foreach(
    # BIOGEOCHEMICAL PARAMETERS
    gamma_sc = param_sets$gamma_sc,                           # Subcatchment scaling exponent
    sigma_w_sc = param_sets$sigma_w_sc,                       # Subcatchment variability parameter
    damkohler_transport = param_sets$damkohler_transport,     # Transport Damkohler number
    gamma_ls = param_sets$gamma_ls,                           # Network-wide scaling exponent  
    sigma_w_ls = param_sets$sigma_w_ls,                       # Network-wide variability parameter
    
    # SIMULATION CONTROL PARAMETERS
    tmax_yrs = param_sets$tmax_yrs,                           # Total simulation time [years]
    warm_up = param_sets$warm_up,                             # Model warm-up period [days]
    overall_mean_concentration = param_sets$overall_mean_concentration, # Target mean concentration
    # HYDROCLIMATIC AND LANDSCAPE PARAMETERS  
    oro_scaling = param_sets$oro_scaling,                                 # Orographic scaling factor
    rho = param_sets$rho,                                 # Degree of concurrency of hydroclimatic forcing
    
    # DAMKOHLER NUMBERS (REACTION-TRANSPORT COUPLING)
    damkohler_interarrival = param_sets$damkohler_interarrival,          # Interarrival timescale coupling
    damkohler_longterm = param_sets$damkohler_longterm,                   # Long-term process coupling
    
    # PRECIPITATION AND EVAPOTRANSPIRATION
    rain_per_year = param_sets$rain_per_year,                            # Annual precipitation [mm/yr]
    interarrival_time_mean = param_sets$interarrival_time_mean,           # Mean precipitation interval [days]
    ET_max = param_sets$ET_max,                                           # Maximum evapotranspiration [mm/day]
    
    # SOIL AND SUBSURFACE PARAMETERS
    z_r = param_sets$z_r,                                                 # Root zone depth [mm]
    z_vz = param_sets$z_vz,                                               # Vadose zone depth [mm] 
    theta_fc = param_sets$theta_fc,                                       # Field capacity soil moisture [-]
    theta_wp = param_sets$theta_wp,                                       # Wilting point soil moisture [-]
    theta_res = param_sets$theta_res,                                     # Residual soil moisture [-]
    R = param_sets$R,                                                     # Recharge rate [mm/day]
    mean_tr = param_sets$mean_tr,                                         # Mean residence time [days]
    theta_sat = param_sets$theta_sat,                                     # Saturated soil moisture [-]
    aquifer_z = param_sets$aquifer_z,                                     # Aquifer depth [mm]
    
    # HYDRAULIC AND TRANSPORT PARAMETERS
    vf = param_sets$vf,                                                   # Settling velocity [m/s]
    mannings_n = param_sets$mannings_n,                                   # Manning's roughness coefficient [-]
    side_slope = param_sets$side_slope,                                   # Channel side slope [-]
    
    # SIMULATION IDENTIFICATION
    simulation_identifier = param_sets$simulation_identifier,             # Unique simulation ID for output files
    # =============================================================================
    # PARALLEL EXECUTION CONFIGURATION
    # =============================================================================
    .options.future = list(
      seed = TRUE,  # Enable reproducible random number generation across workers
      packages = c("OCNet", "progress", "reticulate", "tidyverse", "arrow", "matrixStats")
    ),
    # .combine = rbind,  # Disabled: simulations save results directly to files
    .errorhandling = "pass"  # Continue execution if individual simulations fail
  ) %dofuture% {
    # =============================================================================
    # INDIVIDUAL SIMULATION EXECUTION
    # =============================================================================
    # Each parallel worker executes this block with unique parameter combination
    
    # Report progress to main process for user feedback
    p()
    
    # Memory management: clean garbage collector for memory efficiency
    # rm(list = ls())  # Disabled: may interfere with parallel variable scoping
    gc()
    
    # Set seed for reproducability
    set.seed(42)
    
    # =============================================================================
    # SLM MODEL EXECUTION
    # =============================================================================
    # Call main SLM function with current parameter combination
    run_SLM(
      # Network and biogeochemical parameters
      network = network,
      gamma_sc = gamma_sc,
      sigma_w_sc = sigma_w_sc,
      damkohler_transport = damkohler_transport,
      gamma_ls = gamma_ls,
      sigma_w_ls = sigma_w_ls,
      
      # Temporal and concentration parameters  
      tmax_yrs = tmax_yrs,
      warm_up = warm_up,
      overall_mean_concentration = overall_mean_concentration,
      
      # Landscape and climate parameters
      oro_scaling = oro_scaling,
      rho = rho,
      
      # Process coupling parameters
      damkohler_interarrival = damkohler_interarrival,
      damkohler_longterm = damkohler_longterm,
      
      # Hydrological parameters
      rain_per_year = rain_per_year,
      interarrival_time_mean = interarrival_time_mean,
      ET_max = ET_max,
      
      # Soil profile parameters
      z_r = z_r,
      z_vz = z_vz,
      theta_fc = theta_fc,
      theta_wp = theta_wp,
      # Additional soil and subsurface parameters
      theta_res = theta_res,
      R = R,
      mean_tr = mean_tr,
      theta_sat = theta_sat,
      aquifer_z = aquifer_z,
      
      # Hydraulic parameters
      vf = vf,
      mannings_n = mannings_n,
      side_slope = param_sets$side_slope,
      
      # Simulation control and identification
      simulation_identifier = simulation_identifier,
      print_progress = FALSE,    # Disable individual progress bars in parallel mode

      # System configuration and logging
      log_file = log_file,
      path_to_python = path_to_python,
      path_to_python_functions = path_to_python_functions,
      path_to_OCN = path_to_OCN,
      path_to_output = path_to_output,
      log_warning_parallel = log_warning_parallel
    )
  }
})

# Report total execution time for performance monitoring
toc()

# =============================================================================
# EXECUTION MONITORING AND LOGGING
# =============================================================================
# Record execution performance and capture any warnings for debugging

# =============================================================================
# EXECUTION TIME REPORTING
# =============================================================================
# Calculate and log total execution time for performance monitoring

# Convert execution time to integer seconds
time_elapsed <- as.integer(time_elapsed)

# Convert seconds to hours for human-readable reporting
time_elapsed <- round(time_elapsed / 3600, 4)

# Log total execution time for performance analysis
log_warning_parallel(warning(glue("Ellapsed time: {time_elapsed} hours")), log_file)

# =============================================================================
# WARNING CAPTURE AND LOGGING
# =============================================================================
# Capture and log any warnings generated during execution for debugging

# Capture all accumulated warnings from R session
warning_message <- capture.output(warnings())

# Write warnings to log file for post-execution debugging and analysis
log_warning_parallel(warning(warning_message), log_file)
