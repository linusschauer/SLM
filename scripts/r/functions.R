# =============================================================================
# Stochastic landscape model (SLM) Model Core Functions
# =============================================================================
#
# @title SLM Model Core Functions
# @description Core computational functions for the SLM stochastic landscape model
#              for simulating solute transport in river networks
# @details This file contains functions for:
#   - River network processing and subcatchment generation
#   - Hydraulic routing and mixing calculations  
#   - Lateral inflow simulation
#   - Results analysis and output processing
#   - Spatial and temporal variability analysis
# 
# @author Linus Schauer
# @date October 2025
# @version 1.0
#
# @note This is part of the SLM model system for analyzing solute transport
#       variability in river networks using Optimal Channel Networks (OCNs)
#
# Dependencies:
#   - OCNet package for river network generation
#   - Python integration via reticulate
#   - Various R packages for data manipulation and analysis
#
# =============================================================================

#' Create lookup dictionary mapping river network nodes to their Strahler stream orders
#' 
#' Builds a lookup table that maps each river network (RN) node to its 
#' corresponding Strahler stream order
#' 
#' @param OCN River network object from OCNet package
#' 
#' @return List with length equal to number of RN nodes, where each element contains
#'         the Strahler stream order for the corresponding node index
get_stream_order_dict <- function(OCN) {
  # Initialize empty list with length equal to total number of river network nodes
  stream_order_dict <- vector("list", OCN$RN$nNodes[1])
  
  # Iterate through each river network (RN) node to build the lookup dictionary
  for (node_index in seq_len(OCN$RN$nNodes[1])) {
    # Get the corresponding aggregated reach ID for this RN node
    reach <- OCN$RN$toAGReach[node_index]  
    
    # Extract Strahler stream order from the aggregated network and store in dictionary
    stream_order_dict[[node_index]] <- OCN$AG$streamOrder[reach]
  }
  
  # Return completed lookup dictionary
  return(stream_order_dict)
}

#' Generates matrix with downstream distances to outlet
#' 
#' @param OCN Optimal Channel Network object from OCNet package
#' 
#' @return Matrix of downstream path lengths
get_downstream_matrix <- function(OCN) {
  
  # Extract cell size from OCN for spatial discretization
  cellsize <- OCN$cellsize
  
  # Calculate grid dimensions based on network extent
  grid_width_cells <- max(OCN$FD$X) / cellsize + 1
  grid_height_cells <- max(OCN$FD$Y) / cellsize + 1
  
  # Create spatial grid matrix to store downstream path lengths to outlet
  df_downstream_path <- matrix(0, grid_width_cells, grid_height_cells)
  
  # Populate grid with downstream path lengths for generation of structured heterogeneity
  for (node in 1:OCN$RN$nNodes[1]) {
    
    # Get downstream distance from node to network outlet
    downstream_path <- OCN$RN$downstreamPathLength[node, OCN$RN$outlet]
    node_x_position <- OCN$RN$X[node]
    node_y_position <- OCN$RN$Y[node]
    
    # Convert spatial coordinates to grid indices
    grid_column_index <- (node_x_position / cellsize) + 1
    grid_row_index <- (node_y_position / cellsize) + 1
    df_downstream_path[grid_column_index, grid_row_index] <- downstream_path
    
    # Special case for outlet node: set small positive value to avoid zero division
    if (node == OCN$RN$outlet) {
      df_downstream_path[grid_column_index, grid_row_index] <- cellsize / 2
    }
  }
  
  # Return completed downstream distance matrix
  return(df_downstream_path)
}

#' Identify upstream boundary nodes for a stream order
#' 
#' Finds source nodes and confluence points that mark stream order reach beginnings
#' 
#' @param OCN Optimal Channel Network object from OCNet package
#' @param stream_order_node_list List of all nodes in the current stream order
#' @param stream_order_dict Dictionary mapping nodes to their stream orders
#' @param stream_order Current stream order being processed
#' 
#' @return List of upstream boundary nodes (sources and confluences)
get_upstream_node_list <- function(OCN, stream_order_node_list, stream_order_dict, stream_order) {
  
  # Initialize storage for upstream boundary nodes
  upstream_node_list <- list()      # Nodes at start of stream order reaches
  upstream_node_counter <- 1
  
  # Examine each node in current stream order to identify upstream boundaries
  for (node in stream_order_node_list) {
    
    # Find all nodes that flow directly into the current node
    upstream_nodes <- which(OCN$RN$downNode == node)
    
    # Case 1: Source node with no upstream connections (first-order streams only)
    if (length(upstream_nodes) == 0) {
      upstream_node_list[upstream_node_counter] <- node
      upstream_node_counter <- upstream_node_counter + 1
    }
    
    # Case 2: Confluence node where all upstream tributaries are from different stream orders
    # This marks the beginning of a new stream order reach after confluence
    upstream_stream_orders <- stream_order_dict[upstream_nodes]
    if (length(upstream_nodes) != 0) {
      if (all(upstream_stream_orders != stream_order)) {
        upstream_node_list[upstream_node_counter] <- node
        upstream_node_counter <- upstream_node_counter + 1
      }
    }
  }
  
  # Return list of identified upstream boundary nodes
  return(upstream_node_list)
}

#' Collect all nodes belonging to a specific stream order
#' 
#' Filters river network nodes by Strahler stream order classification
#' 
#' @param OCN Optimal Channel Network object from OCNet package
#' @param stream_order_dict Dictionary mapping nodes to their stream orders
#' @param stream_order Target stream order to filter for
#' 
#' @return List of node IDs belonging to the specified stream order
get_stream_order_node_list <- function(OCN, stream_order_dict, stream_order) {
  
  # Initialize storage for nodes in target stream order
  stream_order_node_counter <- 1
  stream_order_node_list <- list()
  
  # Iterate through all river network nodes to find matches
  for (node in 1:OCN$RN$nNodes[1]) {
    stream_order_node <- stream_order_dict[[node]]
    
    # Add node to list if it matches target stream order
    if (stream_order_node == stream_order) {
      stream_order_node_list[stream_order_node_counter] <- node
      stream_order_node_counter <- stream_order_node_counter + 1
    }
  }
  
  # Return list of nodes belonging to specified stream order
  return(stream_order_node_list)
}


#' Recursively find all upstream flow direction (FD) nodes
#' 
#' Traverses network upstream from given node to collect all contributing nodes
#' 
#' @param OCN River network object
#' @param FD_node Starting flow direction node ID
#' @param all_upstream_FD_nodes Accumulator vector for upstream node IDs
#' 
#' @return Vector of all upstream FD node IDs
get_all_upstream_FD_nodes <- function(OCN,
                                      FD_node,
                                      all_upstream_FD_nodes) {
  
  # ==================================================================
  # DIRECT UPSTREAM NODE IDENTIFICATION
  # ==================================================================
  
  # Find all Flow Direction (FD) nodes that drain directly into the current FD node
  # This identifies the immediate upstream neighbors in the flow direction grid
  # Each FD cell points to its downstream neighbor, so we search for cells pointing to current cell
  upstream_FD_nodes <- which(OCN$FD$downNode == FD_node)
  
  # ==================================================================
  # ACCUMULATION OF UPSTREAM NODES
  # ==================================================================
  
  # Add newly discovered upstream nodes to the accumulator list
  # This builds up the complete set of all contributing FD cells progressively
  # The accumulator carries the results through the recursive calls
  all_upstream_FD_nodes <- c(all_upstream_FD_nodes, upstream_FD_nodes)
  
  # ==================================================================
  # RECURSIVE TERMINATION CHECK
  # ==================================================================
  
  # Check if we've reached the headwater cells (no further upstream connections)
  # This serves as the termination condition for the recursive algorithm
  if (length(upstream_FD_nodes) == 0) {
    # Base case: no more upstream nodes found, return complete drainage network
    return(all_upstream_FD_nodes)
  } else {
    
    # ==================================================================
    # RECURSIVE UPSTREAM TRAVERSAL
    # ==================================================================
    
    # Continue recursive search through each upstream branch
    # This implements a depth-first search through the drainage network tree
    # Each branch is followed to its headwater source before backtracking
    for (node_index in 1:length(upstream_FD_nodes)) {
      
      # Process current upstream node as the new starting point for recursive search
      current_node <- upstream_FD_nodes[node_index]
      
      # Recursively collect all nodes upstream of the current node
      # The accumulator carries forward all previously discovered nodes
      # This ensures no nodes are lost during the recursive traversal
      all_upstream_FD_nodes <- get_all_upstream_FD_nodes(
        OCN = OCN,                              # River network object (unchanged)
        FD_node = current_node,                 # New starting point for upstream search
        all_upstream_FD_nodes = all_upstream_FD_nodes  # Accumulated results from all branches
      )
    }
  }
  
  # ==================================================================
  # RETURN COMPLETE UPSTREAM NETWORK
  # ==================================================================
  
  # Return the complete set of all upstream FD nodes that contribute flow to the original FD node
  # This represents the entire drainage area upstream of the starting point
  # Includes nodes from all tributary branches at all levels of the network hierarchy
  return(all_upstream_FD_nodes)
}

#' Generate lateral inflow nodes for OCN river network
#' 
#' Creates subcatchments above area threshold and assigns parameters to patches
#' 
#' @param OCN River network object from OCNet package
#' @param cellsize Grid cell size in meters
#' @param thresh Area threshold for subcatchment creation
#' @param landscape_configuration_matrix Matrix of landscape heterogeneity (mean immobile concentration of subcatchments)
#' @param precip_timing Precipitation timing parameter or "random"
#' @param precip_intensity Precipitation intensity parameter or "random" 
#' @param traveltime_seed Travel time seed parameter or "random"
#' @param fixed_params Fixed model parameters
#' 
#' @return List containing lateral_inflow_nodes, modified OCN, and random state params
generate_lateral_inflow_nodes <- function(OCN,
                                          cellsize,
                                          thresh,
                                          landscape_configuration_matrix,
                                          precip_timing,
                                          precip_intensity,
                                          traveltime_seed,
                                          fixed_params) {
  
  # ==================================================================
  # SUBCATCHMENT PATCH INITIALIZATION
  # ==================================================================
  
  # Initialize Flow Direction (FD) grid patches for subcatchment assignment
  OCN$FD$patches <- matrix(0, OCN$FD$nNodes, 1)                      # Patch assignments
  OCN$FD$patches_concentration <- matrix(0, OCN$FD$nNodes, 1)        # Mean immobile concentrations
  
  # Get lookup dictionary mapping river network (RN) nodes to their Strahler stream orders
  stream_order_dict <- get_stream_order_dict(OCN)
  
  # Initialize container for nodes that will receive lateral inflow
  # These nodes represent locations where landscape drainage enters the river network
  lateral_inflow_nodes <- list()
  
  # Calculate grid dimensions for spatial parameter assignment
  grid_width_cells <- max(OCN$FD$X) / cellsize + 1    # Number of grid columns 
  grid_height_cells <- max(OCN$FD$Y) / cellsize + 1   # Number of grid rows
  
  # Identify headwater nodes; these represent the most upstream points in the river network
  upstream <- OCN$RN$upstream
  upstream <- upstream[sapply(upstream, length) == 1]  # Keep only self-draining nodes
  
  # Track Flow Direction (FD) nodes already assigned to subcatchments
  # This prevents double-counting landscape areas in multiple subcatchments
  accounted_for_FD_nodes <- list()
  
  # ==================================================================
  # STREAM ORDER PROCESSING LOOP
  # ==================================================================
  
  # Process river network by stream order (1st order → higher orders)
  for (stream_order in c(1:max(OCN$AG$streamOrder))) {
    
    # ==================================================================
    # STREAM ORDER NODE IDENTIFICATION
    # ==================================================================
    
    # Identify all nodes belonging to the current stream order
    stream_order_node_list <- get_stream_order_node_list(OCN, stream_order_dict, stream_order)
    
    # Identify upstream boundary nodes for the current stream order
    upstream_node_list <- get_upstream_node_list(OCN, stream_order_node_list, stream_order_dict, stream_order)
    
    # ==================================================================
    # SUBCATCHMENT GENERATION ALONG STREAM ORDER REACHES
    # ==================================================================
    
    stream_order_progress_bar <- progress_bar$new(
      format = glue("  Lateral inflow nodes stream order: {stream_order} [:bar] :percent eta: :eta"),
      total = length(upstream_node_list), clear = FALSE, width = 60
    )
    
    for (node in upstream_node_list) {
      stream_order_progress_bar$tick()
      
      # ==================================================================
      # DOWNSTREAM TRAVERSAL WITHIN STREAM ORDER REACH
      # ==================================================================
      
      # Initialize variables to traverse downstream within current stream order
      # This ensures we only process nodes within the same stream order reach
      next_down_strahler_order <- stream_order  # Track stream order consistency
      next_down <- node                         # Current position in downstream traversal
      
      # Move downstream through nodes of same stream order, creating subcatchments
      while (next_down_strahler_order == stream_order) {
        
        # ==================================================================
        # DRAINAGE AREA CALCULATION FOR SUBCATCHMENT
        # ==================================================================
        
        # Find corresponding Flow Direction (FD) grid cell for current river node
        FD_node <- OCN$RN$toFD[next_down]
        upstream_FD_nodes <- list()
        
        # Recursively collect all FD cells that drain to this river node
        # This represents the total landscape area contributing flow to this point
        upstream_FD_nodes <- get_all_upstream_FD_nodes(
          OCN = OCN,
          FD_node = FD_node,
          all_upstream_FD_nodes = upstream_FD_nodes
        )
        
        # Subtract FD cells already assigned to other upstream subcatchments
        # This prevents double-counting landscape areas in overlapping subcatchments
        upstream_nodes_minus_accounted_for <- upstream_FD_nodes[!upstream_FD_nodes %in% accounted_for_FD_nodes]
        
        # Calculate drainage area for potential new subcatchment
        # Area = number of unassigned FD cells × cell area
        area <- cellsize * cellsize * length(upstream_nodes_minus_accounted_for)
        
        # ==================================================================
        # AREA THRESHOLD CHECK AND SUBCATCHMENT CREATION
        # ==================================================================
        
        # Create subcatchment only if drainage area exceeds minimum threshold
        if (area > thresh) {
          
          # Register this river node as a lateral inflow point
          # Lateral inflow nodes will receive discharge/concentration from landscape drainage
          lateral_inflow_nodes <- c(lateral_inflow_nodes, next_down)
          
          # Mark all FD cells in this subcatchment as assigned to prevent double-counting
          accounted_for_FD_nodes <- c(accounted_for_FD_nodes, upstream_nodes_minus_accounted_for)
          
          # Get spatial coordinates of the outlet node for this subcatchment
          node_x_position <- OCN$RN$X[next_down]         # X coordinate [m]
          node_y_position <- OCN$RN$Y[next_down]         # Y coordinate [m]
          
          # Convert spatial coordinates to parameter grid indices
          grid_column_index <- (node_x_position / cellsize) + 1    # Grid column index
          grid_row_index <- (node_y_position / cellsize) + 1       # Grid row index
          
          # Extract subcatchment immobile concentration 
          upstream_area <- OCN$RN$A[next_down]                                    # Total drainage area [m2]
          mean_c_im_number <- landscape_configuration_matrix[grid_column_index, grid_row_index]        # Mean concentration parameter
          
          # Assign subcatchment parameters to all FD cells draining to this river node
          # This links landscape-scale parameters to the appropriate outlet location
          for (FD_node_acc in upstream_nodes_minus_accounted_for) {
            OCN$FD$patches[FD_node_acc] <- next_down                              # Assign outlet node ID
            OCN$FD$patches_concentration[FD_node_acc] <- mean_c_im_number          # Mean concentration
          }
        }
        
        # ==================================================================
        # DOWNSTREAM ADVANCEMENT CONTROL
        # ==================================================================
        
        # Advance to the next downstream node or handle network outlet
        if (OCN$RN$downNode[next_down] != 0) {
          # Move to next downstream node within the river network
          next_down <- OCN$RN$downNode[next_down]
          next_down_strahler_order <- stream_order_dict[[next_down]]  # Update stream order for loop control
        } else {
          # ==================================================================
          # SPECIAL HANDLING FOR NETWORK OUTLET NODE
          # ==================================================================
          
          # Reached network outlet - create final subcatchment and exit while loop
          # The outlet always receives lateral inflow regardless of area threshold
          lateral_inflow_nodes <- c(lateral_inflow_nodes, next_down)
          accounted_for_FD_nodes <- c(accounted_for_FD_nodes, upstream_nodes_minus_accounted_for)
          
          # ==================================================================
          # OUTLET SUBCATCHMENT PARAMETER EXTRACTION
          # ==================================================================
          
          # Extract spatial parameters for the outlet subcatchment
          node_x_position <- OCN$RN$X[next_down]         # Outlet X coordinate [m]
          node_y_position <- OCN$RN$Y[next_down]         # Outlet Y coordinate [m]
          
          # Convert to parameter grid indices
          grid_column_index <- (node_x_position / cellsize) + 1
          grid_row_index <- (node_y_position / cellsize) + 1
          
          # Extract landscape parameters for outlet subcatchment
          upstream_area <- OCN$RN$A[next_down]                                    # Total drainage area [m2]
          mean_c_im_number <- landscape_configuration_matrix[grid_column_index, grid_row_index]        # Mean concentration parameter
          
          # ==================================================================
          # OUTLET FD CELL PARAMETER ASSIGNMENT
          # ==================================================================
          
          for (FD_node_acc in upstream_nodes_minus_accounted_for) {
            OCN$FD$patches[FD_node_acc] <- next_down                              # Assign to outlet node
            OCN$FD$patches_concentration[FD_node_acc] <- mean_c_im_number          # Mean concentration
          }
          
          # ==================================================================
          # OUTLET NODE FD CELL ASSIGNMENT
          # ==================================================================
          
          # Assign parameters to the FD cell corresponding to the outlet node itself
          # This ensures the outlet location has proper parameter values
          FD_outlet <- OCN$RN$toFD[next_down]                                  # Get FD cell for outlet node
          OCN$FD$patches[FD_outlet] <- next_down                               # Assign to itself
          OCN$FD$patches_concentration[FD_outlet] <- mean_c_im_number          # Mean concentration
          
          # Exit the downstream traversal loop since we've reached the outlet
          break
        }
      }
    }
  }
  
  # ==================================================================
  # RETURN RESULTS
  # ==================================================================
  
  # Package all results for return to calling function
  # Returns subcatchment configuration and updated OCN with patch assignments
  return(list(
    lateral_inflow_nodes,        # [1] Vector of river nodes receiving lateral inflow from subcatchments
    OCN                          # [2] Updated OCN object with FD patch assignments and parameters
  ))
}

#' Route discharge and concentration from upstream nodes to downstream node
#' 
#' Applies hydraulic routing with travel time and dispersion effects
#' 
#' @param OCN River network object
#' @param upstream_nodes Vector of source node IDs
#' @param downstream_node Target node ID
#' @param df_nodes_discharge Time series dataframe of discharge values
#' @param df_nodes_conc Time series dataframe of concentration values
#' @param cellsize Grid cell size in meters
#' @param vf Settling velocity parameter
#' @param mannings_n Manning's roughness coefficient
#' @param side_slope Channel side slope parameter
#' @param routing_function_ocn Python routing function
#' 
#' @return List with routed discharge/concentration dataframes and hydraulic properties
routing_routine <- function(OCN,
                            upstream_nodes,
                            downstream_node,
                            df_nodes_discharge,
                            df_nodes_conc,
                            cellsize,
                            vf,
                            mannings_n,
                            side_slope,
                            routing_function_ocn) {
  
  # ==================================================================
  # INPUT TIME SERIES EXTRACTION AND VALIDATION
  # ==================================================================
  
  # Extract time series for all upstream nodes that flow into the downstream target
  # These represent the discharge and concentration inputs that need to be routed
  df_discharge_upstream_routed <- df_nodes_discharge[, upstream_nodes]  # Discharge time series [m3/s]
  df_conc_upstream_routed <- df_nodes_conc[, upstream_nodes]            # Concentration time series [mass/volume]
  
  # ==================================================================
  # DATA STRUCTURE VALIDATION AND CONVERSION
  # ==================================================================
  
  # Ensure input data is in proper dataframe format for consistent processing
  # Single column inputs may be converted to vectors, so force dataframe structure
  if (!is.data.frame(df_discharge_upstream_routed)) {
    df_discharge_upstream_routed <- as.data.frame(df_discharge_upstream_routed)
  }
  
  if (!is.data.frame(df_conc_upstream_routed)) {
    df_conc_upstream_routed <- as.data.frame(df_conc_upstream_routed)
  }
  
  # ==================================================================
  # INDIVIDUAL TRIBUTARY ROUTING LOOP
  # ==================================================================
  
  # Process each upstream tributary separately through hydraulic routing
  # before they are mixed at the downstream confluence
  for (node_index in 1:length(upstream_nodes)) {
    
    # ==================================================================
    # TRIBUTARY-SPECIFIC INPUT PREPARATION
    # ==================================================================
    
    # Extract current tributary information for routing calculations
    upstream_node <- upstream_nodes[node_index]                    # Current tributary node ID
    inflow <- df_discharge_upstream_routed[, node_index]          # Tributary discharge time series [m3/s]
    conc <- df_conc_upstream_routed[, node_index]                 # Tributary concentration time series [mass/volume]
    
    # Calculate routing distance from tributary to confluence point
    reach_length <- OCN$RN$downstreamPathLength[upstream_node, downstream_node]  # Network path length [m]
    routing_distance <- reach_length                                             # Distance for routing calculations [m]
    
    # ==================================================================
    # HYDRAULIC ROUTING CALCULATION VIA PYTHON FUNCTION
    # ==================================================================
    
    # Call Python-based hydraulic routing function to simulate water and solute transport
    routing_results <- routing_function_ocn(
      inflow = inflow,                                # Input discharge time series [m3/s]
      conc = conc,                                   # Input concentration time series [mass/volume]
      dx = routing_distance,                         # Routing distance [m]
      dt_ref = 24,                                   # Reference time step [hours] - daily time step
      bottom_slope = OCN$RN$slope[upstream_node],    # Channel bottom slope [dimensionless]
      mannings_n = mannings_n,                       # Manning's roughness coefficient [s/m^(1/3)]
      reach_length = reach_length,                   # Physical reach length [m]
      bottom_width_int = OCN$RN$width[upstream_node], # Channel bottom width [m]
      side_slope = side_slope,                       # Channel side slope [horizontal:vertical]
      vf = vf                                        # Settling velocity parameter [m/s]
    )
    
    # ==================================================================
    # ROUTING RESULTS EXTRACTION AND PROCESSING
    # ==================================================================
    
    # Extract routing outputs from Python function results
    # The Python function returns multiple components of the hydraulic solution
    routed_discharge <- routing_results[[1]]       # Time-shifted and attenuated discharge [m3/s]
    routed_concentration <- routing_results[[2]]   # Time-shifted and attenuated concentration [mg/l]
    
    # Extract hydraulic properties for channel characterization and analysis
    # These represent time-averaged hydraulic conditions in the reach
    median_water_depth <- routing_results[[3]]     # Median hydraulic depth over simulation [m]
    median_flow_celerity <- routing_results[[4]]   # Median flow celerity (wave speed) [m/s]
    
    # ==================================================================
    # ROUTED TIME SERIES STORAGE
    # ==================================================================
    
    # Store routed time series back into tributary-specific columns
    # Convert Python list outputs to R vectors for consistent data handling
    df_discharge_upstream_routed[, node_index] <- unlist(routed_discharge)    # Routed discharge time series
    df_conc_upstream_routed[, node_index] <- unlist(routed_concentration)     # Routed concentration time series
  }
  
  # ==================================================================
  # RETURN RESULTS
  # ==================================================================
  
  # Package routing results for return to calling function
  # Returns both routed time series and hydraulic characterization
  return(list(
    df_discharge_upstream_routed,   # [1] Routed discharge time series for all tributaries [m3/s]
    df_conc_upstream_routed,        # [2] Routed concentration time series for all tributaries [mass/volume]
    median_water_depth,             # [3] Median hydraulic depth for channel characterization [m]
    median_flow_celerity            # [4] Median flow celerity for velocity analysis [m/s]
  ))
}

#' Mix flows and concentrations at confluence points
#' 
#' @param downstream_node Target node ID
#' @param df_discharge_upstream_routed Dataframe of routed discharge inputs
#' @param df_conc_upstream_routed Dataframe of routed concentration inputs
#' @param df_nodes_discharge Output discharge time series dataframe
#' @param df_nodes_conc Output concentration time series dataframe
#' 
#' @return List with updated discharge and concentration dataframes
mixing_routine <- function(downstream_node,
                           df_discharge_upstream_routed,
                           df_conc_upstream_routed,
                           df_nodes_discharge,
                           df_nodes_conc) {
  
  # Sum flow of all upstream nodes
  df_nodes_discharge[, downstream_node] <- base::rowSums(df_discharge_upstream_routed)
  
  # Calculate solute load (mass flux) for each tributary: Load = Discharge × Concentration
  df_nodes_upstream_loads <- df_discharge_upstream_routed * df_conc_upstream_routed
  
  # Sum load of all upstream nodes
  df_nodes_load <- base::rowSums(df_nodes_upstream_loads)
  
  # Calcualte concetration based on load and discharge
  df_nodes_conc[, downstream_node] <- df_nodes_load / df_nodes_discharge[, downstream_node]
  
  # ==================================================================
  # RETURN UPDATED TIME SERIES
  # ==================================================================
  
  # Return updated network-wide time series with mixed flows at confluence
  # Both discharge and concentration matrices now include the downstream mixed values
  return(list(
    df_nodes_discharge,   # [1] Updated discharge time series matrix [m3/s]
    df_nodes_conc        # [2] Updated concentration time series matrix [mass/volume]
  ))
}

#' Generate lateral inflow time series for a subcatchment
#' 
#' Uses subcatchment-scale model to simulate discharge and concentration inputs
#' 
#' @param OCN River network object
#' @param df_nodes_discharge Discharge time series dataframe
#' @param df_nodes_conc Concentration time series dataframe
#' @param node Target node ID
#' @param cellsize Grid cell size in meters
#' @param tmax_yrs Maximum simulation time in years
#' @param warm_up Warm-up period length
#' @param fixed_params Fixed model parameters
#' @param volume_in Volume mass balance input
#' @param mass_in Mass balance input
#' @param subcatchment_scale_module Python subcatchment modeling function
#' @param precip_timing Precipitation timing parameter or "random"
#' @param precip_intensity Precipitation intensity parameter or "random"
#' @param traveltime_seed Travel time seed parameter or "random"
#' 
#' @return List with updated discharge/concentration and mass balance terms
lateral_inflow_routine <- function(OCN,
                                   df_nodes_discharge,
                                   df_nodes_conc,
                                   node,
                                   cellsize,
                                   tmax_yrs,
                                   warm_up,
                                   fixed_params,
                                   volume_in,
                                   mass_in,
                                   subcatchment_scale_module,
                                   precip_timing,
                                   precip_intensity,
                                   traveltime_seed) {
  
  # ==================================================================
  # SUBCATCHMENT IDENTIFICATION AND AREA CALCULATION
  # ==================================================================
  
  # Get the corresponding Flow Direction (FD) node for the river network node
  # This links the river node to the landscape grid for subcatchment delineation
  FD_node <- OCN$RN$toFD[node]
  
  # Identify the subcatchment (patch) ID assigned to this FD node
  # All FD cells with the same patch ID drain to this river node
  patch_id <- OCN$FD$patches[FD_node]
  
  # Find all FD cells that belong to this subcatchment
  # This represents the complete landscape area draining to the target river node
  all_upstream_FD_nodes <- which(OCN$FD$patches == patch_id)
  
  # ==================================================================
  # TIME SERIES INITIALIZATION
  # ==================================================================
  
  # Initialize arrays to store lateral inflow contributions
  # These will accumulate discharge and load inputs from subcatchment drainage
  lateral_inflow <- rep(0, length(df_nodes_discharge[, node]))    # Lateral discharge time series [m3/s]
  lateral_load <- rep(0, length(df_nodes_discharge[, node]))      # Lateral load time series [mass/time]
  
  # Calculate total subcatchment drainage area
  # Area = number of FD cells × area per cell
  patch_area <- length(all_upstream_FD_nodes) * cellsize * cellsize  # Total drainage area [m2]
  
  # ==================================================================
  # SUBCATCHMENT PARAMETER EXTRACTION
  # ==================================================================
  
  # Extract biogeochemical and hydrological parameters for this subcatchment
  # These parameters were assigned during subcatchment generation based on spatial location
  mean_c_im_number <- OCN$FD$patches_concentration[FD_node]        # Mean concentration parameter [mass/volume]
  elevation <- OCN$FD$Z[FD_node]                                   # Elevation at subcatchment outlet [m]
  
  # ==================================================================
  # ADDITIONAL STOCHASTIC SEED GENERATION
  # ==================================================================
  
  # Generate additional random seed for asynchrony parameter
  # This provides independent stochastic variability for landscape heterogeneity
  rho_seed <- sample(2^16, 1)  # Random seed for asynchrony processes
  
  if (precip_timing == "random") {
    random_state_times <- sample(2^16, 1)
  } else {
    random_state_times <- precip_timing
  }
  
  if (precip_intensity == "random") {
    random_state_rain <- sample(2^16, 1)
  } else {
    random_state_rain <- precip_intensity
  }
  
  if (traveltime_seed == "random") {
    real_random <- sample(2^16, 1)
  } else {
    real_random <- traveltime_seed
  }
  
  # ==================================================================
  # PARAMETER COMPILATION FOR SUBCATCHMENT MODEL
  # ==================================================================
  
  # Compile variable parameters specific to this subcatchment
  # These control biogeochemical processes, stochastic timing, and scaling behavior
  params <- list(
    mean_c_im_number = mean_c_im_number,   # Mean concentration parameter
    random_state_times = random_state_times, # Precipitation timing seed
    random_state_rain = random_state_rain,   # Precipitation intensity seed
    real_random = real_random,               # Travel time seed
    elevation = elevation,                   # Subcatchment elevation
    rho_seed = rho_seed      # Asynchrony random seed
  )
  
  # Combine subcatchment-specific parameters with fixed hydrological parameters
  # Fixed parameters include soil properties, climate drivers, and model constants
  params <- c(params, fixed_params)
  
  # ==================================================================
  # SUBCATCHMENT-SCALE MODEL EXECUTION VIA PYTHON
  # ==================================================================
  
  # Call Python-based subcatchment model to generate discharge and concentration time series
  # This model simulates landscape-scale hydrological and biogeochemical processes including:
  # - Precipitation-driven runoff generation
  # - Soil moisture dynamics and evapotranspiration
  # - Groundwater-surface water interactions
  # - Biogeochemical transformations and solute mobilization
  result <- subcatchment_scale_module(param_dict = params)
  
  # ==================================================================
  # DISCHARGE TIME SERIES PROCESSING AND UNIT CONVERSION
  # ==================================================================
  
  # Extract and process discharge time series with proper unit conversion
  discharge <- result[[1]]                    # Raw discharge from Python model [cm/day]
  discharge <- discharge / 100                # Convert cm/day to m/day
  discharge <- discharge / (24 * 60 * 60)    # Convert m/day to m/s per unit area [m/s per m2]
  discharge <- discharge * patch_area         # Scale by subcatchment area to get total discharge [m3/s]
  
  # ==================================================================
  # VOLUME MASS BALANCE ACCUMULATION
  # ==================================================================
  
  # Accumulate volume input for system-wide mass balance (excluding warm-up period)
  # This tracks total water input from all subcatchments for model validation
  volume_in <- volume_in + sum(discharge[warm_up:length(discharge)] * 3600 * 24)  # Convert m3/s to m3/day and sum
  
  # ==================================================================
  # CONCENTRATION TIME SERIES EXTRACTION
  # ==================================================================
  
  # Extract concentration time series from subcatchment model results
  # This represents the biogeochemical signature of water draining from the landscape
  concentration <- result[[2]]  # Concentration time series [mass/volume]
  
  # ==================================================================
  # LATERAL INFLOW TIME SERIES CONSTRUCTION
  # ==================================================================
  
  # Add subcatchment discharge to lateral inflow time series
  # This accumulates contributions from all subcatchments draining to this river node
  lateral_inflow <- unlist(lapply(1:length(lateral_inflow), function(i) lateral_inflow[i] + discharge[i]))

  # Calculate lateral solute load time series: Load = Discharge × Concentration
  # This converts concentration to mass transport rate for proper mass balance
  lateral_load <- unlist(lapply(1:length(lateral_load), function(i) lateral_load[i] + (discharge[i] * concentration[i])))
  
  # ==================================================================
  # MASS BALANCE ACCUMULATION
  # ==================================================================
  
  # Accumulate solute mass input for system-wide mass balance (excluding warm-up period)
  # This tracks total solute mass input from all subcatchments for model validation
  mass_in <- mass_in + sum(lateral_load[warm_up:length(lateral_load)] * 3600 * 24)  # Convert mass/s to mass/day and sum
  
  # Calculate flow-weighted average concentration for lateral inflow
  # This represents the mixed concentration from all subcatchment contributions
  lateral_concentration <- lateral_load / lateral_inflow
  
  # ==================================================================
  # RIVER NETWORK MIXING WITH LATERAL INFLOW
  # ==================================================================
  
  # Calculate current solute load in the river node before lateral inflow addition
  # Load = Discharge × Concentration (converts concentration to mass transport rate)
  df_nodes_load <- df_nodes_discharge[, node] * df_nodes_conc[, node]
  
  # ==================================================================
  # DISCHARGE MIXING AT RIVER NODE
  # ==================================================================
  
  # Add lateral inflow to existing river discharge using mass conservation
  # New discharge = River discharge + Lateral inflow discharge
  df_nodes_discharge[, node] <- df_nodes_discharge[, node] + lateral_inflow
  
  # ==================================================================
  # SOLUTE LOAD MIXING AT RIVER NODE
  # ==================================================================
  
  # Add lateral solute load to existing river load using mass conservation
  # New load = River load + Lateral load
  df_nodes_load <- df_nodes_load + lateral_load
  
  # ==================================================================
  # MIXED CONCENTRATION CALCULATION
  # ==================================================================
  
  # Calculate mixed concentration after lateral inflow addition
  # Mixed concentration = Total load / Total discharge
  # This implements mass conservation for solutes at the river node
  df_nodes_conc[, node] <- df_nodes_load / df_nodes_discharge[, node]
  
  # ==================================================================
  # RETURN UPDATED RESULTS
  # ==================================================================
  
  # Return updated time series and mass balance components
  # All results now include the lateral contributions from subcatchment drainage
  return(list(
    df_nodes_discharge,   # [1] Updated discharge time series with lateral inflow [m3/s]
    df_nodes_conc,        # [2] Updated concentration time series with lateral mixing [mass/volume]
    volume_in,            # [3] Cumulative volume input for mass balance [m3]
    mass_in               # [4] Cumulative mass input for mass balance [mass units]
  ))
}

#' Recursively populate downstream nodes with routed discharge and concentration
#' 
#' Processes river network downstream, applying routing and mixing at each step
#' 
#' @param OCN River network object
#' @param node Current node ID
#' @param df_nodes_discharge Discharge time series dataframe
#' @param df_nodes_conc Concentration time series dataframe
#' @param cellsize Grid cell size in meters
#' @param landscape_configuration_matrix Matrix of landscape heterogeneity (mean immobile concentration of subcatchments)
#' @param lateral_inflow_nodes Vector of lateral inflow node IDs
#' @param tmax_yrs Maximum simulation time in years
#' @param warm_up Warm-up period length
#' @param vf Settling velocity parameter
#' @param mannings_n Manning's roughness coefficient
#' @param side_slope Channel side slope parameter
#' @param fixed_params Fixed model parameters
#' @param volume_in Volume mass balance input
#' @param mass_in Mass balance input
#' @param precip_timing Precipitation timing parameter or "random"
#' @param precip_intensity Precipitation intensity parameter or "random"
#' @param traveltime_seed Travel time seed parameter or "random"
#' @param subcatchment_scale_module Python subcatchment modeling function
#' @param routing_function_ocn Python routing function
#' @param depth_dict Dictionary for storing hydraulic depth properties
#' @param celerity_dict Dictionary for storing hydraulic celerity properties
#' 
#' @return List with updated time series, mass balance, and hydraulic dictionaries
populate_downstream <- function(OCN,
                                node,
                                df_nodes_discharge,
                                df_nodes_conc,
                                cellsize,
                                landscape_configuration_matrix,
                                lateral_inflow_nodes,
                                tmax_yrs,
                                warm_up,
                                vf,
                                mannings_n,
                                side_slope,
                                fixed_params,
                                volume_in,
                                mass_in,
                                precip_timing,
                                precip_intensity,
                                traveltime_seed,
                                subcatchment_scale_module,
                                routing_function_ocn,
                                depth_dict,
                                celerity_dict) {
  # ==================================================================
  # INITIALIZATION AND STREAM ORDER LOOKUP
  # ==================================================================
  
  # Get lookup dictionary mapping river network (RN) nodes to their Strahler stream orders
  stream_order_dict <- get_stream_order_dict(OCN)
  
  # Identify the immediate downstream node from current position
  downstream_node <- OCN$RN$downNode[node]
  
  # ==================================================================
  # RECURSIVE TERMINATION CONDITIONS
  # ==================================================================
  
  # Check if we've reached the network outlet (no further downstream node)
  if (downstream_node != 0) {
    
    # Check stream order compatibility for continued processing
    # Only process nodes within the same stream order reach
    stream_order_node <- stream_order_dict[[node]]           # Current node's stream order
    stream_order_downstream <- stream_order_dict[[downstream_node]]  # Downstream node's stream order
    
    # Terminate recursion if downstream node belongs to different stream order
    # This prevents crossing stream order boundaries during within-reach processing
    if (stream_order_node != stream_order_downstream) {
      result_list <- list(df_nodes_discharge, df_nodes_conc, volume_in, mass_in, depth_dict, celerity_dict)
      return(result_list)
    }
    
    # ==================================================================
    # UPSTREAM NODE IDENTIFICATION AND FILTERING
    # ==================================================================
    
    # Identify all nodes that flow directly into the downstream target node
    # This includes the current node and any other tributaries joining at this point
    upstream_nodes <- which(OCN$RN$downNode == downstream_node)
    
    # Filter upstream nodes to include only those of equal or lower stream order
    # This ensures we don't include higher-order tributaries that shouldn't be processed yet
    upstream_nodes <- upstream_nodes[stream_order_dict[upstream_nodes] <= stream_order_node]
    
    # ==================================================================
    # HYDRAULIC ROUTING TO DOWNSTREAM NODE
    # ==================================================================
    
    # Apply hydraulic routing to transport discharge and concentration
    # from all upstream nodes to the downstream node with travel time and dispersion
    result_routing <- routing_routine(
      OCN = OCN,
      upstream_nodes = upstream_nodes,
      downstream_node = downstream_node,
      df_nodes_discharge = df_nodes_discharge,
      df_nodes_conc = df_nodes_conc,
      cellsize = cellsize,
      vf = vf,
      mannings_n = mannings_n,
      side_slope = side_slope,
      routing_function_ocn = routing_function_ocn
    )
    
    # Extract hydraulic routing results
    df_discharge_upstream_routed <- result_routing[[1]]  # Routed discharge time series from all tributaries
    df_conc_upstream_routed <- result_routing[[2]]       # Routed concentration time series from all tributaries
    median_depth <- result_routing[[3]]                  # Median hydraulic depth for channel characterization
    median_celerity <- result_routing[[4]]               # Median flow celerity for velocity analysis
    
    # Store hydraulic properties for post-simulation analysis
    # These values characterize the hydraulic conditions at each node
    depth_dict[downstream_node] <- median_depth
    celerity_dict[downstream_node] <- median_celerity
    
    # ==================================================================
    # FLOW MIXING AT DOWNSTREAM NODE
    # ==================================================================
    
    # Combine all routed upstream flows using mass conservation principles
    # This calculates the mixed discharge and concentration at the confluence
    result_mixing <- mixing_routine(
      downstream_node = downstream_node,
      df_discharge_upstream_routed = df_discharge_upstream_routed,
      df_conc_upstream_routed = df_conc_upstream_routed,
      df_nodes_discharge = df_nodes_discharge,
      df_nodes_conc = df_nodes_conc
    )
    
    # Update network-wide time series with mixed flow results
    df_nodes_discharge <- result_mixing[[1]]  # Updated discharge matrix
    df_nodes_conc <- result_mixing[[2]]       # Updated concentration matrix
    
    # ==================================================================
    # LATERAL INFLOW ADDITION (if applicable)
    # ==================================================================
    
    # Check if downstream node receives lateral inflow from surrounding landscape
    if (downstream_node %in% lateral_inflow_nodes) {
      
      # Generate additional discharge and concentration from subcatchment drainage
      # This represents water and solutes entering from the landscape surrounding this node
      result_lateral <- lateral_inflow_routine(
        OCN = OCN,
        df_nodes_discharge = df_nodes_discharge,
        df_nodes_conc = df_nodes_conc,
        node = downstream_node,
        cellsize = cellsize,
        tmax_yrs = tmax_yrs,
        warm_up = warm_up,
        fixed_params = fixed_params,
        volume_in = volume_in,
        mass_in = mass_in,
        subcatchment_scale_module = subcatchment_scale_module,
        precip_timing,
        precip_intensity,
        traveltime_seed
      )

      # Update time series and mass balance with lateral contributions
      df_nodes_discharge <- result_lateral[[1]]  # Discharge including lateral inflow
      df_nodes_conc <- result_lateral[[2]]       # Mixed concentration with lateral inputs
      volume_in <- result_lateral[[3]]           # Updated volume mass balance
      mass_in <- result_lateral[[4]]             # Updated mass balance
    }
    
    # ==================================================================
    # RECURSIVE DOWNSTREAM PROPAGATION
    # ==================================================================
    
    # Continue processing further downstream if not at network outlet
    if (OCN$RN$downNode[downstream_node] != 0) {
      
      # Recursively call this function to continue processing downstream
      # This enables depth-first traversal of the river network, processing
      # each reach completely before moving to the next confluence
      result_pop_down <- populate_downstream(
        OCN = OCN,
        node = downstream_node,                    # Start from current downstream node
        df_nodes_discharge = df_nodes_discharge,
        df_nodes_conc = df_nodes_conc,
        cellsize = cellsize,
        landscape_configuration_matrix = landscape_configuration_matrix,
        lateral_inflow_nodes = lateral_inflow_nodes,
        tmax_yrs = tmax_yrs,
        warm_up = warm_up,
        vf = vf,
        mannings_n = mannings_n,
        side_slope = side_slope,
        fixed_params = fixed_params,
        volume_in = volume_in,                     # Pass updated mass balance
        mass_in = mass_in,
        precip_timing = precip_timing,
        precip_intensity = precip_intensity,
        traveltime_seed = traveltime_seed,
        subcatchment_scale_module = subcatchment_scale_module,
        routing_function_ocn = routing_function_ocn,
        depth_dict = depth_dict,                   # Pass hydraulic property storage
        celerity_dict = celerity_dict
      )
      
      # Extract results from recursive call
      df_nodes_discharge <- result_pop_down[[1]]   # Updated discharge time series
      df_nodes_conc <- result_pop_down[[2]]        # Updated concentration time series
      volume_in <- result_pop_down[[3]]            # Cumulative volume input
      mass_in <- result_pop_down[[4]]              # Cumulative mass input
      depth_dict <- result_pop_down[[5]]           # Hydraulic depth properties
      celerity_dict <- result_pop_down[[6]]        # Hydraulic celerity properties
    }
  }
  
  # ==================================================================
  # RETURN RESULTS
  # ==================================================================
  
  # Package all updated results for return to calling function
  # This includes both simulation outputs and hydraulic characterization
  result_list <- list(
    df_nodes_discharge,   # [1] Updated discharge time series matrix
    df_nodes_conc,        # [2] Updated concentration time series matrix  
    volume_in,            # [3] Cumulative volume input for mass balance
    mass_in,              # [4] Cumulative mass input for mass balance
    depth_dict,           # [5] Hydraulic depth properties for all processed nodes
    celerity_dict         # [6] Hydraulic celerity properties for all processed nodes
  )
  
  return(result_list)
}

#' Create matrix of the subcatchment immobile concentration based on gamma_ls and sigma_w_ls
#' 
#' @param df_downstream_path Network topology dataframe
#' @param gamma_ls Network-wide scaling parameter
#' @param sigma_w_ls Network-wide variability parameter
#' @param overall_mean_concentration Overall mean concentration value
#' 
#' @return Matrix of immobile concentration
generate_landscape_configuration <- function(
    df_downstream_path,
    gamma_ls,
    sigma_w_ls,
    overall_mean_concentration) {
  # ===============================================
  # GRID DIMENSIONS AND INITIALIZATION
  # ===============================================
  # Extract spatial grid dimensions from network topology matrix
  grid_width_cells <- nrow(df_downstream_path)
  grid_height_cells <- ncol(df_downstream_path)
  
  # ===============================================
  # PARAMETER MATRIX INITIALIZATION
  # ===============================================
  # Create base parameter matrix template with grid dimensions
  landscape_configuration_matrix <- matrix(0, nrow = grid_width_cells, ncol = grid_height_cells)
  
  # ===============================================
  # SPATIALLY CORRELATED CONCENTRATION GENERATION
  # ===============================================
  # Initialize list to track concentration values for mean calculation
  mean_conc_list <- c()
  
  # Generate spatially correlated concentration based on dwonstream flow distance
  for (row_index in 1:grid_width_cells) {
    for (col_index in 1:grid_height_cells) {
      # Only generate parameters for cells that are part of the drainage network
      if (df_downstream_path[row_index, col_index] != 0) {
        # Calculate concentration
        mean_c <- calc_c_immobile(df_downstream_path[row_index, col_index], gamma_ls, sigma_w_ls)
        landscape_configuration_matrix[row_index, col_index] <- mean_c
        mean_conc_list <- c(mean_conc_list, landscape_configuration_matrix[row_index, col_index])
      }
    }
  }
  
  # Calculate scaling factor to achieve target mean concentration
  calc_mean_conc <- mean(mean_conc_list)
  scaling_factor <- overall_mean_concentration / calc_mean_conc
  
  # Apply scaling factor to all concentration parameters
  for (row_index in 1:grid_width_cells) {
    for (col_index in 1:grid_height_cells) {
      if (df_downstream_path[row_index, col_index] != 0) {
        # Scale concentration to match target overall mean
        landscape_configuration_matrix[row_index, col_index] <- landscape_configuration_matrix[row_index, col_index] * scaling_factor
      }
    }
  }
  
  # ===============================================
  # RETURN PARAMETER MATRIX COLLECTION
  # ===============================================
  return(landscape_configuration_matrix)
}
# ===== END OF generate_landscape_configuration FUNCTION =====


#' Calculate immobile zone concentration based on distance scaling
#' 
#' Uses power law with log-normal variability to generate concentration values
#' 
#' @param distance Downstream path length
#' @param gamma_ls Scaling exponent parameter
#' @param sigma_w_ls Log-normal standard deviation parameter
#' 
#' @return Numeric concentration value
calc_c_immobile <- function(distance, gamma_ls, sigma_w_ls) {
  
  # Generate random variability from log-normal distribution
  random_variable <- rlnorm(1, meanlog = 0, sdlog = sigma_w_ls)
  
  # Physical interpretation:
  # - distance: downstream path length from outlet [L]
  # - gamma: scaling exponent controlling distance dependence
  #   * gamma > 0: concentrations increase with distance from outlet
  #   * gamma < 0: concentrations decrease with distance from outlet
  #   * gamma = 0: no distance dependence (spatially uniform mean)
  # - random_variable: log-normal multiplicative noise factor
  
  c_im <- (distance^gamma_ls) * random_variable
  
  # Return concentration value with both systematic (power law)
  # and stochastic (log-normal) components
  return(c_im)
}
# ===== END OF calc_c_immobile FUNCTION =====


#' Main function to run the stochastic landscape model
#' 
#' Orchestrates full simulation: loads network, generates subcatchments, runs routing
#' 
#' @param network OCN file path identifier
#' @param gamma_sc Subcatchment biogeochemical scaling parameter
#' @param sigma_w_sc Subcatchment biogeochemical variability parameter
#' @param gamma_ls Network-scale scaling parameter
#' @param sigma_w_ls Network-scale variability parameter
#' @param interarrival_time_mean Mean interarrival time
#' @param rho Asynchrony parameter for landscape heterogeneity
#' @param vf Settling velocity
#' @param damkohler_transport Transport Damkohler number
#' @param tmax_yrs Maximum simulation time in years
#' @param warm_up Warm-up period length
#' @param overall_mean_concentration Overall mean concentration value
#' @param oro_scaling Orographic scaling parameter
#' @param damkohler_interarrival Interarrival Damkohler number
#' @param damkohler_longterm Long-term Damkohler number
#' @param rain_per_year Annual rainfall amount
#' @param ET_max Maximum evapotranspiration rate
#' @param z_r Root zone depth
#' @param z_vz Vadose zone depth
#' @param theta_fc Field capacity soil moisture
#' @param theta_wp Wilting point soil moisture
#' @param theta_res Residual soil moisture
#' @param R Recharge rate
#' @param mean_tr Mean residence time
#' @param theta_sat Saturated soil moisture
#' @param aquifer_z Aquifer depth
#' @param mannings_n Manning's roughness coefficient
#' @param side_slope Channel side slope
#' @param simulation_identifier Unique simulation identifier
#' @param print_progress Logical flag for progress output
#' @param log_file Log file path
#' @param path_to_python Python executable path
#' @param path_to_python_functions Python functions file path
#' @param path_to_OCN OCN data file path
#' @param path_to_output Output directory path
#' @param log_warning_parallel Parallel processing warning log
#' 
#' @return Void (saves results to files)
run_SLM <- function(network,
                    gamma_sc,
                    sigma_w_sc,
                    gamma_ls,
                    sigma_w_ls,
                    interarrival_time_mean,
                    rho,
                    vf,
                    damkohler_transport,
                    tmax_yrs,
                    warm_up,
                    overall_mean_concentration,
                    oro_scaling,
                    damkohler_interarrival,
                    damkohler_longterm,
                    rain_per_year,
                    ET_max,
                    z_r,
                    z_vz,
                    theta_fc,
                    theta_wp,
                    theta_res,
                    R,
                    mean_tr,
                    theta_sat,
                    aquifer_z,
                    mannings_n,
                    side_slope,
                    simulation_identifier,
                    print_progress,
                    log_file,
                    path_to_python,
                    path_to_python_functions,
                    path_to_OCN,
                    path_to_output,
                    log_warning_parallel) {
  
  # ==================================================================
  # INITIALIZE SIMULATION  
  # ==================================================================
  
  # Initialize mass and volume balance tracking variables
  volume_in <- 0    # Total water volume entering the system [m3]
  volume_out <- 0   # Total water volume leaving the system [m3]
  mass_in <- 0      # Total solute mass entering the system [mass units]
  mass_out <- 0     # Total solute mass leaving the system [mass units]
  
  # Configure progress bar display based on user preference
  options(progress_enabled = print_progress)
  
  # Setup Python environment
  use_python(path_to_python)
  source_python(path_to_python_functions)
  
  # Load the Optimal Channel Network (OCN) object from file
  load(glue(path_to_OCN))
  
  # Extract cell size
  cellsize <- OCN$cellsize
  
  # Extract subcatchment area threshold from network identifier string
  # Network name format: "X_Y_threshold_identifier"
  subcatchment_area_treshold <- as.numeric(strsplit(network, "_")[[1]][3])
  
  # Convert threshold from km2 to m2 for area comparison
  thresh <- subcatchment_area_treshold * cellsize^2
  
  # Get lookup dictionary mapping river network (RN) nodes to their Strahler stream orders
  stream_order_dict <- get_stream_order_dict(OCN)
  
  # ==================================================================
  # SIMULATION
  # ==================================================================
  
  # Create spatial grid matrix to store downstream path lengths to outlet
  # This is used for generating spatially correlated parameter fields
  df_downstream_path <- get_downstream_matrix(OCN)
  
  # Generate matrix with immobile concentrations based on parameters gamma_ls and sigma_w_ls
  landscape_configuration_matrix <- generate_landscape_configuration(
    df_downstream_path = df_downstream_path,
    gamma_ls = gamma_ls,
    sigma_w_ls = sigma_w_ls,
    overall_mean_concentration = overall_mean_concentration
  )
  
  # Initialize storage for hydraulic properties computed during routing
  depth_dict <- matrix(0, OCN$RN$nNodes, 1)     # Median water depths [m]
  celerity_dict <- matrix(0, OCN$RN$nNodes, 1)  # Median flow celerities [m/s]
  
  # Configure random seed usage based on hydroclimatic scenario:
  precip_timing <- sample(2^16, 1)     # Fixed seed for precipitation timing
  precip_intensity <- sample(2^16, 1)  # Fixed seed for precipitation intensity
  traveltime_seed <- "random"
  
  # SUBCATCHMENT GENERATION AND PARAMETER ASSIGNMENT
  
  # Generate lateral inflow nodes based on area threshold
  # This creates subcatchments of approximately equal size to headwater catchments
  # and assigns spatially distributed parameters (immobile concentration of subcatchment)
  # to each subcatchment based on the landscape_configuration_matrix
  result <- generate_lateral_inflow_nodes(
    OCN = OCN,
    cellsize = cellsize,
    thresh = thresh,
    landscape_configuration_matrix = landscape_configuration_matrix,
    precip_timing = precip_timing,
    precip_intensity = precip_intensity,
    traveltime_seed = traveltime_seed,
    fixed_params = fixed_params
  )
  lateral_inflow_nodes <- result[[1]]  # List of nodes receiving lateral inflow
  OCN <- result[[2]] # Updated OCN with immobile concentration assigned to subcatchments
  
  # ==================================================================
  # HEADWATER NODE IDENTIFICATION AND TIME SERIES INITIALIZATION
  # ==================================================================
  
  # Identify headwater (source) nodes: nodes with no upstream connections
  # These nodes represent the most upstream points in the network
  upstream <- OCN$RN$upstream
  upstream <- upstream[sapply(upstream, length) == 1]  # Select nodes with only themselves as upstream
  
  # Remove headwater nodes from lateral inflow nodes to avoid double-counting
  # Headwaters get direct subcatchment input, not lateral inflow
  lateral_inflow_nodes <- setdiff(lateral_inflow_nodes, upstream)
  
  # Initialize time series matrices for concentration and discharge
  # Rows = daily time steps, Columns = network nodes
  df_nodes_conc <- data.frame(matrix(NA, nrow = (tmax_yrs * 365), ncol = OCN$RN$nNodes[1]))
  df_nodes_discharge <- data.frame(matrix(NA, nrow = (tmax_yrs * 365), ncol = OCN$RN$nNodes[1]))
  
  # Set column names to correspond to node IDs
  colnames(df_nodes_conc) <- 1:OCN$RN$nNodes[1]
  colnames(df_nodes_discharge) <- 1:OCN$RN$nNodes[1]
  
  # ==================================================================
  # HEADWATER TIME SERIES GENERATION
  # ==================================================================
  
  # Generate concentration and discharge time series for all headwater nodes
  pb <- progress_bar$new(
    format = "  Creating headwater time series [:bar] :percent eta: :eta",
    total = length(upstream), clear = FALSE, width = 60
  )
  
  for (upstream_index in 1:length(upstream)) {
    pb$tick()
    
    # Get spatial coordinates and grid indices for current headwater node
    upstream_node <- upstream[[upstream_index]]
    node_x_position <- OCN$RN$X[upstream_node]
    node_y_position <- OCN$RN$Y[upstream_node]
    
    # Convert spatial coordinates to parameter grid indices  
    grid_column_index <- (node_x_position / cellsize) + 1
    grid_row_index <- (node_y_position / cellsize) + 1
    
    # Extract subcatchment parameters for this headwater node
    FD_node <- OCN$RN$toFD[upstream_node]  # Corresponding flow direction node
    upstream_area <- OCN$RN$A[upstream_node]                         # Upstream drainage area [m2]
    mean_c_im_number <- OCN$FD$patches_concentration[FD_node]         # Mean immobile concentration parameter
    elevation <- OCN$FD$Z[FD_node]                                   # Elevation [m]
    
    # Compile fixed hydrological and biogeochemical parameters
    fixed_params <- list(
      gamma_sc = gamma_sc,
      sigma_w_sc = sigma_w_sc,
      oro_scaling = oro_scaling,
      tmax_yrs = tmax_yrs,
      rho = rho,
      damkohler_transport = damkohler_transport, 
      damkohler_interarrival = damkohler_interarrival,              # Damkohler number for interarrival time
      damkohler_longterm = damkohler_longterm,                      # Damkohler number for long-term processes
      interarrival_time_mean = interarrival_time_mean,              # Mean time between precipitation events [days]
      ET_max = ET_max, z_r = z_r, z_vz = z_vz,                     # Evapotranspiration and depth parameters
      theta_fc = theta_fc, theta_wp = theta_wp, theta_res = theta_res, # Soil moisture parameters
      R = R, mean_tr = mean_tr, theta_sat = theta_sat,              # Recharge and saturation parameters
      aquifer_z = aquifer_z, rain_per_year = rain_per_year         # Aquifer depth and annual rainfall
    )
    
    # Generate unique random seed for asynchrony parameter
    rho_seed <- sample(2^16, 1)
    
    # Random seed for precipitation timing patterns
    if (precip_timing == "random") {
      random_state_times <- sample(2^16, 1)   # Generate new random seed
    } else {
      random_state_times <- precip_timing            # Use fixed seed
    }
    
    # Random seed for precipitation intensity patterns  
    if (precip_intensity == "random") {
      random_state_rain <- sample(2^16, 1)     # Generate new random seed
    } else {
      random_state_rain <- precip_intensity          # Use fixed seed
    }
    
    # Random seed for travel time variability
    if (traveltime_seed == "random") {
      real_random <- sample(2^16, 1)           # Generate new random seed
    } else {
      real_random <- traveltime_seed                 # Use fixed seed
    }
    
    # Combine variable and fixed parameters for subcatchment model
    params <- list(
      mean_c_im_number = mean_c_im_number,   # Mean concentration parameter
      random_state_times = random_state_times, # Precipitation timing seed
      random_state_rain = random_state_rain,   # Precipitation intensity seed
      real_random = real_random,               # Travel time seed
      elevation = elevation,                   # Subcatchment elevation
      rho_seed = rho_seed      # Asynchrony random seed
    )
    
    params <- c(params, fixed_params)
    
    # Call Python subcatchment-scale model to generate time series
    result <- subcatchment_scale_module(param_dict = params)
    
    # Process discharge time series and convert units
    discharge <- result[[1]]                        # Raw discharge [cm/day]  
    discharge <- discharge / 100                    # Convert cm/day to m/day
    discharge <- discharge / (24 * 60 * 60)        # Convert m/day to m/s per unit area  
    discharge <- discharge * upstream_area          # Scale by drainage area to get total m3/s
    
    # Process concentration time series
    concentration <- result[[2]]                    # Concentration time series
    
    # Handle NaN values by replacing with mean concentration
    concentration[is.nan(concentration)] <- mean(concentration, na.rm = TRUE)
    
    
    # Accumulate volume input for mass balance (excluding warm-up period)
    volume_in <- volume_in + sum(discharge[warm_up:length(discharge)] * 3600 * 24)
    
    # Accumulate mass input for mass balance (excluding warm-up period) 
    mass_in <- mass_in + sum(discharge[warm_up:length(discharge)] * 3600 * 24 * concentration[warm_up:length(concentration)])
    
    # Store time series in network-wide matrices
    df_nodes_conc[, upstream_node] <- concentration
    df_nodes_discharge[, upstream_node] <- discharge
  }
  
  # ==================================================================
  # STREAM ORDER PROCESSING AND NETWORK ROUTING
  # ==================================================================
  
  # Process network by stream order (1st order → higher orders)
  # This ensures upstream nodes are processed before downstream nodes
  for (stream_order in c(1:max(OCN$AG$streamOrder))) {
    
    # Identify all nodes belonging to the current stream order
    stream_order_node_list <- get_stream_order_node_list(OCN, stream_order_dict, stream_order)
    
    # Identify upstream boundary nodes for the current stream order
    upstream_node_list <- get_upstream_node_list(OCN, stream_order_node_list, stream_order_dict, stream_order)
    
    # ==================================================================  
    # CONFLUENCE PROCESSING (Stream orders > 1 only)
    # ==================================================================
    
    # Process confluence points where multiple streams of lower order join
    # First-order streams don't have confluences (they are headwaters)
    if (stream_order > 1) {
      for (node in upstream_node_list) {
        # Identify all streams flowing into this confluence
        upstream_nodes <- which(OCN$RN$downNode == node)
        
        # Apply hydraulic routing to transport discharge and concentration
        # from upstream nodes to confluence, accounting for travel time and dispersion
        result_routing <- routing_routine(
          OCN = OCN,
          upstream_nodes = upstream_nodes,
          downstream_node = node,
          df_nodes_discharge = df_nodes_discharge,
          df_nodes_conc = df_nodes_conc,
          cellsize = cellsize,
          vf = vf,
          mannings_n = mannings_n,
          side_slope = side_slope,
          routing_function_ocn = routing_function_ocn
        )
        
        # Extract routing results
        df_discharge_upstream_routed <- result_routing[[1]]  # Routed discharge time series
        df_conc_upstream_routed <- result_routing[[2]]       # Routed concentration time series
        median_depth <- result_routing[[3]]                  # Median water depth for hydraulics
        median_celerity <- result_routing[[4]]               # Median flow celerity for hydraulics
        
        # Store hydraulic properties
        depth_dict[node] <- median_depth
        celerity_dict[node] <- median_celerity
        
        # Apply mixing routine to combine routed flows at confluence
        result_mixing <- mixing_routine(
          downstream_node = node,
          df_discharge_upstream_routed = df_discharge_upstream_routed,
          df_conc_upstream_routed = df_conc_upstream_routed,
          df_nodes_discharge = df_nodes_discharge,
          df_nodes_conc = df_nodes_conc
        )
        
        # Update network-wide time series with mixed results
        df_nodes_discharge <- result_mixing[[1]]
        df_nodes_conc <- result_mixing[[2]]
        
        # Add lateral inflow if this confluence receives subcatchment drainage
        if (node %in% lateral_inflow_nodes) {
          # Generate lateral inflow from subcatchment using stochastic model
          # This adds discharge and concentration from the surrounding landscape
          result_lateral <- lateral_inflow_routine(
            OCN = OCN,
            df_nodes_discharge = df_nodes_discharge,
            df_nodes_conc = df_nodes_conc,
            node = node,
            cellsize = cellsize,
            tmax_yrs = tmax_yrs,
            warm_up = warm_up,
            fixed_params = fixed_params,
            volume_in = volume_in,
            mass_in = mass_in,
            subcatchment_scale_module = subcatchment_scale_module,
            precip_timing,
            precip_intensity,
            traveltime_seed
          )

          # Update time series and mass balance with lateral inputs
          df_nodes_discharge <- result_lateral[[1]]
          df_nodes_conc <- result_lateral[[2]]
          volume_in <- result_lateral[[3]]
          mass_in <- result_lateral[[4]]
        }
      }
    }
    
    # For every upstream node of the repsective Strahler order, recursively iterate 
    # through the reach and populate downstream nodes
    pb3 <- progress_bar$new(
      format = glue("  Calculating stream order: {stream_order} [:bar] :percent eta: :eta"),
      total = length(upstream_node_list), clear = FALSE, width = 60
    )
    for (node in upstream_node_list) {
      pb3$tick()
      
      result_pop_down <- populate_downstream(
        OCN = OCN,
        node = node,
        df_nodes_discharge = df_nodes_discharge,
        df_nodes_conc = df_nodes_conc,
        cellsize = cellsize,
        landscape_configuration_matrix = landscape_configuration_matrix,
        lateral_inflow_nodes = lateral_inflow_nodes,
        tmax_yrs = tmax_yrs,
        warm_up = warm_up,
        vf = vf,
        mannings_n = mannings_n,
        side_slope = side_slope,
        fixed_params = fixed_params,
        volume_in = volume_in,
        mass_in = mass_in,
        precip_timing = precip_timing,
        precip_intensity = precip_intensity,
        traveltime_seed = traveltime_seed,
        subcatchment_scale_module = subcatchment_scale_module,
        routing_function_ocn = routing_function_ocn,
        depth_dict = depth_dict,
        celerity_dict = celerity_dict
      )
      df_nodes_discharge <- result_pop_down[[1]]
      df_nodes_conc <- result_pop_down[[2]]
      volume_in <- result_pop_down[[3]]
      mass_in <- result_pop_down[[4]]
      
      depth_dict <- result_pop_down[[5]]
      celerity_dict <- result_pop_down[[6]]
    }
  }
  
  # ==================================================================
  # MASS BALANCE CALCULATION AND LOGGING
  # ==================================================================
  
  # Identify network outlet for mass balance calculations
  outlet_node <- which(OCN$RN$downNode == 0)            # Main outlet node
  outlet_node2 <- which(OCN$RN$downNode == outlet_node) # Node downstream of outlet
  
  # Log simulation diagnostics
  log_warning_parallel(warning(glue("Length of time series: {nrow(df_nodes_conc)}")), log_file)
  
  # Extract outlet time series (excluding warm-up period) for mass balance
  outlet_node_conc <- df_nodes_conc[warm_up:nrow(df_nodes_conc), outlet_node]
  outlet_node_discharge <- df_nodes_discharge[warm_up:nrow(df_nodes_discharge), outlet_node]
  
  # Calculate total volume and mass outputs at network outlet
  volume_out <- sum(outlet_node_discharge * 3600 * 24)                        # Convert m3/s to m3/day, sum over time
  mass_out <- sum(outlet_node_discharge * outlet_node_conc * 3600 * 24)       # Mass flux = discharge × concentration
  
  # Compute mass balance errors as percentage of inputs
  percentage_error_mass <- mass_out / mass_in * 100
  percentage_error_volume <- volume_out / volume_in * 100
  
  # Log mass balance results for model validation
  log_warning_parallel(warning(glue("mass in: {mass_in}, mass out: {mass_out}, percent error: {percentage_error_mass}")), log_file)
  log_warning_parallel(warning(glue("\nvolume in: {volume_in}, volume out: {volume_out}, percent error: {percentage_error_volume}")), log_file)
  log_warning_parallel(warning(glue("Simulation Identifier: {simulation_identifier}")), log_file)
  log_warning_parallel("", log_file)
  
  # Remove warm-up period from final time series for analysis
  df_nodes_conc <- df_nodes_conc[warm_up:(tmax_yrs * 365), ]
  df_nodes_discharge <- df_nodes_discharge[warm_up:(tmax_yrs * 365), ]
  
  # Create dataframe with paramters of this simulation
  df_params <- data.frame(
    simulation_identifier = simulation_identifier,
    gamma_sc = gamma_sc,
    sigma_w_sc = sigma_w_sc,
    gamma_ls = gamma_ls,
    sigma_w_ls = sigma_w_ls,
    interarrival_time_mean = interarrival_time_mean,
    rho = rho,
    vf = vf,
    damkohler_transport = damkohler_transport,
    tmax_yrs = tmax_yrs,
    warm_up = warm_up,
    overall_mean_concentration = overall_mean_concentration,
    oro_scaling = oro_scaling,
    damkohler_interarrival = damkohler_interarrival,
    damkohler_longterm = damkohler_longterm,
    rain_per_year = rain_per_year,
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
    mannings_n = mannings_n,
    side_slope = side_slope,
    volume_in = volume_in,
    volume_out = volume_out,
    mass_in = mass_in,
    mass_out = mass_out,
    percentage_error_mass = percentage_error_mass,
    percentage_error_volume = percentage_error_volume,
    outlet_node = outlet_node2[[1]],
    network = network
  )
  
  # Ensure dataframes are in correct format for analysis
  df_nodes_conc <- as.data.frame(df_nodes_conc)
  df_nodes_discharge <- as.data.frame(df_nodes_discharge)
  
  # ==================================================================
  # SAVE RESULTS AS PARQUET and RDATA FILES
  # ==================================================================
  
  # save df_params
  write_parquet(df_params, file.path(path_to_output, glue("data/{simulation_identifier}_params.parquet")))
  
  # save df_nodes_conc
  write_parquet(df_nodes_conc, file.path(path_to_output, glue("data/{simulation_identifier}_conc_wide.parquet")))
  
  # save df_nodes_discharge
  write_parquet(df_nodes_discharge, file.path(path_to_output, glue("data/{simulation_identifier}_discharge_wide.parquet")))
  
  # Save OCN object
  save(OCN, file = file.path(path_to_output, glue("data/{simulation_identifier}_OCN.RData")))
  
  # save depth_dict as pickle
  saveRDS(depth_dict, file.path(path_to_output, glue("data/{simulation_identifier}_depth_dict.rds")))
  
  # save celerity_dict as pickle
  saveRDS(celerity_dict, file.path(path_to_output, glue("data/{simulation_identifier}_celerity_dict.rds")))
  
  # Return NULL to prevent doFuture from printing last variable
  NULL
}
# ===== END OF run_SLM FUNCTION =====