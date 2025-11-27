# =============================================================================
# Stochastic landscape model (SLM)
# =============================================================================


#' Create lookup dictionary mapping river network nodes to their Strahler stream orders
#'
#' @param OCN River network object
#'
#' @return List where each element contains stream order for the corresponding node index
get_stream_order_dict <- function(OCN) {
  # Initialize empty list
  stream_order_dict <- vector("list", OCN$RN$nNodes[1])

  # Iterate through each river network (RN) node to build the lookup dictionary
  for (node_index in seq_len(OCN$RN$nNodes[1])) {
    reach <- OCN$RN$toAGReach[node_index]
    stream_order_dict[[node_index]] <- OCN$AG$streamOrder[reach]
  }
  return(stream_order_dict)
}


#' Generates matrix with downstream distances to outlet
#'
#' @param OCN River network object
#'
#' @return Matrix of downstream path lengths
get_downstream_matrix <- function(OCN) {
  cellsize <- OCN$cellsize

  grid_width_cells <- max(OCN$FD$X) / cellsize + 1
  grid_height_cells <- max(OCN$FD$Y) / cellsize + 1

  # Create matrix to store downstream path lengths to outlet
  df_downstream_path <- matrix(0, grid_width_cells, grid_height_cells)

  # Populate grid with downstream path lengths
  for (node in 1:OCN$RN$nNodes[1]) {
    downstream_path <- OCN$RN$downstreamPathLength[node, OCN$RN$outlet]
    node_x_position <- OCN$RN$X[node]
    node_y_position <- OCN$RN$Y[node]

    grid_column_index <- (node_x_position / cellsize) + 1
    grid_row_index <- (node_y_position / cellsize) + 1
    df_downstream_path[grid_column_index, grid_row_index] <- downstream_path

    # Special case for outlet node: set small positive value to avoid zero division
    if (node == OCN$RN$outlet) {
      df_downstream_path[grid_column_index, grid_row_index] <- cellsize / 2
    }
  }
  return(df_downstream_path)
}


#' Finds source nodes and confluence points that mark stream order reach beginnings
#'
#' @param OCN River network object
#' @param stream_order_node_list List of all nodes in the current stream order
#' @param stream_order_dict Dictionary mapping nodes to their stream orders
#' @param stream_order Current stream order being processed
#'
#' @return List of upstream boundary nodes (sources and confluences)
get_upstream_node_list <- function(OCN, stream_order_node_list, stream_order_dict, stream_order) {
  upstream_node_list <- list() # Nodes at start of stream order reaches
  upstream_node_counter <- 1

  for (node in stream_order_node_list) {
    # Find all nodes that flow directly into the current node
    upstream_nodes <- which(OCN$RN$downNode == node)

    # Case 1: Source node with no upstream connections (first-order streams only)
    if (length(upstream_nodes) == 0) {
      upstream_node_list[upstream_node_counter] <- node
      upstream_node_counter <- upstream_node_counter + 1
    }

    # Case 2: Confluence node where all upstream tributaries are from different stream orders
    upstream_stream_orders <- stream_order_dict[upstream_nodes]
    if (length(upstream_nodes) != 0) {
      if (all(upstream_stream_orders != stream_order)) {
        upstream_node_list[upstream_node_counter] <- node
        upstream_node_counter <- upstream_node_counter + 1
      }
    }
  }
  return(upstream_node_list)
}


#' Collect all nodes belonging to a specific stream order
#'
#' @param OCN River network object
#' @param stream_order_dict Dictionary mapping nodes to their stream orders
#' @param stream_order Current stream order being processed
#'
#' @return List of node IDs belonging to the specified stream order
get_stream_order_node_list <- function(OCN, stream_order_dict, stream_order) {
  stream_order_node_counter <- 1
  stream_order_node_list <- list()

  # Iterate through all river network nodes
  for (node in 1:OCN$RN$nNodes[1]) {
    stream_order_node <- stream_order_dict[[node]]

    # Add node to list if it matches target stream order
    if (stream_order_node == stream_order) {
      stream_order_node_list[stream_order_node_counter] <- node
      stream_order_node_counter <- stream_order_node_counter + 1
    }
  }
  return(stream_order_node_list)
}


#' Recursively find all upstream flow direction (FD) nodes
#'
#' @param OCN River network object
#' @param FD_node Starting flow direction node ID
#' @param all_upstream_FD_nodes Accumulator vector for upstream node IDs
#'
#' @return Vector of all upstream FD node IDs
get_all_upstream_FD_nodes <- function(OCN,
                                      FD_node,
                                      all_upstream_FD_nodes) {
  # Find all Flow Direction (FD) nodes that drain directly into the current FD node
  upstream_FD_nodes <- which(OCN$FD$downNode == FD_node)

  all_upstream_FD_nodes <- c(all_upstream_FD_nodes, upstream_FD_nodes)

  if (length(upstream_FD_nodes) == 0) {
    # Base case: no more upstream nodes found, return complete drainage network
    return(all_upstream_FD_nodes)
  } else {
    for (node_index in seq_along(upstream_FD_nodes)) {
      current_node <- upstream_FD_nodes[node_index]
      # Recursively collect all nodes upstream of the current node
      all_upstream_FD_nodes <- get_all_upstream_FD_nodes(
        OCN = OCN,
        FD_node = current_node,
        all_upstream_FD_nodes = all_upstream_FD_nodes
      )
    }
  }
  return(all_upstream_FD_nodes)
}


#' Generate lateral inflow nodes for OCN river network
#'
#' @param OCN River network object
#' @param cellsize Grid cell size in meters
#' @param thresh Area threshold for subcatchment creation
#' @param landscape_configuration_matrix Matrix of landscape heterogeneity (mean immobile concentration of subcatchments)
#' @param precip_timing Seed for precipitation timing parameter or "random"
#' @param precip_intensity Seed for precipitation intensity parameter or "random"
#' @param traveltime_seed Seed for travel time parameter or "random"
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
  # Initialize Flow Direction (FD) grid patches for subcatchment assignment
  OCN$FD$patches <- matrix(0, OCN$FD$nNodes, 1)
  OCN$FD$patches_concentration <- matrix(0, OCN$FD$nNodes, 1)

  stream_order_dict <- get_stream_order_dict(OCN)

  # Initialize container for nodes that will receive lateral inflow
  lateral_inflow_nodes <- list()

  upstream <- OCN$RN$upstream
  upstream <- upstream[sapply(upstream, length) == 1]

  # Track FD nodes already assigned to subcatchments
  accounted_for_FD_nodes <- list()

  # Process river network by stream order
  for (stream_order in c(1:max(OCN$AG$streamOrder))) {
    # Identify all nodes belonging to the current stream order
    stream_order_node_list <- get_stream_order_node_list(OCN, stream_order_dict, stream_order)

    # Identify upstream boundary nodes for the current stream order
    upstream_node_list <- get_upstream_node_list(OCN, stream_order_node_list, stream_order_dict, stream_order)

    stream_order_progress_bar <- progress_bar$new(
      format = glue("  Lateral inflow nodes stream order: {stream_order} [:bar] :percent eta: :eta"),
      total = length(upstream_node_list), clear = FALSE, width = 60
    )

    for (node in upstream_node_list) {
      stream_order_progress_bar$tick()

      next_down_strahler_order <- stream_order
      next_down <- node

      # Move downstream through nodes of same stream order, creating subcatchments
      while (next_down_strahler_order == stream_order) {
        # Find corresponding Flow Direction (FD) grid cell for current river node
        FD_node <- OCN$RN$toFD[next_down]
        upstream_FD_nodes <- list()

        # Collect all FD cells that drain to this river node
        upstream_FD_nodes <- get_all_upstream_FD_nodes(
          OCN = OCN,
          FD_node = FD_node,
          all_upstream_FD_nodes = upstream_FD_nodes
        )

        # Subtract FD cells already assigned to other upstream subcatchments
        upstream_nodes_minus_accounted_for <- upstream_FD_nodes[!upstream_FD_nodes %in% accounted_for_FD_nodes]

        # Calculate drainage area for potential new subcatchment
        area <- cellsize * cellsize * length(upstream_nodes_minus_accounted_for)

        # Create subcatchment only if drainage area exceeds minimum threshold
        if (area > thresh) {
          # Register this river node as a lateral inflow point
          lateral_inflow_nodes <- c(lateral_inflow_nodes, next_down)

          # Mark all FD cells in this subcatchment as assigned
          accounted_for_FD_nodes <- c(accounted_for_FD_nodes, upstream_nodes_minus_accounted_for)

          # Get spatial coordinates of the outlet node for this subcatchment
          node_x_position <- OCN$RN$X[next_down] # X coordinate [m]
          node_y_position <- OCN$RN$Y[next_down] # Y coordinate [m]

          # Convert spatial coordinates to parameter grid indices
          grid_column_index <- (node_x_position / cellsize) + 1 # Grid column index
          grid_row_index <- (node_y_position / cellsize) + 1 # Grid row index

          # Extract subcatchment immobile concentration
          mean_c_im_number <- landscape_configuration_matrix[grid_column_index, grid_row_index] # Mean concentration parameter

          # Assign subcatchment parameters to all FD cells draining to this river node
          for (FD_node_acc in upstream_nodes_minus_accounted_for) {
            OCN$FD$patches[FD_node_acc] <- next_down # Assign outlet node ID
            OCN$FD$patches_concentration[FD_node_acc] <- mean_c_im_number # Mean concentration
          }
        }

        # Advance to the next downstream node or handle network outlet
        if (OCN$RN$downNode[next_down] != 0) {
          # Move to next downstream node within the river network
          next_down <- OCN$RN$downNode[next_down]
          next_down_strahler_order <- stream_order_dict[[next_down]] # Update stream order for loop control
        } else {
          # Reached network outlet - create final subcatchment and exit while loop
          lateral_inflow_nodes <- c(lateral_inflow_nodes, next_down)
          accounted_for_FD_nodes <- c(accounted_for_FD_nodes, upstream_nodes_minus_accounted_for)

          # Extract spatial parameters for the outlet subcatchment
          node_x_position <- OCN$RN$X[next_down] # Outlet X coordinate [m]
          node_y_position <- OCN$RN$Y[next_down] # Outlet Y coordinate [m]

          # Convert to parameter grid indices
          grid_column_index <- (node_x_position / cellsize) + 1
          grid_row_index <- (node_y_position / cellsize) + 1

          # Extract landscape parameters for outlet subcatchment
          mean_c_im_number <- landscape_configuration_matrix[grid_column_index, grid_row_index] # Mean concentration parameter

          for (FD_node_acc in upstream_nodes_minus_accounted_for) {
            OCN$FD$patches[FD_node_acc] <- next_down # Assign to outlet node
            OCN$FD$patches_concentration[FD_node_acc] <- mean_c_im_number # Mean concentration
          }

          # Assign parameters to the FD cell corresponding to the outlet node itself
          FD_outlet <- OCN$RN$toFD[next_down] # Get FD cell for outlet node
          OCN$FD$patches[FD_outlet] <- next_down # Assign to itself
          OCN$FD$patches_concentration[FD_outlet] <- mean_c_im_number # Mean concentration

          # Exit the downstream traversal loop since outlet is reached
          break
        }
      }
    }
  }
  return(list(
    lateral_inflow_nodes, # [1] Vector of river nodes receiving lateral inflow from subcatchments
    OCN # [2] Updated OCN object with FD patch assignments and parameters
  ))
}

#' Route discharge and concentration from upstream nodes to downstream node
#'
#' @param OCN River network object
#' @param upstream_nodes Vector of source node IDs
#' @param downstream_node Target node ID
#' @param df_nodes_discharge Discharge time series dataframe
#' @param df_nodes_conc Concentration time series dataframe
#' @param cellsize Grid cell size in meters
#' @param vf Mass transfer parameter
#' @param mannings_n Manning's roughness coefficient
#' @param side_slope Channel side slope
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
  # Extract time series for all upstream nodes that flow into the downstream target
  df_discharge_upstream_routed <- df_nodes_discharge[, upstream_nodes]
  df_conc_upstream_routed <- df_nodes_conc[, upstream_nodes]

  if (!is.data.frame(df_discharge_upstream_routed)) {
    df_discharge_upstream_routed <- as.data.frame(df_discharge_upstream_routed)
  }

  if (!is.data.frame(df_conc_upstream_routed)) {
    df_conc_upstream_routed <- as.data.frame(df_conc_upstream_routed)
  }

  for (node_index in seq_along(upstream_nodes)) {
    # Extract current tributary information for routing calculations
    upstream_node <- upstream_nodes[node_index] # Current tributary node ID
    inflow <- df_discharge_upstream_routed[, node_index] # Tributary discharge time series
    conc <- df_conc_upstream_routed[, node_index] # Tributary concentration time series

    # Calculate routing distance from tributary to confluence point
    reach_length <- OCN$RN$downstreamPathLength[upstream_node, downstream_node]
    routing_distance <- reach_length

    # Call Python-based hydraulic routing function to simulate water and solute transport
    routing_results <- routing_function_ocn(
      inflow = inflow, # Input discharge time series [m3/s]
      conc = conc, # Input concentration time series [mg/l]
      dx = routing_distance, # Routing distance [m]
      dt_ref = 24, # Reference time step [hours]
      bottom_slope = OCN$RN$slope[upstream_node], # Channel bottom slope [dimensionless]
      mannings_n = mannings_n, # Manning's roughness coefficient [s/m^(1/3)]
      reach_length = reach_length, # reach length [m]
      bottom_width_int = OCN$RN$width[upstream_node], # Channel bottom width [m]
      side_slope = side_slope, # Channel side slope [horizontal:vertical]
      vf = vf # mass transfer parameter [m/d]
    )

    # Extract routing outputs from Python function results
    routed_discharge <- routing_results[[1]]
    routed_concentration <- routing_results[[2]]

    median_water_depth <- routing_results[[3]]
    median_flow_celerity <- routing_results[[4]]

    # Store routed time series back into tributary-specific columns
    df_discharge_upstream_routed[, node_index] <- unlist(routed_discharge)
    df_conc_upstream_routed[, node_index] <- unlist(routed_concentration)
  }
  return(list(
    df_discharge_upstream_routed, # [1] Routed discharge time series for all tributaries [m3/s]
    df_conc_upstream_routed, # [2] Routed concentration time series for all tributaries [mg/l]
    median_water_depth, # [3] Median depth
    median_flow_celerity # [4] Median celerity
  ))
}

#' Mix flows and concentrations at confluence
#'
#' @param downstream_node Target node ID
#' @param df_discharge_upstream_routed Dataframe of routed discharge inputs
#' @param df_conc_upstream_routed Dataframe of routed concentration inputs
#' @param df_nodes_discharge Discharge time series dataframe
#' @param df_nodes_conc Concentration time series dataframe
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

  return(list(
    df_nodes_discharge, # [1] Updated discharge time series matrix [m3/s]
    df_nodes_conc # [2] Updated concentration time series matrix [mg/l]
  ))
}

#' Generate lateral inflow time series for a subcatchment
#'
#' @param OCN River network object
#' @param df_nodes_discharge Discharge time series dataframe
#' @param df_nodes_conc Concentration time series dataframe
#' @param node Target node ID
#' @param cellsize Grid cell size in meters
#' @param tmax_yrs Simulation time in years
#' @param warm_up Warm-up period
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
  # Get the corresponding Flow Direction (FD) node for the river network node
  FD_node <- OCN$RN$toFD[node]

  # Identify the subcatchment (patch) ID assigned to this FD node
  patch_id <- OCN$FD$patches[FD_node]

  # Find all FD cells that belong to this subcatchment
  all_upstream_FD_nodes <- which(OCN$FD$patches == patch_id)

  lateral_inflow <- rep(0, length(df_nodes_discharge[, node])) # Lateral discharge time series
  lateral_load <- rep(0, length(df_nodes_discharge[, node])) # Lateral load time series

  # Calculate total subcatchment drainage area
  patch_area <- length(all_upstream_FD_nodes) * cellsize * cellsize

  mean_c_im_number <- OCN$FD$patches_concentration[FD_node] # Mean concentration parameter [mg/l]
  elevation <- OCN$FD$Z[FD_node] # Elevation at subcatchment outlet [m]

  # Generate random seed for asynchrony parameter
  rho_seed <- sample(2^16, 1) # Random seed for asynchrony processes

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

  params <- list(
    mean_c_im_number = mean_c_im_number, # Mean concentration parameter
    random_state_times = random_state_times, # Precipitation timing seed
    random_state_rain = random_state_rain, # Precipitation intensity seed
    real_random = real_random, # Travel time seed
    elevation = elevation, # Subcatchment elevation
    rho_seed = rho_seed # Asynchrony random seed
  )

  # Combine subcatchment-specific parameters with fixed hydrological parameters
  params <- c(params, fixed_params)

  # Call Python-based subcatchment model to generate discharge and concentration time series
  result <- subcatchment_scale_module(param_dict = params)

  # Extract discharge time series
  discharge <- result[[1]]
  discharge <- discharge / 100
  discharge <- discharge / (24 * 60 * 60)
  discharge <- discharge * patch_area

  # Accumulate volume input for system-wide mass balance (excluding warm-up period)
  volume_in <- volume_in + sum(discharge[warm_up:length(discharge)] * 3600 * 24) # Convert m3/s to m3/day and sum

  concentration <- result[[2]] # Concentration time series [mg/l]

  lateral_inflow <- unlist(lapply(seq_along(lateral_inflow), function(i) lateral_inflow[i] + discharge[i]))

  lateral_load <- unlist(lapply(seq_along(lateral_load), function(i) lateral_load[i] + (discharge[i] * concentration[i])))

  mass_in <- mass_in + sum(lateral_load[warm_up:length(lateral_load)] * 3600 * 24)

  df_nodes_load <- df_nodes_discharge[, node] * df_nodes_conc[, node]

  df_nodes_discharge[, node] <- df_nodes_discharge[, node] + lateral_inflow

  df_nodes_load <- df_nodes_load + lateral_load

  df_nodes_conc[, node] <- df_nodes_load / df_nodes_discharge[, node]

  return(list(
    df_nodes_discharge,
    df_nodes_conc,
    volume_in,
    mass_in
  ))
}

#' Recursively populate downstream nodes with routed discharge and concentration
#'
#' @param OCN River network object
#' @param node Current node ID
#' @param df_nodes_discharge Discharge time series dataframe
#' @param df_nodes_conc Concentration time series dataframe
#' @param cellsize Grid cell size in meters
#' @param landscape_configuration_matrix Matrix of landscape heterogeneity (mean immobile concentration of subcatchments)
#' @param lateral_inflow_nodes Vector of lateral inflow node IDs
#' @param tmax_yrs Simulation time in years
#' @param warm_up Warm-up period
#' @param vf Mass transfer parameter
#' @param mannings_n Manning's roughness coefficient
#' @param side_slope Channel side slope
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
  # Get lookup dictionary mapping river network (RN) nodes to their Strahler stream orders
  stream_order_dict <- get_stream_order_dict(OCN)

  # Identify the immediate downstream node from current position
  downstream_node <- OCN$RN$downNode[node]

  if (downstream_node != 0) {
    stream_order_node <- stream_order_dict[[node]]
    stream_order_downstream <- stream_order_dict[[downstream_node]]

    # Terminate recursion if downstream node belongs to different stream order
    if (stream_order_node != stream_order_downstream) {
      result_list <- list(df_nodes_discharge, df_nodes_conc, volume_in, mass_in, depth_dict, celerity_dict)
      return(result_list)
    }

    # Identify all nodes that flow directly into the downstream target node
    upstream_nodes <- which(OCN$RN$downNode == downstream_node)

    # Filter upstream nodes to include only those of equal or lower stream order
    upstream_nodes <- upstream_nodes[stream_order_dict[upstream_nodes] <= stream_order_node]

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
    df_discharge_upstream_routed <- result_routing[[1]]
    df_conc_upstream_routed <- result_routing[[2]]
    median_depth <- result_routing[[3]]
    median_celerity <- result_routing[[4]]

    depth_dict[downstream_node] <- median_depth
    celerity_dict[downstream_node] <- median_celerity

    result_mixing <- mixing_routine(
      downstream_node = downstream_node,
      df_discharge_upstream_routed = df_discharge_upstream_routed,
      df_conc_upstream_routed = df_conc_upstream_routed,
      df_nodes_discharge = df_nodes_discharge,
      df_nodes_conc = df_nodes_conc
    )

    df_nodes_discharge <- result_mixing[[1]] # Updated discharge matrix
    df_nodes_conc <- result_mixing[[2]] # Updated concentration matrix

    # Check if downstream node receives lateral inflow from surrounding landscape
    if (downstream_node %in% lateral_inflow_nodes) {
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
      df_nodes_discharge <- result_lateral[[1]]
      df_nodes_conc <- result_lateral[[2]]
      volume_in <- result_lateral[[3]]
      mass_in <- result_lateral[[4]]
    }

    # Continue processing further downstream if not at network outlet
    if (OCN$RN$downNode[downstream_node] != 0) {
      # Recursively call this function to continue processing downstream
      result_pop_down <- populate_downstream(
        OCN = OCN,
        node = downstream_node, # Start from current downstream node
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
        volume_in = volume_in, # Pass updated mass balance
        mass_in = mass_in,
        precip_timing = precip_timing,
        precip_intensity = precip_intensity,
        traveltime_seed = traveltime_seed,
        subcatchment_scale_module = subcatchment_scale_module,
        routing_function_ocn = routing_function_ocn,
        depth_dict = depth_dict, # Pass hydraulic property storage
        celerity_dict = celerity_dict
      )

      # Extract results from recursive call
      df_nodes_discharge <- result_pop_down[[1]]
      df_nodes_conc <- result_pop_down[[2]]
      volume_in <- result_pop_down[[3]]
      mass_in <- result_pop_down[[4]]
      depth_dict <- result_pop_down[[5]]
      celerity_dict <- result_pop_down[[6]]
    }
  }

  result_list <- list(
    df_nodes_discharge,
    df_nodes_conc,
    volume_in,
    mass_in,
    depth_dict,
    celerity_dict
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
  grid_width_cells <- nrow(df_downstream_path)
  grid_height_cells <- ncol(df_downstream_path)

  landscape_configuration_matrix <- matrix(0, nrow = grid_width_cells, ncol = grid_height_cells)

  mean_conc_list <- c()

  for (row_index in 1:grid_width_cells) {
    for (col_index in 1:grid_height_cells) {
      if (df_downstream_path[row_index, col_index] != 0) {
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
  return(landscape_configuration_matrix)
}
# ===== END OF generate_landscape_configuration FUNCTION =====


#' Calculate immobile zone concentration based on flow distance
#'
#' @param distance Downstream path length
#' @param gamma_ls Scaling exponent parameter
#' @param sigma_w_ls Log-normal standard deviation parameter
#'
#' @return Numeric concentration value
calc_c_immobile <- function(distance, gamma_ls, sigma_w_ls) {
  # Generate random variability from log-normal distribution
  random_variable <- rlnorm(1, meanlog = 0, sdlog = sigma_w_ls)

  c_im <- (distance^gamma_ls) * random_variable

  return(c_im)
}
# ===== END OF calc_c_immobile FUNCTION =====


#' Main function to run the stochastic landscape model
#'
#' @param network OCN file path identifier
#' @param gamma_sc Subcatchment biogeochemical scaling parameter
#' @param sigma_w_sc Subcatchment biogeochemical variability parameter
#' @param gamma_ls Network-scale scaling parameter
#' @param sigma_w_ls Network-scale variability parameter
#' @param interarrival_time_mean Mean interarrival time
#' @param rho Asynchrony parameter for landscape heterogeneity
#' @param vf Mass transfer parameter
#' @param damkohler_transport Transport Damkohler number
#' @param tmax_yrs Simulation time in years
#' @param warm_up Warm-up period
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
  # Initialize mass and volume balance tracking variables
  volume_in <- 0
  volume_out <- 0
  mass_in <- 0
  mass_out <- 0

  # Configure progress bar display based on user preference
  options(progress_enabled = print_progress)

  # Setup Python environment
  use_python(path_to_python)
  source_python(path_to_python_functions)

  # Load the Optimal Channel Network (OCN) object from file
  load(glue(path_to_OCN))

  cellsize <- OCN$cellsize

  # Network name format: "X_Y_threshold_identifier"
  subcatchment_area_treshold <- as.numeric(strsplit(network, "_")[[1]][3])

  thresh <- subcatchment_area_treshold * cellsize^2

  stream_order_dict <- get_stream_order_dict(OCN)

  # Spatial grid matrix to store downstream path lengths to outlet
  df_downstream_path <- get_downstream_matrix(OCN)

  # Generate matrix with immobile concentrations based on parameters gamma_ls and sigma_w_ls
  landscape_configuration_matrix <- generate_landscape_configuration(
    df_downstream_path = df_downstream_path,
    gamma_ls = gamma_ls,
    sigma_w_ls = sigma_w_ls,
    overall_mean_concentration = overall_mean_concentration
  )

  depth_dict <- matrix(0, OCN$RN$nNodes, 1) # Median water depths [m]
  celerity_dict <- matrix(0, OCN$RN$nNodes, 1) # Median flow celerities [m/s]

  # Configure random seed usage based on hydroclimatic scenario:
  precip_timing <- sample(2^16, 1) # Fixed seed for precipitation timing
  precip_intensity <- sample(2^16, 1) # Fixed seed for precipitation intensity
  traveltime_seed <- "random"


  # Generate lateral inflow nodes based on area threshold
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
  lateral_inflow_nodes <- result[[1]] # List of nodes receiving lateral inflow
  OCN <- result[[2]] # Updated OCN with immobile concentration assigned to subcatchments


  # Identify headwater (source) nodes: nodes with no upstream connections
  upstream <- OCN$RN$upstream
  upstream <- upstream[sapply(upstream, length) == 1]

  # Remove headwater nodes from lateral inflow nodes to avoid double-counting
  lateral_inflow_nodes <- setdiff(lateral_inflow_nodes, upstream)

  df_nodes_conc <- data.frame(matrix(NA, nrow = (tmax_yrs * 365), ncol = OCN$RN$nNodes[1]))
  df_nodes_discharge <- data.frame(matrix(NA, nrow = (tmax_yrs * 365), ncol = OCN$RN$nNodes[1]))

  # Set column names to correspond to node IDs
  colnames(df_nodes_conc) <- 1:OCN$RN$nNodes[1]
  colnames(df_nodes_discharge) <- 1:OCN$RN$nNodes[1]


  # Generate concentration and discharge time series for all headwater nodes
  pb <- progress_bar$new(
    format = "  Creating headwater time series [:bar] :percent eta: :eta",
    total = length(upstream), clear = FALSE, width = 60
  )

  for (upstream_index in seq_along(upstream)) {
    pb$tick()

    upstream_node <- upstream[[upstream_index]]

    # Extract subcatchment parameters for this headwater node
    FD_node <- OCN$RN$toFD[upstream_node] # Corresponding flow direction node
    upstream_area <- OCN$RN$A[upstream_node] # Upstream drainage area [m2]
    mean_c_im_number <- OCN$FD$patches_concentration[FD_node] # Mean immobile concentration parameter
    elevation <- OCN$FD$Z[FD_node] # Elevation [m]

    fixed_params <- list(
      gamma_sc = gamma_sc,
      sigma_w_sc = sigma_w_sc,
      oro_scaling = oro_scaling,
      tmax_yrs = tmax_yrs,
      rho = rho,
      damkohler_transport = damkohler_transport,
      damkohler_interarrival = damkohler_interarrival, # Damkohler number for interarrival time
      damkohler_longterm = damkohler_longterm, # Damkohler number for long-term processes
      interarrival_time_mean = interarrival_time_mean, # Mean time between precipitation events [days]
      ET_max = ET_max, z_r = z_r, z_vz = z_vz, # Evapotranspiration and depth parameters
      theta_fc = theta_fc, theta_wp = theta_wp, theta_res = theta_res, # Soil moisture parameters
      R = R, mean_tr = mean_tr, theta_sat = theta_sat, # Recharge and saturation parameters
      aquifer_z = aquifer_z, rain_per_year = rain_per_year # Aquifer depth and annual rainfall
    )

    rho_seed <- sample(2^16, 1)

    if (precip_timing == "random") {
      random_state_times <- sample(2^16, 1) # Generate new random seed
    } else {
      random_state_times <- precip_timing # Use fixed seed
    }

    if (precip_intensity == "random") {
      random_state_rain <- sample(2^16, 1) # Generate new random seed
    } else {
      random_state_rain <- precip_intensity # Use fixed seed
    }

    if (traveltime_seed == "random") {
      real_random <- sample(2^16, 1) # Generate new random seed
    } else {
      real_random <- traveltime_seed # Use fixed seed
    }

    # Combine variable and fixed parameters for subcatchment model
    params <- list(
      mean_c_im_number = mean_c_im_number, # Mean concentration parameter
      random_state_times = random_state_times, # Precipitation timing seed
      random_state_rain = random_state_rain, # Precipitation intensity seed
      real_random = real_random, # Travel time seed
      elevation = elevation, # Subcatchment elevation
      rho_seed = rho_seed # Asynchrony random seed
    )

    params <- c(params, fixed_params)

    # Call Python subcatchment-scale model to generate time series
    result <- subcatchment_scale_module(param_dict = params)

    discharge <- result[[1]]
    discharge <- discharge / 100
    discharge <- discharge / (24 * 60 * 60)
    discharge <- discharge * upstream_area

    concentration <- result[[2]] # Concentration time series

    # Handle NaN values by replacing with mean concentration
    concentration[is.nan(concentration)] <- mean(concentration, na.rm = TRUE)

    volume_in <- volume_in + sum(discharge[warm_up:length(discharge)] * 3600 * 24)

    mass_in <- mass_in + sum(discharge[warm_up:length(discharge)] * 3600 * 24 * concentration[warm_up:length(concentration)])

    df_nodes_conc[, upstream_node] <- concentration
    df_nodes_discharge[, upstream_node] <- discharge
  }

  # Process network by stream order (1st order → higher orders)
  for (stream_order in c(1:max(OCN$AG$streamOrder))) {
    stream_order_node_list <- get_stream_order_node_list(OCN, stream_order_dict, stream_order)

    upstream_node_list <- get_upstream_node_list(OCN, stream_order_node_list, stream_order_dict, stream_order)

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
        df_discharge_upstream_routed <- result_routing[[1]] # Routed discharge time series
        df_conc_upstream_routed <- result_routing[[2]] # Routed concentration time series
        median_depth <- result_routing[[3]] # Median water depth for hydraulics
        median_celerity <- result_routing[[4]] # Median flow celerity for hydraulics

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

        # Add lateral inflow
        if (node %in% lateral_inflow_nodes) {
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

  outlet_node <- which(OCN$RN$downNode == 0) # Main outlet node
  outlet_node2 <- which(OCN$RN$downNode == outlet_node) # Node downstream of outlet

  # Log simulation diagnostics
  log_warning_parallel(warning(glue("Length of time series: {nrow(df_nodes_conc)}")), log_file)

  outlet_node_conc <- df_nodes_conc[warm_up:nrow(df_nodes_conc), outlet_node]
  outlet_node_discharge <- df_nodes_discharge[warm_up:nrow(df_nodes_discharge), outlet_node]

  volume_out <- sum(outlet_node_discharge * 3600 * 24) # Convert m3/s to m3/day, sum over time
  mass_out <- sum(outlet_node_discharge * outlet_node_conc * 3600 * 24) # Mass flux = discharge × concentration

  percentage_error_mass <- mass_out / mass_in * 100
  percentage_error_volume <- volume_out / volume_in * 100

  log_warning_parallel(warning(glue("mass in: {mass_in}, mass out: {mass_out}, percent error: {percentage_error_mass}")), log_file)
  log_warning_parallel(warning(glue("\nvolume in: {volume_in}, volume out: {volume_out}, percent error: {percentage_error_volume}")), log_file)
  log_warning_parallel(warning(glue("Simulation Identifier: {simulation_identifier}")), log_file)
  log_warning_parallel("", log_file)

  df_nodes_conc <- df_nodes_conc[warm_up:(tmax_yrs * 365), ]
  df_nodes_discharge <- df_nodes_discharge[warm_up:(tmax_yrs * 365), ]

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

  df_nodes_conc <- as.data.frame(df_nodes_conc)
  df_nodes_discharge <- as.data.frame(df_nodes_discharge)

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
