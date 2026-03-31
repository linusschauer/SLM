# =============================================================================
# Stochastic landscape model (SLM)
# =============================================================================


#' Create lookup vector mapping river network nodes to their Strahler stream orders
#'
#' @param OCN River network object
#'
#' @return Integer vector where element i is the stream order of node i
get_stream_order_dict <- function(OCN) {
  OCN$AG$streamOrder[OCN$RN$toAGReach]
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
#' A reach-start is any node of the target stream order that has no same-order
#' node flowing into it (i.e. it is either a source or a confluence fed only
#' by lower-order tributaries).
#'
#' @param OCN River network object
#' @param stream_order_node_list Integer vector of all nodes in the current stream order
#' @param stream_order_dict Integer vector mapping nodes to their stream orders
#' @param stream_order Current stream order being processed
#'
#' @return Integer vector of upstream boundary nodes (sources and confluences)
get_upstream_node_list <- function(OCN, stream_order_node_list, stream_order_dict, stream_order) {
  # Downstream targets of all same-order nodes
  same_order_targets <- OCN$RN$downNode[stream_order_node_list]

  # Reach-starts are same-order nodes that no same-order node flows into
  stream_order_node_list[!(stream_order_node_list %in% same_order_targets)]
}


#' Collect all nodes belonging to a specific stream order
#'
#' @param OCN River network object (unused, retained for API compatibility)
#' @param stream_order_dict Integer vector mapping nodes to their stream orders
#' @param stream_order Current stream order being processed
#'
#' @return Integer vector of node IDs belonging to the specified stream order
get_stream_order_node_list <- function(OCN, stream_order_dict, stream_order) {
  which(stream_order_dict == stream_order)
}


#' Build upstream adjacency list from flow direction downstream node vector
#'
#' Pre-computes which FD nodes drain into each FD node, avoiding repeated
#' linear scans of OCN$FD$downNode during tree traversal.
#'
#' @param downNode Integer vector where downNode[i] is the downstream node of node i
#'
#' @return List where element [[i]] is an integer vector of nodes draining into node i
build_upstream_adjacency <- function(downNode) {
  n <- length(downNode)
  adj <- vector("list", n)
  for (i in seq_len(n)) {
    dn <- downNode[i]
    if (dn > 0 && dn <= n) {
      adj[[dn]] <- c(adj[[dn]], i)
    }
  }
  return(adj)
}


#' Resolve a seed parameter: return the fixed value, or draw a new random seed
#'
#' @param seed_value Either "random" (draw a new seed) or a numeric seed to reuse
#'
#' @return Integer seed value
resolve_seed <- function(seed_value) {
  if (identical(seed_value, "random")) sample(2^16, 1) else seed_value
}


#' Find all upstream flow direction (FD) nodes using iterative traversal
#'
#' @param FD_node Starting flow direction node ID
#' @param upstream_adj Upstream adjacency list (from build_upstream_adjacency)
#'
#' @return Integer vector of all upstream FD node IDs
get_all_upstream_FD_nodes <- function(FD_node, upstream_adj) {
  # Iterative depth-first traversal using explicit stack
  stack <- upstream_adj[[FD_node]]
  if (is.null(stack) || length(stack) == 0) {
    return(integer(0))
  }

  result <- vector("list", length(stack))
  result_count <- 0
  stack_pos <- length(stack)

  while (stack_pos > 0) {
    # Pop from stack
    current <- stack[stack_pos]
    stack_pos <- stack_pos - 1

    # Collect node
    result_count <- result_count + 1
    if (result_count > length(result)) {
      # Double capacity if needed
      result <- c(result, vector("list", length(result)))
    }
    result[[result_count]] <- current

    # Push upstream neighbors onto stack
    neighbors <- upstream_adj[[current]]
    if (!is.null(neighbors) && length(neighbors) > 0) {
      n_new <- length(neighbors)
      # Grow stack if needed
      if (stack_pos + n_new > length(stack)) {
        stack <- c(stack, integer(n_new * 2))
      }
      stack[(stack_pos + 1):(stack_pos + n_new)] <- neighbors
      stack_pos <- stack_pos + n_new
    }
  }

  return(unlist(result[seq_len(result_count)]))
}


#' Generate lateral inflow nodes for OCN river network
#'
#' @param OCN River network object
#' @param cellsize Grid cell size in meters
#' @param thresh Area threshold for subcatchment creation
#' @param landscape_configuration_matrix Matrix of landscape heterogeneity (mean immobile concentration of subcatchments)
#' @param stream_order_dict List mapping node indices to their Strahler stream orders
#'
#' @return List containing lateral_inflow_nodes and modified OCN
generate_lateral_inflow_nodes <- function(OCN,
                                          cellsize,
                                          thresh,
                                          landscape_configuration_matrix,
                                          stream_order_dict) {
  # Initialize Flow Direction (FD) grid patches for subcatchment assignment
  OCN$FD$patches <- matrix(0, OCN$FD$nNodes, 1)
  OCN$FD$patches_concentration <- matrix(0, OCN$FD$nNodes, 1)

  # Initialize container for nodes that will receive lateral inflow
  # Use list with counter to avoid O(N²) from repeated c() calls
  lateral_inflow_collector <- vector("list", OCN$RN$nNodes[1])
  lateral_inflow_count <- 0

  # Track FD nodes already assigned to subcatchments using logical index
  # for O(1) lookup and assignment instead of O(N) with %in%
  is_accounted_for <- logical(OCN$FD$nNodes)

  # Pre-compute upstream adjacency list for FD tree traversal (O(N) once,
  # avoids O(N) which() scan per node during traversal)
  upstream_adj <- build_upstream_adjacency(OCN$FD$downNode)

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

        # Collect all FD cells that drain to this river node
        upstream_FD_nodes_vec <- get_all_upstream_FD_nodes(
          FD_node = FD_node,
          upstream_adj = upstream_adj
        )

        # Subtract FD cells already assigned to other upstream subcatchments
        upstream_nodes_minus_accounted_for <- upstream_FD_nodes_vec[!is_accounted_for[upstream_FD_nodes_vec]]

        # Calculate drainage area for potential new subcatchment
        area <- cellsize * cellsize * length(upstream_nodes_minus_accounted_for)

        # Create subcatchment only if drainage area exceeds minimum threshold
        if (area > thresh) {
          # Register this river node as a lateral inflow point
          lateral_inflow_count <- lateral_inflow_count + 1
          lateral_inflow_collector[[lateral_inflow_count]] <- next_down

          # Mark all FD cells in this subcatchment as assigned
          is_accounted_for[upstream_nodes_minus_accounted_for] <- TRUE

          # Get spatial coordinates of the outlet node for this subcatchment
          node_x_position <- OCN$RN$X[next_down] # X coordinate [m]
          node_y_position <- OCN$RN$Y[next_down] # Y coordinate [m]

          # Convert spatial coordinates to parameter grid indices
          grid_column_index <- (node_x_position / cellsize) + 1 # Grid column index
          grid_row_index <- (node_y_position / cellsize) + 1 # Grid row index

          # Extract subcatchment immobile concentration
          mean_c_im_number <- landscape_configuration_matrix[grid_column_index, grid_row_index] # Mean concentration parameter

          # Assign subcatchment parameters to all FD cells draining to this river node
          OCN$FD$patches[upstream_nodes_minus_accounted_for] <- next_down
          OCN$FD$patches_concentration[upstream_nodes_minus_accounted_for] <- mean_c_im_number
        }

        # Advance to the next downstream node or handle network outlet
        if (OCN$RN$downNode[next_down] != 0) {
          # Move to next downstream node within the river network
          next_down <- OCN$RN$downNode[next_down]
          next_down_strahler_order <- stream_order_dict[[next_down]] # Update stream order for loop control
        } else {
          # Reached network outlet - create final subcatchment and exit while loop
          lateral_inflow_count <- lateral_inflow_count + 1
          lateral_inflow_collector[[lateral_inflow_count]] <- next_down
          is_accounted_for[upstream_nodes_minus_accounted_for] <- TRUE

          # Extract spatial parameters for the outlet subcatchment
          node_x_position <- OCN$RN$X[next_down] # Outlet X coordinate [m]
          node_y_position <- OCN$RN$Y[next_down] # Outlet Y coordinate [m]

          # Convert to parameter grid indices
          grid_column_index <- (node_x_position / cellsize) + 1
          grid_row_index <- (node_y_position / cellsize) + 1

          # Extract landscape parameters for outlet subcatchment
          mean_c_im_number <- landscape_configuration_matrix[grid_column_index, grid_row_index] # Mean concentration parameter

          OCN$FD$patches[upstream_nodes_minus_accounted_for] <- next_down
          OCN$FD$patches_concentration[upstream_nodes_minus_accounted_for] <- mean_c_im_number

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
  # Flatten collected lateral inflow nodes to a vector
  lateral_inflow_nodes <- unlist(lateral_inflow_collector[seq_len(lateral_inflow_count)])

  return(list(
    lateral_inflow_nodes = lateral_inflow_nodes,
    OCN = OCN
  ))
}

#' Route discharge and concentration from upstream nodes to downstream node
#'
#' @param OCN River network object
#' @param upstream_nodes Vector of source node IDs
#' @param downstream_node Target node ID
#' @param df_nodes_discharge Discharge time series dataframe
#' @param df_nodes_conc Concentration time series dataframe
#' @param channel_params List with vf, mannings_n, side_slope
#' @param routing_function_ocn Python routing function
#'
#' @return Named list: discharge, concentration, depth, celerity
routing_routine <- function(OCN,
                            upstream_nodes,
                            downstream_node,
                            df_nodes_discharge,
                            df_nodes_conc,
                            channel_params,
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

  # Accumulators for flow-weighted hydraulic properties across tributaries
  weighted_depth_sum <- 0
  weighted_celerity_sum <- 0
  total_discharge_weight <- 0

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
      mannings_n = channel_params$mannings_n, # Manning's roughness coefficient [s/m^(1/3)]
      reach_length = reach_length, # reach length [m]
      bottom_width_int = OCN$RN$width[upstream_node], # Channel bottom width [m]
      side_slope = channel_params$side_slope, # Channel side slope [horizontal:vertical]
      vf = channel_params$vf # mass transfer parameter [m/d]
    )

    # Extract routing outputs from Python function results
    routed_discharge <- routing_results[[1]]
    routed_concentration <- routing_results[[2]]

    # Accumulate flow-weighted hydraulic properties
    trib_mean_Q <- mean(unlist(routed_discharge))
    weighted_depth_sum <- weighted_depth_sum + routing_results[[3]] * trib_mean_Q
    weighted_celerity_sum <- weighted_celerity_sum + routing_results[[4]] * trib_mean_Q
    total_discharge_weight <- total_discharge_weight + trib_mean_Q

    # Store routed time series back into tributary-specific columns
    df_discharge_upstream_routed[, node_index] <- unlist(routed_discharge)
    df_conc_upstream_routed[, node_index] <- unlist(routed_concentration)
  }

  # Compute flow-weighted average hydraulic properties across all tributaries
  if (total_discharge_weight > 0) {
    median_water_depth <- weighted_depth_sum / total_discharge_weight
    median_flow_celerity <- weighted_celerity_sum / total_discharge_weight
  } else {
    median_water_depth <- 0
    median_flow_celerity <- 0
  }

  return(list(
    discharge = df_discharge_upstream_routed,
    concentration = df_conc_upstream_routed,
    depth = median_water_depth,
    celerity = median_flow_celerity
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

  # Calculate concentration from load and discharge, guarding against zero-flow timesteps
  mixed_Q <- df_nodes_discharge[, downstream_node]
  df_nodes_conc[, downstream_node] <- ifelse(mixed_Q == 0, 0, df_nodes_load / mixed_Q)

  return(list(
    discharge = df_nodes_discharge,
    concentration = df_nodes_conc
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
#' @param stochastic_headwater_model Python stochastic headwater model function
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
                                   stochastic_headwater_model,
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

  # Assemble subcatchment-specific parameters
  params <- c(list(
    mean_c_im_number = mean_c_im_number,
    random_state_times = resolve_seed(precip_timing),
    random_state_rain = resolve_seed(precip_intensity),
    real_random = resolve_seed(traveltime_seed),
    elevation = elevation,
    rho_seed = sample(2^16, 1)
  ), fixed_params)

  # Call Python-based stochastic headwater model to generate discharge and concentration time series
  result <- stochastic_headwater_model(param_dict = params)

  # Extract and convert discharge time series from SHM output
  discharge <- convert_shm_discharge(result[[1]], patch_area)

  # Accumulate volume input for system-wide mass balance (excluding warm-up period)
  volume_in <- volume_in + sum(discharge[warm_up:length(discharge)] * 3600 * 24) # Convert m3/s to m3/day and sum

  concentration <- result[[2]] # Concentration time series [mg/l]

  lateral_inflow <- discharge

  lateral_load <- discharge * concentration

  mass_in <- mass_in + sum(lateral_load[warm_up:length(lateral_load)] * 3600 * 24)

  df_nodes_load <- df_nodes_discharge[, node] * df_nodes_conc[, node]

  df_nodes_discharge[, node] <- df_nodes_discharge[, node] + lateral_inflow

  df_nodes_load <- df_nodes_load + lateral_load

  total_Q <- df_nodes_discharge[, node]
  df_nodes_conc[, node] <- ifelse(total_Q == 0, 0, df_nodes_load / total_Q)

  return(list(
    discharge = df_nodes_discharge,
    concentration = df_nodes_conc,
    volume_in = volume_in,
    mass_in = mass_in
  ))
}

#' Iteratively populate downstream nodes with routed discharge and concentration
#'
#' Walks downstream along a single stream-order reach from the starting node,
#' applying routing, mixing, and lateral inflow at each step. Terminates when
#' the reach ends (stream order changes, outlet reached, or downNode == 0).
#'
#' @param OCN River network object
#' @param node Starting node ID
#' @param df_nodes_discharge Discharge time series dataframe
#' @param df_nodes_conc Concentration time series dataframe
#' @param cellsize Grid cell size in meters
#' @param is_lateral_inflow Logical vector where TRUE indicates node receives lateral inflow
#' @param tmax_yrs Simulation time in years
#' @param warm_up Warm-up period
#' @param channel_params List with vf, mannings_n, side_slope
#' @param fixed_params Fixed model parameters
#' @param volume_in Volume mass balance input
#' @param mass_in Mass balance input
#' @param precip_timing Precipitation timing parameter or "random"
#' @param precip_intensity Precipitation intensity parameter or "random"
#' @param traveltime_seed Travel time seed parameter or "random"
#' @param stochastic_headwater_model Python stochastic headwater model function
#' @param routing_function_ocn Python routing function
#' @param depth_dict Dictionary for storing hydraulic depth properties
#' @param celerity_dict Dictionary for storing hydraulic celerity properties
#' @param stream_order_dict Integer vector mapping node indices to their Strahler stream orders
#'
#' @return Named list: discharge, concentration, volume_in, mass_in, depth_dict, celerity_dict
populate_downstream <- function(OCN,
                                node,
                                df_nodes_discharge,
                                df_nodes_conc,
                                cellsize,
                                is_lateral_inflow,
                                tmax_yrs,
                                warm_up,
                                channel_params,
                                fixed_params,
                                volume_in,
                                mass_in,
                                precip_timing,
                                precip_intensity,
                                traveltime_seed,
                                stochastic_headwater_model,
                                routing_function_ocn,
                                depth_dict,
                                celerity_dict,
                                stream_order_dict) {
  current_node <- node
  stream_order_current <- stream_order_dict[[current_node]]

  while (TRUE) {
    downstream_node <- OCN$RN$downNode[current_node]

    # Stop at network outlet
    if (downstream_node == 0) break

    # Stop at stream order boundary
    if (stream_order_dict[[downstream_node]] != stream_order_current) break

    # Identify all nodes that flow directly into the downstream target node
    upstream_nodes <- which(OCN$RN$downNode == downstream_node)

    # Filter upstream nodes to include only those of equal or lower stream order
    upstream_nodes <- upstream_nodes[stream_order_dict[upstream_nodes] <= stream_order_current]

    result_routing <- routing_routine(
      OCN = OCN,
      upstream_nodes = upstream_nodes,
      downstream_node = downstream_node,
      df_nodes_discharge = df_nodes_discharge,
      df_nodes_conc = df_nodes_conc,
      channel_params = channel_params,
      routing_function_ocn = routing_function_ocn
    )

    # Extract hydraulic routing results
    df_discharge_upstream_routed <- result_routing$discharge
    df_conc_upstream_routed <- result_routing$concentration

    depth_dict[downstream_node] <- result_routing$depth
    celerity_dict[downstream_node] <- result_routing$celerity

    result_mixing <- mixing_routine(
      downstream_node = downstream_node,
      df_discharge_upstream_routed = df_discharge_upstream_routed,
      df_conc_upstream_routed = df_conc_upstream_routed,
      df_nodes_discharge = df_nodes_discharge,
      df_nodes_conc = df_nodes_conc
    )

    df_nodes_discharge <- result_mixing$discharge
    df_nodes_conc <- result_mixing$concentration

    # Check if downstream node receives lateral inflow from surrounding landscape
    if (is_lateral_inflow[downstream_node]) {
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
        stochastic_headwater_model = stochastic_headwater_model,
        precip_timing,
        precip_intensity,
        traveltime_seed
      )

      df_nodes_discharge <- result_lateral$discharge
      df_nodes_conc <- result_lateral$concentration
      volume_in <- result_lateral$volume_in
      mass_in <- result_lateral$mass_in
    }

    # Advance to next node in reach
    current_node <- downstream_node
  }

  return(list(
    discharge = df_nodes_discharge,
    concentration = df_nodes_conc,
    volume_in = volume_in,
    mass_in = mass_in,
    depth_dict = depth_dict,
    celerity_dict = celerity_dict
  ))
}


#' Create matrix of the subcatchment immobile concentration based on gamma_ls and sigma_w_ls
#'
#' For each grid cell with a nonzero downstream path length, computes
#' C_im = d_eff^|gamma_ls| * W, where W ~ Lognormal with E[W] = 1.
#'
#' When gamma_ls >= 0, d_eff = distance (concentration increases with path length).
#' When gamma_ls <  0, d_eff = d_max + d_min - distance (concentration decreases
#' with path length, mirroring the positive-gamma curvature).
#'
#' The resulting matrix is then rescaled so the mean of nonzero cells
#' equals overall_mean_concentration.
#'
#' @param df_downstream_path Network topology dataframe
#' @param gamma_ls Network-wide scaling parameter (sign controls direction,
#'   magnitude controls curvature)
#' @param sigma_w_ls Network-wide variability parameter
#' @param overall_mean_concentration Overall mean concentration value
#'
#' @return Matrix of immobile concentration
generate_landscape_configuration <- function(
    df_downstream_path,
    gamma_ls,
    sigma_w_ls,
    overall_mean_concentration) {
  landscape_configuration_matrix <- matrix(0, nrow = nrow(df_downstream_path), ncol = ncol(df_downstream_path))

  # Identify nonzero cells (cells that belong to the network)
  nonzero_mask <- df_downstream_path != 0
  n_nonzero <- sum(nonzero_mask)

  # Compute effective distance: reverse axis for negative gamma
  distances <- df_downstream_path[nonzero_mask]
  if (gamma_ls < 0) {
    d_eff <- max(distances) + min(distances) - distances
  } else {
    d_eff <- distances
  }

  # C_im = d_eff^|gamma| * W, where W is unit-mean lognormal
  W <- rlnorm(n_nonzero, meanlog = -sigma_w_ls^2 / 2, sdlog = sigma_w_ls)
  landscape_configuration_matrix[nonzero_mask] <- (d_eff^abs(gamma_ls)) * W

  # Rescale to achieve target mean concentration
  calc_mean_conc <- mean(landscape_configuration_matrix[nonzero_mask])
  scaling_factor <- overall_mean_concentration / calc_mean_conc
  landscape_configuration_matrix[nonzero_mask] <- landscape_configuration_matrix[nonzero_mask] * scaling_factor

  return(landscape_configuration_matrix)
}
# ===== END OF generate_landscape_configuration FUNCTION =====


#' Convert SHM discharge output to volumetric discharge
#'
#' The stochastic headwater model returns discharge as specific discharge
#' in [cm/day]. This function converts to volumetric discharge [m³/s].
#'
#' @param discharge_cm_day Specific discharge from SHM [cm/day]
#' @param area Drainage area [m²]
#'
#' @return Volumetric discharge [m³/s]
convert_shm_discharge <- function(discharge_cm_day, area) {
  discharge_m_day <- discharge_cm_day / 100       # cm/day → m/day
  discharge_m_s <- discharge_m_day / (24 * 60 * 60) # m/day → m/s
  discharge_m3_s <- discharge_m_s * area           # m/s × m² → m³/s
  return(discharge_m3_s)
}


#' Initialize simulation: load OCN, build spatial lookups, delineate subcatchments,
#' allocate time series matrices.
#'
#' @param network OCN file path identifier (format: "X_Y_threshold_identifier")
#' @param path_to_OCN Path to OCN .RData file
#' @param gamma_ls Network-scale scaling parameter
#' @param sigma_w_ls Network-scale variability parameter
#' @param overall_mean_concentration Target mean immobile concentration
#' @param tmax_yrs Simulation time in years
#'
#' @return Named list with OCN, spatial lookups, time series matrices, and seed config
setup_simulation <- function(network, path_to_OCN, gamma_ls, sigma_w_ls,
                             overall_mean_concentration, tmax_yrs) {
  load(glue(path_to_OCN))
  cellsize <- OCN$cellsize

  subcatchment_area_treshold <- as.numeric(strsplit(network, "_")[[1]][3])
  thresh <- subcatchment_area_treshold * cellsize^2

  stream_order_dict <- get_stream_order_dict(OCN)
  df_downstream_path <- get_downstream_matrix(OCN)

  landscape_configuration_matrix <- generate_landscape_configuration(
    df_downstream_path = df_downstream_path,
    gamma_ls = gamma_ls,
    sigma_w_ls = sigma_w_ls,
    overall_mean_concentration = overall_mean_concentration
  )

  result <- generate_lateral_inflow_nodes(
    OCN = OCN,
    cellsize = cellsize,
    thresh = thresh,
    landscape_configuration_matrix = landscape_configuration_matrix,
    stream_order_dict = stream_order_dict
  )
  OCN <- result$OCN

  # Identify headwater (source) nodes: nodes with no upstream connections
  upstream <- OCN$RN$upstream
  upstream <- upstream[sapply(upstream, length) == 1]

  # Remove headwater nodes from lateral inflow to avoid double-counting
  lateral_inflow_nodes <- setdiff(result$lateral_inflow_nodes, upstream)
  is_lateral_inflow <- logical(OCN$RN$nNodes[1])
  is_lateral_inflow[lateral_inflow_nodes] <- TRUE

  # Allocate time series matrices
  n_nodes <- OCN$RN$nNodes[1]
  n_days <- tmax_yrs * 365
  df_nodes_conc <- data.frame(matrix(NA, nrow = n_days, ncol = n_nodes))
  df_nodes_discharge <- data.frame(matrix(NA, nrow = n_days, ncol = n_nodes))
  colnames(df_nodes_conc) <- seq_len(n_nodes)
  colnames(df_nodes_discharge) <- seq_len(n_nodes)

  return(list(
    OCN = OCN,
    cellsize = cellsize,
    stream_order_dict = stream_order_dict,
    is_lateral_inflow = is_lateral_inflow,
    upstream = upstream,
    depth_dict = matrix(0, n_nodes, 1),
    celerity_dict = matrix(0, n_nodes, 1),
    precip_timing = sample(2^16, 1),
    precip_intensity = sample(2^16, 1),
    traveltime_seed = "random",
    df_nodes_conc = df_nodes_conc,
    df_nodes_discharge = df_nodes_discharge
  ))
}


#' Generate discharge and concentration time series for all headwater nodes
#'
#' @param OCN River network object
#' @param upstream List of headwater node IDs
#' @param df_nodes_discharge Discharge time series matrix (modified in place)
#' @param df_nodes_conc Concentration time series matrix (modified in place)
#' @param fixed_params Fixed hydrological parameters for the SHM
#' @param precip_timing Precipitation timing seed or "random"
#' @param precip_intensity Precipitation intensity seed or "random"
#' @param traveltime_seed Travel time seed or "random"
#' @param warm_up Warm-up period in days
#' @param stochastic_headwater_model Python SHM function
#'
#' @return Named list: discharge, concentration, volume_in, mass_in
generate_headwater_timeseries <- function(OCN, upstream,
                                          df_nodes_discharge, df_nodes_conc,
                                          fixed_params,
                                          precip_timing, precip_intensity, traveltime_seed,
                                          warm_up,
                                          stochastic_headwater_model) {
  volume_in <- 0
  mass_in <- 0

  pb <- progress_bar$new(
    format = "  Creating headwater time series [:bar] :percent eta: :eta",
    total = length(upstream), clear = FALSE, width = 60
  )

  for (upstream_index in seq_along(upstream)) {
    pb$tick()
    upstream_node <- upstream[[upstream_index]]

    FD_node <- OCN$RN$toFD[upstream_node]
    upstream_area <- OCN$RN$A[upstream_node]
    mean_c_im_number <- OCN$FD$patches_concentration[FD_node]
    elevation <- OCN$FD$Z[FD_node]

    params <- c(list(
      mean_c_im_number = mean_c_im_number,
      random_state_times = resolve_seed(precip_timing),
      random_state_rain = resolve_seed(precip_intensity),
      real_random = resolve_seed(traveltime_seed),
      elevation = elevation,
      rho_seed = sample(2^16, 1)
    ), fixed_params)

    result <- stochastic_headwater_model(param_dict = params)

    discharge <- convert_shm_discharge(result[[1]], upstream_area)
    concentration <- result[[2]]
    concentration[is.nan(concentration)] <- mean(concentration, na.rm = TRUE)

    volume_in <- volume_in + sum(discharge[warm_up:length(discharge)] * 3600 * 24)
    mass_in <- mass_in + sum(discharge[warm_up:length(discharge)] * 3600 * 24 *
                               concentration[warm_up:length(concentration)])

    df_nodes_conc[, upstream_node] <- concentration
    df_nodes_discharge[, upstream_node] <- discharge
  }

  return(list(
    discharge = df_nodes_discharge,
    concentration = df_nodes_conc,
    volume_in = volume_in,
    mass_in = mass_in
  ))
}


#' Route flows through the network by Strahler stream order
#'
#' Processes confluences (route + mix + lateral inflow) then propagates
#' each reach downstream via populate_downstream, for each stream order
#' from 1 to max.
#'
#' @param OCN River network object
#' @param stream_order_dict Integer vector mapping nodes to Strahler orders
#' @param df_nodes_discharge Discharge time series matrix
#' @param df_nodes_conc Concentration time series matrix
#' @param is_lateral_inflow Logical vector for lateral inflow nodes
#' @param depth_dict Hydraulic depth storage matrix
#' @param celerity_dict Hydraulic celerity storage matrix
#' @param volume_in Cumulative volume input (from headwater phase)
#' @param mass_in Cumulative mass input (from headwater phase)
#' @param cellsize Grid cell size in meters
#' @param tmax_yrs Simulation time in years
#' @param warm_up Warm-up period in days
#' @param channel_params List with vf, mannings_n, side_slope
#' @param fixed_params Fixed hydrological parameters for the SHM
#' @param precip_timing Precipitation timing seed or "random"
#' @param precip_intensity Precipitation intensity seed or "random"
#' @param traveltime_seed Travel time seed or "random"
#' @param stochastic_headwater_model Python SHM function
#' @param routing_function_ocn Python routing function
#'
#' @return Named list: discharge, concentration, volume_in, mass_in, depth_dict, celerity_dict
route_network <- function(OCN, stream_order_dict,
                          df_nodes_discharge, df_nodes_conc,
                          is_lateral_inflow, depth_dict, celerity_dict,
                          volume_in, mass_in,
                          cellsize, tmax_yrs, warm_up,
                          channel_params, fixed_params,
                          precip_timing, precip_intensity, traveltime_seed,
                          stochastic_headwater_model, routing_function_ocn) {
  for (stream_order in seq_len(max(OCN$AG$streamOrder))) {
    stream_order_node_list <- get_stream_order_node_list(OCN, stream_order_dict, stream_order)
    upstream_node_list <- get_upstream_node_list(OCN, stream_order_node_list, stream_order_dict, stream_order)

    # Process confluences for higher-order streams
    if (stream_order > 1) {
      for (node in upstream_node_list) {
        upstream_nodes <- which(OCN$RN$downNode == node)

        result_routing <- routing_routine(
          OCN = OCN,
          upstream_nodes = upstream_nodes,
          downstream_node = node,
          df_nodes_discharge = df_nodes_discharge,
          df_nodes_conc = df_nodes_conc,
          channel_params = channel_params,
          routing_function_ocn = routing_function_ocn
        )

        depth_dict[node] <- result_routing$depth
        celerity_dict[node] <- result_routing$celerity

        result_mixing <- mixing_routine(
          downstream_node = node,
          df_discharge_upstream_routed = result_routing$discharge,
          df_conc_upstream_routed = result_routing$concentration,
          df_nodes_discharge = df_nodes_discharge,
          df_nodes_conc = df_nodes_conc
        )

        df_nodes_discharge <- result_mixing$discharge
        df_nodes_conc <- result_mixing$concentration

        if (is_lateral_inflow[node]) {
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
            stochastic_headwater_model = stochastic_headwater_model,
            precip_timing,
            precip_intensity,
            traveltime_seed
          )

          df_nodes_discharge <- result_lateral$discharge
          df_nodes_conc <- result_lateral$concentration
          volume_in <- result_lateral$volume_in
          mass_in <- result_lateral$mass_in
        }
      }
    }

    # Propagate each reach downstream
    pb <- progress_bar$new(
      format = glue("  Calculating stream order: {stream_order} [:bar] :percent eta: :eta"),
      total = length(upstream_node_list), clear = FALSE, width = 60
    )
    for (node in upstream_node_list) {
      pb$tick()

      result_pop_down <- populate_downstream(
        OCN = OCN,
        node = node,
        df_nodes_discharge = df_nodes_discharge,
        df_nodes_conc = df_nodes_conc,
        cellsize = cellsize,
        is_lateral_inflow = is_lateral_inflow,
        tmax_yrs = tmax_yrs,
        warm_up = warm_up,
        channel_params = channel_params,
        fixed_params = fixed_params,
        volume_in = volume_in,
        mass_in = mass_in,
        precip_timing = precip_timing,
        precip_intensity = precip_intensity,
        traveltime_seed = traveltime_seed,
        stochastic_headwater_model = stochastic_headwater_model,
        routing_function_ocn = routing_function_ocn,
        depth_dict = depth_dict,
        celerity_dict = celerity_dict,
        stream_order_dict = stream_order_dict
      )
      df_nodes_discharge <- result_pop_down$discharge
      df_nodes_conc <- result_pop_down$concentration
      volume_in <- result_pop_down$volume_in
      mass_in <- result_pop_down$mass_in
      depth_dict <- result_pop_down$depth_dict
      celerity_dict <- result_pop_down$celerity_dict
    }
  }

  return(list(
    discharge = df_nodes_discharge,
    concentration = df_nodes_conc,
    volume_in = volume_in,
    mass_in = mass_in,
    depth_dict = depth_dict,
    celerity_dict = celerity_dict
  ))
}


#' Write simulation outputs to disk
#'
#' Trims warm-up period from time series and writes all output files.
#'
#' @param df_nodes_conc Concentration time series matrix
#' @param df_nodes_discharge Discharge time series matrix
#' @param OCN River network object
#' @param depth_dict Hydraulic depth matrix
#' @param celerity_dict Hydraulic celerity matrix
#' @param df_params Parameter record dataframe
#' @param warm_up Warm-up period in days
#' @param tmax_yrs Simulation time in years
#' @param path_to_output Output directory path
#' @param simulation_identifier Unique simulation ID
write_simulation_output <- function(df_nodes_conc, df_nodes_discharge,
                                    OCN, depth_dict, celerity_dict,
                                    df_params, warm_up, tmax_yrs,
                                    path_to_output, simulation_identifier) {
  # Trim warm-up period
  df_nodes_conc <- as.data.frame(df_nodes_conc[warm_up:(tmax_yrs * 365)-1, ])
  df_nodes_discharge <- as.data.frame(df_nodes_discharge[warm_up:(tmax_yrs * 365)-1, ])

  write_parquet(df_params, file.path(path_to_output, glue("data/{simulation_identifier}_params.parquet")))
  write_parquet(df_nodes_conc, file.path(path_to_output, glue("data/{simulation_identifier}_conc_wide.parquet")))
  write_parquet(df_nodes_discharge, file.path(path_to_output, glue("data/{simulation_identifier}_discharge_wide.parquet")))
  save(OCN, file = file.path(path_to_output, glue("data/{simulation_identifier}_OCN.RData")))
  saveRDS(depth_dict, file.path(path_to_output, glue("data/{simulation_identifier}_depth_dict.rds")))
  saveRDS(celerity_dict, file.path(path_to_output, glue("data/{simulation_identifier}_celerity_dict.rds")))
}


#' Main orchestrator for the stochastic landscape model
#'
#' Coordinates setup, headwater generation, network routing, diagnostics,
#' and output writing. Called from parallel dispatch (doFuture).
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
#' @return NULL (saves results to files)
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
  options(progress_enabled = print_progress)
  use_python(path_to_python)
  source_python(path_to_python_functions)

  # --- Setup ---
  sim <- setup_simulation(network, path_to_OCN, gamma_ls, sigma_w_ls,
                          overall_mean_concentration, tmax_yrs)

  channel_params <- list(vf = vf, mannings_n = mannings_n, side_slope = side_slope)

  fixed_params <- list(
    gamma_sc = gamma_sc, sigma_w_sc = sigma_w_sc, oro_scaling = oro_scaling,
    tmax_yrs = tmax_yrs, rho = rho, damkohler_transport = damkohler_transport,
    damkohler_interarrival = damkohler_interarrival,
    damkohler_longterm = damkohler_longterm,
    interarrival_time_mean = interarrival_time_mean,
    ET_max = ET_max, z_r = z_r, z_vz = z_vz,
    theta_fc = theta_fc, theta_wp = theta_wp, theta_res = theta_res,
    R = R, mean_tr = mean_tr, theta_sat = theta_sat,
    aquifer_z = aquifer_z, rain_per_year = rain_per_year
  )

  # --- Generate headwater time series ---
  hw <- generate_headwater_timeseries(
    OCN = sim$OCN, upstream = sim$upstream,
    df_nodes_discharge = sim$df_nodes_discharge,
    df_nodes_conc = sim$df_nodes_conc,
    fixed_params = fixed_params,
    precip_timing = sim$precip_timing,
    precip_intensity = sim$precip_intensity,
    traveltime_seed = sim$traveltime_seed,
    warm_up = warm_up,
    stochastic_headwater_model = stochastic_headwater_model
  )

  # --- Route through network ---
  net <- route_network(
    OCN = sim$OCN, stream_order_dict = sim$stream_order_dict,
    df_nodes_discharge = hw$discharge, df_nodes_conc = hw$concentration,
    is_lateral_inflow = sim$is_lateral_inflow,
    depth_dict = sim$depth_dict, celerity_dict = sim$celerity_dict,
    volume_in = hw$volume_in, mass_in = hw$mass_in,
    cellsize = sim$cellsize, tmax_yrs = tmax_yrs, warm_up = warm_up,
    channel_params = channel_params, fixed_params = fixed_params,
    precip_timing = sim$precip_timing,
    precip_intensity = sim$precip_intensity,
    traveltime_seed = sim$traveltime_seed,
    stochastic_headwater_model = stochastic_headwater_model,
    routing_function_ocn = routing_function_ocn
  )

  # --- Mass balance diagnostics ---
  outlet_node <- which(sim$OCN$RN$downNode == 0)
  outlet_node2 <- which(sim$OCN$RN$downNode == outlet_node)

  log_warning_parallel(warning(glue("Length of time series: {nrow(net$concentration)}")), log_file)

  outlet_conc <- net$concentration[warm_up:nrow(net$concentration)-1, outlet_node]
  outlet_discharge <- net$discharge[warm_up:nrow(net$discharge)-1, outlet_node]

  volume_out <- sum(outlet_discharge * 3600 * 24)
  mass_out <- sum(outlet_discharge * outlet_conc * 3600 * 24)

  pct_mass <- mass_out / net$mass_in * 100
  pct_volume <- volume_out / net$volume_in * 100

  log_warning_parallel(warning(glue("mass in: {net$mass_in}, mass out: {mass_out}, percent error: {pct_mass}")), log_file)
  log_warning_parallel(warning(glue("\nvolume in: {net$volume_in}, volume out: {volume_out}, percent error: {pct_volume}")), log_file)
  log_warning_parallel(warning(glue("Simulation Identifier: {simulation_identifier}")), log_file)
  log_warning_parallel("", log_file)

  # --- Build parameter record ---
  df_params <- data.frame(
    simulation_identifier = simulation_identifier,
    gamma_sc = gamma_sc, sigma_w_sc = sigma_w_sc,
    gamma_ls = gamma_ls, sigma_w_ls = sigma_w_ls,
    interarrival_time_mean = interarrival_time_mean, rho = rho, vf = vf,
    damkohler_transport = damkohler_transport,
    tmax_yrs = tmax_yrs, warm_up = warm_up,
    overall_mean_concentration = overall_mean_concentration,
    oro_scaling = oro_scaling,
    damkohler_interarrival = damkohler_interarrival,
    damkohler_longterm = damkohler_longterm,
    rain_per_year = rain_per_year, ET_max = ET_max,
    z_r = z_r, z_vz = z_vz,
    theta_fc = theta_fc, theta_wp = theta_wp, theta_res = theta_res,
    R = R, mean_tr = mean_tr, theta_sat = theta_sat,
    aquifer_z = aquifer_z, mannings_n = mannings_n, side_slope = side_slope,
    volume_in = net$volume_in, volume_out = volume_out,
    mass_in = net$mass_in, mass_out = mass_out,
    percentage_error_mass = pct_mass,
    percentage_error_volume = pct_volume,
    outlet_node = outlet_node2[[1]],
    network = network
  )

  # --- Write output ---
  write_simulation_output(
    df_nodes_conc = net$concentration, df_nodes_discharge = net$discharge,
    OCN = sim$OCN, depth_dict = net$depth_dict, celerity_dict = net$celerity_dict,
    df_params = df_params, warm_up = warm_up, tmax_yrs = tmax_yrs,
    path_to_output = path_to_output, simulation_identifier = simulation_identifier
  )

  # Return NULL to prevent doFuture from printing last variable
  NULL
}
# ===== END OF run_SLM FUNCTION =====