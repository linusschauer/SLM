# -*- coding: utf-8 -*-
"""
Stochastic Landscape Model (SLM) - MCT Routing Functions

Implements the Muskingum-Cunge-Todini (MCT) routing algorithm
(https://doi.org/10.5194/hess-11-1645-2007) for channel routing of
discharge and solute transport through river networks.

Model implementation adapted after:
- https://doi.org/10.5281/zenodo.13128294
- https://github.com/hydpy-dev/hydpy

Author: schauer
Created: Fri Feb 23 09:40:07 2024
"""

from typing import List, Tuple

import numpy as np
import pandas as pd
from numba import jit

# Constants
MINIMUM_DISCHARGE = 1e-6  # Minimum discharge to avoid numerical issues [m³/s]
COURANT_NUMBER_THRESHOLD_HIGH = 0.95  # Courant number threshold for stability
COURANT_NUMBER_THRESHOLD_LOW = 0.5  # Courant number threshold for stability
MAX_ITER = 50  # Maximum iterations for MCT
SECONDS_PER_DAY = 86400  # Seconds per day
INFLOW_RATIO_THRESHOLD = 5.0  # Max allowed inflow ratio between timesteps for concentration stability


# ============================================================================
# HYDRAULIC EQUATIONS AND NUMERICAL METHODS
# ============================================================================


@jit(nopython=True, fastmath=True)
def stage_discharge_relationship(
        stage_y: float,
        bottom_width_int: float,
        mannings_n: float,
        bottom_slope: float,
        sinus_alpha: float,
        cot_alpha: float,
) -> float:
    """
    Calculate discharge for a given stage (Equation B2/C4 in Todini et al., 2007).

    Parameters
    ----------
    stage_y : float
        Water stage (depth) above channel bottom [m].
    bottom_width_int : float
        Bottom width of the trapezoidal channel [m].
    mannings_n : float
        Manning's roughness coefficient [-].
    bottom_slope : float
        Longitudinal slope of the channel bottom [-].
    sinus_alpha : float
        Sine of the bank slope angle [-].
    cot_alpha : float
        Cotangent of the bank slope angle (horizontal/vertical) [-].

    Returns
    -------
    float
        Discharge through the channel cross-section [m³/s].
    """
    # Calculate wetted area: A(y) = (bottom_width + stage * cot_alpha) * stage
    wetted_area = (bottom_width_int + stage_y * cot_alpha) * stage_y

    # Calculate wetted perimeter: P(y) = bottom_width + 2*stage/sin_alpha
    wetted_perimeter = bottom_width_int + ((2 * stage_y) / sinus_alpha)

    # Discharge based on Equation B2/C4
    discharge = (
            (1 / mannings_n)
            * (bottom_slope ** 0.5)
            * (wetted_area ** (5 / 3))
            * (1 / (wetted_perimeter ** (2 / 3)))
    )

    return discharge


@jit(nopython=True, fastmath=True)
def d_stage_discharge_relationship(
        stage_y: float,
        bottom_width_int: float,
        mannings_n: float,
        bottom_slope: float,
        sinus_alpha: float,
        cot_alpha: float,
) -> float:
    """
    Calculate the derivative of discharge with respect to stage (Equation C6 for c(y)
    and B6 for the derivative of Q(y) in Todini et al., 2007).

    This function computes the derivative dQ/dy of the stage-discharge relationship,
    which is used in the Newton-Raphson iteration for solving the stage given discharge.

    Parameters
    ----------
    stage_y : float
        Water stage (depth) above channel bottom [m].
    bottom_width_int : float
        Bottom width of the trapezoidal channel [m].
    mannings_n : float
        Manning's roughness coefficient [-].
    bottom_slope : float
        Longitudinal slope of the channel bottom [-].
    sinus_alpha : float
        Sine of the bank slope angle [-].
    cot_alpha : float
        Cotangent of the bank slope angle (horizontal/vertical) [-].

    Returns
    -------
    float
        Derivative of discharge with respect to stage [m²/s].
    """
    # Calculate hydraulic geometry components
    # A(y): wetted area
    wetted_area = (bottom_width_int + stage_y * cot_alpha) * stage_y

    # B(y): surface width
    surface_width = bottom_width_int + (2 * stage_y * cot_alpha)

    # P(y): wetted perimeter
    wetted_perimeter = bottom_width_int + ((2 * stage_y) / sinus_alpha)

    # Equation C6: Calculate celerity
    celerity = (
            (5 / 3)
            * (bottom_slope ** 0.5 / mannings_n)
            * (wetted_area ** (2 / 3) / wetted_perimeter ** (2 / 3))
            * (
                    1
                    - (4 / 5) * (wetted_area / (surface_width * wetted_perimeter * sinus_alpha))
            )
    )

    # Equation B6: Derivative of discharge with respect to stage
    return surface_width * celerity


@jit(nopython=True, fastmath=True)
def newton_raphson(
        q_target: float,
        initial_guess: float,
        bottom_width_int: float,
        mannings_n: float,
        bottom_slope: float,
        sinus_alpha: float,
        cot_alpha: float,
        tolerance: float = 0.0001,
        max_iterations: int = 100,
) -> float:
    """
    Solve for water stage given target discharge using Newton-Raphson method
    based on Equation B3 in Todini et al. (2007).

    Parameters
    ----------
    q_target : float
        Target discharge for which to find the corresponding stage [m³/s].
    initial_guess : float
        Initial guess for the water stage [m].
    bottom_width_int : float
        Bottom width of the trapezoidal channel [m].
    mannings_n : float
        Manning's roughness coefficient [-].
    bottom_slope : float
        Longitudinal slope of the channel bottom [-].
    sinus_alpha : float
        Sine of the bank slope angle [-].
    cot_alpha : float
        Cotangent of the bank slope angle (horizontal/vertical) [-].
    tolerance : float, optional
        Convergence tolerance for stage difference [m]. Default is 0.0001 m (0.1 mm).
    max_iterations : int, optional
        Maximum number of iterations before stopping. Default is 100.

    Returns
    -------
    float
        Water stage corresponding to target discharge [m].
    """
    stage_estimate = initial_guess

    for iteration in range(max_iterations):
        # Calculate residual: f(h) = Q(h) - Q_target
        discharge_residual = (
                stage_discharge_relationship(
                    stage_estimate,
                    bottom_width_int,
                    mannings_n,
                    bottom_slope,
                    sinus_alpha,
                    cot_alpha,
                )
                - q_target
        )

        # Calculate derivative: f'(h) = dQ/dh
        discharge_derivative = d_stage_discharge_relationship(
            stage_estimate,
            bottom_width_int,
            mannings_n,
            bottom_slope,
            sinus_alpha,
            cot_alpha,
        )

        # Newton-Raphson update: h_new = h - f(h)/f'(h)
        stage_new = stage_estimate - discharge_residual / discharge_derivative

        # Check convergence criterion
        if abs(stage_new - stage_estimate) < tolerance:
            return stage_new

        # Update estimate for next iteration
        stage_estimate = stage_new

    # Return best estimate if maximum iterations reached
    return stage_estimate


@jit(nopython=True, fastmath=True)
def compute_hydraulic_properties(
        stage_y: float,
        bottom_width_int: float,
        mannings_n: float,
        bottom_slope: float,
        sinus_alpha: float,
        cot_alpha: float,
) -> Tuple[float, float, float, float, float, float]:
    """
    Compute hydraulic geometry properties for a trapezoidal channel.

    Given a water stage and channel geometry, computes wetted area, surface
    width, wetted perimeter, velocity, celerity, and beta (Eqs C1-C7 in
    Todini et al. (2007)).

    Parameters
    ----------
    stage_y : float
        Water stage (depth) above channel bottom [m].
    bottom_width_int : float
        Bottom width of the trapezoidal channel [m].
    mannings_n : float
        Manning's roughness coefficient [-].
    bottom_slope : float
        Longitudinal slope of the channel bottom [-].
    sinus_alpha : float
        Sine of the bank slope angle [-].
    cot_alpha : float
        Cotangent of the bank slope angle (horizontal/vertical) [-].

    Returns
    -------
    wetted_area : float
        Wetted cross-sectional area [m²].
    surface_width : float
        Water surface width [m].
    wetted_perimeter : float
        Wetted perimeter [m].
    velocity : float
        Mean flow velocity [m/s].
    celerity : float
        Flood wave celerity [m/s].
    beta : float
        Ratio of celerity to velocity [-].
    """
    # Equation C1: A(y) - wetted area
    wetted_area = (bottom_width_int + stage_y * cot_alpha) * stage_y

    # Equation C2: B(y) - surface width
    surface_width = bottom_width_int + (2 * stage_y * cot_alpha)

    # Equation C3: P(y) - wetted perimeter
    wetted_perimeter = bottom_width_int + ((2 * stage_y) / sinus_alpha)

    # Equation C5: velocity
    velocity = (bottom_slope ** 0.5 / mannings_n) * (
            (wetted_area ** (2 / 3)) / wetted_perimeter ** (2 / 3)
    )

    # Equation C6: celerity
    celerity = (
            (5 / 3)
            * (bottom_slope ** 0.5 / mannings_n)
            * (wetted_area ** (2 / 3) / wetted_perimeter ** (2 / 3))
            * (
                    1
                    - (4 / 5) * (wetted_area / (surface_width * wetted_perimeter * sinus_alpha))
            )
    )

    # Equation C7: beta
    beta = celerity / velocity

    return wetted_area, surface_width, wetted_perimeter, velocity, celerity, beta


# ============================================================================
# REACH ROUTING FUNCTIONS
# ============================================================================


@jit(nopython=True, fastmath=True)
def reach_routing(
        inflow: np.ndarray,
        conc: np.ndarray,
        dt_ref: float,
        dx: float,
        bottom_slope: float,
        mannings_n: float,
        bottom_width_int: float,
        side_slope: float,
        vf: float,
        optimized_calc: bool = False,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Perform routing using Muskingum-Cunge-Todini method.

    This function routes water and solute flux through a trapezoidal channel reach using
    the Muskingum-Cunge-Todini approach based on Todini et al. (2007) (https://doi.org/10.5194/hess-11-1645-2007).
    Model implementation is adapted after https://doi.org/10.5281/zenodo.13128294 and https://github.com/hydpy-dev/hydpy.

    Parameters
    ----------
    inflow : np.ndarray
        Time series of inflow discharge values [m³/s].
    conc : np.ndarray
        Time series of inflow concentration values [mg/L].
    dt_ref : float
        Reference time step for the simulation [h].
    dx : float
        Spatial discretization step (reach length) [m].
    bottom_slope : float
        Longitudinal slope of the channel bottom [-].
    mannings_n : float
        Manning's roughness coefficient [-].
    bottom_width_int : float
        Bottom width of the trapezoidal channel [m].
    side_slope : float
        Side slope of the trapezoidal channel (horizontal:vertical) [-].
    vf : float
        Settling velocity for in-stream first-order loss [m/day].
        Loss rate is computed as ke = vf / depth / 86400 [1/s].
    optimized_calc : bool, optional
        If True, returns only essential outputs for computational efficiency.
        If False, returns all diagnostic variables. Default is False.

    Returns
    -------
    out_list : np.ndarray
        Time series of outflow discharge values [m³/s].
    conc_routed_list : np.ndarray
        Time series of outflow concentration values [mg/L].
    t_list : np.ndarray
        Time series of time steps [h].
    storage_list : np.ndarray
        Time series of storage values [m³] (NaN array if optimized_calc is True).
    in_list : np.ndarray
        Time series of inflow discharge values [m³/s] (NaN array if optimized_calc is True).
    conc_list : np.ndarray
        Time series of inflow concentration values [mg/L] (NaN array if optimized_calc is True).
    load_storage_list : np.ndarray
        Time series of storage load values [mg] (NaN array if optimized_calc is True).
    depth_list : np.ndarray
        Time series of flow depth values [m] (NaN array if optimized_calc is True).
    celerity_list : np.ndarray
        Time series of wave celerity values [m/s] (NaN array if optimized_calc is True).

    Notes
    -----
    The trapezoidal channel cross-section is defined by:
    - Bottom width: bottom_width_int
    - Side slopes: specified by side_slope parameter
    """

    # Array of possible time steps [h]
    dt_array = np.array(
        [
            1 / 60,
            5 / 60,
            6 / 60,
            10 / 60,
            12 / 60,
            15 / 60,
            20 / 60,
            30 / 60,
            1,
            2,
            3,
            4,
            6,
            8,
            12,
            24,
        ]
    )

    # Hard-coded channel properties; trapezoidal cross section
    sinus_alpha = 1 / (1 + side_slope ** 2) ** 0.5
    cot_alpha = side_slope / 1
    safety_factor = 0.9

    # Set initial dt to reference dt
    dt = dt_ref

    # Initialize counters
    i = 0.0
    counter = 0

    # Pre-allocate arrays (oversized, will be trimmed later)
    # Sized for 10 years of minute-resolution output.
    # Sufficient for up to 50 years of daily input at the minimum adaptive
    # timestep (5 min). Exceeding this raises an IndexError.
    len_empty_array = 10 * 365 * 24 * 60

    # Essential output arrays (always calculated)
    out_list = np.zeros(len_empty_array)
    conc_routed_list = np.zeros(len_empty_array)
    t_list = np.zeros(len_empty_array)

    # Additional diagnostic arrays (only if not optimized)
    if not optimized_calc:
        in_list = np.zeros(len_empty_array)
        storage_list = np.zeros(len_empty_array)
        conc_list = np.zeros(len_empty_array)
        load_storage_list = np.zeros(len_empty_array)
        depth_list = np.zeros(len_empty_array)
        celerity_list = np.zeros(len_empty_array)

    # Initialize first value of time series
    # Set outflow to inflow (initial condition)
    out_list[0] = inflow[0]
    in_t_minus_1 = inflow[0]
    conc_routed_list[0] = conc[0]
    t_list[0] = i

    if not optimized_calc:
        in_list[0] = inflow[0]
        conc_list[0] = conc[0]

    # Equation 46a: calculate reference discharge at time t
    qref1 = (inflow[0] + out_list[0]) / 2

    # Equation 47a: calculate stage at reference discharge using Newton-Raphson approach
    stageref1 = newton_raphson(
        qref1, 1, bottom_width_int, mannings_n, bottom_slope, sinus_alpha, cot_alpha
    )

    if not optimized_calc:
        depth_list[0] = stageref1

    # Compute hydraulic properties at initial conditions (Eqs C1-C7)
    wetted_arearef1, surface_widthref1, _, vref1, cref1, betaref1 = (
        compute_hydraulic_properties(
            stageref1, bottom_width_int, mannings_n, bottom_slope, sinus_alpha, cot_alpha
        )
    )

    # Calculate initial storage load based on wetted area
    load_storage_t_minus_1 = wetted_arearef1 * dx * conc[0]
    if not optimized_calc:
        storage_list[0] = wetted_arearef1 * dx
        load_storage_list[0] = storage_list[0] * conc[0]

    if not optimized_calc:
        celerity_list[0] = cref1

    # Equation 51a: Reynolds number
    reynoldref1 = qref1 / (betaref1 * surface_widthref1 * bottom_slope * cref1 * dx)

    # Calculate optimal dt based on dx and wave celerity
    initial_dt = safety_factor * ((dx / cref1) / 3600)

    # Select next smallest dt from dt_array
    initial_dt = max(min(dt_array), initial_dt)
    initial_dt = max([x for x in dt_array if x <= initial_dt])
    dt_next = initial_dt

    # Main routing loop: calculate time steps until end of simulation
    while i < (((len(inflow)) * dt_ref) - dt_ref):
        dt = dt_next
        dt_in_sec = dt * 3600

        # Ensure last time step is calculated correctly
        if i + dt > (((len(inflow)) * dt_ref) - dt_ref) - (1 / 60):
            dt = (((len(inflow)) * dt_ref) - dt_ref) - i

        in_t = in_t_minus_1
        out_t = out_list[counter]

        # Interpolate inflow and concentration at time t+dt
        in_t_1 = np.interp(((i + dt) / dt_ref), np.arange(0, len(inflow)), inflow)
        conc_in_t_1 = np.interp(((i + dt) / dt_ref), np.arange(0, len(conc)), conc)

        # Reduce dt if inflow changes too rapidly between timesteps.
        inflow_ratio = max(in_t_1, in_t) / max(min(in_t_1, in_t), MINIMUM_DISCHARGE)
        while inflow_ratio > INFLOW_RATIO_THRESHOLD and dt > min(dt_array):
            dt = max([x for x in dt_array if x < dt][-1], min(dt_array))
            dt_in_sec = dt * 3600
            in_t_1 = np.interp(((i + dt) / dt_ref), np.arange(0, len(inflow)), inflow)
            conc_in_t_1 = np.interp(((i + dt) / dt_ref), np.arange(0, len(conc)), conc)
            inflow_ratio = max(in_t_1, in_t) / max(min(in_t_1, in_t), MINIMUM_DISCHARGE)

        # Call function that calculates the next time step
        (
            out_t_1,
            dt_next,
            qref2,
            vref2,
            surface_widthref2,
            cref2,
            betaref2,
            courantref2,
            reynoldref2,
        ) = reach_routing_single_time_step(
            in_t=in_t,
            in_t_1=in_t_1,
            out_t=out_t,
            dt=dt,
            dt_ref=dt_ref,
            dt_array=dt_array,
            safety_factor=safety_factor,
            dx=dx,
            bottom_slope=bottom_slope,
            bottom_width_int=bottom_width_int,
            mannings_n=mannings_n,
            sinus_alpha=sinus_alpha,
            cot_alpha=cot_alpha,
            reynoldref1=reynoldref1,
            cref1=cref1,
            betaref1=betaref1,
        )

        stageref2 = newton_raphson(
            out_t_1,
            1,
            bottom_width_int,
            mannings_n,
            bottom_slope,
            sinus_alpha,
            cot_alpha,
        )

        # Equation 54: storage calculations
        storage_t_1 = (((1 - reynoldref2) * dt_in_sec) / (2 * courantref2)) * in_t_1 + (
                ((1 + reynoldref2) * dt_in_sec) / (2 * courantref2)
        ) * out_t_1

        volume_in_t_1 = in_t_1 * dt_in_sec
        volume_out_t_1 = out_t_1 * dt_in_sec

        load_storage_t = load_storage_t_minus_1  # Load in storage at time t

        # Concentration at (t + dt): implicit scheme
        conc_storage_t_1 = (
                (load_storage_t + volume_in_t_1 * conc_in_t_1)
                / (storage_t_1 + volume_out_t_1)
        )

        # In-stream 1st order loss (settling/decay)
        ke = vf / stageref2 / SECONDS_PER_DAY  # vf [m/day] / depth [m] → [1/s]
        conc_storage_t_1 = conc_storage_t_1 * np.exp(-ke * dt_in_sec)

        # Update load in storage from new concentration
        load_storage_t_1_loss = storage_t_1 * conc_storage_t_1

        # Store results for current time step
        conc_routed_list[counter + 1] = conc_storage_t_1
        out_list[counter + 1] = out_t_1
        t_list[counter + 1] = i + dt

        # Store additional diagnostics if requested
        if not optimized_calc:
            conc_list[counter + 1] = conc_in_t_1
            storage_list[counter + 1] = storage_t_1
            in_list[counter + 1] = in_t_1
            load_storage_list[counter + 1] = load_storage_t_1_loss
            depth_list[counter + 1] = stageref2
            celerity_list[counter + 1] = cref2

        # Update variables for next iteration
        in_t_minus_1 = in_t_1
        load_storage_t_minus_1 = load_storage_t_1_loss
        qref1 = qref2
        vref1 = vref2
        surface_widthref1 = surface_widthref2
        cref1 = cref2
        betaref1 = betaref2
        reynoldref1 = reynoldref2

        counter += 1
        i = i + dt

    # Trim arrays to actual length
    out_list = out_list[: counter + 1]
    conc_routed_list = conc_routed_list[: counter + 1]
    t_list = t_list[: counter + 1]

    if not optimized_calc:
        storage_list = storage_list[: counter + 1]
        in_list = in_list[: counter + 1]
        conc_list = conc_list[: counter + 1]
        load_storage_list = load_storage_list[: counter + 1]
        depth_list = depth_list[: counter + 1]
        celerity_list = celerity_list[: counter + 1]
    else:
        # Return NaN arrays for optimized calculation
        storage_list = np.array([np.nan])
        in_list = np.array([np.nan])
        conc_list = np.array([np.nan])
        load_storage_list = np.array([np.nan])
        depth_list = np.array([np.nan])
        celerity_list = np.array([np.nan])

    return (
        out_list,
        conc_routed_list,
        t_list,
        storage_list,
        in_list,
        conc_list,
        load_storage_list,
        depth_list,
        celerity_list,
    )


@jit(nopython=True, fastmath=True)
def reach_routing_single_time_step(
        in_t: float,
        in_t_1: float,
        out_t: float,
        dt: float,
        dt_ref: float,
        dt_array: np.ndarray,
        safety_factor: float,
        dx: float,
        bottom_slope: float,
        bottom_width_int: float,
        mannings_n: float,
        sinus_alpha: float,
        cot_alpha: float,
        reynoldref1: float,
        cref1: float,
        betaref1: float,
) -> Tuple[float, float, float, float, float, float, float, float, float]:
    """
    This function performs a single time step of the reach routing and is called by reach_routing function.

    Parameters
    ----------
    in_t : float
        inflow at time t.
    in_t_1 : float
        inflow at time t+dt.
    out_t : float
        outflow at time t.
    dt : float
        time step.
    dt_ref : float
        reference time step.
    dt_array : numpy.array
        array of time steps.
    safety_factor : float
        safety factor.
    dx : float
        spatial step.
    bottom_slope : float
        bottom slope.
    bottom_width_int : float
        bottom width.
    mannings_n : float
        Manning's n.
    sinus_alpha : float
        sinus of alpha.
    cot_alpha : float
        cotangent of alpha.
    reynoldref1 : float
        Reynolds number.
    cref1 : float
        celerity.
    betaref1 : float
        beta.

    Returns
    -------
    out_t_i : float
        outflow at time t+dt.
    dt_next : float
        next time step.
    qref2 : float
        reference discharge at t+dt.
    vref2 : float
        velocity at t+dt.
    surface_widthref2 : float
        surface width at t+dt.
    cref2 : float
        celerity at t+dt.
    betaref2 : float
        beta at t+dt.
    courantref2 : float
        Courant number at t+dt.
    reynoldref2 : float
        Reynolds number at t+dt.
    """

    # Set maximum iteration limit and initialize counters
    counter = 1
    # courant_acceptable: True while Courant number is within bounds.
    # Set to False when Courant triggers a dt adjustment, which exits
    # the convergence loop and returns the adjusted dt to the caller.
    courant_acceptable = True

    # Equation 45: initial guess estimate
    q_guess = out_t + (in_t_1 - in_t)
    last_guess = q_guess * 2  # Ensure at least one iteration

    # Convergence threshold and safeguards
    tresh = 0.003
    if abs(last_guess - q_guess) <= tresh:
        last_guess = last_guess + abs(last_guess - q_guess) + tresh

    # Iterative solution loop
    while (abs(last_guess - q_guess) > tresh) and courant_acceptable:
        last_guess = q_guess
        dt_in_sec = dt * 3600

        # Equation 50a: Courant number
        courantref1 = (cref1 / betaref1) * (dt_in_sec / dx)

        # repeat the computation of equations (46b), (47b), (48b), (49b), (50b), (51b), (52) and (53)
        # twice to eliminate influence of the first guess (equation 45)
        repeat = 1
        while repeat < 3:
            # Equation 46b: calculate reference discharge at t+dt
            qref2 = (in_t_1 + q_guess) / 2

            # Ensure discharge doesn't approach zero (numerical stability)
            if qref2 < MINIMUM_DISCHARGE:
                qref2 = MINIMUM_DISCHARGE

            # Calculate hydraulic properties for reference discharge
            # Equation 47b: calculate stage using Newton-Raphson approach
            stageref2 = newton_raphson(
                qref2, 1, bottom_width_int, mannings_n, bottom_slope, sinus_alpha, cot_alpha
            )

            # Compute hydraulic properties at reference discharge (Eqs C1-C7)
            _, surface_widthref2, _, vref2, cref2, betaref2 = (
                compute_hydraulic_properties(
                    stageref2, bottom_width_int, mannings_n, bottom_slope, sinus_alpha, cot_alpha
                )
            )

            # Equation 50b: Courant number
            courantref2 = (cref2 / betaref2) * (dt_in_sec / dx)

            # Equation 51b: Reynolds number
            reynoldref2 = qref2 / (betaref2 * surface_widthref2 * bottom_slope * cref2 * dx)

            # Equation 52: Muskingum-Cunge-Todini (MCT) parameters
            c0 = (-1 + courantref2 + reynoldref2) / (1 + courantref2 + reynoldref2)
            c1 = ((1 + courantref1 - reynoldref1) / (1 + courantref2 + reynoldref2)) * (
                    courantref2 / courantref1
            )
            c2 = ((1 - courantref1 + reynoldref1) / (1 + courantref2 + reynoldref2)) * (
                    courantref2 / courantref1
            )

            # Equation 53: calculate new discharge estimate
            q_guess = (c0 * in_t_1) + (c1 * in_t) + (c2 * out_t)

            repeat += 1

        counter += 1
        courant_acceptable = True

        # Ensure discharge doesn't approach zero
        if q_guess < MINIMUM_DISCHARGE:
            q_guess = MINIMUM_DISCHARGE

        # Force convergence if maximum iterations reached
        # in exemplary tests, this threshold was never reached
        if counter == MAX_ITER:
            last_guess = q_guess

        # Adaptive time stepping based on Courant number
        if courantref2 < COURANT_NUMBER_THRESHOLD_LOW:  # Courant too small, increase time step
            dt = safety_factor * ((dx / cref2) / 3600)
            courant_acceptable = False
        elif courantref2 > COURANT_NUMBER_THRESHOLD_HIGH:  # Courant too large, decrease time step
            dt = safety_factor * ((dx / cref2) / 3600)
            courant_acceptable = False

        # Select appropriate dt from available options
        if not courant_acceptable:
            dt = max(min(dt_array), dt)
            dt = max([n for n in dt_array if n <= dt])
            dt = min(dt_ref, dt)

        # Check convergence percentage
        percent = abs(100 - (q_guess / last_guess) * 100)
        if percent <= 1:
            last_guess = q_guess

    # Final results
    out_t_i = q_guess
    dt_next = dt

    return (
        out_t_i,
        dt_next,
        qref2,
        vref2,
        surface_widthref2,
        cref2,
        betaref2,
        courantref2,
        reynoldref2,
    )


# ============================================================================
# RESAMPLING AND NETWORK ROUTING
# ============================================================================


def resample_routed_lists(
        in_list: np.ndarray,
        t_list: np.ndarray,
        dt_ref: float) -> pd.DataFrame:
    """
    Resample the routed arrays to original frequency.

    This function interpolates routing results from variable time steps back to
    a regular time grid matching the original simulation frequency.

    Parameters
    ----------
    in_list : np.ndarray
        Array of values to resample (e.g., discharge, concentration).
    t_list : np.ndarray
        Array of time values corresponding to in_list [hours].
    dt_ref : float
        Reference time step for resampling [hours].

    Returns
    -------
    pd.DataFrame
        DataFrame with resampled values and datetime index starting from 2000-01-01.
    """

    # Generate original timestamp grid
    original_timestamps = np.arange(0, max(t_list) + dt_ref, dt_ref)

    # Interpolate values to original frequency
    in_list = np.interp(original_timestamps, t_list, in_list)

    # Create DataFrame with datetime index
    df = pd.DataFrame(in_list)

    df.index = pd.to_datetime(original_timestamps, unit="h", origin="2000-01-01")

    # Create regular date range and filter to match
    if dt_ref < 1:
        dt_min = int(dt_ref * 60)
        freq = f"{int(dt_min)}min"
    else:
        freq = f"{int(dt_ref)}h"

    daterange = pd.date_range(
        start=df.index[0], end=df.index[-1], freq=freq, inclusive="both"
    )
    df = df[df.index.isin(daterange)]

    return df


def routing_function_ocn(
        inflow: np.ndarray,
        conc: np.ndarray,
        dx: float,
        dt_ref: float,
        bottom_slope: float,
        mannings_n: float,
        reach_length: float,
        bottom_width_int: float,
        side_slope: float,
        vf: float,
) -> List:
    """
    This function is called in R and performs the routing using reach_routing function.

    Parameters
    ----------
    inflow : np.array
        Array of inflow values [m^3/s]
    conc : np.array
        Array of concentration values [mg/L]
    dx : float
        Spatial step [m]
    dt_ref : float
        Time step [h]
    bottom_slope : float
        Slope of the channel bottom
    mannings_n : float
        Manning's roughness coefficient
    reach_length : float
        Length of the reach [m]
    bottom_width_int : float
        Bottom width of the channel at the start of the reach
    side_slope : float
        Side slope of the channel
    vf : float
        Settling velocity for in-stream first-order loss [m/day]

    Returns
    -------
    outflow : np.array
        Array of outflow values [m^3/s]
    conc_routed : np.array
        Array of routed concentration values [mg/L]
    median_depth : float
        Median depth of the reach [m]
    median_celerity : float
        Median celerity of the reach [m/s]
    """
    # Convert inputs to numpy arrays and flatten
    inflow = np.array(inflow).flatten()
    conc = np.array(conc).flatten()

    # Calculate number of sub-reaches
    n_reaches = int(reach_length / dx)

    if n_reaches < 1:
        raise ValueError(
            f"reach_length ({reach_length} m) must be >= dx ({dx} m), "
            f"got n_reaches = {n_reaches}"
        )

    # Check for reach length discretization consistency
    if n_reaches * dx != reach_length:
        print(f"reach_length of {reach_length} is not {n_reaches * dx}")

    # Route through each sub-reach sequentially
    for i in range(n_reaches):
        # Perform reach routing for current sub-reach
        (
            out_list,
            conc_routed_list,
            t_list,
            storage_list,
            in_list,
            conc_list,
            load_storage_list,
            depth_list,
            celerity_list,
        ) = reach_routing(
            inflow=inflow,
            conc=conc,
            dt_ref=dt_ref,
            dx=dx,
            bottom_slope=bottom_slope,
            mannings_n=mannings_n,
            bottom_width_int=bottom_width_int,
            side_slope=side_slope,
            vf=vf,
            optimized_calc=False,
        )

        # Resample routing results to regular time grid
        outflow = resample_routed_lists(out_list, t_list, dt_ref)
        conc_routed = resample_routed_lists(conc_routed_list, t_list, dt_ref)
        depth = resample_routed_lists(depth_list, t_list, dt_ref)
        celerity = resample_routed_lists(celerity_list, t_list, dt_ref)

        # Calculate median hydraulic properties for reach
        median_depth = np.median(depth)
        median_celerity = np.median(celerity)

        # Update inflow for next sub-reach (if not last reach)
        if i < (n_reaches - 1):
            inflow = outflow.values.flatten()
            conc = conc_routed.values.flatten()

    # Return results as list for R interface
    return [outflow, conc_routed, median_depth, median_celerity]