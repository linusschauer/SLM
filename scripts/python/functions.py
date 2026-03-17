# -*- coding: utf-8 -*-
"""
Stochastic Landscape Model (SLM) - Core Functions

This script contains:
(1) Subcatchment-scale stochastic model of Musolff et al. (2017) (https://doi.org/10.1002/2017GL072630)
(2) MCT routing algorithm (https://doi.org/10.5194/hess-11-1645-2007)

Author: schauer
Created: Fri Feb 23 09:40:07 2024
"""

from typing import Dict, Tuple, Union, List

import numpy as np
import pandas as pd
import scipy.special
import scipy.stats
from numba import jit
from scipy.stats import expon, lognorm

# Constants
SPIN_UP_YEARS = 7  # Default spin-up time for subcatchment_scale_module in years
DAYS_PER_YEAR = 365  # Days per year
MINIMUM_DISCHARGE = 1e-6  # Minimum discharge to avoid numerical issues [m³/s]
COURANT_NUMBER_THRESHOLD_HIGH = 0.95  # Courant number threshold for stability
COURANT_NUMBER_THRESHOLD_LOW = 0.5  # Courant number threshold for stability
MAX_ITER = 50  # Maximum iterations for MCT

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def calculate_cumulative_rainfall_for_event(
    event_index: int, cumulative_rainfall: np.ndarray, max_time: int
) -> np.ndarray:
    """
    Calculate cumulative rainfall array for a specific event.

    This function extracts the relevant portion of cumulative rainfall for
    distributing solute masses among future rainfall events.

    Parameters
    ----------
    event_index : int
        Index of the rainfall event (0-based).
    cumulative_rainfall : np.ndarray
        Array of cumulative rainfall amounts over time.
    max_time : int
        Maximum simulation time steps.

    Returns
    -------
    np.ndarray
        Cumulative rainfall array starting from the event time.

    Notes
    -----
    For the first event (index 0), returns the full cumulative rainfall array.
    For subsequent events, returns the incremental rainfall from the event onwards.
    """
    if event_index == 0:
        rainfall_array = cumulative_rainfall
    else:
        # Calculate incremental rainfall from this event to max_time
        end_index = event_index + max_time
        rainfall_array = (
            cumulative_rainfall[event_index:end_index]
            - cumulative_rainfall[event_index - 1]
        )
    return rainfall_array


def calculate_event_mass_fractions(
    cumulative_rainfall: np.ndarray, mu_ln_streamtube: float, sigma_ln_streamtube: float
) -> np.ndarray:
    """
    Calculate mass fractions for distributing solute among future events.

    This function calculates how mass from one event is distributed among
    subsequent events based on the lognormal travel time distribution.

    Parameters
    ----------
    cumulative_rainfall : np.ndarray
        Cumulative rainfall amounts for the event sequence.
    mu_ln_streamtube : float
        Mean of log-normal streamtube length distribution.
    sigma_ln_streamtube : float
        Standard deviation of log-normal streamtube length distribution.

    Returns
    -------
    np.ndarray
        Array of mass fractions for each subsequent event.

    Notes
    -----
    Uses the cumulative distribution function (CDF) of a lognormal distribution
    to determine what fraction of mass from an event reaches each future event.
    The first element represents the fraction remaining at the source event.
    """
    # Calculate cumulative fractions using lognormal CDF
    cumulative_fractions = np.asarray(lognorm.cdf(
        x=cumulative_rainfall,
        s=sigma_ln_streamtube,
        scale=np.exp(mu_ln_streamtube),
    ))

    # Convert cumulative to incremental fractions
    incremental_fractions = cumulative_fractions.copy()
    incremental_fractions[1:] = incremental_fractions[1:] - cumulative_fractions[:-1]

    return incremental_fractions


# ============================================================================
# MAIN SIMULATION FUNCTIONS
# ============================================================================


def subcatchment_scale_module(
    param_dict: Dict[str, Union[float, int]],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    This function implements the stochastic model described in Musolff et al. (2017)
    (https://doi.org/10.1002/2017GL072630) to simulate discharge and solute
    transport in catchments with different solute source distributions.

    Parameters
    ----------
    param_dict : dict
        Dictionary containing all model parameters including:

        Random seeds:
            - random_state_times : int
                Seed for timing of rainfall events
            - random_state_rain : int
                Seed for rainfall intensity
            - real_random : int
                General random seed
            - rho_seed : int
                Seed for hydroclimatic asynchrony

        Model parameters:
            - gamma_sc : float
                Structured heterogeneity (correlation of solute sources to travel time)
            - sigma_w_sc : float
                Random heterogeneity parameter
            - damkohler_interarrival : float
                Damköhler number for inter-arrival times (Da1)
            - damkohler_transport : float
                Damköhler number for transport (Da2)
            - damkohler_longterm : float
                Damköhler number for long-term changes (Da3)
            - oro_scaling : float
                Orographic scaling factor for precipitation
            - tmax_yrs : float
                Simulation time in years
            - rain_per_year : float
                Mean annual rainfall [cm/yr]
            - interarrival_time_mean : float
                Mean rainfall inter-arrival time [days]
            - ET_max : float
                Maximum evapotranspiration [cm/day]
            - elevation : float
                Catchment elevation [m]
            - rho : float
                Asynchrony parameter (0-1)
            - z_r : float
                Root zone thickness [cm]
            - z_vz : float
                Vadose zone thickness [cm]
            - theta_fc : float
                Water content at field capacity [-]
            - theta_wp : float
                Water content at wilting point [-]
            - theta_res : float
                Residual water content [-]
            - theta_sat : float
                Saturated water content [-]
            - R : float
                Retardation factor
            - mean_tr : float
                Hydraulic response time [days]
            - aquifer_z : float
                Aquifer thickness [cm]
            - mean_c_im_number : float
                Mean immobile concentration [mg/L]

    Returns
    -------
    tuple of np.ndarray
        discharge : np.ndarray
            Daily discharge time series [m³/s or mm/day]
        concentration : np.ndarray
            Daily concentration time series [mg/L]
        mass_flux : np.ndarray
            Daily mass flux time series [mg/L × m³/s]

    Notes
    -----
    The simulation includes a 7-year spin-up period which is removed from outputs.
    """

    # Set random number generator seeds
    r_nonzero_times = np.random.RandomState(
        int(param_dict["random_state_times"])
    )  # Timing of rain events
    r_nonzero_rain = np.random.RandomState(
        int(param_dict["random_state_rain"])
    )  # Intensity of rain events
    r = np.random.RandomState(int(param_dict["real_random"]))  # General seed
    r_async = np.random.RandomState(
        int(param_dict["rho_seed"])
    )  # Hydroclimatic asynchrony

    # Extract biogeochemical parameters
    gamma = param_dict[
        "gamma_sc"
    ]  # Structured heterogeneity/correlation of solute sources to travel time
    sigma_w = param_dict["sigma_w_sc"]  # Random heterogeneity
    damkohler_interarrival = param_dict["damkohler_interarrival"]  # Da1
    damkohler_transport = param_dict["damkohler_transport"]  # Da2
    damkohler_longterm = param_dict["damkohler_longterm"]  # Da3

    oro_scaling = param_dict["oro_scaling"]

    # Simulation time setup
    spin_up_days = SPIN_UP_YEARS * DAYS_PER_YEAR  # [days] spin-up time
    simulation_years = int(param_dict["tmax_yrs"])  # simulation time [years]
    total_simulation_days = (simulation_years * DAYS_PER_YEAR) + spin_up_days

    # Climate parameters
    rain_per_year = param_dict["rain_per_year"]  # [cm/yr] mean rainfall amount
    interarrival_time_mean = param_dict[
        "interarrival_time_mean"
    ]  # [d] mean rainfall interarrival time
    lambda_time = (
        1 / interarrival_time_mean
    )  # [1/days] mean of the rainfall arrival time exponential PDF
    rainfall_mean = rain_per_year / 365 / lambda_time  # [cm/day]
    lambda_intensity = (
        1 / rainfall_mean
    )  # [1/cm] mean of the rainfall depth exponential PDF
    ET_max = param_dict["ET_max"]  # [cm/day] potential evapotranspiration

    # Soil and hydrology parameters
    z_r = param_dict["z_r"]  # [cm] root zone thickness
    z_vz = param_dict[
        "z_vz"
    ]  # [cm] thickness of the entire vadose zone including root zone
    theta_fc = param_dict["theta_fc"]  # [-] water content at field capacity
    theta_wp = param_dict["theta_wp"]  # [-] water content at permanent wilting point
    theta_res = param_dict["theta_res"]  # [-] residual water content
    soil_storage = z_r * (theta_fc - theta_wp)  # [cm] soil water storage root zone
    soil_storage_vz = z_vz * (
        theta_fc - theta_wp
    )  # [cm] soil water storage entire vadose zone
    R = param_dict["R"]  # Retardation factor (Harman et al. 2011)
    mean_tr = param_dict["mean_tr"]  # [days] hydraulic response time

    theta_sat = param_dict["theta_sat"]  # [-] saturated water content aquifer
    aquifer_z = param_dict["aquifer_z"]  # [cm] aquifer thickness

    # Dimensionless ratios
    gamma_harman = (
        soil_storage / rainfall_mean
    )  # Harman et al. [2011] ratio for root zone
    gamma_harman_vz = (
        soil_storage_vz / rainfall_mean
    )  # Harman et al. [2011] ratio for entire vadose zone
    bigF_harman = (theta_fc - theta_res) / (
        theta_fc - theta_wp
    )  # Harman et al. 2011 parameter
    phi = ET_max / rainfall_mean / lambda_time  # Harman et al. [2011] dryness index

    # [1/days] mean of the drainage arrival time exponential PDF root zone
    lambda_d = (
        lambda_time
        / (
            scipy.special.gamma(gamma_harman / phi)
            - scipy.stats.gamma.sf(x=gamma_harman, a=gamma_harman / phi)
            * scipy.special.gamma(gamma_harman / phi)
        )
        * np.exp(-gamma_harman)
        * gamma_harman ** (gamma_harman / phi)
        * phi
        / gamma_harman
    )

    # [1/days] mean of the drainage arrival time exponential PDF entire vadose zone
    lambda_d_vz = (
        lambda_time
        / (
            scipy.special.gamma(gamma_harman_vz / phi)
            - scipy.stats.gamma.sf(x=gamma_harman_vz, a=gamma_harman_vz / phi)
            * scipy.special.gamma(gamma_harman_vz / phi)
        )
        * np.exp(-gamma_harman_vz)
        * gamma_harman_vz ** (gamma_harman_vz / phi)
        * phi
        / gamma_harman_vz
    )

    # "Effective gamma" from Harman et al. [2011] entire vadose zone
    lambda_e_vz = lambda_d_vz
    gamma_harman_e_vz = gamma_harman_vz * (lambda_d_vz / lambda_time) ** 0.5

    # Travel time parameters vadose zone (=waiting time)
    # These use "effective gamma" from Harman et al. [2011]
    mean_waiting = (
        R * bigF_harman * gamma_harman_e_vz / lambda_e_vz
    )  # Equation (6), mean travel time vadose zone

    # Linear reservoir recession constant, tr (k)
    lambda_tr = 1 / mean_tr  # [1/day]

    # SOLUTE RESPONSE - TRAVEL TIME
    # Travel time as sum of vadose zone travel (waiting) time and groundwater travel time
    # lambda_intensity does not change from vadose filtering
    # total drainage does change = lambda_d/lambda_intensity

    drainage_mean = lambda_d / lambda_intensity  # [cm/day] mean drainage from root zone
    mean_g = aquifer_z * theta_sat / drainage_mean  # [days], equation (7)

    # Vadose and groundwater travel times: exponential
    lambda_waiting = 1 / mean_waiting
    lambda_g = 1 / mean_g

    # m1 and m2 based on lambda time vadose zone and lambda time groundwater
    m1_tp = 1 / lambda_g + 1 / lambda_waiting  # Equation (9)
    m2_tp = (
        2 / lambda_g**2 + 2 / lambda_waiting**2 + 2 / (lambda_g * lambda_waiting)
    )  # Equation (10)

    # Hyperexponential fit of travel time with LN
    mu_ln_time = 2 * np.log(m1_tp) - np.log(m2_tp) / 2
    sigma_ln_time = (np.log(m2_tp) - 2 * np.log(m1_tp)) ** 0.5
    mean_travel_time = np.exp(
        mu_ln_time + sigma_ln_time**2 / 2
    )  # Mean total travel time

    # Times mean drainage rate, "mean log(streamtube length)": needed for
    # redistribution of mobilized masses to following events
    mu_ln_tube = mu_ln_time + np.log(drainage_mean)
    sigma_ln_tube = sigma_ln_time

    # SOURCE CONCENTRATION PARAMETERS
    # Immobile concentration correlation to travel time - random heterogeneity
    c_c = 1  # See equation (14), 1 for simplicity
    mu_w = -sigma_w**2 / 2  # Unit-mean log-normally distributed random variable (equation 14)

    # Mean immobile concentration
    mean_c_im_number = param_dict["mean_c_im_number"]  # [mg/L]
    mean_c_im = np.zeros(
        (total_simulation_days, 1)
    )  # Matrix to account for temporal variance
    mean_c_im[:] = mean_c_im_number

    # Damkohler 3: logistic long-term change of source concentration
    max_c_im = 10 * mean_c_im_number
    longterm_change_time = mean_travel_time / damkohler_longterm
    logistic_shape = 2 / longterm_change_time
    mean_c_im = np.array(
        [
            max_c_im
            / (
                1
                + np.exp(
                    -logistic_shape
                    * (np.arange(1, total_simulation_days + 1, 1) - longterm_change_time)
                )
            )
        ]
    )
    mean_c_im = mean_c_im.T

    # Rate constants for reactions
    # Damkohler 1: drainage event inter-arrival, initial mobile concentration
    rate_interarrival = damkohler_interarrival * lambda_d  # Equation (19)

    # Damkohler 2: mobile concentration (low value = little degradation)
    rate_transport = damkohler_transport / mean_travel_time  # Equation (21)

    # RAINFALL EVENT GENERATION
    # Generate sufficient events to cover simulation period
    # Oversized to avoid time-consuming loops [approach from Gavan McGrath]
    N = int(100 * total_simulation_days * lambda_d)

    interval_l = r_nonzero_times.exponential(1 / lambda_d, int(N)).round(0).astype(int)
    nonzero_rain = r_nonzero_rain.exponential(
        1 / lambda_intensity, int(N)
    )  # lambda_intensity does not change from filtering (Botter 2007)

    # Apply binomial thinning based on asynchrony parameter
    rho = param_dict["rho"]
    keep = r_async.rand(int(N)) < rho
    keep_index = np.where(keep)[0]
    new_index = np.where(~keep)[0]

    # Create 2D array for hydrologic events
    df_hydro = np.zeros((int(N), 6))
    df_hydro[:, 0] = interval_l

    interval_l_new = r_async.exponential(1 / lambda_d, int(N)).round(0).astype(int)
    df_hydro[:, 1] = interval_l_new
    df_hydro[:, 2] = nonzero_rain

    nonzero_rain_new = r_async.exponential(1 / lambda_intensity, int(N))
    df_hydro[:, 3] = nonzero_rain_new

    df_hydro[:, 0] = np.cumsum(df_hydro[:, 0])
    df_hydro[:, 1] = np.cumsum(df_hydro[:, 1])

    df_hydro[keep_index, 4] = df_hydro[keep_index, 0]
    df_hydro[new_index, 4] = df_hydro[new_index, 1]

    df_hydro[keep_index, 5] = df_hydro[keep_index, 2]
    df_hydro[new_index, 5] = df_hydro[new_index, 3]

    # Sort df_hydro based on column 4
    df_hydro = df_hydro[np.argsort(df_hydro[:, 4])]

    all_times = df_hydro[:, 4]
    all_times = np.unique(all_times)  # Remove duplicate events on same day
    all_times = all_times[1:]  # Remove first day as this can be zero

    number_events = len(all_times[all_times < total_simulation_days])

    # Cut vector to correct length and calculate intervals between events
    interval = np.diff(all_times)  # forward-looking
    interval = np.insert(interval, 0, all_times[0])  # prepend gap before first event
    interval = interval[0:number_events].astype(int)

    nonzero_times = all_times[0:number_events].astype(int)
    nonzero_rain = df_hydro[0:number_events, 5]

    # Scale rain intensity based on orographic scaling factor
    elevation = param_dict["elevation"]
    nonzero_rain = nonzero_rain + (elevation / 1000 * oro_scaling) * nonzero_rain

    # Initialize output arrays
    discharge = np.zeros((number_events, 1))
    mass_flux = np.zeros((number_events, 1))
    event_c_m = np.zeros((number_events, 1))  # Immobile and flux-avg concentrations

    # vector of running sums of nonzero rain (effective)
    rain_sum = np.cumsum(nonzero_rain)

    # DISCHARGE CALCULATION
    # Initialize matrices for event-based discharge and flux tracking
    event_discharge = np.zeros((total_simulation_days, number_events))
    event_flux = np.zeros((total_simulation_days, number_events))

    # Pre-calculate unit hydrograph (exponential recession)
    unit_hydrograph = np.asarray(expon.pdf(
        np.arange(0, total_simulation_days), scale=1 / lambda_tr
    ))

    for events_c in np.arange(0, number_events):
        # Calculate discharge hydrograph for each event
        event_start_time = int(nonzero_times[events_c])

        if event_start_time == total_simulation_days:
            continue
        else:
            event_discharge[event_start_time:total_simulation_days, events_c] = (
                float(nonzero_rain[events_c])
                * unit_hydrograph[1 : total_simulation_days - event_start_time + 1]
            )

    # =============================================================================
    # SOLUTE FLUX
    # =============================================================================
    # Distribute solute among events
    # Initialize with zero values so can add to each other

    # Event initial mobile concentration after immobile-mobile exchange equation (18)
    event_c_m = mean_c_im[nonzero_times - 1, 0] * (
        1 - np.exp(-rate_interarrival * interval)
    )

    # Mean mobile concentration scaling parameter
    a_c = event_c_m * np.exp(
        -0.5 * (gamma**2 * sigma_ln_time**2 + c_c**2 * sigma_w**2)
        - gamma * mu_ln_time
        - c_c * mu_w
    ) # equation (17)

    event_timers = nonzero_times + 1

    if event_timers[-1] == total_simulation_days:
        event_timers[-1] = nonzero_times[-1]

    tau_total = np.array(np.arange(1, total_simulation_days + 1))

    # Calculate cumulative rain for each event
    rain_sums = [
        calculate_cumulative_rainfall_for_event(event_index=xi, cumulative_rainfall=rain_sum, max_time=total_simulation_days)
        for xi in np.arange(0, number_events)
    ]

    # Calculate event fractions for future events for each event
    event_fractions_list = [
        calculate_event_mass_fractions(cumulative_rainfall=xi, mu_ln_streamtube=mu_ln_tube, sigma_ln_streamtube=sigma_ln_tube)
        for xi in rain_sums
    ]

    # Create scaling matrix; corrected from original version
    df_scaling = pd.DataFrame(
        np.zeros((number_events, number_events)),
        index=np.arange(0, number_events),
        columns=np.arange(0, number_events),
    )

    for events_c in np.arange(0, number_events):
        event_fractions_shifted = event_fractions_list[events_c]

        k_count = 0
        for k in np.arange(events_c, number_events):
            df_scaling.loc[events_c, k] = event_fractions_shifted[k_count]
            k_count += 1

    df_scaling_summed = df_scaling.sum(axis=0)

    # Iterate over events
    for events_c in np.arange(0, number_events):
        # Track allocation of "events_c" to later events
        k_count = 0
        event_fractions_shifted = event_fractions_list[events_c]

        for k in np.arange(events_c, number_events):
            # From time of event k to end of simulation
            tau = tau_total[0 : (total_simulation_days + 1) - event_timers[k]]

            random_lognorm = r.lognormal(mean=mu_w, sigma=sigma_w, size=len(tau))
            c_m_tau_vector = (
                a_c[events_c]  # Scaling parameter
                * tau**gamma  # Structured heterogeneity
                * random_lognorm**c_c
            )  # Unstructured heterogeneity (equation 14)

            # MOBILE CONCENTRATION
            # Exponential decay along streamtube
            c_m_decay_tau_vector = c_m_tau_vector * np.exp(
                -rate_transport * tau
            )  # Equation (20)

            # Add solute flux to the matrix
            event_flux[event_timers[k] - 1 : total_simulation_days, events_c] = (
                event_flux[event_timers[k] - 1 : total_simulation_days, events_c]
                + c_m_decay_tau_vector
                * event_discharge[event_timers[k] - 1 : total_simulation_days, k]
                * event_fractions_shifted[k_count]
                / df_scaling_summed[k]
            )

            k_count += 1

    # OUTPUT DISCHARGE, LOAD AND CONCENTRATION [daily]
    # Sum of rows
    discharge = event_discharge.sum(axis=1)
    mass_flux = event_flux.sum(axis=1)

    # Get number of discharge == 0 from initiation to avoid divide by 0 warning
    len_init = len(discharge[discharge == 0])

    # Set discharge to -9999 if discharge == 0
    discharge[discharge == 0] = -9999
    concentration = mass_flux / discharge

    # Set concentration to np.nan for the initiation phase
    concentration[0:len_init] = np.nan
    discharge[0:len_init] = np.nan

    # Remove spin-up period from output time series
    discharge = discharge[spin_up_days:]
    concentration = concentration[spin_up_days:]

    return discharge, concentration, mass_flux


# ============================================================================
# ROUTING FUNCTIONS
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
        Reference time step for the simulation [s].
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
        Velocity factor or parameter (usage depends on implementation).
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
        Time series of time steps [s].
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
            # 1 / 60,
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
    sinus_alpha = 1 / (1 + side_slope**2) ** 0.5
    cot_alpha = side_slope / 1
    safety_factor = 0.9

    # Set initial dt to reference dt
    dt = dt_ref

    # Initialize counters
    i = 0.0
    counter = 0

    # Pre-allocate arrays (oversized, will be trimmed later)
    len_empty_array = 10 * 365 * 24 * 60  # 10 years in minutes

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

    # Equation C1: A(y) - wetted area
    wetted_arearef1 = (bottom_width_int + stageref1 * cot_alpha) * stageref1

    # Calculate initial storage load based on wetted area
    load_storage_t_minus_1 = wetted_arearef1 * dx * conc[0]
    if not optimized_calc:
        storage_list[0] = wetted_arearef1 * dx
        load_storage_list[0] = storage_list[0] * conc[0]

    # Equation C2: B(y) - surface width
    surface_widthref1 = bottom_width_int + (2 * stageref1 * cot_alpha)

    # Equation C3: P(y) - wetted perimeter
    wetted_perimeterref1 = bottom_width_int + ((2 * stageref1) / sinus_alpha)

    # Equation C5: velocity
    vref1 = (bottom_slope**0.5 / mannings_n) * (
        (wetted_arearef1 ** (2 / 3)) / wetted_perimeterref1 ** (2 / 3)
    )

    # Equation C6: celerity
    cref1 = (
        (5 / 3)
        * (bottom_slope**0.5 / mannings_n)
        * (wetted_arearef1 ** (2 / 3) / wetted_perimeterref1 ** (2 / 3))
        * (
            1
            - (4 / 5)
            * (
                wetted_arearef1
                / (surface_widthref1 * wetted_perimeterref1 * sinus_alpha)
            )
        )
    )

    if not optimized_calc:
        celerity_list[0] = cref1

    # Equation C7: equivalent to Equation 49a
    betaref1 = cref1 / vref1

    # Equation 51a: Reynolds number
    reynoldref1 = qref1 / (betaref1 * surface_widthref1 * bottom_slope * cref1 * dx)

    # Calculate optimal dt based on dx and wave celerity
    initial_dt = safety_factor * ((dx / cref1) / 3600)

    # Select next smallest dt from dt_array
    initial_dt = max(min(dt_array), initial_dt)
    initial_dt = max([i for i in dt_array if i <= initial_dt])
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

        # Concentration routing and mixing routines
        volume_in_t_1 = in_t_1 * dt_in_sec
        volume_out_t_1 = out_t_1 * dt_in_sec
        conc_storage_t = conc_routed_list[counter]

        # Mass balance calculations
        load_storage_t = load_storage_t_minus_1  # Load in storage at time t
        load_storage_t_1 = (
            load_storage_t
            - (volume_out_t_1 * conc_storage_t)
            + (volume_in_t_1 * conc_in_t_1)
        )  # Load in storage at time t+1

        # Check for negative load (numerical instability)
        # Occurs when concentration decreases rapidly - outgoing > incoming + stored
        if load_storage_t_1 < 0:
            if dt != min(dt_array):
                # Select smaller dt and recalculate
                dt_next = [dt_next for dt_next in dt_array if dt_next < dt][0]
                continue

        # In-stream 1st order loss (settling/decay)
        ke = vf / stageref2 / 86400  # Convert vf to loss rate [1/s]
        load_storage_t_1_loss = load_storage_t_1 * np.exp(-ke * dt_in_sec)

        # Calculate concentration in storage at t+1
        conc_storage_t_1 = load_storage_t_1_loss / storage_t_1

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
    
    # Trim arrays to actual length (remove pre-allocated zeros)
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
    flag = True

    # Equation 45: initial guess estimate
    q_guess = out_t + (in_t_1 - in_t)
    last_guess = q_guess * 2  # Ensure at least one iteration

    # Convergence threshold and safeguards
    tresh = 0.003
    if abs(last_guess - q_guess) <= tresh:
        last_guess = last_guess + abs(last_guess - q_guess) + tresh

    percent = abs(100 - (q_guess / last_guess) * 100)

    # Iterative solution loop
    while (abs(last_guess - q_guess) > tresh) and flag:
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
    
            # Equation C1: A(y) - wetted area
            wetted_arearef2 = (bottom_width_int + stageref2 * cot_alpha) * stageref2
    
            # Equation C2: B(y) - surface width
            surface_widthref2 = bottom_width_int + (2 * stageref2 * cot_alpha)
    
            # Equation C3: P(y) - wetted perimeter
            wetted_perimeterref2 = bottom_width_int + ((2 * stageref2) / sinus_alpha)
    
            # Equation C5: velocity
            vref2 = (bottom_slope**0.5 / mannings_n) * (
                (wetted_arearef2 ** (2 / 3)) / wetted_perimeterref2 ** (2 / 3)
            )
    
            # Equation C6: celerity
            cref2 = (
                (5 / 3)
                * (bottom_slope**0.5 / mannings_n)
                * (wetted_arearef2 ** (2 / 3) / wetted_perimeterref2 ** (2 / 3))
                * (
                    1
                    - (4 / 5)
                    * (
                        wetted_arearef2
                        / (surface_widthref2 * wetted_perimeterref2 * sinus_alpha)
                    )
                )
            )
    
            # Equation C7: equivalent to Equation 49b
            betaref2 = cref2 / vref2
            
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
        flag = True

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
            flag = False
        elif courantref2 > COURANT_NUMBER_THRESHOLD_HIGH:  # Courant too large, decrease time step
            dt = safety_factor * ((dx / cref2) / 3600)
            flag = False

        # Select appropriate dt from available options
        if not flag:
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

# ============================================================================
# NETWORK ROUTING FUNCTIONS
# ============================================================================


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
        mass transfer coefficient

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
        * (bottom_slope**0.5)
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
        * (bottom_slope**0.5 / mannings_n)
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
    based on Equation B3 in Todini et al., (2007).

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
