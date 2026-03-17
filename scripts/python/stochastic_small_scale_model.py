# -*- coding: utf-8 -*-
"""
Stochastic Landscape Model (SLM) - Subcatchment-Scale Stochastic Model

Implements the stochastic model of Musolff et al. (2017)
(https://doi.org/10.1002/2017GL072630) for simulating coupled hydrologic
and biogeochemical responses in catchments.

Author: schauer
Created: Fri Feb 23 09:40:07 2024
"""

from typing import Dict, Tuple, Union

import numpy as np
import pandas as pd
import scipy.special
import scipy.stats
from scipy.stats import expon, lognorm

# Constants
SPIN_UP_YEARS = 7  # Default spin-up time for subcatchment_scale_module in years
DAYS_PER_YEAR = 365  # Days per year


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
# MAIN SIMULATION FUNCTION
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
    # Original equation (17) included gamma-dependent terms:
    #   a_c = event_c_m * np.exp(
    #       -0.5 * (gamma**2 * sigma_ln_time**2 + c_c**2 * sigma_w**2)
    #       - gamma * mu_ln_time
    #       - c_c * mu_w
    #   )
    # These terms are no longer needed because tau^gamma is normalized to have
    # discharge-weighted mean of 1 within each recession (see weighted_tau_gamma_mean
    # below), so a_c only needs to correct for the unstructured heterogeneity W.
    a_c = event_c_m * np.exp(
        -0.5 * c_c**2 * sigma_w**2
        - c_c * mu_w
    )

    event_timers = nonzero_times + 1

    if event_timers[-1] == total_simulation_days:
        event_timers[-1] = nonzero_times[-1]

    tau_total = np.array(np.arange(1, total_simulation_days + 1))
    
    # --- Correction for tau^gamma scale mismatch ---
    # The scaling parameter a_c (equation 17) is calibrated for tau drawn from the
    # lognormal travel time distribution (mean ~hundreds to thousands of days).
    # But in the inner loop, tau^gamma is evaluated at recession indices 1, 2, 3, ...
    # weighted by exponential discharge decay. Without correction, negative gamma
    # produces systematically higher concentrations than positive gamma because
    # a_c overcompensates at short recession times.
    # Fix: normalize tau^gamma within each recession to have discharge-weighted
    # mean of 1, so that tau^gamma only controls the *shape* of within-event
    # concentration (enrichment vs dilution) without affecting the *level*.
    h_weights = expon.pdf(tau_total, scale=1 / lambda_tr)
    cumulative_h = np.cumsum(h_weights)
    cumulative_tau_gamma_h = np.cumsum(tau_total**gamma * h_weights)
    weighted_tau_gamma_mean = cumulative_tau_gamma_h / cumulative_h

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
                a_c[events_c]
                * tau**gamma
                / weighted_tau_gamma_mean[len(tau) - 1]
                * random_lognorm**c_c
            )

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
