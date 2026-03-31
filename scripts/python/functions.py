# -*- coding: utf-8 -*-
"""
Stochastic Landscape Model (SLM) - Core Functions

Author: schauer
Created: Fri Feb 23 09:40:07 2024
"""

from stochastic_headwater_model import (
    stochastic_headwater_model,
    calculate_cumulative_rainfall_for_event,
    calculate_event_mass_fractions,
    SPIN_UP_YEARS,
    DAYS_PER_YEAR,
)

from routing import (
    routing_function_ocn,
    reach_routing,
    reach_routing_single_time_step,
    resample_routed_lists,
    compute_hydraulic_properties,
    stage_discharge_relationship,
    d_stage_discharge_relationship,
    newton_raphson,
    MINIMUM_DISCHARGE,
    COURANT_NUMBER_THRESHOLD_HIGH,
    COURANT_NUMBER_THRESHOLD_LOW,
    MAX_ITER,
    SECONDS_PER_DAY,
    INFLOW_RATIO_THRESHOLD,
)