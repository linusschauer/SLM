# -*- coding: utf-8 -*-
"""
Stochastic Landscape Model (SLM) - Core Functions

Backward-compatible re-export module. Import from stochastic_model or routing
directly for new code.

This script contains:
(1) Subcatchment-scale stochastic model of Musolff et al. (2017) (https://doi.org/10.1002/2017GL072630)
(2) MCT routing algorithm (https://doi.org/10.5194/hess-11-1645-2007)

Author: schauer
Created: Fri Feb 23 09:40:07 2024
"""

from stochastic_small_scale_model import (
    subcatchment_scale_module,
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
    stage_discharge_relationship,
    d_stage_discharge_relationship,
    newton_raphson,
    MINIMUM_DISCHARGE,
    COURANT_NUMBER_THRESHOLD_HIGH,
    COURANT_NUMBER_THRESHOLD_LOW,
    MAX_ITER,
)
