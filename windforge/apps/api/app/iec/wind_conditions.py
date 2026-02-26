"""
IEC 61400-1 Ed.4 wind condition calculations.

Implements all standard wind condition formulas from IEC 61400-1:2019
including turbulence models, extreme wind speeds, and deterministic
gust/direction-change amplitudes.

Reference values for turbulence categories:
  Class A: I_ref = 0.16
  Class B: I_ref = 0.14
  Class C: I_ref = 0.12

Reference wind speeds by turbine class:
  Class I:   V_ref = 50 m/s
  Class II:  V_ref = 42.5 m/s
  Class III: V_ref = 37.5 m/s
  Class S:   site-specific

References
----------
IEC 61400-1:2019 (Ed.4), Sections 6.3 and 6.4
"""

from __future__ import annotations

import math
from enum import Enum
from typing import Optional


# ---------------------------------------------------------------------------
# Constants from IEC 61400-1 Ed.4
# ---------------------------------------------------------------------------

# Turbulence intensity reference values by category (Table 1)
IREF_VALUES: dict[str, float] = {
    "A": 0.16,
    "B": 0.14,
    "C": 0.12,
}

# Reference wind speeds by turbine class (Table 1)
VREF_VALUES: dict[int, float] = {
    1: 50.0,    # Class I
    2: 42.5,    # Class II
    3: 37.5,    # Class III
}

# IEC standard Lambda_1 (turbulence scale parameter)
# Lambda_1 = 0.7 * min(hub_height, 60) for hub_ht <= 60m
# Lambda_1 = 42 for hub_ht > 60m
LAMBDA1_CUTOFF_HEIGHT = 60.0  # meters


def get_iref(turbulence_class: str) -> float:
    """Get the reference turbulence intensity for a given class.

    Parameters
    ----------
    turbulence_class : str
        IEC turbulence category: "A", "B", or "C".

    Returns
    -------
    float
        Reference turbulence intensity I_ref.

    Raises
    ------
    ValueError
        If turbulence_class is not "A", "B", or "C".
    """
    tc = turbulence_class.upper()
    if tc not in IREF_VALUES:
        raise ValueError(
            f"Invalid turbulence class '{turbulence_class}'. "
            f"Must be one of: A, B, C"
        )
    return IREF_VALUES[tc]


def get_vref(turbine_class: int) -> float:
    """Get the reference wind speed for a given turbine class.

    Parameters
    ----------
    turbine_class : int
        IEC wind turbine class: 1, 2, or 3.

    Returns
    -------
    float
        Reference wind speed V_ref (m/s).

    Raises
    ------
    ValueError
        If turbine_class is not 1, 2, or 3.
    """
    if turbine_class not in VREF_VALUES:
        raise ValueError(
            f"Invalid turbine class {turbine_class}. Must be 1, 2, or 3."
        )
    return VREF_VALUES[turbine_class]


def get_lambda1(hub_height: float) -> float:
    """Calculate the longitudinal turbulence scale parameter Lambda_1.

    Per IEC 61400-1 Ed.4 Eq. 5:
      Lambda_1 = 0.7 * min(z_hub, 60)

    Parameters
    ----------
    hub_height : float
        Hub height above ground (m).

    Returns
    -------
    float
        Turbulence scale parameter Lambda_1 (m).
    """
    return 0.7 * min(hub_height, LAMBDA1_CUTOFF_HEIGHT)


# ---------------------------------------------------------------------------
# Normal Turbulence Model (NTM) -- Eq. 11
# ---------------------------------------------------------------------------

def ntm_sigma(v_hub: float, i_ref: float) -> float:
    """Calculate the NTM representative turbulence standard deviation.

    IEC 61400-1 Ed.4, Eq. 11:
      sigma_1 = I_ref * (0.75 * V_hub + 5.6)

    This gives the representative value of the hub-height longitudinal
    wind velocity standard deviation.

    Parameters
    ----------
    v_hub : float
        Hub-height mean wind speed (m/s).
    i_ref : float
        Reference turbulence intensity (-).

    Returns
    -------
    float
        Representative longitudinal turbulence standard deviation sigma_1 (m/s).
    """
    return i_ref * (0.75 * v_hub + 5.6)


# ---------------------------------------------------------------------------
# Extreme Turbulence Model (ETM) -- Eq. 19
# ---------------------------------------------------------------------------

def etm_sigma(v_hub: float, i_ref: float, v_ref: float) -> float:
    """Calculate the ETM representative turbulence standard deviation.

    IEC 61400-1 Ed.4, Eq. 19:
      c = 2.0  (default)
      sigma_1 = c * I_ref * (0.072 * (V_ave / c + 3) * (V_hub / c - 4) + 10)
    where V_ave = 0.2 * V_ref

    Simplified per Ed.4:
      sigma_1 = I_ref * (c1 * V_hub / V_ave + c2)
    with appropriate constants.

    Parameters
    ----------
    v_hub : float
        Hub-height mean wind speed (m/s).
    i_ref : float
        Reference turbulence intensity (-).
    v_ref : float
        Reference wind speed for the turbine class (m/s).

    Returns
    -------
    float
        ETM longitudinal turbulence standard deviation sigma_1 (m/s).
    """
    c = 2.0
    v_ave = 0.2 * v_ref
    sigma_1 = c * i_ref * (
        0.072 * (v_ave / c + 3.0) * (v_hub / c - 4.0) + 10.0
    )
    # Ensure non-negative
    return max(sigma_1, 0.0)


# ---------------------------------------------------------------------------
# Extreme Wind Model (EWM) -- Eqs. 20-21
# ---------------------------------------------------------------------------

def ewm_ve50(v_ref: float) -> float:
    """Calculate the 50-year extreme wind speed at hub height.

    IEC 61400-1 Ed.4, Eq. 20:
      V_e50 = 1.4 * V_ref

    Parameters
    ----------
    v_ref : float
        Reference wind speed for the turbine class (m/s).

    Returns
    -------
    float
        50-year recurrence extreme wind speed V_e50 (m/s).
    """
    return 1.4 * v_ref


def ewm_ve1(v_ref: float) -> float:
    """Calculate the 1-year extreme wind speed at hub height.

    IEC 61400-1 Ed.4, Eq. 21:
      V_e1 = 0.8 * V_e50 = 0.8 * 1.4 * V_ref = 1.12 * V_ref

    Parameters
    ----------
    v_ref : float
        Reference wind speed for the turbine class (m/s).

    Returns
    -------
    float
        1-year recurrence extreme wind speed V_e1 (m/s).
    """
    return 0.8 * ewm_ve50(v_ref)


# ---------------------------------------------------------------------------
# Extreme Operating Gust (EOG) -- Eq. 15
# ---------------------------------------------------------------------------

def eog_amplitude(
    v_hub: float,
    v_ref: float,
    i_ref: float,
    rotor_diameter: float,
) -> float:
    """Calculate the EOG gust amplitude.

    IEC 61400-1 Ed.4, Eq. 15:
      V_gust = min(
        1.35 * (V_e1 - V_hub),
        3.3 * sigma_1 / (1 + 0.1 * D / Lambda_1)
      )

    where sigma_1 is from the NTM model and Lambda_1 is the
    turbulence scale parameter.

    Parameters
    ----------
    v_hub : float
        Hub-height mean wind speed (m/s).
    v_ref : float
        Reference wind speed (m/s).
    i_ref : float
        Reference turbulence intensity (-).
    rotor_diameter : float
        Rotor diameter (m).

    Returns
    -------
    float
        EOG gust amplitude V_gust (m/s).
    """
    v_e1 = ewm_ve1(v_ref)
    sigma_1 = ntm_sigma(v_hub, i_ref)

    # Use hub height ~ rotor center; approximate as D/2 + tower base
    # For the Lambda_1 calculation we need hub height, approximate it
    hub_ht_approx = rotor_diameter  # reasonable approximation
    lambda_1 = get_lambda1(hub_ht_approx)

    term1 = 1.35 * (v_e1 - v_hub)
    term2 = 3.3 * sigma_1 / (1.0 + 0.1 * rotor_diameter / lambda_1)

    return min(term1, term2)


# ---------------------------------------------------------------------------
# Extreme Direction Change (EDC) -- Eq. 16
# ---------------------------------------------------------------------------

def edc_amplitude(
    v_hub: float,
    v_ref: float,
    i_ref: float,
    rotor_diameter: float,
) -> float:
    """Calculate the EDC extreme direction change magnitude.

    IEC 61400-1 Ed.4, Eq. 16:
      theta_e = +/- arctan(
        sigma_1 / (V_hub * (1 + 0.1 * D / Lambda_1))
      )

    The result is limited to +/- 180 degrees.

    Parameters
    ----------
    v_hub : float
        Hub-height mean wind speed (m/s).
    v_ref : float
        Reference wind speed (m/s).
    i_ref : float
        Reference turbulence intensity (-).
    rotor_diameter : float
        Rotor diameter (m).

    Returns
    -------
    float
        EDC direction change magnitude theta_e (degrees).
    """
    sigma_1 = ntm_sigma(v_hub, i_ref)

    hub_ht_approx = rotor_diameter
    lambda_1 = get_lambda1(hub_ht_approx)

    denominator = v_hub * (1.0 + 0.1 * rotor_diameter / lambda_1)
    if denominator <= 0.0:
        return 180.0

    theta_e_rad = math.atan(sigma_1 / denominator)
    theta_e_deg = math.degrees(theta_e_rad)

    # Limited to +/- 180 degrees
    return min(theta_e_deg, 180.0)


# ---------------------------------------------------------------------------
# Extreme Coherent gust with Direction change (ECD) -- Eq. 14
# ---------------------------------------------------------------------------

def ecd_amplitude(v_hub: float) -> float:
    """Calculate the ECD coherent gust amplitude.

    IEC 61400-1 Ed.4, Eq. 14:
      V_cg = 15 m/s (constant coherent gust amplitude)

    The direction change component theta_cg:
      theta_cg = 180 deg              for V_hub < 4 m/s
      theta_cg = 720 deg / V_hub      for 4 <= V_hub <= V_ref

    Parameters
    ----------
    v_hub : float
        Hub-height mean wind speed (m/s).

    Returns
    -------
    float
        ECD coherent gust magnitude V_cg (m/s). The constant value
        is 15 m/s per the standard.
    """
    return 15.0


def ecd_direction_change(v_hub: float) -> float:
    """Calculate the ECD direction change component.

    Parameters
    ----------
    v_hub : float
        Hub-height mean wind speed (m/s).

    Returns
    -------
    float
        Direction change theta_cg (degrees).
    """
    if v_hub < 4.0:
        return 180.0
    else:
        return min(720.0 / v_hub, 180.0)


# ---------------------------------------------------------------------------
# Extreme Wind Shear (EWS) -- Eq. 17-18
# ---------------------------------------------------------------------------

def ews_amplitude(
    v_hub: float,
    i_ref: float,
    rotor_diameter: float,
) -> float:
    """Calculate the EWS extreme wind shear magnitude.

    IEC 61400-1 Ed.4, Eqs. 17-18. The EWS consists of both vertical
    and horizontal transient shear events. This function returns
    the maximum shear parameter.

    Vertical shear:
      V(z,t) = V_hub * (z/z_hub)^alpha +/- (
        (2.5 + 0.2 * beta * sigma_1 * (D / Lambda_1)^0.25)
        * (z - z_hub) / D
      ) * (1 - cos(2*pi*t/T))

    where beta = 6.4 (vertical) or 6.4 (horizontal),
    and T = 12 s (transient period).

    Parameters
    ----------
    v_hub : float
        Hub-height mean wind speed (m/s).
    i_ref : float
        Reference turbulence intensity (-).
    rotor_diameter : float
        Rotor diameter (m).

    Returns
    -------
    float
        EWS shear amplitude parameter (m/s). This is the peak
        additional velocity at the blade tip due to the shear event.
    """
    sigma_1 = ntm_sigma(v_hub, i_ref)

    hub_ht_approx = rotor_diameter
    lambda_1 = get_lambda1(hub_ht_approx)

    beta = 6.4
    radius = rotor_diameter / 2.0

    # Shear velocity at blade tip (z - z_hub = R)
    shear_amp = (
        2.5 + 0.2 * beta * sigma_1 * (rotor_diameter / lambda_1) ** 0.25
    ) * (radius / rotor_diameter)

    # The factor of 2 accounts for the (1 - cos) term peak value
    return shear_amp * 2.0


# ---------------------------------------------------------------------------
# Weibull distribution parameters
# ---------------------------------------------------------------------------

def weibull_probability(
    v: float,
    v_ave: float,
    k: float = 2.0,
) -> float:
    """Calculate the Weibull probability density for a wind speed.

    Parameters
    ----------
    v : float
        Wind speed (m/s).
    v_ave : float
        Annual average wind speed = 0.2 * V_ref (m/s).
    k : float
        Weibull shape parameter (default 2.0 = Rayleigh).

    Returns
    -------
    float
        Probability density f(v).
    """
    c = v_ave / math.gamma(1.0 + 1.0 / k)  # Scale parameter
    if c <= 0 or v < 0:
        return 0.0
    return (k / c) * (v / c) ** (k - 1) * math.exp(-(v / c) ** k)


def weibull_cdf(
    v: float,
    v_ave: float,
    k: float = 2.0,
) -> float:
    """Calculate the Weibull CDF for a wind speed.

    Parameters
    ----------
    v : float
        Wind speed (m/s).
    v_ave : float
        Annual average wind speed (m/s).
    k : float
        Weibull shape parameter (default 2.0).

    Returns
    -------
    float
        Cumulative probability F(v).
    """
    c = v_ave / math.gamma(1.0 + 1.0 / k)
    if c <= 0 or v < 0:
        return 0.0
    return 1.0 - math.exp(-(v / c) ** k)


def wind_speed_bin_probability(
    v_low: float,
    v_high: float,
    v_ave: float,
    k: float = 2.0,
) -> float:
    """Calculate the probability of wind speed falling in a bin.

    Parameters
    ----------
    v_low : float
        Lower edge of wind speed bin (m/s).
    v_high : float
        Upper edge of wind speed bin (m/s).
    v_ave : float
        Annual average wind speed (m/s).
    k : float
        Weibull shape parameter.

    Returns
    -------
    float
        Probability of wind speed in [v_low, v_high].
    """
    return weibull_cdf(v_high, v_ave, k) - weibull_cdf(v_low, v_ave, k)
