"""Electric propulsion design tools.

.. autosummary::

  m_Xe
  m_Kr
  thrust
  jet_power
  double_ion_thrust_correction
  specific_impulse
  total_efficiency
  thrust_per_power
  stuhlinger_velocity
  optimal_isp_thrust_time
  optimal_isp_delta_v
"""

from proptools.constants import amu_kg

from proptools.electric.generic import (
    thrust,
    jet_power,
    double_ion_thrust_correction,
    specific_impulse,
    total_efficiency,
    thrust_per_power,
    stuhlinger_velocity,
    optimal_isp_thrust_time,
    optimal_isp_delta_v
)

m_Xe = 131.29 * amu_kg
"""Xenon atom mass [units: kilogram]."""

m_Kr = 83.798 * amu_kg
"""Krypton atom mass [units: kilogram]."""
