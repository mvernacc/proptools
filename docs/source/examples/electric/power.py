"""Example electric propulsion power calculation."""
import numpy as np
from proptools import electric

F = 52.2e-3    # Thrust [units: newton].
m_dot = 1.59e-6    # Mass flow [units: kilogram second**-1].

# Compute the jet power [units: watt].
P_jet = electric.jet_power(F, m_dot)

# Compute the total efficiency [units: dimensionless].
eta_T = electric.total_efficiency(
    divergence_correction=np.cos(np.deg2rad(10)),
    double_fraction=0.1,
    mass_utilization=0.9,
    electrical_efficiency=0.85)

# Compute the input power [units: watt].
P_in = P_jet / eta_T

print 'Jet power = {:.0f} W'.format(P_jet)
print 'Total efficiency = {:.3f}'.format(eta_T)
print 'Input power = {:.0f} W'.format(P_in)
