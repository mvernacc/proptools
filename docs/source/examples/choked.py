"""Check that the nozzle is choked and find the mass flow."""
from math import pi
from proptools import nozzle

# Declare engine design parameters
p_c = 10e6    # Chamber pressure [units: pascal]
p_e = 100e3    # Exit pressure [units: pascal]
gamma = 1.2    # Exhaust heat capacity ratio [units: dimensionless]
m_molar = 20e-3    # Exhaust molar mass [units: kilogram mole**1]
T_c = 3000.    # Chamber temperature [units: kelvin]
A_t = pi * (0.1 / 2)**2    # Throat area [units: meter**2]

# Check choking
if nozzle.is_choked(p_c, p_e, gamma):
    print 'The flow is choked'

# Compute the mass flow [units: kilogram second**-1]
m_dot = nozzle.mass_flow(A_t, p_c, T_c, gamma, m_molar)

print 'Mass flow = {:.1f} kg s**-1'.format(m_dot)
