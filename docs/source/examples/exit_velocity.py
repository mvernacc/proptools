"""Estimate exit velocity."""
from proptools import isentropic

# Declare engine design parameters
p_c = 10e6    # Chamber pressure [units: pascal]
p_e = 100e3    # Exit pressure [units: pascal]
gamma = 1.2    # Exhaust heat capacity ratio [units: dimensionless]
m_molar = 20e-3    # Exhaust molar mass [units: kilogram mole**1]
T_c = 3000.    # Chamber temperature [units: kelvin]

# Compute the exit velocity
v_e = isentropic.velocity(v_1=0, p_1=p_c, T_1=T_c, p_2=p_e, gamma=gamma, m_molar=m_molar)

print 'Exit velocity = {:.0f} m s**-1'.format(v_e)
