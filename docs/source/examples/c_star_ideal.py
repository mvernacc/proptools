"""Ideal characteristic velocity."""
from proptools import nozzle

gamma = 1.2    # Exhaust heat capacity ratio [units: dimensionless]
m_molar = 20e-3    # Exhaust molar mass [units: kilogram mole**1]
T_c = 3000.    # Chamber temperature [units: kelvin]

# Compute the characteristic velocity [units: meter second**-1]
c_star = nozzle.c_star(gamma, m_molar, T_c)

print 'Ideal characteristic velocity = {:.0f} m s**-1'.format(c_star)
