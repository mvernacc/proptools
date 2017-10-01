"""Estimate specific impulse, thrust and mass flow."""
from math import pi
from proptools import nozzle

# Declare engine design parameters
p_c = 10e6    # Chamber pressure [units: pascal]
p_e = 100e3    # Exit pressure [units: pascal]
gamma = 1.2    # Exhaust heat capacity ratio [units: dimensionless]
m_molar = 20e-3    # Exhaust molar mass [units: kilogram mole**1]
T_c = 3000.    # Chamber temperature [units: kelvin]
A_t = pi * (0.1 / 2)**2    # Throat area [units: meter**2]

# Predict engine performance
C_f = nozzle.thrust_coef(p_c, p_e, gamma)    # Thrust coefficient [units: dimensionless]
c_star = nozzle.c_star(gamma, m_molar, T_c)    # Characteristic velocity [units: meter second**-1]
I_sp = C_f * c_star / nozzle.g    # Specific impulse [units: second]
F = A_t * p_c * C_f    # Thrust [units: newton]
m_dot = A_t * p_c / c_star    # Propellant mass flow [units: kilogram second**-1]

print 'Specific impulse = {:.1f} s'.format(I_sp)
print 'Thrust = {:.1f} kN'.format(F * 1e-3)
print 'Mass flow = {:.1f} kg s**-1'.format(m_dot)
