"""Compute the pressure ratio from a given expansion ratio."""
from scipy.optimize import fsolve
from proptools import nozzle

p_c = 10e6    # Chamber pressure [units: pascal]
gamma = 1.2    # Exhaust heat capacity ratio [units: dimensionless]
exp_ratio = 11.9    # Expansion ratio [units: dimensionless]

# Solve for the exit pressure [units: pascal].
p_e = fsolve(lambda p_e: exp_ratio - nozzle.er_from_p(p_c, p_e, gamma),
             x0=2.)[0]

print 'Exit pressure = {:.0f} kPa'.format(p_e * 1e-3)
