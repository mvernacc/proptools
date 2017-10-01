"""Compute the expansion ratio for a given pressure ratio."""
from proptools import nozzle

p_c = 10e6    # Chamber pressure [units: pascal]
p_e = 100e3    # Exit pressure [units: pascal]
gamma = 1.2    # Exhaust heat capacity ratio [units: dimensionless]

# Solve for the expansion ratio [units: dimensionless]
exp_ratio = nozzle.er_from_p(p_c, p_e, gamma)

print 'Expansion ratio = {:.1f}'.format(exp_ratio)
