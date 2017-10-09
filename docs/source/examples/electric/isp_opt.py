"""Example calculation of optimal specific impulse for an electrically propelled spacecraft."""
from proptools import electric

dv = 2e3    # Delta-v [units: meter second**-1].
t_m = 100 * 24 * 60 * 60    # Thrust duration [units: second].
eta_T = 0.7    # Total efficiency [units: dimensionless].
specific_mass = 50e-3    # Specific mass of poower supply [units: kilogram watt**-1].
I_sp = electric.optimal_isp_delta_v(dv, eta_T, t_m, specific_mass)

print 'Optimal specific impulse = {:.0f} s'.format(I_sp)
