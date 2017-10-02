"""Plot thrust vs ambient pressure."""
from scipy.optimize import fsolve
import numpy as np
from matplotlib import pyplot as plt
import skaero.atmosphere.coesa as atmo
from proptools import nozzle

p_c = 10e6   # Chamber pressure [units: pascal]
p_e = 100e3    # Exit pressure [units: pascal]
p_a = np.linspace(0, 100e3)    # Ambient pressure [units: pascal]
gamma = 1.2    # Exhaust heat capacity ratio [units: dimensionless]
A_t = np.pi * (0.1 / 2)**2    # Throat area [units: meter**2]

# Compute thrust [units: newton]
F = nozzle.thrust(A_t, p_c, p_e, gamma,
                  p_a=p_a, er=nozzle.er_from_p(p_c, p_e, gamma))
   
ax1 = plt.subplot(111)
plt.plot(p_a * 1e-3, F * 1e-3)
plt.xlabel('Ambient pressure $p_a$ [kPa]')
plt.ylabel('Thrust $F$ [kN]')
plt.suptitle('Thrust vs ambient pressure at $p_c = {:.0f}$ MPa, $p_e = {:.0f}$ kPa'.format(
    p_c *1e-6, p_e * 1e-3))

# Add altitude on second axis
ylim = plt.ylim()
ax2 = ax1.twiny()
new_tick_locations = np.array([100, 75, 50, 25, 1])
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(new_tick_locations)


def tick_function(p):
    """Map atmospheric pressure [units: kilopascal] to altitude [units: kilometer]."""
    h_table = np.linspace(84e3, 0)    # altitude [units: meter]
    p_table = atmo.pressure(h_table)    # atmo pressure [units: pascal]
    return np.interp(p * 1e3, p_table, h_table) * 1e-3


ax2.set_xticklabels(['{:.0f}'.format(h) for h in tick_function(new_tick_locations)])
ax2.set_xlabel('Altitude [km]')
ax2.tick_params(axis='y', direction='in', pad=-25)
plt.subplots_adjust(top=0.8)
plt.show()
