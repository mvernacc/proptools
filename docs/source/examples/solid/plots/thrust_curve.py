"""Plot the thrust curve of a solid rocket motor with a cylindrical propellant grain."""

from matplotlib import pyplot as plt
import numpy as np
from proptools import solid


# Grain geometry (Clinder with circular port)
r_in = 0.15    # Grain inner radius [units: meter].
r_out = 0.20    # Grain outer radius [units: meter].
length = 1.0    # Grain length [units: meter].

# Propellant properties
gamma = 1.26    # Exhaust gas ratio of specific heats [units: dimensionless].
rho_solid = 1510.    # Solid propellant density [units: kilogram meter**-3].
n = 0.5    # Propellant burn rate exponent [units: dimensionless].
a = 2.54e-3 * (6.9e6)**(-n)    # Burn rate coefficient, such that the propellant
# burns at 2.54 mm s**-1 at 6.9 MPa [units: meter second**-1 pascal**-n].
c_star = 1209.    # Characteristic velocity [units: meter second**-1].

# Nozzle geometry
A_t = 839e-6    # Throat area [units: meter**2].
A_e = 8 * A_t    # Exit area [units: meter**2].
p_a = 101e3    # Ambeint pressure during motor firing [units: pascal].

# Burning surface evolution
x = np.linspace(0, r_out - r_in)    # Flame front progress steps [units: meter].
A_b = 2 * np.pi * (r_in + x) * length    # Burn area at each flame progress step [units: meter**2].

# Compute thrust curve.
t, p_c, F = solid.thrust_curve(A_b, x, A_t, A_e, p_a, a, n, rho_solid, c_star, gamma)

# Plot results.
ax1 = plt.subplot(2, 1, 1)
plt.plot(t, p_c * 1e-6)
plt.ylabel('Chamber pressure [MPa]')

ax2 = plt.subplot(2, 1, 2)
plt.plot(t, F * 1e-3)
plt.ylabel('Thrust, sea level [kN]')
plt.xlabel('Time [s]')
plt.setp(ax1.get_xticklabels(), visible=False)

plt.tight_layout()
plt.subplots_adjust(hspace=0)
plt.show()
