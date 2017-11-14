"""Illustrate the chamber pressure equilibrium of a solid rocket motor."""

from matplotlib import pyplot as plt
import numpy as np

p_c = np.linspace(1e6, 10e6)    # Chamber pressure [units: pascal].

# Propellant properties
gamma = 1.26    # Exhaust gas ratio of specific heats [units: dimensionless].
rho_solid = 1510.    # Solid propellant density [units: kilogram meter**-3].
n = 0.5    # Propellant burn rate exponent [units: dimensionless].
a = 2.54e-3 * (6.9e6)**(-n)    # Burn rate coefficient, such that the propellant
# burns at 2.54 mm s**-1 at 6.9 MPa [units: meter second**-1 pascal**-n].
c_star = 1209.    # Characteristic velocity [units: meter second**-1].

# Motor geometry
A_t = 839e-6    # Throat area [units: meter**2].
A_b = 1.25    # Burn area [units: meter**2].

# Compute the nozzle mass flow rate at each chamber pressure.
# [units: kilogram second**-1].
m_dot_nozzle = p_c * A_t / c_star

# Compute the combustion mass addition rate at each chamber pressure.
# [units: kilogram second**-1].
m_dot_combustion = A_b * rho_solid * a * p_c**n

# Plot the mass rates
plt.plot(p_c * 1e-6, m_dot_nozzle, label='Nozzle')
plt.plot(p_c * 1e-6, m_dot_combustion, label='Combustion')
plt.xlabel('Chamber pressure [MPa]')
plt.ylabel('Mass rate [kg / s]')

# Find where the mass rates are equal (e.g. the equilibrium).
i_equil = np.argmin(abs(m_dot_combustion - m_dot_nozzle))
m_dot_equil = m_dot_nozzle[i_equil]
p_c_equil = p_c[i_equil]

# Plot the equilibrium point.
plt.scatter(p_c_equil * 1e-6, m_dot_equil, marker='o', color='black', label='Equilibrium')
plt.axvline(x=p_c_equil * 1e-6, color='grey', linestyle='--')

plt.title('Chamber pressure: stable equilibrium, $n =$ {:.1f}'.format(n))
plt.legend()
plt.show()
