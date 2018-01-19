"""Film cooting example"""
import numpy as np
from matplotlib import pyplot as plt

from proptools import convection

T_aw = 3200   # Adiabatic wall temperature of core flow [units: kelvin].
T_f = 1600    # Film temperature [units: kelvin].
T_w = 700    # Wall temperature [units: kelvin].
x = np.linspace(0, 1)    # Distance downstream  [units: meter].
D = 0.5    # Diameter [units: meter].
m_dot_core = np.pi / 4 * D**2 * 5.33 * 253    # Core mass flow [units: kilogram second**-1].
m_dot_film = (1./99) * m_dot_core
mu_core = 2e-5 / 0.66   # Dynamic viscosity of the core fluid [units: pascal second].
Pr_film = 0.8
film_param = 1.265
cp_ratio = 0.8

eta = np.array([
    convection.film_efficiency(x_, D, m_dot_core, m_dot_film,
                               mu_core, Pr_film, film_param, cp_ratio)
    for x_ in x])
T_aw_f = convection.film_adiabatic_wall_temperature(eta, T_aw, T_f)
q_ratio = (T_aw_f - T_w) / (T_aw - T_w)

plt.plot(x, eta, label='film eff. $\eta$')
plt.plot(x, q_ratio, label='heat flux reduction $q_{w,f} / q_{w,0}$')
plt.xlabel('Distane from injection [m]')
plt.ylabel('[-]')
plt.title('Film cooling in a typical liquid rocket engine')
plt.legend()
plt.grid(True)

plt.show()
