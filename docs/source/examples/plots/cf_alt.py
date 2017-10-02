"""Plot C_F vs altitude."""
import numpy as np
from matplotlib import pyplot as plt
import skaero.atmosphere.coesa as atmo
from proptools import nozzle

p_c = 10e6   # Chamber pressure [units: pascal]
gamma = 1.2    # Exhaust heat capacity ratio [units: dimensionless]
p_e_1 = 100e3    # Nozzle exit pressure, 1st stage [units: pascal]
exp_ratio_1 = nozzle.er_from_p(p_c, p_e_1, gamma)    # Nozzle expansion ratio [units: dimensionless]
p_e_2 = 15e3    # Nozzle exit pressure, 2nd stage [units: pascal]
exp_ratio_2 = nozzle.er_from_p(p_c, p_e_2, gamma)    # Nozzle expansion ratio [units: dimensionless]

alt = np.linspace(0, 84e3)    # Altitude [units: meter]
p_a = atmo.pressure(alt)    # Ambient pressure [units: pascal]

# Compute the thrust coeffieicient of the fixed-area nozzle, 1st stage [units: dimensionless]
C_F_fixed_1 = nozzle.thrust_coef(p_c, p_e_1, gamma, p_a=p_a, er=exp_ratio_1)

# Compute the thrust coeffieicient of the fixed-area nozzle, 2nd stage [units: dimensionless]
C_F_fixed_2 = nozzle.thrust_coef(p_c, p_e_2, gamma, p_a=p_a, er=exp_ratio_2)

# Compute the thrust coeffieicient of a variable-area matched nozzle [units: dimensionless]
C_F_matched = nozzle.thrust_coef(p_c, p_a, gamma)

plt.plot(alt * 1e-3, C_F_fixed_1, label='1st stage $\\epsilon_1 = {:.1f}$'.format(exp_ratio_1))
plt.plot(alt[0.4 * p_a < p_e_2] * 1e-3, C_F_fixed_2[0.4 * p_a < p_e_2],
    label='2nd stage $\\epsilon_2 = {:.1f}$'.format(exp_ratio_2))
plt.plot(alt * 1e-3, C_F_matched, label='matched', color='grey', linestyle=':')
plt.xlabel('Altitude [km]')
plt.ylabel('Thrust coefficient $C_F$ [-]')
plt.title('Effect of altitude on nozzle performance')
plt.legend()
plt.show()
