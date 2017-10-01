"""Plot thrust vs chmaber pressure."""

import numpy as np
from matplotlib import pyplot as plt
from proptools import nozzle

p_c = 10e6   # Chamber pressure [units: pascal]
p_e = 100e3    # Exit pressure [units: pascal]
p_a = np.linspace(0, 100e3)    # Ambient pressure [units: pascal]
gamma = 1.2    # Exhaust heat capacity ratio [units: dimensionless]
A_t = np.pi * (0.1 / 2)**2    # Throat area [units: meter**2]

C_f = nozzle.thrust_coef(p_c, p_e, gamma,    # Thrust coefficient [units: dimensionless]
                         p_a=p_a, er=nozzle.er_from_p(p_c, p_e, gamma))
F = A_t * p_c * C_f    # Thrust [units: newton]

plt.plot(p_a * 1e-3, F * 1e-3)
plt.xlabel('Ambient pressure $p_a$ [kPa]')
plt.ylabel('Thrust $F$ [kN]')
plt.title('Thrust vs ambient pressure at $p_c = {:.0f}$ MPa, $p_e = {:.0f}$ kPa'.format(
    p_c *1e-6, p_e * 1e-3))
plt.show()
