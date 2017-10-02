"""Plot thrust vs chmaber pressure."""

import numpy as np
from matplotlib import pyplot as plt
from proptools import nozzle

p_c = np.linspace(1e6, 20e6)    # Chamber pressure [units: pascal]
p_e = 100e3    # Exit pressure [units: pascal]
p_a = 100e3    # Ambient pressure [units: pascal]
gamma = 1.2    # Exhaust heat capacity ratio [units: dimensionless]
A_t = np.pi * (0.1 / 2)**2    # Throat area [units: meter**2]

# Compute thrust [units: newton]
F = nozzle.thrust(A_t, p_c, p_e, gamma)

plt.plot(p_c * 1e-6, F * 1e-3)
plt.xlabel('Chamber pressure $p_c$ [MPa]')
plt.ylabel('Thrust $F$ [kN]')
plt.title('Thrust vs chamber pressure at $p_e = p_a = {:.0f}$ kPa'.format(p_e * 1e-3))
plt.show()
