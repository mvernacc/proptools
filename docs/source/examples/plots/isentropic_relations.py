"""Plot isentropic relations."""
import numpy as np
from matplotlib import pyplot as plt
from proptools import isentropic

M = np.logspace(-1, 1)
gamma = 1.2

plt.loglog(M, 1 / isentropic.stag_temperature_ratio(M, gamma), label='$T / T_0$')
plt.loglog(M, 1 / isentropic.stag_pressure_ratio(M, gamma), label='$p / p_0$')
plt.loglog(M, 1 / isentropic.stag_density_ratio(M, gamma), label='$\\rho / \\rho_0$')
plt.xlabel('Mach number [-]')
plt.ylabel('Static / stagnation [-]')
plt.title('Isentropic flow relations for $\gamma={:.2f}$'.format(gamma))
plt.xlim(0.1, 10)
plt.ylim([1e-3, 1.1])
plt.axvline(x=1, color='grey', linestyle=':')
plt.legend(loc='lower left')
plt.show()