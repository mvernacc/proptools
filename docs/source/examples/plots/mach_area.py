"""Plot the Mach-Area relation."""
import numpy as np
from matplotlib import pyplot as plt
from proptools import nozzle

M_1 = np.linspace(0.1, 10)    # Mach number [units: dimensionless]
gamma = 1.2    # Exhaust heat capacity ratio [units: dimensionless]

area_ratio = nozzle.area_from_mach(M_1, gamma)

plt.loglog(M_1, area_ratio)
plt.xlabel('Mach number $M_1$ [-]')
plt.ylabel('Area ratio $A_1 / A_2$ [-]')
plt.title('Mach-Area relation for $M_2 = 1$')
plt.ylim([1, 1e3])
plt.show()
