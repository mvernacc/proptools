"""Effect of expansion ratio on thrust coefficient."""
import numpy as np
from matplotlib import pyplot as plt
from proptools import nozzle

p_c = 10e6    # Chamber pressure [units: pascal]
p_a = 100e3    # Ambient pressure [units: pascal]
gamma = 1.2    # Exhaust heat capacity ratio [units: dimensionless]
p_e = np.linspace(0.4 * p_a, 2 * p_a)    # Exit pressure [units: pascal]

# Compute the expansion ratio and thrust coefficient for each p_e
exp_ratio = nozzle.er_from_p(p_c, p_e, gamma)
C_F = nozzle.thrust_coef(p_c, p_e, gamma, p_a=p_a, er=exp_ratio)

# Compute the matched (p_e = p_a) expansion ratio
exp_ratio_matched = nozzle.er_from_p(p_c, p_a, gamma)

plt.plot(exp_ratio, C_F)
plt.axvline(x=exp_ratio_matched, color='grey')
plt.annotate('matched $p_e = p_a$, $\epsilon = {:.1f}$'.format(exp_ratio_matched),
            xy=(exp_ratio_matched - 0.7, 1.62),
            xytext=(exp_ratio_matched - 0.7, 1.62),
            color='black',
            fontsize=10,
            rotation=90
            )
plt.xlabel('Expansion ratio $\\epsilon = A_e / A_t$ [-]')
plt.ylabel('Thrust coefficient $C_F$ [-]')
plt.title('$C_F$ vs expansion ratio at $p_c = {:.0f}$ MPa, $p_a = {:.0f}$ kPa'.format(
    p_c *1e-6, p_a * 1e-3))
plt.show()
