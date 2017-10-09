"""Plot the thrust, power, Isp trade-off."""
from matplotlib import pyplot as plt
import numpy as np
from proptools import electric

eta_T = 0.7    # Total efficiency [units: dimensionless].
I_sp = np.linspace(1e3, 7e3)    # Specific impulse [units: second].

ax1 = plt.subplot(111)
for P_in in [2e3, 1e3, 500]:
    T = P_in * electric.thrust_per_power(I_sp, eta_T)
    plt.plot(I_sp, T * 1e3, label='$P_{{in}} = {:.1f}$ kW'.format(P_in * 1e-3))

plt.xlim([1e3, 7e3])
plt.ylim([0, 290])
plt.xlabel('Specific impulse [s]')
plt.ylabel('Thrust [mN]')
plt.legend()
plt.suptitle('Thrust, Power and $I_{sp}$ (70% efficient thruster)')
plt.grid(True)

# Add thrust/power on second axis.
ax2 = ax1.twiny()
new_tick_locations = ax1.get_xticks()
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(['{:.0f}'.format(tp * 1e6)
                     for tp in electric.thrust_per_power(new_tick_locations, eta_T)])
ax2.set_xlabel('Thrust/power [$\\mathrm{mN \\, kW^{-1}}$]')
ax2.tick_params(axis='y', direction='in', pad=-25)
plt.subplots_adjust(top=0.8)

plt.show()
