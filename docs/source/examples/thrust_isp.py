from math import pi
from proptools import nozzle

p_c = 10e6
p_e = 100e3
gamma = 1.2
m_molar = 20e-3
T_c = 3000.
A_t = pi * (0.1 / 2)**2

C_f = nozzle.thrust_coef(p_c, p_e, gamma)
c_star = nozzle.c_star(gamma, m_molar, T_c)

I_sp = C_f * c_star / nozzle.g
F = A_t * p_c * C_f
m_dot = A_t * p_c / c_star

print 'Specific impulse = {:.1f} s'.format(I_sp)
print 'Thrust = {:.1f} kN'.format(F * 1e-3)
print 'Mass flow = {:.1f} kg s**-1'.format(m_dot)
