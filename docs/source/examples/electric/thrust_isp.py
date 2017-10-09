"""Example electric propulsion thrust and specific impulse calculation."""
from proptools import electric

V_b = 1000.    # Beam voltage [units: volt].
I_b = 1.    # Beam current [units: ampere].
F = electric.thrust(I_b, V_b, electric.m_Xe)
I_sp = electric.specific_impulse(V_b, electric.m_Xe)
print 'Thrust = {:.1f} mN'.format(F * 1e3)
print 'Specific Impulse = {:.0f} s'.format(I_sp)
