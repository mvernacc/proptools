"""Find the chamber pressure and thrust of a solid rocket motor."""
from proptools import solid, nozzle

# Propellant properties
gamma = 1.26    # Exhaust gas ratio of specific heats [units: dimensionless].
rho_solid = 1510.    # Solid propellant density [units: kilogram meter**-3].
n = 0.5    # Propellant burn rate exponent [units: dimensionless].
a = 2.54e-3 * (6.9e6)**(-n)    # Burn rate coefficient, such that the propellant
# burns at 2.54 mm s**-1 at 6.9 MPa [units: meter second**-1 pascal**-n].
c_star = 1209.    # Characteristic velocity [units: meter second**-1].

# Motor geometry
A_t = 839e-6    # Throat area [units: meter**2].
A_b = 1.25    # Burn area [units: meter**2].

# Nozzle exit pressure [units: pascal].
p_e = 101e3

# Compute the chamber pressure [units: pascal].
p_c = solid.chamber_pressure(A_b / A_t, a, n, rho_solid, c_star)

# Compute the sea level thrust [units: newton].
F = nozzle.thrust(A_t, p_c, p_e, gamma)

print 'Chamber pressure = {:.1f} MPa'.format(p_c * 1e-6)
print 'Thrust (sea level) = {:.1f} kN'.format(F * 1e-3)
