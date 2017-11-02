""" Solid rocket motor equations.
"""

def chamber_pressure(K_n, a, n, rho_solid, c_star):
    """ Chamber pressure due to solid propellant combustion.

    See equation 12-6 in Rocket Propulsion Elements 8th edition.

    Args:
        K (scalar): Ratio of burning area to throat area, A_b/A_t [units: dimensionless].
        a (scalar): Propellant burn rate coefficient [units: meter second**-1 pascal**-n].
        n (scalar): Propellant burn rate exponent [units: dimensionless].
        rho_solid (scalar): Solid propellant density [units: kilogram meter**-3].
        c_star (scalar): Propellant combustion charateristic velocity [units: meter second**-1].

    Returns:
        Chamber pressure [units: pascal].
    """
    return (K_n * rho_solid * a * c_star) ** (1 / (1 - n))


def burn_area_ratio(p_c, a, n, rho_solid, c_star):
    """ Get the burn area ratio, given chamber pressure and propellant properties.

    Reference: Equation 12-6 in Rocket Propulsion Elements 8th edition.

    Arguments:
        p_c (scalar): Chamber pressure [units: pascal].
        a (scalar): Propellant burn rate coefficient [units: meter second**-1 pascal**-n].
        n (scalar): Propellant burn rate exponent [units: none].
        rho_solid (scalar): Solid propellant density [units: kilogram meter**-3].
        c_star (scalar): Propellant combustion charateristic velocity [units: meter second**-1].

    Returns:
        scalar: Ratio of burning area to throat area, :math:`K = A_b/A_t` [units: dimensionless].
    """
    return p_c**(1 - n) / (rho_solid * a * c_star)
