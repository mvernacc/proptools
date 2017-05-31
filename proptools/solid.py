''' Solid rocket motor equations.

Matt Vernacchia
proptools
2016 Aug 22
'''

def chamber_pressure(K_n, a, n, rho_solid, c_star):
    ''' Chamber pressure due to solid propellant combustion.

    See equation 12-6 in Rocket Propulsion Elements 8th edition.

    Args:
        K_n (scalar): Ratio of burning area to throat area, A_b/A_t [units: none].
        a (scalar): Propellant burn rate coefficient [units: meter second**-1 pascal**-n].
        n (scalar): Propellant burn rate exponent [units: none].
        rho_solid (scalar): Solid propellant density [units: kilogram meter**-3].
        c_star (scalar): Propellant combustion charateristic velocity [units: meter second**-1].

    Returns:
        Chamber pressure [units: pascal].
    '''
    return (K_n * rho_solid * a * c_star) ** (1 / (1 - n))
