""" Solid rocket motor equations.

.. autosummary::

  chamber_pressure
  burn_area_ratio
  burn_and_throat_area
  thrust_curve
"""

from scipy.integrate import cumtrapz

from proptools import nozzle


def chamber_pressure(K, a, n, rho_solid, c_star):
    """Chamber pressure due to solid propellant combustion.

    Reference: Equation 12-6 in Rocket Propulsion Elements, 8th edition.

    Args:
        K (scalar): Ratio of burning area to throat area, :math:`A_b/A_t` [units: dimensionless].
        a (scalar): Propellant burn rate coefficient [units: meter second**-1 pascal**-n].
        n (scalar): Propellant burn rate exponent [units: dimensionless].
        rho_solid (scalar): Solid propellant density [units: kilogram meter**-3].
        c_star (scalar): Propellant combustion characteristic velocity [units: meter second**-1].

    Returns:
        scalar: Chamber pressure [units: pascal].
    """
    return (K * rho_solid * a * c_star) ** (1 / (1 - n))


def burn_area_ratio(p_c, a, n, rho_solid, c_star):
    """Get the burn area ratio, given chamber pressure and propellant properties.

    Reference: Equation 12-6 in Rocket Propulsion Elements, 8th edition.

    Arguments:
        p_c (scalar): Chamber pressure [units: pascal].
        a (scalar): Propellant burn rate coefficient [units: meter second**-1 pascal**-n].
        n (scalar): Propellant burn rate exponent [units: none].
        rho_solid (scalar): Solid propellant density [units: kilogram meter**-3].
        c_star (scalar): Propellant combustion characteristic velocity [units: meter second**-1].

    Returns:
        scalar: Ratio of burning area to throat area, :math:`K = A_b/A_t` [units: dimensionless].
    """
    return p_c**(1 - n) / (rho_solid * a * c_star)


def burn_and_throat_area(F, p_c, p_e, a, n, rho_solid, c_star, gamma):
    """Given thrust and chamber pressure, and propellant properties, find the burn area and throat area.

    Assumes that the exit pressure is matched (:math:`p_e = p_a`).

    Arguments:
        F (scalar): Thrust force [units: newton].
        p_c (scalar): Chamber pressure [units: pascal].
        p_e (scalar): Nozzle exit pressure [units: pascal].
        a (scalar): Propellant burn rate coefficient [units: meter second**-1 pascal**-n].
        n (scalar): Propellant burn rate exponent [units: none].
        rho_solid (scalar): Solid propellant density [units: kilogram meter**-3].
        c_star (scalar): Propellant combustion characteristic velocity [units: meter second**-1].
        gamma (scalar): Exhaust gas ratio of specific heats [units: dimensionless].

    Returns:
        tuple: tuple containing:

            A_b (scalar): Burn area [units: meter**2].

            A_t (scalar): Throat area [units: meter**2].
    """
    C_F = nozzle.thrust_coef(p_c, p_e, gamma)
    A_t = F / (C_F * p_c)
    A_b = A_t * burn_area_ratio(p_c, a, n, rho_solid, c_star)
    return (A_b, A_t)


def thrust_curve(A_b, x, A_t, A_e, p_a, a, n, rho_solid, c_star, gamma):
    """Thrust vs time curve for a solid rocket motor.

    Given information about the evolution of the burning surface of the propellant grain,
    this function predicts the time-varying thrust of a solid rocket motor.

    The evolution of the burning surface is described by two lists, ``A_b`` and ``x``.
    Each element in the lists describes a step in the (discretized) evolution of the burning
    surface. ``x[i]`` is the distance which the flame front must progress (normal to the burning
    surface) to reach step ``i``. ``A_b[i]`` is the burn area at step ``i``.

    Arguments:
        A_b (list): Burn area at each step [units: meter**-2].
        x (list): flame front progress distance at each step [units: meter].
        A_t (scalar): Nozzle throat area [units: meter**2].
        A_e (scalar): Nozzle exit area [units: meter**2].
        p_a (scalar): Ambient pressure during motor firing [units: pascal].
        a (scalar): Propellant burn rate coefficient [units: meter second**-1 pascal**-n].
        n (scalar): Propellant burn rate exponent [units: none].
        rho_solid (scalar): Solid propellant density [units: kilogram meter**-3].
        c_star (scalar): Propellant combustion characteristic velocity [units: meter second**-1].
        gamma (scalar): Exhaust gas ratio of specific heats [units: dimensionless].

    Returns:
        tuple: tuple containing:

            t (list): time at each step [units: second].

            p_c (list): Chamber pressure at each step [units: pascal].

            F (list): Thrust at each step [units: newton].
    """
    # Compute chamber pressure and exit pressure each flame progress distance x
    # [units: pascal].
    p_c = chamber_pressure(A_b / A_t, a, n, rho_solid, c_star)
    p_e = p_c * nozzle.pressure_from_er(A_e / A_t, gamma)

    # Compute the burn rate for each flame progress distance x [units: meter second**-1]
    r = a * p_c**n

    # Compute the thrust for each flame progress distance x [units: newton]
    F = nozzle.thrust(A_t, p_c, p_e, gamma, p_a, A_e / A_t)

    # Compute the time to reach each flame progress distance x [units: second]
    t = cumtrapz(1 / r, x, initial=0)

    return (t, p_c, F)
