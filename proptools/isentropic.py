"""Isentropic relations.

See :ref:`isentropic-relations-tutorial-label` for a
physical explaination of the isentropic relations.

.. autosummary::

  stag_temperature_ratio
  stag_pressure_ratio
  stag_density_ratio
  velocity
"""
from proptools.constants import R_univ


def stag_temperature_ratio(M, gamma):
    """Stagnation temperature / static temperature ratio.

    Arguments:
        M (scalar): Mach number [units: dimensionless].
        gamma (scalar): Heat capacity ratio [units: dimensionless].

    Returns:
        scalar: the stagnation temperature ratio :math:`T_0 / T` [units: dimensionless].
    """
    return 1 + (gamma - 1) / 2 * M**2


def stag_pressure_ratio(M, gamma):
    """Stagnation pressure / static pressure ratio.

    Arguments:
        M (scalar): Mach number [units: dimensionless].
        gamma (scalar): Heat capacity ratio [units: dimensionless].

    Returns:
        scalar: the stagnation pressure ratio :math:`p_0 / p` [units: dimensionless].
    """
    return (1 + (gamma - 1) / 2 * M**2)**(gamma / (gamma - 1))


def stag_density_ratio(M, gamma):
    """Stagnation density / static density ratio.

    Arguments:
        M (scalar): Mach number [units: dimensionless].
        gamma (scalar): Heat capacity ratio [units: dimensionless].

    Returns:
        scalar: the stagnation density ratio :math:`\\rho_0 / \\rho` [units: dimensionless].
    """
    return (1 + (gamma - 1) / 2 * M**2)**(1 / (gamma - 1))


def velocity(v_1, p_1, T_1, p_2, gamma, m_molar):    # pylint: disable=too-many-arguments
    """Velocity relation between two points in an isentropic flow.

    Given the velocity, pressure, and temperature at station 1 and the pressure at station 2,
    find the velocity at station 2. See Rocket Propulsion Elements, 8th edition, equation 3-15b.

    Arguments:
        v_1 (scalar): Velocity at station 1 [units: meter second**-1].
        p_1 (scalar): Pressure at station 1 [units: pascal].
        T_1 (scalar): Temperature at station 1 [units kelvin].
        p_2 (scalar): Pressure at station 2 [units: pascal].
        gamma (scalar): Gas ratio of specific heats [units: dimensionless].
        m_molar (scalar): Gas mean molar mass [units: kilogram mole**-1].

    Returns:
        scalar: velocity at station 2 [units: meter second**-1].
    """
    return ((2 * gamma) / (gamma - 1) * R_univ * T_1 / m_molar
            * (1 - (p_2 / p_1)**((gamma - 1) / gamma)) + v_1**2)**0.5
