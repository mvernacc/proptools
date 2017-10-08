"""Generic electric propulsion design equations."""

from proptools.constants import charge, amu_kg


# Xenon atom mass [units: kilogram].
m_Xe = 131.29 * amu_kg

# Krypton atom mass [units: kilogram].
m_Kr = 83.798 * amu_kg


def thrust(I_b, V_b, m_ion):
    """Thrust of an electric thruster.

    Compute the ideal thrust of an electric thruster from the beam current and voltage,
    assuming singly charged ions and no beam divergence.

    Reference: Goebel and Katz, equation 2.3-8.

    Arguments:
        I_b (scalar): Beam current [units: ampere].
        V_b (scalar): Beam voltage [units: volt].
        m_ion (scalar): Ion mass [units: kilogram].

    Returns:
        scalar: Thrust force [units: newton].
    """
    return (2 * (m_ion / charge) * V_b)**0.5 * I_b


def jet_power(F, m_dot):
    """

    Reference: Goebel and Katz, equation 2.3-4.
    """
    pass


def double_ion_thrust_correction(double_fraction):
    """

    Reference: Goebel and Katz, equation 2.3-14.
    """
    pass


def specific_impulse(V_b, m_ion, divergence_correction=1, double_fraction=1, mass_utilization=1):
    """

    Reference: Goebel and Katz, equation 2.4-8.
    """


def total_efficiency(divergence_correction=1, double_fraction=1, mass_utilization=1,
                     electrical_efficiency=1):
    """

    Reference: Goebel and Katz, equation 2.5-7. 
    """
    pass


def thrust_per_power(I_sp, divergence_correction=1, double_fraction=1, mass_utilization=1):
    """

    Reference: Goebel and Katz, equation 2.5-9
    """
    pass


def stuhlinger_velocity(total_efficiency, t_m, specific_mass):
    """

    Reference: Lozano, equation 3-7
    """
    pass


def optimal_isp_thrust_time(total_efficiency, t_m, specific_mass):
    """

    Reference: Lozano, equation 3-5.
    """


def optimal_isp_delta_v(dv, total_efficiency, t_m, specific_mass,
                        discharge_loss=None, m_ion=None):
    """
    """
    pass
