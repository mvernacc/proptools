"""Generic electric propulsion design equations."""
from __future__ import division
import numpy as np
from scipy.optimize import minimize
from proptools.constants import charge, amu_kg, g


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
    """Jet power of a rocket propulsion device.

    Compute the kinetic power of the jet for a given thrust level and mass flow.

    Reference: Goebel and Katz, equation 2.3-4.

    Arguments:
        F (scalar): Thrust force [units: newton].
        m_dot (scalar): Jet mass flow rate [units: kilogram second**-1].

    Returns:
        scalar: jet power [units: watt].
    """
    return F**2 / (2 *  m_dot)


def double_ion_thrust_correction(double_fraction):
    """Doubly-charged ion thrust correction factor.

    Compute the thrust correction factor for the presence of doubly charged ions in the beam.
    This factor is denoted as :math:`\\alpha` in Goebel and Katz.

    Reference: Goebel and Katz, equation 2.3-14.

    Arguments:
        double_fraction (scalar in [0, 1]): The doubly-charged ion current over the singly-charged
            ion current, :math:`I^{++} / I^+` [units: dimensionless].

    Returns:
        scalar in (0, 1]: The thrust correction factor, :math:`\\alpha` [units: dimensionless].
    """
    if double_fraction < 0 or double_fraction > 1:
        raise ValueError('double_fraction {:.f} is not in [0, 1]'.format(double_fraction))

    return (1  + (0.5)**0.5 * double_fraction) / (1 + double_fraction)


def specific_impulse(V_b, m_ion, divergence_correction=1, double_fraction=1, mass_utilization=1):
    """Specific impulse of an electric thruster.

    If only ``V_b`` and ``m_ion`` are provided, the ideal specific impulse will be computed.
    If ``divergence_correction``, ``double_fraction``, or ``mass_utilization`` are provided,
    the specific impulse will be reduced by the corresponding efficiency factors.

    Reference: Goebel and Katz, equation 2.4-8.

    Arguments:
        V_b (scalar): Beam voltage [units: volt].
        m_ion (scalar): Ion mass [units: kilogram].
        divergence_correction (scalar in (0, 1])): Thrust correction factor for beam divergence
            [units: dimensionless].
        double_fraction (scalar in [0, 1]): The doubly-charged ion current over the singly-charged
            ion current, :math:`I^{++} / I^+` [units: dimensionless].
        mass_utilization (scalar in (0, 1])): Mass utilization efficiency [units: dimensionless].

    Returns:
        scalar: the specific impulse [units: second].
    """
    # Check inputs
    if divergence_correction < 0 or divergence_correction > 1:
        raise ValueError('divergence_correction {:.f} is not in [0, 1]'.format(
            divergence_correction))
    if mass_utilization < 0 or mass_utilization > 1:
        raise ValueError('mass_utilization {:.f} is not in [0, 1]'.format(mass_utilization))

    # Compute the efficiency factor
    efficiency = (divergence_correction * double_ion_thrust_correction(double_fraction)
                  * mass_utilization)

    # Compute the ideal specific impulse
    I_sp_ideal = 1 / g * (2 * (charge / m_ion) * V_b)**0.5

    return efficiency * I_sp_ideal


def total_efficiency(divergence_correction=1, double_fraction=1, mass_utilization=1,
                     electrical_efficiency=1):
    """Total efficiency of an electric thruster.

    The total efficiency is defined as the ratio of jet power to input power:

    :math:`\\eta_T \\equiv \\frac{P_{jet}}{P_{in}}`

    Reference: Goebel and Katz, equation 2.5-7.

    Arguments:
        divergence_correction (scalar in (0, 1])): Thrust correction factor for beam divergence
            [units: dimensionless].
        double_fraction (scalar in [0, 1]): The doubly-charged ion current over the singly-charged ion
            current, :math:`I^{++} / I^+` [units: dimensionless].
        mass_utilization (scalar in (0, 1])): Mass utilization efficiency [units: dimensionless].
        electrical_efficiency (scalar in (0, 1])): Electrical efficiency [units: dimensionless].

    Returns:
        scalar: Total efficiency [units: dimensionless].
    """
    # Check inputs
    if divergence_correction < 0 or divergence_correction > 1:
        raise ValueError('divergence_correction {:.f} is not in [0, 1]'.format(
            divergence_correction))
    if mass_utilization < 0 or mass_utilization > 1:
        raise ValueError('mass_utilization {:.f} is not in [0, 1]'.format(mass_utilization))
    if electrical_efficiency < 0 or electrical_efficiency > 1:
        raise ValueError('electrical_efficiency {:.f} is not in [0, 1]'.format(
            electrical_efficiency))

    gamma = divergence_correction * double_ion_thrust_correction(double_fraction)
    return gamma**2 * mass_utilization * electrical_efficiency


def thrust_per_power(I_sp, total_efficiency=1):
    """Thrust/power ratio of an electric thruster.

    Reference: Goebel and Katz, equation 2.5-9.

    Arguments:
        I_sp (scalar): Specific impulse [units:second].
        total_efficiency (scalar in (0, 1]): The total efficiency of the thruster
            [units: dimensionless].

    Returns:
        scalar: Thrust force per unit power input [units: newton watt**-1].
    """
    # Check inputs
    if total_efficiency < 0 or total_efficiency > 1:
        raise ValueError('total_efficiency {:.f} is not in [0, 1]'.format(total_efficiency))

    return (2 / g) * (total_efficiency / I_sp)


def stuhlinger_velocity(total_efficiency, t_m, specific_mass):
    """Stuhlinger velocity, a characteristic velocity for electric propulsion missions.

    The Stuhlinger velocity is a characteristic velocity for electric propulsion missions,
    and is used in calculation of optimal specific impulse. It is defined as:

    :math:`v_{ch} \\equiv \\sqrt{ \\frac{2 \\eta_T t_m}{\\alpha} }`

    where :math:`\\eta_T` is the thruster total efficiency,
    :math:`t_m` is the mission thrust duration, and
    :math:`\\alpha` is the propulsion and power system specific mass.

    The Stuhlinger velocity is also an approximate upper limit for the :math:`\\Delta v`
    capability of an electric propulsion spacecraft.

    Reference: `Lozano 16.522 notes`_, equation 3-7.

    Arguments:
        total_efficiency (scalar in (0, 1]): The total efficiency of the thruster
            [units: dimensionless].
        t_m: mission thrust duration [units: second].
        specific_mass: specific mass (per unit power) of the thruster and power system
            [units: kilogram watt**1].

    Returns:
        scalar: Stuhlinger velocity [units: meter second**-1].

    .. _Lozano 16.522 notes: https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-522-space-propulsion-spring-2015/lecture-notes/MIT16_522S15_Lecture3-4.pdf
    """
    # Check inputs
    if total_efficiency < 0 or total_efficiency > 1:
        raise ValueError('total_efficiency {:.f} is not in [0, 1]'.format(total_efficiency))

    return (2 * total_efficiency * t_m / specific_mass)**0.5


def optimal_isp_thrust_time(total_efficiency, t_m, specific_mass):
    """Optimal specific impulse for constant thrust / fixed time mission. 

    Reference: `Lozano 16.522 notes`_, equation 3-5.

    Arguments:
        total_efficiency (scalar in (0, 1]): The total efficiency of the thruster
            [units: dimensionless].
        t_m: mission thrust duration [units: second].
        specific_mass (scalar): specific mass (per unit power) of the thruster and power system
            [units: kilogram watt**1].

    Returns:
        scalar: Optimal specific impulse [units: second].

    """
    return stuhlinger_velocity(total_efficiency, t_m, specific_mass) / g


def optimal_isp_delta_v(dv, total_efficiency, t_m, specific_mass,
                        discharge_loss=None, m_ion=None):
    """Optimal specific impulse for a fixed :math:`\\Delta v` mission.

    The dependence of efficiency on specific impulse can optionally be included in the
    optimization.
    If ``discharge_loss`` and ``m_ion`` are not provided, then the total efficiency
    :math:`\\eta_T` has a constant value of the given ``total_efficiency``. If ``discharge_loss``
    and ``m_ion`` are provided, :math:`\\eta_T` is assumed to vary with specific impulse as:

    :math:`\\eta_T = \\frac{\\eta_{T,0}}{1 + (v_L / I_{sp} g)^2}`

    where :math:`\\eta_{T,0}` is the given value of ``total_efficiency``, and :math:`v_L` is
    the loss velocity:

    :math:`v_L = \\sqrt{2 q V_d / m_i}`

    where :math:`V_d` is the discharge loss and :math:`m_i` is the ion mass.

    Reference: `Lozano 16.522 notes`_, equation 3-13.

    Arguments:
        total_efficiency (scalar in (0, 1]): The total efficiency of the thruster
            [units: dimensionless].
        t_m: mission thrust duration [units: second].
        specific_mass: specific mass (per unit power) of the thruster and power system
            [units: kilogram watt**1].
        discharge_loss (scalar, optional): The discharge loss per ion [units: volt].
        m_ion (scalar, optional): Ion mass [units: kilogram].

    Returns:
        scalar: Optimal specific impulse [units: second].
    """
    v_ch = stuhlinger_velocity(total_efficiency, t_m, specific_mass)
    v_L = 0
    if discharge_loss is not None and m_ion is not None:
        v_L = (2 * (charge / m_ion) * discharge_loss)**0.5

    def mass_fraction(c):
        """Optimization objective function.

        Payload + structure mass fraction.

        Reference: Lozano equation 3-13

        Arguments:
            c : effective exhaust velocity, I_sp * g.
        """
        return np.exp(-dv / c) - (c**2 + v_L**2) / (v_ch**2) * (1 - np.exp(-dv / c))

    # Scale factor to avoid numerical problems with minimize.
    scale_factor = 10e3

    # Find the c which maximizes the mass fraction.
    result = minimize(lambda c: -mass_fraction(c * scale_factor), 1.)
    c = result.x[0] * scale_factor

    return c / g


# Put functions in the electric module.
thrust.__module__ = 'proptools.electric'
jet_power.__module__ = 'proptools.electric'
double_ion_thrust_correction.__module__ = 'proptools.electric'
specific_impulse.__module__ = 'proptools.electric'
total_efficiency.__module__ = 'proptools.electric'
thrust_per_power.__module__ = 'proptools.electric'
stuhlinger_velocity.__module__ = 'proptools.electric'
optimal_isp_thrust_time.__module__ = 'proptools.electric'
optimal_isp_delta_v.__module__ = 'proptools.electric'
