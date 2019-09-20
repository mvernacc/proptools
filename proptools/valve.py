"""Valve flow calculations.


References:
  [1]  "Flow Calculation for Gases", Ideal Valve, Inc.
        http://www.idealvalve.com/pdf/Flow-Calculation-for-Gases.pdf
"""
import numpy as np
import proptools.constants


def m_dot_to_scfh(m_dot, m_molar):
    """Convert mass flow rate to standard cubic feet per hour.

    Arguments:
        m_dot (scalar): Mass flow rate [units: kilogram second**-1].
        m_molar (scalar): Gas molar mass [units: kilogram mole**-1].
    """
    flow_mole = m_dot / m_molar    # units: mole second**-1
    flow_scfs = flow_mole / proptools.constants.scf_mole    # units: scf second**-1
    flow_scfh = flow_scfs * 3600
    return flow_scfh


def scfh_to_m_dot(flow_scfh, m_molar):
    """Convert standard cubic feet per hour to mass flow rate."""
    flow_scfs = flow_scfh / 3600
    flow_mole = flow_scfs * proptools.constants.scf_mole
    m_dot = flow_mole * m_molar
    return m_dot


def valve_gas_cv(m_dot, p_1, p_2, m_molar, T):
    """Find the required valve Cv for a given mass flow and pressure drop.
    Assumes that a compressible gas is flowing through the valve.
 
    Arguments:
        m_dot (scalar): Mass flow rate [units: kilogram second**-1].
        p_1 (scalar): Inlet pressure [units: pascal].
        p_2 (scalar): Outlet pressure [units: pascal].
        m_molar (scalar): Gas molar mass [units: kilogram mole**-1].
        T (scalar): Gas temperature [units: kelvin].

    Returns:
        scalar: Valve flow coefficient Cv [units: gallon minute**-1 psi**-1].
    """
    # Specific gravity of the gas [units: dimensionless]:
    spec_grav = m_molar / proptools.constants.m_molar_air

    # Convert gas flow to standard cubic feet per hour
    flow_scfh = m_dot_to_scfh(m_dot, m_molar)

    # Determine if the flow is choked.
    # Checking if `p_1 >= 2 * p_2` is suggested by [1].
    # There is a more accurate choked flow criterion which depends
    # on the ratio of specific heats.
    choked = p_1 >= 2 * p_2

    if choked:
        cv = flow_scfh / 0.08821 * (spec_grav * T)**0.5 / p_1
    else:
        cv = flow_scfh / 0.1040 * (spec_grav * T / (p_1**2 - p_2**2))**0.5
    return cv


def valve_gas_pressure(cv, m_dot, p_1, m_molar, T):
    # TODO, use fsolve on valve gas mass flow
    pass


def valve_gas_mass_flow(cv, p_1, p_2, m_molar, T):
    # Specific gravity of the gas [units: dimensionless]:
    spec_grav = m_molar / proptools.constants.m_molar_air

    # Determine if the flow is choked.
    # Checking if `p_1 >= 2 * p_2` is suggested by [1].
    # There is a more accurate choked flow criterion which depends
    # on the ratio of specific heats.
    choked = p_1 >= 2 * p_2

    if choked:
        flow_scfh = cv * 0.08821 * p_1 / (spec_grav * T)**0.5
    else:
        flow_scfh = cv * 0.1040 * ((p_1**2 - p_2**2) / (spec_grav * T))**0.5

    m_dot = scfh_to_m_dot(flow_scfh, m_molar)
    return m_dot


def demo_plots():
    """Demonstrate by plotting mass flow vs pressure."""
    cv = 0.1
    p_1 = np.linspace(101e3, 500e3)    # units: pascal
    p_2 = 101e3    # units: pascal
    m_molar = proptools.constants.m_molar_air
    T = 300    # units: kelvin
    m_dot = np.array([valve_gas_mass_flow(cv, p, p_2, m_molar, T) for p in p_1])
    plt.plot(p_1 * 1e-3, m_dot, label='Cv = {:.2f}'.format(cv))
    plt.xlabel('Inlet pressure [kPa]')
    plt.ylabel('Mass flow [kg/s]')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    demo_plots()
