''' Nozzle flow calculations.

Matt Vernacchia
proptools
2016 Apr 3
'''

import nozzle

def pump_power(dp, m_dot, rho, eta):
    ''' Get the input drive power for a pump.

    Arguments:
        dp: Pump pressure rise [units: pascal].
        m_dot: Pump mass flow [units: kilogram second**-1].
        rho: Density of pumped fluid [units: kilogram meter**-3].
        eta: Pump efficiency [units: none].

    Returns:
        The shaft power required to drive the pump [units: watt].
    '''
    return 1 / eta * dp * m_dot / rho


def trubine_power(p_o, p_e, m_dot, T_o, eta, gamma, c_p):
    ''' Get the output drive power for a turbine.

    Arguments:
        p_o: Turbine inlet stagnation pressure [units: same as p_e].
        p_e: Turbine exit pressure [units: same as p_o].
        m_dot: Turbine working gas mass flow [units: kilogram second**-1].
        T_o: Turbine inlet stagnation temperature [units: kelvin].
        gamma: Turbine working gas ratio of specific heats [units: none].
        c_p: working gas heat capacity at const pressure
            [units: joule kilogram**-1 kelvin**-1].

    Returns:
        The shaft power produced by the turbine [units: watt].
    '''
    # Turbine specific enthalpy drop [units: joule kilogram**-1]
    # Eqn. 6-13 in Huzel and Huang, assumes ideal and calorically perfect gas.
    dh_turb = c_p * T_o * (1 - (p_e / p_o)**((gamma -1) / gamma))

    # Turbine output power [units: watt]
    power_turb = eta * m_dot * dh_turb
    return power_turb

def gg_dump_isp(p_o, p_te, p_ne, T_o, eta, gamma, c_p, m_molar):
    '''Get the specific impulse of a Gas Generator turbine exhaust dump.
    
    Arguments:
        p_o: turbine inlet stagnation pressure [units: pascal].
        p_te: turbine exit pressure [units: pascal].
        p_ne: Dump nozzle exit pressure [units: pascal].
        T_o: turbine inlet stagnation temperature [units: kelvin].
        eta: turbine efficiency.
        gamma: working gas ratio of specific heats [units: none].
        c_p: working gas heat capacity at const pressure
            [units: joule kilogram**-1 kelvin**-1].
        m_molar: working gas molar mass [units: kilogram mole**-1].
    '''
    # Turbine specific enthalpy drop [units: joule kilogram**-1]
    # Eqn. 6-13 in Huzel and Huang, assumes ideal and calorically perfect gas.
    dh_turb_ideal = c_p * T_o * (1 - (p_te / p_o)**((gamma -1) / gamma))
    dh_turb = eta * dh_turb_ideal
    # Turbine exit temperature
    T_te = T_o - (dh_turb / c_p)

    # Dump nozzle thrust coefficient and characteristic velocity.
    # Assume optimal expansion.
    C_f = nozzle.thrust_coef(p_c=p_te, p_e=p_ne, gamma=gamma)
    c_star = nozzle.c_star(gamma=gamma, m_molar=m_molar, T_c=T_te)

    # Dump nozzle specific impulse [units: second]
    Isp = c_star * C_f / 9.81

    return Isp
