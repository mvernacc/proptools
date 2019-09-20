''' Nozzle flow calculations.

Matt Vernacchia
proptools
2016 Apr 3
'''
from __future__ import division
from proptools import nozzle
import math
import numpy as np
from scipy.interpolate import RectBivariateSpline, interp1d


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


def turbine_enthalpy(p_o, p_e, T_o, gamma, c_p):
    ''' Get the specific enthalpy drop for a turbine.

    Arguments:
        p_o: Turbine inlet stagnation pressure [units: same as p_e].
        p_e: Turbine exit pressure [units: same as p_o].
        T_o: Turbine inlet stagnation temperature [units: kelvin].
        gamma: Turbine working gas ratio of specific heats [units: none].
        c_p: working gas heat capacity at const pressure
            [units: joule kilogram**-1 kelvin**-1].

    Returns:
        The specific enthalpy drop across the turbine [units: joule kilogram**-1].
    '''
    # Turbine specific enthalpy drop [units: joule kilogram**-1]
    # Eqn. 6-13 in Huzel and Huang, assumes ideal and calorically perfect gas.
    return c_p * T_o * (1 - (p_e / p_o)**((gamma -1) / gamma))


def turbine_spout_velocity(p_o, p_e, T_o, gamma, c_p):
    ''' Get the theoretical spouting velocity for a turbine.

    Arguments:
        p_o: Turbine inlet stagnation pressure [units: same as p_e].
        p_e: Turbine exit pressure [units: same as p_o].
        T_o: Turbine inlet stagnation temperature [units: kelvin].
        gamma: Turbine working gas ratio of specific heats [units: none].
        c_p: working gas heat capacity at const pressure
            [units: joule kilogram**-1 kelvin**-1].

    Returns:
        The theoretical spouting velocity c_o of the turbine [units: mater second**-1].
    '''
    dh_turb = turbine_enthalpy(p_o, p_e, T_o, gamma, c_p)
    c_o = (2 * dh_turb)**0.5
    return c_o


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
    dh_turb = turbine_enthalpy(p_o, p_e, T_o, gamma, c_p)

    # Turbine output power [units: watt]
    power_turb = eta * m_dot * dh_turb
    return power_turb


def turbine_exit_temperature(p_o, p_te, T_o, eta, gamma, c_p):
    '''Get the turbine exit temperature.
    
    Arguments:
        p_o: turbine inlet stagnation pressure [units: pascal].
        p_te: turbine exit pressure [units: pascal].
        T_o: turbine inlet stagnation temperature [units: kelvin].
        eta: turbine efficiency.
        gamma: working gas ratio of specific heats [units: none].
        c_p: working gas heat capacity at const pressure
            [units: joule kilogram**-1 kelvin**-1].
    '''
    # Turbine specific enthalpy drop [units: joule kilogram**-1]
    dh_turb_ideal = turbine_enthalpy(p_o, p_te, T_o, gamma, c_p)
    dh_turb = eta * dh_turb_ideal
    # Turbine exit temperature
    T_te = T_o - (dh_turb / c_p)
    return T_te


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
    T_te = turbine_exit_temperature(p_o, p_te, T_o, eta, gamma, c_p)

    # Dump nozzle thrust coefficient and characteristic velocity.
    # Assume optimal expansion.
    C_f = nozzle.thrust_coef(p_c=p_te, p_e=p_ne, gamma=gamma)
    c_star = nozzle.c_star(gamma=gamma, m_molar=m_molar, T_c=T_te)

    # Dump nozzle specific impulse [units: second]
    Isp = c_star * C_f / 9.81

    return Isp



def m_dot2gpm(m_dot, rho):
    ''' Convert mass flow to gallons per minute.

    Arguments:
        m_dot: mass flow [units: kilogram second**-1].
        rho: density [units: kilogram meter**-3].

    Returns:
        Volume flow [gallon minute**-1].
    '''
    m3_per_minute = m_dot / rho * 60
    liter_per_minute = m3_per_minute * 1000
    gpm = liter_per_minute * 0.2642
    return gpm


def gpm2m_dot(gpm, rho):
    ''' Convert gallons per minute to mass flow.

    Arguments:
        gpm: Volume flow [gallon minute**-1].
        rho: density [units: kilogram meter**-3].

    Returns:
        m_dot: mass flow [units: kilogram second**-1].
    '''

    liter_per_minute = gpm / 0.2642
    m3_per_minute = liter_per_minute / 1000
    m_dot = m3_per_minute * rho / 60
    return m_dot


def dp2head(dp, rho):
    ''' Convert pump pressure rise to US-units head.

    Arguments:
        dp: Pump pressure rise [units: pascal].
        rho: density [units: kilogram meter**-3].

    Returns:
        pump head [units: feet].
    '''
    return dp / (rho * nozzle.g) / 0.0254 / 12


def radsec2rpm(radsec):
    ''' Convert radian second**-1 to rpm.
    '''
    return radsec / (2 * math.pi) * 60


def rpm2radsec(rpm):
    ''' Convert rpm to radian second**-1.
    '''
    return rpm * (2 * math.pi) / 60


def pump_specific_speed_us(dp, m_dot, rho, N):
    ''' Pump specific speed N_s in US units.

    Arguments:
        dp: Pump pressure rise [units: pascal].
        m_dot: Pump mass flow [units: kilogram second**-1].
        rho: Density of pumped fluid [units: kilogram meter**-3].
        N: Pump rotation speed [radian second**-1]. 

    Returns:
        N_s [units: rpm gallon**0.5 minute**-0.5 feet**-0.75)].
    '''
    N_rpm = radsec2rpm(N)
    # Volume flow rate [units: gallon minute**-1].
    Q = m_dot2gpm(m_dot, rho)
    # Pump head [units: feet]
    H = dp2head(dp, rho)
    N_s = N_rpm * Q**0.5 / H**0.75
    return N_s


# Pump efficiency interpolation.
# Based on figure 6-23 in Huzel and Huang.
gpm = np.array([50, 150, 350, 750])
Ns_us = np.array([500, 700, 900, 1500, 3000, 5000, 7000, 9000, 15000])
eta = np.array([
    [46.51, 53.74, 58.40, 64.56, 65.69, 63.91, 62.23, 60.72, 56.58],
    [52.58, 58.67, 62.86, 69.25, 71.00, 69.09, 67.25, 65.61, 61.32],
    [62.07, 66.41, 69.93, 75.66, 76.74, 74.50, 72.47, 70.69, 66.16],
    [67.72, 72.77, 75.82, 81.07, 82.34, 79.62, 77.22, 75.00, 69.86],
    ]) / 100.0

pump_eff_interpolator = RectBivariateSpline(gpm, Ns_us, eta)

def pump_efficiency(dp, m_dot, rho, N):
    ''' Pump efficiency estimate.

    Based on figure 6-23 in Huzel and Huang.

    Arguments:
        dp: Pump pressure rise [units: pascal].
        m_dot: Pump mass flow [units: kilogram second**-1].
        rho: Density of pumped fluid [units: kilogram meter**-3].
        N: Pump rotation speed [radian second**-1]. 

    Returns:
        pump efficiency  [units: none].
    '''
    # Volume flow rate [units: gallon minute**-1].
    Q = m_dot2gpm(m_dot, rho)
    # Specific speed in US units
    Ns_us = pump_specific_speed_us(dp, m_dot, rho, N)
    eta = pump_eff_interpolator(Q, Ns_us)
    eta = eta[0][0]
    return eta


def pump_efficiency_demo():
    plt.figure()
    rho = 1000
    m_dot = gpm2m_dot(np.array([50, 100, 150, 350]), rho)
    N = rpm2radsec(np.linspace(10e3, 400e3))
    dp = 10e6
    print(dp2head(dp, rho))
    for m in m_dot:
        eta = np.zeros(len(N))
        Ns_us = np.zeros(len(N))
        for i in range(len(N)):
            eta[i] = pump_efficiency(dp, m, rho, N[i])
            Ns_us[i] = pump_specific_speed_us(dp, m, rho, N[i])
        plt.semilogx(Ns_us, eta, label='Q = {0:.0f} gpm'.format(m_dot2gpm(m, rho)))
    plt.xlabel('Ns [US units]')
    plt.ylabel('$\eta$ [-]')
    plt.legend()
    plt.grid(True)


# Turbine efficiency interpolation.
# From figure 10-9 in Rocket Propulsion Elements.
# Single-stage impulse turbine.
uco = np.array([0.03, 0.04, 0.05, 0.07, 0.09, 0.11, 0.14, 0.17, 0.20,
    0.25, 0.28, 0.33, 0.39, 0.45, 0.49, 0.54, 0.56])
eta = np.array([10.10, 13.22, 16.68, 20.10, 27.36, 32.01, 39.90,
    47.01, 53.80, 62.62, 68.06, 74.53, 77.71, 77.53, 75.85, 72.59,
    69.81]) / 100.0
ssi_turbine_eff_interpolator = interp1d(uco, eta)

def ssi_turbine_efficiency(uco):
    ''' Efficiency of a single-stage impulse turbine.

    Data from 10-9 in Rocket Propulsion Elements.

    Arguments:
        uco: Velocity ratio u / c_o [units: none].

    Returns:
        turbine efficiency [units: none].
    '''
    return float(ssi_turbine_eff_interpolator(uco))


def turbine_efficiency_demo():
    plt.figure()
    uco = np.logspace(np.log10(0.04), np.log10(0.5))
    eta = [ssi_turbine_efficiency(u) for u in uco]
    plt.loglog(uco, eta)
    plt.xlabel('$u / c_o$')
    plt.ylabel('$\eta$')
    plt.grid(True)


if __name__ == '__main__':
    from matplotlib import pyplot as plt
    pump_efficiency_demo()
    turbine_efficiency_demo()
    plt.show()
