''' Non-simple compressible flow.

Calculate quasi-1D compressible flow properties with varying area, friction, and heat addition.
"One-dimensional compressible flows of calorically perfect gases in which only a single driving potential is
present are called simple flows" [1]. This module implements a numerical solution for non-simple flows,
i.e. flows with multiple driving potentials.

References:
    [1] L. Pekker, "One-Dimensional Compressible Flow in Variable Area Duct with Heat Addition,"
        Air Force Research Laboratory, Edwards, CA, Rep. AFRL-RZ-ED-JA-2010-303, 2010.
        Online: http://www.dtic.mil/dtic/tr/fulltext/u2/a524450.pdf.

    [2] A. Bandyopadhyay and A. Majumdar, "Modeling of Compressible Flow with Friction and Heat
        Transfer using the Generalized Fluid System Simulation Program (GFSSP),"
        Thermal Fluid Analysis Workshop, Cleveland, OH, 2007.
        Online: https://tfaws.nasa.gov/TFAWS07/Proceedings/TFAWS07-1016.pdf

Matt Vernacchia
proptools
2016 Oct 3
'''

from scipy.misc import derivative
from scipy.integrate import odeint
import numpy as np

from proptools import nozzle


def differential(state, x, mdot, c_p, gamma, f_f, f_q, f_A):
    ''' Differential equation for Mach number in non-simple duct flow.

    Arguments:
        state (2-vector): Stagnation temperature [units: kelvin],
            Mach number [units: none].
        x (scalar): Distance from the duct inlet [units: meter].
        mdot (scalar): The mass flow through the duct [units: kilogram second**-1].
        c_p (scalar): Fluid heat capacity at constant pressure [units: joule kilogram**-1 kelvin**-1].
        gamma (scalar): Fluid ratio of specific heats [units: none].
        f_f (function mapping scalar->scalar): The Fanning friction factor
            as a function of distance from the inlet [units: none].
        f_q (function mapping scalar->scalar): The heat transfer into the fluid
            per unit wall area as a function of distance from the inlet
            [units: joule meter**-2].
        f_A (function mapping scalar->scalar): The duct area as a function of
            distance from the inlet [units: meter**2].

    Returns:
        d state / dx
    '''
    T_o = state[0]
    M = state[1]

    # Duct area
    A = f_A(x)
    # Duct diameter
    D = (A / np.pi)**2
    # Duct area derivative
    dA_dx = derivative(f_A, x, dx=1e-6)

    # Use Equation 4 from Bandyopadhyay to find dT_o/dx.
    # dT_o_dx = f_q(x) * np.pi * D / (mdot * c_p)
    dT_o_dx = 0

    # Use Equation 3 from Bandyopadhyay to find dM/dx.
    # Note: There is an error in Bandyopadhyay's equation. The area
    # term should not be multiplied by gamma * M**2. Pekker's equation 10
    # shows the correct area term; but lacks a friction term.
    dM_dx = M * (1 + (gamma - 1) / 2 * M**2) / (1 - M**2) \
        * ( \
            # gamma * M**2 * f_f(x) / D \
            # + (1 + gamma * M**2) / (2 * T_o) * dT_o_dx \
            # - gamma * M**2 / A * dA_dx
             - dA_dx / A
            )
    return [dT_o_dx, dM_dx]


def solve_nonsimple(x, M_in, T_o_in, mdot, c_p, gamma, f_f, f_q, f_A):
    ''' Solve a non-simple flow case

    Arguments:
        state (2-vector): Stagnation temperature [units: kelvin],
            Mach number [units: none].
        x (array): Distances from the duct inlet at which to return solution [units: meter].
        T_o_in (scalar): Inlet stagnation temperature [units: kelvin].
        M_in (scalar): Inlet Mach number [units: none].
        mdot (scalar): The mass flow through the duct [units: kilogram second**-1].
        c_p (scalar): Fluid heat capacity at constant pressure [units: joule kilogram**-1 kelvin**-1].
        gamma (scalar): Fluid ratio of specific heats [units: none].
        f_f (function mapping scalar->scalar): The Fanning friction factor
            as a function of distance from the inlet [units: none].
        f_q (function mapping scalar->scalar): The heat transfer into the fluid
            per unit wall area as a function of distance from the inlet
            [units: joule meter**-2].
        f_A (function mapping scalar->scalar): The duct area as a function of
            distance from the inlet [units: meter**2].

    Returns:
        T_o (array of length len(x)): The stagnation temperature at each station in x [units: none].
        M (array of length len(x)): The Mach number at each station in x [units: none].
    '''
    y = odeint(differential, [T_o_in, M_in], x, args=(mdot, c_p, gamma, f_f, f_q, f_A))
    T_o = y[:,0]
    M = y[:,1]
    return (T_o, M)


def main():
    from matplotlib import pyplot as plt
    def f_f(x):
        return 0
    def f_q(x):
        return 0
    def f_A(x):
        return 1 + x
    mdot = 100
    c_p = 2000
    gamma = 1.4
    R = c_p * (1 - 1 / gamma)
    M_in = 0.2
    T_o_in = 300
    T_in = T_o_in * (1 + (gamma - 1) / 2 * M_in**2)**-1
    v_in = M_in * (gamma * R * T_in)**0.5
    rho_in = mdot / (v_in * f_A(0))
    p_in = rho_in * R * T_in
    p_o_in = p_in * (T_o_in / T_in)**(gamma / (gamma -1))
    x = np.linspace(0, 1)

    A_t = nozzle.throat_area(mdot, p_o_in, T_o_in, gamma, nozzle.R_univ / R)

    if M_in < 1:
        M_area = [nozzle.mach_from_area_subsonic(f_A(xi) / A_t, gamma) for xi in x]
    else:
        M_area = [nozzle.mach_from_er(f_A(xi) / A_t, gamma) for xi in x]

    T_o, M = solve_nonsimple(x, M_in, T_o_in, mdot, c_p, gamma, f_f, f_q, f_A)

    plt.subplot(2,1,1)
    plt.plot(x, T_o)
    plt.xlabel('x [m]')
    plt.ylabel('T_o [K]')

    plt.subplot(2,1,2)
    plt.plot(x, M_area, color='blue', label='Simple area solution', marker='x')
    plt.plot(x, M, color='red', label='Non-simple solution', marker='+')
    plt.xlabel('x [m]')
    plt.ylabel('M [-]')
    plt.legend()

    plt.show()



if __name__ == '__main__':
    main()
