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

    [3] J. D. Anderson, Modern Compressible Flow with Historical Perspective, 2nd ed.
        New York, NY: McGraw-Hill, 1990.

Matt Vernacchia
proptools
2016 Oct 3
'''

from scipy.misc import derivative
from scipy.integrate import ode
import numpy as np

from proptools import nozzle


def differential(x, state, mdot, c_p, gamma, f_f, f_q, f_A):
    ''' Differential equation for Mach number in non-simple duct flow.

    Note: This method will not be accurate (and may divide by zero) for flows which
    contain a region at Mach 1, e.g. a choked convergent-divergent nozzle.

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
    D = (A / np.pi)**0.5
    # Duct area derivative
    dA_dx = derivative(f_A, x, dx=1e-6)

    # Use Equation 4 from Bandyopadhyay to find dT_o/dx.
    dT_o_dx = f_q(x) * np.pi * D / (mdot * c_p)

    # Use Equation 3 from Bandyopadhyay to find dM/dx.
    # Note: There is an error in Bandyopadhyay's equation. The area
    # term should not be multiplied by gamma * M**2. Pekker's equation 10
    # shows the correct area term; but lacks a friction term.
    # Note: Compared to Anderson's Eqn 3.96, Bandyopadhyay's equation
    # is missing a factor of 2 on the friction term. I use Anderson's
    # friction term here.
    dM_dx = M * (1 + (gamma - 1) / 2 * M**2) / (1 - M**2) \
        * ( \
            gamma * M**2 * 2 * f_f(x) / D \
            + (1 + gamma * M**2) / (2 * T_o) * dT_o_dx \
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
        choked (boolean): True if the flow chokes at M=1 in the duct. M and T_o for x past the
            choke point will be nan. Choking can cause shocks or upstream effects which this model
            does not capture; therefore results for choked scenarios may not be accurate.
    '''
    T_o = np.full(len(x), np.nan)
    M = np.full(len(x), np.nan)
    solver = ode(differential).set_f_params(mdot, c_p, gamma, f_f, f_q, f_A)
    solver.set_initial_value([T_o_in, M_in], 0)
    T_o[0] = T_o_in
    M[0] = M_in
    i = 1
    while solver.successful() and i < len(x):
        y = solver.integrate(x[i])
        T_o[i] = y[0]
        M[i] = y[1]
        i += 1
    choked = False
    if not solver.successful() and abs(solver.y[1] - 1) < 1e-3:
        choked = True
    return (T_o, M, choked)


def main():
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
    x = np.linspace(0, 1)
    T_o_in = 300

    for M_in in [0.2, 0.8, 1.2, 2]:
        T_in = T_o_in * (1 + (gamma - 1) / 2 * M_in**2)**-1
        v_in = M_in * (gamma * R * T_in)**0.5
        rho_in = mdot / (v_in * f_A(0))
        p_in = rho_in * R * T_in
        p_o_in = p_in * (T_o_in / T_in)**(gamma / (gamma -1))
        

        A_t = nozzle.throat_area(mdot, p_o_in, T_o_in, gamma, nozzle.R_univ / R)

        if M_in < 1:
            M_area = [nozzle.mach_from_area_subsonic(f_A(xi) / A_t, gamma) for xi in x]
        else:
            M_area = [nozzle.mach_from_er(f_A(xi) / A_t, gamma) for xi in x]

        T_o, M, choked = solve_nonsimple(x, M_in, T_o_in, mdot, c_p, gamma, f_f, f_q, f_A)

        plt.subplot(2,1,1)
        plt.plot(x, T_o)
        plt.xlabel('x [m]')
        plt.ylabel('T_o [K]')

        plt.subplot(2,1,2)
        plt.plot(x, M_area, label='Simple area solution $M_{{in}}$={:.1f}'.format(M_in),
            marker='x')
        plt.plot(x, M, label='Non-simple solution $M_{{in}}$={:.1f}'.format(M_in),
            marker='+')
        plt.xlabel('x [m]')
        plt.ylabel('M [-]')
    plt.legend()
    plt.suptitle('Demonstration with linear area variation, no heat addition, no friction')

    # Fanno Flow demo
    plt.figure()
    f = 0.005
    def f_f(x):
        return f
    def f_q(x):
        return 0
    def f_A(x):
        return 1

    for M_in in [0.2, 0.8, 1.2, 2]:
        # Duct diameter
        D = (f_A(0) / np.pi)**0.5
        # Choking length
        # Anderson Modern Compressible Flow Equation 3.107
        L = D / (4 * f) * ((1 - M_in**2) / (gamma * M_in**2) + (gamma + 1) / (2 * gamma) \
            * np.log((gamma + 1) * M_in**2 / (2 + (gamma - 1) * M_in**2)))
        x = np.linspace(0, L)
        T_in = T_o_in * (1 + (gamma - 1) / 2 * M_in**2)**-1
        v_in = M_in * (gamma * R * T_in)**0.5
        rho_in = mdot / (v_in * f_A(0))
        p_in = rho_in * R * T_in
        p_o_in = p_in * (T_o_in / T_in)**(gamma / (gamma -1))

        T_o, M, choked = solve_nonsimple(x, M_in, T_o_in, mdot, c_p, gamma, f_f, f_q, f_A)

        plt.subplot(2,1,1)
        plt.plot(x, T_o)
        plt.xlabel('x [m]')
        plt.ylabel('T_o [K]')

        plt.subplot(2,1,2)
        plt.plot(x, M, label='$M_{{in}}$={:.1f}, {:s}'.format(M_in, 'choked' if choked else 'not choked'),
            marker='+')
        plt.xlabel('x [m]')
        plt.ylabel('M [-]')
    plt.axhline(y=1, color='black')
    plt.legend()
    plt.suptitle('Demonstration of Fanno Flow')

    # Rayleigh Flow demo
    plt.figure()
    q = 1e5
    def f_f(x):
        return 0
    def f_q(x):
        return q
    def f_A(x):
        return 1

    for M_in in [0.2, 0.8, 1.2, 2]:
        x = np.linspace(0, 100)

        T_o, M, choked = solve_nonsimple(x, M_in, T_o_in, mdot, c_p, gamma, f_f, f_q, f_A)

        plt.subplot(2,1,1)
        plt.plot(x, T_o)
        plt.xlabel('x [m]')
        plt.ylabel('T_o [K]')

        plt.subplot(2,1,2)
        plt.plot(x, M, label='$M_{{in}}$={:.1f}, {:s}'.format(M_in, 'choked' if choked else 'not choked'),
            marker='+')
        plt.xlabel('x [m]')
        plt.ylabel('M [-]')
    plt.axhline(y=1, color='black')
    plt.legend()
    plt.suptitle('Demonstration of Rayleigh Flow')
    plt.show()



if __name__ == '__main__':
    from matplotlib import pyplot as plt
    main()
