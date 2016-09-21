''' Nozzle flow calculations.

Matt Vernacchia
proptools
2016 Apr 3
'''
import numpy as np

# Warning logging.
import logging
formatter = logging.Formatter('%(levelname)s {%(pathname)s:%(lineno)d}: %(message)s')
ch = logging.StreamHandler()
ch.setFormatter(formatter)
logger = logging.getLogger('simple_example')
logger.addHandler(ch)

import constants

R_univ = constants.R_univ
g = constants.g

def thrust_coef(p_c, p_e, gamma, p_a=None, er=None):
    '''Nozzle thrust coefficient, C_f.

    Equation 1-33a in Huzel and Huang.

    Arguments:
        p_c: Nozzle stagnation chamber pressure [units: pascal].
        p_e: Nozzle exit static pressure [units: pascal].
        gamma: Exhaust gas ratio of specific heats [units: none].
        p_a (optional): Ambient pressure [units: pascal]. If None,
            then p_a = p_e.
        er (optional): Ambient pressure [units: pascal]. If None,
            then p_a = p_e.


    Returns:
        C_f [units: none].
    '''
    if (p_a is None and er is not None) or (er is None and p_a is not None):
        raise ValueError('Both p_a and er must be provided.')
    C_f = (2 * gamma**2 / (gamma - 1) \
        * (2 / (gamma + 1))**((gamma + 1) / (gamma - 1)) \
        * (1 - (p_e / p_c)**((gamma - 1) / gamma))
        )**0.5
    if p_a is not None and er is not None:
        C_f += er * (p_e - p_a) / p_c
    return C_f


def c_star(gamma, m_molar, T_c):
    '''Characteristic velocity, c*.

    Equation 1-32a in Huzel and Huang. Note that the g in Huzel is removed here,
    because Huzel uses US units while this function uses SI.

    Arguments:
        gamma: Exhaust gas ratio of specific heats [units: none].
        m_molar: Exhaust gas mean molar mass [units: kilogram mole**-1].
        T_c: Nozzle stagnation temperature [units: kelvin].

    Returns:
        c_star [units: meter second**-1].
    '''
    return (gamma * (R_univ / m_molar) * T_c)**0.5 \
        / gamma \
        / (2 / (gamma + 1))**((gamma + 1) / (2 * (gamma - 1)))


def er_from_p(p_c, p_e, gamma):
    ''' Find the nozzle expansion ratio from thechamber and exit pressures.
    
    Arguments:
        p_c: Nozzle stagnation chamber pressure [units: pascal].
        p_e: Nozzle exit static pressure [units: pascal].
        gamma: Exhaust gas ratio of specific heats [units: none].

    Returns:
        Expansion ratio [units: none]
    '''
    # Rocket Propulsion Elements 7th Ed, Equation 3-25
    AtAe = ((gamma + 1) / 2)**(1 / (gamma - 1)) \
        * (p_e / p_c)**(1 / gamma) \
        * ((gamma + 1) / (gamma - 1)*( 1 - (p_e / p_c)**((gamma -1) / gamma)))**0.5
    er = 1/AtAe
    return er


def throat_area(m_dot, p_c, T_c, gamma, m_molar):
    ''' Find the nozzle throat area.

    Arguments:
        p_c: Nozzle stagnation chamber pressure [units: pascal].
        T_c: Nozzle stagnation temperature [units: kelvin].
        gamma: Exhaust gas ratio of specific heats [units: none].
        m_molar: Exhaust gas mean molar mass [units: kilogram mole**-1].

    Returns:
        Throat area [units: meter**2].
    '''
    R = R_univ / m_molar
    # Find the Throat Area require for the specified mass flow, using
    # Eocket Propulsion Equations  7th Ed, Equation 3-24
    A_t = m_dot / ( p_c * gamma \
        * (2 / (gamma + 1))**((gamma + 1) / (2*gamma - 2)) \
        / (gamma * R  * T_c)**0.5)
    return A_t


def mach_from_er(er, gamma):
    '''
    Find the exit Mach number from the area expansion ratio.
    
    Explicit Inversion of Stodola's Area-Mach Equation
    Source: J. Majdalani and B. A. Maickie
    http://maji.utsi.edu/publications/pdf/HT02_11.pdf

    Arguments:
        er: Diverging nozzle area expansion ratio, A_e / A_t [units: none].
        gamma: Exhaust gas ratio of specific heats [units: none].

    Returns:
        The exit Mach number [units: none].
    '''    
    n = 5 # order of the aproximation
    X = np.zeros((n,))
    M = np.zeros((n,))

    e = 1/float(er) # expansion ratio
    y = gamma # ratio of specific heats
    B = (y+1)/(y-1)
    k = np.sqrt( 0.5*(y-1) )
    u = e**(1/B) / np.sqrt( 1+k**2 )
    X[0] = (u*k)**(B/(1-B))
    M[0] = X[0]

    for i in xrange(1,n):
        lamb = 1/( 2*M[i-1]**(2/B)*(B-2) + M[i-1]**2 *B**2*k**2*u**2 )
        X[i] = lamb*M[i-1]*B*( M[i-1]**(2/B) - M[i-1]**2*B*k**2*u**2 \
            + ( M[i-1]**(2+2/B)*k**2*u**2*(B**2-4*B+4) \
            - M[i-1]**2*B**2*k**2*u**4 + M[i-1]**(4/B)*(2*B-3) \
            + 2*M[i-1]**(2/B)*u**2*(2-B) )**0.5 )
        M[i] = M[i-1] + X[i]
    if abs( np.imag( M[n-1] ) ) > 1e-5:
        logger.warning('Exit Mach Number has nonzero imaginary part!')
    Me = float(np.real(M[n-1]))
    return Me


def mach_from_pr(p_c, p_e, gamma):
    ''' Find the exit Mach number from the pressure ratio.

    Arugments:
        p_c: Nozzle stagnation chamber pressure [units: pascal].
        p_e: Nozzle exit static pressure [units: pascal].
        gamma: Exhaust gas ratio of specific heats [units: none].

    Returns:
        scalar: Exit Mach number [units: none].
    '''
    return (2 / (gamma - 1) * ((p_e / p_c)**((1 - gamma) / gamma) -1))**0.5


def is_choked(p_c, p_e, gamma):
    ''' Determine whether the nozzle flow is choked.

    See https://en.wikipedia.org/wiki/Choked_flow#Choking_in_change_of_cross_section_flow

    Arguments:
        p_c: Nozzle stagnation chamber pressure [units: pascal].
        p_e: Nozzle exit static pressure [units: pascal].
        gamma: Exhaust gas ratio of specific heats [units: none].

    Returns:
        True if flow is choked, false otherwise.
    '''
    return p_e/p_c < (2 / (gamma + 1))**(gamma / (gamma - 1))

def main():
    # Do sample problem 1-3 from Huzel and Huang.
    # Inputs
    T_c = 3633 # = 6540 R
    p_c = 6.895e6 # = 1000 psia
    p_e = 67.9e3 # = 9.85 psia
    p_a = 101e3 # = 14.7 psia
    m_molar = 22.67e-3
    gamma = 1.2
    er = 12

    # COrrect results from Huzel and Huang
    C_f_hh = 1.5918
    c_hh = 1777 # = 5830 ft s**-1

    C_f = thrust_coef(p_c, p_e, gamma, p_a, er)
    c = c_star(gamma, m_molar, T_c)
    
    print 'C_f  = {0:.4f}, should be {1:.4f}'.format(C_f, C_f_hh)
    print 'c*  = {0:.0f} m/s, should be {1:.0f} m/s'.format(c, c_hh)
    assert(abs(C_f - C_f_hh) < 0.01)
    assert(abs(c - c_hh) < 5)

    # Test mach_from_er against H&H figure 1-12
    gamma = 1.4
    er = 4
    M_hh = 3.0
    M = mach_from_er(er, gamma)
    print 'M  = {0:.3f}, should be {1:.1f}'.format(M, M_hh)
    assert(abs(M - M_hh) < 0.1)

    # Test mach_from_pr against RPE figure 3-1.
    M = mach_from_pr(1, 1, 1.3)
    M_RPE = 0.
    assert(abs(M - M_RPE) < 0.05)
    M = mach_from_pr(1, 0.8, 1.3)
    M_RPE = 0.6
    assert(abs(M - M_RPE) < 0.05)
    M = mach_from_pr(1, 0.1, 1.3)
    M_RPE = 2.2
    assert(abs(M - M_RPE) < 0.05)

    # Test choked flow
    assert(is_choked(p_c, p_e, gamma))



if __name__ == '__main__':
    main()