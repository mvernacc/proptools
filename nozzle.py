''' Nozzle flow calculations.

Matt Vernacchia
proptools
2016 Apr 3
'''

# Universal gas constant [units: joule kelvin**-1 mole**-1].
R_univ = 8.314

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


if __name__ == '__main__':
    main()