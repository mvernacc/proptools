"""Models for convective heat transfer."""
from math import pi


def adiabatic_wall_temperature(T_c, M, gamma, r=0.9, Pr=None):
    """Adiabatic wall temperature: the driving temperature for boundary layer convection.

    Reference: Huzel and Huang, equation 4-10-a
    Arguments:
        T_c (scalar): Stagnation temperature of the fluid external to the boundary layer
            [units: kelvin].
        M (scalar): Mach number of the flow [units: dimensionless].
        gamma (scalar): Fluid's ratio of specific heats [units: dimensionless].
        r (scalar): Recovery factor [units: dimensionless]. Value of 0.9
            recommended for turbulent flow.
        Pr (scalar): Prandtl number [units: dimensionless]. If Pr is provided,
            Pr will be used to estimate r (instead of the value of r given).

    Returns:
        scalar: The adiabatic wall temperature :math:`T_{aw}` [units: kelvin].
    """
    if Pr is not None:
        r = Pr**0.33    # Correlation from H&H for turbulent flow.

    return T_c * (1 + r * (gamma - 1) / 2 * M**2) / (1 + (gamma - 1) / 2 * M**2)


def bartz(p_c, c_star, D_t, D, c_p, mu_e, Pr, sigma=1.):
    """Bartz equation for estimation of the convection coefficient.

    The Bartz equation is an empirical estimate of the convection coefficient
    :math:`h_g` of the exhaust flow in a rocket nozzle. The derivation of the
    Bartz equation starts from the modified Reynolds analogy, and uses an
    empirical correlation for turbulent pipe flow to estimate the friction
    factor.

    References:
      [1] Huzel and Huang Equation 4-13.
      [2] M. Martinez-Sanchez, "Convective Heat Transfer: Reynolds Analogy,"
        MIT 16.512 Lecture 7. https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-512-rocket-propulsion-fall-2005/lecture-notes/lecture_7.pdf

    Arguments:
        p_c (scalar): Chamber pressure [units: pascal].
        c_star (scalar): Characteristic velocity [units: meter second**-1].
        D_t (scalar): Throat diameter [units: meter].
        D (scalar): Nozzle diameter where convection is to be computed [units: meter].
        c_p (scalar): Heat capacity at constant pressure of the fluid,
            at stagnation conditions [units: joule kilogram**-1 kelvin**-1].
        mu_e (scalar): dynamic viscosity of the fluid,
            at stagnation conditions [units: pascal second].
        Pr (scalar): Prandtl number of the fluid, at stagnation
            conditions [units: dimensionless].
        sigma (scalar): Correction factor [units: dimensionless].

    Returns:
        scalar: The convection coefficient :math:`h_g` [units: watt meter**-2 kelvin**-1].
    """
    # Note: this neglects the factor (D_T/R)**0.1 present in Huzel's version of the
    # equation, but not present in Sanchez
    # Note: this includes a factor of Pr**0.6, which is present in Huzel but
    # assumed to be unity in Sanchez.
    return (0.026 / D_t**0.2
            * (p_c / c_star)**0.8
            * (D_t / D)**1.8
            * c_p * mu_e**0.2 / Pr**0.6
            * sigma)


def bartz_sigma_sanchez(T_e, T_avg, w=0.6):
    """Correction factor for the Bartz equation.

    Reference:
      [1] M. Martinez-Sanchez, "Convective Heat Transfer: Reynolds Analogy,"
        MIT 16.512 Lecture 7. https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-512-rocket-propulsion-fall-2005/lecture-notes/lecture_7.pdf

    Arguments:
        T_e (scalar): Static temperature of the fluid external to the boundary layer
            [units: kelvin].
        T_avg (scalar): Average temperature through the boundary layer, used to adjust
            the effective density and viscosity of the boundary layer. Martinez-Sanchez
            suggests :math:`T_{avg} = (T_e - T_w)/2` for flow around Mach 1.
        w (scalar): Viscosity is assumed to vary as :math:`T^w`. Martinez-Sanchez
            suggests :math:`w = 0.6` [units: dimensionless].

    Returns:
        scalar: The Bartz equation correction factor :math:`\sigma`.
    """
    return (T_e / T_avg)**(0.8 - 0.2 * w)


def bartz_sigma_huzel(T_c, T_w, M, gamma):
    """Correction factor for the Bartz equation.

    Reference:
      [1] Huzel and Huang Equation 4-14.

    Arguments:
        T_c (scalar): Chamber temperature (flow stagnation temperature)
            [units: kelvin].
        T_w (scalar): Wall temperature [units: kelvin].
        M (scalar): Mach number of the flow [units: dimensionless].
        gamma (scalar): Fluid's ratio of specific heats [units: dimensionless].

    Returns:
        scalar: The Bartz equation correction factor :math:`\sigma`.
    """
    stag = 1 + (gamma - 1) / 2 * M**2    # Stagnation temperature ratio
    return ((0.5 * T_w / T_c * stag + 0.5)**0.68
            * (stag)**0.12)**-1


def film_adiabatic_wall_temperature(eta_film, T_aw, T_f):
    """The effective adiabatic wall temperature in the presence of film cooling.

    The heat flux with film cooling can be approximated as:

    .. math::
      q_w = h_g (T_{aw}^f - T_w)

    Where :math:`T_{aw}^f`is effective adiabatic wall temperature in the presence of film cooling,
    and :math:`h_g` is the convection coefficient computed as if there were no film.

    Reference:
      [1] M. Martinez-Sanchez, "Ablative Cooling, Film Cooling,"
        MIT 16.512 Lecture 10.
    """
    T_aw_f = T_aw - eta_film * (T_aw - T_f)
    return T_aw_f


def film_efficiency(x, D, m_dot_core, m_dot_film, mu_core, Pr_film=1, film_param=1, cp_ratio=1):
    """The film efficiency for a film cooling layer.

    Arguments:
        x (scalar): Distance downstream from the film injection location [units: meter].
        D (scalar): Passage diameter [units: meter].
        m_dot_core (scalar): Core mass flow [units: kilogram second**-1].
        m_dot_film (scalar): Film mass flow [units: kilogram second**-1].
        mu_core (scalar): Dynamic viscosity of the core fluid [units: pascal second].
        Pr_film (scalar, optional): Prandtl number of the film fluid [units: dimensionless].
        film_param (scalar, optional): The film parameter, which is the ratio of the
            film density to the core density: :math:`M_f = \rho_f / \rho_c`
            [units: dimensionless].
        cp_ratio (scalar, optional): The ratio of the core specific heat capacity to the film
            specific heat capacity, :math:`c_{p,c}/c_{p,f}` [units: dimensionless].

    References:
      [1] M. Martinez-Sanchez, "Ablative Cooling, Film Cooling,"
        MIT 16.512 Lecture 10.
      [2] W. Rohsenow, J Hartnett, and Y Cho, "Handbook of Heat Transfer,"
        Ch 17
    """
    # Estimate the film thickness [units: meter].
    s = D / 4 * m_dot_film / m_dot_core / film_param

    # Film flow area [units: meter**2].
    A_film = pi * D * s

    # Film flow mass flux [units: kilogram meter**-2 second**-1].
    flux_film = m_dot_film / A_film

    # Zeta parameter, dimensionless distance downstream from the film injection.
    zeta = x / (film_param * s) * (flux_film * s / mu_core)**(-0.25)

    # Find the film efficiency by an empirical correlation with the zeta distance and
    # the Prandtl number. Martinez-Sanchez recommends this correlation, which is
    # from Rohsenow.
    eta = 1.9 * Pr_film**(0.66) / (1 + 0.329 * cp_ratio * zeta**0.8)
    eta = min([1.0, eta])

    return eta
