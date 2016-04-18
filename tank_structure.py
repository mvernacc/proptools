''' Tank structure calculations.

Matt Vernacchia
proptools
2015 April 17
'''

from math import pi
import numpy as np

def crown_thickness(p_t, R, stress, weld_eff):
    ''' Crown thickness of a spherical or ellipsoidal tank end.

    Implements eqn 8-16 from Huzel and Huang. The crown is the center of the
    tank end, see figure 8-6 in Huzel and Huang.
    
    Arguments:
        p_t: Tank internal pressure (less the ambient pressure) [units: pascal].
        R: Crown radius [units: meter].
        stress: Max allowable stress in thank wall material [units: pascal].
        weld_eff: Weld efficiency, scalar in [0, 1] [units: none].

    Returns:
        crown thickness [units: meter].
    '''
    return p_t * R / (2 * stress * weld_eff)


def knuckle_thickness(p_t, a, b, stress, weld_eff):
    ''' Knuckle thickness of a ellipsoidal tank end.

    Implements eqn 8-15 from Huzel and Huang. The knuckle is the transition from
    the cylindrical section to the tank end, see figure 8-6 in Huzel and Huang.
    
    Arguments:
        p_t: Tank internal pressure (less the ambient pressure) [units: pascal].
        a: Tank radius [units: meter].
        a: Semiminor axis [units: meter].
        stress: Max allowable stress in thank wall material [units: pascal].
        weld_eff: Weld efficiency, scalar in [0, 1] [units: none].

    Returns:
        knuckle thickness [units: meter].
    '''
    # K is the stress factor, see figure 8-7 in Huzel and Huang. Its value is
    # 0.67 for spherical ends, and increases for ellipsoidal ends.
    K = knuckle_factor(a / b)
    return K * p_t * a / (stress * weld_eff)


def cylinder_thickness(p_t, a, stress, weld_eff):
    ''' Thickness of a cylindrical tank section.

    Implements eqn 8-28 from Huzel and Huang.
    
    Arguments:
        p_t: Tank internal pressure (less the ambient pressure) [units: pascal].
        a: Tank radius [units: meter].
        stress: Max allowable stress in thank wall material [units: pascal].
        weld_eff: Weld efficiency, scalar in [0, 1] [units: none].

    Returns:
        cylinder thickness [units: meter].
    '''
    return p_t * a / (stress * weld_eff)


def sphere_thickness(p_t, a, stress, weld_eff):
    ''' Thickness of a spherical tank.

    Implements eqn 8-9 from Huzel and Huang.
    
    Arguments:
        p_t: Tank internal pressure (less the ambient pressure) [units: pascal].
        a: Tank radius [units: meter].
        stress: Max allowable stress in thank wall material [units: pascal].
        weld_eff: Weld efficiency, scalar in [0, 1] [units: none].

    Returns:
        sphere thickness [units: meter].
    '''
    return p_t * a / (2 * stress * weld_eff)


def max_axial_load(p_t, a, t_c, l_c, E):
    '''Maximum compressive axial load that a cylindrical section can support.

    Implements eqn 8-33 from Huzel and Huang.

    Arguments:
        p_t: Tank internal pressure (less the ambient pressure) [units: pascal].
        a: Tank radius [units: meter].
        t_c: Cylinder wall thickness [units: meter].
        l_c: Cylinder length [units: meter].
        E: Wall material modulus of elasticity [units: pascal].

    Returns:
        Critical compressive axial load [units: newtons].
    '''
    # Critical axial stress (unpressurized) [units: pascal].
    S_c = (9 * (t_c / a)**1.6 + 0.16 * (t_c / l_c)**1.3) * E
    # Axial cross section area [units: meter**2].
    A_section = 2 * pi * a * t_c
    # Critical axial compression load (unpressurized) [units: newton].
    F_c = S_c *  A_section
    # Axial tension load due to internal pressure [units: newton].
    F_a = pi * a**2 * p_t
    # The total critical compressive load is the load supported by the tank
    # walls plus the load supported by the internal pressure.
    return F_c + F_a


def cylinder_mass(a, t_c, l_c, rho):
    ''' Mass of a cylindrical tank section.

    Arguments:
        a: Tank radius [units: meter].
        t_c: Tank wall thickness [units: meter].
        l_c: Cylinder length [units: meter].
        rho: Tank material density [units: kilogram meter**-3].

    Returns:
        mass of spherical tank [kilogram].
    '''
    return 2 * pi * a * l_c * t_c * rho


def sphere_mass(a, t, rho):
    ''' Mass of a spherical tank.

    Arguments:
        a: Tank radius [units: meter].
        t: Tank wall thickness [units: meter].
        rho: Tank material density [units: kilogram meter**-3].

    Returns:
        mass of spherical tank [kilogram].
    '''
    return 4 * pi * a**2 * t * rho


def ellipse_mass(a, b, t, rho):
    ''' Mass of a ellipsoidal tank.

    Arguments:
        a: Tank radius [units: meter].
        b: Tank semimajor axis [units: meter].
        t: Tank wall thickness [units: meter].
        rho: Tank material density [units: kilogram meter**-3].

    Returns:
        mass of ellipsoidal tank [kilogram].
    '''
    k = a / b
    return pi * a**2 * t * rho * ellipse_design_factor(k) / k


def cr_ex_press_sphere(a, t, E, v):
    ''' Critical external pressure difference to buckle a spherical tank.

    Implements eqn 8-12 in Huzel and Huang.

    Arguments:
        a: Tank radius [units: meter].
        t: Wall thickness [units: meter].
        E: Wall material modulus of elasticity [units: pascal].
        v: Wall material Poisson's ratio [units: none].

    Returns:
        Critical external pressure for buckling [units: pascal].
    '''
    return 2 * E * t**2 / a**2 * (3 * (1 - v**2))**0.5


def cr_ex_press_sphere_end(a, t, E):
    ''' Critical external pressure difference to buckle a spherical tank end.

    Implements eqn 8-26 in Huzel and Huang.

    Arguments:
        a: Tank radius [units: meter].
        t: Wall thickness [units: meter].
        E: Wall material modulus of elasticity [units: pascal].

    Returns:
        Critical external pressure for buckling [units: pascal].
    '''
    return 0.342 * E * t**2 / a**2


def cr_ex_press_ellipse_end(a, b, t, E, C_b=0.05):
    ''' Critical external pressure difference to buckle a ellipsoidal tank end.

    Implements eqn 8-25 in Huzel and Huang.

    Arguments:
        a: Tank radius [units: meter].
        a: Semiminor axis [units: meter].
        t: Wall thickness [units: meter].
        E: Wall material modulus of elasticity [units: pascal].
        Cb: Buckling coefficient [units: none].

    Returns:
        Critical external pressure for buckling [units: pascal].
    '''
    return C_b * 2 * E * t**2 / (a / b * a)**2


def cr_ex_press_cylinder(a, t_c, l_c, E, v):
    ''' Critical external pressure difference to buckle a cylindrical tank section.

    Implements eqn 8-12 in Huzel and Huang.

    Arguments:
        a: Tank radius [units: meter].
        t_c: Tank wall thickness [units: meter].
        l_c: Cylinder length [units: meter].
        E: Wall material modulus of elasticity [units: pascal].
        v: Wall material Poisson's ratio [units: none].

    Returns:
        Critical external pressure for buckling [units: pascal].
    '''
    if l_c < 4.9 * a * (a / t_c)**0.5:
        # Short tank.
        return 0.807 * E * t_c**2 / (l_c * a) \
            * ( (1 / (1 - v**2))**3 * t_c**2 / a**2)**0.25
    else:
        # Long tank.
        return E * t_c**3 / (4 * (1 - v**2) * a**3)


def sphere_volume(a):
    '''Volume enclosed by a spherical tank.

    Arguments:
        a: Tank radius [units: meter].

    Returns:
        tank volume [units: meter**3].
    '''
    return 4.0 / 3 * pi * a**3


def ellipse_volume(a, b):
    '''Volume enclosed by a ellipsoidal tank.

    Arguments:
        a: Tank radius (semimajor axis) [units: meter].
        a: Semiminor axis [units: meter].

    Returns:
        tank volume [units: meter**3].
    '''
    return 4.0 / 3 * pi * a**2 * b


def cylinder_volume(a, l_c):
    '''Volume enclosed by a cylindrical tank section.

    Arguments:
        a: Tank radius [units: meter].
        l_c: Cylinder length [units: meter].

    Returns:
        tank volume [units: meter**3].
    '''
    return pi * a**2 * l_c


def knuckle_factor(ellipse_ratio):
    '''Get the knuckle factor K for an ellipse ratio.

    Implements the "Envelope Curve for K for Combined Stress" curve from 
    figure 8-7 in Huzel and Huang.

    Arguments: 
        ellipse_ratio: Ratio of major and minor axes of ellipse end [units: none].

    Returns:
        knuckle factor K for ellipsoidal end stress calculations [units: none].
    '''
    # Linear fit to:
    # 1.0 -> 0.67
    # 1.75 -. 0.95
    return 0.67 + (ellipse_ratio - 1) * (0.95 - 0.67) / 0.75


def ellipse_design_factor(ellipse_ratio):
    '''Get the ellipse design factor K for an ellipse ratio.

    Implements eqn bs-16 in Huzel and Huang.

    Arguments: 
        ellipse_ratio: Ratio of major and minor axes of ellipse end [units: none].

    Returns:
        ellipse design factor K for ellipsoidal end stress calculations [units: none].
    '''
    k1 = ellipse_ratio
    k2 = (ellipse_ratio**2 - 1)**0.5
    return 2 * k1 + 1 / k2 * np.log((k1 + k2) / (k1 - k2))