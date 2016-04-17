''' Tank structure calculations.

Matt Vernacchia
proptools
2015 April 17
'''

from math import pi

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


def knuckle_thickness(p_t, a, stress, weld_eff):
    ''' Knuckle thickness of a spherical tank end.

    Implements eqn 8-15 from Huzel and Huang. The knuckle is the transition from
    the cylindrical section to the tank end, see figure 8-6 in Huzel and Huang.
    
    Arguments:
        p_t: Tank internal pressure (less the ambient pressure) [units: pascal].
        a: Tank radius [units: meter].
        stress: Max allowable stress in thank wall material [units: pascal].
        weld_eff: Weld efficiency, scalar in [0, 1] [units: none].

    Returns:
        knuckle thickness [units: meter].
    '''
    # K is the stress factor, see figure 8-7 in Huzel and Huang. Its value is
    # 0.67 for spherical ends, and increases for ellipsoidal ends.
    K = 0.67
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

