"""Generic electric propulsion design equations."""


def thrust(I_b, V_b, m_molar):
    """Thrust of an electric thruster.

    Compute the ideal thrust of an electric thruster from the beam current and voltage,
    assuming singly charged ions and no beam divergence.

    Reference: Goebel and Katz, equation 2.3-8.

    Arguments:
        I_b (scalar): Beam current [units: ampere].
        V_b (scalar): Beam voltage [units: volt].
        m_molar (scalar): Ion molar mass [units: kilogram mole**-1].

    Returns:
        scalar: Thrust force [units: newton].
    """
    pass

