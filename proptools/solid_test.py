"""Unit tests for solid rocket motor equations."""
import unittest
import numpy as np
from proptools import solid


class TestBurnAreaRatio(unittest.TestCase):
    """Unit tests for solid.burn_area_ratio."""

    def test_ex_12_3(self):
        """Test against example problem 12-3 from Rocket Propulsion Elements."""
        gamma = 1.26    # Given ratio of specific heats [units: dimensionless].
        rho_solid = 1510.    # Given solid propellant density [units: kilogram meter**-3].
        p_c = 6.895e6    # Given chamber pressure [units: pascal].
        n = 0.5    # Burn rate exponent, guess
        a = 2.54e-3 * (p_c)**(-n)    # Burn rate coefficient, from given burn rate of
        # 0.1 inch second**-1 at 1000 psi [units: meter second**-1 pascal**-n].
        c_star = 1209.    # Given characteristic velocity [units: meter second**-1].
        K_RPE = 1933. / 1.30    # Given burn area ratio [units: dimensionless].

        K = solid.burn_area_ratio(p_c, a, n, rho_solid, c_star)

        self.assertTrue(abs(K - K_RPE) < 0.1)


class TestBurnAndThroatArea(unittest.TestCase):
    """Unit tests for solid.burn_and_throat_area."""

    def test_ex_12_3(self):
        """Test against example problem 12-3 from Rocket Propulsion Elements."""
        gamma = 1.26    # Given ratio of specific heats [units: dimensionless].
        rho_solid = 1510.    # Given solid propellant density [units: kilogram meter**-3].
        p_c = 6.895e6    # Given chamber pressure [units: pascal].
        n = 0.5    # Burn rate exponent, guess
        a = 2.54e-3 * (p_c)**(-n)    # Burn rate coefficient, from given burn rate of
        # 0.1 inch second**-1 at 1000 psi [units: meter second**-1 pascal**-n].
        c_star = 1209.    # Given characteristic velocity [units: meter second**-1].
        F = 8.9e3    # Given thrust [units: newton].
        p_e = 101e3    # Exit pressure [units: pascal].

        A_t_RPE = 839e-6   # Given throat area [units: meter**2].
        A_b_RPE = 1.25    # Given burn area [units: meter**2].

        A_b, A_t = solid.burn_and_throat_area(F, p_c, p_e, a, n, rho_solid, c_star, gamma)

        self.assertTrue(np.isclose(A_t, A_t_RPE, rtol=5e-2))
        self.assertTrue(np.isclose(A_b, A_b_RPE, rtol=5e-2))


class TestThrustCurve(unittest.TestCase):
    """Unit tests for solid.thrust_curve."""

    def test_ex_12_3_const(self):
        """Test against example problem 12-3 from Rocket Propulsion Elements.

        Constant burn area and thrust."""
        gamma = 1.26    # Given ratio of specific heats [units: dimensionless].
        rho_solid = 1510.    # Given solid propellant density [units: kilogram meter**-3].
        p_c = 6.895e6    # Given chamber pressure [units: pascal].
        n = 0.5    # Burn rate exponent, guess
        a = 2.54e-3 * (p_c)**(-n)    # Burn rate coefficient, from given burn rate of
        # 0.1 inch second**-1 at 1000 psi [units: meter second**-1 pascal**-n].
        c_star = 1209.    # Given characteristic velocity [units: meter second**-1].
        A_t = 839e-6    # Throat area [units: meter**2].
        A_e = 8 * A_t    # Exit area [units: meter**2].
        p_a = 101e3    # Ambeint pressure during motor firing [units: pascal].
        x = np.array([0, 2.54e-3, 2 * 2.54e-3])    # Flame front progress distance.
        A_b = 1.25 * np.ones(3)    # Burn area [units: meter**2].

        F_RPE = 8.9e3    # Given thrust [units: newton].
        p_c_RPE = 6.9e6    # Given chamber pressure [units: pascal].
        t_expected = np.array([0, 1, 2])    # Expected times to reach each step [units: second].

        t, p_c, F = solid.thrust_curve(A_b, x, A_t, A_e, p_a, a, n, rho_solid, c_star, gamma)

        for i in range(3):
            self.assertTrue(abs(p_c[i] - p_c_RPE) < 0.1e6)
            self.assertTrue(abs(F[i] - F_RPE) < 0.3e3)
            self.assertTrue(abs(t[i] - t_expected[i]) < 1e-2)


if __name__ == '__main__':
    unittest.main()
