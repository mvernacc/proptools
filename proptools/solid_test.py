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

        print K_RPE
        print K

        self.assertTrue(abs(K - K_RPE) < 0.1)


if __name__ == '__main__':
    unittest.main()
