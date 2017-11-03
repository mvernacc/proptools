"""Unit tests for nozzle flow."""
import unittest
import numpy as np
from proptools import nozzle
from proptools.constants import R_univ


class TestMassFlow(unittest.TestCase):
    """Unit tests for nozzle.mass_flow."""

    def test_rpe_3_3(self):
        """Test against example problem 3-3 from Rocket Propulsion Elements."""
        T_c = 2800.
        gamma = 1.2
        m_molar = R_univ / 360.
        p_c = 2.039e6
        A_t = 13.32e-4

        m_dot = nozzle.mass_flow(A_t, p_c, T_c, gamma, m_molar)
        self.assertAlmostEqual(1.754, m_dot, places=3)


class TestThrust(unittest.TestCase):
    """Unit tests for nozzle.thrust."""

    def test_rpe_3_3(self):
        """Test against example problem 3-3 from Rocket Propulsion Elements."""
        T_c = 2800.
        gamma = 1.2
        m_molar = R_univ / 360.
        p_c = 2.039e6
        A_t = 13.32e-4
        p_e = 2.549e3

        F = nozzle.thrust(A_t, p_c, p_e, gamma)
        self.assertAlmostEqual(5001., F, places=0)


class TestThrustCoef(unittest.TestCase):
    """Unit tests for nozzle.thrust_coef."""

    def test_hh_1_3(self):
        """Test against example problem 1-3 from Huzel and Huang."""
         # Inputs
        T_c = 3633 # = 6540 R
        p_c = 6.895e6 # = 1000 psia
        p_e = 67.9e3 # = 9.85 psia
        p_a = 101e3 # = 14.7 psia
        m_molar = 22.67e-3
        gamma = 1.2
        er = 12

        # Correct results from Huzel and Huang
        C_F_hh = 1.5918

        C_F = nozzle.thrust_coef(p_c, p_e, gamma, p_a, er)

        self.assertTrue(abs(C_F - C_F_hh) < 0.01)


class TestCStar(unittest.TestCase):
    """Unit tests for nozzle.c_star."""

    def test_hh_1_3(self):
        """Test against example problem 1-3 from Huzel and Huang."""
         # Inputs
        T_c = 3633 # = 6540 R
        p_c = 6.895e6 # = 1000 psia
        p_e = 67.9e3 # = 9.85 psia
        p_a = 101e3 # = 14.7 psia
        m_molar = 22.67e-3
        gamma = 1.2
        er = 12

        # Correct results from Huzel and Huang
        c_star_hh = 1777 # = 5830 ft s**-1

        c_star = nozzle.c_star(gamma, m_molar, T_c)

        self.assertTrue(abs(c_star - c_star_hh) < 5)


class TestMachFormEr(unittest.TestCase):
    """Unit tests for nozzle.mach_from_er."""

    def test_hh_1_12(self):
        """Test mach_from_er against H&H figure 1-12"""
        gamma = 1.4
        er = 4
        M_hh = 3.0
        M = nozzle.mach_from_er(er, gamma)
        self.assertTrue(abs(M - M_hh) < 0.1)


class TestMachFromPr(unittest.TestCase):
    """Unit tests for nozzle.mach_from_pr."""

    def test_rpe_3_1(self):
        """Test mach_from_pr against RPE figure 3-1."""
        M = nozzle.mach_from_pr(1, 1, 1.3)
        M_RPE = 0.
        self.assertTrue(abs(M - M_RPE) < 0.05)
        M = nozzle.mach_from_pr(1, 0.8, 1.3)
        M_RPE = 0.6
        self.assertTrue(abs(M - M_RPE) < 0.05)
        M = nozzle.mach_from_pr(1, 0.1, 1.3)
        M_RPE = 2.2
        self.assertTrue(abs(M - M_RPE) < 0.05)


class TestMachArea(unittest.TestCase):
    """Unit tests for nozzle.area_from_mach."""

    def test_gamma_sonic(self):
        """Test that the sonic area ratio is 1, across range of gamma."""
        for gamma in np.linspace(1.1, 1.6, 10):
           self.assertTrue(np.isclose(nozzle.area_from_mach(1, gamma), 1))

    def test_inverse(self):
        """Test that mach_from_area_subsonic is the inverse of area_from_mach."""
        for gamma in np.linspace(1.1, 1.6, 10):
            for Ar in np.linspace(1.1, 10, 10):
                M = nozzle.mach_from_area_subsonic(Ar, gamma)
                self.assertTrue(np.isclose(nozzle.area_from_mach(M, gamma), Ar))


class PressureFromEr(unittest.TestCase):
    """Unit tests for nozzle.pressure_from_er and nozzle.er_from_p"""

    def test_rpe_3_3(self):
        """Test against example problem 3-3 from Rocket Propulsion Elements."""
        gamma = 1.20    # Given ratio of specific heats [units: dimensionless].
        er_rpe = 60.    # Given area expansion ratio [units: dimensionless].
        p_c_rpe = 2.039e6    # Given chamber pressure [units: pascal].
        p_e_rpe = 2.549e3    # Given exit pressure [units: pascal].

        er = nozzle.er_from_p(p_c_rpe, p_e_rpe, gamma)
        self.assertTrue(np.isclose(er, er_rpe, rtol=1e-2))

        pr = nozzle.pressure_from_er(er_rpe, gamma)
        self.assertTrue(np.isclose(pr, p_e_rpe / p_c_rpe, rtol=1e-2))

    def test_inverse(self):
        """Test that pressure_from_er is the inverse of area_from_mach."""
        for gamma in np.linspace(1.1, 1.6, 5):
            for er in np.linspace(2, 400, 10):
                pr = nozzle.pressure_from_er(er, gamma)
                p_e = 1.
                p_c = p_e / pr
                er_calc = nozzle.er_from_p(p_c, p_e, gamma)
                self.assertTrue(np.isclose(er_calc, er, rtol=1e-3))


if __name__ == '__main__':
    unittest.main()
