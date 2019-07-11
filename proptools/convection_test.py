"""Unit tests for convection models."""
import unittest
from math import pi
from proptools import convection, nozzle


class TestTaw(unittest.TestCase):
    """Unit tests for convection.adiabatic_wall_temperature"""

    def test_ssme(self):
        """Test against the SSME example from 16.512 Lecture 7."""
        # The T_aw given in the example (3400 K) does not seem to be correct.
        # For now, I am not going to implement this test, for want of a reference
        # example to compare to.
        pass

class TestLongTube(unittest.TestCase):
    """Unit tests for convection.long_tube_coeff."""

    def test_1(self):
        """No examples are given in Hill & Peterson, I had to come up with this."""
        # Setup
        Re = 1e6
        Pr = 0.7
        mass_flux = 1.
        mu = 2e-5
        c_p = 1000.
        k = mu * c_p / Pr
        D = Re * mu / mass_flux

        # Action
        h = convection.long_tube_coeff(mass_flux, D, c_p, mu, k)

        # Verification
        # Lefthand side of equation 11.35 from H & P
        #  = h / (c_p mass_flux)
        # units: dimensionless
        lhs = 0.023 * Re**-0.2 * Pr**-0.67
        self.assertTrue(abs(lhs - h / (c_p * mass_flux)) < 1e-6)


class TestBartz(unittest.TestCase):
    """Unit tests for convection.bartz."""

    def test_ssme(self):
        """Test agaisnt the SSME example from 16.512 Lecture 7."""
        T_e = 3200.
        T_w = 1000.
        T_avg = (T_e + T_w) / 2
        sigma = convection.bartz_sigma_sanchez(T_e, T_avg)
        p_c = 22e6
        c_star = 2600.
        c_p = 2800.
        mu_e = 3e-5
        Pr = 1
        D_t = 0.25    # Not given in example, taken from another reference.

        # Note: the value for h_g given in the lecture notes, 160e3 W m**-2 K**-1
        # appears to be incorrect. Kelly Mathesius and Matt Vernacchia re-did
        # the calculation, and instead got 22e3 W m**-2 K**-1.
        h_g_exp = 22e3    # Answer from the example

        h_g = convection.bartz(p_c, c_star, D_t, D_t, c_p, mu_e, Pr, sigma)

        self.assertTrue(abs(h_g - h_g_exp) < 0.1e3)

    def test_huzel_43_a1_throat(self):
        """Test against Huzel and Huang example problem 4-3, for the A-1 engine
        at the throat."""
        # Values given in the problem statement
        p_c = 6.9e6    # 1000 psi
        c_star = 1725.    # 5660 ft second**-1
        D_t = 0.633    # 24.9 inch
        c_p = 2031.    # 0.485 BTU lbm**-1 rankine**-1
        mu_e =  7.47e-5   # 4.18e-6 lbm inch**-1 second**-1
        Pr = 0.816    # Prandtl number
        sigma = 1.0    # Correction factor

        # Answer given for the convection coefficient
        h_g_huzel = 7950.    # 0.0027 BTU inch**-2 second**-1 fahrenheit**-1

        h_g = convection.bartz(p_c, c_star, D_t, D_t, c_p, mu_e, Pr, sigma)

        self.assertTrue(abs(h_g - h_g_huzel) < 1e3)

    def test_huzel_43_a1_exit(self):
        """Test against Huzel and Huang example problem 4-3, for the A-1 engine
        at the nozzle exit."""
        # Values given in the problem statement
        p_c = 6.9e6    # 1000 psi
        c_star = 1725.    # 5660 ft second**-1
        D_t = 0.633    # 24.9 inch
        D = 5**0.5 * D_t    # Expansion ratio of 5
        c_p = 2031.    # 0.485 BTU lbm**-1 rankine**-1
        mu_e =  7.47e-5   # 4.18e-6 lbm inch**-1 second**-1
        Pr = 0.816    # Prandtl number
        sigma = 0.8    # Correction factor

        # Answer given for the convection coefficient
        h_g_huzel = 1492.    # 0.000507 BTU inch**-2 second**-1 fahrenheit**-1

        h_g = convection.bartz(p_c, c_star, D_t, D, c_p, mu_e, Pr, sigma)

        self.assertTrue(abs(h_g - h_g_huzel) < 0.1e3)


class TestBartzSigmaHuzel(unittest.TestCase):
    """Unit tests for convection.bartz_sigma_huzel."""

    def test_huzel_43_a1_throat(self):
        """Test against Huzel and Huang example problem 4-3, for the A-1 engine
        at the throat."""
        # Values given in the problem statement
        T_c = 3411.    # 6140 rankine
        T_w = 0.8 * T_c   # "Since the carbon-deposit approached to the gas temperature,
        # a T_w/T_c value of 0.8 is used to determine the sigma values."
        M = 1.    # Mach=1 at throat
        gamma = 1.2

        # Answer given for the correction factor
        sigma_huzel = 1.

        sigma = convection.bartz_sigma_huzel(T_c, T_w, M, gamma)

        self.assertTrue(abs(sigma - sigma_huzel) < 0.05)

    def test_huzel_43_a1_exit(self):
        """Test against Huzel and Huang example problem 4-3, for the A-1 engine
        at the exit."""
        # Values given in the problem statement
        T_c = 3411.    # 6140 rankine
        T_w = 0.8 * T_c   # "Since the carbon-deposit approached to the gas temperature,
        # a T_w/T_c value of 0.8 is used to determine the sigma values."
        gamma = 1.2
        M = nozzle.mach_from_er(5., gamma)    # Expansion ratio of 5

        # Answer given for the correction factor
        sigma_huzel = 0.8

        sigma = convection.bartz_sigma_huzel(T_c, T_w, M, gamma)

        self.assertTrue(abs(sigma - sigma_huzel) < 0.05)


class TestFilmCooling(unittest.TestCase):
    """Test the film cooling functions."""

    def test_ms_lec10_1(self):
        """Test against the example from Martinez-Sanchez Lecture 10,
        with 1% film mass flow."""
        x = 0.5    # Distance downstream  [units: meter].
        D = 0.5    # Diameter [units: meter].
        m_dot_core = pi / 4 * D**2 * 5.33 * 253    # Core mass flow [units: kilogram second**-1].
        m_dot_film = (1./99) * m_dot_core
        mu_core = 2e-5 / 0.66   # Dynamic viscosity of the core fluid [units: pascal second].
        Pr_film = 0.8
        film_param = 1.265
        cp_ratio = 0.8

        # Answer from lecture notes
        eta_ms = 0.361

        eta = convection.film_efficiency(x, D, m_dot_core, m_dot_film,
                                         mu_core, Pr_film, film_param, cp_ratio)
        self.assertTrue(abs(eta - eta_ms) < 0.05)

    def test_ms_lec10_10(self):
        """Test against the example from Martinez-Sanchez Lecture 10,
        with 10% film mass flow."""
        x = 0.5    # Distance downstream  [units: meter].
        D = 0.5    # Diameter [units: meter].
        m_dot_core = pi / 4 * D**2 * 5.33 * 253    # Core mass flow [units: kilogram second**-1].
        m_dot_film = (1./9) * m_dot_core
        mu_core = 2e-5 / 0.66   # Dynamic viscosity of the core fluid [units: pascal second].
        Pr_film = 0.8
        film_param = 1.265
        cp_ratio = 0.8

        # Answer from lecture notes
        eta_ms = 1.0

        eta = convection.film_efficiency(x, D, m_dot_core, m_dot_film,
                                         mu_core, Pr_film, film_param, cp_ratio)
        self.assertTrue(abs(eta - eta_ms) < 0.05)

    def test_film_adiabatic_wall_temperature(self):
        """Test against the example from Martinez-Sanchez Lecture 10"""
        T_aw_f = convection.film_adiabatic_wall_temperature(1, 3200, 1600)
        self.assertEqual(T_aw_f, 1600)

        T_aw_f = convection.film_adiabatic_wall_temperature(0.361, 3200, 1600)
        # Error in notes, correct value is 2622 K, not 2296 K (In the notes,
        # 700 K, the wall temp, is incorrectly used in place of 1600 K, the film temp)
        self.assertTrue(abs(2622 - T_aw_f) < 1)


if __name__ == '__main__':
    unittest.main()
