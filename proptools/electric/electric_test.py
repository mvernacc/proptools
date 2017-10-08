"""Unit tests for proptools.electric."""

import unittest
import random
import numpy as np
from proptools import electric

class TestThrust(unittest.TestCase):
    """Unit tests for electric.thrust."""

    def test_gk_239(self):
        """Test against the rule-of-thumb equation for Xe given in
        Goebel and Katz equation 2.3-9."""
        random.seed(357) # Seed random for test repeatability.
        for i in range(10):
            I_b = random.random()   # Beam current [units: ampere].
            V_b = 100 * random.random()    # Beam voltage [units: volt].
            F_gk_239 = 1.65 * I_b * V_b**0.5 * 1e-3
            F = electric.thrust(I_b, V_b, electric.m_Xe)
            self.assertAlmostEqual(F, F_gk_239, places=1)

class TestJetPower(unittest.TestCase):
    """Unit tests for electric.jet_power."""

    def test_unity(self):
        """Power should be 1/2 if F=1 and m_dot=1."""
        self.assertEqual(0.5, electric.jet_power(1, 1))

    def test_ten(self):
        """Power should be 50 if F=10 and m_dot=1."""
        self.assertEqual(50, electric.jet_power(10, 1))

class TestThrustCorrection(unittest.TestCase):
    """Unit test for electric.double_ion_thrust_correction."""

    def test_gk_example_gamma(self):
        """Test against the example value of gamma given in Goebel and Katz page 24.

        "For example, assuming an ion thruster with a 10-deg half-angle beam divergence and
        # a 10% doubles-to-singles ratio results in gamma=0.958."
        """
        gamma_gk = 0.958
        alpha = electric.double_ion_thrust_correction(0.1)
        self.assertAlmostEqual(gamma_gk, alpha * np.cos(np.deg2rad(10.)), places=2)

    def test_exception(self):
        """Check that the function raises a value error for invalid double_fraction."""
        with self.assertRaises(ValueError):
            electric.double_ion_thrust_correction(-1)


class TestIsp(unittest.TestCase):
    """Unit test for electric.specific_impusle."""

    def test_gk_example(self):
        """Test against the example values given in Goebel and Katz page 27.

        "Using our previous example of a 10-deg half-angle beam divergence and a 10%
        doubles-to-singles ratio with a 90% propellant utilization of xenon at [...] 1500 V,
        the Isp is [...] 4127 s"
        """
        I_sp_gk = 4127    # Correct I_sp value from Goebel and Katz [units: second].
        V_b = 1500    # Beam voltage [units: volt].
        divergence_correction = np.cos(np.deg2rad(10.))    # Beam divergence correction factor. 
        I_sp = electric.specific_impulse(V_b, electric.m_Xe,
                                         divergence_correction=divergence_correction,
                                         double_fraction=0.1,
                                         mass_utilization=0.9)
        self.assertLess(abs(I_sp - I_sp_gk), 3)


    def test_exception(self):
        """Check that the function raises a value error for invalid inputs."""
        with self.assertRaises(ValueError):
            electric.specific_impulse(1, 1, divergence_correction=-1)
        with self.assertRaises(ValueError):
            electric.specific_impulse(1, 1, mass_utilization=-1)


class TestTotalEffciency(unittest.TestCase):
    """Unit tests for electric.total_efficiency"""

    def test_gk_example(self):
        """Test against the example values given in Goebel and Katz page 29.

        "Using our previous example of an ion thruster with 10-deg half-angle divergence,
        10% double ion current, 90% mass utilization efficiency and [...] beam at 1500 V
        [...] the electrical efficiency is [...] 0.857 [...] and the total efficiency is
        [...] 0.708"
        """
        eta_t_gk = 0.708    # Correct total efficiency value from Goebel and Katz [units: dimensionless].
        V_b = 1500    # Beam voltage [units: volt].
        divergence_correction = np.cos(np.deg2rad(10.))    # Beam divergence correction factor.
        eta_t = electric.total_efficiency(divergence_correction=divergence_correction,
                                          double_fraction=0.1,
                                          mass_utilization=0.9,
                                          electrical_efficiency=0.857)
        self.assertAlmostEqual(eta_t_gk, eta_t, places=2)

    def test_exception(self):
        """Check that the function raises a value error for invalid inputs."""
        with self.assertRaises(ValueError):
            electric.total_efficiency(divergence_correction=-1)
        with self.assertRaises(ValueError):
            electric.total_efficiency(mass_utilization=-1)
        with self.assertRaises(ValueError):
            electric.total_efficiency(electrical_efficiency=-1)


class TestThrustPerPower(unittest.TestCase):
    """Unit tests for electric.thrust_per_power."""

    def test_gk_example(self):
        """Test against the example values given in Goebel and Katz.
        "thrust produced is 122.4 mN [...] 2 A beam at 1500 V [...] dissipated power of 528.3 W"
        """
        eta_t = 0.708    # Total efficiency [units: dimensionless].
        I_sp = 4127    # Specific impulse [units: seconds].
        P_in = 2 * 1500 + 528.3    # Input power [units: watt].
        F = 122.4e-3    # Thrust force [units: watt].
        fp = electric.thrust_per_power(I_sp, eta_t)
        self.assertAlmostEqual(fp, F / P_in, places=6)

    def test_exception(self):
        """Check that the function raises a value error for invalid inputs."""
        with self.assertRaises(ValueError):
            electric.total_efficiency(divergence_correction=-1)
        with self.assertRaises(ValueError):
            electric.total_efficiency(mass_utilization=-1)
        with self.assertRaises(ValueError):
            electric.total_efficiency(electrical_efficiency=-1)

if __name__ == '__main__':
    unittest.main()
