"""Unit test for isentropic relations."""
import unittest
from proptools import isentropic
from proptools.constants import R_univ


class TestStagTemperatureRatio(unittest.TestCase):
    """Unit tests for isentropic.stag_temperature_ratio"""

    def test_still(self):
        """Check that the ratio is 1 when Mach=0."""
        self.assertEqual(1, isentropic.stag_temperature_ratio(0, 1.2))

    def test_sonic(self):
        """Check the ratio when Mach=1."""
        # Rocket Propulsion Elements, 8th edition, page 59.
        self.assertAlmostEqual(1 / 0.91, isentropic.stag_temperature_ratio(1, 1.2), places=2)


class TestStagPressureRatio(unittest.TestCase):
    """Unit tests for isentropic.stag_pressure_ratio"""

    def test_still(self):
        """Check that the ratio is 1 when Mach=0."""
        self.assertEqual(1, isentropic.stag_pressure_ratio(0, 1.2))

    def test_sonic(self):
        """Check the ratio when Mach=1."""
        # Rocket Propulsion Elements, 8th edition, page 59.
        self.assertAlmostEqual(1 / 0.56, isentropic.stag_pressure_ratio(1, 1.2), places=1)


class TestStagDesnityRatio(unittest.TestCase):
    """Unit tests for isentropic.stag_density_ratio"""

    def test_still(self):
        """Check that the ratio is 1 when Mach=0."""
        self.assertEqual(1, isentropic.stag_density_ratio(0, 1.2))

    def test_sonic(self):
        """Check the ratio when Mach=1."""
        # Rocket Propulsion Elements, 8th edition, page 59.
        self.assertAlmostEqual(1.61, isentropic.stag_density_ratio(1, 1.2), places=2)


class TestVelocity(unittest.TestCase):
    """Unit tests for isentropic.velocity."""

    def test_equal_pressure(self):
        """Test that the velocity is the same if the pressure is the same."""
        v_1 = 100.
        p = 1e6
        v_2 = isentropic.velocity(v_1, p_1=p, T_1=300, p_2=p, gamma=1.2, m_molar=20e-3)
        self.assertEqual(v_1, v_2)

    def test_rpe_3_2(self):
        """Test against example problem 3-2 from Rocket Propulsion Elements."""
        v_1 = 0
        p_1 = 2.068e6
        T_1 = 2222.
        gamma = 1.3
        m_molar = R_univ / 345.7
        p_2 = 101e3
        v_2 = isentropic.velocity(v_1, p_1, T_1, p_2, gamma, m_molar)
        self.assertAlmostEqual(v_2, 1828., places=0)


if __name__ == '__main__':
    unittest.main()
