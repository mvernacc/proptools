"""Unit tests for nozzle flow."""
import unittest
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


if __name__ == '__main__':
    unittest.main()
