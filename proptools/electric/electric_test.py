"""Unit tests for proptools.electric."""

import unittest
import random
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


if __name__ == '__main__':
    unittest.main()
