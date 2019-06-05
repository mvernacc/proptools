"""Unit tests for valve flow."""

import unittest
from proptools import valve

class TestValveGasCv(unittest.TestCase):
    def test_methane_example(self):
        """Test the methane example from
        http://www.idealvalve.com/pdf/Flow-Calculation-for-Gases.pdf
        """
        # Setup
        p_1 = 790.83e3    # Inlet pressure, 100 psig [units: pascal]
        p_2 = 101e3    # Outlet pressure [units: pascal]
        T = 294.3    # gas temperature, 70 F [units: kelvin]
        m_molar = 0.01604    # Methane molar mass [units: kilogram mole**-1]
        flow_scfh = 600
        m_dot = valve.scfh_to_m_dot(flow_scfh, m_molar)
        cv_given = 0.1098    # Flow coefficient Cv given in example.

        # Action
        cv_calc = valve.valve_gas_cv(m_dot, p_1, p_2, m_molar, T)

        # Verification
        self.assertAlmostEqual(cv_given, cv_calc, places=3)


class TestValveGasMassFlow(unittest.TestCase):
    def test_methane_example(self):
        """Test the methane example from
        http://www.idealvalve.com/pdf/Flow-Calculation-for-Gases.pdf
        """
        # Setup
        p_1 = 790.83e3    # Inlet pressure, 100 psig [units: pascal]
        p_2 = 101e3    # Outlet pressure [units: pascal]
        T = 294.3    # gas temperature, 70 F [units: kelvin]
        m_molar = 0.01604    # Methane molar mass [units: kilogram mole**-1]
        flow_scfh_given = 600    # Flow in scfh given in example
        m_dot_given = valve.scfh_to_m_dot(flow_scfh_given, m_molar)
        cv = 0.1098    # Flow coefficient Cv

        # Action
        m_dot_calc = valve.valve_gas_mass_flow(cv, p_1, p_2, m_molar, T)

        # Verification
        self.assertAlmostEqual(m_dot_given, m_dot_calc, places=3)


if __name__ == '__main__':
    unittest.main()
