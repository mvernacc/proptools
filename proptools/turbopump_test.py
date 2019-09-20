import unittest
from proptools import turbopump

class TestStringMethods(unittest.TestCase):

    def test_m_dot2gpm(self):
        # 100 kg s**-1, 1000 kg m**-3 --> 1584 gal min**-1
        self.assertAlmostEqual(1584, turbopump.m_dot2gpm(100, 1000), delta=2)


    def test_gpm_m_dot_inverse(self):
        for m_dot in [1, 100, 1000, 2434]:
            for rho in [1, 100, 1000, 2000]:
                self.assertAlmostEqual(m_dot, turbopump.gpm2m_dot(
                    turbopump.m_dot2gpm(m_dot, rho), rho))


    def test_rpm_radsec_inverse(self):
        for rpm in [1, 1e3, 1e6]:
            self.assertAlmostEqual(rpm, turbopump.radsec2rpm(
                turbopump.rpm2radsec(rpm)))

    def test_sample61(self):
        # Check against sample problem 6-1 from Huzel and Huang.
        rho_ox = 1143 # 71.38 lbm ft**-3
        rho_fuel = 808.1 # 50.45 lbm ft**-3
        dp_ox = 9.997e6 # 1450 psi
        dp_fuel = 11.55e6 # 1675 psi
        m_dot_ox = 894 # 1971 lbm s**-1
        m_dot_fuel = 405 # 892 lbm s**-1
        N = turbopump.rpm2radsec(7000)

        # Check speed
        self.assertAlmostEqual(7000, turbopump.radsec2rpm(N))

        # Check pump head
        self.assertAlmostEqual(2930, turbopump.dp2head(dp_ox, rho_ox), delta=10)
        self.assertAlmostEqual(4790, turbopump.dp2head(dp_fuel, rho_fuel), delta=10)

        # Check volume flow rate
        self.assertAlmostEqual(12420, turbopump.m_dot2gpm(m_dot_ox, rho_ox), delta=40)
        self.assertAlmostEqual(7960, turbopump.m_dot2gpm(m_dot_fuel, rho_fuel), delta=20)

        # Check specific speed
        self.assertAlmostEqual(1980, turbopump.pump_specific_speed_us(
            dp_ox, m_dot_ox, rho_ox, N),
            delta=30)
        self.assertAlmostEqual(1083, turbopump.pump_specific_speed_us(
            dp_fuel, m_dot_fuel, rho_fuel, N),
            delta=20)


if __name__ == '__main__':
    unittest.main()