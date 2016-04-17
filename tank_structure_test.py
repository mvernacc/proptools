import unittest
import tank_structure as ts
from units import inch2meter, psi2pascal


class TestStringMethods(unittest.TestCase):

    def test_sample_8_3(self):
        # Do sample problem 8-3 from Huzel and Huang.
        stress = psi2pascal(38e3)
        a = inch2meter(41.0)
        l_c = inch2meter(46.9)
        E = psi2pascal(10.4e6)
        v = 0.36
        weld_eff = 1.0
        p_to = psi2pascal(180) # Oxidizer max pressure 180 psi
        p_tf = psi2pascal(170) # fuel max pressure 170 psi

        # Ox tank crown 0.135 inch thick.
        self.assertAlmostEqual(inch2meter(0.135), ts.crown_thickness(
            p_to, 1.395*a, stress, weld_eff), delta=inch2meter(0.005))
        # Fuel tank cylinder 0.183 inch thick.
        self.assertAlmostEqual(inch2meter(0.183), ts.cylinder_thickness(
            p_tf, a, stress, weld_eff), delta=inch2meter(0.005))

        # Critical external loading for fuel tank cylinder 10.8 psi
        self.assertAlmostEqual(psi2pascal(10.8), ts.cr_ex_press_cylinder(
            a, inch2meter(0.183), l_c, E, v), delta=1e3)


if __name__ == '__main__':
    unittest.main()