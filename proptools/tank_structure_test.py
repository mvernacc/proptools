import unittest
from proptools import tank_structure as ts
from proptools.units import inch2meter, psi2pascal, lbf2newton, lbm2kilogram


class TestStringMethods(unittest.TestCase):

    def test_sample_8_3(self):
        # Do sample problem 8-3 from Huzel and Huang.
        stress = psi2pascal(38e3)
        a = inch2meter(41.0)
        b = inch2meter(29.4)
        l_c = inch2meter(46.9)
        E = psi2pascal(10.4e6)
        v = 0.36
        weld_eff = 1.0
        p_to = psi2pascal(180) # Oxidizer max pressure 180 psi
        p_tf = psi2pascal(170) # fuel max pressure 170 psi
        rho = 0.101 * lbm2kilogram(1) / inch2meter(1)**3

        # knuckle factor K = 0.80
        self.assertAlmostEqual(0.8, ts.knuckle_factor(a / b), delta=0.02)

        # Ox tank knuckle 0.155 inch thick.
        self.assertAlmostEqual(inch2meter(0.155), ts.knuckle_thickness(
            p_to, a, b, stress, weld_eff), delta=inch2meter(0.005))

        # Ox tank crown 0.135 inch thick.
        self.assertAlmostEqual(inch2meter(0.135), ts.crown_thickness(
            p_to, (a / b) * a, stress, weld_eff), delta=inch2meter(0.005))
        
        # Fuel tank cylinder 0.183 inch thick.
        self.assertAlmostEqual(inch2meter(0.183), ts.cylinder_thickness(
            p_tf, a, stress, weld_eff), delta=inch2meter(0.005))

        # Ellipse design factor E' = 4.56
        self.assertAlmostEqual(4.56, ts.ellipse_design_factor(a / b), delta = 0.005)

        # Oxidizer tank end weighs 126.4 lmb
        self.assertAlmostEqual(2 * lbm2kilogram(126.4), ts.ellipse_mass(
            a, b, inch2meter(0.145), rho), delta=lbm2kilogram(0.2))

        # Cylindrical section weighs 223.3 lbm
        self.assertAlmostEqual(lbm2kilogram(223.3), ts.cylinder_mass(
            a, inch2meter(0.183), l_c, rho), delta=lbm2kilogram(0.2))

        # Critical external pressure for ox tank ends = 13.4 psi.
        self.assertAlmostEqual(psi2pascal(13.4), ts.cr_ex_press_ellipse_end(
            a, b, inch2meter(0.145), E, C_b=0.10), delta=1e3)

        # Critical external loading for fuel tank cylinder 10.8 psi
        self.assertAlmostEqual(psi2pascal(10.8), ts.cr_ex_press_cylinder(
            a, inch2meter(0.183), l_c, E, v), delta=1e3)


    def test_sample_8_4(self):
        # Do sample problem 8-4 from Huzel and Huang.            
        a = inch2meter(41.0)
        l_c = inch2meter(46.9)
        E = psi2pascal(10.4e6)
        t_c = inch2meter(0.183)

        self.assertAlmostEqual(lbf2newton(823900), ts.max_axial_load(
            0, a, t_c, l_c, E), delta=100)

if __name__ == '__main__':
    unittest.main()
