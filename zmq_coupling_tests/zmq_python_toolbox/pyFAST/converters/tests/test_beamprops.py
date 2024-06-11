import unittest
import os
import numpy as np
from pyFAST.converters.beam import *


class Test(unittest.TestCase):

    def test_MM_translate(self):

        # Inertia param
        m       =  262.
        I_x     =  163.
        I_y     =  54.
        I_p     =  217.
        theta_i = -0.1
        x_G     =  0.175
        y_G     = -0.1
        x_A     = -1
        y_A     =  1
        # Referenec values - Mass and stiffness matrix for simple beam
        M       = MM(m, I_x, I_y, I_p, x_G, y_G, theta_i) # NOTE: theta_i in rad
        M_G_ref = MM(m, I_x, I_y, I_p, 0,   0  , theta_i) 
        M_A_ref = MM(m, I_x, I_y, I_p, x_G-x_A, y_G-y_A, theta_i) # NOTE: theta_i in rad
        J_G_ref = M_G_ref[3:6, 3:6] # Inertia at COG

        # Functions to compute inertial properties
        inertia            = ComputeInertiaProps()
        xm , ym            = inertia.ComputeMassCenter(M)
        mass, J_G, Ref2COG = identifyRigidBodyMM(M)

        np.testing.assert_almost_equal(xm,  x_G, 6)
        np.testing.assert_almost_equal(ym,  y_G, 6)
        np.testing.assert_almost_equal(mass,  m, 6)
        np.testing.assert_almost_equal(Ref2COG[0],  x_G, 6)
        np.testing.assert_almost_equal(Ref2COG[1],  y_G, 6)
        np.testing.assert_almost_equal(Ref2COG[2],  0, 6)
        np.testing.assert_almost_equal(J_G,  J_G_ref, 6)


        # Translating to the COG
        M_G = TranslateSectionMassMatrix(M, x_G, y_G)
        np.testing.assert_almost_equal(M_G,  M_G_ref, 6)

        # Translating to another point
        M_A = TranslateSectionMassMatrix(M, x_A, y_A)
        np.testing.assert_almost_equal(M_A,  M_A_ref, 6)


if __name__ == '__main__':

    unittest.main()
