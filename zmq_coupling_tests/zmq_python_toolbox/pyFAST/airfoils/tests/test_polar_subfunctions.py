
import unittest
import numpy as np
import os
import matplotlib.pyplot as plt
MyDir=os.path.dirname(__file__)
from pyFAST.airfoils.Polar import _intersections

# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{


class TestPolarSubFunctions(unittest.TestCase):


    def test_intersection1(self):
        a, b = 1, 2
        phi = np.linspace(3, 10, 100)
        x1 = a*phi - b*np.sin(phi)
        y1 = a - b*np.cos(phi)

        x2=phi
        y2=np.sin(phi)+2
        x,y=_intersections(x1,y1,x2,y2, plot=False)

        np.testing.assert_almost_equal(x,[6.10766,8.3648], decimal=4)
        np.testing.assert_almost_equal(y,[1.825397,2.87209], decimal=4)

    def test_intersection2(self):

        alpha  = np.linspace(-180,180,121)
        Cl     = alpha*np.pi/180 * 2*np.pi
        Cn     = Cl*np.cos(alpha*np.pi/180)
        Cn_lin =  0.10876153755016482* (alpha) #*np.pi/180)
        Cn_f   = Cn_lin * ((1 + np.sqrt(0.7)) / 2) ** 2

        x1 = alpha
        y1 = Cn

        x2 = alpha
        y2 = Cn_f

        x,y=_intersections(x1,y1,x2,y2, plot=False, verbose=False)

        np.testing.assert_almost_equal(x,[-33.2116, 0, 33.2116 ], decimal=3)
        np.testing.assert_almost_equal(y,[-3.046  , 0 , 3.04623], decimal=3)


if __name__ == '__main__':
    unittest.main()
