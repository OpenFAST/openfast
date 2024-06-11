import unittest
import numpy as np
import os
MyDir=os.path.dirname(__file__)
from pyFAST.airfoils.Polar import * 

# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
class TestPolarInterp(unittest.TestCase):

    def test_interp(self):
        # --- Interpolation of self is self
        P1=Polar(os.path.join(MyDir,'../data/FFA-W3-241-Re12M.dat'))
        P2=Polar(os.path.join(MyDir,'../data/FFA-W3-241-Re12M.dat'))
        P3= blend(P1,P2,0.5)
        np.testing.assert_equal(P3.alpha,P1.alpha)
        np.testing.assert_equal(P3.cl,P1.cl)
        np.testing.assert_equal(P3.cd,P1.cd)
        np.testing.assert_equal(P3.cm,P1.cm)
        P2.cl=P2.cl+0.3
        P2.cd=P2.cd+0.3
        P2.cm=P2.cm+0.3

        # --- Interpolation with weight 0 is first
        P3= blend(P1,P2,0.0)
        np.testing.assert_equal(P3.cl,P1.cl)
        np.testing.assert_equal(P3.cd,P1.cd)
        np.testing.assert_equal(P3.cm,P1.cm)
        # --- Interpolation with weight 1 is second
        P3= blend(P1,P2,1.0)
        np.testing.assert_equal(P3.cl,P2.cl)
        np.testing.assert_equal(P3.cd,P2.cd)
        np.testing.assert_equal(P3.cm,P2.cm)

        # --- Interpolation, weight=0.5, same alpha on both
        P1.cl=P1.cl*0+1.0
        P1.cd=P1.cl*0+1.0
        P1.cm=P1.cl*0+1.0
        P2.cl=P1.cl*0+2.0
        P2.cd=P1.cl*0+2.0
        P2.cm=P1.cl*0+2.0
        P3= blend(P1,P2,0.5)
        np.testing.assert_equal(P3.cl,P2.cl*0+1.5)
        np.testing.assert_equal(P3.cd,P2.cl*0+1.5)
        np.testing.assert_equal(P3.cm,P2.cl*0+1.5)

