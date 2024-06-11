import unittest
import numpy as np
import os
MyDir=os.path.dirname(__file__)
from pyFAST.airfoils.Polar import * 

# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
class TestPolarManip(unittest.TestCase):
    def assertNaN(self,x):
        self.assertTrue(np.isnan(x))

    def test_read(self):
        P=Polar(os.path.join(MyDir,'../data/FFA-W3-241-Re12M.dat'))
        self.assertEqual(P.alpha[-1],180)
        self.assertEqual(P.cl[-1],0)

        P=Polar(os.path.join(MyDir,'../data/Cylinder.dat'))
        self.assertEqual(P.cl.size,3)

if __name__ == '__main__':
    unittest.main()
