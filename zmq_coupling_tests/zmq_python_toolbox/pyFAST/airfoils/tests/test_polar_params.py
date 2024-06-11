import unittest
import numpy as np
import os
MyDir=os.path.dirname(__file__)
from pyFAST.airfoils.Polar import * 

# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
class TestPolarParams(unittest.TestCase):
    def setUp(self):
        self.P235 = Polar(os.path.join(MyDir,'../data/63-235.csv'))
        self.PFFA = Polar(os.path.join(MyDir,'../data/FFA-W3-241-Re12M.dat'))
        self.PCyl = Polar(os.path.join(MyDir,'../data/Cylinder.csv'))
        self.PFFA_rad = Polar(os.path.join(MyDir,'../data/FFA-W3-241-Re12M.dat'), radians=True)

    def assertNaN(self,x):
        self.assertTrue(np.isnan(x))

    def test_alpha0(self):
        # --- Polar with one Cl value
        # non zero cl, alpha0 is nan
        self.assertNaN  (Polar(None,[100],[0.1],[],[]).alpha0())
        # cl value is 0, alpha0 is arbitrarily 0
        self.assertEqual(Polar(None,[100],[0.0],[],[]).alpha0(), 0)

        # --- Polar with one zero crossing
        P=Polar(None,[-10,10],[-0.1,0.1],[],[])
        # Alpha0 is found as long as the window holds it
        self.assertEqual(P.alpha0(window=[-50,50]),0.0)
        self.assertEqual(P.alpha0(window=[-10,10]),0.0)
        self.assertEqual(P.alpha0(window=[ -2, 2]),0.0)
        # Error when window outside, no crossing found
        self.assertRaises(Exception,P.alpha0, window=[-100,-50])

        # --- Polar with many zero crossing
        P=Polar(None,[-10,-5,0,5,10],[-0.1,0.1,-0.1,0.1,0.2],[],[])
        self.assertEqual(P.alpha0(window=[-10,-5]), -7.5)
        # Error when several zero crossing are found
        #self.assertRaises(Exception,P.alpha0, window=[-10,10])
        print('\n\n>>>>> TODO alpha0 with no zero crossing was commented!!\n\n')

        # --- Polar with constant values 
        # non zero cl, alpha0 is nan
        self.assertNaN  (Polar(None,[-10,10],[0.1,0.1],[],[]).alpha0())
        # cl is 0, alpha0 is arbitrarily 0
        self.assertEqual(Polar(None,[-10,10],[0.0,0.0],[],[]).alpha0(), 0)

        # --- Real Polars
        np.testing.assert_almost_equal(self.P235.alpha0(),-1.26, decimal=2)
        np.testing.assert_almost_equal(self.PFFA.alpha0(),-2.68, decimal=2)
        np.testing.assert_almost_equal(self.PCyl.alpha0(),0.00, decimal=2)
 

    def test_slope(self):
        alpha0 = 10
        # --- Polar with two points
        P=Polar(None,np.array([-1,1])+alpha0,[-1,1],[],[])
        sl,a0=P.cl_linear_slope(method='optim')
        np.testing.assert_almost_equal(sl,1.0)
        np.testing.assert_almost_equal(a0, alpha0)
        sl,a0=P.cl_linear_slope(method='max')
        np.testing.assert_almost_equal(sl,1.0)
        np.testing.assert_almost_equal(a0,alpha0)
        # --- Polar three points lin
        P=Polar(None,np.array([-1,0,1])+alpha0,[-1,0,1],[],[])
        sl,a0=P.cl_linear_slope(method='optim')
        np.testing.assert_almost_equal(sl,1.0)
        np.testing.assert_almost_equal(a0,alpha0)
        # --- Polar three points cst
        P=Polar(None,np.array([-1,0,2])+alpha0,[1,1,1],[],[])
        sl,a0=P.cl_linear_slope(method='optim')
        np.testing.assert_almost_equal(sl,0.0)
        np.testing.assert_almost_equal(a0,np.nan) # for constant Cl/=0, we return nan
        # --- Polar with sine shape
        P=Polar(None,[-3,-2,-1,0,1,2,3],[-1,-2,-1,0,1,0,0],[],[])
        sl,a0=P.cl_linear_slope(method='optim')
        np.testing.assert_almost_equal(sl,1.0)
        np.testing.assert_almost_equal(a0,0.0)
        sl,a0=P.cl_linear_slope(method='max')
        np.testing.assert_almost_equal(sl,1.0)
        np.testing.assert_almost_equal(a0,0.0)
        # --- Polar sine with plateaux 
        P=Polar(None,[-3,-2,-1,0,1,2,3],[-1,-2,-2,-1,0,1,1],[],[],radians=False)
        P.alpha0()
        sl,a0=P.cl_linear_slope(method='optim')
        np.testing.assert_almost_equal(sl,1.0)
        sl,a0=P.cl_linear_slope(method='max')
        np.testing.assert_almost_equal(sl,1.0)
        # --- Polar sine-line  -  Difficult to evaluate
        P=Polar(None,[-3,-2,-1,0,1,2,3],[-1,-2.1,-2,-1.1,0,1.1,1.2],[],[])
        sl,a0=P.cl_linear_slope(method='optim')
        np.testing.assert_almost_equal(sl,1.0,decimal=1)
        sl,a0=P.cl_linear_slope(method='max')
        np.testing.assert_almost_equal(sl,1.0,decimal=1)
        # --- Polar with a kink - Difficult
        P=Polar(None,[-3,-2,-1,0,1,2,3],[-1,-2,-2,-2,0,1,1],[],[])
        sl,a0=P.cl_linear_slope(method='optim')
        np.testing.assert_almost_equal(sl,1.5,decimal=1)
        sl,a0=P.cl_linear_slope(method='max')
        np.testing.assert_almost_equal(sl,2.0)
        # --- Polar step function
        P=Polar(None,np.array([-3,-2,-1,0,1,2,3])+alpha0,[-.5,-.5,-.5,-.5,.5,.5,.5],[],[])
        sl,a0=P.cl_linear_slope(method='optim')
        np.testing.assert_almost_equal(a0,10.5)
        np.testing.assert_almost_equal(sl,1.0)
        sl,a0=P.cl_linear_slope(method='max')
        np.testing.assert_almost_equal(sl,1.0)
        np.testing.assert_almost_equal(a0,10.5)
        # --- Sine
        alpha = np.linspace(-50,50,100) 
        Cl = np.sin(alpha*np.pi/180.)*180/np.pi
        P=Polar(None,alpha,Cl,[],[])
        sl,a0=P.cl_linear_slope(window=[-10,10])
        np.testing.assert_almost_equal(sl,1.0, decimal=2)
        sl,a0=P.cl_linear_slope(window=[-10,10],method='max')
        np.testing.assert_almost_equal(sl,1.0, decimal=2)
        # --- Real Polars
        P=self.PFFA
        sl,a0=P.cl_linear_slope(method='optim', radians=True) # Requesting radians outputs
        np.testing.assert_almost_equal(sl,7.091, decimal=3)
        np.testing.assert_almost_equal(a0,-0.04682, decimal=3)
        sl,a0=P.cl_linear_slope(method='optim')
        np.testing.assert_almost_equal(sl,0.123, decimal=3)
        np.testing.assert_almost_equal(a0,-2.683, decimal=3)
        sl,a0=P.cl_linear_slope(method='max')
        np.testing.assert_almost_equal(sl,0.13, decimal=3)
        # --- Real Polars, with alpha already in radians
        P=self.PFFA_rad
        sl,a0=P.cl_linear_slope(method='optim')
        np.testing.assert_almost_equal(sl,7.091, decimal=3) # This polar is already in radians
        np.testing.assert_almost_equal(a0,-0.04682, decimal=3)
        sl,a0=P.cl_linear_slope(method='optim', radians=True)
        np.testing.assert_almost_equal(sl,7.091, decimal=3) # This polar is already in radians
        np.testing.assert_almost_equal(a0,-0.04682, decimal=3)
        # --- Cylinder
        P=self.PCyl
        sl,a0=P.cl_linear_slope(method='optim')
        self.assertEqual(sl,0.0)
        # --- Real Polars
        P=self.P235
        sl,a0=P.cl_linear_slope(method='optim')
        np.testing.assert_almost_equal(sl,0.102, decimal=3)
        sl,a0=P.cl_linear_slope(method='max')
        np.testing.assert_almost_equal(sl,0.113, decimal=3)
        # --- Default method (NOTE: Might change in the future!)
        sl,a0=P.cl_linear_slope()
        np.testing.assert_almost_equal(sl,0.102, decimal=3)

        ##print(sl,a0)
        #import matplotlib.pyplot as plt
        #fig=plt.figure()
        #ax = fig.add_subplot(111)
        ## deg
        #ax.plot(P.alpha, P.cl)
        #ax.plot(WinLin,(np.array(WinLin)-a0)*sl,'o')
        #ax.plot(a0,0,'ko')
        #ax.plot(WinSearch[0],0,'ko')
        #ax.plot(WinSearch[1],0,'ko')
        #ax.plot(P.alpha,(P.alpha-a0)*sl,'--')
        #ax.set_xlim(np.array(WinLin)+[-20,+20])
        ## rad
        ##ax.plot(np.deg2rad(P.alpha), P.cl)
        ##ax.plot(WinLin,(np.array(WinLin)-a0)*sl,'o')
        ##ax.plot(np.deg2rad(P.alpha),(np.deg2rad(P.alpha)-a0)*sl,'--')
        #ax.set_ylim([-3,2.0])
        #plt.show()

    def test_linear_region(self):

        P=self.PCyl
        a_lin,cl_lin,slope,alpha0 = P.linear_region()
        np.testing.assert_almost_equal(a_lin[0],-180, decimal=1)
        np.testing.assert_almost_equal(a_lin[1], 180, decimal=1)

        P=self.P235
        a_lin,cl_lin,slope,alpha0 = P.linear_region()
        np.testing.assert_almost_equal(a_lin[0],-6.9, decimal=1)
        np.testing.assert_almost_equal(a_lin[1], 6.2, decimal=1)

        P=self.PFFA
        a_lin,cl_lin,slope,alpha0 = P.linear_region()
        np.testing.assert_almost_equal(a_lin[0],-9.8, decimal=1)
        np.testing.assert_almost_equal(a_lin[1], 9.0, decimal=1)


        #import matplotlib.pyplot as plt
        #fig=plt.figure()
        #ax = fig.add_subplot(111)
        #ax.plot(P.alpha, P.cl)
        #ax.plot(P.alpha, slope*(P.alpha-alpha0),'--')
        #ax.plot(a_lin,cl_lin)
        #ax.plot(a0,0,'ko')
        #ax.plot(WinSearch[0],0,'ko')
        #ax.plot(WinSearch[1],0,'ko')
        #ax.plot(P.alpha,(P.alpha-a0)*sl,'--')
        #ax.set_xlim(np.array(WinLin)+[-20,+20])
        #ax.set_ylim([-3,3.0])
        #plt.show()

    def test_fully_sep(self):
        # --- 63-235, for that polar
        # merging occurs at i=31 and i=120
        # at i=63 there is a singularity (f_st==1, cl_fs=cl/2)
        P=self.P235
        cl_fs,f_st=P.cl_fully_separated()
        # Below and above merging, fully sep polar is the same as original
        np.testing.assert_almost_equal(cl_fs[30] ,P.cl[30])
        np.testing.assert_almost_equal(cl_fs[121],P.cl[121])
        # Singularity at i=63
        np.testing.assert_almost_equal(cl_fs[63],P.cl[63]/2)
        np.testing.assert_almost_equal(f_st[63],1.0)
        #np.testing.assert_almost_equal(f_st[63],1.0023744)
        # Close to singularity, should be not far from cl/2
        np.testing.assert_almost_equal(cl_fs[64],P.cl[64]/2*1.004,decimal=4)
        np.testing.assert_almost_equal(cl_fs[62],P.cl[62]/2*0.998,decimal=4)
    
        # TODO TODO TODO Ensure harmony between f_st if computed with slope that is not max
        #P=Polar(os.path.join(MyDir,'../data/FFA-W3-241-Re12M.dat'),compute_params=True)
        #cl_fs,f_st0=P.cl_fully_separated()
        #f_st=(P.cl-cl_fs)/(P.cl_inv-cl_fs);

        #import matplotlib.pyplot as plt
        #fig=plt.figure()
        #ax = fig.add_subplot(111)
        #ax.plot(P.alpha  , f_st0 )
        #ax.plot(P.alpha  , f_st )
        #plt.show()

        #import matplotlib.pyplot as plt
        #fig=plt.figure()
        #ax = fig.add_subplot(111)
        #ax.plot(P.alpha,P.cl,label='Cl')
        #ax.plot(P.alpha,cl_fs,'--',label='Cl_fs')
        #ax.plot(P.alpha,f_st,label='f_st')
        #plt.xlim([-50,50])
        #plt.ylim([-3,3])
        #plt.legend()
        #plt.show()
        #print(f_st)
# 
#         P=Polar(os.path.join(MyDir,'../data/Cylinder.dat'))
#         sl,offset=P.cl_linear_slope()
# 
#         plt.show()

#         P=Polar(os.path.join(MyDir,'../data/FFA-W3-241-Re12M.dat'))
#         P=Polar(os.path.join(MyDir,'../data/63-235.csv'))
#         P=Polar(os.path.join(MyDir,'../data/Cylinder.csv'))
#         P=Polar(os.path.join(MyDir,'../data/AD_3-63-224_mod.csv'))
#         P=Polar(os.path.join(MyDir,'../data/AD_4-63-218_mod.csv'))
#         P=Polar(os.path.join(MyDir,'../data/AD_5_63-214_mod.csv'))


if __name__ == '__main__':
    unittest.main()
