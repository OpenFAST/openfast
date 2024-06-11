import unittest
import os
import numpy as np
from pyFAST.converters.beam import ComputeStiffnessProps, ComputeInertiaProps, TransformCrossSectionMatrix
from pyFAST.converters.beam import MM, KK
from pyFAST.converters.beam import K66toPropsDecoupled, M66toPropsDecoupled


class Test(unittest.TestCase):

    def test_BD_H2(self):
        
        K = np.array([[ 3.22114734e+09, 3.25671541e+07, 0.00000000e+00, 0.00000000e+00,  0.00000000e+00, -4.99668423e+07],
              [ 3.25671541e+07, 2.43648638e+09, 0.00000000e+00, 0.00000000e+00,  0.00000000e+00,  3.71784304e+07],
              [ 0.00000000e+00, 0.00000000e+00, 2.44325680e+10, 1.90703565e+08,  1.78086138e+09,  0.00000000e+00],
              [ 0.00000000e+00, 0.00000000e+00, 1.90703565e+08, 4.84065153e+10,  4.60061848e+09,  0.00000000e+00],
              [ 0.00000000e+00, 0.00000000e+00, 1.78086138e+09, 4.60061848e+09,  5.23688767e+10,  0.00000000e+00],
              [-4.99668423e+07, 3.71784304e+07, 0.00000000e+00, 0.00000000e+00,  0.00000000e+00,  2.33148898e+10]])

        I = np.array([[1439.88989372 ,    0.         ,    0.         ,    0.         ,    0.         ,   -5.92729626],
                    [   0.         , 1439.88989372 ,    0.         ,    0.         ,    0.         ,   -68.63720816],
                    [   0.         ,    0.         , 1439.88989372 ,    5.92729626 ,   68.63720816 ,    0.        ],
                    [   0.         ,    0.         ,    5.92729626 , 2594.38760678 ,   47.71380715 ,    0.        ],
                    [   0.         ,    0.         ,   68.63720816 ,   47.71380715 , 3434.04030538 ,    0.        ],
                    [  -5.92729626 ,  -68.63720816 ,    0.         ,    0.         ,    0.         ,   6028.42791216]])
        
        stiff = ComputeStiffnessProps()
        inertia = ComputeInertiaProps()
        transform = TransformCrossSectionMatrix()
        xs , ys = stiff.ComputeShearCenter(K)
        xt , yt = stiff.ComputeTensionCenter(K)
        xm , ym = inertia.ComputeMassCenter(I)

        # Approach BECAS
        # Translate to tension (elastic) center
        Kel = transform.CrossSectionRotoTranslationMatrix(K, xt, yt, 0.)
        # Find delta
        DeltaBecas = stiff.OrientationPrincipalAxesBecas(Kel)
        # Rotate by delta
        Kel_DeltaBecas = transform.CrossSectionRotoTranslationMatrix(Kel, 0., 0., DeltaBecas)

        # Approach ANBA4
        # Translate to tension (elastic) center
        Kdec = stiff.DecoupleStiffness(K)
        # Find delta
        DeltaANBA4 = stiff.PrincipalAxesRotationAngle(Kdec)
        # Rotate by delta
        Kel_DeltaANBA4 = transform.CrossSectionRotoTranslationMatrix(Kdec, 0., 0., -DeltaANBA4)

        #print(stiff.ComputeTensionCenter(Kel_DeltaBecas))
        #print(stiff.ComputeTensionCenter(Kel_DeltaANBA4))
        #print(stiff.OrientationPrincipalAxesBecas(Kel_DeltaBecas))
        #print(stiff.OrientationPrincipalAxesBecas(Kel_DeltaANBA4))
        #print(stiff.PrincipalAxesRotationAngle(Kel_DeltaBecas))
        #print(stiff.PrincipalAxesRotationAngle(Kel_DeltaANBA4))

        np.testing.assert_almost_equal(xs,  0.015468467117843322, 6)
        np.testing.assert_almost_equal(ys,  0.015668518364738194, 6)
        np.testing.assert_almost_equal(xt, -0.07288883346195946, 6)
        np.testing.assert_almost_equal(yt,  0.0078053017185913485, 6)
        np.testing.assert_almost_equal(xm, -0.04766837274110846, 6)
        np.testing.assert_almost_equal(ym,  0.004116492716458095, 6)
        np.testing.assert_almost_equal(DeltaBecas, -0.5874557755802033, 6)
        np.testing.assert_almost_equal(DeltaANBA4,  0.5874557755802033, 6)

        # --- Using decoupled methods
        EA2, EIx2, EIy2, kxsGA2, kysGA2, GKt2, x_C2, y_C2, x_S2, y_S2, theta_p2, theta_s2 = K66toPropsDecoupled(K)
        m2, Ixi2, Iyi2, Ip2, x_G2, y_G2, theta_i2 = M66toPropsDecoupled(I)

        np.testing.assert_almost_equal(x_C2, xt, 6)
        np.testing.assert_almost_equal(y_C2, yt, 6)
        np.testing.assert_almost_equal(x_S2, xs, 6)
        np.testing.assert_almost_equal(y_S2, ys, 6)
        np.testing.assert_almost_equal(DeltaANBA4, -theta_p2, 6)

        


    def test_decoupleMat(self):
        # Test the detection of main points (center of mass, tension and shear center)
        # for a simple (decoupled) beam cross section matrix in BeamDyn coordinates

        # Main points
        x_G     =  0.175
        y_G     = -0.1
        x_C     =  0.164
        y_C     = -0.13
        x_S     =  0.18
        y_S     = -0.136
        # Inertia param
        m       =  262.
        Ixi     =  163.
        Iyi     =  54.
        I_p     =  217.
        theta_i = -0.1
        # Stiffness param
        EA      = 4.1e10
        EIxp    = 2.1e9
        EIyp    = 0.9e9
        GKt     = 2.3e8
        GA      = 3.75e9
        kxs     = 0.5
        kys     = 0.5
        theta_p = -0.1
        theta_s = -0.1
        # Mass and stiffness matrix for simple beam
        M =  MM(m, Ixi, Iyi, I_p, x_G, y_G, theta_i) # NOTE: theta_i in rad
        K =  KK(EA, EIxp, EIyp, GKt, GA, kxs, kys, x_C, y_C, theta_p, x_S, y_S, theta_s) # Note theta_p/s in rad

        # --- Using general methods
        stiff = ComputeStiffnessProps()
        inertia = ComputeInertiaProps()
        xs , ys = stiff.ComputeShearCenter(K)
        xt , yt = stiff.ComputeTensionCenter(K)
        xm , ym = inertia.ComputeMassCenter(M)

        np.testing.assert_almost_equal(xs,  x_S, 6)
        np.testing.assert_almost_equal(ys,  y_S, 6)
        np.testing.assert_almost_equal(xm,  x_G, 6)
        np.testing.assert_almost_equal(ym,  y_G, 6)
        np.testing.assert_almost_equal(xt,  x_C, 6)
        np.testing.assert_almost_equal(yt,  y_C, 6)

        # --- Using decoupled methods
        EA2, EIx2, EIy2, kxsGA2, kysGA2, GKt2, x_C2, y_C2, x_S2, y_S2, theta_p2, theta_s2 = K66toPropsDecoupled(K)
        m2, Ixi2, Iyi2, Ip2, x_G2, y_G2, theta_i2 = M66toPropsDecoupled(M)

        np.testing.assert_almost_equal(x_S2    ,x_S     ,6)
        np.testing.assert_almost_equal(y_S2    ,y_S     ,6)
        np.testing.assert_almost_equal(x_G2    ,x_G     ,6)
        np.testing.assert_almost_equal(y_G2    ,y_G     ,6)
        np.testing.assert_almost_equal(x_C2    ,x_C     ,6)
        np.testing.assert_almost_equal(y_C2    ,y_C     ,6)
        np.testing.assert_almost_equal(EA2     ,EA      ,6)
        np.testing.assert_almost_equal(EIx2/1e9,EIxp/1e9,6)
        np.testing.assert_almost_equal(EIy2/1e9,EIyp/1e9,6)
        np.testing.assert_almost_equal(kxsGA2  ,kxs*GA  ,6)
        np.testing.assert_almost_equal(kysGA2  ,kys*GA  ,6)
        np.testing.assert_almost_equal(GKt2/1e9,GKt/1e9 ,6)
        np.testing.assert_almost_equal(theta_p2,theta_p ,6)
        #np.testing.assert_almost_equal(theta_s2,theta_s ,6) # ....
        np.testing.assert_almost_equal(m2      ,m       ,6)
        np.testing.assert_almost_equal(Ixi2    ,Ixi     ,6)
        np.testing.assert_almost_equal(Iyi2    ,Iyi     ,6)
        np.testing.assert_almost_equal(Ip2     ,I_p     ,6)




if __name__ == '__main__':

    unittest.main()
