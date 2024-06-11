"""
Set of tools for stiffness and mass matrix of a beam section
"""

import numpy as np
from numpy import cos, sin

# --------------------------------------------------------------------------------}
# --- Functions for 6x6 rigid body mass matrices 
# --------------------------------------------------------------------------------{
def identifyRigidBodyMM(MM):
    """ 
    Based on a 6x6 mass matrix at a reference point:
     - Identify the position of the center of mass
     - Compute the inertia at the center of mass
    """
    mass = MM[0,0]
    # Using average of two coeffs to get estimate of COG
    xCM = 0.5*( MM[1,5]-MM[2,4])/mass
    zCM = 0.5*( MM[0,4]-MM[1,3])/mass
    yCM = 0.5*(-MM[0,5]+MM[2,3])/mass
    # Destance from refopint to COG
    Ref2COG=np.array((xCM,yCM,zCM))
    # Inertia at ref oint
    J_P = MM[3:6,3:6].copy()
    # Inertia at COG
    J_G = translateInertiaMatrixToCOG(J_P, mass, r_PG=Ref2COG ) 
    return mass, J_G, Ref2COG

def translateInertiaMatrixToCOG(I_P, Mass, r_PG): 
    """ Transform inertia matrix with respect to point P to the inertia matrix with respect to the COG
    NOTE: the vectors and the inertia matrix needs to be expressed in the same coordinate system.
    
    INPUTS:
      I_G  : Inertia matrix 3x3 with respect to COG
      Mass : Mass of the body
      r_PG: vector from P to COG 
    """
    I_G = I_P + Mass * np.dot(skew(r_PG), skew(r_PG))
    return I_G

def translateInertiaMatrixFromCOG(I_G, Mass, r_GP): 
    """
    Transform inertia matrix with respect to COG to the inertia matrix with respect to point P
    NOTE: the vectors and the inertia matrix needs to be expressed in the same coordinate system.
    INPUTS:
       I_G  : Inertia matrix 3x3 with respect to COG
       Mass : Mass of the body
       r_GP: vector from COG of the body to point P
    """
    I_P = I_G - Mass * np.dot(skew(r_GP),skew(r_GP))
    return I_P

def rigidBodyMassMatrixAtP(m=None, J_G=None, Ref2COG=None):
    """ 
    Rigid body mass matrix (6x6) at a given reference point: 
      the center of gravity (if Ref2COG is None) 


    INPUTS:
     - m/tip: (scalar) body mass 
                     default: None, no mass
     - J_G: (3-vector or 3x3 matrix), diagonal coefficients or full inertia matrix
                     with respect to COG of body! 
                     The inertia is transferred to the reference point if Ref2COG is not None
                     default: None 
     - Ref2COG: (3-vector) x,y,z position of center of gravity (COG) with respect to a reference point
                     default: None, at first/last node.
    OUTPUTS:
      - M66 (6x6) : rigid body mass matrix at COG or given point 
    """
    # Default values
    if m is None: m=0
    if Ref2COG is None: Ref2COG=(0,0,0)
    if J_G is None: J_G=np.zeros((3,3))
    if len(J_G.flatten()==3): J_G = np.eye(3).dot(J_G)

    M66 = np.zeros((6,6))
    x,y,z = Ref2COG
    Jxx,Jxy,Jxz = J_G[0,:]
    _  ,Jyy,Jyz = J_G[1,:]
    _  ,_  ,Jzz = J_G[2,:]
    M66[0, :] =[   m     ,   0     ,   0     ,   0                 ,  z*m                , -y*m                 ]
    M66[1, :] =[   0     ,   m     ,   0     , -z*m                ,   0                 ,  x*m                 ]
    M66[2, :] =[   0     ,   0     ,   m     ,  y*m                , -x*m                ,   0                  ]
    M66[3, :] =[   0     , -z*m    ,  y*m    , Jxx + m*(y**2+z**2) , Jxy - m*x*y         , Jxz  - m*x*z         ]
    M66[4, :] =[  z*m    ,   0     , -x*m    , Jxy - m*x*y         , Jyy + m*(x**2+z**2) , Jyz  - m*y*z         ]
    M66[5, :] =[ -y*m    , x*m     ,   0     , Jxz - m*x*z         , Jyz - m*y*z         , Jzz  + m*(x**2+y**2) ]
    return M66

def skew(x):
    x=np.asarray(x).ravel()
    """ Returns the skew symmetric matrix M, such that: cross(x,v) = M v """
    return np.array([[0, -x[2], x[1]],[x[2],0,-x[0]],[-x[1],x[0],0]])


# --------------------------------------------------------------------------------}
# --- Bauchau 
# --------------------------------------------------------------------------------{
def K_sheartorsion_xbeam(J,K22,K23,K33,x2,x3):
    """ Returns Shear-torsion stiffness matrix.  See Eq.(13) of DyMore manual """
    return np.array( [
        [J + K22*x3**2 + 2*K23*x2*x3 + K33*x2**2, -K22*x3 - K23*x2, K23*x3 + K33*x2],
        [-K22*x3 - K23*x2, K22, -K23],
        [K23*x3 + K33*x2, -K23, K33]])

def K_axialbending_xbeam(S,H22,H23,H33,x2,x3):
    """ Returns Axial-Bending stiffness matrix. See Eq.(20) of DyMore manual """
    return np.array([
        [S, S*x3, -S*x2],
        [S*x3, H22 + S*x3**2, -H23 - S*x2*x3],
        [-S*x2, -H23 - S*x2*x3, H33 + S*x2**2]])

# --------------------------------------------------------------------------------}
# --- BeamDyn 
# --------------------------------------------------------------------------------{
def K_axialbending(EA, EI_x, EI_y, x_C=0, y_C=0, theta_p=0):
    """
    Axial bending problem. See KK for notations.
    """
    H_xx = EI_x*cos(theta_p)**2 + EI_y*sin(theta_p)**2 
    H_yy = EI_x*sin(theta_p)**2 + EI_y*cos(theta_p)**2
    H_xy = (EI_y-EI_x)*sin(theta_p)*cos(theta_p)
    return np.array([
        [EA      , EA*y_C             , -EA*x_C            ] ,
        [EA*y_C  , H_xx + EA*y_C**2   , -H_xy - EA*x_C*y_C ] ,
        [-EA*x_C , -H_xy - EA*x_C*y_C , H_yy + EA*x_C**2   ] 
        ])

def K_sheartorsion(GKt, GA, kxs, kys, x_S=0, y_S=0, theta_s=0):
    """
    Shear torsion problem. See KK for notations.
    """
    K_xx = GA * ( kxs*cos(theta_s)**2 + kys*sin(theta_s)**2   ) 
    K_yy = GA * ( kxs*sin(theta_s)**2 + kys*cos(theta_s)**2   )
    K_xy = GA * ( (kys-kxs)*sin(theta_s)*cos(theta_s)         )
    return np.array([
        [K_xx                 , -K_xy               , -K_xx*y_S - K_xy*x_S                             ] ,
        [-K_xy                , K_yy                , K_xy*y_S + K_yy*x_S                              ] ,
        [-K_xx*y_S - K_xy*x_S , K_xy*y_S + K_yy*x_S , GKt + K_xx*y_S**2 + 2*K_xy*x_S*y_S + K_yy*x_S**2 ]
        ])


# --------------------------------------------------------------------------------}
# --- Beam stiffness and mass matrix as function of "beam" properties
# --------------------------------------------------------------------------------{
def KK(EA, EI_x, EI_y, GKt, GA, kxs, kys, x_C=0, y_C=0, theta_p=0, x_S=0, y_S=0, theta_s=0):
    """ 
    Returns 6x6 stiffness matrix at the cross section origin O, based on inputs at centroid and shear center.
    INPUTS:
        - EA, EI_x, EI_y: diagonal terms for the axial bending expressed at the centroid and in the principal axis frame
        - GKt, GA*kxs, GA*kys: diagonal terms for the shear/torsion expressed at the shear center and in the princial shear direction frame
        - kxs, kys: dimensionless shear parameters
        - x_C, y_C: coordinates of the centroid (elastic center/ neutral axis), expressed FROM the origin of the cross section O
        - x_S, y_S:       "            shear center            "                  "                                             
        - theta_p : angle (around z) FROM the reference axes to the principal axes [rad]
        - theta_s :       "            "             "              principal shear axes [rad]
    """
    H_xx = EI_x*cos(theta_p)**2 + EI_y*sin(theta_p)**2 
    H_yy = EI_x*sin(theta_p)**2 + EI_y*cos(theta_p)**2
    H_xy = (EI_y-EI_x)*sin(theta_p)*cos(theta_p)
    K_xx = GA * ( kxs*cos(theta_s)**2 + kys*sin(theta_s)**2   ) 
    K_yy = GA * ( kxs*sin(theta_s)**2 + kys*cos(theta_s)**2   )
    K_xy = GA * ( (kys-kxs)*sin(theta_s)*cos(theta_s)         )
    return np.array([
        [K_xx                 , -K_xy               , 0*EA    , 0*EA               , 0*EA               , -K_xx*y_S - K_xy*x_S                             ]    , 
        [-K_xy                , K_yy                , 0*EA    , 0*EA               , 0*EA               , K_xy*y_S + K_yy*x_S                              ]    , 
        [0*EA                 , 0*EA                , EA      , EA*y_C             , -EA*x_C            , 0*EA                                                ] , 
        [0*EA                 , 0*EA                , EA*y_C  , H_xx + EA*y_C**2   , -H_xy - EA*x_C*y_C , 0*EA                                                ] , 
        [0*EA                 , 0*EA                , -EA*x_C , -H_xy - EA*x_C*y_C , H_yy + EA*x_C**2   , 0*EA                                                ] , 
        [-K_xx*y_S - K_xy*x_S , K_xy*y_S + K_yy*x_S , 0*EA    , 0*EA               , 0*EA               , GKt + K_xx*y_S**2 + 2*K_xy*x_S*y_S + K_yy*x_S**2 ]
        ])

def MM(m,I_x,I_y,I_p,x_G=0,y_G=0,theta_i=0):
    """ 
    Returns the mass matrix at a given point O and with respect to given orientation axes based 
    on the values at the center of gravity and in the inertia axis frame.
    The convention is such that:
      - x_G,y_G      : the distaances FROM point O to point G
      - theta_i      : angle (around z) FROM the reference axes to the inertial axes
      - I_x, I_y, I_p: "diagonal" inertias for the body expressed in the inertial frame and at point G
    """
    Ixx = I_x*cos(theta_i)**2 + I_y*sin(theta_i)**2 
    Iyy = I_x*sin(theta_i)**2 + I_y*cos(theta_i)**2
    Ixy = (I_y-I_x)*sin(theta_i)*cos(theta_i)

    return np.array([
        [m      , 0*m     , 0*m      , 0*m                , 0*m                , -m*y_G]                    , 
        [0*m      , m     , 0*m      , 0*m                , 0*m                , m*x_G]                     , 
        [0*m      , 0*m     , m      , m*y_G            , -m*x_G           , 0*m]                         , 
        [0*m      , 0*m     , m*y_G  , Ixx + m*y_G**2   , -Ixy - m*x_G*y_G , 0*m]                         , 
        [0*m      , 0*m     , -m*x_G , -Ixy - m*x_G*y_G , Iyy + m*x_G**2   , 0*m]                         , 
        [-m*y_G , m*x_G , 0*m      , 0*m                , 0*m                , I_p + m*x_G**2 + m*y_G**2]
        ])

# --------------------------------------------------------------------------------}
# --- Tools to manipulate 6x6 stiffness and mass matrix
# --------------------------------------------------------------------------------{
class TransformCrossSectionMatrix(object):

    def CrossSectionTranslationMatrix(self, x, y):
        T = np.eye(6)
        T[0,5] = y
        T[1,5] = -x
        T[2,3] = -y
        T[2,4] = x
        return T

    def CrossSectionRotationMatrix(self, alpha):
        c=np.cos(alpha)
        s=np.sin(alpha)
        R1=[[c,s,0],
            [-s,c,0],
            [0,0,1]]
        R=np.vstack((np.hstack((R1, np.zeros((3,3)))),
           np.hstack((np.zeros((3,3)), R1))))
        return R

    def CrossSectionRotoTranslationMatrix(self, M1, x, y, alpha):
        # Translation
        T = self.CrossSectionTranslationMatrix(x, y)
        M2 = T.T @ M1 @ T 
        # Rotation 
        R = self.CrossSectionRotationMatrix(alpha)
        M3 = R @ M2 @ R.T
        return M3

class ComputeStiffnessProps(object):

    def ComputeShearCenter(self, K):   # shear center equiv. to elastic axes
        """ 
        Shear center location for a 6x6 cross section matrix of a beam oriented along "z"
        NOTE: works for BeamDyn or Hawc2 convention
        """
        K1 = np.array([[K[i, j] for j in range(3)] for i in range(3)])
        K3 = np.array([[K[i, j+3] for j in range(3)] for i in range(3)])
        Y = np.linalg.solve(K1, -K3)
        return [-Y[1,2], Y[0,2]]

    def ComputeTensionCenter(self, K):  # tension center equiv. to neutral axes
        """ 
        Tension center (also called elastic center) location for a 6x6 cross section matrix of a beam oriented along "z"
        NOTE: works for BeamDyn or Hawc2 convention
        """
        K1 = np.array([[K[i, j] for j in range(3)] for i in range(3)])
        K3 = np.array([[K[i, j+3] for j in range(3)] for i in range(3)])
        Y = np.linalg.solve(K1, -K3)
        return [Y[2,1], -Y[2,0]]

    def OrientationPrincipalAxesBecas(self, K):
        ksub=K[3:5,3:5]
        [ val, mod ] = np.linalg.eig(ksub)
        val = np.sort(np.diag(val))
        ind = np.argsort(np.diag(val))
        mod = mod[:,ind]
        Delta = np.arctan(mod[1,0]/mod[0,0])
        return Delta

    def DecoupleStiffness(self, K):
        K1 = np.array([[K[i, j] for j in range(3)] for i in range(3)])
        K3 = np.array([[K[i, j+3] for j in range(3)] for i in range(3)])
        Y = np.linalg.solve(K1, -K3)
        I3 = np.eye(3)
        Z3 = np.zeros((3,3))
        TL = np.block([[I3, Z3], [Y.T, I3]])
        TR = np.block([[I3, Y], [Z3, I3]])
        return TL @ K @ TR

    def PrincipalAxesRotationAngle(self, decoupledK):
        K1 = np.array([[decoupledK[i, j] for j in range(3)] for i in range(3)])
        K3 = np.array([[decoupledK[i+3, j+3] for j in range(3)] for i in range(3)])
        (w1, v1) = np.linalg.eig(K1)
        (w3, v3) = np.linalg.eig(K3)
        if np.abs(v3[0,0]) < np.abs(v3[0,1]):
            angle = np.arccos(v3[0,0])
        else:
            angle = -np.arcsin(v3[0,1])
        return angle

class ComputeInertiaProps(object):     
    
    def ComputeMassCenter(self, M):
        """ 
        Mass center location for a 6x6 cross section matrix of a beam oriented along "z"
        NOTE: works for BeamDyn or Hawc2 convention
        """
        M1 = np.array([[M[i, j] for j in range(3)] for i in range(3)])
        M3 = np.array([[M[i, j+3] for j in range(3)] for i in range(3)])
        Y = np.linalg.solve(M1, -M3)
        return [Y[2,1], -Y[2,0]]

    def Translate(self, M, x, y):
        return TranslateSectionMassMatrix(M, x ,y)


def TranslateSectionMassMatrix(M, x, y):
    """ 
    Translate a 6x6 section mass matrix to a different point
    """
    # First identify main properties (mass, inertia, location of center of mass from previous ref point)
    mass, J_G, Ref2COG = identifyRigidBodyMM(M)
    # New COG location from new (x,y) ref point
    Ref2COG[0] -= x
    Ref2COG[1] -= y
    Ref2COG[2]  = 0 # for sanity
    # Compute mass matrix (NOTE: using 6x6 of rigid body here, could use MM but this interface is simpler)
    M_new =  rigidBodyMassMatrixAtP(mass, J_G, Ref2COG)
    return M_new


# --------------------------------------------------------------------------------}
# --- Functions to extract "beam" properties from 6x6 matrices
# --------------------------------------------------------------------------------{
def solvexytheta(Hxx,Hyy,Hxy):
    """ 
     Solve for a system of three unknown, used to get:
       - EIy, EIx and thetap given Hxx,Hyy,Hxy
       - kxs*GA, kys*GA and thetas given Kxx,Kyy,Kxy
       - I_x, I_y and theta_is given Ixx,Iyy,Ixy
    """
    from scipy.optimize import fsolve
    def residual(x):
        EI_x, EI_y, theta_p =x
        res=np.array([
            Hxx - EI_x*np.cos(theta_p)**2 - EI_y*np.sin(theta_p)**2 ,
            Hyy - EI_x*np.sin(theta_p)**2 - EI_y*np.cos(theta_p)**2,
            Hxy - (EI_y-EI_x)*np.sin(theta_p)*np.cos(theta_p)]
                ).astype(float)
        return res
    x0 = [Hxx,Hyy,0]
    x_opt   = fsolve(residual, x0)
    EI_x,EI_y,theta_p = x_opt
    theta_p = np.mod(theta_p,2*np.pi)
    return EI_x, EI_y, theta_p


def M66toPropsDecoupled(M, convention='BeamDyn'):
    """ 
    Convert mass properties of a 6x6 section to beam properties
    This assumes that the axial and bending loads are decoupled.
    Has been tested for BeamDyn coordinate system

    INPUTS:
     - M : 6x6 array of mass elements. Each element may be an array (e.g. for all spanwise values)
     OUTPUTS:
      - m: section mass
      - Ixx, Iyy, Ixy: area moment of inertia
      - x_G, y_G
    """
    M11=M[0,0]
    M44=M[3,3]
    M55=M[4,4]
    M66=M[5,5]
    M16=M[0,5]
    M26=M[1,5]
    M45=M[3,4]

    Mlist = flat2matlist(M)
    inertia = ComputeInertiaProps()

    m=M11
    if convention=='BeamDyn':
        if np.all(np.abs(m)<1e-16):
            Ixi = Iyi = Ipp = x_G = y_G = theta_i = m*0
            return m, Ixi, Iyi, Ipp, x_G, y_G, theta_i
        y_G= -M16/m
        x_G=  M26/m
        # sanity
        np.testing.assert_array_almost_equal([M[0,3],M[0,4]],[0*M[0,3],0*M[0,3]])
        np.testing.assert_array_almost_equal([M[1,3],M[1,4]],[0*M[0,3],0*M[0,3]])
        np.testing.assert_array_almost_equal([M[3,5],M[4,5]],[0*M[0,3],0*M[0,3]])
        
        Ixx =  M44-m*y_G**2
        Iyy =  M55-m*x_G**2
        Ixy = -M45-m*x_G*y_G
        Ipp =  M66-m*(x_G**2 + y_G**2)

        if np.all(np.abs(Ixy)<1e-16):
            # NOTE: Assumes theta_i ==0
            #print('>>> Assume theta_i 0')
            Ixi     = Ixx
            Iyi     = Iyy
            theta_i = Ixx*0
        else:
            if len(Mlist)==1:
                Ixi,Iyi,theta_i = solvexytheta(Ixx,Iyy,Ixy)
            else:
                #print('>>> Minimize theta_i')
                Ixi = np.zeros(Ixx.shape)
                Iyi = np.zeros(Ixx.shape)
                theta_i = np.zeros(Ixx.shape)
                for i, (hxx, hyy, hxy) in enumerate(zip(Ixx,Iyy,Ixy)):
                    Ixi[i],Iyi[i],theta_i[i] = solvexytheta(hxx,hyy,hxy)

                    MM2= MM(m[i],Ixi[i],Iyi[i],Ipp[i],x_G[i],y_G[i],theta_i[i])

                    np.testing.assert_allclose(MM2[3,3], M[3,3][i], rtol=1e-3)
                    np.testing.assert_allclose(MM2[4,4], M[4,4][i], rtol=1e-3)
                    np.testing.assert_allclose(MM2[5,5], M[5,5][i], rtol=1e-3)
                    np.testing.assert_allclose(MM2[3,4], M[3,4][i], rtol=1e-3)

        np.testing.assert_array_almost_equal(Ipp, Ixx+Iyy, 2)
        np.testing.assert_array_almost_equal(Ipp, Ixi+Iyi, 2)
         
    else:
        raise NotImplementedError()

    return m, Ixi, Iyi, Ipp, x_G, y_G, theta_i



def K66toPropsDecoupled(K, theta_p_in=None, convention='BeamDyn'):
    """ 
    Convert stiffness properties of a 6x6 section to beam properties
    This assumes that the axial and bending loads are decoupled.
    Has been tested for BeamDyn coordinate system

    INPUTS:
     - K : 6x6 array of stiffness elements. Each element may be an array (e.g. for all spanwise values)
    INPUTS OPTIONAL:
     - theta_p_in : angle from section to principal axis [rad], positive around z
     - convention : to change coordinate systems in the future
    OUTPUTS:
     - EA, EIx, EIy: axial and bending stiffnesses
     - kxGA, kyGA, GKt: shear and torsional stiffness
     - xC,yC : centroid
     - xS,yS : shear center
     - theta_p, theta_s: angle to principal axes and shear axes [rad]
    """
    K11=K[0,0]
    K22=K[1,1]
    K33=K[2,2]
    K44=K[3,3]
    K55=K[4,4]
    K66=K[5,5]

    K12=K[0,1]
    K16=K[0,5]
    K26=K[1,5]
    K34=K[2,3]
    K35=K[2,4]
    K45=K[3,4]

    Klist = flat2matlist(K)
    stiff = ComputeStiffnessProps()

    if convention=='BeamDyn':
        # --------------------------------------------------------------------------------}
        # --- Axial/bending problem
        # --------------------------------------------------------------------------------{
        # Find: EA, EI, position of axial strain centroid (tension center, elastic center), principal axes
        EA =  K33
        # xC yC - method 1
        yC =  K34/EA
        xC = -K35/EA
        # xC yC - method 2, more general
        if len(Klist)==1:
            xC2, yC2 = stiff.ComputeTensionCenter(Klist[0])   # tension center
        else:
            xC2   = np.zeros(xC.shape)
            yC2   = np.zeros(yC.shape)
            for i,Kloc in enumerate(Klist):
                xs, ys = stiff.ComputeTensionCenter(Kloc)   # tension center
                xC2[i] = xs
                yC2[i] = ys
        # Sanity checks, comment in future
        np.testing.assert_allclose(xC2, xC, rtol=1e-3)
        np.testing.assert_allclose(yC2, yC, rtol=1e-3)

        # --- Find EI and theta_p (principal axes)
        Hxx=  K44-EA*yC**2
        Hyy=  K55-EA*xC**2 
        Hxy= -K45-EA*xC*yC
        if theta_p_in is not None:
            theta_p=theta_p_in
            print('>>> theta_p given')
            C2=np.cos(theta_p)**2
            S2=np.sin(theta_p)**2
            C4=np.cos(theta_p)**4
            S4=np.sin(theta_p)**4
            EIxp = (Hxx*C2 - Hyy*S2)/(C4-S4)
            EIyp = (Hxx*S2 - Hyy*C2)/(S4-C4)
            Hxyb = (EIyp-EIxp)*np.sin(theta_p)*np.cos(theta_p)

            bNZ=np.logical_and(Hxy!=0, Hxyb!=0)
            np.testing.assert_allclose(Hxy[bNZ], Hxyb[bNZ], rtol=1e-3)
            np.testing.assert_allclose(EIxp+EIyp, Hxx+Hyy, rtol=1e-3)

        else:
            if np.all(np.abs(Hxy)<1e-16):
                #print('>>>> assume theta_p=0')
                # NOTE: Assumes theta_p ==0
                EIxp = Hxx
                EIyp = Hyy
                theta_p=0*EA
            else:
                #print('>>> Minimization for theta_p')
                if len(Klist)==1:
                    EIxp,EIyp,theta_p = solvexytheta(Hxx,Hyy,Hxy)
                    theta_p=np.asarray(theta_p)
                    if theta_p>np.pi:
                        theta_p-=2*np.pi
                else:
                    EIxp= np.zeros(Hxx.shape)
                    EIyp= np.zeros(Hxx.shape)
                    theta_p = np.zeros(Hxx.shape)
                    for i, (hxx, hyy, hxy) in enumerate(zip(Hxx,Hyy,Hxy)):
                        EIxp[i],EIyp[i],theta_p[i] = solvexytheta(hxx,hyy,hxy)

                    theta_p[theta_p>np.pi]=theta_p[theta_p>np.pi]-2*np.pi

        # --------------------------------------------------------------------------------}
        # ---Torsion/shear problem
        # --------------------------------------------------------------------------------{
        # Find: Torsion, shear terms, shear center
        Kxx =  K11
        Kxy = -K12
        Kyy =  K22
        # Method 1
        yS  = (Kyy*K16+Kxy*K26)/(-Kyy*Kxx + Kxy**2)
        xS  = (Kxy*K16+Kxx*K26)/( Kyy*Kxx - Kxy**2)
        # Method 2, more general
        if len(Klist)==1:
            xS2, yS2 = stiff.ComputeShearCenter(Klist[0])  # shear center
        else:
            xS2 = np.zeros(xC.shape)
            yS2 = np.zeros(yC.shape)
            for i,Kloc in enumerate(Klist):
                xS2[i], yS2[i] = stiff.ComputeShearCenter(Kloc)    # shear center
        # Sanity check, comment in future
        np.testing.assert_allclose(xS2, xS, rtol=1e-3)
        np.testing.assert_allclose(yS2, yS, rtol=1e-3)

        # --- Find shear coefficients and main direction
        GKt = K66 - Kxx*yS**2 -2*Kxy*xS*yS - Kyy*xS**2
        if np.all(np.abs(Kxy)<1e-16):
            # Assumes theta_s=0
            kxsGA = Kxx # Kxx = kxs*GA
            kysGA = Kyy
            theta_s=0*EA
        else:
            if len(Klist)==1:
                kxsGA,kysGA,theta_s = solvexytheta(Kxx,Kyy,Kxy)
                if theta_s>np.pi:
                    theta_s-=2*np.pi
            else:
                kxsGA = np.zeros(Kxx.shape)
                kysGA = np.zeros(Kxx.shape)
                theta_s = np.zeros(Hxx.shape)
                for i, (kxx, kyy, kxy) in enumerate(zip(Kxx,Kyy,Kxy)):
                    kxsGA[i],kysGA[i],theta_s[i] = solvexytheta(kxx,kyy,kxy)

                theta_s[theta_s>np.pi]=theta_s[theta_s>np.pi]-2*np.pi

        # Sanity checks, comment in future
        KK2= KK(EA, EIxp, EIyp, GKt, EA*0+1, kxsGA, kysGA, xC, yC, theta_p, xS, yS, theta_s)
        np.testing.assert_allclose(KK2[0,0], K[0,0], rtol=1e-2)
        np.testing.assert_allclose(KK2[1,1], K[1,1], rtol=1e-2)
        np.testing.assert_allclose(KK2[2,2], K[2,2], rtol=1e-2)
        np.testing.assert_allclose(KK2[3,3], K[3,3], rtol=1e-1)
#         np.testing.assert_allclose(KK2[4,4], K[4,4], rtol=1e-2)
        np.testing.assert_allclose(KK2[5,5], K[5,5], rtol=1e-1)
        np.testing.assert_allclose(KK2[2,3], K[2,3], rtol=1e-2)
        np.testing.assert_allclose(KK2[2,4], K[2,4], rtol=1e-2)

        np.testing.assert_allclose(K16, -Kxx*yS-Kxy*xS)
#         np.testing.assert_allclose(KK2[0,5], K[0,5],rtol=1e-3)
#         np.testing.assert_allclose(KK2[1,5], K[1,5],rtol=5e-2) # Kxy harder to get
        #np.testing.assert_allclose(KK2[3,4], K[3,4]) # <<< hard to match

    else:
        raise NotImplementedError()

    return EA, EIxp, EIyp, kxsGA, kysGA, GKt, xC, yC, xS, yS, theta_p, theta_s

def flat2matlist(M):
    """
    Convert a matrix of size 6x6 where each element is an array of size n
    into a list of n 6x6 matrices
    """
    Mlist=[]
    if not hasattr(M[0,0],'__len__'):
        return [M]
    for iSpan in np.arange(len(M[0,0])):
        Mloc = np.zeros((6,6))
        for i in np.arange(6):
            for j in np.arange(6):
                Mloc[i,j] = M[i,j][iSpan]
        Mlist.append(Mloc)
    return Mlist
