import numpy as np
from numpy import cos, sin
import pandas as pd
import os
# from weio.hawc2_htc_file import HAWC2HTCFile
# from weio.csv_file import CSVFile
# from weio.fast_input_file import FASTInputFile
from pyFAST.input_output.hawc2_htc_file import HAWC2HTCFile
from pyFAST.input_output.csv_file import CSVFile
from pyFAST.input_output.fast_input_file import FASTInputFile
from pyFAST.converters.beam import ComputeStiffnessProps, TransformCrossSectionMatrix
from pyFAST.input_output.fast_input_deck import FASTInputDeck

from .beam import *
from .hawc2 import dfstructure2stfile


def arc_length(points):
    """
    Compute the distances between points along a curve and return those
    cumulative distances as a flat array.

    This function works for 2D, 3D, and N-D arrays.

    Parameters
    ----------
    points : numpy array[n_points, n_dimensions]
        Array of coordinate points that we compute the arc distances for.

    Returns
    -------
    arc_distances : numpy array[n_points]
        Array, starting at 0, with the cumulative distance from the first
        point in the points array along the arc.

    See Also
    --------
    wisdem.commonse.utilities.arc_length_deriv : computes derivatives for
    the arc_length function

    Examples
    --------
    Here is a simple example of how to use this function to find the cumulative
    distances between points on a 2D curve.

    >>> x_values = numpy.linspace(0., 5., 10)
    >>> y_values = numpy.linspace(2., 4., 10)
    >>> points = numpy.vstack((x_values, y_values)).T
    >>> arc_length(points)
    array([0.        , 0.59835165, 1.19670329, 1.79505494, 2.39340658,
           2.99175823, 3.59010987, 4.18846152, 4.78681316, 5.38516481])
    """
    cartesian_distances = np.sqrt(np.sum(np.diff(points, axis=0) ** 2, axis=1))
    arc_distances = np.r_[0.0, np.cumsum(cartesian_distances)]

    return arc_distances


# --------------------------------------------------------------------------------}
# ---beamDynToHawc2
# --------------------------------------------------------------------------------{
def beamDynToHawc2(BD_mainfile, BD_bladefile, H2_htcfile=None, H2_stfile=None, bodyname=None, A=None, E=None, G=None, fstIn=None, theta_p_in=None, FPM=False, verbose=False):
    """ 
    
     FPM: fully populated matrix, if True, use the FPM format of hawc2
    """
    # --- Read BeamDyn files
    if isinstance(BD_mainfile, str):
        BD_mainfile = FASTInputFile(BD_mainfile)
    if isinstance(BD_bladefile, str):
        BD_bladefile = FASTInputFile(BD_bladefile)
    bdLine = BD_mainfile.toDataFrame()
    prop   = BD_bladefile.toDataFrame()

    # Define BeamDyn reference axis
    r_bar = prop['Span'].values
    
    kp_x_raw  = bdLine['kp_xr_[m]'].values
    kp_y_raw  = bdLine['kp_yr_[m]'].values
    kp_z_raw  = bdLine['kp_zr_[m]'].values
    twist_raw = bdLine['initial_twist_[deg]'].values*np.pi/180 # BeamDyn convention

    myarc = arc_length(np.array([kp_x_raw, kp_y_raw, kp_z_raw]).T)
    s_coord = myarc/myarc[-1]

    kp_x = np.interp(r_bar, s_coord, kp_x_raw)
    kp_y = np.interp(r_bar, s_coord, kp_y_raw)
    kp_z = np.interp(r_bar, s_coord, kp_z_raw)
    twist = np.interp(r_bar, s_coord, twist_raw)

    # Create K and M matrices
    K = np.zeros((6,6),dtype='object')
    M = np.zeros((6,6),dtype='object')
    for i in np.arange(6):
        for j in np.arange(6):
            K[i,j]=prop['K{}{}'.format(i+1,j+1)].values
            M[i,j]=prop['M{}{}'.format(i+1,j+1)].values

    if fstIn != None:
        fst = FASTInputDeck(fstIn, verbose=True)
        Bld = fst.fst_vt['AeroDynBlade']
        AF = fst.fst_vt['af_data']
        BlSpn = Bld['BldAeroNodes'][:,0]
        n_span = len(BlSpn)
        le2ac_raw = np.zeros(n_span)
        for iSpan in range(n_span):
            le2ac_raw[iSpan] = fst.fst_vt['ac_data'][iSpan]['AirfoilRefPoint'][0]
        # Define axis of aero centers
        BlCrvAC = Bld['BldAeroNodes'][:,1]
        BlSwpAC = Bld['BldAeroNodes'][:,2]
        chord = Bld['BldAeroNodes'][:,5]
        s_aero = BlSpn/BlSpn[-1]
        ac_x = np.interp(r_bar, s_aero, BlCrvAC)
        ac_y = np.interp(r_bar, s_aero, BlSwpAC)
        le2ac = np.interp(r_bar, s_aero, le2ac_raw) # Leading edge to aerodynamic center (in chord)

        # Get x and y coordinates of c2 axis
        ac2c2 = (0.5 - le2ac) * chord
        c2_x = ac_x + ac2c2 * np.sin(twist)
        c2_y = ac_y + ac2c2 * np.cos(twist)
        # Get offsets from BD axis to c2 axis along the twisted frame of reference
        c2BD_y = np.sqrt( (c2_y - kp_y)**2 + (c2_x - kp_x)**2 )
        c2BD_x = np.zeros_like(c2BD_y) # no x translation, we should be translating along the twisted chord

        # Translate matrices from BD to c2 axis (translate along chord, x and twist are 0)
        transform = TransformCrossSectionMatrix()
        for iSpan in np.arange(len(K[0,0])):
            K_bd_temp = np.zeros((6,6))
            M_bd_temp = np.zeros((6,6))
            for i in np.arange(6):
                for j in np.arange(6):
                    K_bd_temp[i,j] = K[i,j][iSpan]
                    M_bd_temp[i,j] = M[i,j][iSpan]
            K_c2_temp = transform.CrossSectionRotoTranslationMatrix(K_bd_temp, 0., c2BD_y[iSpan], 0.)
            M_c2_temp = transform.CrossSectionRotoTranslationMatrix(M_bd_temp, 0., c2BD_y[iSpan], 0.)
            for i in np.arange(6):
                for j in np.arange(6):
                    K[i,j][iSpan]=K_c2_temp[i,j]
                    M[i,j][iSpan]=M_c2_temp[i,j]

        # Update BeamDyn axis to c2 axis
        kp_x = c2_x
        kp_y = c2_y

    # Map 6x6 data to "beam" data
    # NOTE: theta_* are in [rad]
    EA, EIx, EIy, kxsGA, kysGA, GKt, x_C, y_C, x_S, y_S, theta_p, theta_s = K66toPropsDecoupled(K, theta_p_in)
    m, Ixi, Iyi, Ip, x_G, y_G, theta_i = M66toPropsDecoupled(M)
#     print('kxGA    {:e}'.format(np.mean(kxsGA)))
#     print('kyGA    {:e}'.format(np.mean(kysGA)))
#     print('EA      {:e}'.format(np.mean(EA)))
#     print('EIx     {:e}'.format(np.mean(EIx)))
#     print('EIy     {:e}'.format(np.mean(EIy)))
#     print('GKt     {:e}'.format(np.mean(GKt)))
#     print('xC    ',np.mean(x_C))
#     print('yC    ',np.mean(y_C))
#     print('xS    ',np.mean(x_S))
#     print('yS    ',np.mean(y_S))
#     print('thetap',np.mean(theta_p))
#     print('thetas',np.mean(theta_s))
#     print('m     ',np.mean(m))
#     print('Ixi   ',np.mean(Ixi))
#     print('Iyi   ',np.mean(Iyi))
#     print('Ip    ',np.mean(Ip))
#     print('x_G   ',np.mean(x_G))
#     print('y_G   ',np.mean(y_G))
#     print('thetai',np.mean(theta_i))

    # Convert to Hawc2 system
    if FPM:
        dfMeanLine , dfStructure = beamDyn2Hawc2FPM_raw(r_bar,
                kp_x, kp_y, kp_z, twist,  # BeamDyn convention, twist around -z [in rad]
                m, Ixi, Iyi, x_G, y_G, theta_i,  # theta_i/p around z (in rad)
                x_C, y_C, theta_p, K, M)

    else:

        dfMeanLine , dfStructure = beamDyn2Hawc2_raw(r_bar,
                kp_x, kp_y, kp_z, twist, 
                m, Ixi, Iyi, x_G, y_G, theta_i,
                EA, EIx, EIy, GKt, kxsGA, kysGA, x_C, y_C, theta_p, x_S, y_S, theta_s, 
                A=A, E=E, G=G)

    # --- Rewrite st file
    if H2_stfile is not None:
        try:
            os.makedirs(os.path.dirname(H2_stfile))
        except:
            pass
        if verbose: 
            print('Writing:   ',H2_stfile)

        dfstructure2stfile(dfStructure, H2_stfile)
        #with open(H2_stfile, 'w') as f:
        #    f.write('%i ; number of sets, Nset\n' % 1)
        #    f.write('-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n')
        #    f.write('#%i ; set number\n' % 1)
        #    if FPM:
        #        cols=['r_[m]','m_[kg/m]','x_cg_[m]','y_cg_[m]','ri_x_[m]','ri_y_[m]','pitch_[deg]','x_e_[m]','y_e_[m]','K11','K12','K13','K14','K15','K16','K22','K23','K24','K25','K26','K33','K34','K35','K36','K44','K45','K46','K55','K56','K66']
        #    else:
        #        cols=['r_[m]','m_[kg/m]','x_cg_[m]','y_cg_[m]','ri_x_[m]','ri_y_[m]', 'x_sh_[m]','y_sh_[m]','E_[N/m^2]','G_[N/m^2]','I_x_[m^4]','I_y_[m^4]','I_p_[m^4]','k_x_[-]','k_y_[-]','A_[m^2]','pitch_[deg]','x_e_[m]','y_e_[m]']
        #    f.write('\t'.join(['{:20s}'.format(s) for s in cols])+'\n')
        #    f.write('$%i %i\n' % (1, dfStructure.shape[0]))
        #    f.write('\n'.join('\t'.join('%19.13e' %x for x in y) for y in dfStructure.values))

    # --- Rewrite htc file
    if H2_htcfile is not None:
        def readToMarker(lines, marker, i, nMax=None, noException=False):
            l_sel=[]
            if nMax is None: nMax=len(lines)
            while i<nMax:
                line=lines[i]
                if line.replace(' ','').lower().find(marker)>=0:
                    break
                l_sel.append(line.strip())
                i+=1
            if line.strip().replace(' ','').lower().find(marker)<0:
                if noException:
                    return None, None, None
                else:
                    raise Exception('Marker not found '+ marker)
            return l_sel, line, i

        with open(H2_htcfile, 'r') as f:
            lines_in = f.readlines()
        lines_out = []
        bodyNotFound=True
        iBodyEnd=0
        nBodies=0
        while bodyNotFound and nBodies<10:
            _, line, iBodyStart = readToMarker(lines_in, 'beginmain_body',iBodyEnd)
            _, line, iBodyEnd = readToMarker(lines_in, 'endmain_body', iBodyStart)
            _, line, iBody = readToMarker(lines_in, 'name'+bodyname, iBodyStart, iBodyEnd, True)
            nBodies+=1
            if line is None:
                iBody=-1
            else:
                #print('Body {} found between lines {} and {} '.format(bodyname, iBodyStart+1, iBodyEnd+1))
                bodyNotFound=False
        if nBodies>=10:
            raise Exception('Body {} not found in file'.format(bodyname))

        _, line, iC2Start = readToMarker(lines_in, 'beginc2_def', iBodyStart, iBodyEnd)
        _, line, iC2End   = readToMarker(lines_in, 'endc2_def'  , iC2Start, iBodyEnd)

        _, line, iTIStart = readToMarker(lines_in, 'begintimoschenko_input', iBodyStart, iBodyEnd)
        _, line, iTIEnd   = readToMarker(lines_in, 'endtimoschenko_input'  , iTIStart, iBodyEnd)


        simdir        = os.path.dirname(H2_htcfile)
        H2_stfile_rel = os.path.relpath(H2_stfile, simdir)

        lines_out  = lines_in[:iTIStart+1]
        lines_out += ['      filename {};\n'.format(H2_stfile_rel)]
        if FPM:
            lines_out += ['      FPM 1 ;\n']
        #    lines_out += ['      FPM 0 ;\n']
        lines_out += ['      set 1 1 ;\n']
        lines_out += lines_in[iTIEnd:iC2Start+1]
        lines_out += ['      nsec {} ;\n'.format(dfMeanLine.shape[0])]
        for i, row in dfMeanLine.iterrows():
            lines_out += ['      sec {:4d}\t{:13.6e}\t{:13.6e}\t{:13.6e}\t{:13.6e};\n'.format(i+1, row['x_[m]'],row['y_[m]'],row['z_[m]'],row['twist_[deg]'])]
        lines_out += lines_in[iC2End:]


        if verbose: 
            print('ReWriting: ',H2_htcfile)
        with open(H2_htcfile, 'w') as f:
            f.write(''.join(lines_out))

    return dfMeanLine, dfStructure


def beamDyn2Hawc2FPM_raw(r_bar, kp_x, kp_y, kp_z, twist,
        m, Ixi, Iyi, x_G, y_G, theta_i,  # TODO remove in the future
        x_C, y_C, theta_p,               # TODO remove in the future
        K, M):
    """
    Convert spanwise quantities from BeamDyn to Hawc2
     - kp_x, kp_y, kp_z, twist: keypoints positions/orientation as defined in BeamDyn "main" file

    NOTE: all angles are in radians
    
    """
    import scipy.linalg
    nSpan = len(K[0,0])
    # --- BeamDyn to Hawc2 Structural data
    # Hawc2 = BeamDyn
    x_g     = -y_G
    y_g     =  x_G
    x_e     = -y_C
    y_e     =  x_C
    pitch   = theta_p*180/np.pi # [deg] NOTE: could use theta_p, theta_i or theta_s
    if np.all(np.abs(m)<1e-16):
        ri_y    = m*0
        ri_x    = m*0
    else:
        ri_y    = np.sqrt(Ixi/m)    # [m]
        ri_x    = np.sqrt(Iyi/m)    # [m]

    # Curvilinear position of keypoints (only used to get max radius...)
    dr= np.sqrt((kp_x[1:]-kp_x[0:-1])**2 +(kp_y[1:]-kp_y[0:-1])**2 +(kp_z[1:]-kp_z[0:-1])**2)
    r_p= np.concatenate(([0],np.cumsum(dr)))
    r=r_bar * r_p[-1]

    # From Hawc2/ANBA4 to BeamDyn
    RotMat_H2_BD=np.array([  
            [ 0, 1, 0],
            [-1, 0, 0],
            [ 0, 0, 1]])
    RR_H2_BD = scipy.linalg.block_diag(RotMat_H2_BD,RotMat_H2_BD)
    # From BeamDyn to Hawc2/ANBA4
    RR_BD_H2 = RR_H2_BD.T

    stiff = ComputeStiffnessProps()
    inertia = ComputeInertiaProps()
    transform = TransformCrossSectionMatrix()

    # --- Spanwise quantities for HAWC2 input file
    x_e_h2  = np.zeros(nSpan) # Point where radial force (z) does not contribute to beanding around x or y
    y_e_h2  = np.zeros(nSpan)
    x_g_h2  = np.zeros(nSpan)
    y_g_h2  = np.zeros(nSpan)
    m_h2    = np.zeros(nSpan)
    ri_x_h2 = np.zeros(nSpan)
    ri_y_h2 = np.zeros(nSpan)
    pitch_h2 = np.zeros(nSpan) # angle between xc2 axis and principal axis xe
    KH2=np.zeros((6,6,nSpan))
    for iSpan in np.arange(nSpan):
        K_bd = np.zeros((6,6))
        M_bd = np.zeros((6,6))
        for i in np.arange(6):
            for j in np.arange(6):
                K_bd[i,j] = K[i,j][iSpan]
                M_bd[i,j] = M[i,j][iSpan]
        
        # Rotate BD stiffness matrix to Hawc2 reference system
        K_h2 = (RR_BD_H2).dot(K_bd).dot(RR_BD_H2.T)
        M_h2 = (RR_BD_H2).dot(M_bd).dot(RR_BD_H2.T)
        # Compute coordinates of main points 
        xt , yt = stiff.ComputeTensionCenter(K_h2) # tension/elastic center
        xs , ys = stiff.ComputeShearCenter(K_h2)   # shear center
        xm , ym = inertia.ComputeMassCenter(M_h2)  # inertia
        x_g_h2[iSpan] = xm
        y_g_h2[iSpan] = ym
        x_e_h2[iSpan] = xt
        y_e_h2[iSpan] = yt
        # Inertia properties
        m_h2[iSpan]   = M_h2[0,0]
        #ri_x_h2[iSpan] = TODO
        #ri_x_h2[iSpan] = TODO

        # Compute stiffness matrix with decoupled forces and moments        
        Kdec = stiff.DecoupleStiffness(K_h2)
        # Compute Delta, the rotation angle of principal axes (rad)   
        #Delta =  stiff.PrincipalAxesRotationAngle(Kdec)
        Delta = - stiff.OrientationPrincipalAxesBecas(Kdec) # appears more robust
        pitch_h2[iSpan]=Delta
        # Translate K matrix into EC and rotate it by -Delta
        Kh2 = transform.CrossSectionRotoTranslationMatrix(K_h2, xt, yt, -Delta)
        
        for i in np.arange(6):
            for j in np.arange(6):
                KH2[i,j][iSpan]=Kh2[i,j]

    pitch_h2 *=180/np.pi # [deg] 

    # sanity checks between previous method and new (general) method
    np.testing.assert_allclose(m  , m_h2, rtol=1e-3)
    np.testing.assert_allclose(x_e, x_e_h2, rtol=1e-3)
    np.testing.assert_allclose(y_e, y_e_h2, rtol=1e-3)
    np.testing.assert_allclose(x_g, x_g_h2, rtol=1e-3)
    np.testing.assert_allclose(y_g, y_g_h2, rtol=1e-3)
    np.testing.assert_allclose(theta_p*180/np.pi, -pitch_h2, rtol=1e-3)

    # Using new values (remove in the future)
    pitch = -pitch_h2
    x_e   = x_e_h2
    y_e   = y_e_h2
    x_g   = x_g_h2
    y_g   = y_g_h2
    m     = m_h2
    if np.mean(pitch)>0:
        print('Pitch (delta) is mostly positive')
    else:
        print('Pitch (delta) is mostly negative')


    K11 = KH2[0,0]
    K22 = KH2[1,1]
    K33 = KH2[2,2]
    K44 = KH2[3,3]
    K55 = KH2[4,4]
    K66 = KH2[5,5]

    K12 = KH2[0,1]
    K13 = KH2[0,2]
    K14 = KH2[0,3]
    K15 = KH2[0,4]
    K16 = KH2[0,5]
    K23 = KH2[1,2]
    K24 = KH2[1,3]
    K25 = KH2[1,4]
    K26 = KH2[1,5]
    K34 = KH2[2,3]
    K35 = KH2[2,4]
    K36 = KH2[2,5]
    K44 = KH2[3,3]
    K45 = KH2[3,4]
    K46 = KH2[3,5]
    K55 = KH2[4,4]
    K56 = KH2[4,5]

    # --- Create a data frame with spanwise beam data
    columns=['r_[m]','m_[kg/m]','x_cg_[m]','y_cg_[m]','ri_x_[m]','ri_y_[m]','pitch_[deg]','x_e_[m]','y_e_[m]','K11','K12','K13','K14','K15','K16','K22','K23','K24','K25','K26','K33','K34','K35','K36','K44','K45','K46','K55','K56','K66']
    data = np.column_stack((r, m, x_g, y_g, ri_x, ri_y, pitch, x_e, y_e, K11,K12,K13,K14,K15,K16,K22,K23,K24,K25,K26,K33,K34,K35,K36,K44,K45,K46,K55,K56,K66))
    dfStructure = pd.DataFrame(data=data, columns=columns)

    # --- BeamDyn to Hawc2 Reference axis
    X_H2     = -kp_y
    Y_H2     =  kp_x
    Z_H2     =  kp_z
    twist_H2 = - twist*180/np.pi # -negative of BeamDyn twist [deg]
    columns=['x_[m]', 'y_[m]', 'z_[m]', 'twist_[deg]']
    data = np.column_stack((X_H2, Y_H2, Z_H2, twist_H2))
    dfMeanLine = pd.DataFrame(data=data, columns=columns)

    return dfMeanLine, dfStructure



def beamDyn2Hawc2_raw(r_bar, kp_x, kp_y, kp_z, twist,
        m, Ixi, Iyi, x_G, y_G, theta_i, 
        EA, EIx, EIy, GKt, kxsGA, kysGA, x_C, y_C, theta_p, x_S, y_S, theta_s, 
        A=None, E=None, G=None):
    """ 
    NOTE: all angles are in radians
    """
    # --- BeamDyn to Hawc2 Structural data
    if A is None: A = np.ones(x_G.shape)
    if E is None: E = EA/A
    if G is None: G = E/2/(1+0.3) # Young modulus
    # Hawc2 = BeamDyn
    x_cg    = -y_G
    y_cg    = x_G
    x_sh    = -y_S
    y_sh    = x_S
    x_e     = -y_C
    y_e     = x_C
    I_y     = EIx/E            # [m^4] Hawc2 Iy is wrt to principal bending ye axis
    I_x     = EIy/E            # [m^4] Hawc2 Ix is wrt to principal bending xe axis
    I_p     = GKt/G            # [m^4]
    k_y     = kxsGA/(G*A)
    k_x     = kysGA/(G*A)
    pitch   = theta_p*180/np.pi # [deg] NOTE: could use theta_p, theta_i or theta_s
    if np.all(np.abs(m)<1e-16):
        ri_y    = m*0
        ri_x    = m*0
    else:
        ri_y    = np.sqrt(Ixi/m)    # [m]
        ri_x    = np.sqrt(Iyi/m)    # [m]
    # Curvilinear position of keypoints (only used to get max radius...)
    dr= np.sqrt((kp_x[1:]-kp_x[0:-1])**2 +(kp_y[1:]-kp_y[0:-1])**2 +(kp_z[1:]-kp_z[0:-1])**2)
    r_p= np.concatenate(([0],np.cumsum(dr)))
    r=r_bar * r_p[-1]

    columns=['r_[m]','m_[kg/m]','x_cg_[m]','y_cg_[m]','ri_x_[m]','ri_y_[m]','x_sh_[m]','y_sh_[m]','E_[N/m^2]','G_[N/m^2]','I_x_[m^4]','I_y_[m^4]','I_p_[m^4]','k_x_[-]','k_y_[-]','A_[m^2]','pitch_[deg]','x_e_[m]','y_e_[m]']
    data = np.column_stack((r,m,x_cg,y_cg,ri_x, ri_y, x_sh, y_sh, E, G, I_x, I_y, I_p, k_x, k_y, A, pitch, x_e, y_e))
    dfStructure = pd.DataFrame(data=data, columns=columns)

    # --- BeamDyn to Hawc2 Reference axis
    X_H2     = -kp_y
    Y_H2     = kp_x
    Z_H2     = kp_z
    twist_H2 = - twist*180/np.pi # -[deg]
    columns=['x_[m]', 'y_[m]', 'z_[m]', 'twist_[deg]']
    data = np.column_stack((X_H2, Y_H2, Z_H2, twist_H2))
    dfMeanLine = pd.DataFrame(data=data, columns=columns)

    return dfMeanLine, dfStructure



if __name__=='__main__':
    np.set_printoptions(linewidth=300)

    # --- BeamDyn 2 Hawc 2
    BD_mainfile  = 'solid_beam_BeamDyn.dat'
    BD_bladefile = '../solid_beam_BeamDyn_Blade.dat'
    H2_htcfile_old  = './_template.htc'
    H2_htcfile_new  = './solid_beam_hawc2.htc'
    H2_stfile    = './solid_beam_st.dat'

    from shutil import copyfile
    copyfile(H2_htcfile_old, H2_htcfile_new)

    beamDyn2Hawc2(BD_mainfile, BD_bladefile, H2_htcfile_new, H2_stfile, 'beam_1', A=None, E=None, G=None)



