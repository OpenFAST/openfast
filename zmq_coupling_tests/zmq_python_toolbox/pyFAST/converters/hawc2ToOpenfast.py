""" 

Set of tools to convert a hawc2 model to OpenFAST

- hawc2toAD: write aerodyn files from HAWC2model

"""
import numpy as np
import pandas as pd
import os
# Local
from pyFAST.input_output.hawc2_htc_file   import HAWC2HTCFile
from pyFAST.input_output.hawc2_ae_file    import HAWC2AEFile
from pyFAST.input_output.hawc2_pc_file    import HAWC2PCFile
from pyFAST.input_output.fast_input_file  import FASTInputFile, ADBladeFile, ADPolarFile
from pyFAST.airfoils.Polar import Polar



# --------------------------------------------------------------------------------}
# --- Helper functions 
# --------------------------------------------------------------------------------{
def multiInterp(x, xp, fp, extrap='bounded'):
    """ 
    Interpolate all the columns of a matrix `fp` based on new values `x`
    INPUTS:
      - x  : array ( n ), new values
      - xp : array ( np ), old values
      - fp : array ( nval, np), matrix values to be interpolated
    """
    # Sanity
    x   = np.asarray(x)
    xp  = np.asarray(xp)
    assert fp.shape[1]==len(xp), 'Second dimension of fp should have the same length as xp'

    j   = np.searchsorted(xp, x) - 1
    dd  = np.zeros(len(x))
    bOK = np.logical_and(j>=0, j< len(xp)-1)
    bLower =j<0
    bUpper =j>=len(xp)-1
    jOK = j[bOK]
    dd[bOK] = (x[bOK] - xp[jOK]) / (xp[jOK + 1] - xp[jOK])
    jBef=j 
    jAft=j+1
    # 
    # Use first and last values for anything beyond xp
    jAft[bUpper] = len(xp)-1
    jBef[bUpper] = len(xp)-1
    jAft[bLower] = 0
    jBef[bLower] = 0
    if extrap=='bounded':
        pass
        # OK
    elif extrap=='nan':
        dd[~bOK] = np.nan
    else:
        raise NotImplementedError()

    return (1 - dd) * fp[:,jBef] + fp[:,jAft] * dd

            
def interpH2Polar(ae, pc, radius, alpha=None, ae_set_nr=1):
    """ 
    Interpolate hawc2 polar to a given radial station
    ae: an AEFile instance from wetb
    pc: an PCFile instance from wetb
    radius: (scalar) radial position where the polar is to be computed
            
    """
    r_ae = ae.ae_sets[ae_set_nr][:,0]
    r_min = np.min(r_ae)
    r_max = np.max(r_ae)
    if radius<r_min or radius>r_max:
        raise Exception('Radius ({}) needs to be between {} and {}'.format(radius, r_min, r_max))

    thickness = ae.thickness(radius, ae_set_nr)
    pc_set_nr = ae.pc_set_nr(radius, ae_set_nr)
    thicknesses, profiles = pc.pc_sets[pc_set_nr]
    index = np.searchsorted(thicknesses, thickness)
    if index == 0:
        index = 1
    Cx0, Cx1 = profiles[index - 1:index + 1]
    alpha0 =  Cx0[:, 0]
    alpha1 =  Cx1[:, 0]
    if alpha is None:
        alpha = alpha1
    Cxi0 = multiInterp(alpha, alpha0, Cx0.T, extrap='bounded').T
    Cxi1 = multiInterp(alpha, alpha1, Cx1.T, extrap='bounded').T
    th0, th1 = thicknesses[index - 1:index + 1]
    M = Cxi0 + (Cxi1 - Cxi0) * (thickness - th0) / (th1 - th0)
    M[:,0] = alpha # to avoid numerics
    return M, thickness, pc_set_nr, [th0, th1]



# --------------------------------------------------------------------------------}
# --- From HAWC2 to AeroDyn  
# --------------------------------------------------------------------------------{
def AE_PC_C2_toAD(ae_filename, pc_filename, blade_c2def, r_AD=None, ADbldFilename_out=None, polarFilebase_out='Polar_', ae_set_nr=1, correction3D=True, tsr=9):
    """ 
    Convert aerodynamic data from a hawc2 model (AE, PC, balde C2def)
    to AeroDyn (AD) files.

    """
    # --- Read AE and PC file
    ae = HAWC2AEFile(ae_filename)
    pc = HAWC2PCFile(pc_filename)
    ae_base = os.path.basename(ae_filename)
    pc_base = os.path.basename(pc_filename)
    #print(ae)
    #print(pc)

    # --- Setting mean line
    # C2def
    if blade_c2def.shape[1]==4:
        iOff=0
    elif blade_c2def.shape[1]==5:
        iOff=1
    else:
        raise Exception('blade_c2_def should have 4 or 5 columns')
    xC2     = blade_c2def[:,0+iOff]
    yC2     = blade_c2def[:,1+iOff]
    zC2     = blade_c2def[:,2+iOff]
    twistC2 = blade_c2def[:,3+iOff] # [deg]

    ae_set = ae.data.ae_sets[ae_set_nr]
    radius_H2 = ae_set[:,0]

    # --- Default radius if user does not provided one
    if r_AD is None:
        r_AD = radius_H2

    # --- Aerodynamic spanwise data
    # Interpolating AE data to user requested positions
    chord_H2 = np.interp(r_AD, ae_set[:,0], ae_set[:,1])
    trel_H2  = np.interp(r_AD, ae_set[:,0], ae_set[:,2])
    pcset_H2 = np.interp(r_AD, ae_set[:,0], ae_set[:,3])
    twist_H2 = np.interp(r_AD, zC2, twistC2)            #[deg]
    x_H2     = np.interp(r_AD, zC2, xC2) # sweep
    y_H2     = np.interp(r_AD, zC2, yC2) # prebend
    # Aerodynamic Center
    dx_H2 = chord_H2/4 * np.cos(twist_H2*np.pi/180)
    dy_H2 = chord_H2/4 * np.sin(twist_H2*np.pi/180)
    xAC_H2=(x_H2+dx_H2)
    yAC_H2=(y_H2+dy_H2)

    # --- AeroDyn nodes
    # BlSpn        BlCrvAC        BlSwpAC        BlCrvAng       BlTwist        BlChord          BlAFID
    aeroNodes = np.zeros((len(r_AD), 7))
    aeroNodes[:,0] = r_AD
    aeroNodes[:,1] =  yAC_H2  # BlCrvAC # NOTE: not c2def but AC
    aeroNodes[:,2] = -xAC_H2  # BlSwpAC (positive along yOF, negative xH2) # NOTE: not c2def but AC
    dr = np.gradient(aeroNodes[:,0])
    dx = np.gradient(aeroNodes[:,1])
    aeroNodes[:,3] = np.degrees(np.arctan2(dx,dr))
    aeroNodes[:,4] = -twist_H2                     # [deg]
    aeroNodes[:,5] =  chord_H2
    aeroNodes[:,6] = (np.arange(len(r_AD))+1).astype(int) # One polar per radius..
    # Write to disk if needed
    if ADbldFilename_out is not None:
        Bld = ADBladeFile()
        Bld.comment='Generated using HAWC2 inputs: AE:{} PC:{}'.format(ae_base, pc_base)
        Bld['BldAeroNodes'] = aeroNodes
        Bld.write(ADbldFilename_out)
    aeroNodes = pd.DataFrame(data=aeroNodes, columns=['BlSpn_[m]', 'BlCrvAC_[m]', 'BlSwpAC_[m]', 'BlCrvAng_[deg]', 'BlTwist_[deg]', 'BlChord_[m]', 'BlAFID_[-]'])


    # --- Write Polars for each radial station, interpolating on thickness
    polarFilenames=[]
    polars=[]
    vAlpha = np.arange(361)-180 # make it an argument
    nAlpha = len(vAlpha)

    writePolars = polarFilebase_out is not None
    if writePolars:
        baseDir = os.path.dirname(polarFilebase_out)
        try:
            os.makedirs(baseDir, exist_ok=True)
        except:
            pass
        # Create a default polar file
        pol = ADPolarFile()

    r_max = np.max(r_AD)
    for ir, (r,c) in enumerate(zip(r_AD, chord_H2)):
        M, t, pc_set, thicknesses = interpH2Polar(ae.data, pc.data, r, alpha=vAlpha, ae_set_nr=ae_set_nr)
        comment = 'Thickness: {} - pc_set:{} - thicknesses:{}\nGenerated using HAWC2 inputs: AE:{} PC:{}'.format(t, pc_set, thicknesses, ae_base, pc_base)
        Re = 1.0 # TODO
        # Ensure that first value match last value
        M[-1,1:] = M[0,1:]
        # Create an instance of Polar class for convenience
        P = Polar(Re=Re, alpha=vAlpha, cl=M[:,1], cd=M[:,2], cm=M[:,3], radians=False)
        # Apply 3D correction
        if r>0 and correction3D:
            try:
                P = P.correction3D(
                    r_over_R     = r/r_max,
                    chord_over_r = c/r,
                    tsr=tsr, 
                    lift_method="DuSelig", drag_method="None", blending_method="linear_25_45",
                    max_cl_corr=0.25)
            except:
                print('3D correction not applied at station ', ir)

        # Store polars
        M = np.column_stack((P.alpha, P.cl, P.cd, P.cm))
        polars.append(M)

        # Write polar to disk
        if writePolars:
            polarFilename_out = polarFilebase_out+'{:03d}.dat'.format(ir+1)
            P.toAeroDyn(filenameOut=polarFilename_out, templateFile=pol, Re=Re, comment=comment)
            #(alpha0,alpha1,alpha2,cnSlope,cn1,cn2,cd0,cm0)=P.unsteadyParams()
            #pol['AFCoeff']= M 
            #pol.comment= comment
            #pol.write(polarFilename_out)
            polarFilenames.append(polarFilename_out)

    return aeroNodes, polars, polarFilenames


def hawc2toAD(htc, r_AD=None, ADbldFilename_out=None, polarFilebase_out='Polar_', bladeBodyName='blade1', correction3D=True, tsr=9):
    """ 
    Convert aerodynamic data from a hawc2 model (taken using an htc file as input)
    to AeroDyn (AD) files.

    INPUTS:
     - htc: filename of main htc file, or instance of HAWC2HTCFile
     - r_AD: radial positions for AeroDyn file. If None, use hawc2 AE stations
     - ADbldFilename_out: filename for AeroDyn blade. If None: nothing is written.
     - polarFilebase_out: filebase, used to write the polar files. If None: nothing is written
     - bladeBodyName: body name of blade in htc file

    TODO:
     - Position of aerodynamic center assuemd to be at c/4 from c2def. TODO
    """
    if isinstance(htc, HAWC2HTCFile):
        pass
    else:
        # Read htc file
        htc = HAWC2HTCFile(htc)
    # Extra c2 def of blade
    bdy = htc.bodyByName(bladeBodyName)
    c2def = htc.bodyC2(bdy)
    # Extract ae and pc filenames
    if 'aero' not in htc.data.keys():
        raise Exception('Aero section not found in htc file: '.format(htcFilename))
    aero = htc.data['aero']
    ae_filename = os.path.join(htc.data.modelpath, aero.ae_filename[0])
    pc_filename = os.path.join(htc.data.modelpath, aero.pc_filename[0])

    # Convert to AeroDyn
    return AE_PC_C2_toAD(ae_filename, pc_filename, c2def, r_AD=r_AD, ADbldFilename_out=ADbldFilename_out, polarFilebase_out=polarFilebase_out, correction3D=correction3D, tsr=9)


