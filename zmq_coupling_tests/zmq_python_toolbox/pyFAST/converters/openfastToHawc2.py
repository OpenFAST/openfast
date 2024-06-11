import os
import numpy as np
import pandas as pd

from pyFAST.input_output.hawc2_htc_file import HAWC2HTCFile
from pyFAST.input_output.hawc2_ae_file import HAWC2AEFile
from pyFAST.input_output.hawc2_pc_file import HAWC2PCFile
from pyFAST.input_output.csv_file import CSVFile
from pyFAST.input_output.fast_input_file import FASTInputFile
from pyFAST.input_output.fast_input_deck import FASTInputDeck

import pyFAST.converters.beamdyn as bd
import pyFAST.converters.elastodyn as ed
import pyFAST.converters.hawc2 as h2


def FAST2Hawc2(fstIn, htcTemplate, htcOut, OPfile=None, TwrFAFreq=0.1, TwrSSFreq=0.1, SftTorFreq=4, FPM = False, Bld_E=None, Bld_G=None, Bld_A=None, Bld_theta_p=None):
    """ 
    Write a Hawc2 and hawcstab2 model from an openfat model
    """

    # --- OpenFAST
    fst   = FASTInputDeck(fstIn, verbose=True)
    IW    = fst.fst_vt['InflowWind']
    ED    = fst.fst_vt['ElastoDyn']
    AD    = fst.fst_vt['AeroDyn15']
    Bld   = fst.fst_vt['AeroDynBlade']
    AF    = fst.fst_vt['af_data']
    twrOF = fst.fst_vt['ElastoDynTower']
    BD    = fst.fst_vt['BeamDyn']
    BDbld = fst.fst_vt['BeamDynBlade']
    EDbld = fst.fst_vt['ElastoDynBlade']

    # print(fst.ED.keys())
    # print(fst.ED['TwrFile'])
    # print(fst.AD.keys())

    # --- Derived parameters
    outDir = os.path.dirname(htcOut)
    if os.path.basename(outDir)=='htc':
        outDir=os.path.dirname(outDir) # we go one more level up
    outDataDir = os.path.join(outDir, 'data')

    try:
        os.mkdir(outDir)
    except:
        pass
    try:
        os.mkdir(outDataDir)
    except:
        pass


    # --- Read template htc file
    htc = HAWC2HTCFile(htcTemplate)
    # --------------------------------------------------------------------------------}
    # --- Tower
    # --------------------------------------------------------------------------------{
    # Damping
    bdy = htc.bodyByName('tower')
    bdy.damping_posdef.values[3] = (2*twrOF['TwrFADmp(1)']/100) / (TwrFAFreq*2*np.pi) # 2 zeta / omega
    bdy.damping_posdef.values[4] = (2*twrOF['TwrSSDmp(1)']/100) / (TwrSSFreq*2*np.pi) # 2 zeta / omega
    bdy.damping_posdef.values[5] = 0.01 # TODO TODO
    # Mean line
    val = htc.bodyC2(bdy)
    val[:,3] = np.linspace(0., -ED['TowerHt'], val.shape[0])
    htc.setBodyC2(bdy, val)

    # --------------------------------------------------------------------------------}
    # --- Towertop/nacelle 
    # --------------------------------------------------------------------------------{
    bdy = htc.bodyByName('towertop')
    # Option 1, at node 2
    # bdy.concentrated_mass.values[0] = 2
    # bdy.concentrated_mass.values[1] = np.around(ED['NacCMyn'],4)
    # bdy.concentrated_mass.values[2] = np.around(ED['NacCMxn'],4)
    # bdy.concentrated_mass.values[3] = np.around(ED['Twr2Shft']-ED['NacCMzn'],4) # z distance from last towertopnode to NacCOG. NOTE: z hawc2 negative
    # bdy.concentrated_mass.values[4] = np.around(ED['NacMass'] ,4)
    # bdy.concentrated_mass.values[7] = np.around(ED['NacYIner']-ED['NacMass']*ED['NacCMxn']**2,2) # Jz at CM
    # bdy.concentrated_mass.values[5] = bdy.concentrated_mass.values[7] # Jx ~ Jz 
    # bdy.concentrated_mass.values[6] = bdy.concentrated_mass.values[7]/10 # Jy ~ Jz /10 ? TODO
    # Option 2, at node 1
    bdy.concentrated_mass.values[0] = 1
    bdy.concentrated_mass.values[1] = np.around(ED['NacCMyn'],4)
    bdy.concentrated_mass.values[2] = np.around(ED['NacCMxn'],4)
    bdy.concentrated_mass.values[3] = np.around(ED['NacCMzn'],4)
    bdy.concentrated_mass.values[4] = np.around(ED['NacMass'] ,4)
    NacInerZ  = ED['NacYIner'] # TODO
    NacInerFA = NacInerZ       # TODO
    NacInerSS = NacInerZ/10    # TODO
    bdy.concentrated_mass.values[5] = NacInerFA
    bdy.concentrated_mass.values[6] = NacInerSS
    bdy.concentrated_mass.values[7] = NacInerZ

    # Mean line
    val = htc.bodyC2(bdy)
    val[:,3] = np.linspace(0, -ED['Twr2Shft'], val.shape[0])
    htc.setBodyC2(bdy, val)
    # print(bdy)

    # --------------------------------------------------------------------------------}
    # --- Shaft 
    # --------------------------------------------------------------------------------{
    bdy = htc.bodyByName('shaft')
    # Torsion
    SftOmega = SftTorFreq * 2 * np.pi
    SftMass = ED['DTTorSpr']/SftOmega                    # Shaft omega = sqrt(k/m)
    if SftMass >0:
        SftZeta = ED['DTTorDmp']*SftOmega/(2*ED['DTTorSpr']) # zeta = c*omega/(2*)
        SftBetaTors = 2*SftZeta / SftOmega                       # beta = 2 zeta/omega
    else:
        SftBetaTors = 1e-3
    SftBetaBending = 4e-4 # TODO
    # Damping
    bdy.damping_posdef.values[3] = SftBetaBending
    bdy.damping_posdef.values[4] = SftBetaBending
    bdy.damping_posdef.values[5] = SftBetaTors
    # Mass and inertia
    bdy.concentrated_mass.values[7] = np.around(ED['GenIner']*ED['GBRatio']**2, 2)
    bdy.concentrated_mass__2.values[7] = np.around(ED['HubIner'],2)
    bdy.concentrated_mass__2.values[4] = np.around(ED['HubMass'],2)
    # Mean line
    val = htc.bodyC2(bdy)
    val[:,3] = np.linspace(0, -ED['OverHang'], val.shape[0])
    htc.setBodyC2(bdy, val)
    # print(bdy)

    # --------------------------------------------------------------------------------}
    # --- Hub 
    # --------------------------------------------------------------------------------{
    bdy = htc.bodyByName('hub1')
    # Mean line
    val = htc.bodyC2(bdy)
    val[:,3] = np.linspace(0, ED['HubRad'], val.shape[0])
    htc.setBodyC2(bdy, val)
    # print(bdy)

    # --------------------------------------------------------------------------------}
    # --- Blade 
    # --------------------------------------------------------------------------------{
    # Mean line
    # val[:,3] = np.linspace(0, ED['HubRad'], val.shape[0])

    st_filefull = os.path.join(outDataDir, 'blade_st.st') # Full path of blade st file
    st_file     = os.path.relpath(st_filefull, outDir)    # Relative to output dir
    dfStructure = None
    dfMeanLine  = None
    damp=[0, 0, 0]
    if BD is not None:
        # --- Convert from BeamDyn to HAWC2 c2def and st data
        # also writes st file
        # TODO A, E, G, theta_p
        dfMeanLine, dfStructure = bd.beamDynToHawc2(BD, BDbld, H2_stfile=st_filefull, A=Bld_A, E=Bld_E, G=Bld_G, theta_p_in = Bld_theta_p, FPM=FPM, verbose=True)

        # --- Damping
        damp[0] = BDbld['DampingCoeffs'][0,1] # TODO  Check order
        damp[1] = BDbld['DampingCoeffs'][0,0] # TODO 
        damp[2] = BDbld['DampingCoeffs'][0,2] # TODO 
            
    elif EDbld is not None:
        print('[WARN] Blade with ElastoDyn is not fully implemented')
        if FPM:
            print('[WARN] Cannot do FPM with ElastoDyn for now')
            FPM = False # Cannot do FPM with ElastoDyn

        # --- ElastoDyn blade properties
        M = EDbld['BldProp']
        dfBld = EDbld.toDataFrame() # r/R, PitchAxis(unused), StructuralTwist, m, EIflap, EIdge
        R           = ED['TipRad'] - ED['HubRad']

        # --- Convert from ElastoDyn to HAWC2 c2def and st data
        dfMeanLine, dfStructure = ed.elastodynBlade2Hawc2_raw(dfBld, R)
        # Write st file
        h2.dfstructure2stfile(dfStructure, st_filefull)

        # --- 
        # TODO TODO TODO Damping is completely wrong here
        print('[WARN] Damping values for blades are wrong when using ElastoDyn blade')
        damp[0] = EDbld['BldFlDmp(1)']/100 # 
        damp[1] = EDbld['BldFlDmp(2)']/100
        damp[2] = EDbld['BldEdDmp(1)']/100


    else:
        raise Exception('No BeamDyn or ElastoDyn blade present')
        # or throw a warning and skip the section below..

    # --- Setup HAWC2 "body" data
    bdy = htc.bodyByName('blade1')
    bdy.nbodies = dfMeanLine.shape[0]-1 # One body per station -1
    bdy.timoschenko_input.ts_filename = st_file
    bdy.timoschenko_input.set.values = [1,1]

    if FPM:
        bdy.timoschenko_input.fpm = 1
        # Damping
        #bdy['damping_aniso'] = bdy.damping_posdef.copy()
        bdy.damping_posdef.name_='damping_aniso'
        bdy.damping_posdef.name_='damping_aniso'
        bdy.damping_posdef.values[0]=0 # no mass proportional damping
        bdy.damping_posdef.values[1]=0 # no mass proportional damping
        bdy.damping_posdef.values[2]=0 # no mass proportional damping
        bdy.damping_posdef.values[3] = damp[0]
        bdy.damping_posdef.values[4] = damp[1]
        bdy.damping_posdef.values[5] = damp[2]
        #raise NotImplementedError('Damping for FPM')
        print('>>>> TODO TODO TODO DAMPING ')
    else:
        bdy.timoschenko_input.fpm = 0
        # Damping
        bdy.damping_posdef.values[3] = damp[0]
        bdy.damping_posdef.values[4] = damp[1]
        bdy.damping_posdef.values[5] = damp[2]

    # Meanline
    c2_def = htc.bodyC2(bdy)
    c2_def = np.column_stack((np.arange(1,len(dfMeanLine)+1), dfMeanLine.values))
    c2_def = np.around(c2_def, 5)
    htc.setBodyC2(bdy, c2_def)

    #print(bdy)


    # --- Orientation
    sec = htc.data.new_htc_structure.orientation 
    # tower top 2 shaft
    #print(sec.relative__2.body2_eulerang__2)
    sec.relative__2.body2_eulerang__2.values[0] = -ED['ShftTilt']
    sec.relative__2.body2_ini_rotvec_d1.values[3] = ED['RotSpeed']/60.*2*np.pi
    # shaft 2 hubs
    sec.relative__3.body2_eulerang__3.values[0] = - ED['PreCone(1)']
    sec.relative__4.body2_eulerang__3.values[0] = - ED['PreCone(2)']
    sec.relative__5.body2_eulerang__3.values[0] = - ED['PreCone(3)']
    # hubs 2 shaft
    # sec.relative__6.body2_eulerang.values[2] = - np.around(ED['BlPitch(1)'],4) # TODO sign
    # sec.relative__7.body2_eulerang.values[2] = - np.around(ED['BlPitch(2)'],4) # TODO sign
    # sec.relative__8.body2_eulerang.values[2] = - np.around(ED['BlPitch(3)'],4) # TODO sign
    #print(sec)

    # --------------------------------------------------------------------------------}
    # --- Wind
    # --------------------------------------------------------------------------------{
    sec = htc.data.wind
    sec.density                = AD['AirDens']
    sec.wsp                    = IW['HWindSpeed']
    sec.center_pos0.values[2]  = -IW['RefHt']
    sec.shear_format.values[1] = IW['PLexp']
    sec.tower_shadow_method    = min(AD['TwrPotent'], 1)
    sec.tower_shadow_potential.radius.values[0:2]      = [0,                   AD['TowProp'][0, 1]/2]
    sec.tower_shadow_potential.radius__2.values[0:2]   = [-AD['TowProp'][-1,0], AD['TowProp'][-1,1]/2]
    sec.tower_shadow_potential_2.radius.values[0:2]    = [0,                   AD['TowProp'][0, 1]/2]
    sec.tower_shadow_potential_2.radius__2.values[0:2] = [ AD['TowProp'][-1,0], AD['TowProp'][-1,1]/2]
    # print(sec)

    # --------------------------------------------------------------------------------}
    # --- AeroDynamics 
    ## --- AE File
    radius    = Bld['BldAeroNodes'][:,0]
    twist     = Bld['BldAeroNodes'][:,4]
    chord     = Bld['BldAeroNodes'][:,5]
    polar_id  = Bld['BldAeroNodes'][:,6].astype(int)
    pc_set_id = np.ones(chord.shape).astype(int) # Only one set
    thickness = np.around((1-(polar_id-np.min(polar_id))/(np.max(polar_id)+1-np.min(polar_id)))*100,2) # fake thickness
    #print(thickness)
    ae = HAWC2AEFile()
    ae.add_set(radius=radius, chord=chord, thickness=thickness, pc_set_id= pc_set_id)
    ae_filenamefull = os.path.join(outDataDir, 'ae_file.ae')
    ae.write(ae_filenamefull)
    ae_filename = os.path.relpath(ae_filenamefull, outDir)

    # --- PC File
    pc = HAWC2PCFile()
    upolar_id, indices = np.unique(polar_id, return_index=True)
    uthickness = np.flip(thickness[indices])                      # NOTE: we flip 
    polars     = [AF[i-1]['AFCoeff'] for i in np.flip(upolar_id)] # NOTE: we flip  

    thicknesses = uthickness
    pc.add_set(1, thicknesses, polars)
    pc_filenamefull = os.path.join(outDataDir, 'pc_file.pc')
    pc.write(pc_filenamefull)
    pc_filename = os.path.relpath(pc_filenamefull, outDir)

    # --- Aero HTC
    sec = htc.data.aero
    sec.pc_filename = pc_filename
    sec.ae_filename = ae_filename
    sec.induction_method 
    sec.induction_method =  AD['WakeMod']
    sec.tiploss_method   =  1 if AD['TipLoss'] else 0
    sec.aerosections = Bld['BldAeroNodes'].shape[0]

    if AD['AFAeroMod']==2:
        sec.dynstall_method= 2
        if AD['UAMod']==6:
            sec.dynstall_method= 1
    else:
        sec.dynstall_method  =  0
    # print(sec)

    # print(htc.data)


    # --------------------------------------------------------------------------------}
    # --- Hawcstab2 
    # --------------------------------------------------------------------------------{
    # Write operating point file
    if OPfile is not None:
        dfOP = CSVFile(OPfile).toDataFrame()
        WS0=dfOP['Wind Speed (m/s)']
        WS0,Pitch0,RPM0,P0,T0 = dfOP['Wind Speed (m/s)'], dfOP['Pitch Angle (deg)'], dfOP['Rotor Speed (rpm)'], dfOP['Aerodynamic Power (W)']/1000, dfOP['Thrust (N)']/1000 

        # interpolating data to less points
        WS=np.arange(np.min(dfOP['Wind Speed (m/s)']),np.max(dfOP['Wind Speed (m/s)'])+0.5,1)
        Pitch = np.interp(WS, WS0, Pitch0)
        RPM   = np.interp(WS, WS0, RPM0)
        P     = np.interp(WS, WS0, P0)
        T     = np.interp(WS, WS0, T0)

        op_filenamefull = os.path.join(outDataDir, 'operational_data.opt')
        op_filename = os.path.relpath(op_filenamefull, outDir)
        with open(op_filenamefull, 'w') as f:
            colnames=['Wind speed [m/s]','Pitch [deg]','Rot. speed [rpm]','Aero power [kW]','Aero thrust [kN]']
            f.write('{}\t{}\n'.format(len(WS), '\t'.join(colnames)))
            for (ws,pi,rpm,p,t) in zip(WS,Pitch,RPM,P,T):
                f.write('{}\t{}\t{}\t{}\t{}\n'.format(ws,pi,rpm,p,t))

        sec=htc.data.hawcstab2
        sec.operational_data_filename = op_filename


    # --- Write htc for hawc2
    htc.write(htcOut)

    # --- Write htc for hawcstab2
    del htc.data.output
    del htc.data.simulation
    del htc.data.dll
    del htc.data.output
    del htc.data.wind.tower_shadow_potential_2  # hawcstab2 can't handle this subblock
    del htc.data.wind.mann
    del htc.data.aerodrag  # can't handle this either fuck annoying
    del htc.data.wind.tower_shadow_potential
    del htc.data.wind.tower_shadow_potential_2
    htc.data.wind.turb_format = 0
    htc.data.wind.tower_shadow_method = 0
    htc.data.wind.shear_format.values[0] = 1

    # interpolate blde to less stations
    bdy = htc.bodyByName('blade1')
    bdy.nbodies = 24
    c2_def = htc.bodyC2(bdy)
    z0 = dfMeanLine['z_[m]'].values
    z2 = np.linspace(z0[0],z0[-1],25)
    dfMeanLine2 = pd.DataFrame(columns=dfMeanLine.columns)
    dfMeanLine2['z_[m]']       = z2
    dfMeanLine2['x_[m]']       = np.interp(z2, z0, dfMeanLine['x_[m]'].values)
    dfMeanLine2['y_[m]']       = np.interp(z2, z0, dfMeanLine['y_[m]'].values)
    dfMeanLine2['twist_[deg]'] = np.interp(z2, z0, dfMeanLine['twist_[deg]'].values)
    c2_def = np.column_stack((np.arange(1,len(dfMeanLine2)+1), dfMeanLine2.values))
    c2_def = np.around(c2_def, 5)
    htc.setBodyC2(bdy, c2_def)

    htc.write(htcOut.replace('.htc','_hs2.htc'))

