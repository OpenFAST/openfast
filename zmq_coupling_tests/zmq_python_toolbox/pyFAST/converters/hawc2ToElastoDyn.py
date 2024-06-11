import numpy as np
from numpy import cos, sin
import pandas as pd
import os
from pyFAST.input_output.hawc2_htc_file import HAWC2HTCFile
from pyFAST.input_output.csv_file import CSVFile
from pyFAST.input_output.fast_input_file import FASTInputFile, EDBladeFile,EDTowerFile
from pyFAST.tools.pandalib import pd_interp1


def htcToElastoDyn(HTCFile, outDir='./', prefix='', suffix='', bladeBodyName='blade1', towerBodyName='tower', rBlade=None, hTower=None):
    """
    Writes ElastoDyn inputs files from a HAWC2 htc file and the blade body name


    NOTE: This is only a preliminary implementation..
          An attempt is maded for the blade and tower file
          A lot of conversion remains (and might be extremely hard to achieve generically):
            - Shape functions for tower and blade (could use cbeam from welib)
            - Damping for tower and blade
            - Shaft
            - RNA masses and inertia
            - Platform

            - Wont work if bodies are defined using FPM

    INPUTS:
      - HTCFile: path to a htc file
      - bladeBodyName: name of blade in htc file
      - towerBodyName: name of blade in htc file

      - rBlade: array of radial positions for the blade (between 0 and 1)

    """
    htc = HAWC2HTCFile(HTCFile)
    dfs = htc.toDataFrame()


    # --- Blade
    H2BMeanLine  = dfs[bladeBodyName+'_c2']
    H2BStructure = dfs[bladeBodyName+'_st']
    H2BMeanLine = H2BMeanLine[['x_[m]','y_[m]','z_[m]','twist_[deg]']] # Changing order
    if rBlade is None:
        rB_new = np.abs(H2BMeanLine['z_[m]'].values) # NOTE: this should be a curvilinear length, but ElastoDyn is straight, so...
    else:
        rB_new = rBlade* np.abs(H2BMeanLine['z_[m]'].values[-1])
    H2BStructure =  pd_interp1(rB_new, 'r_[m]', H2BStructure)

    # --- Tower
    H2TMeanLine  = dfs[towerBodyName+'_c2']
    H2TStructure = dfs[towerBodyName+'_st']
    H2TMeanLine = H2TMeanLine[['x_[m]','y_[m]','z_[m]','twist_[deg]']] # Changing order
    if hTower is None:
        hT_new = np.abs(H2BMeanLine['z_[m]'].values)
    else:
        hT_new = hBlade* np.abs(H2TMeanLine['z_[m]'].values[-1])
    H2TStructure =  pd_interp1(hT_new, 'r_[m]', H2TStructure)

    # --- ElastoDyn Blade
    bld = EDBladeFile()
    M = np.zeros( (len(rB_new) ,6) )
    M[:,0] = rB_new / rB_new[-1] # BlFract
    #M[:,1] = PitchAxis          # PitchAxis
    M[:,2] = - H2BMeanLine['twist_[deg]'] # StrcTwst [deg]
    M[:,3] =   H2BStructure['m_[kg/m]']   # BMassDen [kg/m]
    M[:,4] =   H2BStructure['E_[N/m^2]'] * H2BStructure['I_x_[m^4]'] # FlpStff [Nm^2] TODO, with twist, make sure that's correct..
    M[:,5] =   H2BStructure['E_[N/m^2]'] * H2BStructure['I_y_[m^4]'] # EdgStff [Nm^2]
    bld['BldProp'] = M
    # TODO find mode shapes
    bldFileName = os.path.join(outDir, '{}ED_blade{}.dat'.format(prefix,suffix))
    bld.write(bldFileName)
    print('Writing ElastoDyn blade file:', bldFileName)

    # --- ElastoDyn Tower
    twr = EDTowerFile()
    M = np.zeros( (len(hT_new), 4) )
    M[:,0] = hT_new / hT_new[-1] # TwFract
    M[:,1] = H2TStructure['m_[kg/m]'] # TMassDen [kg/m]
    M[:,2] =   H2BStructure['E_[N/m^2]'] * H2BStructure['I_x_[m^4]'] # TwFAStif [Nm^2] TODO TODO TODO VERIFY xy
    M[:,3] =   H2BStructure['E_[N/m^2]'] * H2BStructure['I_y_[m^4]'] # TwSSStif [Nm^2]
    twr['TowProp'] = M
    # TODO find mode shapes
    twrFileName = os.path.join(outDir, '{}ED_tower{}.dat'.format(prefix,suffix))
    twr.write(twrFileName)
    print('Writing ElastoDyn tower file:', twrFileName)
