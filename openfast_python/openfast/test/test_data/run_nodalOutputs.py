"""

Example script to compute the steady-state performance in OpenFAST with Nodal outputs for AeroDyn, ElastoDyn and BeamDyn

"""
import weis
from weis.aeroelasticse.runFAST_pywrapper import runFAST_pywrapper_batch
from weis.aeroelasticse.CaseGen_General import CaseGen_General
import numpy as np
import os
import platform
import shutil


if __name__ == '__main__':

    # Paths calling the standard modules of WEIS
    fastBatch = runFAST_pywrapper_batch()
    run_dir1                    = os.path.dirname( os.path.dirname( os.path.realpath(weis.__file__) ) )
    fastBatch.FAST_exe          = shutil.which( 'openfast' )   # Path to executable
    run_dir2                    = os.path.dirname( os.path.realpath(__file__) ) + os.sep
    fastBatch.FAST_directory    = os.path.join(run_dir2, 'OpenFAST_models','IEA-15-240-RWT','IEA-15-240-RWT-NodalOutputs')   # Path to fst directory files
    fastBatch.FAST_InputFile    = 'IEA-15-240-RWT-NodalOutputs.fst'   # FAST input file (ext=.fst)
    fastBatch.FAST_runDirectory = 'steady_state-nodal/iea15mw'
    fastBatch.debug_level       = 2

    # User settings
    n_cores     = 2     # Number of available cores
    TMax        = 0.01  # Length of wind grids and OpenFAST simulations, suggested 720 s
    cut_in      = 3.    # Cut in wind speed
    cut_out     = 25.   # Cut out wind speed
    wind_speeds = [10] #np.linspace(cut_in, cut_out, int(n_ws)) # Wind speeds to run OpenFAST at
    Ttrans      = min([0., TMax - 0.])  # Start of the transient for DLC with a transient, e.g. DLC 1.4
    TStart      = min([0., TMax - 0.]) # Start of the recording of the channels of OpenFAST

    # Initial conditions for ElastoDyn
    u_ref       = np.arange(3.,26.) # Wind speed vector to specify the initial conditions
    pitch_ref   = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.5058525323662666, 5.253759185225932, 7.50413344606208, 9.310153958810268, 10.8972969450052, 12.412247669440042, 13.883219268525659, 15.252012626933068, 16.53735488246438, 17.76456777500061, 18.953261878035104, 20.11055307762722, 21.238680277668898, 22.30705111326602, 23.455462501156205] # Pitch values in deg
    omega_ref   = [2.019140272160114, 2.8047214918577925, 3.594541645994511, 4.359025795823625, 5.1123509774611025, 5.855691196288371, 6.589281196735111, 7.312788026081227, 7.514186181824161, 7.54665511646938, 7.573823812448151, 7.600476033113538, 7.630243938880304, 7.638301051122195, 7.622050377183605, 7.612285710588359, 7.60743945212863, 7.605865650155881, 7.605792924227456, 7.6062185247519825, 7.607153933765292, 7.613179734210654, 7.606737845170748] # Rotor speeds in rpm
    pitch_init = np.interp(wind_speeds, u_ref, pitch_ref)
    omega_init = np.interp(wind_speeds, u_ref, omega_ref)

    # Settings passed to OpenFAST
    case_inputs = {}
    case_inputs[("Fst","TMax")]             = {'vals':[TMax], 'group':0}
    case_inputs[("Fst","DT")]               = {'vals':[0.0002], 'group':0}
    case_inputs[("Fst","CompInflow")]       = {'vals':[1], 'group':0}
    case_inputs[("Fst","CompServo")]        = {'vals':[1], 'group':0}
    case_inputs[("Fst","OutFileFmt")]       = {'vals':[1], 'group':0}
    case_inputs[("Fst","DT_Out")]           = {'vals':[0.01], 'group':0}
    case_inputs[("ElastoDyn","GenDOF")]     = {'vals':['True'], 'group':0}
    case_inputs[("ServoDyn","PCMode")]      = {'vals':[5], 'group':0}
    case_inputs[("ServoDyn","VSContrl")]    = {'vals':[5], 'group':0}
    case_inputs[("InflowWind","WindType")]  = {'vals':[1], 'group':0}
    case_inputs[("InflowWind","HWindSpeed")]= {'vals': wind_speeds, 'group': 0}
    case_inputs[("Fst","OutFileFmt")]       = {'vals':[0], 'group':0}
    case_inputs[("ElastoDyn","RotSpeed")]   = {'vals': omega_init, 'group': 0}
    case_inputs[("ElastoDyn","BlPitch1")]   = {'vals': pitch_init, 'group': 0}
    case_inputs[("ElastoDyn","BlPitch2")]   = case_inputs[("ElastoDyn","BlPitch1")]
    case_inputs[("ElastoDyn","BlPitch3")]   = case_inputs[("ElastoDyn","BlPitch1")]

    case_inputs[("Fst","CompElast")]       = {'vals':[1,2], 'group':1}


    # Find the controller
    if platform.system() == 'Windows':
        sfx = 'dll'
    elif platform.system() == 'Darwin':
        sfx = 'dylib'
    else:
        sfx = 'so'
    path2dll = os.path.join(run_dir1, 'local','lib','libdiscon.'+sfx)

    # case_inputs[("ServoDyn","DLL_FileName")] = {'vals':[path2dll], 'group':0}

    # Generate the matrix of cases
    case_list, case_name_list = CaseGen_General(case_inputs, dir_matrix=fastBatch.FAST_runDirectory, namebase='iea15mw-nodal')

    fastBatch.case_list = case_list
    fastBatch.case_name_list = case_name_list

    # Run OpenFAST, either serially or sequentially
    if n_cores == 1:
        fastBatch.run_serial()
    else:
        fastBatch.run_multi(n_cores)
