import os, re, copy
import numpy as np
from functools import reduce
import operator
from openfast_io.FAST_vars_out import FstOutput

try:
    from rosco.toolbox.utilities import read_DISCON, load_from_txt
    from rosco.toolbox import turbine as ROSCO_turbine
    ROSCO = True
except:
    ROSCO = False


def readline_filterComments(f):
            read = True
            while read:
                line = f.readline().strip()
                if len(line)>0:
                    if line[0] != '!':
                        read = False
            return line

def fix_path(name):
    """ split a path, then reconstruct it using os.path.join """
    name = re.split("\\\|/", name)
    new = name[0]
    for i in range(1,len(name)):
        new = os.path.join(new, name[i])
    return new

def read_array(f,len,array_type=str):
    strings = re.split(',| ',f.readline().strip())
    while '' in strings:    # remove empties
        strings.remove('')

    arr = strings[:len]    # select len strings

    if array_type==str:
        arr = [ar.replace('"','') for ar in arr]  # remove quotes and commas
    elif array_type==float:
        arr = [float_read(ar) for ar in arr]
    elif array_type==int:
        arr = [int_read(ar) for ar in arr]
    elif array_type==bool:
        arr = [bool_read(ar) for ar in arr]
    else:
        raise Exception(f"read_array with type {str(array_type)} not currently supported")

    return arr

def bool_read(text):
    # convert true/false strings to boolean
    if 'default' in text.lower():
        return str(text)
    else:
        if text.lower() == 'true':
            return True
        else:
            return False

def float_read(text):
    # return float with error handing for "default" values
    if 'default' in text.lower():
        return str(text)
    else:
        try:
            return float(text)
        except:
            return str(text)

def int_read(text):
    # return int with error handing for "default" values
    if 'default' in text.lower():
        return str(text)
    else:
        try:
            return int(text)
        except:
            return str(text)

class InputReader_OpenFAST(object):
    """ OpenFAST input file reader """

    def __init__(self):

        self.FAST_InputFile = None   # FAST input file (ext=.fst)
        self.FAST_directory = None   # Path to fst directory files
        self.path2dll       = None   # Path to dll file
        self.fst_vt         = {}
        self.fst_vt['Fst']  = {}
        self.fst_vt['outlist']  = FstOutput
        self.fst_vt['ElastoDyn'] = {}
        self.fst_vt['ElastoDynBlade'] = {}
        self.fst_vt['ElastoDynTower'] = {}
        self.fst_vt['InflowWind'] = {}
        self.fst_vt['AeroDyn15'] = {}
        self.fst_vt['AeroDyn14'] = {}
        self.fst_vt['AeroDynBlade'] = {}
        self.fst_vt['AeroDynTower'] = {}
        self.fst_vt['AeroDynPolar'] = {}
        self.fst_vt['ServoDyn'] = {}
        self.fst_vt['DISCON_in'] = {}
        self.fst_vt['HydroDyn'] = {}
        self.fst_vt['MoorDyn'] = {}
        self.fst_vt['SubDyn'] = {}
        self.fst_vt['MAP'] = {}
        self.fst_vt['BeamDyn'] = {}
        self.fst_vt['BeamDynBlade'] = {}

    def set_outlist(self, vartree_head, channel_list):
        """ Loop through a list of output channel names, recursively set them to True in the nested outlist dict """

        # given a list of nested dictionary keys, return the dict at that point
        def get_dict(vartree, branch):
            return reduce(operator.getitem, branch, vartree_head)
        # given a list of nested dictionary keys, set the value of the dict at that point
        def set_dict(vartree, branch, val):
            get_dict(vartree, branch[:-1])[branch[-1]] = val
        # recursively loop through outlist dictionaries to set output channels
        def loop_dict(vartree, search_var, branch):
            for var in vartree.keys():
                branch_i = copy.copy(branch)
                branch_i.append(var)
                if type(vartree[var]) is dict:
                    loop_dict(vartree[var], search_var, branch_i)
                else:
                    if var == search_var:
                        set_dict(vartree_head, branch_i, True)

        # loop through outchannels on this line, loop through outlist dicts to set to True
        for var in channel_list:
            var = var.replace(' ', '')
            loop_dict(vartree_head, var, [])

    def read_MainInput(self):
        # Main FAST v8.16-v8.17 Input File
        # Currently no differences between FASTv8.16 and OpenFAST.
        fst_file = os.path.join(self.FAST_directory, self.FAST_InputFile)
        f = open(fst_file)

        # Header of .fst file
        f.readline()
        self.fst_vt['description'] = f.readline().rstrip()

        # Simulation Control (fst_sim_ctrl)
        f.readline()
        self.fst_vt['Fst']['Echo'] = bool_read(f.readline().split()[0])
        self.fst_vt['Fst']['AbortLevel'] = f.readline().split()[0][1:-1]
        self.fst_vt['Fst']['TMax'] = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['DT']  = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['InterpOrder']  = int(f.readline().split()[0])
        self.fst_vt['Fst']['NumCrctn']  = int(f.readline().split()[0])
        self.fst_vt['Fst']['DT_UJac']  = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['UJacSclFact']  = float_read(f.readline().split()[0])

        # Feature Switches and Flags (ftr_swtchs_flgs)
        f.readline()
        self.fst_vt['Fst']['CompElast'] = int(f.readline().split()[0])
        self.fst_vt['Fst']['CompInflow'] = int(f.readline().split()[0])
        self.fst_vt['Fst']['CompAero'] = int(f.readline().split()[0])
        self.fst_vt['Fst']['CompServo'] = int(f.readline().split()[0])
        self.fst_vt['Fst']['CompHydro'] = int(f.readline().split()[0])
        self.fst_vt['Fst']['CompSub'] = int(f.readline().split()[0])
        self.fst_vt['Fst']['CompMooring'] = int(f.readline().split()[0])
        self.fst_vt['Fst']['CompIce'] = int(f.readline().split()[0])
        self.fst_vt['Fst']['MHK'] = int(f.readline().split()[0])

        # Environmental conditions
        f.readline()
        self.fst_vt['Fst']['Gravity']   = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['AirDens']   = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['WtrDens']   = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['KinVisc']   = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['SpdSound']  = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['Patm']      = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['Pvap']      = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['WtrDpth']   = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['MSL2SWL']   = float_read(f.readline().split()[0])

        # Input Files (input_files)
        f.readline()
        self.fst_vt['Fst']['EDFile'] = f.readline().split()[0][1:-1]
        self.fst_vt['Fst']['BDBldFile(1)'] = f.readline().split()[0][1:-1]
        self.fst_vt['Fst']['BDBldFile(2)'] = f.readline().split()[0][1:-1]
        self.fst_vt['Fst']['BDBldFile(3)'] = f.readline().split()[0][1:-1]
        self.fst_vt['Fst']['InflowFile'] = f.readline().split()[0][1:-1]
        self.fst_vt['Fst']['AeroFile'] = f.readline().split()[0][1:-1]
        self.fst_vt['Fst']['ServoFile'] = f.readline().split()[0][1:-1]
        self.fst_vt['Fst']['HydroFile'] = f.readline().split()[0][1:-1]
        self.fst_vt['Fst']['SubFile'] = f.readline().split()[0][1:-1]
        self.fst_vt['Fst']['MooringFile'] = f.readline().split()[0][1:-1]
        self.fst_vt['Fst']['IceFile'] = f.readline().split()[0][1:-1]

        # FAST Output Parameters (fst_output_params)
        f.readline()
        self.fst_vt['Fst']['SumPrint'] = bool_read(f.readline().split()[0])
        self.fst_vt['Fst']['SttsTime'] = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['ChkptTime'] = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['DT_Out'] = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['TStart'] = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['OutFileFmt'] = int(f.readline().split()[0])
        self.fst_vt['Fst']['TabDelim'] = bool_read(f.readline().split()[0])
        self.fst_vt['Fst']['OutFmt'] = f.readline().split()[0][1:-1]

        # Fst
        f.readline()
        self.fst_vt['Fst']['Linearize']  = f.readline().split()[0]
        self.fst_vt['Fst']['CalcSteady'] = f.readline().split()[0]
        self.fst_vt['Fst']['TrimCase']   = f.readline().split()[0]
        self.fst_vt['Fst']['TrimTol']    = f.readline().split()[0]
        self.fst_vt['Fst']['TrimGain']   = f.readline().split()[0]
        self.fst_vt['Fst']['Twr_Kdmp']   = f.readline().split()[0]
        self.fst_vt['Fst']['Bld_Kdmp']   = f.readline().split()[0]
        self.fst_vt['Fst']['NLinTimes']  = f.readline().split()[0]
        self.fst_vt['Fst']['LinTimes']   = re.findall(r'[^,\s]+', f.readline())[0:2]
        self.fst_vt['Fst']['LinInputs']  = f.readline().split()[0]
        self.fst_vt['Fst']['LinOutputs'] = f.readline().split()[0]
        self.fst_vt['Fst']['LinOutJac']  = f.readline().split()[0]
        self.fst_vt['Fst']['LinOutMod']  = f.readline().split()[0]

        # Visualization ()
        f.readline()
        self.fst_vt['Fst']['WrVTK'] = int(f.readline().split()[0])
        self.fst_vt['Fst']['VTK_type'] = int(f.readline().split()[0])
        self.fst_vt['Fst']['VTK_fields'] = bool_read(f.readline().split()[0])
        self.fst_vt['Fst']['VTK_fps'] = float_read(f.readline().split()[0])
        
        f.close()
        
        # File paths
        self.fst_vt['Fst']['EDFile_path']       = os.path.split(self.fst_vt['Fst']['EDFile'])[0]
        self.fst_vt['Fst']['BDBldFile(1_path)'] = os.path.split(self.fst_vt['Fst']['BDBldFile(1)'])[0]
        self.fst_vt['Fst']['BDBldFile(2_path)'] = os.path.split(self.fst_vt['Fst']['BDBldFile(2)'])[0]
        self.fst_vt['Fst']['BDBldFile(3_path)'] = os.path.split(self.fst_vt['Fst']['BDBldFile(3)'])[0]
        self.fst_vt['Fst']['InflowFile_path']   = os.path.split(self.fst_vt['Fst']['InflowFile'])[0]
        self.fst_vt['Fst']['AeroFile_path']     = os.path.split(self.fst_vt['Fst']['AeroFile'])[0]
        self.fst_vt['Fst']['ServoFile_path']    = os.path.split(self.fst_vt['Fst']['ServoFile'])[0]
        self.fst_vt['Fst']['HydroFile_path']    = os.path.split(self.fst_vt['Fst']['HydroFile'])[0]
        self.fst_vt['Fst']['SubFile_path']      = os.path.split(self.fst_vt['Fst']['SubFile'])[0]
        self.fst_vt['Fst']['MooringFile_path']  = os.path.split(self.fst_vt['Fst']['MooringFile'])[0]
        self.fst_vt['Fst']['IceFile_path']      = os.path.split(self.fst_vt['Fst']['IceFile'])[0]

    def read_ElastoDyn(self, ed_file):
        # ElastoDyn v1.03 Input File
        # Currently no differences between FASTv8.16 and OpenFAST.

        f = open(ed_file)

        f.readline()
        f.readline()

        # Simulation Control (ed_sim_ctrl)
        f.readline()
        self.fst_vt['ElastoDyn']['Echo'] = bool_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['Method']  = int(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['DT'] = float_read(f.readline().split()[0])

        # Degrees of Freedom (dof)
        f.readline()
        self.fst_vt['ElastoDyn']['FlapDOF1'] = bool_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['FlapDOF2'] = bool_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['EdgeDOF'] = bool_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['TeetDOF'] = bool_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['DrTrDOF'] = bool_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['GenDOF'] = bool_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['YawDOF'] = bool_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['TwFADOF1'] = bool_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['TwFADOF2'] = bool_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['TwSSDOF1'] = bool_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['TwSSDOF2'] = bool_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['PtfmSgDOF'] = bool_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['PtfmSwDOF'] = bool_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['PtfmHvDOF'] = bool_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['PtfmRDOF'] = bool_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['PtfmPDOF'] = bool_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['PtfmYDOF'] = bool_read(f.readline().split()[0])

        # Initial Conditions (init_conds)
        f.readline()
        self.fst_vt['ElastoDyn']['OoPDefl']    = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['IPDefl']     = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['BlPitch1']   = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['BlPitch2']   = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['BlPitch3']   = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['TeetDefl']   = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['Azimuth']    = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['RotSpeed']   = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['NacYaw']     = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['TTDspFA']    = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['TTDspSS']    = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['PtfmSurge']  = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['PtfmSway']   = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['PtfmHeave']  = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['PtfmRoll']   = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['PtfmPitch']  = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['PtfmYaw']    = float_read(f.readline().split()[0])


        # Turbine Configuration (turb_config)
        f.readline()
        self.fst_vt['ElastoDyn']['NumBl']      = int(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['TipRad']     = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['HubRad']     = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['PreCone(1)']   = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['PreCone(2)']   = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['PreCone(3)']   = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['HubCM']      = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['UndSling']   = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['Delta3']     = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['AzimB1Up']   = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['OverHang']   = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['ShftGagL']   = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['ShftTilt']   = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['NacCMxn']    = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['NacCMyn']    = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['NacCMzn']    = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['NcIMUxn']    = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['NcIMUyn']    = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['NcIMUzn']    = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['Twr2Shft']   = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['TowerHt']    = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['TowerBsHt']  = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['PtfmCMxt']   = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['PtfmCMyt']   = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['PtfmCMzt']   = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['PtfmRefzt']  = float_read(f.readline().split()[0])

        # Mass and Inertia (mass_inertia)
        f.readline()
        self.fst_vt['ElastoDyn']['TipMass(1)']   = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['TipMass(2)']   = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['TipMass(3)']   = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['HubMass']    = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['HubIner']    = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['GenIner']    = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['NacMass']    = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['NacYIner']   = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['YawBrMass']  = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['PtfmMass']   = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['PtfmRIner']  = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['PtfmPIner']  = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['PtfmYIner']  = float_read(f.readline().split()[0])

        # ElastoDyn Blade (blade_struc)
        f.readline()
        self.fst_vt['ElastoDyn']['BldNodes'] = int(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['BldFile1'] = f.readline().split()[0][1:-1]
        self.fst_vt['ElastoDyn']['BldFile2'] = f.readline().split()[0][1:-1]
        self.fst_vt['ElastoDyn']['BldFile3'] = f.readline().split()[0][1:-1]

        # Rotor-Teeter (rotor_teeter)
        f.readline()
        self.fst_vt['ElastoDyn']['TeetMod']  = int(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['TeetDmpP'] = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['TeetDmp']  = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['TeetCDmp'] = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['TeetSStP'] = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['TeetHStP'] = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['TeetSSSp'] = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['TeetHSSp'] = float_read(f.readline().split()[0])

        # Drivetrain (drivetrain)
        f.readline()
        self.fst_vt['ElastoDyn']['GBoxEff']  = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['GBRatio']  = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['DTTorSpr'] = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['DTTorDmp'] = float_read(f.readline().split()[0])

        # Furling (furling)
        f.readline()
        self.fst_vt['ElastoDyn']['Furling'] = bool_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['FurlFile'] = f.readline().split()[0][1:-1]

        # Tower (tower)
        f.readline()
        self.fst_vt['ElastoDyn']['TwrNodes'] = int(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['TwrFile'] = f.readline().split()[0][1:-1]

        # ED Output Parameters (ed_out_params)
        f.readline()
        self.fst_vt['ElastoDyn']['SumPrint'] = bool_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['OutFile']  = int(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['TabDelim'] = bool_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['OutFmt']   = f.readline().split()[0][1:-1]
        self.fst_vt['ElastoDyn']['TStart']   = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['DecFact']  = int(f.readline().split()[0])
        self.fst_vt['ElastoDyn']['NTwGages'] = int(f.readline().split()[0])
        if self.fst_vt['ElastoDyn']['NTwGages'] != 0: #loop over elements if there are gauges to be added, otherwise assign directly
            self.fst_vt['ElastoDyn']['TwrGagNd'] = f.readline().strip().split()[:self.fst_vt['ElastoDyn']['NTwGages']]
            for i, bldgag in enumerate(self.fst_vt['ElastoDyn']['TwrGagNd']):
                self.fst_vt['ElastoDyn']['TwrGagNd'][i] = int(bldgag.strip(','))
        else:
            self.fst_vt['ElastoDyn']['TwrGagNd'] = 0
            f.readline()
        self.fst_vt['ElastoDyn']['NBlGages'] = int(f.readline().split()[0])
        if self.fst_vt['ElastoDyn']['NBlGages'] != 0:
            self.fst_vt['ElastoDyn']['BldGagNd'] = f.readline().strip().split()[:self.fst_vt['ElastoDyn']['NBlGages']]
            for i, bldgag in enumerate(self.fst_vt['ElastoDyn']['BldGagNd']):
                self.fst_vt['ElastoDyn']['BldGagNd'][i] = int(bldgag.strip(','))
        else:
            self.fst_vt['ElastoDyn']['BldGagNd'] = 0
            f.readline()

        # Loop through output channel lines
        f.readline()
        data = f.readline()
        if data != '':
            while data.split()[0] != 'END':
                channels = data.split('"')
                channel_list = channels[1].split(',')
                self.set_outlist(self.fst_vt['outlist']['ElastoDyn'], channel_list)

                data = f.readline()
            
        # ElastoDyn optional outlist
        try:
            f.readline()
            self.fst_vt['ElastoDyn']['BldNd_BladesOut']    = int(f.readline().split()[0])
            self.fst_vt['ElastoDyn']['BldNd_BlOutNd']    = f.readline().split()[0]

            f.readline()
            data =  f.readline()
            while data.split()[0] != 'END':
                if data.find('"')>=0:
                    channels = data.split('"')
                    opt_channel_list = channels[1].split(',')
                else:
                    row_string = data.split(',')
                    if len(row_string)==1:
                        opt_channel_list = row_string[0].split('\n')[0]
                    else:
                        opt_channel_list = row_string
                self.set_outlist(self.fst_vt['outlist']['ElastoDyn_Nodes'], opt_channel_list)
                data = f.readline()
        except:
            # The optinal outlist does not exist.
            None

        f.close()

    def read_ElastoDynBlade(self, blade_file):
        # ElastoDyn v1.00 Blade Input File
        # Currently no differences between FASTv8.16 and OpenFAST.

        f = open(blade_file)
        # print blade_file
        f.readline()
        f.readline()
        f.readline()
        
        # Blade Parameters
        self.fst_vt['ElastoDynBlade']['NBlInpSt'] = int(f.readline().split()[0])
        self.fst_vt['ElastoDynBlade']['BldFlDmp1'] = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDynBlade']['BldFlDmp2'] = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDynBlade']['BldEdDmp1'] = float_read(f.readline().split()[0])
        
        # Blade Adjustment Factors
        f.readline()
        self.fst_vt['ElastoDynBlade']['FlStTunr1'] = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDynBlade']['FlStTunr2'] = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDynBlade']['AdjBlMs'] = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDynBlade']['AdjFlSt'] = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDynBlade']['AdjEdSt'] = float_read(f.readline().split()[0])
        
        # Distrilbuted Blade Properties
        f.readline()
        f.readline()
        f.readline()
        self.fst_vt['ElastoDynBlade']['BlFract'] = [None] * self.fst_vt['ElastoDynBlade']['NBlInpSt']
        self.fst_vt['ElastoDynBlade']['PitchAxis'] = [None] * self.fst_vt['ElastoDynBlade']['NBlInpSt']
        self.fst_vt['ElastoDynBlade']['StrcTwst'] = [None] * self.fst_vt['ElastoDynBlade']['NBlInpSt']
        self.fst_vt['ElastoDynBlade']['BMassDen'] = [None] * self.fst_vt['ElastoDynBlade']['NBlInpSt']
        self.fst_vt['ElastoDynBlade']['FlpStff'] = [None] * self.fst_vt['ElastoDynBlade']['NBlInpSt']
        self.fst_vt['ElastoDynBlade']['EdgStff'] = [None] * self.fst_vt['ElastoDynBlade']['NBlInpSt']

        for i in range(self.fst_vt['ElastoDynBlade']['NBlInpSt']):
            data = f.readline().split()          
            self.fst_vt['ElastoDynBlade']['BlFract'][i]  = float_read(data[0])
            self.fst_vt['ElastoDynBlade']['PitchAxis'][i]  = float_read(data[1])
            self.fst_vt['ElastoDynBlade']['StrcTwst'][i]  = float_read(data[2])
            self.fst_vt['ElastoDynBlade']['BMassDen'][i]  = float_read(data[3])
            self.fst_vt['ElastoDynBlade']['FlpStff'][i]  = float_read(data[4])
            self.fst_vt['ElastoDynBlade']['EdgStff'][i]  = float_read(data[5])

        f.readline()
        self.fst_vt['ElastoDynBlade']['BldFl1Sh'] = [None] * 5
        self.fst_vt['ElastoDynBlade']['BldFl2Sh'] = [None] * 5        
        self.fst_vt['ElastoDynBlade']['BldEdgSh'] = [None] * 5
        for i in range(5):
            self.fst_vt['ElastoDynBlade']['BldFl1Sh'][i]  = float_read(f.readline().split()[0])
        for i in range(5):
            self.fst_vt['ElastoDynBlade']['BldFl2Sh'][i]  = float_read(f.readline().split()[0])            
        for i in range(5):
            self.fst_vt['ElastoDynBlade']['BldEdgSh'][i]  = float_read(f.readline().split()[0])        

        f.close()

    def read_ElastoDynTower(self, tower_file):
        # ElastoDyn v1.00 Tower Input Files
        # Currently no differences between FASTv8.16 and OpenFAST.
        
        f = open(tower_file)

        f.readline()
        f.readline()

        # General Tower Paramters
        f.readline()
        self.fst_vt['ElastoDynTower']['NTwInpSt'] = int(f.readline().split()[0])
        self.fst_vt['ElastoDynTower']['TwrFADmp1'] = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDynTower']['TwrFADmp2'] = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDynTower']['TwrSSDmp1'] = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDynTower']['TwrSSDmp2'] = float_read(f.readline().split()[0])
    
        # Tower Adjustment Factors
        f.readline()
        self.fst_vt['ElastoDynTower']['FAStTunr1'] = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDynTower']['FAStTunr2'] = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDynTower']['SSStTunr1'] = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDynTower']['SSStTunr2'] = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDynTower']['AdjTwMa'] = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDynTower']['AdjFASt'] = float_read(f.readline().split()[0])
        self.fst_vt['ElastoDynTower']['AdjSSSt'] = float_read(f.readline().split()[0])
     
        # Distributed Tower Properties   
        f.readline()
        f.readline()
        f.readline()
        self.fst_vt['ElastoDynTower']['HtFract'] = [None] * self.fst_vt['ElastoDynTower']['NTwInpSt']
        self.fst_vt['ElastoDynTower']['TMassDen'] = [None] * self.fst_vt['ElastoDynTower']['NTwInpSt']
        self.fst_vt['ElastoDynTower']['TwFAStif'] = [None] * self.fst_vt['ElastoDynTower']['NTwInpSt']
        self.fst_vt['ElastoDynTower']['TwSSStif'] = [None] * self.fst_vt['ElastoDynTower']['NTwInpSt']

        for i in range(self.fst_vt['ElastoDynTower']['NTwInpSt']):
            data = f.readline().split()
            self.fst_vt['ElastoDynTower']['HtFract'][i]  = float_read(data[0])
            self.fst_vt['ElastoDynTower']['TMassDen'][i]  = float_read(data[1])
            self.fst_vt['ElastoDynTower']['TwFAStif'][i]  = float_read(data[2])
            self.fst_vt['ElastoDynTower']['TwSSStif'][i]  = float_read(data[3])       
        
        # Tower Mode Shapes
        f.readline()
        self.fst_vt['ElastoDynTower']['TwFAM1Sh'] = [None] * 5
        self.fst_vt['ElastoDynTower']['TwFAM2Sh'] = [None] * 5
        for i in range(5):
            self.fst_vt['ElastoDynTower']['TwFAM1Sh'][i]  = float_read(f.readline().split()[0])
        for i in range(5):
            self.fst_vt['ElastoDynTower']['TwFAM2Sh'][i]  = float_read(f.readline().split()[0])        
        f.readline()
        self.fst_vt['ElastoDynTower']['TwSSM1Sh'] = [None] * 5
        self.fst_vt['ElastoDynTower']['TwSSM2Sh'] = [None] * 5          
        for i in range(5):
            self.fst_vt['ElastoDynTower']['TwSSM1Sh'][i]  = float_read(f.readline().split()[0])
        for i in range(5):
            self.fst_vt['ElastoDynTower']['TwSSM2Sh'][i]  = float_read(f.readline().split()[0]) 

        f.close()

    def read_BeamDyn(self, bd_file):
        # BeamDyn Input File
        f = open(bd_file)
        f.readline()
        f.readline()
        f.readline()
        # ---------------------- SIMULATION CONTROL --------------------------------------
        self.fst_vt['BeamDyn']['Echo']             = bool_read(f.readline().split()[0])
        self.fst_vt['BeamDyn']['QuasiStaticInit']  = bool_read(f.readline().split()[0])
        self.fst_vt['BeamDyn']['rhoinf']           = float_read(f.readline().split()[0])
        self.fst_vt['BeamDyn']['quadrature']       = int_read(f.readline().split()[0])
        self.fst_vt['BeamDyn']['refine']           = int_read(f.readline().split()[0])
        self.fst_vt['BeamDyn']['n_fact']           = int_read(f.readline().split()[0])
        self.fst_vt['BeamDyn']['DTBeam']           = float_read(f.readline().split()[0])
        self.fst_vt['BeamDyn']['load_retries']     = int_read(f.readline().split()[0])
        self.fst_vt['BeamDyn']['NRMax']            = int_read(f.readline().split()[0])
        self.fst_vt['BeamDyn']['stop_tol']         = float_read(f.readline().split()[0])
        self.fst_vt['BeamDyn']['tngt_stf_fd']      = bool_read(f.readline().split()[0])
        self.fst_vt['BeamDyn']['tngt_stf_comp']    = bool_read(f.readline().split()[0])
        self.fst_vt['BeamDyn']['tngt_stf_pert']    = float_read(f.readline().split()[0])
        self.fst_vt['BeamDyn']['tngt_stf_difftol'] = float_read(f.readline().split()[0])
        self.fst_vt['BeamDyn']['RotStates']        = bool_read(f.readline().split()[0])
        f.readline()
        #---------------------- GEOMETRY PARAMETER --------------------------------------
        self.fst_vt['BeamDyn']['member_total']     = int_read(f.readline().split()[0])
        self.fst_vt['BeamDyn']['kp_total']         = int_read(f.readline().split()[0])
        self.fst_vt['BeamDyn']['members']          = []
        for i in range(self.fst_vt['BeamDyn']['member_total']):
            ln = f.readline().split()
            n_pts_i                   = int(ln[1])
            member_i                  = {}
            member_i['kp_xr']         = [None]*n_pts_i
            member_i['kp_yr']         = [None]*n_pts_i
            member_i['kp_zr']         = [None]*n_pts_i
            member_i['initial_twist'] = [None]*n_pts_i
            f.readline()
            f.readline()
            for j in range(n_pts_i):
                ln = f.readline().split()
                member_i['kp_xr'][j]          = float(ln[0])
                member_i['kp_yr'][j]          = float(ln[1])
                member_i['kp_zr'][j]          = float(ln[2])
                member_i['initial_twist'][j]  = float(ln[3])

            self.fst_vt['BeamDyn']['members'].append(member_i)
        #---------------------- MESH PARAMETER ------------------------------------------
        f.readline()
        self.fst_vt['BeamDyn']['order_elem']  = int_read(f.readline().split()[0])
        #---------------------- MATERIAL PARAMETER --------------------------------------
        f.readline()
        self.fst_vt['BeamDyn']['BldFile']     = f.readline().split()[0].replace('"','').replace("'",'')
        #---------------------- PITCH ACTUATOR PARAMETERS -------------------------------
        f.readline()
        self.fst_vt['BeamDyn']['UsePitchAct'] = bool_read(f.readline().split()[0])
        self.fst_vt['BeamDyn']['PitchJ']      = float_read(f.readline().split()[0])
        self.fst_vt['BeamDyn']['PitchK']      = float_read(f.readline().split()[0])
        self.fst_vt['BeamDyn']['PitchC']      = float_read(f.readline().split()[0])
        #---------------------- OUTPUTS -------------------------------------------------
        f.readline()
        self.fst_vt['BeamDyn']['SumPrint']    = bool_read(f.readline().split()[0])
        self.fst_vt['BeamDyn']['OutFmt']      = f.readline().split()[0][1:-1]
        self.fst_vt['BeamDyn']['NNodeOuts']   = int_read(f.readline().split()[0])
        self.fst_vt['BeamDyn']['OutNd']       = [idx.strip() for idx in f.readline().split('NNodeOuts')[0].split(',')]
        # BeamDyn Outlist
        f.readline()
        data = f.readline()
        while data.split()[0] != 'END':
            channels = data.split('"')
            channel_list = channels[1].split(',')
            self.set_outlist(self.fst_vt['outlist']['BeamDyn'], channel_list)
            data = f.readline()
            
        # BeamDyn optional outlist
        try:
            f.readline()
            # self.fst_vt['BeamDyn']['BldNd_BladesOut']    = int(f.readline().split()[0])
            self.fst_vt['BeamDyn']['BldNd_BlOutNd']    = f.readline().split()[0]

            f.readline()
            data =  f.readline()
            while data.split()[0] != 'END':
                if data.find('"')>=0:
                    channels = data.split('"')
                    opt_channel_list = channels[1].split(',')
                else:
                    row_string = data.split(',')
                    if len(row_string)==1:
                        opt_channel_list = row_string[0].split('\n')[0]
                    else:
                        opt_channel_list = row_string
                self.set_outlist(self.fst_vt['outlist']['BeamDyn_Nodes'], opt_channel_list)
                data = f.readline()
        except:
            # The optinal outlist does not exist.
            None

        f.close()

        beamdyn_blade_file = os.path.join(os.path.dirname(bd_file), self.fst_vt['BeamDyn']['BldFile'])
        self.read_BeamDynBlade(beamdyn_blade_file)

    def read_BeamDynBlade(self, beamdyn_blade_file):
        # BeamDyn Blade

        f = open(beamdyn_blade_file)
        
        f.readline()
        f.readline()
        f.readline()
        #---------------------- BLADE PARAMETERS --------------------------------------
        self.fst_vt['BeamDynBlade']['station_total'] = int_read(f.readline().split()[0])
        self.fst_vt['BeamDynBlade']['damp_type']     = int_read(f.readline().split()[0])
        f.readline()
        f.readline()
        f.readline()
        #---------------------- DAMPING COEFFICIENT------------------------------------
        ln = f.readline().split()
        self.fst_vt['BeamDynBlade']['mu1']           = float(ln[0])
        self.fst_vt['BeamDynBlade']['mu2']           = float(ln[1])
        self.fst_vt['BeamDynBlade']['mu3']           = float(ln[2])
        self.fst_vt['BeamDynBlade']['mu4']           = float(ln[3])
        self.fst_vt['BeamDynBlade']['mu5']           = float(ln[4])
        self.fst_vt['BeamDynBlade']['mu6']           = float(ln[5])
        f.readline()
        #---------------------- DISTRIBUTED PROPERTIES---------------------------------
        
        self.fst_vt['BeamDynBlade']['radial_stations'] = np.zeros((self.fst_vt['BeamDynBlade']['station_total']))
        self.fst_vt['BeamDynBlade']['beam_stiff']      = np.zeros((self.fst_vt['BeamDynBlade']['station_total'], 6, 6))
        self.fst_vt['BeamDynBlade']['beam_inertia']    = np.zeros((self.fst_vt['BeamDynBlade']['station_total'], 6, 6))
        for i in range(self.fst_vt['BeamDynBlade']['station_total']):
            self.fst_vt['BeamDynBlade']['radial_stations'][i]  = float_read(f.readline().split()[0])
            for j in range(6):
                self.fst_vt['BeamDynBlade']['beam_stiff'][i,j,:] = np.array([float(val) for val in f.readline().strip().split()])
            f.readline()
            for j in range(6):
                self.fst_vt['BeamDynBlade']['beam_inertia'][i,j,:] = np.array([float(val) for val in f.readline().strip().split()])
            f.readline()

        f.close()

    def read_InflowWind(self):
        # InflowWind v3.01
        # Currently no differences between FASTv8.16 and OpenFAST.
        inflow_file = os.path.normpath(os.path.join(self.FAST_directory, self.fst_vt['Fst']['InflowFile']))
        f = open(inflow_file)
        
        f.readline()
        f.readline()
        f.readline()

        # Inflow wind header parameters (inflow_wind)
        self.fst_vt['InflowWind']['Echo']           = bool_read(f.readline().split()[0])
        self.fst_vt['InflowWind']['WindType']       = int(f.readline().split()[0])
        self.fst_vt['InflowWind']['PropagationDir'] = float_read(f.readline().split()[0])
        self.fst_vt['InflowWind']['VFlowAng']       = float_read(f.readline().split()[0])
        self.fst_vt['InflowWind']['VelInterpCubic'] = bool_read(f.readline().split()[0])
        self.fst_vt['InflowWind']['NWindVel']       = int(f.readline().split()[0])
        self.fst_vt['InflowWind']['WindVxiList']    = float_read(f.readline().split()[0])
        self.fst_vt['InflowWind']['WindVyiList']    = float_read(f.readline().split()[0])
        self.fst_vt['InflowWind']['WindVziList']    = float_read(f.readline().split()[0])

        # Parameters for Steady Wind Conditions [used only for WindType = 1] (steady_wind_params)
        f.readline()
        self.fst_vt['InflowWind']['HWindSpeed'] = float_read(f.readline().split()[0])
        self.fst_vt['InflowWind']['RefHt'] = float_read(f.readline().split()[0])
        self.fst_vt['InflowWind']['PLexp'] = float_read(f.readline().split()[0])

        # Parameters for Uniform wind file   [used only for WindType = 2] (uniform_wind_params)
        f.readline()
        self.fst_vt['InflowWind']['Filename_Uni'] = os.path.join(os.path.split(inflow_file)[0], f.readline().split()[0][1:-1])
        self.fst_vt['InflowWind']['RefHt_Uni'] = float_read(f.readline().split()[0])
        self.fst_vt['InflowWind']['RefLength'] = float_read(f.readline().split()[0])

        # Parameters for Binary TurbSim Full-Field files   [used only for WindType = 3] (turbsim_wind_params)
        f.readline()
        self.fst_vt['InflowWind']['FileName_BTS'] = os.path.join(os.path.split(inflow_file)[0], f.readline().split()[0][1:-1])
        # Parameters for Binary Bladed-style Full-Field files   [used only for WindType = 4] (bladed_wind_params)
        f.readline()
        self.fst_vt['InflowWind']['FilenameRoot'] = f.readline().split()[0][1:-1]       
        self.fst_vt['InflowWind']['TowerFile'] = bool_read(f.readline().split()[0])

        # Parameters for HAWC-format binary files  [Only used with WindType = 5] (hawc_wind_params)
        f.readline()
        self.fst_vt['InflowWind']['FileName_u'] = os.path.normpath(os.path.join(os.path.split(inflow_file)[0], f.readline().split()[0][1:-1]))
        self.fst_vt['InflowWind']['FileName_v'] = os.path.normpath(os.path.join(os.path.split(inflow_file)[0], f.readline().split()[0][1:-1]))
        self.fst_vt['InflowWind']['FileName_w'] = os.path.normpath(os.path.join(os.path.split(inflow_file)[0], f.readline().split()[0][1:-1]))
        self.fst_vt['InflowWind']['nx']    = int(f.readline().split()[0])
        self.fst_vt['InflowWind']['ny']    = int(f.readline().split()[0])
        self.fst_vt['InflowWind']['nz']    = int(f.readline().split()[0])
        self.fst_vt['InflowWind']['dx']    = float_read(f.readline().split()[0])
        self.fst_vt['InflowWind']['dy']    = float_read(f.readline().split()[0])
        self.fst_vt['InflowWind']['dz']    = float_read(f.readline().split()[0])
        self.fst_vt['InflowWind']['RefHt_Hawc'] = float_read(f.readline().split()[0])

        # Scaling parameters for turbulence (still hawc_wind_params)
        f.readline()
        self.fst_vt['InflowWind']['ScaleMethod'] = int(f.readline().split()[0])
        self.fst_vt['InflowWind']['SFx']         = float_read(f.readline().split()[0])
        self.fst_vt['InflowWind']['SFy']         = float_read(f.readline().split()[0])
        self.fst_vt['InflowWind']['SFz']         = float_read(f.readline().split()[0])
        self.fst_vt['InflowWind']['SigmaFx']     = float_read(f.readline().split()[0])
        self.fst_vt['InflowWind']['SigmaFy']     = float_read(f.readline().split()[0])
        self.fst_vt['InflowWind']['SigmaFz']     = float_read(f.readline().split()[0])

        # Mean wind profile parameters (added to HAWC-format files) (still hawc_wind_params)
        f.readline()
        self.fst_vt['InflowWind']['URef']        = float_read(f.readline().split()[0])
        self.fst_vt['InflowWind']['WindProfile'] = int(f.readline().split()[0])
        self.fst_vt['InflowWind']['PLExp_Hawc']  = float_read(f.readline().split()[0])
        self.fst_vt['InflowWind']['Z0']          = float_read(f.readline().split()[0])
        self.fst_vt['InflowWind']['XOffset']     = float_read(f.readline().split()[0])

        # LIDAR parameters
        f.readline()
        self.fst_vt['InflowWind']['SensorType'] = int(f.readline().split()[0])
        self.fst_vt['InflowWind']['NumPulseGate'] = int(f.readline().split()[0])
        self.fst_vt['InflowWind']['PulseSpacing'] = float_read(f.readline().split()[0])
        self.fst_vt['InflowWind']['NumBeam'] = int(f.readline().split()[0])
        self.fst_vt['InflowWind']['FocalDistanceX'] = float_read(f.readline().split()[0])
        self.fst_vt['InflowWind']['FocalDistanceY'] = float_read(f.readline().split()[0])
        self.fst_vt['InflowWind']['FocalDistanceZ'] = float_read(f.readline().split()[0])
        self.fst_vt['InflowWind']['RotorApexOffsetPos'] = [idx.strip() for idx in f.readline().split('RotorApexOffsetPos')[0].split(',')]
        self.fst_vt['InflowWind']['URefLid'] = float_read(f.readline().split()[0])
        self.fst_vt['InflowWind']['MeasurementInterval'] = float_read(f.readline().split()[0])
        self.fst_vt['InflowWind']['LidRadialVel'] = bool_read(f.readline().split()[0])
        self.fst_vt['InflowWind']['ConsiderHubMotion'] = int(f.readline().split()[0])

        # Inflow Wind Output Parameters (inflow_out_params)
        f.readline()
        self.fst_vt['InflowWind']['SumPrint'] = bool_read(f.readline().split()[0])
        
        # InflowWind Outlist
        f.readline()
        data = f.readline()
        while data.split()[0] != 'END':
            if data.find('"')>=0:
                channels = data.split('"')
                channel_list = channels[1].split(',')
            else:
                row_string = data.split(',')
                if len(row_string)==1:
                    channel_list = row_string[0].split('\n')[0]
                else:
                    channel_list = row_string
            self.set_outlist(self.fst_vt['outlist']['InflowWind'], channel_list)
            data = f.readline()

        f.close()
                
    def read_AeroDyn15(self):
        # AeroDyn v15.03

        ad_file = os.path.join(self.FAST_directory, self.fst_vt['Fst']['AeroFile'])
        f = open(ad_file)

        # General Option
        f.readline()
        f.readline()
        f.readline()
        self.fst_vt['AeroDyn15']['Echo']          = bool_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['DTAero']        = float_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['WakeMod']       = int(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['AFAeroMod']     = int(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['TwrPotent']     = int(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['TwrShadow']     = int(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['TwrAero']       = bool_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['FrozenWake']    = bool_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['CavitCheck']    = bool_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['Buoyancy']      = bool_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['CompAA']        = bool_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['AA_InputFile']  = f.readline().split()[0]

        # Environmental Conditions
        f.readline()
        self.fst_vt['AeroDyn15']['AirDens']        = float_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['KinVisc']        = float_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['SpdSound']       = float_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['Patm']           = float_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['Pvap']           = float_read(f.readline().split()[0])
        #self.fst_vt['AeroDyn15']['FluidDepth']           = float_read(f.readline().split()[0])

        # Blade-Element/Momentum Theory Options
        f.readline()
        self.fst_vt['AeroDyn15']['SkewMod']               = int_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['SkewModFactor']         = float_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['TipLoss']               = bool_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['HubLoss']               = bool_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['TanInd']                = bool_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['AIDrag']                = bool_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['TIDrag']                = bool_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['IndToler']              = float_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['MaxIter']               = int(f.readline().split()[0])


        # Dynamic Blade-Element/Momentum Theory Options 
        f.readline()
        self.fst_vt['AeroDyn15']['DBEMT_Mod']          = int(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['tau1_const']         = float_read(f.readline().split()[0])

        # Olaf -- cOnvecting LAgrangian Filaments (Free Vortex Wake) Theory Options
        f.readline()
        self.fst_vt['AeroDyn15']['OLAFInputFileName']  = f.readline().split()[0]
        
        # Beddoes-Leishman Unsteady Airfoil Aerodynamics Options
        f.readline()
        self.fst_vt['AeroDyn15']['UAMod']                  = int(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['FLookup']                = bool_read(f.readline().split()[0])
        
        file_pos = f.tell()
        line = f.readline()
        if 'UAStartRad' in line:
            self.fst_vt['AeroDyn15']['UAStartRad']             = float_read(line.split()[0])
        else:
            f.seek(file_pos)
        
        file_pos = f.tell()
        line = f.readline()
        if 'UAEndRad' in line:
            self.fst_vt['AeroDyn15']['UAEndRad']             = float_read(line.split()[0])
        else:
            f.seek(file_pos)


        # Airfoil Information
        f.readline()
        self.fst_vt['AeroDyn15']['AFTabMod']         = int(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['InCol_Alfa']       = int(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['InCol_Cl']         = int(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['InCol_Cd']         = int(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['InCol_Cm']         = int(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['InCol_Cpmin']      = int(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['NumAFfiles']       = int(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['AFNames']          = [None] * self.fst_vt['AeroDyn15']['NumAFfiles']
        for i in range(self.fst_vt['AeroDyn15']['NumAFfiles']):
            af_filename = fix_path(f.readline().split()[0])[1:-1]
            self.fst_vt['AeroDyn15']['AFNames'][i]   = os.path.abspath(os.path.join(self.FAST_directory, self.fst_vt['Fst']['AeroFile_path'], af_filename))

        # Rotor/Blade Properties
        f.readline()
        self.fst_vt['AeroDyn15']['UseBlCm']        = bool_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['ADBlFile1']      = f.readline().split()[0][1:-1]
        self.fst_vt['AeroDyn15']['ADBlFile2']      = f.readline().split()[0][1:-1]
        self.fst_vt['AeroDyn15']['ADBlFile3']      = f.readline().split()[0][1:-1]

        # Hub, nacelle, and tail fin aerodynamics
        f.readline()
        self.fst_vt['AeroDyn15']['VolHub'] = float_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['HubCenBx'] = float_read(f.readline().split()[0])
        f.readline()
        self.fst_vt['AeroDyn15']['VolNac'] = float_read(f.readline().split()[0])
        # data = [float(val) for val in f.readline().split(',')]
        self.fst_vt['AeroDyn15']['NacCenB'] = [idx.strip() for idx in f.readline().split('NacCenB')[0].split(',')]
        f.readline()
        self.fst_vt['AeroDyn15']['TFinAero'] = bool_read(f.readline().split()[0])
        tfa_filename = fix_path(f.readline().split()[0])[1:-1]
        self.fst_vt['AeroDyn15']['TFinFile'] = os.path.abspath(os.path.join(self.FAST_directory, tfa_filename))

        # Tower Influence and Aerodynamics
        f.readline()
        self.fst_vt['AeroDyn15']['NumTwrNds']      = int(f.readline().split()[0])
        f.readline()
        f.readline()
        self.fst_vt['AeroDyn15']['TwrElev'] = [None]*self.fst_vt['AeroDyn15']['NumTwrNds']
        self.fst_vt['AeroDyn15']['TwrDiam'] = [None]*self.fst_vt['AeroDyn15']['NumTwrNds']
        self.fst_vt['AeroDyn15']['TwrCd'] = [None]*self.fst_vt['AeroDyn15']['NumTwrNds']
        self.fst_vt['AeroDyn15']['TwrTI'] = [None]*self.fst_vt['AeroDyn15']['NumTwrNds']
        self.fst_vt['AeroDyn15']['TwrCb'] = [None]*self.fst_vt['AeroDyn15']['NumTwrNds']
        for i in range(self.fst_vt['AeroDyn15']['NumTwrNds']):
            data = [float(val) for val in f.readline().split()]
            self.fst_vt['AeroDyn15']['TwrElev'][i] = data[0] 
            self.fst_vt['AeroDyn15']['TwrDiam'][i] = data[1] 
            self.fst_vt['AeroDyn15']['TwrCd'][i]   = data[2]
            self.fst_vt['AeroDyn15']['TwrTI'][i]   = data[3]
            self.fst_vt['AeroDyn15']['TwrCb'][i]   = data[4]

        # Outputs
        f.readline()
        self.fst_vt['AeroDyn15']['SumPrint']    = bool_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['NBlOuts']     = int(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['BlOutNd']     = [idx.strip() for idx in f.readline().split('BlOutNd')[0].split(',')]
        self.fst_vt['AeroDyn15']['NTwOuts']     = int(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['TwOutNd']     = [idx.strip() for idx in f.readline().split('TwOutNd')[0].split(',')]

        # AeroDyn15 Outlist
        f.readline()
        data = f.readline()
        while data.split()[0] != 'END':
            if data.find('"')>=0:
                channels = data.split('"')
                channel_list = channels[1].split(',')
            else:
                row_string = data.split(',')
                if len(row_string)==1:
                    channel_list = row_string[0].split('\n')[0]
                else:
                    channel_list = row_string
            self.set_outlist(self.fst_vt['outlist']['AeroDyn'], channel_list)
            data = f.readline()

        # AeroDyn15 optional outlist
        try:
            f.readline()
            self.fst_vt['AeroDyn15']['BldNd_BladesOut']    = int(f.readline().split()[0])
            self.fst_vt['AeroDyn15']['BldNd_BlOutNd']    = f.readline().split()[0]

            f.readline()
            data =  f.readline()
            while data.split()[0] != 'END':
                if data.find('"')>=0:
                    channels = data.split('"')
                    opt_channel_list = channels[1].split(',')
                else:
                    row_string = data.split(',')
                    if len(row_string)==1:
                        opt_channel_list = row_string[0].split('\n')[0]
                    else:
                        opt_channel_list = row_string
                self.set_outlist(self.fst_vt['outlist']['AeroDyn_Nodes'], opt_channel_list)
                data = f.readline()
        except:
            # The optinal outlist does not exist.
            None

        f.close()

        self.read_AeroDyn15Blade()
        self.read_AeroDyn15Polar()
        self.read_AeroDyn15Coord()
        olaf_filename = os.path.join(self.FAST_directory, self.fst_vt['AeroDyn15']['OLAFInputFileName'])
        if os.path.isfile(olaf_filename):
            self.read_AeroDyn15OLAF(olaf_filename)

    def read_AeroDyn15Blade(self):
        # AeroDyn v5.00 Blade Definition File

        ad_blade_file = os.path.join(self.FAST_directory, self.fst_vt['Fst']['AeroFile_path'], self.fst_vt['AeroDyn15']['ADBlFile1'])
        f = open(ad_blade_file)

        f.readline()
        f.readline()
        f.readline()
        # Blade Properties
        self.fst_vt['AeroDynBlade']['NumBlNds']       = int(f.readline().split()[0])
        f.readline()
        f.readline()
        self.fst_vt['AeroDynBlade']['BlSpn']          = [None]*self.fst_vt['AeroDynBlade']['NumBlNds']
        self.fst_vt['AeroDynBlade']['BlCrvAC']        = [None]*self.fst_vt['AeroDynBlade']['NumBlNds']
        self.fst_vt['AeroDynBlade']['BlSwpAC']        = [None]*self.fst_vt['AeroDynBlade']['NumBlNds']
        self.fst_vt['AeroDynBlade']['BlCrvAng']       = [None]*self.fst_vt['AeroDynBlade']['NumBlNds']
        self.fst_vt['AeroDynBlade']['BlTwist']        = [None]*self.fst_vt['AeroDynBlade']['NumBlNds']
        self.fst_vt['AeroDynBlade']['BlChord']        = [None]*self.fst_vt['AeroDynBlade']['NumBlNds']
        self.fst_vt['AeroDynBlade']['BlAFID']         = [None]*self.fst_vt['AeroDynBlade']['NumBlNds']
        for i in range(self.fst_vt['AeroDynBlade']['NumBlNds']):
            data = [float(val) for val in f.readline().split()]
            self.fst_vt['AeroDynBlade']['BlSpn'][i]   = data[0] 
            self.fst_vt['AeroDynBlade']['BlCrvAC'][i] = data[1] 
            self.fst_vt['AeroDynBlade']['BlSwpAC'][i] = data[2]
            self.fst_vt['AeroDynBlade']['BlCrvAng'][i]= data[3]
            self.fst_vt['AeroDynBlade']['BlTwist'][i] = data[4]
            self.fst_vt['AeroDynBlade']['BlChord'][i] = data[5]
            self.fst_vt['AeroDynBlade']['BlAFID'][i]  = data[6]
        
        f.close()

    def read_AeroDyn15Polar(self):
        # AirfoilInfo v1.01

        
        self.fst_vt['AeroDyn15']['af_data'] = [None]*self.fst_vt['AeroDyn15']['NumAFfiles']

        for afi, af_filename in enumerate(self.fst_vt['AeroDyn15']['AFNames']):
            f = open(af_filename)
            # print af_filename

            polar = {}

            polar['InterpOrd']      = int_read(readline_filterComments(f).split()[0])
            polar['NonDimArea']     = float_read(readline_filterComments(f).split()[0])
            polar['NumCoords']      = readline_filterComments(f).split()[0]
            polar['BL_file']        = readline_filterComments(f).split()[0]
            polar['NumTabs']        = int_read(readline_filterComments(f).split()[0])
            self.fst_vt['AeroDyn15']['af_data'][afi] = [None]*polar['NumTabs']

            for tab in range(polar['NumTabs']): # For multiple tables
                polar['Re']             = float_read(readline_filterComments(f).split()[0]) * 1.e+6
                polar['Ctrl']           = int_read(readline_filterComments(f).split()[0])
                polar['InclUAdata']     = bool_read(readline_filterComments(f).split()[0])

                # Unsteady Aero Data
                if polar['InclUAdata']:
                    polar['alpha0']     = float_read(readline_filterComments(f).split()[0])
                    polar['alpha1']     = float_read(readline_filterComments(f).split()[0])
                    polar['alpha2']     = float_read(readline_filterComments(f).split()[0])
                    # polar['alphaUpper']     = float_read(readline_filterComments(f).split()[0])
                    # polar['alphaLower']     = float_read(readline_filterComments(f).split()[0])
                    polar['eta_e']      = float_read(readline_filterComments(f).split()[0])
                    polar['C_nalpha']   = float_read(readline_filterComments(f).split()[0])
                    polar['T_f0']       = float_read(readline_filterComments(f).split()[0])
                    polar['T_V0']       = float_read(readline_filterComments(f).split()[0])
                    polar['T_p']        = float_read(readline_filterComments(f).split()[0])
                    polar['T_VL']       = float_read(readline_filterComments(f).split()[0])
                    polar['b1']         = float_read(readline_filterComments(f).split()[0])
                    polar['b2']         = float_read(readline_filterComments(f).split()[0])
                    polar['b5']         = float_read(readline_filterComments(f).split()[0])
                    polar['A1']         = float_read(readline_filterComments(f).split()[0])
                    polar['A2']         = float_read(readline_filterComments(f).split()[0])
                    polar['A5']         = float_read(readline_filterComments(f).split()[0])
                    polar['S1']         = float_read(readline_filterComments(f).split()[0])
                    polar['S2']         = float_read(readline_filterComments(f).split()[0])
                    polar['S3']         = float_read(readline_filterComments(f).split()[0])
                    polar['S4']         = float_read(readline_filterComments(f).split()[0])
                    polar['Cn1']        = float_read(readline_filterComments(f).split()[0])
                    polar['Cn2']        = float_read(readline_filterComments(f).split()[0])
                    polar['St_sh']      = float_read(readline_filterComments(f).split()[0])
                    polar['Cd0']        = float_read(readline_filterComments(f).split()[0])
                    polar['Cm0']        = float_read(readline_filterComments(f).split()[0])
                    polar['k0']         = float_read(readline_filterComments(f).split()[0])
                    polar['k1']         = float_read(readline_filterComments(f).split()[0])
                    polar['k2']         = float_read(readline_filterComments(f).split()[0])
                    polar['k3']         = float_read(readline_filterComments(f).split()[0])
                    polar['k1_hat']     = float_read(readline_filterComments(f).split()[0])
                    polar['x_cp_bar']   = float_read(readline_filterComments(f).split()[0])
                    polar['UACutout']   = float_read(readline_filterComments(f).split()[0])
                    polar['filtCutOff'] = float_read(readline_filterComments(f).split()[0])

                # Polar Data
                polar['NumAlf']         = int_read(readline_filterComments(f).split()[0])
                polar['Alpha']          = [None]*polar['NumAlf']
                polar['Cl']             = [None]*polar['NumAlf']
                polar['Cd']             = [None]*polar['NumAlf']
                polar['Cm']             = [None]*polar['NumAlf']
                polar['Cpmin']          = [None]*polar['NumAlf']
                for i in range(polar['NumAlf']):
                    data = [float(val) for val in readline_filterComments(f).split()]
                    if self.fst_vt['AeroDyn15']['InCol_Alfa'] > 0:
                        polar['Alpha'][i] = data[self.fst_vt['AeroDyn15']['InCol_Alfa']-1]
                    if self.fst_vt['AeroDyn15']['InCol_Cl'] > 0:
                        polar['Cl'][i]    = data[self.fst_vt['AeroDyn15']['InCol_Cl']-1]
                    if self.fst_vt['AeroDyn15']['InCol_Cd'] > 0:
                        polar['Cd'][i]    = data[self.fst_vt['AeroDyn15']['InCol_Cd']-1]
                    if self.fst_vt['AeroDyn15']['InCol_Cm'] > 0:
                        polar['Cm'][i]    = data[self.fst_vt['AeroDyn15']['InCol_Cm']-1]
                    if self.fst_vt['AeroDyn15']['InCol_Cpmin'] > 0:
                        polar['Cpmin'][i] = data[self.fst_vt['AeroDyn15']['InCol_Cpmin']-1]

                self.fst_vt['AeroDyn15']['af_data'][afi][tab] = copy.copy(polar) # For multiple tables
            
            f.close()

    def read_AeroDyn15Coord(self):
        
        self.fst_vt['AeroDyn15']['af_coord'] = []
        self.fst_vt['AeroDyn15']['ac']   = np.zeros(len(self.fst_vt['AeroDyn15']['AFNames']))

        for afi, af_filename in enumerate(self.fst_vt['AeroDyn15']['AFNames']):
            self.fst_vt['AeroDyn15']['af_coord'].append({})
            if not (self.fst_vt['AeroDyn15']['af_data'][afi][0]['NumCoords'] == 0 or self.fst_vt['AeroDyn15']['af_data'][afi][0]['NumCoords'] == '0'):
                coord_filename = af_filename[0:af_filename.rfind(os.sep)] + os.sep + self.fst_vt['AeroDyn15']['af_data'][afi][0]['NumCoords'][2:-1]
                f = open(coord_filename)
                n_coords = int_read(readline_filterComments(f).split()[0])
                x = np.zeros(n_coords)
                y = np.zeros(n_coords)
                f.readline()
                f.readline()
                f.readline()
                self.fst_vt['AeroDyn15']['ac'][afi] = float(f.readline().split()[0])
                f.readline()
                f.readline()
                f.readline()
                for j in range(n_coords - 1):
                    x[j], y[j] = f.readline().split()

                self.fst_vt['AeroDyn15']['af_coord'][afi]['x'] = x
                self.fst_vt['AeroDyn15']['af_coord'][afi]['y'] = y

                f.close()

    def read_AeroDyn15OLAF(self, olaf_filename):
        
        self.fst_vt['AeroDyn15']['OLAF'] = {}
        f = open(olaf_filename)
        f.readline()
        f.readline()
        f.readline()
        self.fst_vt['AeroDyn15']['OLAF']['IntMethod']       = int_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['OLAF']['DTfvw']           = float_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['OLAF']['FreeWakeStart']   = float_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['OLAF']['FullCircStart']   = float_read(f.readline().split()[0])
        f.readline()
        self.fst_vt['AeroDyn15']['OLAF']['CircSolvMethod']   = int_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['OLAF']['CircSolvConvCrit']    = float_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['OLAF']['CircSolvRelaxation']  = float_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['OLAF']['CircSolvMaxIter']     = int_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['OLAF']['PrescribedCircFile']  = f.readline().split()[0]
        f.readline()
        f.readline()
        f.readline()
        self.fst_vt['AeroDyn15']['OLAF']['nNWPanels']   = int_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['OLAF']['nNWPanelsFree'] = int_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['OLAF']['nFWPanels']  = int_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['OLAF']['nFWPanelsFree']  = int_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['OLAF']['FWShedVorticity'] = bool_read(f.readline().split()[0])
        f.readline()
        self.fst_vt['AeroDyn15']['OLAF']['DiffusionMethod'] = int_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['OLAF']['RegDeterMethod']  = int_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['OLAF']['RegFunction']     = int_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['OLAF']['WakeRegMethod']   = int_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['OLAF']['WakeRegFactor']   = float(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['OLAF']['WingRegFactor']   = float(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['OLAF']['CoreSpreadEddyVisc'] = int(f.readline().split()[0])
        f.readline()
        self.fst_vt['AeroDyn15']['OLAF']['TwrShadowOnWake'] = bool_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['OLAF']['ShearModel']      = int_read(f.readline().split()[0])
        f.readline()
        self.fst_vt['AeroDyn15']['OLAF']['VelocityMethod']  = int_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['OLAF']['TreeBranchFactor']= float_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['OLAF']['PartPerSegment']  = int_read(f.readline().split()[0])
        f.readline()
        f.readline()
        self.fst_vt['AeroDyn15']['OLAF']['WrVTk']       = int_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['OLAF']['nVTKBlades']  = int_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['OLAF']['VTKCoord']    = int_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['OLAF']['VTK_fps']     = float_read(f.readline().split()[0])
        self.fst_vt['AeroDyn15']['OLAF']['nGridOut']    = int_read(f.readline().split()[0])
        f.readline()
        f.close()

    def read_AeroDyn14(self):
        # AeroDyn v14.04

        ad_file = os.path.join(self.FAST_directory, self.fst_vt['Fst']['AeroFile'])
        f = open(ad_file)
        # AeroDyn file header (aerodyn)
        f.readline()
        f.readline()
        self.fst_vt['AeroDyn14']['StallMod'] = f.readline().split()[0]
        self.fst_vt['AeroDyn14']['UseCm'] = f.readline().split()[0]
        self.fst_vt['AeroDyn14']['InfModel'] = f.readline().split()[0]
        self.fst_vt['AeroDyn14']['IndModel'] = f.readline().split()[0]
        self.fst_vt['AeroDyn14']['AToler'] = float_read(f.readline().split()[0])
        self.fst_vt['AeroDyn14']['TLModel'] = f.readline().split()[0]
        self.fst_vt['AeroDyn14']['HLModel'] = f.readline().split()[0]
        self.fst_vt['AeroDyn14']['TwrShad'] = int(f.readline().split()[0])
        if self.fst_vt['AeroDyn14']['TwrShad'] > 0:
            self.fst_vt['AeroDyn14']['TwrPotent'] = bool_read(f.readline().split()[0])
            self.fst_vt['AeroDyn14']['TwrShadow'] = bool_read(f.readline().split()[0])
            self.fst_vt['AeroDyn14']['TwrFile'] = f.readline().split()[0].replace('"','').replace("'",'')
            self.fst_vt['AeroDyn14']['CalcTwrAero'] = bool_read(f.readline().split()[0])
        else:
            self.fst_vt['AeroDyn14']['ShadHWid'] = float_read(f.readline().split()[0])
            self.fst_vt['AeroDyn14']['T_Shad_Refpt'] = float_read(f.readline().split()[0])
        self.fst_vt['AeroDyn14']['AirDens'] = float_read(f.readline().split()[0])
        self.fst_vt['AeroDyn14']['KinVisc'] = float_read(f.readline().split()[0])
        self.fst_vt['AeroDyn14']['DTAero'] = float_read(f.readline().split()[0])

        # AeroDyn Blade Properties (blade_aero)
        self.fst_vt['AeroDyn14']['NumFoil'] = int(f.readline().split()[0])
        self.fst_vt['AeroDyn14']['FoilNm'] = [None] * self.fst_vt['AeroDyn14']['NumFoil']
        for i in range(self.fst_vt['AeroDyn14']['NumFoil']):
            af_filename = f.readline().split()[0]
            af_filename = fix_path(af_filename)
            self.fst_vt['AeroDyn14']['FoilNm'][i]  = af_filename[1:-1]
        
        self.fst_vt['AeroDynBlade']['BldNodes'] = int(f.readline().split()[0])
        f.readline()
        self.fst_vt['AeroDynBlade']['RNodes'] = [None] * self.fst_vt['AeroDynBlade']['BldNodes']
        self.fst_vt['AeroDynBlade']['AeroTwst'] = [None] * self.fst_vt['AeroDynBlade']['BldNodes']
        self.fst_vt['AeroDynBlade']['DRNodes'] = [None] * self.fst_vt['AeroDynBlade']['BldNodes']
        self.fst_vt['AeroDynBlade']['Chord'] = [None] * self.fst_vt['AeroDynBlade']['BldNodes']
        self.fst_vt['AeroDynBlade']['NFoil'] = [None] * self.fst_vt['AeroDynBlade']['BldNodes']
        self.fst_vt['AeroDynBlade']['PrnElm'] = [None] * self.fst_vt['AeroDynBlade']['BldNodes']       
        for i in range(self.fst_vt['AeroDynBlade']['BldNodes']):
            data = f.readline().split()
            self.fst_vt['AeroDynBlade']['RNodes'][i]  = float_read(data[0])
            self.fst_vt['AeroDynBlade']['AeroTwst'][i]  = float_read(data[1])
            self.fst_vt['AeroDynBlade']['DRNodes'][i]  = float_read(data[2])
            self.fst_vt['AeroDynBlade']['Chord'][i]  = float_read(data[3])
            self.fst_vt['AeroDynBlade']['NFoil'][i]  = int(data[4])
            self.fst_vt['AeroDynBlade']['PrnElm'][i]  = data[5]

        f.close()

        # create airfoil objects
        self.fst_vt['AeroDynBlade']['af_data'] = []
        for i in range(self.fst_vt['AeroDyn14']['NumFoil']):
             self.fst_vt['AeroDynBlade']['af_data'].append(self.read_AeroDyn14Polar(os.path.join(self.FAST_directory,self.fst_vt['AeroDyn14']['FoilNm'][i])))

        # tower
        if self.fst_vt['AeroDyn14']['TwrShad'] > 0:
            self.read_AeroDyn14Tower()

    def read_AeroDyn14Tower(self):
        # AeroDyn v14.04 Tower

        ad_tower_file = os.path.join(self.FAST_directory, self.fst_vt['AeroDyn14']['TwrFile'])
        f = open(ad_tower_file)

        f.readline()
        f.readline()
        self.fst_vt['AeroDynTower']['NTwrHt'] = int(f.readline().split()[0])
        self.fst_vt['AeroDynTower']['NTwrRe'] = int(f.readline().split()[0])
        self.fst_vt['AeroDynTower']['NTwrCD'] = int(f.readline().split()[0])
        self.fst_vt['AeroDynTower']['Tower_Wake_Constant'] = float_read(f.readline().split()[0])
        
        f.readline()
        f.readline()
        self.fst_vt['AeroDynTower']['TwrHtFr'] = [None]*self.fst_vt['AeroDynTower']['NTwrHt']
        self.fst_vt['AeroDynTower']['TwrWid'] = [None]*self.fst_vt['AeroDynTower']['NTwrHt']
        self.fst_vt['AeroDynTower']['NTwrCDCol'] = [None]*self.fst_vt['AeroDynTower']['NTwrHt']
        for i in range(self.fst_vt['AeroDynTower']['NTwrHt']):
            data = [float(val) for val in f.readline().split()]
            self.fst_vt['AeroDynTower']['TwrHtFr'][i]  = data[0] 
            self.fst_vt['AeroDynTower']['TwrWid'][i]  = data[1]
            self.fst_vt['AeroDynTower']['NTwrCDCol'][i]  = data[2]

        f.readline()
        f.readline()
        self.fst_vt['AeroDynTower']['TwrRe'] = [None]*self.fst_vt['AeroDynTower']['NTwrRe']
        self.fst_vt['AeroDynTower']['TwrCD'] = np.zeros((self.fst_vt['AeroDynTower']['NTwrRe'], self.fst_vt['AeroDynTower']['NTwrCD']))
        for i in range(self.fst_vt['AeroDynTower']['NTwrRe']):
            data = [float(val) for val in f.readline().split()]
            self.fst_vt['AeroDynTower']['TwrRe'][i]  = data[0]
            self.fst_vt['AeroDynTower']['TwrCD'][i,:]  = data[1:]

        f.close()

    def read_AeroDyn14Polar(self, aerodynFile):
        # AeroDyn v14 Airfoil Polar Input File

        # open aerodyn file
        f = open(aerodynFile, 'r')
                
        airfoil = copy.copy(self.fst_vt['AeroDynPolar'])

        # skip through header
        airfoil['description'] = f.readline().rstrip()  # remove newline
        f.readline()
        airfoil['number_tables'] = int(f.readline().split()[0])

        IDParam = [float_read(val) for val in f.readline().split()[0:airfoil['number_tables']]]
        StallAngle = [float_read(val) for val in f.readline().split()[0:airfoil['number_tables']]]
        f.readline()
        f.readline()
        f.readline()
        ZeroCn = [float_read(val) for val in f.readline().split()[0:airfoil['number_tables']]]
        CnSlope = [float_read(val) for val in f.readline().split()[0:airfoil['number_tables']]]
        CnPosStall = [float_read(val) for val in f.readline().split()[0:airfoil['number_tables']]]
        CnNegStall = [float_read(val) for val in f.readline().split()[0:airfoil['number_tables']]]
        alphaCdMin = [float_read(val) for val in f.readline().split()[0:airfoil['number_tables']]]
        CdMin = [float_read(val) for val in f.readline().split()[0:airfoil['number_tables']]]

        data = []
        airfoil['af_tables'] = []
        while True:
            line = f.readline()
            if 'EOT' in line:
                break
            line = [float_read(s) for s in line.split()]
            if len(line) < 1:
                break
            data.append(line)

        # loop through tables
        for i in range(airfoil['number_tables']):
            polar = {}
            polar['IDParam'] = IDParam[i]
            polar['StallAngle'] = StallAngle[i]
            polar['ZeroCn'] = ZeroCn[i]
            polar['CnSlope'] = CnSlope[i]
            polar['CnPosStall'] = CnPosStall[i]
            polar['CnNegStall'] = CnNegStall[i]
            polar['alphaCdMin'] = alphaCdMin[i]
            polar['CdMin'] = CdMin[i]

            alpha = []
            cl = []
            cd = []
            cm = []
            # read polar information line by line
            for datai in data:
                if len(datai) == airfoil['number_tables']*3+1:
                    alpha.append(datai[0])
                    cl.append(datai[1 + 3*i])
                    cd.append(datai[2 + 3*i])
                    cm.append(datai[3 + 3*i])
                elif len(datai) == airfoil['number_tables']*2+1:
                    alpha.append(datai[0])
                    cl.append(datai[1 + 2*i])
                    cd.append(datai[2 + 2*i])

            polar['alpha'] = alpha
            polar['cl'] = cl
            polar['cd'] = cd
            polar['cm'] = cm
            airfoil['af_tables'].append(polar)

        f.close()

        return airfoil

    def read_ServoDyn(self):
        # ServoDyn v1.05 Input File
        # Currently no differences between FASTv8.16 and OpenFAST.


        sd_file = os.path.normpath(os.path.join(self.FAST_directory, self.fst_vt['Fst']['ServoFile']))
        f = open(sd_file)

        f.readline()
        f.readline()

        # Simulation Control (sd_sim_ctrl)
        f.readline()
        self.fst_vt['ServoDyn']['Echo'] = bool_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['DT'] = float_read(f.readline().split()[0])

        # Pitch Control (pitch_ctrl)
        f.readline()
        self.fst_vt['ServoDyn']['PCMode']       = int(f.readline().split()[0])
        self.fst_vt['ServoDyn']['TPCOn']        = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['TPitManS1']    = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['TPitManS2']    = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['TPitManS3']    = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['PitManRat(1)'] = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['PitManRat(2)'] = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['PitManRat(3)'] = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['BlPitchF(1)']  = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['BlPitchF(2)']  = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['BlPitchF(3)']  = float_read(f.readline().split()[0])

        # Geneartor and Torque Control (gen_torq_ctrl)
        f.readline()
        self.fst_vt['ServoDyn']['VSContrl'] = int(f.readline().split()[0])
        self.fst_vt['ServoDyn']['GenModel'] = int(f.readline().split()[0])
        self.fst_vt['ServoDyn']['GenEff']   = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['GenTiStr'] = bool_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['GenTiStp'] = bool_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['SpdGenOn'] = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['TimGenOn'] = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['TimGenOf'] = float_read(f.readline().split()[0])

        # Simple Variable-Speed Torque Control (var_speed_torq_ctrl)
        f.readline()
        self.fst_vt['ServoDyn']['VS_RtGnSp'] = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['VS_RtTq']   = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['VS_Rgn2K']  = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['VS_SlPc']   = float_read(f.readline().split()[0])

        # Simple Induction Generator (induct_gen)
        f.readline()
        self.fst_vt['ServoDyn']['SIG_SlPc'] = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['SIG_SySp'] = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['SIG_RtTq'] = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['SIG_PORt'] = float_read(f.readline().split()[0])

        # Thevenin-Equivalent Induction Generator (theveq_induct_gen)
        f.readline()
        self.fst_vt['ServoDyn']['TEC_Freq'] = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['TEC_NPol'] = int(f.readline().split()[0])
        self.fst_vt['ServoDyn']['TEC_SRes'] = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['TEC_RRes'] = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['TEC_VLL']  = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['TEC_SLR']  = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['TEC_RLR']  = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['TEC_MR']   = float_read(f.readline().split()[0])

        # High-Speed Shaft Brake (shaft_brake)
        f.readline()
        self.fst_vt['ServoDyn']['HSSBrMode'] = int(f.readline().split()[0])
        self.fst_vt['ServoDyn']['THSSBrDp']  = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['HSSBrDT']   = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['HSSBrTqF']  = float_read(f.readline().split()[0])

        # Nacelle-Yaw Control (nac_yaw_ctrl)
        f.readline()
        self.fst_vt['ServoDyn']['YCMode']    = int(f.readline().split()[0])
        self.fst_vt['ServoDyn']['TYCOn']     = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['YawNeut']   = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['YawSpr']    = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['YawDamp']   = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['TYawManS']  = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['YawManRat'] = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['NacYawF']   = float_read(f.readline().split()[0])

        # Aero flow control
        f.readline()
        self.fst_vt['ServoDyn']['AfCmode']      = int(f.readline().split()[0])
        self.fst_vt['ServoDyn']['AfC_Mean']     = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['AfC_Amp']      = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['AfC_Phase']    = float_read(f.readline().split()[0])

        # Structural Control
        f.readline()
        self.fst_vt['ServoDyn']['NumBStC']  = int(f.readline().split()[0])
        self.fst_vt['ServoDyn']['BStCfiles'] = read_array(f,self.fst_vt['ServoDyn']['NumBStC'],str)
        self.fst_vt['ServoDyn']['NumNStC']  = int(f.readline().split()[0])
        self.fst_vt['ServoDyn']['NStCfiles'] = read_array(f,self.fst_vt['ServoDyn']['NumNStC'],str)
        self.fst_vt['ServoDyn']['NumTStC'] = int(f.readline().split()[0])
        self.fst_vt['ServoDyn']['TStCfiles'] = read_array(f,self.fst_vt['ServoDyn']['NumTStC'],str)
        self.fst_vt['ServoDyn']['NumSStC'] = int(f.readline().split()[0])
        self.fst_vt['ServoDyn']['SStCfiles'] = read_array(f,self.fst_vt['ServoDyn']['NumSStC'],str)

        # Initialize Struct Control trees
        self.fst_vt['BStC'] = [] * self.fst_vt['ServoDyn']['NumBStC']
        self.fst_vt['NStC'] = [] * self.fst_vt['ServoDyn']['NumNStC']
        self.fst_vt['TStC'] = [] * self.fst_vt['ServoDyn']['NumTStC']
        self.fst_vt['SStC'] = [] * self.fst_vt['ServoDyn']['NumSStC']

        # Cable control
        f.readline()
        self.fst_vt['ServoDyn']['CCmode']  = int(f.readline().split()[0])

        # Bladed Interface and Torque-Speed Look-Up Table (bladed_interface)
        f.readline()
        if self.path2dll == '' or self.path2dll == None:
            self.fst_vt['ServoDyn']['DLL_FileName'] = os.path.abspath(os.path.normpath(os.path.join(os.path.split(sd_file)[0], f.readline().split()[0][1:-1])))
        else:
            f.readline()
            self.fst_vt['ServoDyn']['DLL_FileName'] = self.path2dll
        self.fst_vt['ServoDyn']['DLL_InFile']   = os.path.abspath(os.path.normpath(os.path.join(os.path.split(sd_file)[0], f.readline().split()[0][1:-1])))
        self.fst_vt['ServoDyn']['DLL_ProcName'] = f.readline().split()[0][1:-1]
        dll_dt_line = f.readline().split()[0]
        try:
            self.fst_vt['ServoDyn']['DLL_DT'] = float_read(dll_dt_line)
        except:
            self.fst_vt['ServoDyn']['DLL_DT'] = dll_dt_line[1:-1]
        self.fst_vt['ServoDyn']['DLL_Ramp']     = bool_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['BPCutoff']     = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['NacYaw_North'] = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['Ptch_Cntrl']   = int(f.readline().split()[0])
        self.fst_vt['ServoDyn']['Ptch_SetPnt']  = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['Ptch_Min']     = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['Ptch_Max']     = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['PtchRate_Min'] = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['PtchRate_Max'] = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['Gain_OM']      = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['GenSpd_MinOM'] = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['GenSpd_MaxOM'] = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['GenSpd_Dem']   = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['GenTrq_Dem']   = float_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['GenPwr_Dem']   = float_read(f.readline().split()[0])

        f.readline()

        self.fst_vt['ServoDyn']['DLL_NumTrq'] = int(f.readline().split()[0])
        f.readline()
        f.readline()
        self.fst_vt['ServoDyn']['GenSpd_TLU'] = [None] * self.fst_vt['ServoDyn']['DLL_NumTrq']
        self.fst_vt['ServoDyn']['GenTrq_TLU'] = [None] * self.fst_vt['ServoDyn']['DLL_NumTrq']
        for i in range(self.fst_vt['ServoDyn']['DLL_NumTrq']):
            data = f.readline().split()
            self.fst_vt['ServoDyn']['GenSpd_TLU'][i]  = float_read(data[0])
            self.fst_vt['ServoDyn']['GenTrq_TLU'][i]  = float_read(data[0])

        # ServoDyn Output Params (sd_out_params)
        f.readline()
        self.fst_vt['ServoDyn']['SumPrint'] = bool_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['OutFile']  = int(f.readline().split()[0])
        self.fst_vt['ServoDyn']['TabDelim'] = bool_read(f.readline().split()[0])
        self.fst_vt['ServoDyn']['OutFmt']   = f.readline().split()[0][1:-1]
        self.fst_vt['ServoDyn']['TStart']   = float_read(f.readline().split()[0])

        # ServoDyn Outlist
        f.readline()
        data = f.readline()
        while data.split()[0] != 'END':
            channels = data.split('"')
            channel_list = channels[1].split(',')
            self.set_outlist(self.fst_vt['outlist']['ServoDyn'], channel_list)
            data = f.readline()

        f.close()

    def read_StC(self,filename):
        '''
        return StC vt so it can be appended to fst_vt['XStC'] list
        '''
        StC_vt = {}

        with open(os.path.join(self.FAST_directory, filename)) as f:

            f.readline()
            f.readline()
            f.readline()
            StC_vt['Echo'] = bool_read(f.readline().split()[0])      #  Echo         - Echo input data to <RootName>.ech (flag)
            f.readline()  # StC DEGREES OF FREEDOM
            StC_vt['StC_DOF_MODE'] = int_read(f.readline().split()[0]) #        4   StC_DOF_MODE - DOF mode (switch) {0: No StC or TLCD DOF; 1: StC_X_DOF, StC_Y_DOF, and/or StC_Z_DOF (three independent StC DOFs); 2: StC_XY_DOF (Omni-Directional StC); 3: TLCD; 4: Prescribed force/moment time series}
            StC_vt['StC_X_DOF'] = bool_read(f.readline().split()[0])  # false         StC_X_DOF    - DOF on or off for StC X (flag) [Used only when StC_DOF_MODE=1]
            StC_vt['StC_Y_DOF'] = bool_read(f.readline().split()[0])  # false         StC_Y_DOF    - DOF on or off for StC Y (flag) [Used only when StC_DOF_MODE=1]
            StC_vt['StC_Z_DOF'] = bool_read(f.readline().split()[0])  # false         StC_Z_DOF    - DOF on or off for StC Z (flag) [Used only when StC_DOF_MODE=1]
            f.readline()    # StC LOCATION 
            StC_vt['StC_P_X'] = float_read(f.readline().split()[0])   #    -51.75   StC_P_X      - At rest X position of StC (m)
            StC_vt['StC_P_Y'] = float_read(f.readline().split()[0])   #        0   StC_P_Y      - At rest Y position of StC (m)
            StC_vt['StC_P_Z'] = float_read(f.readline().split()[0])   #        -10   StC_P_Z      - At rest Z position of StC (m)
            f.readline()    # StC INITIAL CONDITIONS 
            StC_vt['StC_X_DSP'] = float_read(f.readline().split()[0])   #        0   StC_X_DSP    - StC X initial displacement (m) [relative to at rest position]
            StC_vt['StC_Y_DSP'] = float_read(f.readline().split()[0])   #        0   StC_Y_DSP    - StC Y initial displacement (m) [relative to at rest position]
            StC_vt['StC_Z_DSP'] = float_read(f.readline().split()[0])   #        0   StC_Z_DSP    - StC Z initial displacement (m) [relative to at rest position; used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
            StC_vt['StC_Z_PreLd'] = f.readline().split()[0] # "none"        StC_Z_PreLd  - StC Z prefloat_read(f.readline().split()[0])  #-load (N) {"gravity" to offset for gravity load; "none" or 0 to turn off} [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
            f.readline()    # StC CONFIGURATION 
            StC_vt['StC_X_NSP'] = float_read(f.readline().split()[0])   #        0   StC_X_NSP    - Negative stop position (minimum X mass displacement) (m)
            StC_vt['StC_X_PSP'] = float_read(f.readline().split()[0])   #        0   StC_X_PSP    - Positive stop position (maximum X mass displacement) (m)
            StC_vt['StC_Y_PSP'] = float_read(f.readline().split()[0])   #        0   StC_Y_PSP    - Positive stop position (maximum Y mass displacement) (m)
            StC_vt['StC_Y_NSP'] = float_read(f.readline().split()[0])   #        0   StC_Y_NSP    - Negative stop position (minimum Y mass displacement) (m)
            StC_vt['StC_Z_PSP'] = float_read(f.readline().split()[0])   #        0   StC_Z_PSP    - Positive stop position (maximum Z mass displacement) (m) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
            StC_vt['StC_Z_NSP'] = float_read(f.readline().split()[0])   #        0   StC_Z_NSP    - Negative stop position (minimum Z mass displacement) (m) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
            f.readline()    # StC MASS, STIFFNESS, & DAMPING
            StC_vt['StC_X_M'] = float_read(f.readline().split()[0])   #        0   StC_X_M      - StC X mass (kg) [must equal StC_Y_M for StC_DOF_MODE = 2]
            StC_vt['StC_Y_M'] = float_read(f.readline().split()[0])   #        50   StC_Y_M      - StC Y mass (kg) [must equal StC_X_M for StC_DOF_MODE = 2]
            StC_vt['StC_Z_M'] = float_read(f.readline().split()[0])   #        0   StC_Z_M      - StC Z mass (kg) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
            StC_vt['StC_XY_M'] = float_read(f.readline().split()[0])   #        0   StC_XY_M     - StC XY mass (kg) [used only when StC_DOF_MODE=2]
            StC_vt['StC_X_K'] = float_read(f.readline().split()[0])   #    2300   StC_X_K      - StC X stiffness (N/m)
            StC_vt['StC_Y_K'] = float_read(f.readline().split()[0])   #    2300   StC_Y_K      - StC Y stiffness (N/m)
            StC_vt['StC_Z_K'] = float_read(f.readline().split()[0])   #        0   StC_Z_K      - StC Z stiffness (N/m) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
            StC_vt['StC_X_C'] = float_read(f.readline().split()[0])   #        35   StC_X_C      - StC X damping (N/(m/s))
            StC_vt['StC_Y_C'] = float_read(f.readline().split()[0])  #float_read(f.readline().split()[0])  #        35   StC_Y_C      - StC Y damping (N/(m/s))
            StC_vt['StC_Z_C'] = float_read(f.readline().split()[0])  #        0   StC_Z_C      - StC Z damping (N/(m/s)) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
            StC_vt['StC_X_KS'] = float_read(f.readline().split()[0])  #        0   StC_X_KS     - Stop spring X stiffness (N/m)
            StC_vt['StC_Y_KS'] = float_read(f.readline().split()[0])  #       0   StC_Y_KS     - Stop spring Y stiffness (N/m)
            StC_vt['StC_Z_KS'] = float_read(f.readline().split()[0])  #        0   StC_Z_KS     - Stop spring Z stiffness (N/m) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
            StC_vt['StC_X_CS'] = float_read(f.readline().split()[0])  #        0   StC_X_CS     - Stop spring X damping (N/(m/s))
            StC_vt['StC_Y_CS'] = float_read(f.readline().split()[0])  #        0   StC_Y_CS     - Stop spring Y damping (N/(m/s))
            StC_vt['StC_Z_CS'] = float_read(f.readline().split()[0])  #        0   StC_Z_CS     - Stop spring Z damping (N/(m/s)) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
            f.readline()    # StC USER-DEFINED SPRING FORCES
            StC_vt['Use_F_TBL'] = bool_read(f.readline().split()[0])  # False         Use_F_TBL    - Use spring force from user-defined table (flag)
            StC_vt['NKInpSt'] = int_read(f.readline().split()[0]) #        17   NKInpSt      - Number of spring force input stations

            StC_vt['SpringForceTable'] = {}
            table = StC_vt['SpringForceTable']
            table['X']      = [None] * StC_vt['NKInpSt']
            table['F_X']    = [None] * StC_vt['NKInpSt']
            table['Y']      = [None] * StC_vt['NKInpSt']
            table['F_Y']    = [None] * StC_vt['NKInpSt']
            table['Z']      = [None] * StC_vt['NKInpSt']
            table['F_Z']    = [None] * StC_vt['NKInpSt']

            f.readline()
            # if StC_vt['Use_F_TBL']:
            f.readline()
            f.readline()
            for i in range(StC_vt['NKInpSt']):
                ln = f.readline().split()
                table['X'][i]      = float(ln[0])
                table['F_X'][i]    = float(ln[1])
                table['Y'][i]      = float(ln[2])
                table['F_Y'][i]    = float(ln[3])
                table['Z'][i]      = float(ln[4])
                table['F_Z'][i]    = float(ln[5])
            # else:
            #     # Skip until next section
            #     data_line = f.readline().strip().split()
            #     while data_line[0][:3] != '---': 
            #         data_line = f.readline().strip().split()

            # StructCtrl CONTROL, skip this readline() because it already happened
            f.readline()
            StC_vt['StC_CMODE'] = int_read(f.readline().split()[0]) #        5   StC_CMODE     - Control mode (switch) {0:none; 1: Semi-Active Control Mode; 4: Active Control Mode through Simulink (not available); 5: Active Control Mode through Bladed interface}
            StC_vt['StC_CChan'] = int_read(f.readline().split()[0]) #        0   StC_CChan     - Control channel group (1:10) for stiffness and damping (StC_[XYZ]_K, StC_[XYZ]_C, and StC_[XYZ]_Brake) (specify additional channels for blade instances of StC active control -- one channel per blade) [used only when StC_DOF_MODE=1 or 2, and StC_CMODE=4 or 5]
            StC_vt['StC_SA_MODE'] = int_read(f.readline().split()[0]) #        1   StC_SA_MODE   - Semi-Active control mode {1: velocity-based ground hook control; 2: Inverse velocity-based ground hook control; 3: displacement-based ground hook control 4: Phase difference Algorithm with Friction Force 5: Phase difference Algorithm with Damping Force} (-)
            StC_vt['StC_X_C_LOW'] = float_read(f.readline().split()[0])  #        0   StC_X_C_LOW   - StC X low damping for ground hook control
            StC_vt['StC_X_C_HIGH'] = float_read(f.readline().split()[0])  #        0   StC_X_C_HIGH  - StC X high damping for ground hook control
            StC_vt['StC_Y_C_HIGH'] = float_read(f.readline().split()[0])  #        0   StC_Y_C_HIGH  - StC Y high damping for ground hook control
            StC_vt['StC_Y_C_LOW'] = float_read(f.readline().split()[0])  #        0   StC_Y_C_LOW   - StC Y low damping for ground hook control
            StC_vt['StC_Z_C_HIGH'] = float_read(f.readline().split()[0])  #        0   StC_Z_C_HIGH  - StC Z high damping for ground hook control [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
            StC_vt['StC_Z_C_LOW'] = float_read(f.readline().split()[0])  #        0   StC_Z_C_LOW   - StC Z low damping for ground hook control  [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
            StC_vt['StC_X_C_BRAKE'] = float_read(f.readline().split()[0])  #        0   StC_X_C_BRAKE - StC X high damping for braking the StC (Don't use it now. should be zero)
            StC_vt['StC_Y_C_BRAKE'] = float_read(f.readline().split()[0])  #        0   StC_Y_C_BRAKE - StC Y high damping for braking the StC (Don't use it now. should be zero)
            StC_vt['StC_Z_C_BRAKE'] = float_read(f.readline().split()[0])  #        0   StC_Z_C_BRAKE - StC Z high damping for braking the StC (Don't use it now. should be zero) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
            f.readline()    # TLCD
            StC_vt['L_X'] = float_read(f.readline().split()[0])  #   7.9325     L_X             - X TLCD total length (m)
            StC_vt['B_X'] = float_read(f.readline().split()[0])  #   6.5929     B_X             - X TLCD horizontal length (m)
            StC_vt['area_X'] = float_read(f.readline().split()[0])  #   2.0217     area_X          - X TLCD cross-sectional area of vertical column (m^2)
            StC_vt['area_ratio_X'] = float_read(f.readline().split()[0])  #   0.913      area_ratio_X    - X TLCD cross-sectional area ratio (vertical column area divided by horizontal column area) (-)
            StC_vt['headLossCoeff_X'] = float_read(f.readline().split()[0])  #   2.5265     headLossCoeff_X - X TLCD head loss coeff (-)
            StC_vt['rho_X'] = float_read(f.readline().split()[0])  #       1000     rho_X           - X TLCD liquid density (kg/m^3)
            StC_vt['L_Y'] = float_read(f.readline().split()[0])  #   3.5767     L_Y             - Y TLCD total length (m)
            StC_vt['B_Y'] = float_read(f.readline().split()[0])  #   2.1788     B_Y             - Y TLCD horizontal length (m)
            StC_vt['area_Y'] = float_read(f.readline().split()[0])  #   1.2252     area_Y          - Y TLCD cross-sectional area of vertical column (m^2)
            StC_vt['area_ratio_Y'] = float_read(f.readline().split()[0])  #   2.7232     area_ratio_Y    - Y TLCD cross-sectional area ratio (vertical column area divided by horizontal column area) (-)
            StC_vt['headLossCoeff_Y'] = float_read(f.readline().split()[0])  #   0.6433     headLossCoeff_Y - Y TLCD head loss coeff (-)
            StC_vt['rho_Y'] = float_read(f.readline().split()[0])  #       1000     rho_Y           - Y TLCD liquid density (kg/m^3)
            f.readline()    # PRESCRIBED TIME SERIES 
            StC_vt['PrescribedForcesCoord'] = int_read(f.readline().split()[0]) #        2   PrescribedForcesCoord- Prescribed forces are in global or local coordinates (switch) {1: global; 2: local}
            StC_vt['PrescribedForcesFile'] = f.readline().split()[0] # "Bld-TimeForceSeries.dat"  PrescribedForcesFile   - Time series force and moment (7 columns of time, FX, FY, FZ, MX, MY, MZ)
            f.readline()

        return StC_vt
    
    def read_DISCON_in(self):
        # Read the Bladed style Interface controller input file, intended for ROSCO https://github.com/NREL/rosco.toolbox

        discon_in_file = os.path.normpath(os.path.join(self.FAST_directory, self.fst_vt['ServoDyn']['DLL_InFile']))

        if os.path.exists(discon_in_file):

            # Read DISCON infiles
            self.fst_vt['DISCON_in'] = read_DISCON(discon_in_file)

            # Some additional filename parsing
            self.fst_vt['DISCON_in']['PerfFileName'] = os.path.abspath(os.path.join(self.FAST_directory, self.fst_vt['DISCON_in']['PerfFileName']))

            # Try to read rotor performance data if it is available
            try:
                pitch_vector, tsr_vector, Cp_table, Ct_table, Cq_table = load_from_txt(self.fst_vt['DISCON_in']['PerfFileName'])

                RotorPerformance = ROSCO_turbine.RotorPerformance
                Cp = RotorPerformance(Cp_table, pitch_vector, tsr_vector)
                Ct = RotorPerformance(Ct_table, pitch_vector, tsr_vector)
                Cq = RotorPerformance(Cq_table, pitch_vector, tsr_vector)

                self.fst_vt['DISCON_in']['Cp'] = Cp
                self.fst_vt['DISCON_in']['Ct'] = Ct
                self.fst_vt['DISCON_in']['Cq'] = Cq
                self.fst_vt['DISCON_in']['Cp_pitch_initial_rad'] = pitch_vector
                self.fst_vt['DISCON_in']['Cp_TSR_initial'] = tsr_vector
                self.fst_vt['DISCON_in']['Cp_table'] = Cp_table
                self.fst_vt['DISCON_in']['Ct_table'] = Ct_table
                self.fst_vt['DISCON_in']['Cq_table'] = Cq_table
            except:
                print('WARNING: Cp table not loaded!')
            
            # Add some DISCON entries that might be needed within WISDEM        
            self.fst_vt['DISCON_in']['v_rated'] = 1.
        
        else:
            del self.fst_vt['DISCON_in']

    def read_HydroDyn(self, hd_file):

        f = open(hd_file)

        f.readline()
        f.readline()

        self.fst_vt['HydroDyn']['Echo'] = bool_read(f.readline().split()[0])
        # ENVIRONMENTAL CONDITIONS
        f.readline()
        self.fst_vt['HydroDyn']['WtrDens'] = float_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['WtrDpth'] = float_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['MSL2SWL'] = float_read(f.readline().split()[0])

        # WAVES
        f.readline()
        self.fst_vt['HydroDyn']['WaveMod']       = int_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['WaveStMod']     = int_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['WaveTMax']      = float_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['WaveDT']        = float_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['WaveHs']        = float_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['WaveTp']        = float_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['WavePkShp']     = float_read(f.readline().split()[0]) # default
        self.fst_vt['HydroDyn']['WvLowCOff']     = float_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['WvHiCOff']      = float_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['WaveDir']       = float_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['WaveDirMod']    = int_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['WaveDirSpread'] = float_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['WaveNDir']      = int_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['WaveDirRange']  = float_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['WaveSeed1']     = int_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['WaveSeed2']     = int_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['WaveNDAmp']     = bool_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['WvKinFile']     = f.readline().split()[0][1:-1]
        self.fst_vt['HydroDyn']['NWaveElev']     = int_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['WaveElevxi']    = [idx.strip() for idx in f.readline().split('WaveElevxi')[0].split(',')]
        self.fst_vt['HydroDyn']['WaveElevyi']    = [idx.strip() for idx in f.readline().split('WaveElevyi')[0].split(',')]

        # 2ND-ORDER WAVES
        f.readline()
        self.fst_vt['HydroDyn']['WvDiffQTF']     = bool_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['WvSumQTF']      = bool_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['WvLowCOffD']    = float_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['WvHiCOffD']     = float_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['WvLowCOffS']    = float_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['WvHiCOffS']     = float_read(f.readline().split()[0])

        # CURRENT
        f.readline()
        self.fst_vt['HydroDyn']['CurrMod']       = int_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['CurrSSV0']      = float_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['CurrSSDir']     = float_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['CurrNSRef']     = float_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['CurrNSV0']      = float_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['CurrNSDir']     = float_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['CurrDIV']       = float_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['CurrDIDir']     = float_read(f.readline().split()[0])

        # FLOATING PLATFORM
        f.readline()
        self.fst_vt['HydroDyn']['PotMod']        = int_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['ExctnMod']      = int_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['RdtnMod']       = int_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['RdtnTMax']      = float_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['RdtnDT']        = float_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['NBody']         = int_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['NBodyMod']      = int_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['PotFile']       = os.path.normpath(os.path.join(os.path.split(hd_file)[0], f.readline().split()[0][1:-1]))
        self.fst_vt['HydroDyn']['WAMITULEN']     = float_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['PtfmRefxt']     = float_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['PtfmRefyt']     = float_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['PtfmRefzt']     = float_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['PtfmRefztRot'] = float_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['PtfmVol0']      = float_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['PtfmCOBxt']     = float_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['PtfmCOByt']     = float_read(f.readline().split()[0])

        # 2ND-ORDER FLOATING PLATFORM FORCES
        f.readline()
        self.fst_vt['HydroDyn']['MnDrift']       = int_read(f.readline().split()[0]) # ?
        self.fst_vt['HydroDyn']['NewmanApp']     = int_read(f.readline().split()[0]) # ?
        self.fst_vt['HydroDyn']['DiffQTF']       = int_read(f.readline().split()[0]) # ?
        self.fst_vt['HydroDyn']['SumQTF']        = int_read(f.readline().split()[0]) # ?

        # PLATFORM ADDITIONAL STIFFNESS AND DAMPING
        f.readline()
        # Get number of F0 terms [If NBodyMod=1, one size 6*NBody x 1 vector; if NBodyMod>1, NBody size 6 x 1 vectors]
        NBody = self.fst_vt['HydroDyn']['NBody']
        if self.fst_vt['HydroDyn']['NBodyMod'] == 1:
            self.fst_vt['HydroDyn']['AddF0']         = [float(f.readline().strip().split()[0]) for i in range(6*NBody)]
        elif self.fst_vt['HydroDyn']['NBodyMod'] > 1:
            self.fst_vt['HydroDyn']['AddF0']         = [[float(idx) for idx in f.readline().strip().split()[:NBody]] for i in range(6)]
        else:
            raise Exception("Invalid value for fst_vt['HydroDyn']['NBodyMod']")

        self.fst_vt['HydroDyn']['AddCLin']       = np.array([[float(idx) for idx in f.readline().strip().split()[:6]] for i in range(6)])
        self.fst_vt['HydroDyn']['AddBLin']       = np.array([[float(idx) for idx in f.readline().strip().split()[:6]] for i in range(6)])
        self.fst_vt['HydroDyn']['AddBQuad']      = np.array([[float(idx) for idx in f.readline().strip().split()[:6]] for i in range(6)])

        #AXIAL COEFFICIENTS
        f.readline()
        self.fst_vt['HydroDyn']['NAxCoef']       = int_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['AxCoefID']      = [None]*self.fst_vt['HydroDyn']['NAxCoef']
        self.fst_vt['HydroDyn']['AxCd']          = [None]*self.fst_vt['HydroDyn']['NAxCoef']
        self.fst_vt['HydroDyn']['AxCa']          = [None]*self.fst_vt['HydroDyn']['NAxCoef']
        self.fst_vt['HydroDyn']['AxCp']          = [None]*self.fst_vt['HydroDyn']['NAxCoef']
        ln = f.readline().split()
        ln = f.readline().split()
        for i in range(self.fst_vt['HydroDyn']['NAxCoef']):
            ln = f.readline().split()
            self.fst_vt['HydroDyn']['AxCoefID'][i] = int(ln[0])
            self.fst_vt['HydroDyn']['AxCd'][i]     = float(ln[1])
            self.fst_vt['HydroDyn']['AxCa'][i]     = float(ln[2])
            self.fst_vt['HydroDyn']['AxCp'][i]     = float(ln[3])

        #MEMBER JOINTS
        f.readline()
        self.fst_vt['HydroDyn']['NJoints']    = int_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['JointID']    = [None]*self.fst_vt['HydroDyn']['NJoints']
        self.fst_vt['HydroDyn']['Jointxi']    = [None]*self.fst_vt['HydroDyn']['NJoints']
        self.fst_vt['HydroDyn']['Jointyi']    = [None]*self.fst_vt['HydroDyn']['NJoints']
        self.fst_vt['HydroDyn']['Jointzi']    = [None]*self.fst_vt['HydroDyn']['NJoints']
        self.fst_vt['HydroDyn']['JointAxID']  = [None]*self.fst_vt['HydroDyn']['NJoints']
        self.fst_vt['HydroDyn']['JointOvrlp'] = [None]*self.fst_vt['HydroDyn']['NJoints']
        ln = f.readline().split()
        ln = f.readline().split()
        for i in range(self.fst_vt['HydroDyn']['NJoints']):
            ln = f.readline().split()
            self.fst_vt['HydroDyn']['JointID'][i]    = int(ln[0])
            self.fst_vt['HydroDyn']['Jointxi'][i]    = float(ln[1])
            self.fst_vt['HydroDyn']['Jointyi'][i]    = float(ln[2])
            self.fst_vt['HydroDyn']['Jointzi'][i]    = float(ln[3])
            self.fst_vt['HydroDyn']['JointAxID'][i]  = int(ln[4])
            self.fst_vt['HydroDyn']['JointOvrlp'][i] = int(ln[5])

        #MEMBER CROSS-SECTION PROPERTIES
        f.readline()
        self.fst_vt['HydroDyn']['NPropSets'] = int_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['PropSetID'] = [None]*self.fst_vt['HydroDyn']['NPropSets']
        self.fst_vt['HydroDyn']['PropD']     = [None]*self.fst_vt['HydroDyn']['NPropSets']
        self.fst_vt['HydroDyn']['PropThck']  = [None]*self.fst_vt['HydroDyn']['NPropSets']
        ln = f.readline().split()
        ln = f.readline().split()
        for i in range(self.fst_vt['HydroDyn']['NPropSets']):
            ln = f.readline().split()
            self.fst_vt['HydroDyn']['PropSetID'][i] = int(ln[0])
            self.fst_vt['HydroDyn']['PropD'][i]     = float(ln[1])
            self.fst_vt['HydroDyn']['PropThck'][i]  = float(ln[2])

        #SIMPLE HYDRODYNAMIC COEFFICIENTS
        f.readline()
        f.readline()
        f.readline()
        ln = f.readline().split()
        self.fst_vt['HydroDyn']['SimplCd']     = float(ln[0])
        self.fst_vt['HydroDyn']['SimplCdMG']   = float(ln[1])
        self.fst_vt['HydroDyn']['SimplCa']     = float(ln[2])
        self.fst_vt['HydroDyn']['SimplCaMG']   = float(ln[3])
        self.fst_vt['HydroDyn']['SimplCp']     = float(ln[4])
        self.fst_vt['HydroDyn']['SimplCpMG']   = float(ln[5])
        self.fst_vt['HydroDyn']['SimplAxCd']   = float(ln[6])
        self.fst_vt['HydroDyn']['SimplAxCdMG'] = float(ln[7])
        self.fst_vt['HydroDyn']['SimplAxCa'] = float(ln[8])
        self.fst_vt['HydroDyn']['SimplAxCaMG'] = float(ln[9])
        self.fst_vt['HydroDyn']['SimplAxCp'] = float(ln[10])
        self.fst_vt['HydroDyn']['SimplAxCpMG'] = float(ln[11])

        #DEPTH-BASED HYDRODYNAMIC COEFFICIENTS
        f.readline()
        self.fst_vt['HydroDyn']['NCoefDpth']  = int_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['Dpth']       = [None]*self.fst_vt['HydroDyn']['NCoefDpth']
        self.fst_vt['HydroDyn']['DpthCd']     = [None]*self.fst_vt['HydroDyn']['NCoefDpth']
        self.fst_vt['HydroDyn']['DpthCdMG']   = [None]*self.fst_vt['HydroDyn']['NCoefDpth']
        self.fst_vt['HydroDyn']['DpthCa']     = [None]*self.fst_vt['HydroDyn']['NCoefDpth']
        self.fst_vt['HydroDyn']['DpthCaMG']   = [None]*self.fst_vt['HydroDyn']['NCoefDpth']
        self.fst_vt['HydroDyn']['DpthCp']     = [None]*self.fst_vt['HydroDyn']['NCoefDpth']
        self.fst_vt['HydroDyn']['DpthCpMG']   = [None]*self.fst_vt['HydroDyn']['NCoefDpth']
        self.fst_vt['HydroDyn']['DpthAxCd']   = [None]*self.fst_vt['HydroDyn']['NCoefDpth']
        self.fst_vt['HydroDyn']['DpthAxCdMG'] = [None]*self.fst_vt['HydroDyn']['NCoefDpth']
        self.fst_vt['HydroDyn']['DpthAxCa']   = [None]*self.fst_vt['HydroDyn']['NCoefDpth']
        self.fst_vt['HydroDyn']['DpthAxCaMG'] = [None]*self.fst_vt['HydroDyn']['NCoefDpth']
        self.fst_vt['HydroDyn']['DpthAxCp']   = [None]*self.fst_vt['HydroDyn']['NCoefDpth']
        self.fst_vt['HydroDyn']['DpthAxCpMG'] = [None]*self.fst_vt['HydroDyn']['NCoefDpth']
        ln = f.readline().split()
        ln = f.readline().split()
        for i in range(self.fst_vt['HydroDyn']['NCoefDpth']):
            ln = f.readline().split()
            self.fst_vt['HydroDyn']['Dpth'][i]       = float(ln[0])
            self.fst_vt['HydroDyn']['DpthCd'][i]     = float(ln[1])
            self.fst_vt['HydroDyn']['DpthCdMG'][i]   = float(ln[2])
            self.fst_vt['HydroDyn']['DpthCa'][i]     = float(ln[3])
            self.fst_vt['HydroDyn']['DpthCaMG'][i]   = float(ln[4])
            self.fst_vt['HydroDyn']['DpthCp'][i]     = float(ln[5])
            self.fst_vt['HydroDyn']['DpthCpMG'][i]   = float(ln[6])
            self.fst_vt['HydroDyn']['DpthAxCd'][i]   = float(ln[7])
            self.fst_vt['HydroDyn']['DpthAxCdMG'][i] = float(ln[8])
            self.fst_vt['HydroDyn']['DpthAxCa'][i]   = float(ln[9])
            self.fst_vt['HydroDyn']['DpthAxCaMG'][i] = float(ln[10])
            self.fst_vt['HydroDyn']['DpthAxCp'][i]   = float(ln[11])
            self.fst_vt['HydroDyn']['DpthAxCpMG'][i] = float(ln[12])

        #MEMBER-BASED HYDRODYNAMIC COEFFICIENTS
        f.readline()
        self.fst_vt['HydroDyn']['NCoefMembers']  = int_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['MemberID_HydC']      = [None]*self.fst_vt['HydroDyn']['NCoefMembers']
        self.fst_vt['HydroDyn']['MemberCd1']     = [None]*self.fst_vt['HydroDyn']['NCoefMembers']
        self.fst_vt['HydroDyn']['MemberCd2']     = [None]*self.fst_vt['HydroDyn']['NCoefMembers']
        self.fst_vt['HydroDyn']['MemberCdMG1']   = [None]*self.fst_vt['HydroDyn']['NCoefMembers']
        self.fst_vt['HydroDyn']['MemberCdMG2']   = [None]*self.fst_vt['HydroDyn']['NCoefMembers']
        self.fst_vt['HydroDyn']['MemberCa1']     = [None]*self.fst_vt['HydroDyn']['NCoefMembers']
        self.fst_vt['HydroDyn']['MemberCa2']     = [None]*self.fst_vt['HydroDyn']['NCoefMembers']
        self.fst_vt['HydroDyn']['MemberCaMG1']   = [None]*self.fst_vt['HydroDyn']['NCoefMembers']
        self.fst_vt['HydroDyn']['MemberCaMG2']   = [None]*self.fst_vt['HydroDyn']['NCoefMembers']
        self.fst_vt['HydroDyn']['MemberCp1']     = [None]*self.fst_vt['HydroDyn']['NCoefMembers']
        self.fst_vt['HydroDyn']['MemberCp2']     = [None]*self.fst_vt['HydroDyn']['NCoefMembers']
        self.fst_vt['HydroDyn']['MemberCpMG1']   = [None]*self.fst_vt['HydroDyn']['NCoefMembers']
        self.fst_vt['HydroDyn']['MemberCpMG2']   = [None]*self.fst_vt['HydroDyn']['NCoefMembers']
        self.fst_vt['HydroDyn']['MemberAxCd1']   = [None]*self.fst_vt['HydroDyn']['NCoefMembers']
        self.fst_vt['HydroDyn']['MemberAxCd2']   = [None]*self.fst_vt['HydroDyn']['NCoefMembers']
        self.fst_vt['HydroDyn']['MemberAxCdMG1'] = [None]*self.fst_vt['HydroDyn']['NCoefMembers']
        self.fst_vt['HydroDyn']['MemberAxCdMG2'] = [None]*self.fst_vt['HydroDyn']['NCoefMembers']
        self.fst_vt['HydroDyn']['MemberAxCa1']   = [None]*self.fst_vt['HydroDyn']['NCoefMembers']
        self.fst_vt['HydroDyn']['MemberAxCa2']   = [None]*self.fst_vt['HydroDyn']['NCoefMembers']
        self.fst_vt['HydroDyn']['MemberAxCaMG1'] = [None]*self.fst_vt['HydroDyn']['NCoefMembers']
        self.fst_vt['HydroDyn']['MemberAxCaMG2'] = [None]*self.fst_vt['HydroDyn']['NCoefMembers']
        self.fst_vt['HydroDyn']['MemberAxCp1']   = [None]*self.fst_vt['HydroDyn']['NCoefMembers']
        self.fst_vt['HydroDyn']['MemberAxCp2']   = [None]*self.fst_vt['HydroDyn']['NCoefMembers']
        self.fst_vt['HydroDyn']['MemberAxCpMG1'] = [None]*self.fst_vt['HydroDyn']['NCoefMembers']
        self.fst_vt['HydroDyn']['MemberAxCpMG2'] = [None]*self.fst_vt['HydroDyn']['NCoefMembers']

        f.readline()
        f.readline()
        for i in range(self.fst_vt['HydroDyn']['NCoefMembers']):
            ln = f.readline().split()
            self.fst_vt['HydroDyn']['MemberID_HydC'][i]      = int(ln[0])
            self.fst_vt['HydroDyn']['MemberCd1'][i]     = float(ln[1])
            self.fst_vt['HydroDyn']['MemberCd2'][i]     = float(ln[2])
            self.fst_vt['HydroDyn']['MemberCdMG1'][i]   = float(ln[3])
            self.fst_vt['HydroDyn']['MemberCdMG2'][i]   = float(ln[4])
            self.fst_vt['HydroDyn']['MemberCa1'][i]     = float(ln[5])
            self.fst_vt['HydroDyn']['MemberCa2'][i]     = float(ln[6])
            self.fst_vt['HydroDyn']['MemberCaMG1'][i]   = float(ln[7])
            self.fst_vt['HydroDyn']['MemberCaMG2'][i]   = float(ln[8])
            self.fst_vt['HydroDyn']['MemberCp1'][i]     = float(ln[9])
            self.fst_vt['HydroDyn']['MemberCp2'][i]     = float(ln[10])
            self.fst_vt['HydroDyn']['MemberCpMG1'][i]   = float(ln[11])
            self.fst_vt['HydroDyn']['MemberCpMG2'][i]   = float(ln[12])
            self.fst_vt['HydroDyn']['MemberAxCd1'][i]   = float(ln[13])
            self.fst_vt['HydroDyn']['MemberAxCd2'][i]   = float(ln[14])
            self.fst_vt['HydroDyn']['MemberAxCdMG1'][i] = float(ln[15])
            self.fst_vt['HydroDyn']['MemberAxCdMG2'][i] = float(ln[16])
            self.fst_vt['HydroDyn']['MemberAxCa1'][i]   = float(ln[17])
            self.fst_vt['HydroDyn']['MemberAxCa2'][i]   = float(ln[18])
            self.fst_vt['HydroDyn']['MemberAxCaMG1'][i] = float(ln[19])
            self.fst_vt['HydroDyn']['MemberAxCaMG2'][i] = float(ln[20])
            self.fst_vt['HydroDyn']['MemberAxCp1'][i]   = float(ln[21])
            self.fst_vt['HydroDyn']['MemberAxCp2'][i]   = float(ln[22])
            self.fst_vt['HydroDyn']['MemberAxCpMG1'][i] = float(ln[23])
            self.fst_vt['HydroDyn']['MemberAxCpMG2'][i] = float(ln[24])

        #MEMBERS
        f.readline()
        self.fst_vt['HydroDyn']['NMembers']    = int_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['MemberID']    = [None]*self.fst_vt['HydroDyn']['NMembers']
        self.fst_vt['HydroDyn']['MJointID1']   = [None]*self.fst_vt['HydroDyn']['NMembers']
        self.fst_vt['HydroDyn']['MJointID2']   = [None]*self.fst_vt['HydroDyn']['NMembers']
        self.fst_vt['HydroDyn']['MPropSetID1'] = [None]*self.fst_vt['HydroDyn']['NMembers']
        self.fst_vt['HydroDyn']['MPropSetID2'] = [None]*self.fst_vt['HydroDyn']['NMembers']
        self.fst_vt['HydroDyn']['MDivSize']    = [None]*self.fst_vt['HydroDyn']['NMembers']
        self.fst_vt['HydroDyn']['MCoefMod']    = [None]*self.fst_vt['HydroDyn']['NMembers']
        self.fst_vt['HydroDyn']['PropPot']     = [None]*self.fst_vt['HydroDyn']['NMembers']
        ln = f.readline().split()
        ln = f.readline().split()
        for i in range(self.fst_vt['HydroDyn']['NMembers']):
            ln = f.readline().split()
            self.fst_vt['HydroDyn']['MemberID'][i]    = int(ln[0])
            self.fst_vt['HydroDyn']['MJointID1'][i]   = int(ln[1])
            self.fst_vt['HydroDyn']['MJointID2'][i]   = int(ln[2])
            self.fst_vt['HydroDyn']['MPropSetID1'][i] = int(ln[3])
            self.fst_vt['HydroDyn']['MPropSetID2'][i] = int(ln[4])
            self.fst_vt['HydroDyn']['MDivSize'][i]    = float(ln[5])
            self.fst_vt['HydroDyn']['MCoefMod'][i]    = int(ln[6])
            self.fst_vt['HydroDyn']['PropPot'][i]     = bool_read(ln[7])

        #FILLED MEMBERS
        f.readline()
        self.fst_vt['HydroDyn']['NFillGroups'] = int_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['FillNumM']    = [None]*self.fst_vt['HydroDyn']['NFillGroups']
        self.fst_vt['HydroDyn']['FillMList']   = [None]*self.fst_vt['HydroDyn']['NFillGroups']
        self.fst_vt['HydroDyn']['FillFSLoc']   = [None]*self.fst_vt['HydroDyn']['NFillGroups']
        self.fst_vt['HydroDyn']['FillDens']    = [None]*self.fst_vt['HydroDyn']['NFillGroups']
        ln = f.readline().split()
        ln = f.readline().split()
        for i in range(self.fst_vt['HydroDyn']['NFillGroups']):
            ln = f.readline().split()
            self.fst_vt['HydroDyn']['FillNumM'][i]  = int(ln[0])
            self.fst_vt['HydroDyn']['FillMList'][i] = [int(j) for j in ln[1:-2]]
            self.fst_vt['HydroDyn']['FillFSLoc'][i] = float(ln[-2])
            self.fst_vt['HydroDyn']['FillDens'][i]  = float(ln[-1])

        #MARINE GROWTH
        f.readline()
        self.fst_vt['HydroDyn']['NMGDepths'] = int_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['MGDpth']    = [None]*self.fst_vt['HydroDyn']['NMGDepths']
        self.fst_vt['HydroDyn']['MGThck']    = [None]*self.fst_vt['HydroDyn']['NMGDepths']
        self.fst_vt['HydroDyn']['MGDens']    = [None]*self.fst_vt['HydroDyn']['NMGDepths']
        ln = f.readline().split()
        ln = f.readline().split()
        for i in range(self.fst_vt['HydroDyn']['NMGDepths']):
            ln = f.readline().split()
            self.fst_vt['HydroDyn']['MGDpth'][i] = float(ln[0])
            self.fst_vt['HydroDyn']['MGThck'][i] = float(ln[1])
            self.fst_vt['HydroDyn']['MGDens'][i] = float(ln[2])

        #MEMBER OUTPUT LIST
        f.readline()
        self.fst_vt['HydroDyn']['NMOutputs'] = int_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['MemberID_out']  = [None]*self.fst_vt['HydroDyn']['NMOutputs']
        self.fst_vt['HydroDyn']['NOutLoc']   = [None]*self.fst_vt['HydroDyn']['NMOutputs']
        self.fst_vt['HydroDyn']['NodeLocs']  = [None]*self.fst_vt['HydroDyn']['NMOutputs']
        ln = f.readline().split()
        ln = f.readline().split()
        for i in range(self.fst_vt['HydroDyn']['NMOutputs']):
            ln = f.readline().split()
            self.fst_vt['HydroDyn']['MemberID_out'][i] = int(ln[0])
            self.fst_vt['HydroDyn']['NOutLoc'][i]  = int(ln[1])
            self.fst_vt['HydroDyn']['NodeLocs'][i] = float(ln[2])

        #JOINT OUTPUT LIST
        f.readline()
        self.fst_vt['HydroDyn']['NJOutputs'] = int_read(f.readline().split()[0])
        if int(self.fst_vt['HydroDyn']['NJOutputs']) > 0:
            self.fst_vt['HydroDyn']['JOutLst']   = [int(idx.strip()) for idx in f.readline().split('JOutLst')[0].split(',')]
        else:
            f.readline()
            self.fst_vt['HydroDyn']['JOutLst']   = [0]

        #OUTPUT
        f.readline()
        self.fst_vt['HydroDyn']['HDSum']     = bool_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['OutAll']    = bool_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['OutSwtch']  = int_read(f.readline().split()[0])
        self.fst_vt['HydroDyn']['OutFmt']    = str(f.readline().split()[0])
        self.fst_vt['HydroDyn']['OutSFmt']   = str(f.readline().split()[0])

        self.fst_vt['HydroDyn']['HDSum']   
        self.fst_vt['HydroDyn']['OutAll']  
        self.fst_vt['HydroDyn']['OutSwtch']
        self.fst_vt['HydroDyn']['OutFmt']  
        self.fst_vt['HydroDyn']['OutSFmt'] 

        # HydroDyn Outlist
        f.readline()
        data = f.readline()
        while data.split()[0] != 'END':
            channels = data.split('"')
            channel_list = channels[1].split(',')
            self.set_outlist(self.fst_vt['outlist']['HydroDyn'], channel_list)
            data = f.readline()

        f.close()

    def read_SubDyn(self, sd_file):

        f = open(sd_file)
        f.readline()
        f.readline()
        f.readline()
        # SIMULATION CONTROL
        self.fst_vt['SubDyn']['Echo']      = bool_read(f.readline().split()[0])
        self.fst_vt['SubDyn']['SDdeltaT']  = float_read(f.readline().split()[0])
        self.fst_vt['SubDyn']['IntMethod'] = int_read(f.readline().split()[0])
        self.fst_vt['SubDyn']['SttcSolve'] = bool_read(f.readline().split()[0])
        self.fst_vt['SubDyn']['GuyanLoadCorrection'] = bool_read(f.readline().split()[0])
        f.readline()
        # FEA and CRAIG-BAMPTON PARAMETERS
        self.fst_vt['SubDyn']['FEMMod']    = int_read(f.readline().split()[0])
        self.fst_vt['SubDyn']['NDiv']      = int_read(f.readline().split()[0])
        self.fst_vt['SubDyn']['CBMod']     = bool_read(f.readline().split()[0])
        self.fst_vt['SubDyn']['Nmodes']    = int_read(f.readline().split()[0])
        self.fst_vt['SubDyn']['JDampings'] = float_read(f.readline().split()[0])
        self.fst_vt['SubDyn']['GuyanDampMod'] = int_read(f.readline().split()[0])
        self.fst_vt['SubDyn']['RayleighDamp'] = read_array(f,2,array_type=float)
        self.fst_vt['SubDyn']['GuyanDampSize'] = int_read(f.readline().split()[0])
        self.fst_vt['SubDyn']['GuyanDamp'] = np.array([[float(idx) for idx in f.readline().strip().split()[:6]] for i in range(self.fst_vt['SubDyn']['GuyanDampSize'])])
        f.readline()
        # STRUCTURE JOINTS
        self.fst_vt['SubDyn']['NJoints']   = int_read(f.readline().split()[0])
        self.fst_vt['SubDyn']['JointID']   = [None]*self.fst_vt['SubDyn']['NJoints']
        self.fst_vt['SubDyn']['JointXss']  = [None]*self.fst_vt['SubDyn']['NJoints']
        self.fst_vt['SubDyn']['JointYss']  = [None]*self.fst_vt['SubDyn']['NJoints']
        self.fst_vt['SubDyn']['JointZss']  = [None]*self.fst_vt['SubDyn']['NJoints']
        self.fst_vt['SubDyn']['JointType']  = [None]*self.fst_vt['SubDyn']['NJoints']
        self.fst_vt['SubDyn']['JointDirX']  = [None]*self.fst_vt['SubDyn']['NJoints']
        self.fst_vt['SubDyn']['JointDirY']  = [None]*self.fst_vt['SubDyn']['NJoints']
        self.fst_vt['SubDyn']['JointDirZ']  = [None]*self.fst_vt['SubDyn']['NJoints']
        self.fst_vt['SubDyn']['JointStiff']  = [None]*self.fst_vt['SubDyn']['NJoints']
        ln = f.readline().split()
        ln = f.readline().split()
        for i in range(self.fst_vt['SubDyn']['NJoints']):
            ln = f.readline().split()
            self.fst_vt['SubDyn']['JointID'][i]    = int(ln[0])
            self.fst_vt['SubDyn']['JointXss'][i]   = float(ln[1])
            self.fst_vt['SubDyn']['JointYss'][i]   = float(ln[2])
            self.fst_vt['SubDyn']['JointZss'][i]   = float(ln[3])
            self.fst_vt['SubDyn']['JointType'][i]   = int(ln[4])
            self.fst_vt['SubDyn']['JointDirX'][i]   = float(ln[5])
            self.fst_vt['SubDyn']['JointDirY'][i]   = float(ln[6])
            self.fst_vt['SubDyn']['JointDirZ'][i]   = float(ln[7])
            self.fst_vt['SubDyn']['JointStiff'][i]   = float(ln[8])
        f.readline()
        # BASE REACTION JOINTS
        self.fst_vt['SubDyn']['NReact']   = int_read(f.readline().split()[0])
        self.fst_vt['SubDyn']['RJointID'] = [None]*self.fst_vt['SubDyn']['NReact']
        self.fst_vt['SubDyn']['RctTDXss'] = [None]*self.fst_vt['SubDyn']['NReact']
        self.fst_vt['SubDyn']['RctTDYss'] = [None]*self.fst_vt['SubDyn']['NReact']
        self.fst_vt['SubDyn']['RctTDZss'] = [None]*self.fst_vt['SubDyn']['NReact']
        self.fst_vt['SubDyn']['RctRDXss'] = [None]*self.fst_vt['SubDyn']['NReact']
        self.fst_vt['SubDyn']['RctRDYss'] = [None]*self.fst_vt['SubDyn']['NReact']
        self.fst_vt['SubDyn']['RctRDZss'] = [None]*self.fst_vt['SubDyn']['NReact']
        self.fst_vt['SubDyn']['Rct_SoilFile'] = [None]*self.fst_vt['SubDyn']['NReact']
        ln = f.readline().split()
        ln = f.readline().split()
        for i in range(self.fst_vt['SubDyn']['NReact']):
            ln = f.readline().split()
            self.fst_vt['SubDyn']['RJointID'][i] = int(ln[0])
            self.fst_vt['SubDyn']['RctTDXss'][i] = int(ln[1])
            self.fst_vt['SubDyn']['RctTDYss'][i] = int(ln[2])
            self.fst_vt['SubDyn']['RctTDZss'][i] = int(ln[3])
            self.fst_vt['SubDyn']['RctRDXss'][i] = int(ln[4])
            self.fst_vt['SubDyn']['RctRDYss'][i] = int(ln[5])
            self.fst_vt['SubDyn']['RctRDZss'][i] = int(ln[6])
            if len(ln) == 8:
                self.fst_vt['SubDyn']['Rct_SoilFile'][i] = ln[7]
            else:
                self.fst_vt['SubDyn']['Rct_SoilFile'][i] = 'None'
        f.readline()
        # INTERFACE JOINTS
        self.fst_vt['SubDyn']['NInterf']   = int_read(f.readline().split()[0])
        self.fst_vt['SubDyn']['IJointID'] = [None]*self.fst_vt['SubDyn']['NInterf']
        self.fst_vt['SubDyn']['ItfTDXss'] = [None]*self.fst_vt['SubDyn']['NInterf']
        self.fst_vt['SubDyn']['ItfTDYss'] = [None]*self.fst_vt['SubDyn']['NInterf']
        self.fst_vt['SubDyn']['ItfTDZss'] = [None]*self.fst_vt['SubDyn']['NInterf']
        self.fst_vt['SubDyn']['ItfRDXss'] = [None]*self.fst_vt['SubDyn']['NInterf']
        self.fst_vt['SubDyn']['ItfRDYss'] = [None]*self.fst_vt['SubDyn']['NInterf']
        self.fst_vt['SubDyn']['ItfRDZss'] = [None]*self.fst_vt['SubDyn']['NInterf']
        ln = f.readline().split()
        ln = f.readline().split()
        for i in range(self.fst_vt['SubDyn']['NInterf']):
            ln = f.readline().split()
            self.fst_vt['SubDyn']['IJointID'][i] = int(ln[0])
            self.fst_vt['SubDyn']['ItfTDXss'][i] = int(ln[1])
            self.fst_vt['SubDyn']['ItfTDYss'][i] = int(ln[2])
            self.fst_vt['SubDyn']['ItfTDZss'][i] = int(ln[3])
            self.fst_vt['SubDyn']['ItfRDXss'][i] = int(ln[4])
            self.fst_vt['SubDyn']['ItfRDYss'][i] = int(ln[5])
            self.fst_vt['SubDyn']['ItfRDZss'][i] = int(ln[6])
        f.readline()
        # MEMBERS
        self.fst_vt['SubDyn']['NMembers']    = int_read(f.readline().split()[0])
        self.fst_vt['SubDyn']['MemberID']    = [None]*self.fst_vt['SubDyn']['NMembers']
        self.fst_vt['SubDyn']['MJointID1']   = [None]*self.fst_vt['SubDyn']['NMembers']
        self.fst_vt['SubDyn']['MJointID2']   = [None]*self.fst_vt['SubDyn']['NMembers']
        self.fst_vt['SubDyn']['MPropSetID1'] = [None]*self.fst_vt['SubDyn']['NMembers']
        self.fst_vt['SubDyn']['MPropSetID2'] = [None]*self.fst_vt['SubDyn']['NMembers']
        self.fst_vt['SubDyn']['MType']       = [None]*self.fst_vt['SubDyn']['NMembers']
        self.fst_vt['SubDyn']['COSMID']      = [None]*self.fst_vt['SubDyn']['NMembers']
        ln = f.readline().split()
        ln = f.readline().split()
        for i in range(self.fst_vt['SubDyn']['NMembers']):
            ln = f.readline().split()
            self.fst_vt['SubDyn']['MemberID'][i]    = int(ln[0])
            self.fst_vt['SubDyn']['MJointID1'][i]   = int(ln[1])
            self.fst_vt['SubDyn']['MJointID2'][i]   = int(ln[2])
            self.fst_vt['SubDyn']['MPropSetID1'][i] = int(ln[3])
            self.fst_vt['SubDyn']['MPropSetID2'][i] = int(ln[4])
            self.fst_vt['SubDyn']['MType'][i]       = int(ln[5])
            if len(ln) > 6:
                self.fst_vt['SubDyn']['COSMID'][i]  = int(ln[6])
        f.readline()
        # MEMBER X-SECTION PROPERTY data 1/2
        self.fst_vt['SubDyn']['NPropSets'] = int_read(f.readline().split()[0])
        self.fst_vt['SubDyn']['PropSetID1'] = [None]*self.fst_vt['SubDyn']['NPropSets']
        self.fst_vt['SubDyn']['YoungE1']    = [None]*self.fst_vt['SubDyn']['NPropSets']
        self.fst_vt['SubDyn']['ShearG1']    = [None]*self.fst_vt['SubDyn']['NPropSets']
        self.fst_vt['SubDyn']['MatDens1']   = [None]*self.fst_vt['SubDyn']['NPropSets']
        self.fst_vt['SubDyn']['XsecD']     = [None]*self.fst_vt['SubDyn']['NPropSets']
        self.fst_vt['SubDyn']['XsecT']     = [None]*self.fst_vt['SubDyn']['NPropSets']
        ln = f.readline().split()
        ln = f.readline().split()
        for i in range(self.fst_vt['SubDyn']['NPropSets']):
            ln = f.readline().split()
            self.fst_vt['SubDyn']['PropSetID1'][i] = int(ln[0])
            self.fst_vt['SubDyn']['YoungE1'][i]    = float(ln[1])
            self.fst_vt['SubDyn']['ShearG1'][i]    = float(ln[2])
            self.fst_vt['SubDyn']['MatDens1'][i]   = float(ln[3])
            self.fst_vt['SubDyn']['XsecD'][i]     = float(ln[4])
            self.fst_vt['SubDyn']['XsecT'][i]     = float(ln[5])
        f.readline()
        # MEMBER X-SECTION PROPERTY data 2/2
        self.fst_vt['SubDyn']['NXPropSets'] = int_read(f.readline().split()[0])
        self.fst_vt['SubDyn']['PropSetID2']  = [None]*self.fst_vt['SubDyn']['NXPropSets']
        self.fst_vt['SubDyn']['YoungE2']     = [None]*self.fst_vt['SubDyn']['NXPropSets']
        self.fst_vt['SubDyn']['ShearG2']     = [None]*self.fst_vt['SubDyn']['NXPropSets']
        self.fst_vt['SubDyn']['MatDens2']    = [None]*self.fst_vt['SubDyn']['NXPropSets']
        self.fst_vt['SubDyn']['XsecA']      = [None]*self.fst_vt['SubDyn']['NXPropSets']
        self.fst_vt['SubDyn']['XsecAsx']    = [None]*self.fst_vt['SubDyn']['NXPropSets']
        self.fst_vt['SubDyn']['XsecAsy']    = [None]*self.fst_vt['SubDyn']['NXPropSets']
        self.fst_vt['SubDyn']['XsecJxx']    = [None]*self.fst_vt['SubDyn']['NXPropSets']
        self.fst_vt['SubDyn']['XsecJyy']    = [None]*self.fst_vt['SubDyn']['NXPropSets']
        self.fst_vt['SubDyn']['XsecJ0']     = [None]*self.fst_vt['SubDyn']['NXPropSets']
        ln = f.readline().split()
        ln = f.readline().split()
        for i in range(self.fst_vt['SubDyn']['NXPropSets']):
            ln = f.readline().split()
            self.fst_vt['SubDyn']['PropSetID2'][i] = int(ln[0])
            self.fst_vt['SubDyn']['YoungE2'][i]    = float(ln[1])
            self.fst_vt['SubDyn']['ShearG2'][i]    = float(ln[2])
            self.fst_vt['SubDyn']['MatDens2'][i]   = float(ln[3])
            self.fst_vt['SubDyn']['XsecA'][i]      = float(ln[4])
            self.fst_vt['SubDyn']['XsecAsx'][i]    = float(ln[5])
            self.fst_vt['SubDyn']['XsecAsy'][i]    = float(ln[6])
            self.fst_vt['SubDyn']['XsecJxx'][i]    = float(ln[7])
            self.fst_vt['SubDyn']['XsecJyy'][i]    = float(ln[8])
            self.fst_vt['SubDyn']['XsecJ0'][i]     = float(ln[9])
        # CABLE PROPERTIES
        f.readline()
        self.fst_vt['SubDyn']['NCablePropSets'] = int_read(f.readline().split()[0])
        self.fst_vt['SubDyn']['CablePropSetID'] = [None]*self.fst_vt['SubDyn']['NCablePropSets']
        self.fst_vt['SubDyn']['CableEA']        = [None]*self.fst_vt['SubDyn']['NCablePropSets']
        self.fst_vt['SubDyn']['CableMatDens']   = [None]*self.fst_vt['SubDyn']['NCablePropSets']
        self.fst_vt['SubDyn']['CableT0']        = [None]*self.fst_vt['SubDyn']['NCablePropSets']
        f.readline()
        f.readline()
        for i in range(self.fst_vt['SubDyn']['NCablePropSets']):
            ln = f.readline().split()
            self.fst_vt['SubDyn']['CablePropSetID'][i] = int(ln[0])
            self.fst_vt['SubDyn']['CableEA'][i]        = float(ln[1])
            self.fst_vt['SubDyn']['CableMatDens'][i]   = float(ln[1])
            self.fst_vt['SubDyn']['CableT0'][i]        = float(ln[1])
        # RIGID LINK PROPERTIES
        f.readline()
        self.fst_vt['SubDyn']['NRigidPropSets'] = int_read(f.readline().split()[0])
        self.fst_vt['SubDyn']['RigidPropSetID'] = [None]*self.fst_vt['SubDyn']['NRigidPropSets']
        self.fst_vt['SubDyn']['RigidMatDens']   = [None]*self.fst_vt['SubDyn']['NRigidPropSets']
        f.readline()
        f.readline()
        for i in range(self.fst_vt['SubDyn']['NRigidPropSets']):
            ln = f.readline().split()
            self.fst_vt['SubDyn']['RigidPropSetID'][i] = int(ln[0])
            self.fst_vt['SubDyn']['RigidMatDens'][i]   = float(ln[1])
        # MEMBER COSINE MATRICES
        f.readline()
        self.fst_vt['SubDyn']['NCOSMs'] = int_read(f.readline().split()[0])
        self.fst_vt['SubDyn']['COSMID'] = [None]*self.fst_vt['SubDyn']['NCOSMs']
        self.fst_vt['SubDyn']['COSM11'] = [None]*self.fst_vt['SubDyn']['NCOSMs']
        self.fst_vt['SubDyn']['COSM12'] = [None]*self.fst_vt['SubDyn']['NCOSMs']
        self.fst_vt['SubDyn']['COSM13'] = [None]*self.fst_vt['SubDyn']['NCOSMs']
        self.fst_vt['SubDyn']['COSM21'] = [None]*self.fst_vt['SubDyn']['NCOSMs']
        self.fst_vt['SubDyn']['COSM22'] = [None]*self.fst_vt['SubDyn']['NCOSMs']
        self.fst_vt['SubDyn']['COSM23'] = [None]*self.fst_vt['SubDyn']['NCOSMs']
        self.fst_vt['SubDyn']['COSM31'] = [None]*self.fst_vt['SubDyn']['NCOSMs']
        self.fst_vt['SubDyn']['COSM32'] = [None]*self.fst_vt['SubDyn']['NCOSMs']
        self.fst_vt['SubDyn']['COSM33'] = [None]*self.fst_vt['SubDyn']['NCOSMs']
        ln = f.readline().split()
        ln = f.readline().split()
        for i in range(self.fst_vt['SubDyn']['NCOSMs']):
            ln = f.readline().split()
            self.fst_vt['SubDyn']['COSMID'][i] = int(ln[0])
            self.fst_vt['SubDyn']['COSM11'][i] = float(ln[1])
            self.fst_vt['SubDyn']['COSM12'][i] = float(ln[2])
            self.fst_vt['SubDyn']['COSM13'][i] = float(ln[3])
            self.fst_vt['SubDyn']['COSM21'][i] = float(ln[4])
            self.fst_vt['SubDyn']['COSM22'][i] = float(ln[5])
            self.fst_vt['SubDyn']['COSM23'][i] = float(ln[6])
            self.fst_vt['SubDyn']['COSM31'][i] = float(ln[7])
            self.fst_vt['SubDyn']['COSM32'][i] = float(ln[8])
            self.fst_vt['SubDyn']['COSM33'][i] = float(ln[9])
        f.readline()
        # JOINT ADDITIONAL CONCENTRATED MASSES
        self.fst_vt['SubDyn']['NCmass']    = int_read(f.readline().split()[0])
        self.fst_vt['SubDyn']['CMJointID'] = [None]*self.fst_vt['SubDyn']['NCmass']
        self.fst_vt['SubDyn']['JMass']     = [None]*self.fst_vt['SubDyn']['NCmass']
        self.fst_vt['SubDyn']['JMXX']      = [None]*self.fst_vt['SubDyn']['NCmass']
        self.fst_vt['SubDyn']['JMYY']      = [None]*self.fst_vt['SubDyn']['NCmass']
        self.fst_vt['SubDyn']['JMZZ']      = [None]*self.fst_vt['SubDyn']['NCmass']
        self.fst_vt['SubDyn']['JMXY']      = [None]*self.fst_vt['SubDyn']['NCmass']
        self.fst_vt['SubDyn']['JMXZ']      = [None]*self.fst_vt['SubDyn']['NCmass']
        self.fst_vt['SubDyn']['JMYZ']      = [None]*self.fst_vt['SubDyn']['NCmass']
        self.fst_vt['SubDyn']['MCGX']      = [None]*self.fst_vt['SubDyn']['NCmass']
        self.fst_vt['SubDyn']['MCGY']      = [None]*self.fst_vt['SubDyn']['NCmass']
        self.fst_vt['SubDyn']['MCGZ']      = [None]*self.fst_vt['SubDyn']['NCmass']
        ln = f.readline().split()
        ln = f.readline().split()
        for i in range(self.fst_vt['SubDyn']['NCmass']):
            ln = f.readline().split()
            self.fst_vt['SubDyn']['CMJointID'][i] = int(ln[0])
            self.fst_vt['SubDyn']['JMass'][i]     = float(ln[1])
            self.fst_vt['SubDyn']['JMXX'][i]      = float(ln[2])
            self.fst_vt['SubDyn']['JMYY'][i]      = float(ln[3])
            self.fst_vt['SubDyn']['JMZZ'][i]      = float(ln[4])
            self.fst_vt['SubDyn']['JMXY'][i]      = float(ln[5])
            self.fst_vt['SubDyn']['JMXZ'][i]      = float(ln[6])
            self.fst_vt['SubDyn']['JMYZ'][i]      = float(ln[7])
            self.fst_vt['SubDyn']['MCGX'][i]      = float(ln[8])
            self.fst_vt['SubDyn']['MCGY'][i]      = float(ln[9])
            self.fst_vt['SubDyn']['MCGZ'][i]      = float(ln[10])
        f.readline()
        # OUTPUT
        self.fst_vt['SubDyn']['SumPrint'] = bool_read(f.readline().split()[0])
        self.fst_vt['SubDyn']['OutCOSM']  = bool_read(f.readline().split()[0])
        self.fst_vt['SubDyn']['OutAll']   = bool_read(f.readline().split()[0])
        self.fst_vt['SubDyn']['OutSwtch'] = int_read(f.readline().split()[0])
        self.fst_vt['SubDyn']['TabDelim'] = bool_read(f.readline().split()[0])
        self.fst_vt['SubDyn']['OutDec']   = int_read(f.readline().split()[0])
        self.fst_vt['SubDyn']['OutFmt']   = f.readline().split()[0]
        self.fst_vt['SubDyn']['OutSFmt']  = f.readline().split()[0]
        f.readline()
        # MEMBER OUTPUT LIST
        self.fst_vt['SubDyn']['NMOutputs']     = int_read(f.readline().split()[0])
        self.fst_vt['SubDyn']['MemberID_out']  = [None]*self.fst_vt['SubDyn']['NMOutputs']
        self.fst_vt['SubDyn']['NOutCnt']       = [None]*self.fst_vt['SubDyn']['NMOutputs']
        self.fst_vt['SubDyn']['NodeCnt']       = [None]*self.fst_vt['SubDyn']['NMOutputs']
        ln = f.readline().split()
        ln = f.readline().split()
        for i in range(self.fst_vt['SubDyn']['NMOutputs']):
            ln = f.readline().split()
            self.fst_vt['SubDyn']['MemberID_out'][i] = int(ln[0])
            self.fst_vt['SubDyn']['NOutCnt'][i]      = int(ln[1])
            self.fst_vt['SubDyn']['NodeCnt'][i]      = int(ln[2])
        f.readline()
        # SSOutList
        data = f.readline()
        while data.split()[0] != 'END':
            channels = data.split('"')
            channel_list = channels[1].split(',')
            self.set_outlist(self.fst_vt['outlist']['SubDyn'], channel_list)
            data = f.readline()
            
        f.close()

        f.close()

    def read_MAP(self, map_file):
        # MAP++

        # TODO: this is likely not robust enough, only tested on the Hywind Spar
        # additional lines in these tables are likely

        f = open(map_file)
        f.readline()
        f.readline()
        f.readline()
        data_line = f.readline().strip().split()
        self.fst_vt['MAP']['LineType']     = str(data_line[0])
        self.fst_vt['MAP']['Diam']         = float(data_line[1])
        self.fst_vt['MAP']['MassDenInAir'] = float(data_line[2])
        self.fst_vt['MAP']['EA']           = float(data_line[3])
        self.fst_vt['MAP']['CB']           = float(data_line[4])
        self.fst_vt['MAP']['CIntDamp']     = float(data_line[5])
        self.fst_vt['MAP']['Ca']           = float(data_line[6])
        self.fst_vt['MAP']['Cdn']          = float(data_line[7])
        self.fst_vt['MAP']['Cdt']          = float(data_line[8])
        f.readline()
        f.readline()
        f.readline()
        for i in range(2):
            data_node = f.readline().strip().split()
            self.fst_vt['MAP']['Node'].append(int(data_node[0]))
            self.fst_vt['MAP']['Type'].append(str(data_node[1]))
            self.fst_vt['MAP']['X'].append(float_read(data_node[2]))
            self.fst_vt['MAP']['Y'].append(float_read(data_node[3]))
            self.fst_vt['MAP']['Z'].append(float_read(data_node[4]))
            self.fst_vt['MAP']['M'].append(float_read(data_node[5]))
            self.fst_vt['MAP']['B'].append(float_read(data_node[6]))
            self.fst_vt['MAP']['FX'].append(float_read(data_node[7]))
            self.fst_vt['MAP']['FY'].append(float_read(data_node[8]))
            self.fst_vt['MAP']['FZ'].append(float_read(data_node[9]))
        f.readline()
        f.readline()
        f.readline()
        data_line_prop = f.readline().strip().split()
        self.fst_vt['MAP']['Line']     = int(data_line_prop[0])
        self.fst_vt['MAP']['LineType'] = str(data_line_prop[1])
        self.fst_vt['MAP']['UnstrLen'] = float(data_line_prop[2])
        self.fst_vt['MAP']['NodeAnch'] = int(data_line_prop[3])
        self.fst_vt['MAP']['NodeFair'] = int(data_line_prop[4])
        self.fst_vt['MAP']['Flags']    = [str(val) for val in data_line_prop[5:]]
        f.readline()
        f.readline()
        f.readline()
        self.fst_vt['MAP']['Option']   = [str(val) for val in f.readline().strip().split()]
        f.close()

    def read_MoorDyn(self, moordyn_file):

        f = open(moordyn_file)

        # MoorDyn
        f.readline()
        f.readline()
        self.fst_vt['MoorDyn']['Echo']     = bool_read(f.readline().split()[0])
        f.readline()
        f.readline()
        f.readline()
        self.fst_vt['MoorDyn']['Name'] = []
        self.fst_vt['MoorDyn']['Diam']     = []
        self.fst_vt['MoorDyn']['MassDen']  = []
        self.fst_vt['MoorDyn']['EA']       = []
        self.fst_vt['MoorDyn']['BA_zeta']  = []
        self.fst_vt['MoorDyn']['EI']  = []
        self.fst_vt['MoorDyn']['Cd']      = []
        self.fst_vt['MoorDyn']['Ca']      = []
        self.fst_vt['MoorDyn']['CdAx']      = []
        self.fst_vt['MoorDyn']['CaAx']      = []
        data_line = f.readline().strip().split()
        while data_line[0][:3] != '---': # OpenFAST searches for ---, so we'll do the same
            self.fst_vt['MoorDyn']['Name'].append(str(data_line[0]))
            self.fst_vt['MoorDyn']['Diam'].append(float(data_line[1]))
            self.fst_vt['MoorDyn']['MassDen'].append(float(data_line[2]))
            self.fst_vt['MoorDyn']['EA'].append(float(data_line[3]))
            self.fst_vt['MoorDyn']['BA_zeta'].append(float(data_line[4]))
            self.fst_vt['MoorDyn']['EI'].append(float(data_line[5]))
            self.fst_vt['MoorDyn']['Cd'].append(float(data_line[6]))
            self.fst_vt['MoorDyn']['Ca'].append(float(data_line[7]))
            self.fst_vt['MoorDyn']['CdAx'].append(float(data_line[8]))
            self.fst_vt['MoorDyn']['CaAx'].append(float(data_line[9]))
            data_line = f.readline().strip().split()
        f.readline()
        f.readline()
        self.fst_vt['MoorDyn']['Point_ID'] = []
        self.fst_vt['MoorDyn']['Attachment'] = []
        self.fst_vt['MoorDyn']['X']    = []
        self.fst_vt['MoorDyn']['Y']    = []
        self.fst_vt['MoorDyn']['Z']    = []
        self.fst_vt['MoorDyn']['M']    = []
        self.fst_vt['MoorDyn']['V']    = []
        self.fst_vt['MoorDyn']['CdA']  = []
        self.fst_vt['MoorDyn']['CA']   = []
        data_line = f.readline().strip().split()
        while data_line[0][:3] != '---': # OpenFAST searches for ---, so we'll do the same
            self.fst_vt['MoorDyn']['Point_ID'].append(int(data_line[0]))
            self.fst_vt['MoorDyn']['Attachment'].append(str(data_line[1]))
            self.fst_vt['MoorDyn']['X'].append(float(data_line[2]))
            self.fst_vt['MoorDyn']['Y'].append(float(data_line[3]))
            self.fst_vt['MoorDyn']['Z'].append(float(data_line[4]))
            self.fst_vt['MoorDyn']['M'].append(float(data_line[5]))
            self.fst_vt['MoorDyn']['V'].append(float(data_line[6]))
            self.fst_vt['MoorDyn']['CdA'].append(float(data_line[7]))
            self.fst_vt['MoorDyn']['CA'].append(float(data_line[8]))
            data_line = f.readline().strip().split()
        f.readline()
        f.readline()
        self.fst_vt['MoorDyn']['Line_ID']          = []
        self.fst_vt['MoorDyn']['LineType']    = []
        self.fst_vt['MoorDyn']['AttachA']     = []
        self.fst_vt['MoorDyn']['AttachB']     = []
        self.fst_vt['MoorDyn']['UnstrLen']    = []
        self.fst_vt['MoorDyn']['NumSegs']     = []
        self.fst_vt['MoorDyn']['Outputs']     = []
        data_line = f.readline().strip().split()
        while data_line[0][:3] != '---': # OpenFAST searches for ---, so we'll do the same
            self.fst_vt['MoorDyn']['Line_ID'].append(int(data_line[0]))
            self.fst_vt['MoorDyn']['LineType'].append(str(data_line[1]))
            self.fst_vt['MoorDyn']['AttachA'].append(int(data_line[2]))
            self.fst_vt['MoorDyn']['AttachB'].append(int(data_line[3]))
            self.fst_vt['MoorDyn']['UnstrLen'].append(float(data_line[4]))
            self.fst_vt['MoorDyn']['NumSegs'].append(int(data_line[5]))
            self.fst_vt['MoorDyn']['Outputs'].append(str(data_line[6]))
            data_line = f.readline().strip().split()

        # read optional control inputs, there are other optional MoorDyn sections/inputs
        self.fst_vt['MoorDyn']['ChannelID'] = []
        self.fst_vt['MoorDyn']['Lines_Control'] = []
        if 'CONTROL' in [dl.upper() for dl in data_line]:
            f.readline()
            f.readline()
            data_line = f.readline().strip().split()
            while data_line[0][:3] != '---': # OpenFAST searches for ---, so we'll do the same
                self.fst_vt['MoorDyn']['ChannelID'].append(int(data_line[0]))
                # Line(s) is a list of mooring lines, spaces are allowed between commas
                control_lines = []
                for lines in data_line[1:]:
                    for line in lines.split(','):
                        control_lines.append(line.strip(','))

                # Spaces show up in control_lines as '', remove them all
                while '' in control_lines:
                    control_lines.remove('')

                self.fst_vt['MoorDyn']['Lines_Control'].append(control_lines)
                data_line = f.readline().strip().split()

        # Solver options, there are a few more optional MoorDyn inputs that can be added line 'CONTROL'
        self.fst_vt['MoorDyn']['dtM']       = float_read(f.readline().split()[0])
        self.fst_vt['MoorDyn']['kbot']      = float_read(f.readline().split()[0])
        self.fst_vt['MoorDyn']['cbot']      = float_read(f.readline().split()[0])
        self.fst_vt['MoorDyn']['dtIC']      = float_read(f.readline().split()[0])
        self.fst_vt['MoorDyn']['TmaxIC']    = float_read(f.readline().split()[0])
        self.fst_vt['MoorDyn']['CdScaleIC'] = float_read(f.readline().split()[0])
        self.fst_vt['MoorDyn']['threshIC']  = float_read(f.readline().split()[0])
        f.readline()

        data = f.readline()
        while data.split()[0] != 'END':
            channels = data.strip().strip('"').strip("'")
            channel_list = channels.split(',')
            self.set_outlist(self.fst_vt['outlist']['MoorDyn'], channel_list)
            data = f.readline()

        f.close()

    def execute(self):
          
        self.read_MainInput()
        ed_file = os.path.join(self.FAST_directory, self.fst_vt['Fst']['EDFile'])
        self.read_ElastoDyn(ed_file)
        if not os.path.isabs(self.fst_vt['ElastoDyn']['BldFile1']):
            ed_blade_file = os.path.join(os.path.dirname(ed_file), self.fst_vt['ElastoDyn']['BldFile1'])
        self.read_ElastoDynBlade(ed_blade_file)
        if not os.path.isabs(self.fst_vt['ElastoDyn']['TwrFile']):
            ed_tower_file = os.path.join(os.path.dirname(ed_file), self.fst_vt['ElastoDyn']['TwrFile'])
        self.read_ElastoDynTower(ed_tower_file)
        self.read_InflowWind()
        # AeroDyn version selection
        if self.fst_vt['Fst']['CompAero'] == 1:
            self.read_AeroDyn14()
        elif self.fst_vt['Fst']['CompAero'] == 2:
            self.read_AeroDyn15()
            
        if self.fst_vt['Fst']['CompServo'] == 1:
            self.read_ServoDyn()
            # Read StC Files
            for StC_file in self.fst_vt['ServoDyn']['BStCfiles']:
                self.fst_vt['BStC'].append(self.read_StC(StC_file))
            for StC_file in self.fst_vt['ServoDyn']['NStCfiles']:
                self.fst_vt['NStC'].append(self.read_StC(StC_file))
            for StC_file in self.fst_vt['ServoDyn']['TStCfiles']:
                self.fst_vt['TStC'].append(self.read_StC(StC_file))
            for StC_file in self.fst_vt['ServoDyn']['SStCfiles']:
                self.fst_vt['SStC'].append(self.read_StC(StC_file))
            if ROSCO:
                self.read_DISCON_in()
        hd_file = os.path.normpath(os.path.join(self.FAST_directory, self.fst_vt['Fst']['HydroFile']))
        if os.path.isfile(hd_file): 
            self.read_HydroDyn(hd_file)
        sd_file = os.path.normpath(os.path.join(self.FAST_directory, self.fst_vt['Fst']['SubFile']))
        if os.path.isfile(sd_file): 
            self.read_SubDyn(sd_file)
        if self.fst_vt['Fst']['CompMooring'] == 1: # only MAP++ implemented for mooring models
            map_file = os.path.normpath(os.path.join(self.FAST_directory, self.fst_vt['Fst']['MooringFile']))
            if os.path.isfile(map_file):
                self.read_MAP(map_file)
        if self.fst_vt['Fst']['CompMooring'] == 3: # MoorDyn implimented
            moordyn_file = os.path.normpath(os.path.join(self.FAST_directory, self.fst_vt['Fst']['MooringFile']))
            if os.path.isfile(moordyn_file):
                self.read_MoorDyn(moordyn_file)
        bd_file = os.path.normpath(os.path.join(self.FAST_directory, self.fst_vt['Fst']['BDBldFile(1)']))
        if os.path.isfile(bd_file):
            self.read_BeamDyn(bd_file)

if __name__=="__main__":
    
    examples_dir = os.path.dirname( os.path.dirname( os.path.dirname( os.path.realpath(__file__) ) ) ) + os.sep

    fast = InputReader_OpenFAST()
    fast.FAST_InputFile = 'IEA-15-240-RWT-UMaineSemi.fst'   # FAST input file (ext=.fst)
    fast.FAST_directory = os.path.join(examples_dir, 'examples', '01_aeroelasticse',
                                                     'OpenFAST_models', 'IEA-15-240-RWT',
                                                     'IEA-15-240-RWT-UMaineSemi')   # Path to fst directory files
    fast.execute()
