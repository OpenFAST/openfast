import os
import copy
import random
import time
import operator
import numpy as np
from functools import reduce


try:
    from rosco.toolbox import utilities as ROSCO_utilities
    ROSCO = True
except:
    ROSCO = False


def auto_format(f, var):
    """
    Error handling for variables with 'Default' options

    args:
    f: file object
    var: variable to write to file

    """
    if isinstance(var, str):
        f.write('{:}\n'.format(var))
    elif isinstance(var, int):
        f.write('{:3}\n'.format(var))
    elif isinstance(var, float):
        f.write('{: 2.15e}\n'.format(var))

def float_default_out(val):
    """
    Formatted float output when 'default' is an option

    args:
    val: value to be formatted

    returns:
    formatted value
    """
    if type(val) is float:
        return '{: 22f}'.format(val)
    else:
        return '{:<22}'.format(val)

def int_default_out(val):
    """
    Formatted int output when 'default' is an option

    args:
    val: value to be formatted

    returns:
    formatted value
    """
    if type(val) is float:
        return '{:<22d}'.format(val)
    else:
        return '{:<22}'.format(val)

def get_dict(vartree, branch):
    """
    Given a list of nested dictionary keys, return the dictionary at that point

    args:
    vartree: dictionary to search
    branch: list of keys to search

    returns:
    dictionary at the specified branch
    """
    return reduce(operator.getitem, branch, vartree)

class InputWriter_OpenFAST(object):
    """ Methods to write OpenFAST input files."""

    def __init__(self):

        self.FAST_namingOut = None    #Master FAST file
        self.FAST_runDirectory = None #Output directory
        self.fst_vt = {}
        self.fst_update = {}

    def update(self, fst_update={}):
        """ Change fast variables based on the user supplied values """
        if fst_update:
            self.fst_update = fst_update

        # recursively loop through fast variable levels and set them to their update values
        def loop_dict(vartree, branch):
            for var in vartree.keys():
                branch_i = copy.copy(branch)
                branch_i.append(var)
                if type(vartree[var]) is dict:
                    loop_dict(vartree[var], branch_i)
                else:
                    # try:
                    get_dict(self.fst_vt, branch_i[:-1])[branch_i[-1]] = get_dict(self.fst_update, branch_i[:-1])[branch_i[-1]]
                    # except:
                        # pass

        # make sure update dictionary is not empty
        if self.fst_update:
            # if update dictionary uses list keys, convert to nested dictionaries
            if type(list(self.fst_update.keys())) in [list, tuple]:
                fst_update = copy.copy(self.fst_update)
                self.fst_update = {}
                for var_list in fst_update.keys():
                    branch = []
                    for i, var in enumerate(var_list[0:-1]):
                        if var not in get_dict(self.fst_update, branch).keys():
                            get_dict(self.fst_update, branch)[var] = {}
                        branch.append(var)

                    get_dict(self.fst_update, branch)[var_list[-1]] = fst_update[var_list]
            else:
                print('WARNING: OpenFAST user settings not correctly applied. Please check the modeling_options.yaml')

            # set fast variables to update values
            loop_dict(self.fst_update, [])

    def get_outlist(self, vartree_head, channel_list=[]):
        """ Loop through a list of output channel names, recursively find values set to True in the nested outlist dict """

        # recursively search nested dictionaries
        def loop_dict(vartree, outlist_i):
            for var in vartree.keys():
                if type(vartree[var]) is dict:
                    loop_dict(vartree[var], outlist_i)
                else:
                    if vartree[var]:
                        outlist_i.append(var)
            return outlist_i

        # if specific outlist branches are not specified, get all
        if not channel_list:
            channel_list = vartree_head.keys()

        # loop through top level of dictionary
        outlist = []
        for var in channel_list:
            var = var.replace(' ', '')
            outlist_i = []
            outlist_i = loop_dict(vartree_head[var], outlist_i)
            if outlist_i:
                outlist.append(sorted(outlist_i))

        return outlist

    def update_outlist(self, channels):
        """ Loop through a list of output channel names, recursively search the nested outlist dict and set to specified value"""
        # 'channels' is a dict of channel names as keys with the boolean value they should be set to

        # given a list of nested dictionary keys, return the dict at that point
        def get_dict(vartree, branch):
            return reduce(operator.getitem, branch, self.fst_vt['outlist'])
        # given a list of nested dictionary keys, set the value of the dict at that point
        def set_dict(vartree, branch, val):
            get_dict(vartree, branch[:-1])[branch[-1]] = val
        # recursively loop through outlist dictionaries to set output channels
        def loop_dict(vartree, search_var, val, branch):
            for var in vartree.keys():
                branch_i = copy.copy(branch)
                branch_i.append(var)
                if type(vartree[var]) is dict:
                    loop_dict(vartree[var], search_var, val, branch_i)
                else:
                    if var == search_var:
                        set_dict(self.fst_vt['outlist'], branch_i, val)

        # loop through outchannels on this line, loop through outlist dicts to set to True
        channel_list = channels.keys()
        for var in channel_list:
            val = channels[var]
            var = var.replace(' ', '')
            loop_dict(self.fst_vt['outlist'], var, val, [])

    def execute(self):
        
        if not os.path.exists(self.FAST_runDirectory):
            os.makedirs(self.FAST_runDirectory)

        if self.fst_vt['Fst']['CompElast'] == 3: # Simplified ElastoDyn
            self.write_SimpleElastoDyn()
        else:
            # If elastodyn blade is being used OR if the blade file exists
            if self.fst_vt['Fst']['CompElast'] == 1 or os.path.isfile(self.fst_vt['ElastoDyn']['BldFile1']):
                
                if isinstance(self.fst_vt['ElastoDynBlade'], list):
                    for i_EDbld, EDbld in enumerate(self.fst_vt['ElastoDynBlade']):
                        self.fst_vt['ElastoDyn']['BldFile%d'%(i_EDbld+1)] = self.FAST_namingOut + '_ElastoDynBlade_%d.dat'%(i_EDbld+1)
                        self.write_ElastoDynBlade(bldInd = i_EDbld)

                elif isinstance(self.fst_vt['ElastoDynBlade'], dict):
                    self.fst_vt['ElastoDyn']['BldFile1'] = self.FAST_namingOut + '_ElastoDynBlade.dat'
                    self.fst_vt['ElastoDyn']['BldFile2'] = self.fst_vt['ElastoDyn']['BldFile1']
                    self.fst_vt['ElastoDyn']['BldFile3'] = self.fst_vt['ElastoDyn']['BldFile1']
                    self.write_ElastoDynBlade()

            self.write_ElastoDynTower()
            self.write_ElastoDyn()
        # self.write_WindWnd()
        if self.fst_vt['Fst']['CompInflow'] == 1:
            self.write_InflowWind()
        if self.fst_vt['Fst']['CompAero'] == 1:
            self.write_AeroDisk()
        elif self.fst_vt['Fst']['CompAero'] == 2:
            self.write_AeroDyn()
        
        if self.fst_vt['Fst']['CompServo'] == 1:
            if 'DISCON_in' in self.fst_vt and ROSCO:
                self.write_DISCON_in()
            self.write_ServoDyn()
            for i_NStC, NStC in enumerate(self.fst_vt['NStC']):
                self.write_StC(NStC,self.fst_vt['ServoDyn']['NStCfiles'][i_NStC])
            for i_BStC, BStC in enumerate(self.fst_vt['BStC']):
                self.write_StC(BStC,self.fst_vt['ServoDyn']['BStCfiles'][i_BStC])
            for i_TStC, TStC in enumerate(self.fst_vt['TStC']):
                self.write_StC(TStC,self.fst_vt['ServoDyn']['TStCfiles'][i_TStC])
            for i_SStC, SStC in enumerate(self.fst_vt['SStC']):
                self.write_StC(SStC,self.fst_vt['ServoDyn']['SStCfiles'][i_SStC])
            if self.fst_vt['ServoDyn']['VSContrl'] == 3: # user-defined from routine UserVSCont refered
                self.write_spd_trq()
        
        if self.fst_vt['Fst']['CompHydro'] == 1:
            self.write_HydroDyn()
        if self.fst_vt['Fst']['CompSeaSt'] == 1:
            self.write_SeaState()
        if self.fst_vt['Fst']['CompSub'] == 1:
            self.write_SubDyn()
        elif self.fst_vt['Fst']['CompSub'] == 2:
            self.write_ExtPtfm()
        if self.fst_vt['Fst']['CompMooring'] == 1:
            self.write_MAP()
        elif self.fst_vt['Fst']['CompMooring'] == 3:
            self.write_MoorDyn()
            if 'options' in self.fst_vt['MoorDyn'] and 'WaterKin' in self.fst_vt['MoorDyn']['options']:
                self.write_WaterKin(os.path.join(self.FAST_runDirectory,self.fst_vt['MoorDyn']['WaterKin_file']))

        if self.fst_vt['Fst']['CompElast'] == 2:

            # look at if the the self.fst_vt['BeamDyn'] is an array, if so, loop through the array
            # if its a dictionary, just write the same BeamDyn file 

            if isinstance(self.fst_vt['BeamDyn'], list):
                for i_BD, BD in enumerate(self.fst_vt['BeamDyn']):

                    self.fst_vt['Fst']['BDBldFile(%d)'%(i_BD+1)] = self.FAST_namingOut + '_BeamDyn_%d.dat'%(i_BD+1)
                    self.write_BeamDyn(bldInd = i_BD)

            elif isinstance(self.fst_vt['BeamDyn'], dict):
                self.fst_vt['Fst']['BDBldFile(1)'] = self.FAST_namingOut + '_BeamDyn.dat'
                self.fst_vt['Fst']['BDBldFile(2)'] = self.fst_vt['Fst']['BDBldFile(1)']
                self.fst_vt['Fst']['BDBldFile(3)'] = self.fst_vt['Fst']['BDBldFile(1)']
                self.write_BeamDyn() 

        self.write_MainInput()

    def write_MainInput(self):
        # Main FAST Input File

        self.FAST_InputFileOut = os.path.join(self.FAST_runDirectory, self.FAST_namingOut+'.fst')

        # Keep simple for now:
        f = open(self.FAST_InputFileOut, 'w')

        # ===== .fst Input File =====

        f.write('------- OpenFAST INPUT FILE -------------------------------------------\n')
        f.write('Generated with OpenFAST_IO\n')
        f.write('---------------------- SIMULATION CONTROL --------------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['Fst']['Echo'], 'Echo', '- Echo input data to <RootName>.ech (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['Fst']['AbortLevel']+'"', 'AbortLevel', '- Error level when simulation should abort (string) {"WARNING", "SEVERE", "FATAL"}\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['TMax'], 'TMax', '- Total run time (s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['DT'], 'DT', '- Recommended module time step (s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['InterpOrder'], 'InterpOrder', '- Interpolation order for input/output time history (-) {1=linear, 2=quadratic}\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['NumCrctn'], 'NumCrctn', '- Number of correction iterations (-) {0=explicit calculation, i.e., no corrections}\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['DT_UJac'], 'DT_UJac', '- Time between calls to get Jacobians (s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['UJacSclFact'], 'UJacSclFact', '- Scaling factor used in Jacobians (-)\n'))
        f.write('---------------------- FEATURE SWITCHES AND FLAGS ------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['CompElast'], 'CompElast', '- Compute structural dynamics (switch) {1=ElastoDyn; 2=ElastoDyn + BeamDyn for blades; 3=Simplified ElastoDyn}\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['CompInflow'], 'CompInflow', '- Compute inflow wind velocities (switch) {0=still air; 1=InflowWind; 2=external from ExtInflow}\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['CompAero'], 'CompAero', '- Compute aerodynamic loads (switch) {0=None; 1=AeroDisk; 2=AeroDyn; 3=ExtLoads}\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['CompServo'], 'CompServo', '- Compute control and electrical-drive dynamics (switch) {0=None; 1=ServoDyn}\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['CompSeaSt'], 'CompSeaSt', '- Compute sea state information (switch) {0=None; 1=SeaState}\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['CompHydro'], 'CompHydro', '- Compute hydrodynamic loads (switch) {0=None; 1=HydroDyn}\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['CompSub'], 'CompSub', '- Compute sub-structural dynamics (switch) {0=None; 1=SubDyn; 2=External Platform MCKF}\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['CompMooring'], 'CompMooring', '- Compute mooring system (switch) {0=None; 1=MAP++; 2=FEAMooring; 3=MoorDyn; 4=OrcaFlex}\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['CompIce'], 'CompIce', '- Compute ice loads (switch) {0=None; 1=IceFloe; 2=IceDyn}\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['MHK'], 'MHK', '- MHK turbine type (switch) {0=Not an MHK turbine; 1=Fixed MHK turbine; 2=Floating MHK turbine}\n'))
        f.write('---------------------- ENVIRONMENTAL CONDITIONS --------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['Gravity'], 'Gravity', '- Gravitational acceleration (m/s^2)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['AirDens'], 'AirDens', '- Air density (kg/m^3)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['WtrDens'], 'WtrDens', '- Water density (kg/m^3)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['KinVisc'], 'KinVisc', '- Kinematic viscosity of working fluid (m^2/s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['SpdSound'], 'SpdSound', '- Speed of sound in working fluid (m/s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['Patm'], 'Patm', '- Atmospheric pressure (Pa) [used only for an MHK turbine cavitation check]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['Pvap'], 'Pvap', '- Vapour pressure of working fluid (Pa) [used only for an MHK turbine cavitation check]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['WtrDpth'], 'WtrDpth', '- Water depth (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['MSL2SWL'], 'MSL2SWL', '- Offset between still-water level and mean sea level (m) [positive upward]\n'))        
        f.write('---------------------- INPUT FILES ---------------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['Fst']['EDFile']+'"', 'EDFile', '- Name of file containing ElastoDyn input parameters (quoted string)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['Fst']['BDBldFile(1)']+'"', 'BDBldFile(1)', '- Name of file containing BeamDyn input parameters for blade 1 (quoted string)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['Fst']['BDBldFile(2)']+'"', 'BDBldFile(2)', '- Name of file containing BeamDyn input parameters for blade 2 (quoted string)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['Fst']['BDBldFile(3)']+'"', 'BDBldFile(3)', '- Name of file containing BeamDyn input parameters for blade 3 (quoted string)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['Fst']['InflowFile']+'"', 'InflowFile', '- Name of file containing inflow wind input parameters (quoted string)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['Fst']['AeroFile']+'"', 'AeroFile', '- Name of file containing aerodynamic input parameters (quoted string)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['Fst']['ServoFile']+'"', 'ServoFile', '- Name of file containing control and electrical-drive input parameters (quoted string)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['Fst']['SeaStFile']+'"', 'SeaStFile', '- Name of file containing sea state input parameters (quoted string)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['Fst']['HydroFile']+'"', 'HydroFile', '- Name of file containing hydrodynamic input parameters (quoted string)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['Fst']['SubFile']+'"', 'SubFile', '- Name of file containing sub-structural input parameters (quoted string)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['Fst']['MooringFile']+'"', 'MooringFile', '- Name of file containing mooring system input parameters (quoted string)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['Fst']['IceFile']+'"', 'IceFile', '- Name of file containing ice input parameters (quoted string)\n'))
        f.write('---------------------- OUTPUT --------------------------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['Fst']['SumPrint'], 'SumPrint', '- Print summary data to "<RootName>.sum" (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['SttsTime'], 'SttsTime', '- Amount of time between screen status messages (s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['ChkptTime'], 'ChkptTime', '- Amount of time between creating checkpoint files for potential restart (s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['DT_Out'], 'DT_Out', '- Time step for tabular output (s) (or "default")\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['TStart'], 'TStart', '- Time to begin tabular output (s)\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['Fst']['OutFileFmt'], 'OutFileFmt', '- Format for tabular (time-marching) output file (switch) {1: text file [<RootName>.out], 2: binary file [<RootName>.outb], 3: both 1 and 2, 4: uncompressed binary [<RootName>.outb], 5: both 1 and 4}\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['Fst']['TabDelim'], 'TabDelim', '- Use tab delimiters in text tabular output file? (flag) {uses spaces if false}\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['Fst']['OutFmt']+'"', 'OutFmt', '- Format used for text tabular output, excluding the time channel.  Resulting field should be 10 characters. (quoted string)\n'))
        f.write('---------------------- LINEARIZATION -------------------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['Fst']['Linearize'],   'Linearize',    '- Linearization analysis (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['Fst']['CalcSteady'],  'CalcSteady',   '- Calculate a steady-state periodic operating point before linearization? [unused if Linearize=False] (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['TrimCase'],      'TrimCase',     '- Controller parameter to be trimmed {1:yaw; 2:torque; 3:pitch} [used only if CalcSteady=True] (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['TrimTol'],       'TrimTol',      '- Tolerance for the rotational speed convergence [used only if CalcSteady=True] (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['TrimGain'],      'TrimGain',     '- Proportional gain for the rotational speed error (>0) [used only if CalcSteady=True] (rad/(rad/s) for yaw or pitch; Nm/(rad/s) for torque)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['Twr_Kdmp'],      'Twr_Kdmp',     '- Damping factor for the tower [used only if CalcSteady=True] (N/(m/s))\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['Bld_Kdmp'],      'Bld_Kdmp',     '- Damping factor for the blades [used only if CalcSteady=True] (N/(m/s))\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['NLinTimes'],     'NLinTimes',    '- Number of times to linearize (-) [>=1] [unused if Linearize=False]\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join(['%f'%i for i in np.array(self.fst_vt['Fst']['LinTimes'], dtype=float)]), 'LinTimes', '- List of times at which to linearize (s) [1 to NLinTimes] [used only when Linearize=True and CalcSteady=False]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['LinInputs'],     'LinInputs',    '- Inputs included in linearization (switch) {0=none; 1=standard; 2=all module inputs (debug)} [unused if Linearize=False]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['LinOutputs'],    'LinOutputs',   '- Outputs included in linearization (switch) {0=none; 1=from OutList(s); 2=all module outputs (debug)} [unused if Linearize=False]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['Fst']['LinOutJac'],   'LinOutJac',    '- Include full Jacobians in linearization output (for debug) (flag) [unused if Linearize=False; used only if LinInputs=LinOutputs=2]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['Fst']['LinOutMod'],   'LinOutMod',    '- Write module-level linearization output files in addition to output for full system? (flag) [unused if Linearize=False]\n'))
        f.write('---------------------- VISUALIZATION ------------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['WrVTK'], 'WrVTK', '- VTK visualization data output: (switch) {0=none; 1=initialization data only; 2=animation; 3=mode shapes}\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['VTK_type'], 'VTK_type', '- Type of VTK visualization data: (switch) {1=surfaces; 2=basic meshes (lines/points); 3=all meshes (debug)} [unused if WrVTK=0]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['Fst']['VTK_fields'], 'VTK_fields', '- Write mesh fields to VTK data files? (flag) {true/false} [unused if WrVTK=0]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['VTK_fps'], 'VTK_fps', '-Frame rate for VTK output (frames per second){will use closest integer multiple of DT} [used only if WrVTK=2 or WrVTK=3]\n'))

        f.flush()
        os.fsync(f)
        f.close()

    def write_ElastoDyn(self):

        self.fst_vt['Fst']['EDFile'] = self.FAST_namingOut + '_ElastoDyn.dat'
        ed_file = os.path.join(self.FAST_runDirectory,self.fst_vt['Fst']['EDFile'])
        f = open(ed_file, 'w')

        f.write('------- ELASTODYN INPUT FILE -------------------------------------------\n')
        f.write('Generated with OpenFAST_IO\n')

        # ElastoDyn Simulation Control (ed_sim_ctrl)
        f.write('---------------------- SIMULATION CONTROL --------------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['Echo'], 'Echo', '- Echo input data to "<RootName>.ech" (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['Method'], 'Method', '- Integration method: {1: RK4, 2: AB4, or 3: ABM4} (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['DT'], 'DT', '- Integration time step (s)\n'))
        f.write('---------------------- DEGREES OF FREEDOM --------------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['FlapDOF1'], 'FlapDOF1', '- First flapwise blade mode DOF (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['FlapDOF2'], 'FlapDOF2', '- Second flapwise blade mode DOF (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['EdgeDOF'], 'EdgeDOF', '- First edgewise blade mode DOF (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['TeetDOF'], 'TeetDOF', '- Rotor-teeter DOF (flag) [unused for 3 blades]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['DrTrDOF'], 'DrTrDOF', '- Drivetrain rotational-flexibility DOF (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['GenDOF'], 'GenDOF', '- Generator DOF (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['YawDOF'], 'YawDOF', '- Yaw DOF (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['TwFADOF1'], 'TwFADOF1', '- First fore-aft tower bending-mode DOF (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['TwFADOF2'], 'TwFADOF2', '- Second fore-aft tower bending-mode DOF (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['TwSSDOF1'], 'TwSSDOF1', '- First side-to-side tower bending-mode DOF (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['TwSSDOF2'], 'TwSSDOF2', '- Second side-to-side tower bending-mode DOF (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['PtfmSgDOF'], 'PtfmSgDOF', '- Platform horizontal surge translation DOF (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['PtfmSwDOF'], 'PtfmSwDOF', '- Platform horizontal sway translation DOF (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['PtfmHvDOF'], 'PtfmHvDOF', '- Platform vertical heave translation DOF (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['PtfmRDOF'], 'PtfmRDOF', '- Platform roll tilt rotation DOF (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['PtfmPDOF'], 'PtfmPDOF', '- Platform pitch tilt rotation DOF (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['PtfmYDOF'], 'PtfmYDOF', '- Platform yaw rotation DOF (flag)\n'))
        f.write('---------------------- INITIAL CONDITIONS --------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['OoPDefl'], 'OoPDefl', '- Initial out-of-plane blade-tip displacement (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['IPDefl'], 'IPDefl', '- Initial in-plane blade-tip deflection (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['BlPitch1'], 'BlPitch(1)', '- Blade 1 initial pitch (degrees)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['BlPitch2'], 'BlPitch(2)', '- Blade 2 initial pitch (degrees)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['BlPitch3'], 'BlPitch(3)', '- Blade 3 initial pitch (degrees) [unused for 2 blades]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['TeetDefl'], 'TeetDefl', '- Initial or fixed teeter angle (degrees) [unused for 3 blades]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['Azimuth'], 'Azimuth', '- Initial azimuth angle for blade 1 (degrees)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['RotSpeed'], 'RotSpeed', '- Initial or fixed rotor speed (rpm)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['NacYaw'], 'NacYaw', '- Initial or fixed nacelle-yaw angle (degrees)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['TTDspFA'], 'TTDspFA', '- Initial fore-aft tower-top displacement (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['TTDspSS'], 'TTDspSS', '- Initial side-to-side tower-top displacement (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['PtfmSurge'], 'PtfmSurge', '- Initial or fixed horizontal surge translational displacement of platform (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['PtfmSway'], 'PtfmSway', '- Initial or fixed horizontal sway translational displacement of platform (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['PtfmHeave'], 'PtfmHeave', '- Initial or fixed vertical heave translational displacement of platform (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['PtfmRoll'], 'PtfmRoll', '- Initial or fixed roll tilt rotational displacement of platform (degrees)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['PtfmPitch'], 'PtfmPitch', '- Initial or fixed pitch tilt rotational displacement of platform (degrees)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['PtfmYaw'], 'PtfmYaw', '- Initial or fixed yaw rotational displacement of platform (degrees)\n'))
        f.write('---------------------- TURBINE CONFIGURATION -----------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['NumBl'], 'NumBl', '- Number of blades (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['TipRad'], 'TipRad', '- The distance from the rotor apex to the blade tip (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['HubRad'], 'HubRad', '- The distance from the rotor apex to the blade root (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['PreCone(1)'], 'PreCone(1)', '- Blade 1 cone angle (degrees)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['PreCone(2)'], 'PreCone(2)', '- Blade 2 cone angle (degrees)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['PreCone(3)'], 'PreCone(3)', '- Blade 3 cone angle (degrees) [unused for 2 blades]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['HubCM'], 'HubCM', '- Distance from rotor apex to hub mass [positive downwind] (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['UndSling'], 'UndSling', '- Undersling length [distance from teeter pin to the rotor apex] (meters) [unused for 3 blades]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['Delta3'], 'Delta3', '- Delta-3 angle for teetering rotors (degrees) [unused for 3 blades]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['AzimB1Up'], 'AzimB1Up', '- Azimuth value to use for I/O when blade 1 points up (degrees)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['OverHang'], 'OverHang', '- Distance from yaw axis to rotor apex [3 blades] or teeter pin [2 blades] (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['ShftGagL'], 'ShftGagL', '- Distance from rotor apex [3 blades] or teeter pin [2 blades] to shaft strain gages [positive for upwind rotors] (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['ShftTilt'], 'ShftTilt', '- Rotor shaft tilt angle (degrees)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['NacCMxn'], 'NacCMxn', '- Downwind distance from the tower-top to the nacelle CM (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['NacCMyn'], 'NacCMyn', '- Lateral  distance from the tower-top to the nacelle CM (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['NacCMzn'], 'NacCMzn', '- Vertical distance from the tower-top to the nacelle CM (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['NcIMUxn'], 'NcIMUxn', '- Downwind distance from the tower-top to the nacelle IMU (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['NcIMUyn'], 'NcIMUyn', '- Lateral  distance from the tower-top to the nacelle IMU (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['NcIMUzn'], 'NcIMUzn', '- Vertical distance from the tower-top to the nacelle IMU (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['Twr2Shft'], 'Twr2Shft', '- Vertical distance from the tower-top to the rotor shaft (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['TowerHt'], 'TowerHt', '- Height of tower relative to ground level [onshore], MSL [offshore wind or floating MHK], or seabed [fixed MHK] (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['TowerBsHt'], 'TowerBsHt', '- Height of tower base relative to ground level [onshore], MSL [offshore wind or floating MHK], or seabed [fixed MHK] (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['PtfmCMxt'], 'PtfmCMxt', '- Downwind distance from the ground level [onshore], MSL [offshore wind or floating MHK], or seabed [fixed MHK] to the platform CM (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['PtfmCMyt'], 'PtfmCMyt', '- Lateral distance from the ground level [onshore], MSL [offshore wind or floating MHK], or seabed [fixed MHK] to the platform CM (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['PtfmCMzt'], 'PtfmCMzt', '- Vertical distance from the ground level [onshore], MSL [offshore wind or floating MHK], or seabed [fixed MHK] to the platform CM (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['PtfmRefzt'], 'PtfmRefzt', '- Vertical distance from the ground level [onshore], MSL [offshore wind or floating MHK], or seabed [fixed MHK] to the platform reference point (meters)\n'))
        f.write('---------------------- MASS AND INERTIA ----------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['TipMass(1)'], 'TipMass(1)', '- Tip-brake mass, blade 1 (kg)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['TipMass(2)'], 'TipMass(2)', '- Tip-brake mass, blade 2 (kg)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['TipMass(3)'], 'TipMass(3)', '- Tip-brake mass, blade 3 (kg) [unused for 2 blades]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['HubMass'], 'HubMass', '- Hub mass (kg)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['HubIner'], 'HubIner', '- Hub inertia about rotor axis [3 blades] or teeter axis [2 blades] (kg m^2)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['GenIner'], 'GenIner', '- Generator inertia about HSS (kg m^2)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['NacMass'], 'NacMass', '- Nacelle mass (kg)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['NacYIner'], 'NacYIner', '- Nacelle inertia about yaw axis (kg m^2)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['YawBrMass'], 'YawBrMass', '- Yaw bearing mass (kg)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['PtfmMass'], 'PtfmMass', '- Platform mass (kg)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['PtfmRIner'], 'PtfmRIner', '- Platform inertia for roll tilt rotation about the platform CM (kg m^2)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['PtfmPIner'], 'PtfmPIner', '- Platform inertia for pitch tilt rotation about the platform CM (kg m^2)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['PtfmYIner'], 'PtfmYIner', '- Platform inertia for yaw rotation about the platform CM (kg m^2)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['PtfmXYIner'], 'PtfmXYIner', '- Platform xy moment of inertia about the platform CM (=-int(xydm)) (kg m^2)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['PtfmYZIner'], 'PtfmYZIner', '- Platform yz moment of inertia about the platform CM (=-int(yzdm)) (kg m^2)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['PtfmXZIner'], 'PtfmXZIner', '- Platform xz moment of inertia about the platform CM (=-int(xzdm)) (kg m^2)\n'))
        f.write('---------------------- BLADE ---------------------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['BldNodes'], 'BldNodes', '- Number of blade nodes (per blade) used for analysis (-)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['ElastoDyn']['BldFile1']+'"', 'BldFile(1)', '- Name of file containing properties for blade 1 (quoted string)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['ElastoDyn']['BldFile2']+'"', 'BldFile(2)', '- Name of file containing properties for blade 2 (quoted string)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['ElastoDyn']['BldFile3']+'"', 'BldFile(3)', '- Name of file containing properties for blade 3 (quoted string) [unused for 2 blades]\n'))
        f.write('---------------------- ROTOR-TEETER --------------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['TeetMod'], 'TeetMod', '- Rotor-teeter spring/damper model {0: none, 1: standard, 2: user-defined from routine UserTeet} (switch) [unused for 3 blades]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['TeetDmpP'], 'TeetDmpP', '- Rotor-teeter damper position (degrees) [used only for 2 blades and when TeetMod=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['TeetDmp'], 'TeetDmp', '- Rotor-teeter damping constant (N-m/(rad/s)) [used only for 2 blades and when TeetMod=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['TeetCDmp'], 'TeetCDmp', '- Rotor-teeter rate-independent Coulomb-damping moment (N-m) [used only for 2 blades and when TeetMod=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['TeetSStP'], 'TeetSStP', '- Rotor-teeter soft-stop position (degrees) [used only for 2 blades and when TeetMod=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['TeetHStP'], 'TeetHStP', '- Rotor-teeter hard-stop position (degrees) [used only for 2 blades and when TeetMod=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['TeetSSSp'], 'TeetSSSp', '- Rotor-teeter soft-stop linear-spring constant (N-m/rad) [used only for 2 blades and when TeetMod=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['TeetHSSp'], 'TeetHSSp', '- Rotor-teeter hard-stop linear-spring constant (N-m/rad) [used only for 2 blades and when TeetMod=1]\n'))
        f.write('---------------------- YAW-FRICTION --------------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['YawFrctMod'], 'YawFrctMod', '- Yaw-friction model {0: none, 1: friction independent of yaw-bearing force and bending moment, 2: friction with Coulomb terms depending on yaw-bearing force and bending moment, 3: user defined model} (switch)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['M_CSmax'], 'M_CSmax', '- Maximum static Coulomb friction torque (N-m) [M_CSmax when YawFrctMod=1; |Fz|*M_CSmax when YawFrctMod=2 and Fz<0]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['M_FCSmax'], 'M_FCSmax', '- Maximum static Coulomb friction torque proportional to yaw bearing shear force (N-m) [sqrt(Fx^2+Fy^2)*M_FCSmax; only used when YawFrctMod=2]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['M_MCSmax'], 'M_MCSmax', '- Maximum static Coulomb friction torque proportional to yaw bearing bending moment (N-m) [sqrt(Mx^2+My^2)*M_MCSmax; only used when YawFrctMod=2]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['M_CD'], 'M_CD', '- Dynamic Coulomb friction moment (N-m) [M_CD when YawFrctMod=1; |Fz|*M_CD when YawFrctMod=2 and Fz<0]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['M_FCD'], 'M_FCD', '- Dynamic Coulomb friction moment proportional to yaw bearing shear force (N-m) [sqrt(Fx^2+Fy^2)*M_FCD; only used when YawFrctMod=2]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['M_MCD'], 'M_MCD', '- Dynamic Coulomb friction moment proportional to yaw bearing bending moment (N-m) [sqrt(Mx^2+My^2)*M_MCD; only used when YawFrctMod=2]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['sig_v'], 'sig_v', '- Linear viscous friction coefficient (N-m/(rad/s))\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['sig_v2'], 'sig_v2', '- Quadratic viscous friction coefficient (N-m/(rad/s)^2)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['OmgCut'], 'OmgCut', '- Yaw angular velocity cutoff below which viscous friction is linearized (rad/s)\n'))
        f.write('---------------------- DRIVETRAIN ----------------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['GBoxEff'], 'GBoxEff', '- Gearbox efficiency (%)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['GBRatio'], 'GBRatio', '- Gearbox ratio (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['DTTorSpr'], 'DTTorSpr', '- Drivetrain torsional spring (N-m/rad)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['DTTorDmp'], 'DTTorDmp', '- Drivetrain torsional damper (N-m/(rad/s))\n'))
        f.write('---------------------- FURLING -------------------------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['Furling'], 'Furling', '- Read in additional model properties for furling turbine (flag) [must currently be FALSE)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['ElastoDyn']['FurlFile']+'"', 'FurlFile', '- Name of file containing furling properties (quoted string) [unused when Furling=False]\n'))
        f.write('---------------------- TOWER ---------------------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['TwrNodes'], 'TwrNodes', '- Number of tower nodes used for analysis (-)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['ElastoDyn']['TwrFile']+'"', 'TwrFile', '- Name of file containing tower properties (quoted string)\n'))
        f.write('---------------------- OUTPUT --------------------------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['SumPrint'], 'SumPrint', '- Print summary data to "<RootName>.sum" (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['OutFile'], 'OutFile', '- Switch to determine where output will be placed: {1: in module output file only; 2: in glue code output file only; 3: both} (currently unused)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['TabDelim'], 'TabDelim', '- Use tab delimiters in text tabular output file? (flag) (currently unused)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['ElastoDyn']['OutFmt']+'"', 'OutFmt', '- Format used for text tabular output (except time).  Resulting field should be 10 characters. (quoted string) (currently unused)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['TStart'], 'TStart', '- Time to begin tabular output (s) (currently unused)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['DecFact'], 'DecFact', '- Decimation factor for tabular output {1: output every time step} (-) (currently unused)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['NTwGages'], 'NTwGages', '- Number of tower nodes that have strain gages for output [0 to 9] (-)\n'))  
        if self.fst_vt['ElastoDyn']['TwrGagNd'] != 0:
            f.write('{:<22} {:<11} {:}'.format(', '.join(['%d'%int(i) for i in self.fst_vt['ElastoDyn']['TwrGagNd']]), 'TwrGagNd', '- List of tower nodes that have strain gages [1 to TwrNodes] (-) [unused if NTwGages=0]\n'))
        else:
            f.write('{:<22} {:<11} {:}'.format('', 'TwrGagNd', '- List of tower nodes that have strain gages [1 to TwrNodes] (-) [unused if NTwGages=0]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['NBlGages'], 'NBlGages', '- Number of blade nodes that have strain gages for output [0 to 9] (-)\n'))
        if self.fst_vt['ElastoDyn']['BldGagNd'] != 0:
            f.write('{:<22} {:<11} {:}'.format(', '.join(['%d'%int(i) for i in self.fst_vt['ElastoDyn']['BldGagNd']]), 'BldGagNd', '- List of blade nodes that have strain gages [1 to BldNodes] (-) [unused if NBlGages=0]\n'))
        else:
            f.write('{:<22} {:<11} {:}'.format('', 'BldGagNd', '- List of blade nodes that have strain gages [1 to BldNodes] (-) [unused if NBlGages=0]\n'))
        f.write('                   OutList             - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)\n')

        outlist = self.get_outlist(self.fst_vt['outlist'], ['ElastoDyn'])
        
        for channel_list in outlist:
            for i in range(len(channel_list)):
                f.write('"' + channel_list[i] + '"\n')
        
        f.write('END of OutList section (the word "END" must appear in the first 3 columns of the last OutList line)\n')
                
        # Optional nodal output section
        if 'BldNd_BladesOut' in self.fst_vt['ElastoDyn']:
            f.write('====== Outputs for all blade stations (same ending as above for Spn1.... =========================== [optional section]\n')
            f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['BldNd_BladesOut'], 'BldNd_BladesOut', '- Number of blades to output all node information at.  Up to number of blades on turbine. (-)\n'))
            f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ElastoDyn']['BldNd_BlOutNd'], 'BldNd_BlOutNd', '- Future feature will allow selecting a portion of the nodes to output.  Not implemented yet. (-)\n'))
            f.write('                   OutList     - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx, ElastoDyn_Nodes tab for a listing of available output channels, (-)\n')
            
            opt_outlist = self.get_outlist(self.fst_vt['outlist'], ['ElastoDyn_Nodes'])      
            for opt_channel_list in opt_outlist:
                for i in range(len(opt_channel_list)):
                    f.write('"' + opt_channel_list[i] + '"\n')
            f.write('END (the word "END" must appear in the first 3 columns of this last OutList line in the optional nodal output section)\n')
        
        f.write('---------------------------------------------------------------------------------------\n')
        f.flush()
        os.fsync(f)
        f.close()

    def write_SimpleElastoDyn(self):
        # Write the simple ElastoDyn file

        self.fst_vt['Fst']['EDFile'] = self.FAST_namingOut + '_SimpleElastoDyn.dat'
        sed_file = os.path.join(self.FAST_runDirectory,self.fst_vt['Fst']['EDFile'])
        f = open(sed_file, 'w')

        f.write('------- SIMPLIFIED ELASTODYN INPUT FILE ----------------------------------------\n')
        f.write('Generated with OpenFAST_IO\n')
        f.write('---------------------- SIMULATION CONTROL --------------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['SimpleElastoDyn']['Echo'], 'Echo', '- Echo input data to "<RootName>.ech" (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SimpleElastoDyn']['IntMethod'], 'IntMethod', '- Integration method: {1: RK4, 2: AB4, or 3: ABM4} (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SimpleElastoDyn']['DT'], 'DT', '- Integration time step (s)\n'))
        f.write('---------------------- DEGREES OF FREEDOM --------------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['SimpleElastoDyn']['GenDOF'], 'GenDOF', '- Generator DOF (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['SimpleElastoDyn']['YawDOF'], 'YawDOF', '- Yaw degree of freedom -- controlled by controller (flag)\n'))
        f.write('---------------------- INITIAL CONDITIONS --------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SimpleElastoDyn']['Azimuth'], 'Azimuth', '- Initial azimuth angle for blades (degrees)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SimpleElastoDyn']['BlPitch'], 'BlPitch', '- Blades initial pitch (degrees)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SimpleElastoDyn']['RotSpeed'], 'RotSpeed', '- Initial or fixed rotor speed (rpm)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SimpleElastoDyn']['NacYaw'], 'NacYaw', '- Initial or fixed nacelle-yaw angle (degrees)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SimpleElastoDyn']['PtfmPitch'], 'PtfmPitch', '- Fixed pitch tilt rotational displacement of platform (degrees)\n'))
        f.write('---------------------- TURBINE CONFIGURATION -----------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SimpleElastoDyn']['NumBl'], 'NumBl', '- Number of blades (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SimpleElastoDyn']['TipRad'], 'TipRad', '- The distance from the rotor apex to the blade tip (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SimpleElastoDyn']['HubRad'], 'HubRad', '- The distance from the rotor apex to the blade root (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SimpleElastoDyn']['PreCone'], 'PreCone', '- Blades cone angle (degrees)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SimpleElastoDyn']['OverHang'], 'OverHang', '- Distance from yaw axis to rotor apex [3 blades] or teeter pin [2 blades] (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SimpleElastoDyn']['ShftTilt'], 'ShftTilt', '- Rotor shaft tilt angle (degrees)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SimpleElastoDyn']['Twr2Shft'], 'Twr2Shft', '- Vertical distance from the tower-top to the rotor shaft (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SimpleElastoDyn']['TowerHt'], 'TowerHt', '- Height of tower above ground level [onshore] or MSL [offshore] (meters)\n'))
        f.write('---------------------- MASS AND INERTIA ----------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SimpleElastoDyn']['RotIner'], 'RotIner', '- Rot inertia about rotor axis [blades + hub] (kg m^2)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SimpleElastoDyn']['GenIner'], 'GenIner', '- Generator inertia about HSS (kg m^2)\n'))
        f.write('---------------------- DRIVETRAIN ----------------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SimpleElastoDyn']['GBoxRatio'], 'GBoxRatio', '- Gearbox ratio (-)\n'))
        f.write('---------------------- OUTPUT --------------------------------------------------\n')
        f.write('                   OutList     - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)\n')

        outlist = self.get_outlist(self.fst_vt['outlist'], ['SimpleElastoDyn'])
        for channel_list in outlist:
            for i in range(len(channel_list)):
                f.write('"' + channel_list[i] + '"\n')

        f.write('END of input file (the word "END" must appear in the first 3 columns of the last OutList line)\n')
        f.write('---------------------------------------------------------------------------------------\n')

        f.flush()
        os.fsync(f)
        f.close()

    def write_ElastoDynBlade(self, bldInd = None):

        if bldInd is None:
            EDbld_dict = self.fst_vt['ElastoDynBlade']
            blade_file = os.path.join(self.FAST_runDirectory,self.fst_vt['ElastoDyn']['BldFile1'])
        else:
            EDbld_dict = self.fst_vt['ElastoDynBlade'][bldInd]
            blade_file = os.path.join(self.FAST_runDirectory,self.fst_vt['ElastoDyn']['BldFile'+bldInd])

        f = open(blade_file, 'w')

        f.write('------- ELASTODYN INDIVIDUAL BLADE INPUT FILE --------------------------\n')
        f.write('Generated with OpenFAST_IO\n')
        f.write('---------------------- BLADE PARAMETERS ----------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(EDbld_dict['NBlInpSt'], 'NBlInpSt', '- Number of blade input stations (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(EDbld_dict['BldFlDmp1'], 'BldFlDmp(1)', '- Blade flap mode #1 structural damping in percent of critical (%)\n'))
        f.write('{:<22} {:<11} {:}'.format(EDbld_dict['BldFlDmp2'], 'BldFlDmp(2)', '- Blade flap mode #2 structural damping in percent of critical (%)\n'))
        f.write('{:<22} {:<11} {:}'.format(EDbld_dict['BldEdDmp1'], 'BldEdDmp(1)', '- Blade edge mode #1 structural damping in percent of critical (%)\n'))
        f.write('---------------------- BLADE ADJUSTMENT FACTORS --------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(EDbld_dict['FlStTunr1'], 'FlStTunr(1)', '- Blade flapwise modal stiffness tuner, 1st mode (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(EDbld_dict['FlStTunr2'], 'FlStTunr(2)', '- Blade flapwise modal stiffness tuner, 2nd mode (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(EDbld_dict['AdjBlMs'], 'AdjBlMs', '- Factor to adjust blade mass density (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(EDbld_dict['AdjFlSt'], 'AdjFlSt', '- Factor to adjust blade flap stiffness (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(EDbld_dict['AdjEdSt'], 'AdjEdSt', '- Factor to adjust blade edge stiffness (-)\n'))
        f.write('---------------------- DISTRIBUTED BLADE PROPERTIES ----------------------------\n')
        f.write('    BlFract      PitchAxis      StrcTwst       BMassDen        FlpStff        EdgStff\n')
        f.write('      (-)           (-)          (deg)          (kg/m)         (Nm^2)         (Nm^2)\n')
        BlFract   = EDbld_dict['BlFract']
        PitchAxis = EDbld_dict['PitchAxis']
        StrcTwst  = EDbld_dict['StrcTwst']
        BMassDen  = EDbld_dict['BMassDen']
        FlpStff   = EDbld_dict['FlpStff']
        EdgStff   = EDbld_dict['EdgStff']
        for BlFracti, PitchAxisi, StrcTwsti, BMassDeni, FlpStffi, EdgStffi in zip(BlFract, PitchAxis, StrcTwst, BMassDen, FlpStff, EdgStff):
            f.write('{: 2.15e} {: 2.15e} {: 2.15e} {: 2.15e} {: 2.15e} {: 2.15e}\n'.format(BlFracti, PitchAxisi, StrcTwsti, BMassDeni, FlpStffi, EdgStffi))
        f.write('---------------------- BLADE MODE SHAPES ---------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(EDbld_dict['BldFl1Sh'][0], 'BldFl1Sh(2)', '- Flap mode 1, coeff of x^2\n'))
        f.write('{:<22} {:<11} {:}'.format(EDbld_dict['BldFl1Sh'][1], 'BldFl1Sh(3)', '-            , coeff of x^3\n'))
        f.write('{:<22} {:<11} {:}'.format(EDbld_dict['BldFl1Sh'][2], 'BldFl1Sh(4)', '-            , coeff of x^4\n'))
        f.write('{:<22} {:<11} {:}'.format(EDbld_dict['BldFl1Sh'][3], 'BldFl1Sh(5)', '-            , coeff of x^5\n'))
        f.write('{:<22} {:<11} {:}'.format(EDbld_dict['BldFl1Sh'][4], 'BldFl1Sh(6)', '-            , coeff of x^6\n'))
        f.write('{:<22} {:<11} {:}'.format(EDbld_dict['BldFl2Sh'][0], 'BldFl2Sh(2)', '- Flap mode 2, coeff of x^2\n'))
        f.write('{:<22} {:<11} {:}'.format(EDbld_dict['BldFl2Sh'][1], 'BldFl2Sh(3)', '-            , coeff of x^3\n'))
        f.write('{:<22} {:<11} {:}'.format(EDbld_dict['BldFl2Sh'][2], 'BldFl2Sh(4)', '-            , coeff of x^4\n'))
        f.write('{:<22} {:<11} {:}'.format(EDbld_dict['BldFl2Sh'][3], 'BldFl2Sh(5)', '-            , coeff of x^5\n'))
        f.write('{:<22} {:<11} {:}'.format(EDbld_dict['BldFl2Sh'][4], 'BldFl2Sh(6)', '-            , coeff of x^6\n'))
        f.write('{:<22} {:<11} {:}'.format(EDbld_dict['BldEdgSh'][0], 'BldEdgSh(2)', '- Edge mode 1, coeff of x^2\n'))
        f.write('{:<22} {:<11} {:}'.format(EDbld_dict['BldEdgSh'][1], 'BldEdgSh(3)', '-            , coeff of x^3\n'))
        f.write('{:<22} {:<11} {:}'.format(EDbld_dict['BldEdgSh'][2], 'BldEdgSh(4)', '-            , coeff of x^4\n'))
        f.write('{:<22} {:<11} {:}'.format(EDbld_dict['BldEdgSh'][3], 'BldEdgSh(5)', '-            , coeff of x^5\n'))
        f.write('{:<22} {:<11} {:}'.format(EDbld_dict['BldEdgSh'][4], 'BldEdgSh(6)', '-            , coeff of x^6\n'))      
         
        f.flush()
        os.fsync(f)
        f.close()

    def write_ElastoDynTower(self):

        self.fst_vt['ElastoDyn']['TwrFile'] = self.FAST_namingOut + '_ElastoDyn_tower.dat'
        tower_file = os.path.join(self.FAST_runDirectory,self.fst_vt['ElastoDyn']['TwrFile'])
        f = open(tower_file, 'w')

        f.write('------- ELASTODYN TOWER INPUT FILE -------------------------------------\n')
        f.write('Generated with OpenFAST_IO\n')
        f.write('---------------------- TOWER PARAMETERS ----------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['NTwInpSt'],  'NTwInpSt', '- Number of input stations to specify tower geometry\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['TwrFADmp1'], 'TwrFADmp(1)', '- Tower 1st fore-aft mode structural damping ratio (%)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['TwrFADmp2'], 'TwrFADmp(2)', '- Tower 2nd fore-aft mode structural damping ratio (%)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['TwrSSDmp1'], 'TwrSSDmp(1)', '- Tower 1st side-to-side mode structural damping ratio (%)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['TwrSSDmp2'], 'TwrSSDmp(2)', '- Tower 2nd side-to-side mode structural damping ratio (%)\n'))
        f.write('---------------------- TOWER ADJUSTMUNT FACTORS --------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['FAStTunr1'], 'FAStTunr(1)', '- Tower fore-aft modal stiffness tuner, 1st mode (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['FAStTunr2'], 'FAStTunr(2)', '- Tower fore-aft modal stiffness tuner, 2nd mode (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['SSStTunr1'], 'SSStTunr(1)', '- Tower side-to-side stiffness tuner, 1st mode (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['SSStTunr2'], 'SSStTunr(2)', '- Tower side-to-side stiffness tuner, 2nd mode (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['AdjTwMa'], 'AdjTwMa', '- Factor to adjust tower mass density (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['AdjFASt'], 'AdjFASt', '- Factor to adjust tower fore-aft stiffness (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['AdjSSSt'], 'AdjSSSt', '- Factor to adjust tower side-to-side stiffness (-)\n'))
        f.write('---------------------- DISTRIBUTED TOWER PROPERTIES ----------------------------\n')
        f.write('  HtFract       TMassDen         TwFAStif       TwSSStif\n')
        f.write('   (-)           (kg/m)           (Nm^2)         (Nm^2)\n')
        HtFract   = self.fst_vt['ElastoDynTower']['HtFract']
        TMassDen  = self.fst_vt['ElastoDynTower']['TMassDen']
        TwFAStif  = self.fst_vt['ElastoDynTower']['TwFAStif']
        TwSSStif  = self.fst_vt['ElastoDynTower']['TwSSStif']
        for HtFracti, TMassDeni, TwFAStifi, TwSSStifi in zip(HtFract, TMassDen, TwFAStif, TwSSStif):
            f.write('{: 2.15e} {: 2.15e} {: 2.15e} {: 2.15e}\n'.format(HtFracti, TMassDeni, TwFAStifi, TwSSStifi))
        f.write('---------------------- TOWER FORE-AFT MODE SHAPES ------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['TwFAM1Sh'][0], 'TwFAM1Sh(2)', '- Mode 1, coefficient of x^2 term\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['TwFAM1Sh'][1], 'TwFAM1Sh(3)', '-       , coefficient of x^3 term\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['TwFAM1Sh'][2], 'TwFAM1Sh(4)', '-       , coefficient of x^4 term\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['TwFAM1Sh'][3], 'TwFAM1Sh(5)', '-       , coefficient of x^5 term\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['TwFAM1Sh'][4], 'TwFAM1Sh(6)', '-       , coefficient of x^6 term\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['TwFAM2Sh'][0], 'TwFAM2Sh(2)', '- Mode 2, coefficient of x^2 term\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['TwFAM2Sh'][1], 'TwFAM2Sh(3)', '-       , coefficient of x^3 term\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['TwFAM2Sh'][2], 'TwFAM2Sh(4)', '-       , coefficient of x^4 term\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['TwFAM2Sh'][3], 'TwFAM2Sh(5)', '-       , coefficient of x^5 term\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['TwFAM2Sh'][4], 'TwFAM2Sh(6)', '-       , coefficient of x^6 term\n'))
        f.write('---------------------- TOWER SIDE-TO-SIDE MODE SHAPES --------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['TwSSM1Sh'][0], 'TwSSM1Sh(2)', '- Mode 1, coefficient of x^2 term\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['TwSSM1Sh'][1], 'TwSSM1Sh(3)', '-       , coefficient of x^3 term\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['TwSSM1Sh'][2], 'TwSSM1Sh(4)', '-       , coefficient of x^4 term\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['TwSSM1Sh'][3], 'TwSSM1Sh(5)', '-       , coefficient of x^5 term\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['TwSSM1Sh'][4], 'TwSSM1Sh(6)', '-       , coefficient of x^6 term\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['TwSSM2Sh'][0], 'TwSSM2Sh(2)', '- Mode 2, coefficient of x^2 term\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['TwSSM2Sh'][1], 'TwSSM2Sh(3)', '-       , coefficient of x^3 term\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['TwSSM2Sh'][2], 'TwSSM2Sh(4)', '-       , coefficient of x^4 term\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['TwSSM2Sh'][3], 'TwSSM2Sh(5)', '-       , coefficient of x^5 term\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ElastoDynTower']['TwSSM2Sh'][4], 'TwSSM2Sh(6)', '-       , coefficient of x^6 term\n'))
        
        f.flush()
        os.fsync(f)
        f.close()

    def write_BeamDyn(self, bldInd = None):

        # if we have bldInd is None, 
        if bldInd is None:

            self.fst_vt['BeamDyn']['BldFile'] = self.FAST_namingOut + '_BeamDyn_Blade.dat'
            bd_dict = self.fst_vt['BeamDyn']
            beamdyn_file = os.path.join(self.FAST_runDirectory,self.fst_vt['Fst']['BDBldFile(1)'])
        else:
            self.fst_vt['BeamDyn'][bldInd]['BldFile'] = self.FAST_namingOut + '_BeamDyn_Blade_%d.dat'%bldInd
            bd_dict = self.fst_vt['BeamDyn'][bldInd]
            beamdyn_file = os.path.join(self.FAST_runDirectory,self.fst_vt['Fst']['BDBldFile(%s)'%bldInd])


        self.write_BeamDynBlade(bldInd = bldInd)

        f            = open(beamdyn_file, 'w')

        f.write('--------- BEAMDYN with OpenFAST INPUT FILE -------------------------------------------\n')
        f.write('Generated with OpenFAST_IO\n')
        f.write('---------------------- SIMULATION CONTROL --------------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(bd_dict['Echo'], 'Echo', '- Echo input data to "<RootName>.ech" (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(bd_dict['QuasiStaticInit'], 'QuasiStaticInit', '- Use quasistatic pre-conditioning with centripetal accelerations in initialization (flag) [dynamic solve only]\n'))
        f.write('{:<22} {:<11} {:}'.format(bd_dict['rhoinf'], 'rhoinf', '- Numerical damping parameter for generalized-alpha integrator\n'))
        f.write('{:<22d} {:<11} {:}'.format(bd_dict['quadrature'], 'quadrature', '- Quadrature method: 1=Gaussian; 2=Trapezoidal (switch)\n'))
        f.write('{:<22} {:<11} {:}'.format(bd_dict['refine'], 'refine', '- Refinement factor for trapezoidal quadrature (-) [DEFAULT = 1; used only when quadrature=2]\n'))
        f.write('{:<22} {:<11} {:}'.format(bd_dict['n_fact'], 'n_fact', '- Factorization frequency for the Jacobian in N-R iteration(-) [DEFAULT = 5]\n'))
        f.write(float_default_out(bd_dict['DTBeam']) + '   {:<11} {:}'.format('DTBeam', '- Time step size (s).\n'))
        f.write(int_default_out(bd_dict['load_retries']) + '   {:<11} {:}'.format('load_retries', '- Number of factored load retries before quitting the aimulation [DEFAULT = 20]\n'))
        f.write(int_default_out(bd_dict['NRMax']) + '   {:<11} {:}'.format('NRMax', '- Max number of iterations in Newton-Raphson algorithm (-). [DEFAULT = 10]\n'))
        f.write(float_default_out(bd_dict['stop_tol']) + '   {:<11} {:}'.format('stop_tol', '- Tolerance for stopping criterion (-) [DEFAULT = 1E-5]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(bd_dict['tngt_stf_fd'], 'tngt_stf_fd', '- Use finite differenced tangent stiffness matrix? (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(bd_dict['tngt_stf_comp'], 'tngt_stf_comp', '- Compare analytical finite differenced tangent stiffness matrix? (flag)\n'))
        f.write(float_default_out(bd_dict['tngt_stf_pert']) + '   {:<11} {:}'.format('tngt_stf_pert', '- Perturbation size for finite differencing (-) [DEFAULT = 1E-6]\n'))
        f.write(float_default_out(bd_dict['tngt_stf_difftol']) + '   {:<11} {:}'.format('tngt_stf_difftol', '- Maximum allowable relative difference between analytical and fd tangent stiffness (-); [DEFAULT = 0.1]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(bd_dict['RotStates'], 'RotStates', '- Orient states in the rotating frame during linearization? (flag) [used only when linearizing]\n'))
        f.write('---------------------- GEOMETRY PARAMETER --------------------------------------\n')
        f.write('{:<22d} {:<11} {:}'.format(bd_dict['member_total'], 'member_total', '- Total number of members (-)\n'))
        f.write('{:<22d} {:<11} {:}'.format(bd_dict['kp_total'], 'kp_total', '- Total number of key points (-) [must be at least 3]\n'))
        for i in range(bd_dict['member_total']):
            mem = bd_dict['members'][i]
            f.write('{:<22} {:<11} {:}'.format(' '.join(['%d'%(i+1),'%d'%len(mem['kp_xr'])]), '', '- Member number; Number of key points in this member\n'))
            f.write(" ".join(['{:^21s}'.format(i) for i in ['kp_xr', 'kp_yr', 'kp_zr', 'initial_twist']])+'\n')
            f.write(" ".join(['{:^21s}'.format(i) for i in ['(m)', '(m)', '(m)', '(deg)']])+'\n')
            for j in range(len(mem['kp_xr'])):
                ln = []
                ln.append('{: 2.14e}'.format(mem['kp_xr'][j]))
                ln.append('{: 2.14e}'.format(mem['kp_yr'][j]))
                ln.append('{: 2.14e}'.format(mem['kp_zr'][j]))
                ln.append('{: 2.14e}'.format(mem['initial_twist'][j]))
                f.write(" ".join(ln) + '\n')
        f.write('---------------------- MESH PARAMETER ------------------------------------------\n')
        f.write('{:<22d} {:<11} {:}'.format(bd_dict['order_elem'], 'order_elem', '- Order of interpolation (basis) function (-)\n'))
        f.write('---------------------- MATERIAL PARAMETER --------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format('"'+bd_dict['BldFile']+'"', 'BldFile', '- Name of file containing properties for blade (quoted string)\n'))
        f.write('---------------------- PITCH ACTUATOR PARAMETERS -------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(bd_dict['UsePitchAct'], 'UsePitchAct', '- Whether a pitch actuator should be used (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format(bd_dict['PitchJ'], 'PitchJ', '- Pitch actuator inertia (kg-m^2) [used only when UsePitchAct is true]\n'))
        f.write('{:<22} {:<11} {:}'.format(bd_dict['PitchK'], 'PitchK', '- Pitch actuator stiffness (kg-m^2/s^2) [used only when UsePitchAct is true]\n'))
        f.write('{:<22} {:<11} {:}'.format(bd_dict['PitchC'], 'PitchC', '- Pitch actuator damping (kg-m^2/s) [used only when UsePitchAct is true]\n'))
        f.write('---------------------- OUTPUTS -------------------------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(bd_dict['SumPrint'], 'SumPrint', '- Print summary data to "<RootName>.sum" (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+bd_dict['OutFmt']+'"', 'OutFmt', '- Format used for text tabular output, excluding the time channel.\n'))
        f.write('{:<22} {:<11} {:}'.format(bd_dict['NNodeOuts'], 'NNodeOuts', '- Number of nodes to output to file [0 - 9] (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join(bd_dict['OutNd']), 'OutNd', '- Nodes whose values will be output  (-)\n'))
        f.write('          OutList            - The next line(s) contains a list of output parameters. See OutListParameters.xlsx for a listing of available output channels, (-)\n')
        outlist = self.get_outlist(self.fst_vt['outlist'], ['BeamDyn'])
        for channel_list in outlist:
            for i in range(len(channel_list)):
                f.write('"' + channel_list[i] + '"\n')
        f.write('END of input file (the word "END" must appear in the first 3 columns of the last OutList line)\n')
                
        # Optional nodal output section
        if 'BldNd_BlOutNd' in bd_dict:
            f.write('====== Outputs for all blade stations (same ending as above for B1N1.... =========================== [optional section]\n')
            # f.write('{:<22d} {:<11} {:}'.format(bd_dict['BldNd_BladesOut'], 'BldNd_BladesOut', '- Number of blades to output all node information at.  Up to number of blades on turbine. (-)\n'))
            f.write('{!s:<22} {:<11} {:}'.format(bd_dict['BldNd_BlOutNd'], 'BldNd_BlOutNd', '- Future feature will allow selecting a portion of the nodes to output.  Not implemented yet. (-)\n'))
            f.write('                   OutList     - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx, BeamDyn_Nodes tab for a listing of available output channels, (-)\n')
            
            opt_outlist = self.get_outlist(self.fst_vt['outlist'], ['BeamDyn_Nodes'])      
            for opt_channel_list in opt_outlist:
                for i in range(len(opt_channel_list)):
                    f.write('"' + opt_channel_list[i] + '"\n')
            f.write('END of input file (the word "END" must appear in the first 3 columns of the last OutList line)\n')
        
        f.write('---------------------------------------------------------------------------------------')
        f.flush()
        os.fsync(f)
        f.close()

    def write_BeamDynBlade(self, bldInd = None):

        if bldInd is None:
            bd_blade_dict = self.fst_vt['BeamDynBlade']
            bd_blade_file = os.path.abspath(os.path.join(self.FAST_runDirectory, self.fst_vt['BeamDyn']['BldFile']))
        else:
            bd_blade_dict = self.fst_vt['BeamDynBlade'][bldInd]
            bd_blade_file = os.path.abspath(os.path.join(self.FAST_runDirectory, self.fst_vt['BeamDyn'][bldInd]['BldFile']))

        f = open(bd_blade_file, 'w')

        f.write('------- BEAMDYN INDIVIDUAL BLADE INPUT FILE --------------------------\n')
        f.write('Generated with OpenFAST_IO\n')
        f.write('---------------------- BLADE PARAMETERS --------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(bd_blade_dict['station_total'], 'station_total', '- Number of blade input stations (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(bd_blade_dict['damp_type'], 'damp_type', '- Damping type: 0: no damping; 1: damped\n'))
        f.write('---------------------- DAMPING COEFFICIENT------------------------------------\n')
        f.write(" ".join(['{:^11s}'.format(i) for i in ['mu1','mu2','mu3','mu4','mu5','mu6']])+'\n')
        f.write(" ".join(['{:^11s}'.format(i) for i in ['(-)','(-)','(-)','(-)','(-)','(-)']])+'\n')
        mu = [bd_blade_dict['mu1'], bd_blade_dict['mu2'], bd_blade_dict['mu3'], bd_blade_dict['mu4'], bd_blade_dict['mu5'], bd_blade_dict['mu6']]
        f.write(" ".join(['{:^11f}'.format(i) for i in mu])+'\n')
        f.write('---------------------- DISTRIBUTED PROPERTIES---------------------------------\n')

        for i in range(len(bd_blade_dict['radial_stations'])):
            f.write('{: 2.15e}\n'.format(bd_blade_dict['radial_stations'][i]))
            for j in range(6):
                f.write(" ".join(['{: 2.15e}'.format(i) for i in bd_blade_dict['beam_stiff'][i,j,:]])+'\n')
            f.write('\n')
            for j in range(6):
                f.write(" ".join(['{: 2.15e}'.format(i) for i in bd_blade_dict['beam_inertia'][i,j,:]])+'\n')
            f.write('\n')

        f.write('\n')

    def write_InflowWind(self):
        self.fst_vt['Fst']['InflowFile'] = self.FAST_namingOut + '_InflowWind.dat'
        inflow_file = os.path.join(self.FAST_runDirectory,self.fst_vt['Fst']['InflowFile'])
        f = open(inflow_file, 'w')

        f.write('------- InflowWind INPUT FILE -------------------------------------------------------------------------\n')
        f.write('Generated with OpenFAST_IO\n')
        f.write('---------------------------------------------------------------------------------------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['Echo'], 'Echo', '- Echo input data to <RootName>.ech (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['WindType'], 'WindType', '- switch for wind file type (1=steady; 2=uniform; 3=binary TurbSim FF; 4=binary Bladed-style FF; 5=HAWC format; 6=User defined; 7=native Bladed FF)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['PropagationDir'], 'PropagationDir', '- Direction of wind propagation (meteoroligical rotation from aligned with X (positive rotates towards -Y) -- degrees) (not used for native Bladed format WindType=7)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['VFlowAng'], 'VFlowAng', '- Upflow angle (degrees) (not used for native Bladed format WindType=7)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['VelInterpCubic'], 'VelInterpCubic', '- Use cubic interpolation for velocity in time (false=linear, true=cubic) [Used with WindType=2,3,4,5,7]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['NWindVel'], 'NWindVel', '- Number of points to output the wind velocity    (0 to 9)\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join(np.array(self.fst_vt['InflowWind']['WindVxiList'], dtype=str)), 'WindVxiList', '- List of coordinates in the inertial X direction (m)\n'))        
        f.write('{:<22} {:<11} {:}'.format(', '.join(np.array(self.fst_vt['InflowWind']['WindVyiList'], dtype=str)), 'WindVyiList', '- List of coordinates in the inertial Y direction (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join(np.array(self.fst_vt['InflowWind']['WindVziList'], dtype=str)), 'WindVziList', '- List of coordinates in the inertial Z direction (m)\n'))
        f.write('================== Parameters for Steady Wind Conditions [used only for WindType = 1] =========================\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['HWindSpeed'], 'HWindSpeed', '- Horizontal wind speed                            (m/s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['RefHt'], 'RefHt', '- Reference height for horizontal wind speed      (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['PLExp'], 'PLExp', '- Power law exponent                              (-)\n'))
        f.write('================== Parameters for Uniform wind file   [used only for WindType = 2] ============================\n')
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['InflowWind']['FileName_Uni']+'"', 'FileName_Uni', '- Filename of time series data for uniform wind field.      (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['RefHt_Uni'], 'RefHt_Uni', '- Reference height for horizontal wind speed                (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['RefLength'], 'RefLength', '- Reference length for linear horizontal and vertical sheer (-)\n'))
        f.write('================== Parameters for Binary TurbSim Full-Field files   [used only for WindType = 3] ==============\n')
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['InflowWind']['FileName_BTS']+'"', 'FileName_BTS', '- Name of the Full field wind file to use (.bts)\n'))
        f.write('================== Parameters for Binary Bladed-style Full-Field files   [used only for WindType = 4 or WindType = 7] =========\n')
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['InflowWind']['FileNameRoot']+'"', 'FileNameRoot', '- WindType=4: Rootname of the full-field wind file to use (.wnd, .sum); WindType=7: name of the intermediate file with wind scaling values\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['TowerFile'], 'TowerFile', '- Have tower file (.twr) (flag) ignored when WindType = 7\n'))
        f.write('================== Parameters for HAWC-format binary files  [Only used with WindType = 5] =====================\n')
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['InflowWind']['FileName_u']+'"', 'FileName_u', '- name of the file containing the u-component fluctuating wind (.bin)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['InflowWind']['FileName_v']+'"', 'FileName_v', '- name of the file containing the v-component fluctuating wind (.bin)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['InflowWind']['FileName_w']+'"', 'FileName_w', '- name of the file containing the w-component fluctuating wind (.bin)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['nx'], 'nx', '- number of grids in the x direction (in the 3 files above) (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['ny'], 'ny', '- number of grids in the y direction (in the 3 files above) (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['nz'], 'nz', '- number of grids in the z direction (in the 3 files above) (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['dx'], 'dx', '- distance (in meters) between points in the x direction    (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['dy'], 'dy', '- distance (in meters) between points in the y direction    (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['dz'], 'dz', '- distance (in meters) between points in the z direction    (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['RefHt_Hawc'], 'RefHt_Hawc', '- reference height; the height (in meters) of the vertical center of the grid (m)\n'))
        f.write('-------------   Scaling parameters for turbulence   ---------------------------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['ScaleMethod'], 'ScaleMethod', '- Turbulence scaling method   [0 = none, 1 = direct scaling, 2 = calculate scaling factor based on a desired standard deviation]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['SFx'], 'SFx', '- Turbulence scaling factor for the x direction (-)   [ScaleMethod=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['SFy'], 'SFy', '- Turbulence scaling factor for the y direction (-)   [ScaleMethod=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['SFz'], 'SFz', '- Turbulence scaling factor for the z direction (-)   [ScaleMethod=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['SigmaFx'], 'SigmaFx', '- Turbulence standard deviation to calculate scaling from in x direction (m/s)    [ScaleMethod=2]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['SigmaFy'], 'SigmaFy', '- Turbulence standard deviation to calculate scaling from in y direction (m/s)    [ScaleMethod=2]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['SigmaFz'], 'SigmaFz', '- Turbulence standard deviation to calculate scaling from in z direction (m/s)    [ScaleMethod=2]\n'))
        f.write('-------------   Mean wind profile parameters (added to HAWC-format files)   ---------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['URef'], 'URef', '- Mean u-component wind speed at the reference height (m/s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['WindProfile'], 'WindProfile', '- Wind profile type (0=constant;1=logarithmic,2=power law)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['PLExp_Hawc'], 'PLExp_Hawc', '- Power law exponent (-) (used for PL wind profile type only)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['Z0'], 'Z0', '- Surface roughness length (m) (used for LG wind profile type only)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['XOffset'], 'XOffset', '- Initial offset in +x direction (shift of wind box) (-)\n'))
        f.write('-------------   LIDAR Parameters   --------------------------------------------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['SensorType'], 'SensorType', '- Switch for lidar configuration (0 = None, 1 = Single Point Beam(s), 2 = Continuous, 3 = Pulsed)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['NumPulseGate'], 'NumPulseGate', '- Number of lidar measurement gates (used when SensorType = 3)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['PulseSpacing'], 'PulseSpacing', '- Distance between range gates (m) (used when SensorType = 3)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['NumBeam'], 'NumBeam', '- Number of lidar measurement beams (0-5)(used when SensorType = 1)\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join(np.array(self.fst_vt['InflowWind']['FocalDistanceX'], dtype=str)), 'FocalDistanceX', '- Focal distance co-ordinates of the lidar beam in the x direction (relative to hub height) (only first coordinate used for SensorType 2 and 3) (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join(np.array(self.fst_vt['InflowWind']['FocalDistanceY'], dtype=str)), 'FocalDistanceY', '- Focal distance co-ordinates of the lidar beam in the y direction (relative to hub height) (only first coordinate used for SensorType 2 and 3) (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join(np.array(self.fst_vt['InflowWind']['FocalDistanceZ'], dtype=str)), 'FocalDistanceZ', '- Focal distance co-ordinates of the lidar beam in the z direction (relative to hub height) (only first coordinate used for SensorType 2 and 3) (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join(np.array(self.fst_vt['InflowWind']['RotorApexOffsetPos'], dtype=str)), 'RotorApexOffsetPos', '- Offset of the lidar from hub height (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['URefLid'], 'URefLid', '- Reference average wind speed for the lidar[m/s]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['MeasurementInterval'], 'MeasurementInterval', '- Time between each measurement [s]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['LidRadialVel'], 'LidRadialVel', '- TRUE => return radial component, FALSE => return \'x\' direction estimate\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['ConsiderHubMotion'], 'ConsiderHubMotion', '- Flag whether to consider the hub motion\'s impact on Lidar measurements\n'))
        f.write('====================== OUTPUT ==================================================\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['InflowWind']['SumPrint'], 'SumPrint', '- Print summary data to <RootName>.IfW.sum (flag)\n'))
        f.write('OutList      - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)\n')
        
        outlist = self.get_outlist(self.fst_vt['outlist'], ['InflowWind'])
        for channel_list in outlist:
            for i in range(len(channel_list)):
                f.write('"' + channel_list[i] + '"\n')

        f.write('END of input file (the word "END" must appear in the first 3 columns of the last OutList line)\n')
        f.write('---------------------------------------------------------------------------------------\n')

        f.flush()
        os.fsync(f)
        f.close()
  
    def write_AeroDyn(self):
        # AeroDyn v15.03

        # Generate AeroDyn v15 blade input file
        if isinstance(self.fst_vt['AeroDynBlade'], list):
            for i_adBld, adBld in enumerate(self.fst_vt['AeroDynBlade']):
                self.fst_vt['AeroDyn']['ADBlFile%d'%(i_adBld+1)] = self.FAST_namingOut + '_AeroDyn_blade_%d.dat'%(i_adBld+1)
                self.write_AeroDynBlade(bldInd = i_adBld)

        elif isinstance(self.fst_vt['AeroDynBlade'], dict):
            self.fst_vt['AeroDyn']['ADBlFile1'] = self.FAST_namingOut + '_AeroDyn_blade.dat'
            self.fst_vt['AeroDyn']['ADBlFile2'] = self.fst_vt['AeroDyn']['ADBlFile1']
            self.fst_vt['AeroDyn']['ADBlFile3'] = self.fst_vt['AeroDyn']['ADBlFile1']
            self.write_AeroDynBlade()


        # self.write_AeroDynBlade()

        # Generate AeroDyn v15 polars
        self.write_AeroDynPolar()
        
        # Generate AeroDyn v15 airfoil coordinates
        # some polars may have airfoil coordinates, need to account for all possible scenarios
        if any([self.fst_vt['AeroDyn']['af_data'][i][0]['NumCoords'] != '0' for i in range(len(self.fst_vt['AeroDyn']['af_data']))]):
            af_coords = [i for i in range(len(self.fst_vt['AeroDyn']['af_data'])) if self.fst_vt['AeroDyn']['af_data'][i][0]['NumCoords'] != '0']
            self.write_AeroDynCoord(af_coords)
        
        if self.fst_vt['AeroDyn']['Wake_Mod'] == 3:
            if self.fst_vt['AeroDyn']['UA_Mod'] > 0:
                raise Exception('OLAF is called with unsteady airfoil aerodynamics, but OLAF currently only supports UA_Mod == 0') #TODO: need to check if this holds true now
            self.write_OLAF()

        # Generate AeroDyn v15.03 input file
        self.fst_vt['Fst']['AeroFile'] = self.FAST_namingOut + '_AeroDyn.dat'
        ad_file = os.path.join(self.FAST_runDirectory, self.fst_vt['Fst']['AeroFile'])
        f = open(ad_file, 'w')

        f.write('------- AERODYN15 INPUT FILE ------------------------------------------------\n')
        f.write('Generated with OpenFAST_IO\n')
        f.write('======  General Options  ============================================================================\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['Echo'], 'Echo', '- Echo the input to "<rootname>.AD.ech"?  (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['DTAero'], 'DTAero', '- Time interval for aerodynamic calculations {or "default"} (s)\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['AeroDyn']['Wake_Mod'], 'Wake_Mod', '- Wake/induction model (switch) {0=none, 1=BEMT, 3=OLAF} [Wake_Mod cannot be 2 or 3 when linearizing]\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['AeroDyn']['TwrPotent'], 'TwrPotent', '- Type tower influence on wind based on potential flow around the tower (switch) {0=none, 1=baseline potential flow, 2=potential flow with Bak correction}\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['AeroDyn']['TwrShadow'], 'TwrShadow', '- Calculate tower influence on wind based on downstream tower shadow (switch) {0=none, 1=Powles model, 2=Eames model}\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['TwrAero'], 'TwrAero', '- Calculate tower aerodynamic loads? (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['CavitCheck'], 'CavitCheck', '- Perform cavitation check? (flag) [UA_Mod must be 0 when CavitCheck=true]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['Buoyancy'], 'Buoyancy', '- Include buoyancy effects? (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['NacelleDrag'], 'NacelleDrag', '- Include Nacelle Drag effects? (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['CompAA'], 'CompAA', '- Flag to compute AeroAcoustics calculation [used only when Wake_Mod = 1 or 2]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['AA_InputFile'], 'AA_InputFile', '- AeroAcoustics input file [used only when CompAA=true]\n'))
        f.write('======  Environmental Conditions  ===================================================================\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['AirDens'], 'AirDens', '- Air density (kg/m^3)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['KinVisc'], 'KinVisc', '- Kinematic viscosity of working fluid (m^2/s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['SpdSound'], 'SpdSound', '- Speed of sound in working fluid (m/s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['Patm'], 'Patm', '- Atmospheric pressure (Pa) [used only when CavitCheck=True]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['Pvap'], 'Pvap', '- Vapour pressure of working fluid (Pa) [used only when CavitCheck=True]\n'))
        f.write('======  Blade-Element/Momentum Theory Options  ====================================================== [unused when Wake_Mod=0 or 3, except for BEM_Mod]\n')
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['AeroDyn']['BEM_Mod'], 'BEM_Mod', '- BEM model {1=legacy NoSweepPitchTwist, 2=polar} (switch) [used for all Wake_Mod to determine output coordinate system]\n'))
        f.write('--- Skew correction\n')
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['AeroDyn']['Skew_Mod'], 'Skew_Mod', '- Skew model {0=No skew model, -1=Remove non-normal component for linearization, 1=skew model active}\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['SkewMomCorr'], 'SkewMomCorr', '- Turn the skew momentum correction on or off [used only when Skew_Mod=1]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['SkewRedistr_Mod'], 'SkewRedistr_Mod', '- Type of skewed-wake correction model (switch) {0=no redistribution, 1=Glauert/Pitt/Peters, default=1} [used only when Skew_Mod=1]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['SkewRedistrFactor'], 'SkewRedistrFactor', '- Constant used in Pitt/Peters skewed wake model {or "default" is 15/32*pi} (-) [used only when Skew_Mod=1 and SkewRedistr_Mod=1]\n'))
        f.write('--- BEM algorithm\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['TipLoss'], 'TipLoss', '- Use the Prandtl tip-loss model? (flag) [unused when Wake_Mod=0 or 3]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['HubLoss'], 'HubLoss', '- Use the Prandtl hub-loss model? (flag) [unused when Wake_Mod=0 or 3]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['TanInd'], 'TanInd', '- Include tangential induction in BEMT calculations? (flag) [unused when Wake_Mod=0 or 3]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['AIDrag'], 'AIDrag', '- Include the drag term in the axial-induction calculation? (flag) [unused when Wake_Mod=0 or 3]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['TIDrag'], 'TIDrag', '- Include the drag term in the tangential-induction calculation? (flag) [unused when Wake_Mod=0,3 or TanInd=FALSE]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['IndToler'], 'IndToler', '- Convergence tolerance for BEMT nonlinear solve residual equation {or "default"} (-) [unused when Wake_Mod=0 or 3]\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['AeroDyn']['MaxIter'], 'MaxIter', '- Maximum number of iteration steps (-) [unused when Wake_Mod=0]\n'))
        f.write('--- Shear correction\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['SectAvg'], 'SectAvg', '- Use sector averaging (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['SectAvgWeighting'], 'SectAvgWeighting', '- Weighting function for sector average {1=Uniform, default=1} within a sector centered on the blade (switch) [used only when SectAvg=True]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['SectAvgNPoints'], 'SectAvgNPoints', '- Number of points per sectors (-) {default=5} [used only when SectAvg=True]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['SectAvgPsiBwd'], 'SectAvgPsiBwd', '- Backward azimuth relative to blade where the sector starts (<=0) {default=-60} (deg) [used only when SectAvg=True]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['SectAvgPsiFwd'], 'SectAvgPsiFwd', '- Forward azimuth relative to blade where the sector ends (>=0) {default=60} (deg) [used only when SectAvg=True]\n'))

        f.write('--- Dynamic wake/inflow\n')
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['AeroDyn']['DBEMT_Mod'], 'DBEMT_Mod', '- Type of dynamic BEMT (DBEMT) model {0=No Dynamic Wake, -1=Frozen Wake for linearization, 1:constant tau1, 2=time-dependent tau1, 3=constant tau1 with continuous formulation} (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['tau1_const'], 'tau1_const', '- Time constant for DBEMT (s) [used only when DBEMT_Mod=1 or 3]\n'))
        f.write('======  OLAF -- cOnvecting LAgrangian Filaments (Free Vortex Wake) Theory Options  ================== [used only when Wake_Mod=3]\n')
        olaf_file = self.FAST_namingOut + '_OLAF.dat'
        f.write('{!s:<22} {:<11} {:}'.format(olaf_file, 'OLAFInputFileName', '- Input file for OLAF [used only when Wake_Mod=3]\n'))  
        f.write('======  Unsteady Airfoil Aerodynamics Options  ===================================== \n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['AoA34'], 'AoA34', "- Sample the angle of attack (AoA) at the 3/4 chord or the AC point {default=True} [always used]\n"))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['AeroDyn']['UA_Mod'], 'UA_Mod', "- Unsteady Aero Model Switch (switch) {0=Quasi-steady (no UA), 2=B-L Gonzalez, 3=B-L Minnema/Pierce, 4=B-L HGM 4-states, 5=B-L HGM+vortex 5 states, 6=Oye, 7=Boeing-Vertol}\n"))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['FLookup'], 'FLookup', "- Flag to indicate whether a lookup for f' will be calculated (TRUE) or whether best-fit exponential equations will be used (FALSE); if FALSE S1-S4 must be provided in airfoil input files (flag) [used only when UA_Mod>0]\n"))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['IntegrationMethod'], 'IntegrationMethod', "- Switch to indicate which integration method UA uses (1=RK4, 2=AB4, 3=ABM4, 4=BDF2)\n"))
        if 'UAStartRad' in self.fst_vt['AeroDyn'] and 'UAEndRad' in self.fst_vt['AeroDyn']:
            f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['UAStartRad'], 'UAStartRad', '- Starting radius for dynamic stall (fraction of rotor radius [0.0,1.0]) [used only when UA_Mod>0; if line is missing UAStartRad=0]\n'))
            f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['UAEndRad'], 'UAEndRad', '- Ending radius for dynamic stall (fraction of rotor radius [0.0,1.0]) [used only when UA_Mod>0; if line is missing UAEndRad=1]\n'))
        f.write('======  Airfoil Information =========================================================================\n')
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['AeroDyn']['AFTabMod'], 'AFTabMod', '- Interpolation method for multiple airfoil tables {1=1D interpolation on AoA (first table only); 2=2D interpolation on AoA and Re; 3=2D interpolation on AoA and UserProp} (-)\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['AeroDyn']['InCol_Alfa'], 'InCol_Alfa', '- The column in the airfoil tables that contains the angle of attack (-)\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['AeroDyn']['InCol_Cl'], 'InCol_Cl', '- The column in the airfoil tables that contains the lift coefficient (-)\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['AeroDyn']['InCol_Cd'], 'InCol_Cd', '- The column in the airfoil tables that contains the drag coefficient (-)\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['AeroDyn']['InCol_Cm'], 'InCol_Cm', '- The column in the airfoil tables that contains the pitching-moment coefficient; use zero if there is no Cm column (-)\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['AeroDyn']['InCol_Cpmin'], 'InCol_Cpmin', '- The column in the airfoil tables that contains the Cpmin coefficient; use zero if there is no Cpmin column (-)\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['AeroDyn']['NumAFfiles'], 'NumAFfiles', '- Number of airfoil files used (-)\n'))
        for i in range(self.fst_vt['AeroDyn']['NumAFfiles']):
            if i == 0:
                f.write('"' + self.fst_vt['AeroDyn']['AFNames'][i] + '"    AFNames            - Airfoil file names (NumAFfiles lines) (quoted strings)\n')
            else:
                f.write('"' + self.fst_vt['AeroDyn']['AFNames'][i] + '"\n')
        f.write('======  Rotor/Blade Properties  =====================================================================\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['UseBlCm'], 'UseBlCm', '- Include aerodynamic pitching moment in calculations?  (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['AeroDyn']['ADBlFile1']+'"', 'ADBlFile(1)', '- Name of file containing distributed aerodynamic properties for Blade #1 (-)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['AeroDyn']['ADBlFile2']+'"', 'ADBlFile(2)', '- Name of file containing distributed aerodynamic properties for Blade #2 (-) [unused if NumBl < 2]\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['AeroDyn']['ADBlFile3']+'"', 'ADBlFile(3)', '- Name of file containing distributed aerodynamic properties for Blade #3 (-) [unused if NumBl < 3]\n'))
        f.write('======  Hub Properties ============================================================================== [used only when Buoyancy=True]\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['VolHub'], 'VolHub', '- Hub volume (m^3)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['HubCenBx'], 'HubCenBx', '- Hub center of buoyancy x direction offset (m)\n'))
        f.write('======  Nacelle Properties ========================================================================== [used only when Buoyancy=True or NacelleDrag=True]\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['VolNac'], 'VolNac', '- Nacelle volume (m^3)\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join(np.array(self.fst_vt['AeroDyn']['NacCenB'], dtype=str)), 'NacCenB', '- Position of nacelle center of buoyancy from yaw bearing in nacelle coordinates (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join(np.array(self.fst_vt['AeroDyn']['NacArea'], dtype=str)), 'NacArea', '- Projected area of the nacelle in X, Y, Z in the nacelle coordinate system (m^2)\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join(np.array(self.fst_vt['AeroDyn']['NacCd'], dtype=str)), 'NacCd', '- Drag coefficient for the nacelle areas defined above (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join(np.array(self.fst_vt['AeroDyn']['NacDragAC'], dtype=str)), 'NacDragAC', '- Position of aerodynamic center of nacelle drag in nacelle coordinates (m)\n'))
        f.write('======  Tail Fin Aerodynamics ========================================================================\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['TFinAero'], 'TFinAero', '- Calculate tail fin aerodynamics model (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['AeroDyn']['TFinFile']+'"', 'TFinFile', '- Input file for tail fin aerodynamics [used only when TFinAero=True]\n'))
        f.write('======  Tower Influence and Aerodynamics ============================================================ [used only when TwrPotent/=0, TwrShadow/=0, TwrAero=True, or Buoyancy=True]\n')
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['AeroDyn']['NumTwrNds'], 'NumTwrNds', '- Number of tower nodes used in the analysis  (-) [used only when TwrPotent/=0, TwrShadow/=0, TwrAero=True, or Buoyancy=True]\n'))
        f.write('TwrElev        TwrDiam        TwrCd          TwrTI          TwrCb !TwrTI used only when TwrShadow=2; TwrCb used only when Buoyancy=True\n')
        f.write('(m)              (m)           (-)            (-)            (-)\n')
        for TwrElev, TwrDiam, TwrCd, TwrTI, TwrCb in zip(self.fst_vt['AeroDyn']['TwrElev'], self.fst_vt['AeroDyn']['TwrDiam'], self.fst_vt['AeroDyn']['TwrCd'], self.fst_vt['AeroDyn']['TwrTI'], self.fst_vt['AeroDyn']['TwrCb']):
            f.write('{: 2.15e} {: 2.15e} {: 2.15e} {: 2.15e} {: 2.15e} \n'.format(TwrElev, TwrDiam, TwrCd, TwrTI, TwrCb))
        f.write('======  Outputs  ====================================================================================\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['SumPrint'], 'SumPrint', '- Generate a summary file listing input options and interpolated properties to "<rootname>.AD.sum"?  (flag)\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['AeroDyn']['NBlOuts'], 'NBlOuts', '- Number of blade node outputs [0 - 9] (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join(self.fst_vt['AeroDyn']['BlOutNd']), 'BlOutNd', '- Blade nodes whose values will be output  (-)\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['AeroDyn']['NTwOuts'], 'NTwOuts', '- Number of tower node outputs [0 - 9]  (-)\n'))
        # if self.fst_vt['AeroDyn']['NTwOuts'] != 0: # TODO its weird that tower nodes is treated differently than blade nodes
        #     f.write('{:<22} {:<11} {:}'.format(', '.join(self.fst_vt['AeroDyn']['TwOutNd']), 'TwOutNd', '- Tower nodes whose values will be output  (-)\n'))
        # else:
        #     f.write('{:<22} {:<11} {:}'.format(0, 'TwOutNd', '- Tower nodes whose values will be output  (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join(np.array(self.fst_vt['AeroDyn']['TwOutNd'], dtype=str)), 'TwOutNd', '- Tower nodes whose values will be output  (-)\n'))
        f.write('                   OutList             - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)\n')

        outlist = self.get_outlist(self.fst_vt['outlist'], ['AeroDyn'])      
        for channel_list in outlist:
            for i in range(len(channel_list)):
                f.write('"' + channel_list[i] + '"\n')
        f.write('END of input file (the word "END" must appear in the first 3 columns of the last OutList line)\n')
        
        # Optional nodal output section
        if 'BldNd_BladesOut' in self.fst_vt['AeroDyn']:
            f.write('====== Outputs for all blade stations (same ending as above for B1N1.... =========================== [optional section]\n')
            f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['AeroDyn']['BldNd_BladesOut'], 'BldNd_BladesOut', '- Number of blades to output all node information at.  Up to number of blades on turbine. (-)\n'))
            f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['BldNd_BlOutNd'], 'BldNd_BlOutNd', '- Future feature will allow selecting a portion of the nodes to output.  Not implemented yet. (-)\n'))
            f.write('                   OutList_Nodal     - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx, AeroDyn_Nodes tab for a listing of available output channels, (-)\n')
            
            opt_outlist = self.get_outlist(self.fst_vt['outlist'], ['AeroDyn_Nodes'])      
            for opt_channel_list in opt_outlist:
                for i in range(len(opt_channel_list)):
                    f.write('"' + opt_channel_list[i] + '"\n')
            f.write('END of input file (the word "END" must appear in the first 3 columns of the last OutList line)\n')
        
        f.write('---------------------------------------------------------------------------------------\n')
        f.flush()
        os.fsync(f)
        f.close()

    def write_AeroDynBlade(self, bldInd = None):

        if bldInd is None:
            filename = os.path.join(self.FAST_runDirectory, self.fst_vt['AeroDyn']['ADBlFile1'])
            adBld_dict = self.fst_vt['AeroDynBlade']
        else:
            filename = os.path.join(self.FAST_runDirectory, self.fst_vt['AeroDyn']['ADBlFile%d'%(bldInd+1)])
            adBld_dict = self.fst_vt['AeroDynBlade'][bldInd]
        
        f = open(filename, 'w')

        f.write('------- AERODYN15 BLADE DEFINITION INPUT FILE -------------------------------------\n')
        f.write('Generated with OpenFAST_IO\n')
        f.write('======  Blade Properties =================================================================\n')
        f.write('{:<11d} {:<11} {:}'.format(adBld_dict['NumBlNds'], 'NumBlNds', '- Number of blade nodes used in the analysis (-)\n'))
        f.write('    BlSpn        BlCrvAC        BlSwpAC        BlCrvAng       BlTwist        BlChord          BlAFID       BlCb        BlCenBn      BlCenBt\n')
        f.write('     (m)           (m)            (m)            (deg)         (deg)           (m)              (-)        (-)         (m)             (m)\n')
        BlSpn    = adBld_dict['BlSpn']
        BlCrvAC  = adBld_dict['BlCrvAC']
        BlSwpAC  = adBld_dict['BlSwpAC']
        BlCrvAng = adBld_dict['BlCrvAng']
        BlTwist  = adBld_dict['BlTwist']
        BlChord  = adBld_dict['BlChord']
        BlAFID   = adBld_dict['BlAFID']
        BlCb     = adBld_dict['BlCb']
        BlCenBn  = adBld_dict['BlCenBn']
        BlCenBt  = adBld_dict['BlCenBt']
        for Spn, CrvAC, SwpAC, CrvAng, Twist, Chord, AFID, BlCb, BlCenBn, BlCenBt in zip(BlSpn, BlCrvAC, BlSwpAC, BlCrvAng, BlTwist, BlChord, BlAFID, BlCb, BlCenBn, BlCenBt):
            f.write('{: 2.15e} {: 2.15e} {: 2.15e} {: 2.15e} {: 2.15e} {: 2.15e} {: 8d} {: 2.15e} {: 2.15e} {: 2.15e}\n'.format(Spn, CrvAC, SwpAC, CrvAng, Twist, Chord, int(AFID), BlCb, BlCenBn, BlCenBt))
        
        f.flush()
        os.fsync(f)
        f.close()
        
    def write_AeroDynPolar(self):
        # Airfoil Info v1.01

        if not os.path.isdir(os.path.join(self.FAST_runDirectory,'Airfoils')):
            try:
                os.makedirs(os.path.join(self.FAST_runDirectory,'Airfoils'))
            except:
                try:
                    time.sleep(random.random())
                    if not os.path.isdir(os.path.join(self.FAST_runDirectory,'Airfoils')):
                        os.makedirs(os.path.join(self.FAST_runDirectory,'Airfoils'))
                except:
                    print("Error tring to make '%s'!"%os.path.join(self.FAST_runDirectory,'Airfoils'))


        self.fst_vt['AeroDyn']['NumAFfiles'] = len(self.fst_vt['AeroDyn']['af_data'])
        self.fst_vt['AeroDyn']['AFNames'] = ['']*self.fst_vt['AeroDyn']['NumAFfiles']

        for afi in range(int(self.fst_vt['AeroDyn']['NumAFfiles'])):

            self.fst_vt['AeroDyn']['AFNames'][afi] = os.path.join('Airfoils', self.FAST_namingOut + '_AeroDyn_Polar_%02d.dat'%afi)
            af_file = os.path.join(self.FAST_runDirectory, self.fst_vt['AeroDyn']['AFNames'][afi])
            f = open(af_file, 'w')

            f.write('! ------------ AirfoilInfo Input File ----------------------------------\n')
            f.write('! Generated with OpenFAST_IO\n')
            f.write('! line\n')
            f.write('! line\n')
            f.write('! ------------------------------------------------------------------------------\n')
            f.write('{:<22}   {:<11} {:}'.format(self.fst_vt['AeroDyn']['af_data'][afi][0]['InterpOrd'], 'InterpOrd', '! Interpolation order to use for quasi-steady table lookup {1=linear; 3=cubic spline; "default"} [default=3]\n'))
            if 'RelThickness' in self.fst_vt['AeroDyn']['af_data'][afi][0]:
                f.write('{:<22f}   {:<11} {:}'.format(self.fst_vt['AeroDyn']['af_data'][afi][0]['RelThickness'], 'RelThickness', '! The non-dimensional thickness of the airfoil (thickness/chord) [only used if UAMod=7] [default=0.2] (-)\n'))
            
            f.write('{:<22}   {:<11} {:}'.format(self.fst_vt['AeroDyn']['af_data'][afi][0]['NonDimArea'], 'NonDimArea', '! The non-dimensional area of the airfoil (area/chord^2) (set to 1.0 if unsure or unneeded)\n'))
            if self.fst_vt['AeroDyn']['af_data'][afi][0]['NumCoords'] != '0':
                f.write('@"{:}_AF{:02d}_Coords.txt"       {:<11} {:}'.format(self.FAST_namingOut, afi, 'NumCoords', '! The number of coordinates in the airfoil shape file. Set to zero if coordinates not included.\n'))
            else:
                f.write('{:<22d}       {:<11} {:}'.format(0, 'NumCoords', '! The number of coordinates in the airfoil shape file. Set to zero if coordinates not included.\n'))
            f.write('AF{:02d}_BL.txt              {:<11} {:}'.format(afi, 'BL_file', '! The file name including the boundary layer characteristics of the profile. Ignored if the aeroacoustic module is not called.\n'))
            # f.write('{:<22d}   {:<11} {:}'.format(self.fst_vt['AeroDyn']['af_data'][afi][0]['NumTabs'], 'NumTabs', '! Number of airfoil tables in this file.  Each table must have lines for Re and UserProp.\n'))


            # Check if we have multiple tables per airfoil
            # if yes, allocate the number of airfoils to the respective radial stations
            if self.fst_vt['AeroDyn']['AFTabMod'] == 2:
                num_tab = len(self.fst_vt['AeroDyn']['af_data'][afi])
            elif self.fst_vt['AeroDyn']['AFTabMod'] == 3:
                # for tab_orig in range(self.fst_vt['AeroDyn']['af_data'][afi][0]['NumTabs'] - 1):
                if len( self.fst_vt['AeroDyn']['af_data'][afi]) == 1 or \
                    self.fst_vt['AeroDyn']['af_data'][afi][0]['UserProp'] == self.fst_vt['AeroDyn']['af_data'][afi][1]['UserProp']:
                    num_tab = 1  # assume that all UserProp angles of the flaps are identical if the first two are -> no flaps!
                else:
                    num_tab = self.fst_vt['AeroDyn']['af_data'][afi][0]['NumTabs']
            else:
                num_tab = 1
            # f.write('{:<22d}   {:<11} {:}'.format(self.fst_vt['AeroDyn']['af_data'][afi][0]['NumTabs'], 'NumTabs','! Number of airfoil tables in this file.  Each table must have lines for Re and UserProp.\n'))
            f.write('{:<22d}   {:<11} {:}'.format(num_tab, 'NumTabs','! Number of airfoil tables in this file.  Each table must have lines for Re and UserProp.\n'))

            # for tab in range(self.fst_vt['AeroDyn']['af_data'][afi][0]['NumTabs']): # For writing multiple tables (different Re or UserProp values)
            for tab in range(num_tab): # For writing multiple tables (different Re or UserProp values)
                f.write('! ------------------------------------------------------------------------------\n')
                f.write("! data for table %i \n" % (tab + 1))
                f.write('! ------------------------------------------------------------------------------\n')
                f.write('{:<22f}   {:<11} {:}'.format(self.fst_vt['AeroDyn']['af_data'][afi][tab]['Re']*1.e-6, 'Re', '! Reynolds number in millions\n'))
                f.write('{:<22d}   {:<11} {:}'.format(int(self.fst_vt['AeroDyn']['af_data'][afi][tab]['UserProp']), 'UserProp', '! User property (control) setting\n'))
                f.write('{!s:<22}   {:<11} {:}'.format(self.fst_vt['AeroDyn']['af_data'][afi][tab]['InclUAdata'], 'InclUAdata', '! Is unsteady aerodynamics data included in this table? If TRUE, then include 30 UA coefficients below this line\n'))
                f.write('!........................................\n')
                if self.fst_vt['AeroDyn']['af_data'][afi][tab]['InclUAdata']:
                    f.write('{:<22f}   {:<11} {:}'.format(self.fst_vt['AeroDyn']['af_data'][afi][tab]['alpha0'], 'alpha0', '! 0-lift angle of attack, depends on airfoil.\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(self.fst_vt['AeroDyn']['af_data'][afi][tab]['alpha1'], 'alpha1', '! Angle of attack at f=0.7, (approximately the stall angle) for AOA>alpha0. (deg)\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(self.fst_vt['AeroDyn']['af_data'][afi][tab]['alpha2'], 'alpha2', '! Angle of attack at f=0.7, (approximately the stall angle) for AOA<alpha0. (deg)\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(self.fst_vt['AeroDyn']['af_data'][afi][tab]['eta_e'], 'eta_e', '! Recovery factor in the range [0.85 - 0.95] used only for UA_Mod=1, it is set to 1 in the code when flookup=True. (-)\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(self.fst_vt['AeroDyn']['af_data'][afi][tab]['C_nalpha'], 'C_nalpha', '! Slope of the 2D normal force coefficient curve. (1/rad)\n'))
                    f.write(float_default_out(self.fst_vt['AeroDyn']['af_data'][afi][tab]['T_f0']) + '   {:<11} {:}'.format('T_f0', '! Initial value of the time constant associated with Df in the expression of Df and f''. [default = 3]\n'))
                    f.write(float_default_out(self.fst_vt['AeroDyn']['af_data'][afi][tab]['T_V0']) + '   {:<11} {:}'.format('T_V0', '! Initial value of the time constant associated with the vortex lift decay process; it is used in the expression of Cvn. It depends on Re,M, and airfoil class. [default = 6]\n'))
                    f.write(float_default_out(self.fst_vt['AeroDyn']['af_data'][afi][tab]['T_p']) + '   {:<11} {:}'.format('T_p', '! Boundary-layer,leading edge pressure gradient time constant in the expression of Dp. It should be tuned based on airfoil experimental data. [default = 1.7]\n'))
                    f.write(float_default_out(self.fst_vt['AeroDyn']['af_data'][afi][tab]['T_VL']) + '   {:<11} {:}'.format('T_VL', '! Initial value of the time constant associated with the vortex advection process; it represents the non-dimensional time in semi-chords, needed for a vortex to travel from LE to trailing edge (TE); it is used in the expression of Cvn. It depends on Re, M (weakly), and airfoil. [valid range = 6 - 13, default = 11]\n'))
                    f.write(float_default_out(self.fst_vt['AeroDyn']['af_data'][afi][tab]['b1']) + '   {:<11} {:}'.format('b1', '! Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils, but may be different for turbine airfoils. [from experimental results, defaults to 0.14]\n'))
                    f.write(float_default_out(self.fst_vt['AeroDyn']['af_data'][afi][tab]['b2']) + '   {:<11} {:}'.format('b2', '! Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils, but may be different for turbine airfoils. [from experimental results, defaults to 0.53]\n'))
                    f.write(float_default_out(self.fst_vt['AeroDyn']['af_data'][afi][tab]['b5']) + '   {:<11} {:}'.format('b5', "! Constant in the expression of K'''_q,Cm_q^nc, and k_m,q.  [from  experimental results, defaults to 5]\n"))
                    f.write(float_default_out(self.fst_vt['AeroDyn']['af_data'][afi][tab]['A1']) + '   {:<11} {:}'.format('A1', '! Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils, but may be different for turbine airfoils. [from experimental results, defaults to 0.3]\n'))
                    f.write(float_default_out(self.fst_vt['AeroDyn']['af_data'][afi][tab]['A2']) + '   {:<11} {:}'.format('A2', '! Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils, but may be different for turbine airfoils. [from experimental results, defaults to 0.7]\n'))
                    f.write(float_default_out(self.fst_vt['AeroDyn']['af_data'][afi][tab]['A5']) + '   {:<11} {:}'.format('A5', "! Constant in the expression of K'''_q,Cm_q^nc, and k_m,q. [from experimental results, defaults to 1]\n"))
                    f.write('{:<22f}   {:<11} {:}'.format(self.fst_vt['AeroDyn']['af_data'][afi][tab]['S1'], 'S1', '! Constant in the f curve best-fit for alpha0<=AOA<=alpha1; by definition it depends on the airfoil. [ignored if UA_Mod<>1]\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(self.fst_vt['AeroDyn']['af_data'][afi][tab]['S2'], 'S2', '! Constant in the f curve best-fit for         AOA> alpha1; by definition it depends on the airfoil. [ignored if UA_Mod<>1]\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(self.fst_vt['AeroDyn']['af_data'][afi][tab]['S3'], 'S3', '! Constant in the f curve best-fit for alpha2<=AOA< alpha0; by definition it depends on the airfoil. [ignored if UA_Mod<>1]\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(self.fst_vt['AeroDyn']['af_data'][afi][tab]['S4'], 'S4', '! Constant in the f curve best-fit for         AOA< alpha2; by definition it depends on the airfoil. [ignored if UA_Mod<>1]\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(self.fst_vt['AeroDyn']['af_data'][afi][tab]['Cn1'], 'Cn1', '! Critical value of C0n at leading edge separation. It should be extracted from airfoil data at a given Mach and Reynolds number. It can be calculated from the static value of Cn at either the break in the pitching moment or the loss of chord force at the onset of stall. It is close to the condition of maximum lift of the airfoil at low Mach numbers.\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(self.fst_vt['AeroDyn']['af_data'][afi][tab]['Cn2'], 'Cn2', '! As Cn1 for negative AOAs.\n'))
                    # f.write('{: 22f}   {:<11} {:}'.format(self.fst_vt['AeroDyn']['af_data'][afi]['St_sh'], 'St_sh', "! Strouhal's shedding frequency constant.  [default = 0.19]\n"))
                    f.write(float_default_out(self.fst_vt['AeroDyn']['af_data'][afi][tab]['St_sh']) + '   {:<11} {:}'.format('St_sh', "! Strouhal's shedding frequency constant.  [default = 0.19]\n"))
                    f.write('{:<22f}   {:<11} {:}'.format(self.fst_vt['AeroDyn']['af_data'][afi][tab]['Cd0'], 'Cd0', '! 2D drag coefficient value at 0-lift.\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(self.fst_vt['AeroDyn']['af_data'][afi][tab]['Cm0'], 'Cm0', '! 2D pitching moment coefficient about 1/4-chord location, at 0-lift, positive if nose up. [If the aerodynamics coefficients table does not include a column for Cm, this needs to be set to 0.0]\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(self.fst_vt['AeroDyn']['af_data'][afi][tab]['k0'], 'k0', '! Constant in the \\hat(x)_cp curve best-fit; = (\\hat(x)_AC-0.25).  [ignored if UA_Mod<>1]\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(self.fst_vt['AeroDyn']['af_data'][afi][tab]['k1'], 'k1', '! Constant in the \\hat(x)_cp curve best-fit.  [ignored if UA_Mod<>1]\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(self.fst_vt['AeroDyn']['af_data'][afi][tab]['k2'], 'k2', '! Constant in the \\hat(x)_cp curve best-fit.  [ignored if UA_Mod<>1]\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(self.fst_vt['AeroDyn']['af_data'][afi][tab]['k3'], 'k3', '! Constant in the \\hat(x)_cp curve best-fit.  [ignored if UA_Mod<>1]\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(self.fst_vt['AeroDyn']['af_data'][afi][tab]['k1_hat'], 'k1_hat', '! Constant in the expression of Cc due to leading edge vortex effects.  [ignored if UA_Mod<>1]\n'))
                    f.write(float_default_out(self.fst_vt['AeroDyn']['af_data'][afi][tab]['x_cp_bar']) + '   {:<11} {:}'.format('x_cp_bar', '! Constant in the expression of \\hat(x)_cp^v. [ignored if UA_Mod<>1, default = 0.2]\n'))
                    f.write(float_default_out(self.fst_vt['AeroDyn']['af_data'][afi][tab]['UACutout']) + '   {:<11} {:}'.format('UACutout', '! Angle of attack above which unsteady aerodynamics are disabled (deg). [Specifying the string "Default" sets UACutout to 45 degrees]\n'))
                    f.write(float_default_out(self.fst_vt['AeroDyn']['af_data'][afi][tab]['filtCutOff']) + '   {:<11} {:}'.format('filtCutOff', '! Reduced frequency cut-off for low-pass filtering the AoA input to UA, as well as the 1st and 2nd derivatives (-) [default = 0.5]\n'))

                f.write('!........................................\n')
                f.write('! Table of aerodynamics coefficients\n')
                f.write('{:<22d}   {:<11} {:}'.format(self.fst_vt['AeroDyn']['af_data'][afi][tab]['NumAlf'], 'NumAlf', '! Number of data lines in the following table\n'))
                f.write('!    Alpha      Cl      Cd        Cm\n')
                f.write('!    (deg)      (-)     (-)       (-)\n')

                polar_map = [self.fst_vt['AeroDyn']['InCol_Alfa'], self.fst_vt['AeroDyn']['InCol_Cl'], self.fst_vt['AeroDyn']['InCol_Cd'], self.fst_vt['AeroDyn']['InCol_Cm'], self.fst_vt['AeroDyn']['InCol_Cpmin']]
                polar_map.remove(0)
                polar_map = [i-1 for i in polar_map]

                alpha = np.asarray(self.fst_vt['AeroDyn']['af_data'][afi][tab]['Alpha'])
                cl = np.asarray(self.fst_vt['AeroDyn']['af_data'][afi][tab]['Cl'])
                cd = np.asarray(self.fst_vt['AeroDyn']['af_data'][afi][tab]['Cd'])
                cm = np.asarray(self.fst_vt['AeroDyn']['af_data'][afi][tab]['Cm'])
                cpmin = np.asarray(self.fst_vt['AeroDyn']['af_data'][afi][tab]['Cpmin'])

                if alpha[0] != -180.:
                    print('Airfoil number ' + str(afi) + ' tab number ' + str(tab) + ' has the min angle of attack different than -180 deg, and equal to ' + str(alpha[0]) + ' deg. This is changed to -180 deg now.')
                    alpha[0] = -180.
                if alpha[-1] != 180.:
                    print('Airfoil number ' + str(afi) + ' tab number ' + str(tab) + ' has the max angle of attack different than 180 deg, and equal to ' + str(alpha[0]) + ' deg. This is changed to 180 deg now.')
                    alpha[-1] = 180.
                if cl[0] != cl[-1]:
                    print('Airfoil number ' + str(afi) + ' tab number ' + str(tab) + ' has the lift coefficient different between +-180 deg. This is changed to be the same now.')
                    cl[0] = cl[-1]
                if cd[0] != cd[-1]:
                    print('Airfoil number ' + str(afi) + ' tab number ' + str(tab) + ' has the drag coefficient different between +-180 deg. This is changed to be the same now.')
                    cd[0] = cd[-1]
                if cm[0] != cm[-1]:
                    print('Airfoil number ' + str(afi) + ' tab number ' + str(tab) + ' has the moment coefficient different between +-180 deg. This is changed to be the same now.')
                    cm[0] = cm[-1]

                if self.fst_vt['AeroDyn']['InCol_Cm'] == 0:
                    cm = np.zeros_like(cl)
                if self.fst_vt['AeroDyn']['InCol_Cpmin'] == 0:
                    cpmin = np.zeros_like(cl)
                polar = np.column_stack((alpha, cl, cd, cm, cpmin))
                polar = polar[:,polar_map]


                for row in polar:
                    f.write(' '.join(['{: 2.14e}'.format(val) for val in row])+'\n')
            
            f.flush()
            os.fsync(f)
            f.close()
            
    def write_AeroDynCoord(self, af_coord):

        self.fst_vt['AeroDyn']['AFNames_coord'] = ['']*self.fst_vt['AeroDyn']['NumAFfiles']
        
        for afi in af_coord: 
            self.fst_vt['AeroDyn']['AFNames_coord'][afi] = os.path.join('Airfoils', self.FAST_namingOut + '_AF%02d_Coords.txt'%afi)
            
            x     = self.fst_vt['AeroDyn']['af_coord'][afi]['x']
            y     = self.fst_vt['AeroDyn']['af_coord'][afi]['y']
            coord = np.vstack((x, y)).T

            af_file = os.path.join(self.FAST_runDirectory, self.fst_vt['AeroDyn']['AFNames_coord'][afi])
            f = open(af_file, 'w')
            
            f.write('{: 22d}   {:<11} {:}'.format(len(x)+1, 'NumCoords', '! The number of coordinates in the airfoil shape file (including an extra coordinate for airfoil reference).  Set to zero if coordinates not included.\n'))
            f.write('! ......... x-y coordinates are next if NumCoords > 0 .............\n')
            f.write('! x-y coordinate of airfoil reference\n')
            f.write('!  x/c        y/c\n')
            f.write('{: 5f}       0\n'.format(self.fst_vt['AeroDyn']['ac'][afi]))
            f.write('! coordinates of airfoil shape\n')
            f.write('! interpolation to 200 points\n')
            f.write('!  x/c        y/c\n')
            for row in coord:
                f.write(' '.join(['{: 2.14e}'.format(val) for val in row])+'\n')
            
            f.flush()
            os.fsync(f)
            f.close()

    def write_OLAF(self):

        olaf_file = os.path.join(self.FAST_runDirectory, self.FAST_namingOut + '_OLAF.dat')
        f = open(olaf_file, 'w')
        
        f.write('--------------------------- OLAF (cOnvecting LAgrangian Filaments) INPUT FILE -----------------\n')
        f.write('Generated by OpenFAST_IO\n')
        f.write('--------------------------- GENERAL OPTIONS ---------------------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['OLAF']['IntMethod'], 'IntMethod', '- Integration method {1: RK4, 5: Forward Euler 1st order, default: 5} (switch)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['OLAF']['DTfvw'], 'DTfvw', '- Time interval for wake propagation. {default: dtaero} (s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['OLAF']['FreeWakeStart'], 'FreeWakeStart', '- Time when wake is free. (-) value = always free. {default: 0.0} (s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['OLAF']['FullCircStart'], 'FullCircStart', '- Time at which full circulation is reached. {default: 0.0} (s)\n'))
        f.write('--------------------------- CIRCULATION SPECIFICATIONS ----------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['OLAF']['CircSolvMethod'], 'CircSolvingMethod', '- Circulation solving method {1: Cl-Based, 2: No-Flow Through, 3: Prescribed, default: 1 }(switch)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['OLAF']['CircSolvConvCrit'], 'CircSolvConvCrit', ' - Convergence criteria {default: 0.001} [only if CircSolvMethod=1] (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['OLAF']['CircSolvRelaxation'], 'CircSolvRelaxation', '- Relaxation factor {default: 0.1} [only if CircSolvMethod=1] (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['OLAF']['CircSolvMaxIter'], 'CircSolvMaxIter', ' - Maximum number of iterations for circulation solving {default: 30} (-)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['AeroDyn']['OLAF']['PrescribedCircFile']+'"', 'PrescribedCircFile','- File containing prescribed circulation [only if CircSolvMethod=3] (quoted string)\n'))
        f.write('===============================================================================================\n')
        f.write('--------------------------- WAKE OPTIONS ------------------------------------------------------\n')
        f.write('------------------- WAKE EXTENT AND DISCRETIZATION --------------------------------------------\n')
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['AeroDyn']['OLAF']['nNWPanels'], 'nNWPanels','- Number of near-wake panels (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['OLAF']['nNWPanelsFree'], 'nNWPanelsFree','- Number of free near-wake panels (-) {default: nNWPanels}\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['AeroDyn']['OLAF']['nFWPanels'], 'nFWPanels','- Number of far-wake panels (-) {default: 0}\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['OLAF']['nFWPanelsFree'], 'nFWPanelsFree','- Number of free far-wake panels (-) {default: nFWPanels}\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['OLAF']['FWShedVorticity'], 'FWShedVorticity','- Include shed vorticity in the far wake {default: False}\n'))
        f.write('------------------- WAKE REGULARIZATIONS AND DIFFUSION -----------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['OLAF']['DiffusionMethod'], 'DiffusionMethod','- Diffusion method to account for viscous effects {0: None, 1: Core Spreading, "default": 0}\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['OLAF']['RegDeterMethod'], 'RegDeterMethod','- Method to determine the regularization parameters {0:  Manual, 1: Optimized, 2: Chord, 3: Span, default: 0 }\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['OLAF']['RegFunction'], 'RegFunction','- Viscous diffusion function {0: None, 1: Rankine, 2: LambOseen, 3: Vatistas, 4: Denominator, "default": 3} (switch)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['OLAF']['WakeRegMethod'], 'WakeRegMethod','- Wake regularization method {1: Constant, 2: Stretching, 3: Age, default: 3} (switch)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['OLAF']['WakeRegFactor'], 'WakeRegFactor','- Wake regularization factor (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['OLAF']['WingRegFactor'], 'WingRegFactor','- Wing regularization factor (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['OLAF']['CoreSpreadEddyVisc'], 'CoreSpreadEddyVisc','- Eddy viscosity in core spreading methods, typical values 1-1000\n'))
        f.write('------------------- WAKE TREATMENT OPTIONS ---------------------------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['OLAF']['TwrShadowOnWake'], 'TwrShadowOnWake','- Include tower flow disturbance effects on wake convection {default:false} [only if TwrPotent or TwrShadow]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['OLAF']['ShearModel'], 'ShearModel','- Shear Model {0: No treatment, 1: Mirrored vorticity, default: 0}\n'))
        f.write('------------------- SPEEDUP OPTIONS -----------------------------------------------------------\n')
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['AeroDyn']['OLAF']['VelocityMethod'], 'VelocityMethod','- Method to determine the velocity {1:Segment N^2, 2:Particle tree, 3:Particle N^2, 4:Segment tree, default: 2}\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['OLAF']['TreeBranchFactor'], 'TreeBranchFactor','- Branch radius fraction above which a multipole calculation is used {default: 1.5} [only if VelocityMethod=2,4]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['OLAF']['PartPerSegment'], 'PartPerSegment','- Number of particles per segment {default: 1} [only if VelocityMethod=2,3]\n'))
        f.write('===============================================================================================\n')
        f.write('--------------------------- OUTPUT OPTIONS  ---------------------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['OLAF']['WrVTk'], 'WrVTk','- Outputs Visualization Toolkit (VTK) (independent of .fst option) {0: NoVTK, 1: Write VTK at VTK_fps, 2: Write VTK at init and final, default: 0} (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['OLAF']['nVTKBlades'], 'nVTKBlades','- Number of blades for which VTK files are exported {0: No VTK per blade, n: VTK for blade 1 to n, default: 0} (-) \n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['OLAF']['VTKCoord'], 'VTKCoord','- Coordinate system used for VTK export. {1: Global, 2: Hub, 3: Both, default: 1} \n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['OLAF']['VTK_fps'], 'VTK_fps','- Frame rate for VTK output (frames per second) {"all" for all glue code timesteps, "default" for all OLAF timesteps} [only if WrVTK=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDyn']['OLAF']['nGridOut'], 'nGridOut','- Number of grid outputs\n'))
        f.write('GridName  GridType  TStart   TEnd    DTGrid   XStart  XEnd   nX   YStart   YEnd   nY   ZStart   ZEnd   nZ\n')
        f.write('(-)         (-)      (s)     (s)      (s)       (m)    (m)   (-)   (m)     (m)    (-)   (m)     (m)    (-)\n')
        f.write('===============================================================================================\n')
        f.write('--------------------------- ADVANCED OPTIONS --------------------------------------------------\n')
        f.write('===============================================================================================\n')

        f.flush()
        os.fsync(f)
        f.close()
    
    def write_AeroDisk(self):
        # Writing the aeroDisk input file
        self.fst_vt['Fst']['AeroFile'] = self.FAST_namingOut + '_AeroDisk.dat'
        adisk_file = os.path.join(self.FAST_runDirectory,self.fst_vt['Fst']['AeroFile'])
        f = open(adisk_file,'w')

        f.write('--- AERO DISK INPUT FILE -------\n')
        f.write('Generated with OpenFAST_IO\n')
        f.write('--- SIMULATION CONTROL ---------\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['AeroDisk']['Echo'], 'Echo', '- Echo input data to "<RootName>.ADsk.ech" (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['AeroDisk']['DT'], 'DT', '- Integration time step (s)\n'))   
        f.write('--- ENVIRONMENTAL CONDITIONS ---\n')
        f.write('{:<22f} {:<11} {:}'.format(self.fst_vt['AeroDisk']['AirDens'], 'AirDens', '- Air density (kg/m^3) (or "default")\n'))
        f.write('--- ACTUATOR DISK PROPERTIES ---\n')
        f.write('{:<22f} {:<11} {:}'.format(self.fst_vt['AeroDisk']['RotorRad'], 'RotorRad', '- Rotor radius (m) (or "default")\n'))
        f.write('"{:<22}" {:<11} {:}'.format(', '.join(['%s'%i for i in self.fst_vt['AeroDisk']['InColNames']]), 'InColNames', '- Input column headers (string) {may include a combination of "TSR, RtSpd, VRel, Pitch, Skew"} (up to 4 columns) [choose TSR or RtSpd,VRel; if Skew is absent, Skew is modeled as (COS(Skew))^2]\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join(['%s'%i for i in self.fst_vt['AeroDisk']['InColDims']]), 'InColDims', '- Number of unique values in each column (-) (must have same number of columns as InColName) [each >=2]\n'))
        self.write_AeroDiskProp()
        f.write('@{:<22} {:}'.format(self.fst_vt['AeroDisk']['actuatorDiskFile'], '\n'))
        f.write('--- OUTPUTS --------------------\n')
        f.write('{:<22} {:<11} {:}'.format('OutList', 'OutList', '- The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)\n'))

        outlist = self.get_outlist(self.fst_vt['outlist'], ['AeroDisk'])
        for channel_list in outlist:
            for i in range(len(channel_list)):
                f.write('"' + channel_list[i] + '"\n')

        f.write('END of input file (the word "END" must appear in the first 3 columns of the last OutList line)\n')
        f.write('---------------------------------------------------------------------------------------\n')

        f.flush()
        os.fsync(f)
        f.close()

    def write_AeroDiskProp(self):
        # Writing the aeroDiskProp input file

        self.fst_vt['AeroDisk']['actuatorDiskFile'] = self.FAST_namingOut + '_AeroDiskProp.csv'
        adiskprop_file = os.path.join(self.FAST_runDirectory,self.fst_vt['AeroDisk']['actuatorDiskFile'])
        f = open(adiskprop_file,'w')

        f.write('{:<22} {:}'.format(self.fst_vt['AeroDisk']['actuatorDiskTable']['dsc'],'\n'))
        f.write('{:<22} {:}'.format(', '.join(['%s'%i for i in self.fst_vt['AeroDisk']['actuatorDiskTable']['attr']]), '\n'))
        f.write('{:<22} {:}'.format(', '.join(['%s'%i for i in self.fst_vt['AeroDisk']['actuatorDiskTable']['units']]), '\n'))
        for idx in range(len(self.fst_vt['AeroDisk']['actuatorDiskTable']['data'])):
            f.write('{:<22} {:}'.format(', '.join(['%.6f'%i for i in self.fst_vt['AeroDisk']['actuatorDiskTable']['data'][idx]]), '\n'))
        
        f.flush()
        os.fsync(f)
        f.close()



    def write_ServoDyn(self):
        # ServoDyn v1.05 Input File

        self.fst_vt['Fst']['ServoFile'] = self.FAST_namingOut + '_ServoDyn.dat'
        sd_file = os.path.join(self.FAST_runDirectory,self.fst_vt['Fst']['ServoFile'])
        f = open(sd_file,'w')

        f.write('------- SERVODYN INPUT FILE --------------------------------------------\n')
        f.write('Generated with OpenFAST_IO\n')
        f.write('---------------------- SIMULATION CONTROL --------------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['Echo'], 'Echo', '- Echo input data to <RootName>.ech (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['DT'], 'DT', '- Communication interval for controllers (s) (or "default")\n'))
        f.write('---------------------- PITCH CONTROL -------------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['PCMode'], 'PCMode', '- Pitch control mode {0: none, 3: user-defined from routine PitchCntrl, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['TPCOn'], 'TPCOn', '- Time to enable active pitch control (s) [unused when PCMode=0]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['TPitManS1'], 'TPitManS(1)', '- Time to start override pitch maneuver for blade 1 and end standard pitch control (s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['TPitManS2'], 'TPitManS(2)', '- Time to start override pitch maneuver for blade 2 and end standard pitch control (s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['TPitManS3'], 'TPitManS(3)', '- Time to start override pitch maneuver for blade 3 and end standard pitch control (s) [unused for 2 blades]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['PitManRat(1)'], 'PitManRat(1)', '- Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 1 (deg/s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['PitManRat(2)'], 'PitManRat(2)', '- Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 2 (deg/s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['PitManRat(3)'], 'PitManRat(3)', '- Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 3 (deg/s) [unused for 2 blades]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['BlPitchF(1)'], 'BlPitchF(1)', '- Blade 1 final pitch for pitch maneuvers (degrees)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['BlPitchF(2)'], 'BlPitchF(2)', '- Blade 2 final pitch for pitch maneuvers (degrees)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['BlPitchF(3)'], 'BlPitchF(3)', '- Blade 3 final pitch for pitch maneuvers (degrees) [unused for 2 blades]\n'))
        f.write('---------------------- GENERATOR AND TORQUE CONTROL ----------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['VSContrl'], 'VSContrl', '- Variable-speed control mode {0: none, 1: simple VS, 3: user-defined from routine UserVSCont, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['GenModel'], 'GenModel', '- Generator model {1: simple, 2: Thevenin, 3: user-defined from routine UserGen} (switch) [used only when VSContrl=0]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['GenEff'], 'GenEff', '- Generator efficiency [ignored by the Thevenin and user-defined generator models] (%)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['GenTiStr'], 'GenTiStr', '- Method to start the generator {T: timed using TimGenOn, F: generator speed using SpdGenOn} (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['GenTiStp'], 'GenTiStp', '- Method to stop the generator {T: timed using TimGenOf, F: when generator power = 0} (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['SpdGenOn'], 'SpdGenOn', '- Generator speed to turn on the generator for a startup (HSS speed) (rpm) [used only when GenTiStr=False]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['TimGenOn'], 'TimGenOn', '- Time to turn on the generator for a startup (s) [used only when GenTiStr=True]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['TimGenOf'], 'TimGenOf', '- Time to turn off the generator (s) [used only when GenTiStp=True]\n'))
        f.write('---------------------- SIMPLE VARIABLE-SPEED TORQUE CONTROL --------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['VS_RtGnSp'], 'VS_RtGnSp', '- Rated generator speed for simple variable-speed generator control (HSS side) (rpm) [used only when VSContrl=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['VS_RtTq'], 'VS_RtTq', '- Rated generator torque/constant generator torque in Region 3 for simple variable-speed generator control (HSS side) (N-m) [used only when VSContrl=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['VS_Rgn2K'], 'VS_Rgn2K', '- Generator torque constant in Region 2 for simple variable-speed generator control (HSS side) (N-m/rpm^2) [used only when VSContrl=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['VS_SlPc'], 'VS_SlPc', '- Rated generator slip percentage in Region 2 1/2 for simple variable-speed generator control (%) [used only when VSContrl=1]\n'))
        f.write('---------------------- SIMPLE INDUCTION GENERATOR ------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['SIG_SlPc'], 'SIG_SlPc', '- Rated generator slip percentage (%) [used only when VSContrl=0 and GenModel=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['SIG_SySp'], 'SIG_SySp', '- Synchronous (zero-torque) generator speed (rpm) [used only when VSContrl=0 and GenModel=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['SIG_RtTq'], 'SIG_RtTq', '- Rated torque (N-m) [used only when VSContrl=0 and GenModel=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['SIG_PORt'], 'SIG_PORt', '- Pull-out ratio (Tpullout/Trated) (-) [used only when VSContrl=0 and GenModel=1]\n'))
        f.write('---------------------- THEVENIN-EQUIVALENT INDUCTION GENERATOR -----------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['TEC_Freq'], 'TEC_Freq', '- Line frequency [50 or 60] (Hz) [used only when VSContrl=0 and GenModel=2]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['TEC_NPol'], 'TEC_NPol', '- Number of poles [even integer > 0] (-) [used only when VSContrl=0 and GenModel=2]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['TEC_SRes'], 'TEC_SRes', '- Stator resistance (ohms) [used only when VSContrl=0 and GenModel=2]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['TEC_RRes'], 'TEC_RRes', '- Rotor resistance (ohms) [used only when VSContrl=0 and GenModel=2]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['TEC_VLL'], 'TEC_VLL', '- Line-to-line RMS voltage (volts) [used only when VSContrl=0 and GenModel=2]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['TEC_SLR'], 'TEC_SLR', '- Stator leakage reactance (ohms) [used only when VSContrl=0 and GenModel=2]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['TEC_RLR'], 'TEC_RLR', '- Rotor leakage reactance (ohms) [used only when VSContrl=0 and GenModel=2]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['TEC_MR'], 'TEC_MR', '- Magnetizing reactance (ohms) [used only when VSContrl=0 and GenModel=2]\n'))
        f.write('---------------------- HIGH-SPEED SHAFT BRAKE ----------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['HSSBrMode'], 'HSSBrMode', '- HSS brake model {0: none, 1: simple, 3: user-defined from routine UserHSSBr, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['THSSBrDp'], 'THSSBrDp', '- Time to initiate deployment of the HSS brake (s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['HSSBrDT'], 'HSSBrDT', '- Time for HSS-brake to reach full deployment once initiated (sec) [used only when HSSBrMode=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['HSSBrTqF'], 'HSSBrTqF', '- Fully deployed HSS-brake torque (N-m)\n'))
        f.write('---------------------- NACELLE-YAW CONTROL -------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['YCMode'], 'YCMode', '- Yaw control mode {0: none, 3: user-defined from routine UserYawCont, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['TYCOn'], 'TYCOn', '- Time to enable active yaw control (s) [unused when YCMode=0]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['YawNeut'], 'YawNeut', '- Neutral yaw position--yaw spring force is zero at this yaw (degrees)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['YawSpr'], 'YawSpr', '- Nacelle-yaw spring constant (N-m/rad)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['YawDamp'], 'YawDamp', '- Nacelle-yaw damping constant (N-m/(rad/s))\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['TYawManS'], 'TYawManS', '- Time to start override yaw maneuver and end standard yaw control (s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['YawManRat'], 'YawManRat', '- Yaw maneuver rate (in absolute value) (deg/s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['NacYawF'], 'NacYawF', '- Final yaw angle for override yaw maneuvers (degrees)\n'))
        f.write('---------------------- Aerodynamic Flow Control -------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['AfCmode'], 'AfCmode', '- Airfoil control mode {0: none, 1: cosine wave cycle, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['AfC_Mean'], 'AfC_Mean', '- Mean level for cosine cycling or steady value (-) [used only with AfCmode==1]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['AfC_Amp'], 'AfC_Amp', '- Amplitude for for cosine cycling of flap signal (AfC = AfC_Amp*cos(Azimuth+phase)+AfC_mean) (-) [used only with AfCmode==1]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['AfC_Phase'], 'AfC_phase', '- Phase relative to the blade azimuth (0 is vertical) for for cosine cycling of flap signal (deg) [used only with AfCmode==1]\n'))
        f.write('---------------------- STRUCTURAL CONTROL ---------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['NumBStC'], 'NumBStC', '- Number of blade structural controllers (integer)\n'))
        f.write('{!s:<22} {:<11} {:}'.format('"' + '" "'.join(self.fst_vt['ServoDyn']['BStCfiles']) + '"', 'BStCfiles', '- Name of the files for blade structural controllers (quoted strings) [unused when NumBStC==0]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['NumNStC'], 'NumNStC', '- Number of nacelle structural controllers (integer)\n'))
        f.write('{!s:<22} {:<11} {:}'.format('"' + '" "'.join(self.fst_vt['ServoDyn']['NStCfiles']) + '"', 'NStCfiles', '- Name of the files for nacelle structural controllers (quoted strings) [unused when NumNStC==0]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['NumTStC'], 'NumTStC', '- Number of tower structural controllers (integer)\n'))
        f.write('{!s:<22} {:<11} {:}'.format('"' + '" "'.join(self.fst_vt['ServoDyn']['TStCfiles']) + '"', 'TStCfiles', '- Name of the files for tower structural controllers (quoted strings) [unused when NumTStC==0]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['NumSStC'], 'NumSStC', '- Number of substructure structural controllers (integer)\n'))
        f.write('{!s:<22} {:<11} {:}'.format('"' + '" "'.join(self.fst_vt['ServoDyn']['SStCfiles']) + '"', 'SStCfiles', '- Name of the files for substructure structural controllers (quoted strings) [unused when NumSStC==0]\n'))
        f.write('---------------------- CABLE CONTROL ---------------------------------------- \n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['CCmode'], 'CCmode', '- Cable control mode {0: none, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n'))
        f.write('---------------------- BLADED INTERFACE ---------------------------------------- [used only with Bladed Interface]\n')
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['ServoDyn']['DLL_FileName']+'"', 'DLL_FileName', '- Name/location of the dynamic library {.dll [Windows] or .so [Linux]} in the Bladed-DLL format (-) [used only with Bladed Interface]\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['ServoDyn']['DLL_InFile']+'"', 'DLL_InFile', '- Name of input file sent to the DLL (-) [used only with Bladed Interface]\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['ServoDyn']['DLL_ProcName']+'"', 'DLL_ProcName', '- Name of procedure in DLL to be called (-) [case sensitive; used only with DLL Interface]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['DLL_DT'], 'DLL_DT', '- Communication interval for dynamic library (s) (or "default") [used only with Bladed Interface]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['DLL_Ramp'], 'DLL_Ramp', '- Whether a linear ramp should be used between DLL_DT time steps [introduces time shift when true] (flag) [used only with Bladed Interface]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['BPCutoff'], 'BPCutoff', '- Cutoff frequency for low-pass filter on blade pitch from DLL (Hz) [used only with Bladed Interface]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['NacYaw_North'], 'NacYaw_North', '- Reference yaw angle of the nacelle when the upwind end points due North (deg) [used only with Bladed Interface]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['Ptch_Cntrl'], 'Ptch_Cntrl', '- Record 28: Use individual pitch control {0: collective pitch; 1: individual pitch control} (switch) [used only with Bladed Interface]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['Ptch_SetPnt'], 'Ptch_SetPnt', '- Record  5: Below-rated pitch angle set-point (deg) [used only with Bladed Interface]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['Ptch_Min'], 'Ptch_Min', '- Record  6: Minimum pitch angle (deg) [used only with Bladed Interface]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['Ptch_Max'], 'Ptch_Max', '- Record  7: Maximum pitch angle (deg) [used only with Bladed Interface]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['PtchRate_Min'], 'PtchRate_Min', '- Record  8: Minimum pitch rate (most negative value allowed) (deg/s) [used only with Bladed Interface]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['PtchRate_Max'], 'PtchRate_Max', '- Record  9: Maximum pitch rate  (deg/s) [used only with Bladed Interface]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['Gain_OM'], 'Gain_OM', '- Record 16: Optimal mode gain (Nm/(rad/s)^2) [used only with Bladed Interface]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['GenSpd_MinOM'], 'GenSpd_MinOM', '- Record 17: Minimum generator speed (rpm) [used only with Bladed Interface]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['GenSpd_MaxOM'], 'GenSpd_MaxOM', '- Record 18: Optimal mode maximum speed (rpm) [used only with Bladed Interface]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['GenSpd_Dem'], 'GenSpd_Dem', '- Record 19: Demanded generator speed above rated (rpm) [used only with Bladed Interface]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['GenTrq_Dem'], 'GenTrq_Dem', '- Record 22: Demanded generator torque above rated (Nm) [used only with Bladed Interface]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['GenPwr_Dem'], 'GenPwr_Dem', '- Record 13: Demanded power (W) [used only with Bladed Interface]\n'))
        f.write('---------------------- BLADED INTERFACE TORQUE-SPEED LOOK-UP TABLE -------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['DLL_NumTrq'], 'DLL_NumTrq', '- Record 26: No. of points in torque-speed look-up table {0 = none and use the optimal mode parameters; nonzero = ignore the optimal mode PARAMETERs by setting Record 16 to 0.0} (-) [used only with Bladed Interface]\n'))
        f.write('{:<22}\t{:<22}\n'.format("GenSpd_TLU", "GenTrq_TLU"))
        f.write('{:<22}\t{:<22}\n'.format("(rpm)", "(Nm)"))
        for i in range(self.fst_vt['ServoDyn']['DLL_NumTrq']):
            a1 = self.fst_vt['ServoDyn']['GenSpd_TLU'][i]
            a2 = self.fst_vt['ServoDyn']['GenTrq_TLU'][i]
            f.write('{:<22}\t{:<22}\n'.format(a1, a2))
        f.write('---------------------- OUTPUT --------------------------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['SumPrint'], 'SumPrint', '- Print summary data to <RootName>.sum (flag) (currently unused)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['OutFile'], 'OutFile', '- Switch to determine where output will be placed: {1: in module output file only; 2: in glue code output file only; 3: both} (currently unused)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['TabDelim'], 'TabDelim', '- Use tab delimiters in text tabular output file? (flag) (currently unused)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['ServoDyn']['OutFmt']+'"', 'OutFmt', '- Format used for text tabular output (except time).  Resulting field should be 10 characters. (quoted string) (currently unused)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['ServoDyn']['TStart'], 'TStart', '- Time to begin tabular output (s) (currently unused)\n'))
        f.write('              OutList      - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)\n')
        
        outlist = self.get_outlist(self.fst_vt['outlist'], ['ServoDyn'])
        for channel_list in outlist:
            for i in range(len(channel_list)):
                f.write('"' + channel_list[i] + '"\n')

        f.write('END of input file (the word "END" must appear in the first 3 columns of the last OutList line)\n')
        f.write('---------------------------------------------------------------------------------------\n')

        f.flush()
        os.fsync(f)
        f.close()

    def write_DISCON_in(self):
        # Generate Bladed style Interface controller input file, intended for ROSCO https://github.com/NREL/ROSCO_toolbox

        # Fill controller and turbine objects for ROSCO 
        # - controller
        controller = type('', (), {})()
        
        turbine = type('', (), {})()
        turbine.Cp = type('', (), {})()
        turbine.Ct = type('', (), {})()
        turbine.Cq = type('', (), {})()
        turbine.v_rated                 = self.fst_vt['DISCON_in']['v_rated']
        turbine.Cp                      = self.fst_vt['DISCON_in']['Cp']
        turbine.Ct                      = self.fst_vt['DISCON_in']['Ct']
        turbine.Cq                      = self.fst_vt['DISCON_in']['Cq']
        turbine.Cp_table                = self.fst_vt['DISCON_in']['Cp_table']
        turbine.Ct_table                = self.fst_vt['DISCON_in']['Ct_table']
        turbine.Cq_table                = self.fst_vt['DISCON_in']['Cq_table']
        turbine.pitch_initial_rad       = self.fst_vt['DISCON_in']['Cp_pitch_initial_rad']
        turbine.TSR_initial             = self.fst_vt['DISCON_in']['Cp_TSR_initial']
        turbine.TurbineName             = 'OpenFAST_IO Turbine'

        # Define DISCON infile paths
        self.fst_vt['ServoDyn']['DLL_InFile'] = self.FAST_namingOut + '_DISCON.IN'
        discon_in_file = os.path.join(self.FAST_runDirectory, self.fst_vt['ServoDyn']['DLL_InFile'])
        self.fst_vt['DISCON_in']['PerfFileName'] = self.FAST_namingOut + '_Cp_Ct_Cq.txt'
        
        # Write DISCON input files
        ROSCO_utilities.write_rotor_performance(
            turbine, 
            txt_filename=os.path.join(self.FAST_runDirectory, self.fst_vt['DISCON_in']['PerfFileName'])
            )
        
        ROSCO_utilities.write_DISCON(
            turbine,
            controller,
            param_file=discon_in_file, 
            txt_filename=self.fst_vt['DISCON_in']['PerfFileName'],
            rosco_vt=self.fst_vt['DISCON_in']
            )

    def write_spd_trq(self):
        # generate the spd_trq.dat file when VSContrl == 3
        spd_trq_file = os.path.join(self.FAST_runDirectory, 'spd_trq.dat')
        f = open(spd_trq_file, 'w')

        f.write('{:}'.format(self.fst_vt['spd_trq']['header'], '\n'))
        for i in range(len(self.fst_vt['spd_trq']['RPM'])):
            f.write('{:<22f} {:<22f} {:}'.format(self.fst_vt['spd_trq']['RPM'][i], self.fst_vt['spd_trq']['Torque'][i], '\n'))

    def write_HydroDyn(self):

        # Generate HydroDyn input file
        self.fst_vt['Fst']['HydroFile'] = self.FAST_namingOut + '_HydroDyn.dat'
        hd_file = os.path.join(self.FAST_runDirectory, self.fst_vt['Fst']['HydroFile'])
        f = open(hd_file, 'w')

        f.write('------- HydroDyn Input File --------------------------------------------\n')
        f.write('Generated with OpenFAST_IO\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['HydroDyn']['Echo'], 'Echo', '- Echo the input file data (flag)\n'))
        f.write('---------------------- FLOATING PLATFORM --------------------------------------- [unused with WaveMod=6]\n')
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['HydroDyn']['PotMod'], 'PotMod', '- Potential-flow model {0: none=no potential flow, 1: frequency-to-time-domain transforms based on WAMIT output, 2: fluid-impulse theory (FIT)} (switch)\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['HydroDyn']['ExctnMod'], 'ExctnMod', '- Wave-excitation model {0: no wave-excitation calculation, 1: DFT, 2: state-space} (switch) [only used when PotMod=1; STATE-SPACE REQUIRES *.ssexctn INPUT FILE]\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['HydroDyn']['ExctnDisp'], 'ExctnDisp','- Method of computing Wave Excitation {0: use undisplaced position, 1: use displaced position, 2: use low-pass filtered displaced position) [only used when PotMod=1 and ExctnMod>0 and SeaState\'s WaveMod>0]} (switch)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['HydroDyn']['ExctnCutOff'], 'ExctnCutOff','- Cutoff (corner) frequency of the low-pass time-filtered displaced position (Hz) [>0.0] [used only when PotMod=1, ExctnMod>0, and ExctnDisp=2]) [only used when PotMod=1 and ExctnMod>0 and SeaState\'s WaveMod>0]} (switch)\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['HydroDyn']['PtfmYMod'], 'PtfmYMod', '- Model for large platform yaw offset {0: Static reference yaw offset based on PtfmRefY, 1: dynamic reference yaw offset based on low-pass filtering the PRP yaw motion with cutoff frequency PtfmYCutOff} (switch)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['HydroDyn']['PtfmRefY'], 'PtfmRefY', '- Constant (if PtfmYMod=0) or initial (if PtfmYMod=1) platform reference yaw offset (deg)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['HydroDyn']['PtfmYCutOff'], 'PtfmYCutOff', '- Cutoff frequency for the low-pass filtering of PRP yaw motion when PtfmYMod=1 [unused when PtfmYMod=0] (Hz)\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['HydroDyn']['NExctnHdg'], 'NExctnHdg', '- Number of evenly distributed platform yaw/heading angles over the range of [-180, 180) deg for which the wave excitation shall be computed [only used when PtfmYMod=1] (-)\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['HydroDyn']['RdtnMod'], 'RdtnMod', '- Radiation memory-effect model {0: no memory-effect calculation, 1: convolution, 2: state-space} (switch) [only used when PotMod=1; STATE-SPACE REQUIRES *.ss INPUT FILE]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['HydroDyn']['RdtnTMax'], 'RdtnTMax', '- Analysis time for wave radiation kernel calculations (sec) [only used when PotMod=1 and RdtnMod>0; determines RdtnDOmega=Pi/RdtnTMax in the cosine transform; MAKE SURE THIS IS LONG ENOUGH FOR THE RADIATION IMPULSE RESPONSE FUNCTIONS TO DECAY TO NEAR-ZERO FOR THE GIVEN PLATFORM!]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['HydroDyn']['RdtnDT'], 'RdtnDT', '- Time step for wave radiation kernel calculations (sec) [only used when PotMod=1 and ExctnMod>0 or RdtnMod>0; DT<=RdtnDT<=0.1 recommended; determines RdtnOmegaMax=Pi/RdtnDT in the cosine transform]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['HydroDyn']['NBody'], 'NBody', '- Number of WAMIT bodies to be used (-) [>=1; only used when PotMod=1. If NBodyMod=1, the WAMIT data contains a vector of size 6*NBody x 1 and matrices of size 6*NBody x 6*NBody; if NBodyMod>1, there are NBody sets of WAMIT data each with a vector of size 6 x 1 and matrices of size 6 x 6]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['HydroDyn']['NBodyMod'], 'NBodyMod', '- Body coupling model {1: include coupling terms between each body and NBody in HydroDyn equals NBODY in WAMIT, 2: neglect coupling terms between each body and NBODY=1 with XBODY=0 in WAMIT, 3: Neglect coupling terms between each body and NBODY=1 with XBODY=/0 in WAMIT} (switch) [only used when PotMod=1]\n'))
        f.write('{:<22} {:<11} {:}'.format('"{}"'.format('", "'.join(self.fst_vt['HydroDyn']['PotFile']) if isinstance(self.fst_vt['HydroDyn']['PotFile'], list) else self.fst_vt['HydroDyn']['PotFile']), 'PotFile', '- Root name of potential-flow model data; WAMIT output files containing the linear, nondimensionalized, hydrostatic restoring matrix (.hst), frequency-dependent hydrodynamic added mass matrix and damping matrix (.1), and frequency- and direction-dependent wave excitation force vector per unit wave amplitude (.3) (quoted string) [1 to NBody if NBodyMod>1] [MAKE SURE THE FREQUENCIES INHERENT IN THESE WAMIT FILES SPAN THE PHYSICALLY-SIGNIFICANT RANGE OF FREQUENCIES FOR THE GIVEN PLATFORM; THEY MUST CONTAIN THE ZERO- AND INFINITE-FREQUENCY LIMITS!]\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join([f'{val}' for val in self.fst_vt['HydroDyn']['WAMITULEN']]), 'WAMITULEN', '- Characteristic body length scale used to redimensionalize WAMIT output (meters) [1 to NBody if NBodyMod>1] [only used when PotMod=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join([f'{val}' for val in self.fst_vt['HydroDyn']['PtfmRefxt']]), 'PtfmRefxt', '- The xt offset of the body reference point(s) from (0,0,0) (meters) [1 to NBody] [only used when PotMod=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join([f'{val}' for val in self.fst_vt['HydroDyn']['PtfmRefyt']]), 'PtfmRefyt', '- The yt offset of the body reference point(s) from (0,0,0) (meters) [1 to NBody] [only used when PotMod=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join([f'{val}' for val in self.fst_vt['HydroDyn']['PtfmRefzt']]), 'PtfmRefzt', '- The zt offset of the body reference point(s) from (0,0,0) (meters) [1 to NBody] [only used when PotMod=1. If NBodyMod=2,PtfmRefzt=0.0]\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join([f'{val}' for val in self.fst_vt['HydroDyn']['PtfmRefztRot']]), 'PtfmRefztRot', '- The rotation about zt of the body reference frame(s) from xt/yt (degrees) [1 to NBody] [only used when PotMod=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join([f'{val}' for val in self.fst_vt['HydroDyn']['PtfmVol0']]), 'PtfmVol0', '- Displaced volume of water when the body is in its undisplaced position (m^3) [1 to NBody] [only used when PotMod=1; USE THE SAME VALUE COMPUTED BY WAMIT AS OUTPUT IN THE .OUT FILE!]\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join([f'{val}' for val in self.fst_vt['HydroDyn']['PtfmCOBxt']]), 'PtfmCOBxt', '- The xt offset of the center of buoyancy (COB) from (0,0) (meters) [1 to NBody] [only used when PotMod=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join([f'{val}' for val in self.fst_vt['HydroDyn']['PtfmCOByt']]), 'PtfmCOByt', '- The yt offset of the center of buoyancy (COB) from (0,0) (meters) [1 to NBody] [only used when PotMod=1]\n'))
        f.write('---------------------- 2ND-ORDER FLOATING PLATFORM FORCES ---------------------- [unused with WaveMod=0 or 6, or PotMod=0 or 2]\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['HydroDyn']['MnDrift'], 'MnDrift', "- Mean-drift 2nd-order forces computed                                       {0: None; [7, 8, 9, 10, 11, or 12]: WAMIT file to use} [Only one of MnDrift, NewmanApp, or DiffQTF can be non-zero]\n"))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['HydroDyn']['NewmanApp'], 'NewmanApp', "- Mean- and slow-drift 2nd-order forces computed with Newman's approximation {0: None; [7, 8, 9, 10, 11, or 12]: WAMIT file to use} [Only one of MnDrift, NewmanApp, or DiffQTF can be non-zero. Used only when WaveDirMod=0]\n"))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['HydroDyn']['DiffQTF'], 'DiffQTF', "- Full difference-frequency 2nd-order forces computed with full QTF          {0: None; [10, 11, or 12]: WAMIT file to use}          [Only one of MnDrift, NewmanApp, or DiffQTF can be non-zero]\n"))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['HydroDyn']['SumQTF'], 'SumQTF', "- Full summation -frequency 2nd-order forces computed with full QTF          {0: None; [10, 11, or 12]: WAMIT file to use}\n"))
        f.write('---------------------- PLATFORM ADDITIONAL STIFFNESS AND DAMPING  -------------- [unused with PotMod=0 or 2]\n')
        for j in range(6):
            if type(self.fst_vt['HydroDyn']['AddF0'][j]) == float:
                ln = '{:14}   '.format(self.fst_vt['HydroDyn']['AddF0'][j])  
            elif type(self.fst_vt['HydroDyn']['AddF0'][j]) in [list, np.ndarray]:
                ln = '{:14}   '.format(' '.join([f'{val}' for val in self.fst_vt['HydroDyn']['AddF0'][j]]))
            else:
                raise Exception("Check type of self.fst_vt['HydroDyn']['AddF0']")
            
            if j == 0:
                ln = ln + 'AddF0    - Additional preload (N, N-m) [If NBodyMod=1, one size 6*NBody x 1 vector; if NBodyMod>1, NBody size 6 x 1 vectors]\n'
            else:
                ln = ln + '\n'
            f.write(ln)
        for j in range(6):
            try:
                ln = " ".join(['{:14}'.format(i) for i in self.fst_vt['HydroDyn']['AddCLin'][j,:]])
            except:
                ln = " ".join(['{:14}'.format(i) for i in self.fst_vt['HydroDyn']['AddCLin'][j]])
            if j == 0:
                ln = ln + "   AddCLin  - Additional linear stiffness (N/m, N/rad, N-m/m, N-m/rad) [If NBodyMod=1, one size 6*NBody x 6*NBody matrix; if NBodyMod>1, NBody size 6 x 6 matrices]\n"
            else:
                ln = ln  + "\n"
            f.write(ln)
        for j in range(6):
            try:
                ln = " ".join(['{:14}'.format(i) for i in self.fst_vt['HydroDyn']['AddBLin'][j,:]])
            except:
                ln = " ".join(['{:14}'.format(i) for i in self.fst_vt['HydroDyn']['AddBLin'][j]])
            if j == 0:
                ln = ln + "   AddBLin  - Additional linear damping(N/(m/s), N/(rad/s), N-m/(m/s), N-m/(rad/s)) [If NBodyMod=1, one size 6*NBody x 6*NBody matrix; if NBodyMod>1, NBody size 6 x 6 matrices]\n"
            else:
                ln = ln  + "\n"
            f.write(ln)
        for j in range(6):
            try:
                ln = " ".join(['{:14}'.format(i) for i in self.fst_vt['HydroDyn']['AddBQuad'][j,:]])
            except:
                ln = " ".join(['{:14}'.format(i) for i in self.fst_vt['HydroDyn']['AddBQuad'][j]])
            if j == 0:
                ln = ln + "   AddBQuad - Additional quadratic drag(N/(m/s)^2, N/(rad/s)^2, N-m(m/s)^2, N-m/(rad/s)^2) [If NBodyMod=1, one size 6*NBody x 6*NBody matrix; if NBodyMod>1, NBody size 6 x 6 matrices]\n"
            else:
                ln = ln  + "\n"
            f.write(ln)

        f.write('---------------------- STRIP THEORY OPTIONS --------------------------------------\n')
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['HydroDyn']['WaveDisp'], 'WaveDisp', '- Method of computing Wave Kinematics {0: use undisplaced position, 1: use displaced position) } (switch)\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['HydroDyn']['AMMod'], 'AMMod', '- Method of computing distributed added-mass force. (0: Only and always on nodes below SWL at the undisplaced position. 2: Up to the instantaneous free surface) [overwrite to 0 when WaveMod = 0 or 6 or when WaveStMod = 0 in SeaState]\n'))
        
        f.write('---------------------- AXIAL COEFFICIENTS --------------------------------------\n')
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['HydroDyn']['NAxCoef'], 'NAxCoef', '- Number of axial coefficients (-)\n'))
        f.write(" ".join(['{:^11s}'.format(i) for i in ['AxCoefID', 'AxCd', 'AxCa', 'AxCp', 'AxFDMod', 'AxVnCOff', 'AxFDLoFSc']])+'\n')
        f.write(" ".join(['{:^11s}'.format(i) for i in ['(-)']*7])+'\n')
        for i in range(self.fst_vt['HydroDyn']['NAxCoef']):
            ln = []
            ln.append('{:^11d}'.format(self.fst_vt['HydroDyn']['AxCoefID'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['AxCd'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['AxCa'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['AxCp'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['AxFDMod'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['AxVnCOff'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['AxFDLoFSc'][i]))
            f.write(" ".join(ln) + '\n')
        f.write('---------------------- MEMBER JOINTS -------------------------------------------\n')
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['HydroDyn']['NJoints'], 'NJoints', '- Number of joints (-)   [must be exactly 0 or at least 2]\n'))
        f.write(" ".join(['{:^11s}'.format(i) for i in ['JointID', 'Jointxi', 'Jointyi', 'Jointzi', 'JointAxID', 'JointOvrlp']])+'  ! [JointOvrlp= 0: do nothing at joint, 1: eliminate overlaps by calculating super member]\n')
        f.write(" ".join(['{:^11s}'.format(i) for i in ['(-)', '(m)', '(m)', '(m)', '(-)', '(switch)']])+'\n')
        for i in range(self.fst_vt['HydroDyn']['NJoints']):
            ln = []
            ln.append('{:^11d}'.format(self.fst_vt['HydroDyn']['JointID'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['Jointxi'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['Jointyi'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['Jointzi'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['HydroDyn']['JointAxID'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['HydroDyn']['JointOvrlp'][i]))
            f.write(" ".join(ln) + '\n')
        f.write('---------------------- MEMBER CROSS-SECTION PROPERTIES -------------------------\n')
        f.write('{:<11d} {:<11} {:}'.format(self.fst_vt['HydroDyn']['NPropSets'], 'NPropSets', '- Number of member property sets (-)\n'))
        f.write(" ".join(['{:^11s}'.format(i) for i in ['PropSetID', 'PropD', 'PropThck']])+'\n')
        f.write(" ".join(['{:^11s}'.format(i) for i in ['(-)', '(m)', '(m)']])+'\n')
        for i in range(self.fst_vt['HydroDyn']['NPropSets']):
            ln = []
            ln.append('{:^11d}'.format(self.fst_vt['HydroDyn']['PropSetID'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['PropD'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['PropThck'][i]))
            f.write(" ".join(ln) + '\n')
        f.write('---------------------- SIMPLE HYDRODYNAMIC COEFFICIENTS (model 1) --------------\n')
        f.write(" ".join(['{:^11s}'.format(i) for i in ['SimplCd', 'SimplCdMG', 'SimplCa', 'SimplCaMG', 'SimplCp', 'SimplCpMG', 'SimplAxCd', 'SimplAxCdMG', 'SimplAxCa', 'SimplAxCaMG', 'SimplAxCp', 'SimplAxCpMG', 'SimplCb', 'SimplCbMG']])+'\n')
        f.write(" ".join(['{:^11s}'.format(i) for i in ['(-)']*14])+'\n')
        ln = []
        ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['SimplCd']))
        ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['SimplCdMG']))
        ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['SimplCa']))
        ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['SimplCaMG']))
        ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['SimplCp']))
        ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['SimplCpMG']))
        ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['SimplAxCd']))
        ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['SimplAxCdMG']))
        ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['SimplAxCa']))
        ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['SimplAxCaMG']))
        ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['SimplAxCp']))
        ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['SimplAxCpMG']))
        ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['SimplCb']))
        ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['SimplCbMG']))
        f.write(" ".join(ln) + '\n')
        f.write('---------------------- DEPTH-BASED HYDRODYNAMIC COEFFICIENTS (model 2) ---------\n')        
        f.write('{:<11d} {:<11} {:}'.format(self.fst_vt['HydroDyn']['NCoefDpth'], 'NCoefDpth', '- Number of depth-dependent coefficients (-)\n'))
        f.write(" ".join(['{:^11s}'.format(i) for i in ['Dpth', 'DpthCd', 'DpthCdMG', 'DpthCa', 'DpthCaMG', 'DpthCp', 'DpthCpMG', 'DpthAxCd', 'DpthAxCdMG', 'DpthAxCa', 'DpthAxCaMG', 'DpthAxCp', 'DpthAxCpMG', 'DpthCb', 'DpthCbMG']])+'\n')
        f.write(" ".join(['{:^11s}'.format(i) for i in ['(m)', '(-)', '(-)', '(-)', '(-)', '(-)', '(-)', '(-)', '(-)', '(-)', '(-)', '(-)', '(-)']])+'\n')
        for i in range(self.fst_vt['HydroDyn']['NCoefDpth']):
            ln = []
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['Dpth'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['DpthCd'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['DpthCdMG'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['DpthCa'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['DpthCaMG'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['DpthCp'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['DpthCpMG'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['DpthAxCd'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['DpthAxCdMG'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['DpthAxCa'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['DpthAxCaMG'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['DpthAxCp'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['DpthAxCpMG'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['DpthCb'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['DpthCbMG'][i]))
            f.write(" ".join(ln) + '\n')
        f.write('---------------------- MEMBER-BASED HYDRODYNAMIC COEFFICIENTS (model 3) --------\n')
        f.write('{:<11d} {:<11} {:}'.format(self.fst_vt['HydroDyn']['NCoefMembers'], 'NCoefMembers', '- Number of member-based coefficients (-)\n'))
        mem_coeff_names = ['MemberID_HydC', 'MemberCd1', 'MemberCd2', 'MemberCdMG1', 'MemberCdMG2', 'MemberCa1', 'MemberCa2', 'MemberCaMG1', 'MemberCaMG2', 'MemberCp1', 'MemberCp2', 'MemberCpMG1', 'MemberCpMG2', 'MemberAxCd1', 'MemberAxCd2', 'MemberAxCdMG1', 'MemberAxCdMG2', 'MemberAxCa1', 'MemberAxCa2', 'MemberAxCaMG1', 'MemberAxCaMG2', 'MemberAxCp1', 'MemberAxCp2', 'MemberAxCpMG1', 'MemberAxCpMG2','MemberCb1','MemberCb2','MemberCbMG1','MemberCbMG2']
        f.write(" ".join(['{:^11s}'.format(i) for i in mem_coeff_names])+'\n')
        f.write(" ".join(['{:^11s}'.format(i) for i in ['(-)']*len(mem_coeff_names)])+'\n')
        for i in range(self.fst_vt['HydroDyn']['NCoefMembers']):
            ln = []
            ln.append('{:^11d}'.format(self.fst_vt['HydroDyn']['MemberID_HydC'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MemberCd1'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MemberCd2'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MemberCdMG1'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MemberCdMG2'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MemberCa1'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MemberCa2'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MemberCaMG1'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MemberCaMG2'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MemberCp1'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MemberCp2'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MemberCpMG1'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MemberCpMG2'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MemberAxCd1'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MemberAxCd2'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MemberAxCdMG1'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MemberAxCdMG2'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MemberAxCa1'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MemberAxCa2'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MemberAxCaMG1'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MemberAxCaMG2'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MemberAxCp1'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MemberAxCp2'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MemberAxCpMG1'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MemberAxCpMG2'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MemberCb1'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MemberCb2'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MemberCbMG1'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MemberCbMG2'][i]))
            f.write(" ".join(ln) + '\n')
        f.write('-------------------- MEMBERS -------------------------------------------------\n')
        f.write('{:<11d} {:<11} {:}'.format(self.fst_vt['HydroDyn']['NMembers'], 'NMembers', '- Number of members (-)\n'))
        f.write(" ".join(['{:^11s}'.format(i) for i in ['MemberID', 'MJointID1', 'MJointID2', 'MPropSetID1', 'MPropSetID2', 'MDivSize', 'MCoefMod', 'MHstLMod', 'PropPot']])+' !  [MCoefMod=1: use simple coeff table, 2: use depth-based coeff table, 3: use member-based coeff table] [ PropPot/=0 if member is modeled with potential-flow theory]\n')
        f.write(" ".join(['{:^11s}'.format(i) for i in ['(-)', '(-)', '(-)', '(-)', '(-)', '(m)', '(switch)', '(switch)', '(flag)']])+'\n')
        for i in range(self.fst_vt['HydroDyn']['NMembers']):
            ln = []
            ln.append('{:^11d}'.format(self.fst_vt['HydroDyn']['MemberID'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['HydroDyn']['MJointID1'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['HydroDyn']['MJointID2'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['HydroDyn']['MPropSetID1'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['HydroDyn']['MPropSetID2'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MDivSize'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['HydroDyn']['MCoefMod'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['HydroDyn']['MHstLMod'][i]))
            ln.append('{!s:^11}'.format(self.fst_vt['HydroDyn']['PropPot'][i]))
            f.write(" ".join(ln) + '\n')
        f.write("---------------------- FILLED MEMBERS ------------------------------------------\n")
        f.write('{:<11d} {:<11} {:}'.format(self.fst_vt['HydroDyn']['NFillGroups'], 'NFillGroups', '- Number of filled member groups (-) [If FillDens = DEFAULT, then FillDens = WtrDens; FillFSLoc is related to MSL2SWL]\n'))
        f.write(" ".join(['{:^11s}'.format(i) for i in ['FillNumM', 'FillMList', 'FillFSLoc', 'FillDens']])+'\n')
        f.write(" ".join(['{:^11s}'.format(i) for i in ['(-)', '(-)', '(m)', '(kg/m^3)']])+'\n')
        for i in range(self.fst_vt['HydroDyn']['NFillGroups']):
            ln = []
            ln.append('{:^11d}'.format(self.fst_vt['HydroDyn']['FillNumM'][i]))
            ln.append(" ".join(['%d'%j for j in self.fst_vt['HydroDyn']['FillMList'][i]]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['FillFSLoc'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['FillDens'][i]))
            f.write(" ".join(ln) + '\n')
        f.write("---------------------- MARINE GROWTH -------------------------------------------\n")
        f.write('{:<11d} {:<11} {:}'.format(self.fst_vt['HydroDyn']['NMGDepths'], 'NMGDepths', '- Number of marine-growth depths specified (-)\n'))
        f.write(" ".join(['{:^11s}'.format(i) for i in ['MGDpth', 'MGThck', 'MGDens']])+'\n')
        f.write(" ".join(['{:^11s}'.format(i) for i in ['(m)', '(m)', '(kg/m^3)']])+'\n')
        for i in range(self.fst_vt['HydroDyn']['NMGDepths']):
            ln = []
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MGDpth'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MGThck'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['MGDens'][i]))
            f.write(" ".join(ln) + '\n')
        f.write("---------------------- MEMBER OUTPUT LIST --------------------------------------\n")
        f.write('{:<11d} {:<11} {:}'.format(self.fst_vt['HydroDyn']['NMOutputs'], 'NMOutputs', '- Number of member outputs (-) [must be <=99]\n'))
        f.write(" ".join(['{:^11s}'.format(i) for i in ['MemberID_out', 'NOutLoc', 'NodeLocs  [NOutLoc < 10; node locations are normalized distance from the start of the member, and must be >=0 and <= 1] [unused if NMOutputs=0]']])+'\n')
        f.write(" ".join(['{:^11s}'.format(i) for i in ['(-)']*3])+'\n')
        for i in range(self.fst_vt['HydroDyn']['NMOutputs']):
            ln = []
            ln.append('{:^11d}'.format(self.fst_vt['HydroDyn']['MemberID_out'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['HydroDyn']['NOutLoc'][i]))
            ln.append('{:^11}'.format(self.fst_vt['HydroDyn']['NodeLocs'][i]))
            f.write(" ".join(ln) + '\n')
        f.write("---------------------- JOINT OUTPUT LIST ---------------------------------------\n")
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['HydroDyn']['NJOutputs'], 'NJOutputs', '- Number of joint outputs [Must be < 10]\n'))
        f.write('{:<22} {:<11} {:}'.format(" ".join(["%d"%i for i in self.fst_vt['HydroDyn']['JOutLst']]), 'JOutLst', '- List of JointIDs which are to be output (-)[unused if NJOutputs=0]\n'))
        f.write("---------------------- OUTPUT --------------------------------------------------\n")
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['HydroDyn']['HDSum'], 'HDSum', '- Output a summary file [flag]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['HydroDyn']['OutAll'], 'OutAll', '- Output all user-specified member and joint loads (only at each member end, not interior locations) [flag]\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['HydroDyn']['OutSwtch'], 'OutSwtch', '- Output requested channels to: [1=Hydrodyn.out, 2=GlueCode.out, 3=both files]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['HydroDyn']['OutFmt'], 'OutFmt', '- Output format for numerical results (quoted string) [not checked for validity!]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['HydroDyn']['OutSFmt'], 'OutSFmt', '- Output format for header strings (quoted string) [not checked for validity!]\n'))
        f.write('---------------------- OUTPUT CHANNELS -----------------------------------------\n')
        outlist = self.get_outlist(self.fst_vt['outlist'], ['HydroDyn'])
        for channel_list in outlist:
            for i in range(len(channel_list)):
                f.write('"' + channel_list[i] + '"\n')
            
        f.write('END of output channels and end of file. (the word "END" must appear in the first 3 columns of this line)\n')
        
        f.close()

    def write_SeaState(self):

        # Generate SeaState input file
        self.fst_vt['Fst']['SeaStFile'] = self.FAST_namingOut + '_SeaState.dat'
        hd_file = os.path.join(self.FAST_runDirectory, self.fst_vt['Fst']['SeaStFile'])
        f = open(hd_file, 'w')

        f.write('------- SeaState Input File --------------------------------------------\n')
        f.write('Generated with OpenFAST_IO\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['Echo'], 'Echo', '- Echo the input file data (flag)\n'))

        f.write('---------------------- ENVIRONMENTAL CONDITIONS --------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['WtrDens'],'WtrDens', '- Water density (kg/m^3)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['WtrDpth'],'WtrDpth', '- Water depth (meters) relative to MSL\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['MSL2SWL'],'MSL2SWL', '- Offset between still-water level and mean sea level (meters) [positive upward; unused when WaveMod = 6; must be zero if PotMod=1 or 2]\n'))


        f.write('---------------------- SPATIAL DISCRETIZATION ---------------------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['X_HalfWidth'], 'X_HalfWidth', '– Half-width of the domain in the X direction (m) [>0, NOTE: X[nX] = nX*dX, where nX = {-NX+1,-NX+2,…,NX-1} and dX = X_HalfWidth/(NX-1)]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['Y_HalfWidth'], 'Y_HalfWidth', '– Half-width of the domain in the Y direction (m) [>0, NOTE: Y[nY] = nY*dY, where nY = {-NY+1,-NY+2,…,NY-1} and dY = Y_HalfWidth/(NY-1)]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['Z_Depth'], 'Z_Depth', '– Depth of the domain the Z direction (m) relative to SWL [0 < Z_Depth <= WtrDpth+MSL2SWL; "default": Z_Depth = WtrDpth+MSL2SWL; Z[nZ] = ( COS( nZ*dthetaZ ) – 1 )*Z_Depth, where nZ = {0,1,…NZ-1} and dthetaZ = pi/( 2*(NZ-1) )]\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SeaState']['NX'], 'NX', '– Number of nodes in half of the X-direction domain (-) [>=2]\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SeaState']['NY'], 'NY', '– Number of nodes in half of the Y-direction domain (-) [>=2]\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SeaState']['NZ'], 'NZ', '– Number of nodes in the Z direction (-) [>=2]\n'))

        f.write('---------------------- WAVES ---------------------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['WaveMod'], 'WaveMod', '- Incident wave kinematics model {0: none=still water, 1: regular (periodic), 1P#: regular with user-specified phase, 2: JONSWAP/Pierson-Moskowitz spectrum (irregular), 3: White noise spectrum (irregular), 4: user-defined spectrum from routine UserWaveSpctrm (irregular), 5: Externally generated wave-elevation time series, 6: Externally generated full wave-kinematics time series [option 6 is invalid for PotMod/=0]} (switch)\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SeaState']['WaveStMod'], 'WaveStMod', '- Model for stretching incident wave kinematics to instantaneous free surface {0: none=no stretching, 1: vertical stretching, 2: extrapolation stretching, 3: Wheeler stretching} (switch) [unused when WaveMod=0 or when PotMod/=0]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['WaveTMax'], 'WaveTMax', '- Analysis time for incident wave calculations (sec) [unused when WaveMod=0; determines WaveDOmega=2Pi/WaveTMax in the IFFT]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['WaveDT'], 'WaveDT', '- Time step for incident wave calculations     (sec) [unused when WaveMod=0 or 7; 0.1<=WaveDT<=1.0 recommended; determines WaveOmegaMax=Pi/WaveDT in the IFFT]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['WaveHs'], 'WaveHs', '- Significant wave height of incident waves (meters) [used only when WaveMod=1, 2, or 3]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['WaveTp'], 'WaveTp', '- Peak-spectral period of incident waves       (sec) [used only when WaveMod=1 or 2]\n'))
        if isinstance(self.fst_vt['SeaState']['WavePkShp'], float):
            if self.fst_vt['SeaState']['WavePkShp'] == 0.:
                WavePkShp = 'Default'
            else: 
                WavePkShp = self.fst_vt['SeaState']['WavePkShp']
        else:
            WavePkShp = self.fst_vt['SeaState']['WavePkShp']
        f.write('{:<22} {:<11} {:}'.format(WavePkShp, 'WavePkShp', '- Peak-shape parameter of incident wave spectrum (-) or DEFAULT (string) [used only when WaveMod=2; use 1.0 for Pierson-Moskowitz]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['WvLowCOff'], 'WvLowCOff', '- Low  cut-off frequency or lower frequency limit of the wave spectrum beyond which the wave spectrum is zeroed (rad/s) [unused when WaveMod=0, 1, or 6]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['WvHiCOff'], 'WvHiCOff', '- High cut-off frequency or upper frequency limit of the wave spectrum beyond which the wave spectrum is zeroed (rad/s) [unused when WaveMod=0, 1, or 6]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['WaveDir'], 'WaveDir', '- Incident wave propagation heading direction                         (degrees) [unused when WaveMod=0, 6 or 7]\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SeaState']['WaveDirMod'], 'WaveDirMod', '- Directional spreading function {0: none, 1: COS2S}                  (-)       [only used when WaveMod=2,3, or 4]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['WaveDirSpread'], 'WaveDirSpread', '- Wave direction spreading coefficient ( > 0 )                        (-)       [only used when WaveMod=2,3, or 4 and WaveDirMod=1]\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SeaState']['WaveNDir'], 'WaveNDir', '- Number of wave directions                                           (-)       [only used when WaveMod=2,3, or 4 and WaveDirMod=1; odd number only]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['WaveDirRange'], 'WaveDirRange', '- Range of wave directions (full range: WaveDir +/- 1/2*WaveDirRange) (degrees) [only used when WaveMod=2,3,or 4 and WaveDirMod=1]\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SeaState']['WaveSeed1'], 'WaveSeed(1)', '- First  random seed of incident waves [-2147483648 to 2147483647]    (-)       [unused when WaveMod=0, 5, or 6]\n'))
        
        try:
            seed2 = int(self.fst_vt['SeaState']['WaveSeed2'])
            f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SeaState']['WaveSeed2'], 'WaveSeed(2)', '- Second random seed of incident waves [-2147483648 to 2147483647]    (-)       [unused when WaveMod=0, 5, or 6]\n'))
        except ValueError:
            f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['WaveSeed2'], 'WaveSeed(2)', '- Second random seed of incident waves [-2147483648 to 2147483647]    (-) for intrinsic pRNG, or an alternative pRNG: "RanLux" (-) [unused when WaveMod=0, 5, or 6]\n'))
            
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['WaveNDAmp'], 'WaveNDAmp', '- Flag for normally distributed amplitudes                            (flag)    [only used when WaveMod=2, 3, or 4]\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['SeaState']['WvKinFile']+'"', 'WvKinFile', '- Root name of externally generated wave data file(s)        (quoted string)    [used only when WaveMod=5, 6 or 7]\n'))
        f.write('---------------------- 2ND-ORDER WAVES ----------------------------------------- [unused with WaveMod=0 or 6]\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['WvDiffQTF'], 'WvDiffQTF', '- Full difference-frequency 2nd-order wave kinematics (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['WvSumQTF'], 'WvSumQTF', '- Full summation-frequency  2nd-order wave kinematics (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['WvLowCOffD'], 'WvLowCOffD', '- Low  frequency cutoff used in the difference-frequencies (rad/s) [Only used with a difference-frequency method]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['WvHiCOffD'], 'WvHiCOffD', '- High frequency cutoff used in the difference-frequencies (rad/s) [Only used with a difference-frequency method]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['WvLowCOffS'], 'WvLowCOffS', '- Low  frequency cutoff used in the summation-frequencies  (rad/s) [Only used with a summation-frequency  method]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['WvHiCOffS'], 'WvHiCOffS', '- High frequency cutoff used in the summation-frequencies  (rad/s) [Only used with a summation-frequency  method]\n'))
        f.write('---------------------- CONSTRAINED WAVES ----------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['ConstWaveMod'], 'ConstWaveMod', '- Constrained wave model: 0=none; 1=Constrained wave with specified crest elevation, alpha; 2=Constrained wave with guaranteed peak-to-trough crest height, HCrest (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['CrestHmax'], 'CrestHmax', '- Crest height (2*alpha for ConstWaveMod=1 or HCrest for ConstWaveMod=2), must be larger than WaveHs (m) [unused when ConstWaveMod=0]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['CrestTime'], 'CrestTime', '- Time at which the crest appears (s) [unused when ConstWaveMod=0]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['CrestXi'], 'CrestXi', '- X-position of the crest.   (m)     [unused when ConstWaveMod=0]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['CrestYi'], 'CrestYi', '- Y-position of the crest.   (m)     [unused when ConstWaveMod=0]\n'))
        f.write('---------------------- CURRENT ------------------------------------------------- [unused with WaveMod=6]\n')
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SeaState']['CurrMod'], 'CurrMod', '- Current profile model {0: none=no current, 1: standard, 2: user-defined from routine UserCurrent} (switch)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['CurrSSV0'], 'CurrSSV0', '- Sub-surface current velocity at still water level  (m/s) [used only when CurrMod=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['CurrSSDir'], 'CurrSSDir', '- Sub-surface current heading direction (degrees) or DEFAULT (string) [used only when CurrMod=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['CurrNSRef'], 'CurrNSRef', '- Near-surface current reference depth            (meters) [used only when CurrMod=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['CurrNSV0'], 'CurrNSV0', '- Near-surface current velocity at still water level (m/s) [used only when CurrMod=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['CurrNSDir'], 'CurrNSDir', '- Near-surface current heading direction         (degrees) [used only when CurrMod=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['CurrDIV'], 'CurrDIV', '- Depth-independent current velocity                 (m/s) [used only when CurrMod=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['CurrDIDir'], 'CurrDIDir', '- Depth-independent current heading direction    (degrees) [used only when CurrMod=1]\n'))
        f.write('---------------------- MacCamy-Fuchs Diffraction Model -------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['MCFD'],'MCFD', '- MacCamy-Fuchs member radius (ignored if radius <= 0) [must be 0 when WaveMod 0 or 6] \n'))
        f.write('---------------------- OUTPUT --------------------------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['SeaStSum'], 'SeaStSum', '- Output a summary file [flag]\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SeaState']['OutSwtch'], 'OutSwtch','- Output requested channels to: [1=SeaState.out, 2=GlueCode.out, 3=both files]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['OutFmt'], 'OutFmt','- Output format for numerical results (quoted string) [not checked for validity!]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['SeaState']['OutSFmt'], 'OutSFmt','- Output format for header strings (quoted string) [not checked for validity!]\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SeaState']['NWaveElev'], 'NWaveElev','- Number of points where the incident wave elevations can be computed (-)       [maximum of 9 output locations]\n'))
        f.write('{:<22} {:<11} {:}'.format(", ".join([f'{float(val):f}' for val in self.fst_vt['SeaState']['WaveElevxi']]), 'WaveElevxi', '- List of xi-coordinates for points where the incident wave elevations can be output (meters) [NWaveElev points, separated by commas or white space; usused if NWaveElev = 0]\n'))
        f.write('{:<22} {:<11} {:}'.format(", ".join([f'{float(val):f}' for val in self.fst_vt['SeaState']['WaveElevyi']]), 'WaveElevyi', '- List of yi-coordinates for points where the incident wave elevations can be output (meters) [NWaveElev points, separated by commas or white space; usused if NWaveElev = 0]\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SeaState']['NWaveKin'], 'NWaveKin','- Number of points where the wave kinematics can be output (-)       [maximum of 9 output locations]\n'))
        
        if self.fst_vt['SeaState']['NWaveKin'] > 0 :
            f.write('{:<22} {:<11} {:}'.format(", ".join([f'{val:f}' for val in self.fst_vt['SeaState']['WaveKinxi']]), 'WaveKinxi', '- List of xi-coordinates for points where the wave kinematics can be output (meters) [NWaveKin points, separated by commas or white space; usused if NWaveKin = 0]\n'))
            f.write('{:<22} {:<11} {:}'.format(", ".join([f'{val:f}' for val in self.fst_vt['SeaState']['WaveKinyi']]), 'WaveKinyi', '- List of yi-coordinates for points where the wave kinematics can be output (meters) [NWaveKin points, separated by commas or white space; usused if NWaveKin = 0]\n'))
            f.write('{:<22} {:<11} {:}'.format(", ".join([f'{val:f}' for val in self.fst_vt['SeaState']['WaveKinzi']]), 'WaveKinzi', '- List of zi-coordinates for points where the wave kinematics can be output (meters) [NWaveKin points, separated by commas or white space; usused if NWaveKin = 0]\n'))
        else:
            f.write('{:<11} {:}'.format('WaveKinxi', '- List of xi-coordinates for points where the wave kinematics can be output (meters) [NWaveKin points, separated by commas or white space; usused if NWaveKin = 0]\n'))
            f.write('{:<11} {:}'.format('WaveKinyi', '- List of yi-coordinates for points where the wave kinematics can be output (meters) [NWaveKin points, separated by commas or white space; usused if NWaveKin = 0]\n'))
            f.write('{:<11} {:}'.format('WaveKinzi', '- List of zi-coordinates for points where the wave kinematics can be output (meters) [NWaveKin points, separated by commas or white space; usused if NWaveKin = 0]\n'))

        f.write('---------------------- OUTPUT CHANNELS -----------------------------------------\n')
        outlist = self.get_outlist(self.fst_vt['outlist'], ['SeaState'])
        for channel_list in outlist:
            for i in range(len(channel_list)):
                f.write('"' + channel_list[i] + '"\n')
            
        f.write('END of output channels and end of file. (the word "END" must appear in the first 3 columns of this line)\n')
        
        f.close()

    def write_SubDyn(self):
        # Generate SubDyn input file
        self.fst_vt['Fst']['SubFile'] = self.FAST_namingOut + '_SubDyn.dat'
        sd_file = os.path.join(self.FAST_runDirectory, self.fst_vt['Fst']['SubFile'])
        f = open(sd_file, 'w')

        f.write('----------- SubDyn MultiMember Support Structure Input File ------------\n')
        f.write('Generated with OpenFAST_IO\n')
        f.write('-------------------------- SIMULATION CONTROL  ---------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['SubDyn']['Echo'], 'Echo', '- Echo input data to "<rootname>.SD.ech" (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SubDyn']['SDdeltaT'], 'SDdeltaT', '- Local Integration Step. If "default", the glue-code integration step will be used.\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SubDyn']['IntMethod'], 'IntMethod', '- Integration Method [1/2/3/4 = RK4/AB4/ABM4/AM2].\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['SubDyn']['SttcSolve'], 'SttcSolve', '- Solve dynamics about static equilibrium point\n'))
        # f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['SubDyn']['GuyanLoadCorrection'], 'GuyanLoadCorrection', '- Include extra moment from lever arm at interface and rotate FEM for floating.\n'))
        f.write('-------------------- FEA and CRAIG-BAMPTON PARAMETERS---------------------------\n')
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SubDyn']['FEMMod'], 'FEMMod', '- FEM switch: element model in the FEM. [1= Euler-Bernoulli(E-B);  2=Tapered E-B (unavailable);  3= 2-node Timoshenko;  4= 2-node tapered Timoshenko (unavailable)]\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SubDyn']['NDiv'], 'NDiv', '- Number of sub-elements per member\n'))
        # f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['SubDyn']['CBMod'], 'CBMod', '- [T/F] If True perform C-B reduction, else full FEM dofs will be retained. If True, select Nmodes to retain in C-B reduced system.\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SubDyn']['Nmodes'], 'Nmodes', '- Number of internal modes to retain. If Nmodes=0 --> Guyan Reduction. If Nmodes<0 --> retain all modes.\n'))
        
        JDampings = self.fst_vt['SubDyn']['JDampings']
        if isinstance(JDampings, float):
            f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SubDyn']['JDampings'], 'JDampings', '- Damping Ratios for each retained mode (% of critical) If Nmodes>0, list Nmodes structural damping ratios for each retained mode (% of critical), or a single damping ratio to be applied to all retained modes. (last entered value will be used for all remaining modes).\n'))
        else: # list of floats
            f.write('{:<22} {:<11} {:}'.format(", ".join([f'{d:f}' for d in self.fst_vt['SubDyn']['JDampings']]), 'JDampings', '- Damping Ratios for each retained mode (% of critical) If Nmodes>0, list Nmodes structural damping ratios for each retained mode (% of critical), or a single damping ratio to be applied to all retained modes. (last entered value will be used for all remaining modes).\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SubDyn']['GuyanDampMod'], 'GuyanDampMod', '- Guyan damping {0=none, 1=Rayleigh Damping, 2=user specified 6x6 matrix}.\n'))
        f.write('{:<10} {:<10} {:<11} {:}'.format(self.fst_vt['SubDyn']['RayleighDamp'][0], self.fst_vt['SubDyn']['RayleighDamp'][1], 'RayleighDamp', '- Mass and stiffness proportional damping  coefficients (Rayleigh Damping) [only if GuyanDampMod=1].\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SubDyn']['GuyanDampSize'], 'GuyanDampSize', '- Guyan damping matrix (6x6) [only if GuyanDampMod=2].\n'))
        for j in range(self.fst_vt['SubDyn']['GuyanDampSize']):
            try:
                ln = " ".join(['{:14}'.format(i) for i in self.fst_vt['SubDyn']['GuyanDamp'][j,:]])
            except:
                ln = " ".join(['{:14}'.format(i) for i in self.fst_vt['SubDyn']['GuyanDamp'][j]])
            ln += "\n"
            f.write(ln)
        
        f.write('---- STRUCTURE JOINTS: joints connect structure members (~Hydrodyn Input File)---\n')
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SubDyn']['NJoints'], 'NJoints', '- Number of joints (-)\n'))
        f.write(" ".join(['{:^11s}'.format(i) for i in ['JointID','JointXss','JointYss','JointZss','JointType','JointDirX','JointDirY','JointDirZ','JointStiff']])+' ![Coordinates of Member joints in SS-Coordinate System][JointType={1:cantilever, 2:universal joint, 3:revolute joint, 4:spherical joint}]\n')
        f.write(" ".join(['{:^11s}'.format(i) for i in ['(-)','(m)','(m)','(m)','(-)','(-)','(-)','(-)','(Nm/rad)']])+'\n')
        for i in range(self.fst_vt['SubDyn']['NJoints']):
            ln = []
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['JointID'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['JointXss'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['JointYss'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['JointZss'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['JointType'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['JointDirX'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['JointDirY'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['JointDirZ'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['JointStiff'][i]))
            f.write(" ".join(ln) + '\n')        
        f.write('------------------- BASE REACTION JOINTS: 1/0 for Locked/Free DOF @ each Reaction Node ---------------------\n')
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SubDyn']['NReact'], 'NReact', '- Number of Joints with reaction forces; be sure to remove all rigid motion DOFs of the structure  (else det([K])=[0])\n'))
        f.write(" ".join(['{:^11s}'.format(i) for i in ['RJointID', 'RctTDXss', 'RctTDYss', 'RctTDZss', 'RctRDXss', 'RctRDYss', 'RctRDZss','SSIfile']])+' ! [Global Coordinate System]\n')
        f.write(" ".join(['{:^11s}'.format(i) for i in ['(-)', '(flag)', '(flag)', '(flag)', '(flag)', '(flag)', '(flag)', '(string)']])+'\n')
        for i in range(self.fst_vt['SubDyn']['NReact']):
            ln = []
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['RJointID'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['RctTDXss'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['RctTDYss'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['RctTDZss'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['RctRDXss'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['RctRDYss'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['RctRDZss'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['Rct_SoilFile'][i]))
            f.write(" ".join(ln) + '\n')
        f.write('------- INTERFACE JOINTS: 1/0 for Locked (to the TP)/Free DOF @each Interface Joint (only Locked-to-TP implemented thus far (=rigid TP)) ---------\n')
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SubDyn']['NInterf'], 'NInterf', '- Number of interface joints locked to the Transition Piece (TP):  be sure to remove all rigid motion dofs\n'))
        f.write(" ".join(['{:^11s}'.format(i) for i in ['IJointID', 'ItfTDXss', 'ItfTDYss', 'ItfTDZss', 'ItfRDXss', 'ItfRDYss', 'ItfRDZss']])+'\n')
        f.write(" ".join(['{:^11s}'.format(i) for i in ['(-)', '(flag)', '(flag)', '(flag)', '(flag)', '(flag)', '(flag)']])+'\n')
        for i in range(self.fst_vt['SubDyn']['NInterf']):
            ln = []
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['IJointID'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['ItfTDXss'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['ItfTDYss'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['ItfTDZss'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['ItfRDXss'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['ItfRDYss'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['ItfRDZss'][i]))
            f.write(" ".join(ln) + '\n')
        f.write('----------------------------------- MEMBERS --------------------------------------\n')
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SubDyn']['NMembers'], 'NMembers', '- Number of frame members\n'))
        f.write(" ".join(['{:^11s}'.format(i) for i in ['MemberID', 'MJointID1', 'MJointID2', 'MPropSetID1', 'MPropSetID2', 'MType', 'COSMID']])+' ![MType={1:beam circ., 2:cable, 3:rigid, 4:beam arb.}. COMSID={-1:none}]\n')
        f.write(" ".join(['{:^11s}'.format(i) for i in ['(-)','(-)','(-)','(-)','(-)','(-)','(-)']])+'\n')
        for i in range(self.fst_vt['SubDyn']['NMembers']):
            ln = []
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['MemberID'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['MJointID1'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['MJointID2'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['MPropSetID1'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['MPropSetID2'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['MType'][i]))
            # Need to change M_COSMID None elements to -1
            self.fst_vt['SubDyn']['M_COSMID'][i] = -1 if self.fst_vt['SubDyn']['M_COSMID'][i] is None else self.fst_vt['SubDyn']['M_COSMID'][i]
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['M_COSMID'][i]))
            f.write(" ".join(ln) + '\n')
        f.write('------------------ CIRCULAR BEAM CROSS-SECTION PROPERTIES -----------------------------\n')
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SubDyn']['NPropSets'], 'NPropSets', '- Number of structurally unique cross-sections\n'))
        f.write(" ".join(['{:^11s}'.format(i) for i in ['PropSetID', 'YoungE', 'ShearG1', 'MatDens', 'XsecD', 'XsecT']])+'\n')
        f.write(" ".join(['{:^11s}'.format(i) for i in ['(-)','(N/m2)','(N/m2)','(kg/m3)','(m)','(m)']])+'\n')
        for i in range(self.fst_vt['SubDyn']['NPropSets']):
            ln = []
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['PropSetID1'][i]))
            ln.append('{:^11e}'.format(self.fst_vt['SubDyn']['YoungE1'][i]))
            ln.append('{:^11e}'.format(self.fst_vt['SubDyn']['ShearG1'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['MatDens1'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['XsecD'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['XsecT'][i]))
            f.write(" ".join(ln) + '\n')
        f.write('----------------- ARBITRARY BEAM CROSS-SECTION PROPERTIES -----------------------------\n')
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SubDyn']['NXPropSets'], 'NXPropSets', '- Number of structurally unique non-circular cross-sections (if 0 the following table is ignored)\n'))
        f.write(" ".join(['{:^11s}'.format(i) for i in ['PropSetID', 'YoungE', 'ShearG2', 'MatDens', 'XsecA', 'XsecAsx', 'XsecAsy', 'XsecJxx', 'XsecJyy', 'XsecJ0']])+'\n')
        f.write(" ".join(['{:^11s}'.format(i) for i in ['(-)','(N/m2)','(N/m2)','(kg/m3)','(m2)','(m2)','(m2)','(m4)','(m4)','(m4)']])+'\n')
        for i in range(self.fst_vt['SubDyn']['NXPropSets']):
            ln = []
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['PropSetID2'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['YoungE2'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['ShearG2'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['MatDens2'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['XsecA'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['XsecAsx'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['XsecAsy'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['XsecJxx'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['XsecJyy'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['XsecJ0'][i]))
            f.write(" ".join(ln) + '\n')
        f.write('-------------------------- CABLE PROPERTIES  -------------------------------------\n')
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SubDyn']['NCablePropSets'], 'NCablePropSets', '- Number of cable cable properties\n'))
        f.write(" ".join(['{:^11s}'.format(i) for i in ['PropSetID', 'EA', 'MatDens', 'T0']])+'\n')
        f.write(" ".join(['{:^11s}'.format(i) for i in ['(-)','(N)','(kg/m)','(N)']])+'\n')
        for i in range(self.fst_vt['SubDyn']['NCablePropSets']):
            ln = []
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['CablePropSetID'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['CableEA'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['CableMatDens'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['CableT0'][i]))
            f.write(" ".join(ln) + '\n')
        f.write('----------------------- RIGID LINK PROPERTIES ------------------------------------\n')
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SubDyn']['NRigidPropSets'], 'NRigidPropSets', '- Number of rigid link properties\n'))
        f.write(" ".join(['{:^11s}'.format(i) for i in ['PropSetID', 'MatDens']])+'\n')
        f.write(" ".join(['{:^11s}'.format(i) for i in ['(-)','(kg/m)']])+'\n')
        for i in range(self.fst_vt['SubDyn']['NRigidPropSets']):
            ln = []
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['RigidPropSetID'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['RigidMatDens'][i]))
            f.write(" ".join(ln) + '\n')
        f.write('----------------------- SPRING ELEMENT PROPERTIES ------------------------------------\n')
        spring_list = ['k11','k12','k13','k14','k15','k16',
                       'k22','k23','k24','k25','k26',
                       'k33','k34','k35','k36',
                       'k44','k45','k46',
                       'k55','k56',
                       'k66']
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SubDyn']['NSpringPropSets'], 'NSpringPropSets', '- Number of spring properties\n'))
        f.write("PropSetID   k11     k12     k13     k14     k15     k16     k22     k23     k24     k25     k26     k33     k34     k35     k36     k44      k45      k46      k55      k56      k66    \n")
        f.write("  (-)      (N/m)   (N/m)   (N/m)  (N/rad) (N/rad) (N/rad)  (N/m)   (N/m)  (N/rad) (N/rad) (N/rad)  (N/m)  (N/rad) (N/rad) (N/rad) (Nm/rad) (Nm/rad) (Nm/rad) (Nm/rad) (Nm/rad) (Nm/rad) \n")
        for i in range(self.fst_vt['SubDyn']['NSpringPropSets']):
            ln = []
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['SpringPropSetID'][i]))
            for sl in spring_list:
                ln.append('{:^11}'.format(self.fst_vt['SubDyn'][sl][i]))
            f.write(" ".join(ln) + '\n')
        f.write('---------------------- MEMBER COSINE MATRICES COSM(i,j) ------------------------\n')
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SubDyn']['NCOSMs'], 'NCOSMs', '- Number of unique cosine matrices (i.e., of unique member alignments including principal axis rotations); ignored if NXPropSets=0   or 9999 in any element below\n'))
        f.write(" ".join(['{:^11s}'.format(i) for i in ['COSMID', 'COSM11', 'COSM12', 'COSM13', 'COSM21', 'COSM22', 'COSM23', 'COSM31', 'COSM32', 'COSM33']])+'\n')
        f.write(" ".join(['{:^11s}'.format(i) for i in ['(-)','(-)','(-)','(-)','(-)','(-)','(-)','(-)','(-)','(-)']])+'\n')
        for i in range(self.fst_vt['SubDyn']['NCOSMs']):
            ln = []
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['COSMID'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['COSM11'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['COSM12'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['COSM13'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['COSM21'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['COSM22'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['COSM23'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['COSM31'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['COSM32'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['COSM33'][i]))
            f.write(" ".join(ln) + '\n')
        f.write('------------------------ JOINT ADDITIONAL CONCENTRATED MASSES--------------------------\n')
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SubDyn']['NCmass'], 'NCmass', '- Number of joints with concentrated masses; Global Coordinate System\n'))
        f.write(" ".join(['{:^11s}'.format(i) for i in ['CMJointID','JMass','JMXX','JMYY','JMZZ','JMXY','JMXZ','JMYZ','MCGX','MCGY','MCGZ']])+'\n')
        f.write(" ".join(['{:^11s}'.format(i) for i in ['(-)','(kg)','(kg*m^2)','(kg*m^2)','(kg*m^2)','(kg*m^2)','(kg*m^2)','(kg*m^2)','(m)','(m)','(m)']])+'\n')
        for i in range(self.fst_vt['SubDyn']['NCmass']):
            ln = []
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['CMJointID'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['JMass'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['JMXX'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['JMYY'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['JMZZ'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['JMXY'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['JMXZ'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['JMYZ'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['MCGX'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['MCGY'][i]))
            ln.append('{:^11}'.format(self.fst_vt['SubDyn']['MCGZ'][i]))
            f.write(" ".join(ln) + '\n')
        f.write('---------------------------- OUTPUT: SUMMARY & OUTFILE ------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['SubDyn']['SumPrint'], 'SumPrint', '- Output a Summary File (flag).It contains: matrices K,M  and C-B reduced M_BB, M-BM, K_BB, K_MM(OMG^2), PHI_R, PHI_L. It can also contain COSMs if requested.\n'))
        if 'OutCBModes' in self.fst_vt['SubDyn']:
            f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SubDyn']['OutCBModes'], 'OutCBModes', '- Output Guyan and Craig-Bampton modes {0 No output, 1 JSON output}, (flag)\n'))
        if 'OutFEMModes' in self.fst_vt['SubDyn']:
            f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SubDyn']['OutFEMModes'], 'OutFEMModes', '-  Output first 30 FEM modes {0 No output, 1 JSON output} (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['SubDyn']['OutCOSM'], 'OutCOSM', '- Output cosine matrices with the selected output member forces (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['SubDyn']['OutAll'], 'OutAll', "- [T/F] Output all members' end forces\n"))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SubDyn']['OutSwtch'], 'OutSwtch', '- [1/2/3] Output requested channels to: 1=<rootname>.SD.out;  2=<rootname>.out (generated by FAST);  3=both files.\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['SubDyn']['TabDelim'], 'TabDelim', '- Generate a tab-delimited output in the <rootname>.SD.out file\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SubDyn']['OutDec'], 'OutDec', '- Decimation of output in the <rootname>.SD.out file\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SubDyn']['OutFmt'], 'OutFmt', '- Output format for numerical results in the <rootname>.SD.out file\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['SubDyn']['OutSFmt'], 'OutSFmt', '- Output format for header strings in the <rootname>.SD.out file\n'))
        f.write('------------------------- MEMBER OUTPUT LIST ------------------------------------------\n')
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['SubDyn']['NMOutputs'], 'NMOutputs', '- Number of members whose forces/displacements/velocities/accelerations will be output (-) [Must be <= 99].\n'))
        f.write(" ".join(['{:^11s}'.format(i) for i in ['MemberID', 'NOutCnt', 'NodeCnt']])+' ! [NOutCnt=how many nodes to get output for [< 10]; NodeCnt are local ordinal numbers from the start of the member, and must be >=1 and <= NDiv+1] If NMOutputs=0 leave blank as well.\n')
        f.write(" ".join(['{:^11s}'.format(i) for i in ['(-)','(-)','(-)']])+'\n')
        for i in range(self.fst_vt['SubDyn']['NMOutputs']):
            ln = []
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['MemberID_out'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['SubDyn']['NOutCnt'][i]))
            ln.append(" ".join(['{:^11d}'.format(node) for node in self.fst_vt['SubDyn']['NodeCnt'][i]]))
            f.write(" ".join(ln) + '\n')
        f.write('------------------------- SDOutList: The next line(s) contains a list of output parameters that will be output in <rootname>.SD.out or <rootname>.out. ------\n')
        outlist = self.get_outlist(self.fst_vt['outlist'], ['SubDyn'])
        for channel_list in outlist:
            for i in range(len(channel_list)):
                f.write('"' + channel_list[i] + '"\n')
        f.write('END of output channels and end of file. (the word "END" must appear in the first 3 columns of this line)\n')
        f.flush()
        os.fsync(f)
        f.close()

    def write_ExtPtfm(self):
        # Generate ExtPtfm input file

        if self.fst_vt['ExtPtfm']['FileFormat'] == 0:
            None
            # self.write_Guyan() # TODO: need to impliment this. An example file not found to test
        elif self.fst_vt['ExtPtfm']['FileFormat'] == 1:
            self.write_Superelement()


        self.fst_vt['Fst']['SubFile'] = self.FAST_namingOut + '_ExtPtfm.dat'
        ep_file = os.path.join(self.FAST_runDirectory, self.fst_vt['Fst']['SubFile'])
        f = open(ep_file, 'w')

        f.write('---------------------- EXTPTFM INPUT FILE --------------------------------------\n')
        f.write('Comment describing the model\n')
        f.write('---------------------- SIMULATION CONTROL --------------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ExtPtfm']['Echo'], 'Echo', '- Echo input data to <RootName>.ech (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ExtPtfm']['DT'], 'DT', '- Communication interval for controllers (s) (or "default")\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['ExtPtfm']['IntMethod'], 'IntMethod', '- Integration Method {1:RK4; 2:AB4, 3:ABM4} (switch)\n'))
        f.write('---------------------- REDUCTION INPUTS ----------------------------------------\n')
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['ExtPtfm']['FileFormat'], 'FileFormat', '- File Format {0:Guyan; 1:FlexASCII} (switch)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ExtPtfm']['Red_FileName'], 'Red_FileName', '- Path of the file containing Guyan/Craig-Bampton inputs (-)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ExtPtfm']['RedCst_FileName'], 'RedCst_FileName', '- Path of the file containing Guyan/Craig-Bampton constant inputs (-) (currently unused)\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['ExtPtfm']['NActiveDOFList'], 'NActiveDOFList', '- Number of active CB mode listed in ActiveDOFList, use -1 for all modes (integer)\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join([f'{val}' for val in self.fst_vt['ExtPtfm']['ActiveDOFList']]), 'ActiveDOFList', '- List of CB modes index that are active, [unused if NActiveDOFList<=0]\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['ExtPtfm']['NInitPosList'], 'NInitPosList', '- Number of initial positions listed in InitPosList, using 0 implies all DOF initialized to 0  (integer)\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join([f'{val}' for val in self.fst_vt['ExtPtfm']['InitPosList']]), 'InitPosList', '- List of initial positions for the CB modes  [unused if NInitPosList<=0 or EquilStart=True]\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['ExtPtfm']['NInitVelList'], 'NInitVelList', '- Number of initial positions listed in InitVelList, using 0 implies all DOF initialized to 0  (integer)\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join([f'{val}' for val in self.fst_vt['ExtPtfm']['InitVelList']]), 'InitVelList', '- List of initial velocities for the CB modes  [unused if NInitVelPosList<=0 or EquilStart=True]\n'))

        f.write('---------------------- OUTPUT --------------------------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ExtPtfm']['SumPrint'], 'SumPrint', '- Print summary data to <RootName>.sum (flag)\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['ExtPtfm']['OutFile'], 'OutFile', '- Switch to determine where output will be placed: {1: in module output file only; 2: in glue code output file only; 3: both} (currently unused)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ExtPtfm']['TabDelim'], 'TabDelim', '- Use tab delimiters in text tabular output file? (flag) (currently unused)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['ExtPtfm']['OutFmt'], 'OutFmt', '- Format used for text tabular output (except time).  Resulting field should be 10 characters. (quoted string) (currently unused)\n'))
        f.write('{:<22f} {:<11} {:}'.format(self.fst_vt['ExtPtfm']['TStart'], 'TStart', '- Time to begin tabular output (s) (currently unused)\n'))
        f.write('                    OutList      - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)\n')
        outlist = self.get_outlist(self.fst_vt['outlist'], ['ExtPtfm'])

        for channel_list in outlist:
            for i in range(len(channel_list)):
                f.write('"' + channel_list[i] + '"\n')
        f.write('END of input file (the word "END" must appear in the first 3 columns of the last OutList line)\n')
        f.flush()
        os.fsync(f)
        f.close()



    def write_Superelement(self):

        def toString(SuperElement):
            # Function based on https://github.com/OpenFAST/openfast_toolbox/blob/353643ed917d113ec8dfd765813fef7d09752757/openfast_toolbox/io/fast_input_file.py#L2034
            # Developed by Emmanuel Branlard (https://github.com/ebranlard)
            s=''
            s+='!Comment\n'
            s+='!Comment Flex 5 Format\n'
            s+='!Dimension: {}\n'.format(SuperElement['nDOF'])
            s+='!Time increment in simulation: {}\n'.format(SuperElement['dt'])
            s+='!Total simulation time in file: {}\n'.format(SuperElement['T'])

            s+='\n!Mass Matrix\n'
            s+='!Dimension: {}\n'.format(SuperElement['nDOF'])
            s+='\n'.join(''.join('{:16.8e}'.format(x) for x in y) for y in SuperElement['MassMatrix'])

            s+='\n\n!Stiffness Matrix\n'
            s+='!Dimension: {}\n'.format(SuperElement['nDOF'])
            s+='\n'.join(''.join('{:16.8e}'.format(x) for x in y) for y in SuperElement['StiffnessMatrix'])

            s+='\n\n!Damping Matrix\n'
            s+='!Dimension: {}\n'.format(SuperElement['nDOF'])
            s+='\n'.join(''.join('{:16.8e}'.format(x) for x in y) for y in SuperElement['DampingMatrix'])

            s+='\n\n!Loading and Wave Elevation\n'
            s+='!Dimension: 1 time column -  {} force columns\n'.format(SuperElement['nDOF'])
            s+='\n'.join(''.join('{:16.8e}'.format(x) for x in y) for y in SuperElement['Loading'])
            return s

        # Generate Superelement input file
        self.fst_vt['ExtPtfm']['Red_FileName'] = self.FAST_namingOut + '_ExtPtfm_SE.dat'
        se_file = os.path.join(self.FAST_runDirectory, self.fst_vt['ExtPtfm']['Red_FileName'])
        f = open(se_file, 'w')

        f.write(toString(self.fst_vt['ExtPtfm']['FlexASCII']))

        f.flush()
        os.fsync(f)
        f.close()

    def write_MAP(self):

        # Generate MAP++ input file
        self.fst_vt['Fst']['MooringFile'] = self.FAST_namingOut + '_MAP.dat'
        map_file = os.path.join(self.FAST_runDirectory, self.fst_vt['Fst']['MooringFile'])
        f = open(map_file, 'w')
        
        f.write('---------------------- LINE DICTIONARY ---------------------------------------\n')
        f.write(" ".join(['{:<11s}'.format(i) for i in ['LineType', 'Diam', 'MassDenInAir', 'EA', 'CB', 'CIntDamp', 'Ca', 'Cdn', 'Cdt']])+'\n')
        f.write(" ".join(['{:<11s}'.format(i) for i in ['(-)', '(m)', '(kg/m)', '(N)', '(-)', '(Pa-s)', '(-)', '(-)', '(-)']])+'\n')
        ln =[]
        for i in range(1):
            ln = []
            ln.append('{:^11}'.format(self.fst_vt['MAP']['LineType'][i]))
            ln.append('{:^11}'.format(self.fst_vt['MAP']['Diam'][i]))
            ln.append('{:^11}'.format(self.fst_vt['MAP']['MassDenInAir'][i]))
            ln.append('{:^11}'.format(self.fst_vt['MAP']['EA'][i]))
            ln.append('{:<11}'.format(self.fst_vt['MAP']['CB'][i]))
            ln.append('{:<11}'.format(self.fst_vt['MAP']['CIntDamp'][i]))
            ln.append('{:<11}'.format(self.fst_vt['MAP']['Ca'][i]))
            ln.append('{:<11}'.format(self.fst_vt['MAP']['Cdn'][i]))
            ln.append('{:<11}'.format(self.fst_vt['MAP']['Cdt'][i]))
        f.write(" ".join(ln) + '\n')
        f.write('---------------------- NODE PROPERTIES ---------------------------------------\n')
        f.write(" ".join(['{:<11s}'.format(i) for i in ['Node', 'Type', 'X', 'Y', 'Z', 'M', 'B', 'FX', 'FY', 'FZ']])+'\n')
        f.write(" ".join(['{:<11s}'.format(i) for i in ['(-)', '(-)', '(m)', '(m)', '(m)', '(kg)', '(m^3)', '(N)', '(N)', '(N)']])+'\n')
        for i, type in enumerate(self.fst_vt['MAP']['Type']):
            ln = []
            ln.append('{:<11}'.format(self.fst_vt['MAP']['Node'][i]))
            ln.append('{:<11}'.format(self.fst_vt['MAP']['Type'][i]))
            ln.append('{:<11}'.format(self.fst_vt['MAP']['X'][i]))
            ln.append('{:<11}'.format(self.fst_vt['MAP']['Y'][i]))
            ln.append('{:<11}'.format(self.fst_vt['MAP']['Z'][i]))
            ln.append('{:<11}'.format(self.fst_vt['MAP']['M'][i]))
            ln.append('{:<11}'.format(self.fst_vt['MAP']['B'][i]))
            ln.append('{:<11}'.format(self.fst_vt['MAP']['FX'][i]))
            ln.append('{:<11}'.format(self.fst_vt['MAP']['FY'][i]))
            ln.append('{:<11}'.format(self.fst_vt['MAP']['FZ'][i]))
            f.write(" ".join(ln) + '\n')
        f.write('---------------------- LINE PROPERTIES ---------------------------------------\n')
        f.write(" ".join(['{:<11s}'.format(i) for i in ['Line', 'LineType', 'UnstrLen', 'NodeAnch', 'NodeFair', 'Flags']])+'\n')
        f.write(" ".join(['{:<11s}'.format(i) for i in ['(-)', '(-)', '(m)', '(-)', '(-)', '(-)']])+'\n')
        for i in range(len(self.fst_vt['MAP']['Line'])):
            ln = []
            ln.append('{:^11d}'.format(self.fst_vt['MAP']['Line'][i]))
            ln.append('{:^11}'.format(self.fst_vt['MAP']['LineType'][i]))
            ln.append('{:^11}'.format(self.fst_vt['MAP']['UnstrLen'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['MAP']['NodeAnch'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['MAP']['NodeFair'][i]))
            # ln.append('{:^11}'.format(self.fst_vt['MAP']['Outputs'][i]))
            ln.append('{:<11}'.format(" ".join(self.fst_vt['MAP']['Flags'][i])))
            f.write(" ".join(ln) + '\n')
        ln =[]
        f.write('---------------------- SOLVER OPTIONS-----------------------------------------\n')
        f.write('{:<11s}'.format('Option'+'\n'))
        f.write('{:<11s}'.format('(-)')+'\n')
        for i in range(len(self.fst_vt['MAP']['Option'])):
            ln = []
            ln.append('{:<11}'.format(" ".join(self.fst_vt['MAP']['Option'][i])))
            f.write("\n".join(ln) + '\n')
        ln = []
        f.write('\n') # adding a blank line after all solver options
        f.flush()
        os.fsync(f)
        f.close()

    def write_MoorDyn(self):

        self.fst_vt['Fst']['MooringFile'] = self.FAST_namingOut + '_MoorDyn.dat'
        moordyn_file = os.path.join(self.FAST_runDirectory, self.fst_vt['Fst']['MooringFile'])
        f = open(moordyn_file, 'w')

        f.write('--------------------- MoorDyn Input File ------------------------------------\n')
        f.write('Generated with OpenFAST_IO\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['MoorDyn']['Echo'], 'Echo', '- echo the input file data (flag)\n'))
        f.write('----------------------- LINE TYPES ------------------------------------------\n')
        f.write(" ".join(['{:^11s}'.format(i) for i in ['Name', 'Diam', 'MassDen', 'EA', 'BA/-zeta', 'EI', 'Cd', 'Ca', 'CdAx', 'CaAx']])+'\n')
        f.write(" ".join(['{:^11s}'.format(i) for i in ['(-)', '(m)', '(kg/m)', '(N)', '(N-s/-)', '(-)', '(-)', '(-)', '(-)', '(-)']])+'\n')
        for i in range(len(self.fst_vt['MoorDyn']['Name'])):
            ln = []
            ln.append('{:^11}'.format(self.fst_vt['MoorDyn']['Name'][i]))
            ln.append('{:^11.4f}'.format(self.fst_vt['MoorDyn']['Diam'][i]))
            ln.append('{:^11.4f}'.format(self.fst_vt['MoorDyn']['MassDen'][i]))
            ln.append('{:^11.4e}'.format(self.fst_vt['MoorDyn']['EA'][i]))
            ln.append('{:^11.4f}'.format(self.fst_vt['MoorDyn']['BA_zeta'][i]))
            ln.append('{:^11.4e}'.format(self.fst_vt['MoorDyn']['EI'][i]))
            ln.append('{:^11.4f}'.format(self.fst_vt['MoorDyn']['Cd'][i]))
            ln.append('{:^11.4f}'.format(self.fst_vt['MoorDyn']['Ca'][i]))
            ln.append('{:^11.4f}'.format(self.fst_vt['MoorDyn']['CdAx'][i]))
            ln.append('{:^11.4f}'.format(self.fst_vt['MoorDyn']['CaAx'][i]))
            f.write(" ".join(ln) + '\n')
        if 'Rod_Name' in self.fst_vt['MoorDyn'] and self.fst_vt['MoorDyn']['Rod_Name']:
            f.write('----------------------- ROD TYPES ------------------------------------------\n')
            f.write(" ".join(['{:^11s}'.format(i) for i in ['TypeName', 'Diam', 'Mass/m', 'Cd', 'Ca', 'CdEnd', 'CaEnd']])+'\n')
            f.write(" ".join(['{:^11s}'.format(i) for i in ['(name)', '(m)', '(kg/m)', '(-)', '(-)', '(-)', '(-)']])+'\n')
            for i in range(len(self.fst_vt['MoorDyn']['Rod_Name'])):
                ln = []
                ln.append('{:^11}'.format(self.fst_vt['MoorDyn']['Rod_Name'][i]))
                ln.append('{:^11.4f}'.format(self.fst_vt['MoorDyn']['Rod_Diam'][i]))
                ln.append('{:^11.4f}'.format(self.fst_vt['MoorDyn']['Rod_MassDen'][i]))
                ln.append('{:^11.4f}'.format(self.fst_vt['MoorDyn']['Rod_Cd'][i]))
                ln.append('{:^11.4f}'.format(self.fst_vt['MoorDyn']['Rod_Ca'][i]))
                ln.append('{:^11.4f}'.format(self.fst_vt['MoorDyn']['Rod_CdEnd'][i]))
                ln.append('{:^11.4f}'.format(self.fst_vt['MoorDyn']['Rod_CaEnd'][i]))
                f.write(" ".join(ln) + '\n')


        if 'Body_ID' in self.fst_vt['MoorDyn'] and self.fst_vt['MoorDyn']['Body_ID']:
            f.write('----------------------- BODIES ------------------------------------------\n')
            f.write(" ".join(['{:^11s}'.format(i) for i in ['ID', 'Attachement', 'X0', 'Y0', 'Z0', 'r0',  'p0','y0','Mass','CG*','I*','Volume','CdA*','Ca*']])+'\n')
            f.write(" ".join(['{:^11s}'.format(i) for i in ['(#)', '(word)', '(m)', '(m)', '(m)', '(deg)',  '(deg)','(deg)','(kg)','(m)','(kg-m^2)','(m^3)','m^2','(kg/m^3)']])+'\n')
            for i in range(len(self.fst_vt['MoorDyn']['Body_ID'])):
                ln = []
                ln.append('{:^11d}'.format(self.fst_vt['MoorDyn']['Body_ID'][i]))
                ln.append('{:^11}'.format(self.fst_vt['MoorDyn']['Body_Attachment'][i]))
                ln.append('{:^11.4f}'.format(self.fst_vt['MoorDyn']['X0'][i]))
                ln.append('{:^11.4f}'.format(self.fst_vt['MoorDyn']['Y0'][i]))
                ln.append('{:^11.4f}'.format(self.fst_vt['MoorDyn']['Z0'][i]))
                ln.append('{:^11.4e}'.format(self.fst_vt['MoorDyn']['r0'][i]))
                ln.append('{:^11.4e}'.format(self.fst_vt['MoorDyn']['p0'][i]))
                ln.append('{:^11.4e}'.format(self.fst_vt['MoorDyn']['y0'][i]))
                ln.append('{:^11.4e}'.format(self.fst_vt['MoorDyn']['Body_Mass'][i]))
                ln.append('|'.join(['{:^.4f}'.format(a) for a in self.fst_vt['MoorDyn']['Body_CG'][i]]))
                ln.append('|'.join(['{:^.4e}'.format(a) for a in self.fst_vt['MoorDyn']['Body_I'][i]]))
                ln.append('{:^11.4e}'.format(self.fst_vt['MoorDyn']['Body_Volume'][i]))
                ln.append('|'.join(['{:^.4f}'.format(a) for a in self.fst_vt['MoorDyn']['Body_CdA'][i]]))
                ln.append('|'.join(['{:^.4f}'.format(a) for a in self.fst_vt['MoorDyn']['Body_Ca'][i]]))
                f.write(" ".join(ln) + '\n')

        if 'Rod_ID' in self.fst_vt['MoorDyn'] and self.fst_vt['MoorDyn']['Rod_ID']:
            f.write('----------------------- RODS ------------------------------------------\n')
            f.write(" ".join(['{:^11s}'.format(i) for i in ['ID', 'RodType', 'Attachment', 'Xa', 'Ya', 'Za', 'Xb','Yb','Zb','NumSegs','RodOutputs']])+'\n')
            f.write(" ".join(['{:^11s}'.format(i) for i in ['(#)', '(name)', '(word/ID)', '(m)', '(m)', '(m)', '(m)','(m)','(m)','(-)','(-)']])+'\n')
            for i in range(len(self.fst_vt['MoorDyn']['Rod_ID'])):
                ln = []
                ln.append('{:^11}'.format(self.fst_vt['MoorDyn']['Rod_ID'][i]))
                ln.append('{:^11}'.format(self.fst_vt['MoorDyn']['Rod_Type'][i]))
                ln.append('{:^11}'.format(self.fst_vt['MoorDyn']['Rod_Attachment'][i]))
                ln.append('{:^11.4f}'.format(self.fst_vt['MoorDyn']['Xa'][i]))
                ln.append('{:^11.4f}'.format(self.fst_vt['MoorDyn']['Ya'][i]))
                ln.append('{:^11.4f}'.format(self.fst_vt['MoorDyn']['Za'][i]))
                ln.append('{:^11.4f}'.format(self.fst_vt['MoorDyn']['Xb'][i]))
                ln.append('{:^11.4f}'.format(self.fst_vt['MoorDyn']['Yb'][i]))
                ln.append('{:^11.4f}'.format(self.fst_vt['MoorDyn']['Zb'][i]))
                ln.append('{:^11d}'.format(self.fst_vt['MoorDyn']['Rod_NumSegs'][i]))
                ln.append('{:^11}'.format(self.fst_vt['MoorDyn']['RodOutputs'][i]))
                f.write(" ".join(ln) + '\n')
        

        f.write('---------------------- POINTS --------------------------------\n')
        f.write(" ".join(['{:^11s}'.format(i) for i in ['ID', 'Attachment', 'X', 'Y', 'Z', 'M', 'V', 'CdA', 'CA']])+'\n')
        f.write(" ".join(['{:^11s}'.format(i) for i in ['(-)', '(-)', '(m)', '(m)', '(m)', '(kg)', '(m^3)', '(m^2)', '(-)']])+'\n')
        for i in range(len(self.fst_vt['MoorDyn']['Point_ID'])):
            ln = []
            ln.append('{:^11d}'.format(self.fst_vt['MoorDyn']['Point_ID'][i]))
            ln.append('{:^11}'.format(self.fst_vt['MoorDyn']['Attachment'][i]))
            ln.append('{:^11.4f}'.format(self.fst_vt['MoorDyn']['X'][i]))
            ln.append('{:^11.4f}'.format(self.fst_vt['MoorDyn']['Y'][i]))
            ln.append('{:^11.4f}'.format(self.fst_vt['MoorDyn']['Z'][i]))
            ln.append('{:^11.4e}'.format(self.fst_vt['MoorDyn']['M'][i]))
            ln.append('{:^11.4e}'.format(self.fst_vt['MoorDyn']['V'][i]))
            ln.append('{:^11.4f}'.format(self.fst_vt['MoorDyn']['CdA'][i]))
            ln.append('{:^11.4f}'.format(self.fst_vt['MoorDyn']['CA'][i]))
            f.write(" ".join(ln) + '\n')
        f.write('---------------------- LINES --------------------------------------\n')
        f.write(" ".join(['{:^11s}'.format(i) for i in ['Line', 'LineType', 'AttachA', 'AttachB', 'UnstrLen', 'NumSegs',  'Outputs']])+'\n')
        f.write(" ".join(['{:^11s}'.format(i) for i in ['(-)', '(-)', '(-)', '(-)', '(m)', '(-)',  '(-)']])+'\n')
        for i in range(len(self.fst_vt['MoorDyn']['Line_ID'])):
            ln = []
            ln.append('{:^11d}'.format(self.fst_vt['MoorDyn']['Line_ID'][i]))
            ln.append('{:^11}'.format(self.fst_vt['MoorDyn']['LineType'][i]))
            ln.append('{:^11}'.format(self.fst_vt['MoorDyn']['AttachA'][i]))
            ln.append('{:^11}'.format(self.fst_vt['MoorDyn']['AttachB'][i]))
            ln.append('{:^11.4f}'.format(self.fst_vt['MoorDyn']['UnstrLen'][i]))
            ln.append('{:^11d}'.format(self.fst_vt['MoorDyn']['NumSegs'][i]))
            ln.append('{:^11}'.format(self.fst_vt['MoorDyn']['Outputs'][i]))
            f.write(" ".join(ln) + '\n')

        if 'ChannelID' in self.fst_vt['MoorDyn'] and self.fst_vt['MoorDyn']['ChannelID']: # There are control inputs
            f.write('---------------------- CONTROL ---------------------------------------\n')
            f.write(" ".join(['{:^11s}'.format(i) for i in ['ChannelID', 'Line(s)']])+'\n')
            f.write(" ".join(['{:^11s}'.format(i) for i in ['()', '(,)']])+'\n')
            for i_line in range(len(self.fst_vt['MoorDyn']['ChannelID'])):
                ln = []
                ln.append('{:^11d}'.format(self.fst_vt['MoorDyn']['ChannelID'][i_line]))
                ln.append(','.join(self.fst_vt['MoorDyn']['Lines_Control'][i_line]))
                f.write(" ".join(ln) + '\n')

        f.write('---------------------- SOLVER OPTIONS ---------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['MoorDyn']['dtM'], 'dtM', '- time step to use in mooring integration (s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['MoorDyn']['kbot'], 'kbot', '- bottom stiffness (Pa/m)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['MoorDyn']['cbot'], 'cbot', '- bottom damping (Pa-s/m)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['MoorDyn']['dtIC'], 'dtIC', '- time interval for analyzing convergence during IC gen (s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['MoorDyn']['TmaxIC'], 'TmaxIC', '- max time for ic gen (s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['MoorDyn']['CdScaleIC'], 'CdScaleIC', '- factor by which to scale drag coefficients during dynamic relaxation (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['MoorDyn']['threshIC'], 'threshIC', '- threshold for IC convergence (-)\n'))
        if 'options' in self.fst_vt['MoorDyn']:
            if 'inertialF' in self.fst_vt['MoorDyn']['options']:
                f.write('{:<22d} {:<11} {:}'.format(int(self.fst_vt['MoorDyn']['inertialF']), 'inertialF', '- Compute the inertial forces (0: no, 1: yes). Switch to 0 if you get: Warning: extreme pitch moment from body-attached Rod.\n'))

            if 'WaterKin' in self.fst_vt['MoorDyn']['options']:
                self.fst_vt['MoorDyn']['WaterKin_file'] = self.FAST_namingOut + '_WaterKin.dat'
                f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['MoorDyn']['WaterKin_file']+'"', 'WaterKin', '- WaterKin input file\n'))

        
        # f.write('{:^11s} {:<11} {:}'.format(self.fst_vt['MoorDyn']['WaterKin'], 'WaterKin', 'Handling of water motion (0=off, 1=on)\n'))
        f.write('------------------------ OUTPUTS --------------------------------------------\n')
        outlist = self.get_outlist(self.fst_vt['outlist'], ['MoorDyn'])
        for channel_list in outlist:
            for i in range(len(channel_list)):
                f.write('"' + channel_list[i] + '"\n')
        f.write('END\n')
        f.write('------------------------- need this line --------------------------------------\n')

        f.flush()
        os.fsync(f)
        f.close()

    def write_WaterKin(self,WaterKin_file):
        f = open(WaterKin_file, 'w')

        f.write('MoorDyn v2 (Feb 2022) Waves and Currents input file set up for USFLOWT\n')
        f.write('Wave kinematics that will have an impact over the cans and the mooring lines.\n')
        f.write('--------------------------- WAVES -------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['WaterKin']['WaveKinMod'], 'WaveKinMod', '- type of wave input {0 no waves; 3 set up grid of wave data based on time series}\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['WaterKin']['WaveKinFile'], 'WaveKinFile', '- file containing wave elevation time series at 0,0,0\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['WaterKin']['dtWave'], 'dtWave', '- time step to use in setting up wave kinematics grid (s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['WaterKin']['WaveDir'], 'WaveDir', '- wave heading (deg)\n'))
        f.write('{:<22} {:}'.format(self.fst_vt['WaterKin']['X_Type'], '- X wave input type (0: not used; 1: list values in ascending order; 2: uniform specified by -xlim, xlim, num)\n'))
        f.write('{:<22} {:}'.format(', '.join(['{:.3f}'.format(i) for i in self.fst_vt['WaterKin']['X_Grid']]), '- X wave grid point data separated by commas\n'))
        f.write('{:<22} {:}'.format(self.fst_vt['WaterKin']['Y_Type'], '- Y wave input type (0: not used; 1: list values in ascending order; 2: uniform specified by -Ylim, Ylim, num)\n'))
        f.write('{:<22} {:}'.format(', '.join(['{:.3f}'.format(i) for i in self.fst_vt['WaterKin']['Y_Grid']]), '- Y wave grid point data separated by commas\n'))
        f.write('{:<22} {:}'.format(self.fst_vt['WaterKin']['Z_Type'], '- Z wave input type (0: not used; 1: list values in ascending order; 2: uniform specified by -Zlim, Zlim, num)\n'))
        f.write('{:<22} {:}'.format(', '.join(['{:.3f}'.format(i) for i in self.fst_vt['WaterKin']['Z_Grid']]), '- Z wave grid point data separated by commas\n'))
        f.write('--------------------------- CURRENT -------------------------------------\n')
        f.write('0                    CurrentMod  - type of current input {0 no current; 1 steady current profile described below} \n')
        f.write('z-depth     x-current      y-current\n')
        f.write('(m)           (m/s)         (m/s)\n')
        f.write('--------------------- need this line ------------------\n')

        f.flush()
        os.fsync(f)
        f.close()


        
        
    def write_StC(self,StC_vt,StC_filename):

        stc_file = os.path.join(self.FAST_runDirectory, StC_filename)
        f = open(stc_file, 'w')
        
        f.write('------- STRUCTURAL CONTROL (StC) INPUT FILE ----------------------------\n')
        f.write('Generated with OpenFAST_IO\n')
        
        f.write('---------------------- SIMULATION CONTROL --------------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(StC_vt['Echo'], 'Echo', '- Echo input data to "<rootname>.SD.ech" (flag)\n'))
        
        f.write('---------------------- StC DEGREES OF FREEDOM ----------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_DOF_MODE'], 'StC_DOF_MODE', '- DOF mode (switch) {0: No StC or TLCD DOF; 1: StC_X_DOF, StC_Y_DOF, and/or StC_Z_DOF (three independent StC DOFs); 2: StC_XY_DOF (Omni-Directional StC); 3: TLCD; 4: Prescribed force/moment time series}\n'))
        f.write('{!s:<22} {:<11} {:}'.format(StC_vt['StC_X_DOF'], 'StC_X_DOF', '- DOF on or off for StC X (flag) [Used only when StC_DOF_MODE=1]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(StC_vt['StC_Y_DOF'], 'StC_Y_DOF', '- DOF on or off for StC Y (flag) [Used only when StC_DOF_MODE=1]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(StC_vt['StC_Z_DOF'], 'StC_Z_DOF', '- DOF on or off for StC Z (flag) [Used only when StC_DOF_MODE=1]\n'))
        
        f.write('---------------------- StC LOCATION ------------------------------------------- [relative to the reference origin of component attached to]\n')
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_P_X'], 'StC_P_X', '- At rest X position of StC (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_P_Y'], 'StC_P_Y', '- At rest Y position of StC (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_P_Z'], 'StC_P_Z', '- At rest Z position of StC (m)\n'))
        
        f.write('---------------------- StC INITIAL CONDITIONS --------------------------------- [used only when StC_DOF_MODE=1 or 2]\n')
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_X_DSP'], 'StC_X_DSP', '- StC X initial displacement (m) [relative to at rest position]\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_Y_DSP'], 'StC_Y_DSP', '- StC Y initial displacement (m) [relative to at rest position]\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_Z_DSP'], 'StC_Z_DSP', '- StC Z initial displacement (m) [relative to at rest position; used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(StC_vt['StC_Z_PreLd'], 'StC_Z_PreLd', '- StC Z pre-load (N) {"gravity" to offset for gravity load; "none" or 0 to turn off} [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]\n'))
        
        f.write('---------------------- StC CONFIGURATION -------------------------------------- [used only when StC_DOF_MODE=1 or 2]\n')
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_X_PSP'], 'StC_X_PSP', '- Positive stop position (maximum X mass displacement) (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_X_NSP'], 'StC_X_NSP', '- Negative stop position (minimum X mass displacement) (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_Y_PSP'], 'StC_Y_PSP', '- Positive stop position (maximum Y mass displacement) (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_Y_NSP'], 'StC_Y_NSP', '- Negative stop position (minimum Y mass displacement) (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_Z_PSP'], 'StC_Z_PSP', '- Positive stop position (maximum Z mass displacement) (m) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_Z_NSP'], 'StC_Z_NSP', '- Negative stop position (minimum Z mass displacement) (m) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]\n'))
        
        f.write('---------------------- StC MASS, STIFFNESS, & DAMPING ------------------------- [used only when StC_DOF_MODE=1 or 2]\n')
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_X_M'], 'StC_X_M', '- StC X mass (kg) [must equal StC_Y_M for StC_DOF_MODE = 2]\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_Y_M'], 'StC_Y_M', '- StC Y mass (kg) [must equal StC_X_M for StC_DOF_MODE = 2]\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_Z_M'], 'StC_Z_M', '- StC Z mass (kg) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_XY_M'], 'StC_XY_M', '- StC Z mass (kg) [used only when StC_DOF_MODE=2]\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_X_K'], 'StC_X_K', '- StC X stiffness (N/m)\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_Y_K'], 'StC_Y_K', '- StC Y stiffness (N/m)\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_Z_K'], 'StC_Z_K', '- StC Z stiffness (N/m) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_X_C'], 'StC_X_C', '- StC X damping (N/(m/s))\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_Y_C'], 'StC_Y_C', '- StC Y damping (N/(m/s))\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_Z_C'], 'StC_Z_C', '- StC Z damping (N/(m/s)) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_X_KS'], 'StC_X_KS', '- Stop spring X stiffness (N/m)\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_Y_KS'], 'StC_Y_KS', '- Stop spring Y stiffness (N/m)\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_Z_KS'], 'StC_Z_KS', '- Stop spring Z stiffness (N/m) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_X_CS'], 'StC_X_CS', '- Stop spring X damping (N/(m/s))\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_Y_CS'], 'StC_Y_CS', '- Stop spring Y damping (N/(m/s))\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_Z_CS'], 'StC_Z_CS', '- Stop spring Z damping (N/(m/s)) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]\n'))
        
        f.write('---------------------- StC USER-DEFINED SPRING FORCES ------------------------- [used only when StC_DOF_MODE=1 or 2]\n')
        f.write('{!s:<22} {:<11} {:}'.format(StC_vt['Use_F_TBL'], 'Use_F_TBL', '- Use spring force from user-defined table (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['NKInpSt'], 'NKInpSt', '- Number of spring force input stations\n'))
        
        f.write('---------------------- StC SPRING FORCES TABLE -------------------------------- [used only when StC_DOF_MODE=1 or 2]\n')
        f.write('X                F_X               Y              F_Y              Z              F_Z\n')
        f.write('(m)               (N)              (m)             (N)             (m)             (N)\n')
        table = StC_vt['SpringForceTable']
        for x, f_x, y, f_y, z, f_z in zip(table['X'],table['F_X'],table['Y'],table['F_Y'],table['Z'],table['F_Z']):
            row = [x, f_x, y, f_y, z, f_z]
            f.write(' '.join(['{: 2.8e}'.format(val) for val in row])+'\n')
        
        f.write('---------------------- StructUserProp CONTROL -------------------------------------------- [used only when StC_DOF_MODE=1 or 2]\n')
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_CMODE'],     'StC_CMODE',        '- Control mode (switch) {0:none; 1: Semi-Active Control Mode; 2: Active Control Mode}\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_CChan'],     'StC_CChan',        '- Control channel group (1:10) for stiffness and damping (StC_[XYZ]_K, StC_[XYZ]_C, and StC_[XYZ]_Brake) (specify additional channels for blade instances of StC active control -- one channel per blade) [used only when StC_DOF_MODE=1 or 2, and StC_CMODE=4 or 5]\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_SA_MODE'],   'StC_SA_MODE',      '- Semi-Active control mode {1: velocity-based ground hook control; 2: Inverse velocity-based ground hook control; 3: displacement-based ground hook control 4: Phase difference Algorithm with Friction Force 5: Phase difference Algorithm with Damping Force} (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_X_C_HIGH'],  'StC_X_C_HIGH',     '- StC X high damping for ground hook control\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_X_C_LOW'],   'StC_X_C_LOW',      '- StC X low damping for ground hook control\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_Y_C_HIGH'],  'StC_Y_C_HIGH',     '- StC Y high damping for ground hook control\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_Y_C_LOW'],   'StC_Y_C_LOW',      '- StC Y low damping for ground hook control\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_Z_C_HIGH'],  'StC_Z_C_HIGH',     '- StC Z high damping for ground hook control [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_Z_C_LOW'],   'StC_Z_C_LOW',      '- StC Z low damping for ground hook control  [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_X_C_BRAKE'], 'StC_X_C_BRAKE',    '- StC X high damping for braking the StC (Don''t use it now. should be zero)\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_Y_C_BRAKE'], 'StC_Y_C_BRAKE',    '- StC Y high damping for braking the StC (Don''t use it now. should be zero)\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['StC_Z_C_BRAKE'], 'StC_Z_C_BRAKE',    '- StC Z high damping for braking the StC (Don''t use it now. should be zero) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]\n'))

        f.write('---------------------- TLCD --------------------------------------------------- [used only when StC_DOF_MODE=3]\n')
        f.write('{:<22} {:<11} {:}'.format(StC_vt['L_X'], 'L_X', '- X TLCD total length (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['B_X'], 'B_X', '- X TLCD horizontal length (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['area_X'], 'area_X', '- X TLCD cross-sectional area of vertical column (m^2)\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['area_ratio_X'], 'area_ratio_X', '- X TLCD cross-sectional area ratio (vertical column area divided by horizontal column area) (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['headLossCoeff_X'], 'headLossCoeff_X', '- X TLCD head loss coeff (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['rho_X'], 'rho_X', '- X TLCD liquid density (kg/m^3)\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['L_Y'], 'L_Y', '- Y TLCD total length (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['B_Y'], 'B_Y', '- Y TLCD horizontal length (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['area_Y'], 'area_Y', '- Y TLCD cross-sectional area of vertical column (m^2)\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['area_ratio_Y'], 'area_ratio_Y', '- Y TLCD cross-sectional area ratio (vertical column area divided by horizontal column area) (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['headLossCoeff_Y'], 'headLossCoeff_Y', '- Y TLCD head loss coeff (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(StC_vt['rho_Y'], 'rho_Y', '- Y TLCD liquid density (kg/m^3)\n'))
        
        f.write('---------------------- PRESCRIBED TIME SERIES --------------------------------- [used only when StC_DOF_MODE=4]\n')
        f.write('{:<22} {:<11} {:}'.format(StC_vt['PrescribedForcesCoord'], 'PrescribedForcesCoord',    '- Prescribed forces are in global or local coordinates (switch) {1: global; 2: local}\n'))
        f.write('{!s:<22} {:<11} {:}'.format(StC_vt['PrescribedForcesFile'], 'PrescribedForcesFile', '- Time series force and moment (7 columns of time, FX, FY, FZ, MX, MY, MZ)\n'))
        f.write('-------------------------------------------------------------------------------\n')

        f.flush()
        os.fsync(f)
        f.close()

if __name__=="__main__":

    from openfast_io.FAST_reader import InputReader_OpenFAST
    from openfast_io.FileTools import check_rtest_cloned
    from pathlib import Path

    fst_update = {}
    fst_update['Fst', 'TMax'] = 20.
    fst_update['AeroDyn', 'TwrAero'] = False

    parent_dir = os.path.dirname( os.path.dirname( os.path.dirname( os.path.realpath(__file__) ) ) ) + os.sep
    build_of_io_dir = os.path.join(parent_dir, 'build_ofio')
    Path(build_of_io_dir).mkdir(parents=True, exist_ok=True)

    # Read the model
    fast = InputReader_OpenFAST()
    fast.FAST_InputFile = '5MW_Land_BD_DLL_WTurb.fst'   # FAST input file (ext=.fst)
    fast.FAST_directory = os.path.join(parent_dir, 'reg_tests', 'r-test', 
                                       'glue-codes', 'openfast', 
                                       '5MW_Land_BD_DLL_WTurb')   # Path to fst directory files
    
    check_rtest_cloned(os.path.join(fast.FAST_directory, fast.FAST_InputFile))
    
    fast.execute()
    
    # Write out the model
    fastout = InputWriter_OpenFAST()
    fastout.fst_vt = fast.fst_vt
    fastout.FAST_runDirectory = os.path.join(build_of_io_dir,'fast_write_main_test')
    fastout.FAST_namingOut = '5MW_Land_BD_DLL_WTurb_write'
    fastout.update(fst_update=fst_update)
    fastout.execute()
