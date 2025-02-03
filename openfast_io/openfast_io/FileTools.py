import os
import copy
import operator
import numpy as np
import yaml
import sys
from functools import reduce
from deepdiff import DeepDiff
try:
    import ruamel_yaml as ry
except Exception:
    try:
        import ruamel.yaml as ry
    except Exception:
        raise ImportError('No module named ruamel.yaml or ruamel_yaml')

"""
Common utilites for handling the text I/O for using AeroelasticSE
"""

def remove_numpy(fst_vt):
    # recursively move through nested dictionary, remove numpy data types
    # for formatting dictionaries before writing to yaml files

    def get_dict(vartree, branch):
        return reduce(operator.getitem, branch, vartree)

    def loop_dict(vartree, branch):
        if type(vartree) is not dict:
            return fst_vt
        for var in vartree.keys():
            branch_i = copy.copy(branch)
            branch_i.append(var)
            if type(vartree[var]) is dict:
                loop_dict(vartree[var], branch_i)
            else:
                data_type = type(get_dict(fst_vt, branch_i[:-1])[branch_i[-1]])

                if data_type in [np.int_, np.intc, np.intp, np.int8, np.int16, np.int32, np.int64, np.uint8, np.uint16, np.uint32, np.uint64]:
                    get_dict(fst_vt, branch_i[:-1])[branch_i[-1]] = int(get_dict(fst_vt, branch_i[:-1])[branch_i[-1]])
                elif data_type in [np.single, np.double, np.longdouble, np.csingle, np.cdouble, np.float16, np.float32, np.float64, np.complex64, np.complex128]:
                    get_dict(fst_vt, branch_i[:-1])[branch_i[-1]] = float(get_dict(fst_vt, branch_i[:-1])[branch_i[-1]])
                elif data_type in [np.bool_]:
                    get_dict(fst_vt, branch_i[:-1])[branch_i[-1]] = bool(get_dict(fst_vt, branch_i[:-1])[branch_i[-1]])
                elif data_type in [np.ndarray]:
                    get_dict(fst_vt, branch_i[:-1])[branch_i[-1]] = get_dict(fst_vt, branch_i[:-1])[branch_i[-1]].tolist()
                elif data_type in [list,tuple]:
                    for item in get_dict(fst_vt, branch_i[:-1])[branch_i[-1]]:
                        remove_numpy(item)

    # set fast variables to update values
    loop_dict(fst_vt, [])

    return fst_vt

def convert_str(val):
    # string parsing tool

    def try_type(val, data_type):
        try:
            data_type(val)
            return True
        except ValueError:
            return False

    if try_type(val, int):
        return int(val)
    elif try_type(val, float):
        return float(val)
    elif val=='True':
        return True
    elif val=='False':
        return False
    else:
        return val


def str_repeats(string):
    # find repeated sequences in a string
    for x in range(1, len(string)):
        substring = string[:x]
        if substring * (len(string)//len(substring))+(substring[:len(string)%len(substring)]) == string:
            return substring


def load_case_matrix(fname):
    # load the AeroelasticSE case_matrix

    fid = open(fname)
    text = fid.readlines()
    fid.close()

    module    = text[0].strip().split()
    var_names = text[1].strip().split()
    if len(module) == len(var_names)+1:
        var_names = ["ID"] + var_names

    out = {}
    for var in var_names:
        out[var] = []

    for ln in text[2:]:
        val_list = ln.strip().split()
        for var, val in zip(var_names, val_list):
            out[var].append(convert_str(val))

    return out

def load_yaml(fname_input, package=0):
    # Import a .yaml file

    if package == 0:
        with open(fname_input) as f:
            data = yaml.safe_load(f)
        return data

    elif package == 1:
        with open(fname_input, 'r') as myfile:
            text_input = myfile.read()
        myfile.close()
        ryaml = ry.YAML()
        return dict(ryaml.load(text_input))


def save_yaml(outdir, fname, data_out):

    if not os.path.isdir(outdir) and outdir!='':
        os.makedirs(outdir)
    fname = os.path.join(outdir, fname)

    data_out = remove_numpy(data_out)

    f = open(fname, "w")
    yaml=ry.YAML()
    yaml.default_flow_style = None
    yaml.width = float("inf")
    yaml.indent(mapping=4, sequence=6, offset=3)
    yaml.dump(data_out, f)
    f.close()

def print_yaml(data_struct):
    data_struct = remove_numpy(data_struct)
    yaml=ry.YAML()
    yaml.indent(mapping=4, sequence=6, offset=3)
    yaml.dump(data_struct,sys.stdout)

def select_cases(cases, var_sel, val_sel):
    # Find a variable value from the AeroelasticSE case_matrix

    n_var = len(var_sel)
    n_cases = len(cases[var_sel[0]])

    truth = [True]*n_cases
    for vari, vali in zip(var_sel, val_sel):
        test = [valj == vali for valj in cases[vari]]
        truth = [truthi and testi for truthi, testi in zip(truth, test)]

    case_idx = [i for i, x in enumerate(truth) if x]
    return case_idx

def get_dlc_label(cases, include_seed=True):
    # Get descriptive string describing IEC DLC cases from the case_matrix

    labels = []
    # from txt
    try:
        for idx in range(len(cases['DLC'])):

            DLC        = cases['DLC'][idx]
            wind_fname = cases['Filename'][idx]

            if DLC == 1.1:
                ntm      = wind_fname.split('NTM')[-1].split('_')
                ntm_U    = float(".".join(ntm[1].strip("U").split('.')[:-1]))
                ntm_Seed = float(".".join(ntm[2].strip("Seed").split('.')[:-1]))
                if include_seed == True:
                    label_i = "DLC 1.1, Normal Turbulence, U=%0.1f m/s, Seed=%d"%(ntm_U, ntm_Seed)
                else:
                    label_i = "DLC 1.1, Normal Turbulence, U=%0.1f m/s"%(ntm_U)

            if DLC == 1.3:
                etm   = wind_fname.split('ETM')[-1].split('_')
                etm_U = float(".".join(etm[1].strip("U").split('.')[:-1]))
                etm_Seed = float(".".join(etm[2].strip("Seed").split('.')[:-1]))
                if include_seed == True:
                    label_i = "DLC 1.3, Extreme Turbulence, U=%0.1f m/s, Seed=%d"%(etm_U, etm_Seed)
                else:
                    label_i = "DLC 1.3, Extreme Turbulence, U=%0.1f m/s"%(etm_U)

            if DLC == 1.4:
                ecd      = wind_fname.split('ECD')[-1].split('_')
                ecd_dir  = ecd[1]
                ecd_U    = float(".".join(ecd[2].strip("U").split('.')[:-1]))
                label_i  = "DLC 1.4, %s ECD, U=%0.1f m/s"%(ecd_dir, ecd_U)

            if DLC == 1.5:
                ews      = wind_fname.split('EWS')[-1].split('_')
                ews_type = ews[1]
                ews_dir  = ews[2]
                ews_U    = float(".".join(ews[3].strip("U").split('.')[:-1]))
                if ews_type == "H":
                    ews_type = "Hor."
                elif ews_type == "V":
                    ews_type = "Vert."
                if ews_dir == "P":
                    ews_dir = "Pos."
                elif ews_dir == "N":
                    ews_dir = "Neg."
                label_i = "DLC 1.5, Extreme %s Shear, %s, U=%0.1f m/s"%(ews_type, ews_dir, ews_U)

            if DLC == 6.1:
                label_i = "DLC 6.1"

            if DLC == 6.3:
                label_i = "DLC 6.3"

            labels.append(label_i)

    # From yaml
    except KeyError:
        for idx in range(len(cases[('IEC','DLC')])):

            DLC        = cases[('IEC','DLC')][idx]
            wind_fname = cases[('InflowWind', 'Filename')][idx]

            if DLC == 1.1:
                ntm      = wind_fname.split('NTM')[-1].split('_')
                ntm_U    = float(".".join(ntm[1].strip("U").split('.')[:-1]))
                ntm_Seed = float(".".join(ntm[2].strip("Seed").split('.')[:-1]))
                if include_seed == True:
                    label_i = "DLC 1.1, Normal Turbulence, U=%0.1f m/s, Seed=%d"%(ntm_U, ntm_Seed)
                else:
                    label_i = "DLC 1.1, Normal Turbulence, U=%0.1f m/s"%(ntm_U)

            if DLC == 1.3:
                etm   = wind_fname.split('ETM')[-1].split('_')
                etm_U = float(".".join(etm[1].strip("U").split('.')[:-1]))
                etm_Seed = float(".".join(etm[2].strip("Seed").split('.')[:-1]))
                if include_seed == True:
                    label_i = "DLC 1.3, Extreme Turbulence, U=%0.1f m/s, Seed=%d"%(etm_U, etm_Seed)
                else:
                    label_i = "DLC 1.3, Extreme Turbulence, U=%0.1f m/s"%(etm_U)

            if DLC == 1.4:
                ecd      = wind_fname.split('ECD')[-1].split('_')
                ecd_dir  = ecd[1]
                ecd_U    = float(".".join(ecd[2].strip("U").split('.')[:-1]))
                label_i  = "DLC 1.4, %s ECD, U=%0.1f m/s"%(ecd_dir, ecd_U)

            if DLC == 1.5:
                ews      = wind_fname.split('EWS')[-1].split('_')
                ews_type = ews[1]
                ews_dir  = ews[2]
                ews_U    = float(".".join(ews[3].strip("U").split('.')[:-1]))
                if ews_type == "H":
                    ews_type = "Hor."
                elif ews_type == "V":
                    ews_type = "Vert."
                if ews_dir == "P":
                    ews_dir = "Pos."
                elif ews_dir == "N":
                    ews_dir = "Neg."
                label_i = "DLC 1.5, Extreme %s Shear, %s, U=%0.1f m/s"%(ews_type, ews_dir, ews_U)

            if DLC == 6.1:
                label_i = "DLC 6.1"

            if DLC == 6.3:
                label_i = "DLC 6.3"

            labels.append(label_i)


    return labels

def load_file_list(fname_flist):
    # load list of filenames from file
    return np.genfromtxt(fname_flist, dtype='str')

def check_rtest_cloned(rtest_dir):
    # check if the rtest directory is cloned
    if not os.path.isdir(rtest_dir):
        raise FileNotFoundError(f"The directory {rtest_dir} does not exist. Please clone the r-test submodule. Try running `git submodule update --init --recursive`")

    return True


def remove_nested_keys(dictionary, keys_to_remove):
    for key in keys_to_remove:
        if key in dictionary:
            del dictionary[key]

    for value in dictionary.values():
        if isinstance(value, dict):
            remove_nested_keys(value, keys_to_remove)

    return dictionary

def removeDeactivatedModules(fst_vt):
    # Mapping of deactivated modules to their corresponding module names
    OFmodules = {
        'CompElast': {
            1: ['ElastoDyn', 'ElastoDynBlade', 'ElastoDynTower'],
            2: ['ElastoDyn', 'ElastoDynTower', 'BeamDyn', 'BeamDynBlade'],
            3: ['SimpleElastoDyn']
        },
        'CompInflow': {
            0: [],
            1: ['InflowWind']
        },
        'CompAero': {
            0: [],
            1: ['AeroDisk'], 
            2: ['AeroDyn', 'AeroDynBlade', 'AeroDynPolar']
        },
        'CompServo': {
            0: [],
            1: ['ServoDyn', 'DISCON_in']
        }, 
        'CompSeaSt': {
            0: [],
            1: ['SeaState']
        },
        'CompHydro': {
            0: [],
            1: ['HydroDyn']
        },
        'CompSub': {
            0: [],
            1: ['SubDyn'],
            2: ['read_ExtPtfm']
        },
        'CompMooring': {
            0: [],
            1: ['MAP'],
            2: [],
            3: ['MoorDyn', 'WaterKin'],
            4: []
        },
        'CompIce': {
            0: [],
            1: [],
            2: []
        }
        # 'MHK': {0:[]}, # no special handling for MHK
    }

    keys2keep = []
    keys2remove = []
    # loop throught the keys of OFmodules, and make two lists, one of the needed ones,
    # and one of the ones to remove, then remove the ones to remove
    for module, active in fst_vt['Fst'].items():
        if module in OFmodules:
            if active in OFmodules[module]:
                # get the list of modules to keep
                keys2keep.extend(OFmodules[module][active])

                # get the list of modules to remove
                for key, value in OFmodules[module].items():
                    if key != active:
                        keys2remove.extend(value)
    
    # remove the keys in keys2remove and NOT in keys2keep
    fst_vt = remove_nested_keys(fst_vt, [key for key in keys2remove if key not in keys2keep])

    return fst_vt

def cleanup_fstvt(fst_vt, ignoreVars=None, removeFileRef=False, removeArrayProps=False,
                    removeDeactivatedModules=False):
    # sanitize the dictionaries from numpy data types
    fst_vt = remove_numpy(fst_vt)

    if ignoreVars is not None:
        fst_vt = remove_nested_keys(fst_vt, ignoreVars)

    if removeFileRef: # not fair to compare file paths
        fileVars = ['af_coord', 'Filename_Uni', 'FileName_BTS', 'FileName_u', 'FileName_v', 'FileName_w', # TODO: orgainze these logically
                    'AFNames', 'ADBlFile1', 'ADBlFile2', 'ADBlFile3', 'NumCoords',
                    'DLL_FileName','DLL_InFile','af_data',
                    'PerfFileName',
                    'InflowFile',
                    'AeroFile',
                    'BldFile1', 'BldFile2', 'BldFile3', 
                    'TwrFile',
                    'EDFile',
                    'ServoFile',
                    'BldFile',
                    'SeaStFile',
                    'HydroFile',
                    'SubFile',
                    'MooringFile',
                    'Red_FileName',
                    'PrescribedForcesFile',
                    'actuatorDiskFile',
                    'BDBldFile(1)', 'BDBldFile(2)', 'BDBldFile(3)',
                    'OLAFInputFileName',
                    'EDFile_path',      
                    'BDBldFile(1_path)',
                    'BDBldFile(2_path)',
                    'BDBldFile(3_path)',
                    'InflowFile_path',  
                    'AeroFile_path',    
                    'ServoFile_path',   
                    'HydroFile_path',   
                    'SubFile_path',     
                    'MooringFile_path', 
                    'IceFile_path',     
                    'description',
                    ]
        fst_vt = remove_nested_keys(fst_vt, fileVars)

    if removeArrayProps: # we can have different array properties, if run through different tools
        arrayVars = ['BlSpn', 'BlCrvAC','BlSwpAC','BlCrvAng','BlTwist','BlChord','BlAFID',
                    'ac','PC_GS_KP','PC_GS_KI','WE_FOPoles','beam_stiff','attr','units']
        fst_vt = remove_nested_keys(fst_vt, arrayVars)


    if removeDeactivatedModules:
        fst_vt = removeDeactivatedModules(fst_vt)


    return fst_vt

def compare_fst_vt(fst_vt1, fst_vt2, ignoreVars = None, removeFileRef=False, removeArrayProps=False, print_diff=False):
    # Compare two FAST variable trees

    # sanitize the dictionaries from numpy data types
    fst_vt1 = cleanup_fstvt(fst_vt1, ignoreVars, removeFileRef, removeArrayProps)
    fst_vt2 = cleanup_fstvt(fst_vt2, ignoreVars, removeFileRef, removeArrayProps)

    diff = DeepDiff(fst_vt1, fst_vt2, ignore_numeric_type_changes=True,
                    # ignore_string_case=True,
                    # verbose_level = 2
                    )

    if diff == {}:
        print('No differences found between the two fst_vt.')

    if print_diff:
        print(diff.pretty())

    return diff
