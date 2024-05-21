import os
import copy
import operator
import numpy as np
import yaml
from functools import reduce
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
                elif data_type in [np.single, np.double, np.longdouble, np.csingle, np.cdouble, np.float_, np.float16, np.float32, np.float64, np.complex64, np.complex128]:
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
