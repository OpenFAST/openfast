import os, itertools
import numpy as np
from weis.aeroelasticse.FileTools import save_yaml

def save_case_matrix(matrix_out, change_vars, dir_matrix):
    # save matrix file
    if type(change_vars[0]) is tuple:
        n_header_lines = len(change_vars[0])
    else:
        change_vars = [(var,) for var in change_vars]
        n_header_lines = 1

    n_cases = np.shape(matrix_out)[0]
    matrix_out = np.hstack((np.asarray([[i] for i in range(n_cases)]), matrix_out))

    change_vars = [('Case_ID',)+('',)*(n_header_lines-1)] + change_vars
    # col_len = [max([len(val) for val in matrix_out[:,j]] + [len(change_vars[j][0]), len(change_vars[j][1])]) for j in range(len(change_vars))]
    col_len = [max([len(str(val)) for val in matrix_out[:,j]] + [len(change_vars[j][header_i]) for header_i in range(n_header_lines)]) for j in range(len(change_vars))]

    text_out = []
    for header_i in range(n_header_lines):
        text_out.append(''.join([val.center(col+2) for val, col in zip([var[header_i] for var in change_vars], col_len)])+'\n')

    for row in matrix_out:
        row_str = ''
        for val, col in zip(row, col_len):
            if val is not str:
                val = str(val)
            row_str += val.center(col+2)
        row_str += '\n'
        text_out.append(row_str)

    if not os.path.exists(dir_matrix):
        os.makedirs(dir_matrix)
    ofh = open(os.path.join(dir_matrix,'case_matrix.txt'),'w')
    for row in text_out:
        ofh.write(row)
    ofh.close()

def save_case_matrix_yaml(matrix_out, change_vars, dir_matrix, case_names):

    matrix_out_yaml = {}
    for var in change_vars:
        matrix_out_yaml[var] = []
    matrix_out_yaml['Case_ID'] = []
    matrix_out_yaml['Case_Name'] = []

    for i, row in enumerate(matrix_out):
        matrix_out_yaml['Case_ID'].append(i)
        matrix_out_yaml['Case_Name'].append(case_names[i])
        for val, var in zip(row, change_vars):
            if type(val) is list:
                if len(val) == 1:
                    val = val[0]
            if type(val) in [np.float32, np.float64, np.single, np.double, np.longdouble]:
                val = float(val)
            elif type(val) in [np.int8, np.int16, np.int32, np.int64, np.uint8, np.uint16, np.uint32, np.uint64, np.intc, np.uintc, np.uint]:
                val = int(val)
            elif type(val) in [np.array, np.ndarray]:
                val = val.tolist()
            elif type(val) in [np.str_]:
                val = str(val)
            # elif len(val) > 0:
            #     val = val.tolist()
            matrix_out_yaml[var].append(val)

    if not os.path.exists(dir_matrix):
        os.makedirs(dir_matrix)

    save_yaml(dir_matrix, 'case_matrix.yaml', matrix_out_yaml)

def case_naming(n_cases, namebase=None):
    # case naming
    case_name = [('%d'%i).zfill(len('%d'%(n_cases-1))) for i in range(n_cases)]
    if namebase:
        case_name = [namebase+'_'+caseid for caseid in case_name]

    return case_name

def convert_str(val):
    def try_type(val, data_type):
        try:
            data_type(val)
            return True
        except:
            return False
#        return isinstance(val, data_type)  ### this doesn't work b/c of numpy data types; they're not instances of base types
    def try_list(val):
        try:
            val[0]
            return True
        except:
            return False

    if try_type(val, int) and int(val) == float(val):
        return int(val)
    elif try_type(val, float):
        return float(val)
    elif val=='True':
        return True
    elif val=='False':
        return False
    # elif type(val)!=str and try_list(val):
    #     return ", ".join(['{:}'.format(i) for i in val])
    else:
        return val

def CaseGen_General(case_inputs, dir_matrix='', namebase='', save_matrix=True):
    """ Cartesian product to enumerate over all combinations of set of variables that are changed together"""

    # put case dict into lists
    change_vars = sorted(case_inputs.keys())
    change_vals = [case_inputs[var]['vals'] for var in change_vars]
    change_group = [case_inputs[var]['group'] for var in change_vars]

    # find number of groups and length of groups
    group_set = list(set(change_group))
    group_len = [len(change_vals[change_group.index(i)]) for i in group_set]

    # case matrix, as indices
    group_idx = [range(n) for n in group_len]
    matrix_idx = list(itertools.product(*group_idx))

    # index of each group
    matrix_group_idx = [np.where([group_i == group_j for group_j in change_group])[0].tolist() for group_i in group_set]

    # build final matrix of variable values
    matrix_out = []
    for i, row in enumerate(matrix_idx):
        row_out = [None]*len(change_vars)
        for j, val in enumerate(row):
            for g in matrix_group_idx[j]:
                row_out[g] = change_vals[g][val]
        matrix_out.append(row_out)
    try:
        matrix_out = np.asarray(matrix_out, dtype=str)
    except:
        matrix_out = np.asarray(matrix_out)
    n_cases = np.shape(matrix_out)[0]

    # case naming
    case_name = case_naming(n_cases, namebase=namebase)
    
    # Save case matrix
    if save_matrix:
        if not dir_matrix:
            dir_matrix = os.getcwd()
        try:
            save_case_matrix(matrix_out, change_vars, dir_matrix)
            save_case_matrix_yaml(matrix_out, change_vars, dir_matrix, case_name)
        except: 
            save_case_matrix_yaml(matrix_out, change_vars, dir_matrix, case_name)

    case_list = []
    for i in range(n_cases):
        case_list_i = {}
        for j, var in enumerate(change_vars):
            case_list_i[var] = convert_str(matrix_out[i,j])
        case_list.append(case_list_i)


    return case_list, case_name
