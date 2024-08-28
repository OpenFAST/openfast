import os
import xlrd

""" 
Script to generate the Output channel dictionaries.  Included so they can be easily updated
for future FAST versions and additional FAST modules.

The script did not work for me with a xlsx worksheet, but it did with a legacy xls file

by Evan Gaertner, maintained by Pietro Bortolotti @ NREL on June 2nd, 2021
"""


def GetOutlistParameters(fname_vars_out, xl_files_in, sheet_list, write_mode, final_dict):

    f = open(fname_vars_out, write_mode)
    f.write('""" Generated from FAST OutListParameters.xlsx files with AeroelasticSE/src/AeroelasticSE/Util/create_output_vars.py """\n')

    for fname, sheet_name in zip(xl_files_in, sheet_list):
        book = xlrd.open_workbook(fname)
        sheet_idx = book.sheet_names().index(sheet_name)
        sheet = book.sheet_by_index(sheet_idx)

        print(sheet_name, '\t:\t', fname)

        # get data from sheet
        outlist = {}
        headers = {}
        var_names = []
        n = sheet.nrows
        n_dup = 0

        for idx in range(1, sheet.nrows):
            row = sheet.row_values(idx)
            if row[1]!='' and row[5]!='':
                outvar = str(row[1])
                var_names.append(str(row[1]))
                outlist[outvar] = {}
                outlist[outvar]['desc'] = ''.join([i if ord(i) < 128 else ' ' for i in row[3]])
                outlist[outvar]['conv'] = ''.join([i if ord(i) < 128 else ' ' for i in row[4]])
                outlist[outvar]['unit'] = ''.join([i if ord(i) < 128 else ' ' for i in row[5]])
                outlist[outvar]['idx'] = idx

                for var in row[2].split(','):
                    if var.strip() != '':
                        outlist[str(var.strip())] = outlist[outvar]
                        var_names.append(str(var.strip()))
                        n += 1
                        n_dup +=1

            else:
                headers[idx+n_dup] = str(row[0])
        
        # write out variables + metadata
        f.write('\n\n""" ' + sheet_name + ' """\n')
        f.write(sheet_name + ' = {}\n')

        var_count = 0
        for idx in range(1, n):
            if idx in headers.keys():
                f.write('\n# ' + headers[idx] + '\n')
            else:
                var = var_names[var_count]
                f.write(("%s['%s']"%(sheet_name, var)).ljust(30) + "= False     # %s; %s; %s\n"%(outlist[var]['unit'], outlist[var]['desc'], outlist[var]['conv']))
                var_count += 1

    # combine in one dict
    f.write('\n\n""" Final Output Dictionary """\n')
    f.write('%s = {}\n'%final_dict)
    for var in sheet_list:
        f.write(("%s['%s']"%(final_dict, var)).ljust(30) + "= %s\n"%var)


    f.close()


if __name__=="__main__":

    try:
        # File Location
        root_dir = os.path.dirname(os.path.dirname( os.path.dirname( ( os.path.realpath(__file__) ) ) )) + os.sep
        outlist_fast_lib = os.path.join(root_dir, 'docs', 'OtherSupporting' , 'OutListParameters.xls')
        # Sheets to grab
        sheet_list = [
                        'AeroDyn', 
                        'BeamDyn', 
                        'ElastoDyn', 
                        'InflowWind', 
                        'ServoDyn', 
                        'HydroDyn', 
                        'Morison', 
                        'SeaState', 
                        'SubDyn',
                        'AeroDyn_Nodes',
                        'BeamDyn_Nodes',
                        'ElastoDyn_Nodes',
                        'AeroDisk',
                        'SimpleElastoDyn',
                        ] 
        xl_files = [outlist_fast_lib]*len(sheet_list)
        # Output naming
        fname_vars_out = 'FAST_vars_out.py'
        # Run
        GetOutlistParameters(fname_vars_out, xl_files, sheet_list, 'w', 'FstOutput')
    except:
        raise Exception("Try saving the OpenFAST xlsx file into a legacy xls file")