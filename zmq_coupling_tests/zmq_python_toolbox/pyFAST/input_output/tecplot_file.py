""" 
Read/Write TecPto ascii files
sea read_tecplot documentation below

Part of weio library: https://github.com/ebranlard/weio

"""
import pandas as pd
import numpy as np
import os
import struct

try:
    from .file import File, EmptyFileError, WrongFormatError, BrokenFormatError
except:
    EmptyFileError = type('EmptyFileError', (Exception,),{})
    WrongFormatError =  type('WrongFormatError', (Exception,),{})
    BrokenFormatError =  type('BrokenFormatError', (Exception,),{})
    File=dict



Keywords=['title','variables','zone','text','geometry','datasetauxdata','customlabels','varauxdata']
# --------------------------------------------------------------------------------}
# --- Helper functions 
# --------------------------------------------------------------------------------{
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
 
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
 
    return False


def _process_merged_line(line, section, dict_out):
    n = len(section)
    line = line[n:].strip()
    if section=='title':
        dict_out[section]=line
    elif section=='variables':
        line = line.replace('=','').strip()
        line = line.replace(',',' ').strip()
        line = line.replace('  ',' ').strip()
        line = line.replace('[','_[').strip()
        line = line.replace('(','_(').strip()
        line = line.replace('__','_').strip()
        if line.find('"')==0:
            line = line.replace('" "',',')
            line = line.replace('"','')
            sp=line.split(',')
        else:
            sp=line.split()
        dict_out[section]=sp
    elif section=='datasetauxdata':
        if section not in dict_out.keys():
            dict_out[section]={} # initialixe an empty directory
        sp    = line.split('=')
        key   = sp[0]
        value = sp[1].replace('"','').strip()
        if is_number(value):
            value=float(value)
        dict_out[section][key]=value

    elif section=='zone':
        if section not in dict_out.keys():
            dict_out[section]={} # initialixe an empty directory
        sp    = line.split('=')
        key   = sp[0]
        value = sp[1].replace('"','').strip()
        if is_number(value):
            value=float(value)
        dict_out[section][key]=value
        
    else:
        print('!!! Reading of section not implemented:')
        print('Processing section {}:'.format(section),line)
        dict_out[section]=line

def read_tecplot(filename, dict_out={}):
    """ Reads a tecplot file
    Limited support:
       - title optional
       - variables mandatory
       - Lines may be continued to next line, stopping when a predefined keyword is detected
    For now, assumes that only one section of numerical data is present
    """

    merged_line=''
    current_section=''
    variables=[]
    with open(filename, "r") as f:
        dfs = []  # list of dataframes
        iline=0
        while True:
            line= f.readline().strip()
            iline+=1
            if not line:
                break
            l=line.lower().strip()
            # Comment
            if l[0]=='#':
                continue
            new_section = [k for k in Keywords if l.find(k)==0 ]

            if len(new_section)==1:
                # --- Start of a new section
                # First, process the previous section
                if len(merged_line)>0: 
                    _process_merged_line(merged_line, current_section, dict_out)
                # Then start the new section
                current_section=new_section[0]
                merged_line =line
            elif len(current_section)==0:
                raise WrongFormatError('No section detected')
            else:
                if current_section=='title' or current_section=='variables':
                    # OK
                    pass
                else:
                    if 'variables' not in dict_out.keys():
                        raise WrongFormatError('The `variables` section should be present')
                sp = l.split()
                if is_number(sp[0]):
                    if len(merged_line)>0: 
                        _process_merged_line(merged_line, current_section, dict_out)
                    # --- Special case of numerical values outside of zone
                    f.close()
                    M = np.loadtxt(filename, skiprows = iline-1)
                    if M.shape[1]!=len(dict_out['variables']):
                        raise BrokenFormatError('Number of columns of data does not match number of variables')
                    dict_out['data']=M
                    break
                else:
                    # --- Continuation of previous section
                    merged_line +=' '+line
    return dict_out


class TecplotFile(File):

    @staticmethod
    def defaultExtensions():
        return ['.dat']

    @staticmethod
    def formatName():
        return 'Tecplot ASCII file'

    def __init__(self,filename=None,**kwargs):
        self.filename = None
        if filename:
            self.read(filename=filename,**kwargs)

    def read(self, filename=None):
        """ read a tecplot ascii file
        sea `read_tecplot` documentation above
        """
        if filename:
            self.filename = filename
        if not self.filename:
            raise Exception('No filename provided')
        if not os.path.isfile(self.filename):
            raise OSError(2,'File not found:',self.filename)
        if os.stat(self.filename).st_size == 0:
            raise EmptyFileError('File is empty:',self.filename)

        try:
            read_tecplot(filename,self)
        except BrokenFormatError:
            raise 
        except WrongFormatError:
            raise 
        except Exception as e:    
            raise WrongFormatError('Tecplot dat File {}: '.format(self.filename)+e.args[0])

    def write(self, filename=None, precision=None):
        """ Write tecplot ascii file """
        if filename:
            self.filename = filename
        if not self.filename:
            raise Exception('No filename provided')

        with open(self.filename, mode='w') as f:            
            if 'title' in self.keys():
                f.write('TITLE = {}\n'.format(self['title']))
            f.write('VARIABLES = ' + ','.join(['"{}"'.format(col) for col in self['variables'] ]) + '\n')
            for k in Keywords[2:]:
                if k in self.keys():
                    f.write('{} = {}\n'.format(k,self[k]))
            # Data
            if 'data' in self.keys():
                for row in self['data']:
                    srow = np.array2string(row, edgeitems=0, separator=' ', precision=precision)
                    f.write(srow[1:-1]+'\n')


    def __repr__(self):
        s='<{} object> with keys:\n'.format(type(self).__name__)
        for k,v in self.items():
            s+=' - {}: {}\n'.format(k,v)
        return s

    def toDataFrame(self):
        return pd.DataFrame(data=self['data'],columns=self['variables'])

if __name__=='__main__':
    mb = MannBoxFile('mann_bin/mini-u.bin', N=(2,4,8))
    F1=mb['field'].ravel()
    mb.write('mann_bin/mini-u-out.bin')

    mb2= MannBoxFile('mann_bin/mini-u-out.bin', N=(2,4,8))
    F2=mb2['field'].ravel()
#     print(F1-F2)
