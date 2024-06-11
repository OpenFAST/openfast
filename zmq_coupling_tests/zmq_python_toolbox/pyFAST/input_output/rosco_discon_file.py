""" 
Input/output class for the fileformat ROSCO DISCON file
"""
import numpy as np
import pandas as pd
import os
from collections import OrderedDict

try:
    from .file import File, WrongFormatError, BrokenFormatError
except:
    File=OrderedDict
    EmptyFileError    = type('EmptyFileError', (Exception,),{})
    WrongFormatError  = type('WrongFormatError', (Exception,),{})
    BrokenFormatError = type('BrokenFormatError', (Exception,),{})

class ROSCODISCONFile(File):
    """ 
    Read/write a ROSCO DISCON file. The object behaves as a dictionary.
    
    Main methods
    ------------
    - read, write, toDataFrame, keys
    
    Examples
    --------
        f = ROSCODISCONFile('DISCON.IN')
        print(f.keys())
        print(f.toDataFrame().columns)  
    
    """

    @staticmethod
    def defaultExtensions():
        """ List of file extensions expected for this fileformat"""
        return ['.in']

    @staticmethod
    def formatName():
        """ Short string (~100 char) identifying the file format"""
        return 'ROSCO DISCON file'

    @staticmethod
    def priority(): return 60 # Priority in weio.read fileformat list between 0=high and 100:low


    def __init__(self, filename=None, **kwargs):
        """ Class constructor. If a `filename` is given, the file is read. """
        self.filename = filename
        if filename:
            self.read(**kwargs)

    def read(self, filename=None, **kwargs):
        """ Reads the file self.filename, or `filename` if provided """
        
        # --- Standard tests and exceptions (generic code)
        if filename:
            self.filename = filename
        if not self.filename:
            raise Exception('No filename provided')
        if not os.path.isfile(self.filename):
            raise OSError(2,'File not found:',self.filename)
        if os.stat(self.filename).st_size == 0:
            raise EmptyFileError('File is empty:',self.filename)
        # --- Calling (children) function to read
        _, comments, lineKeys = read_DISCON(self.filename, self)
        self.comments=comments
        self.lineKeys=lineKeys

    def write(self, filename=None):
        """ Rewrite object to file, or write object to `filename` if provided """
        if filename:
            self.filename = filename
        if not self.filename:
            raise Exception('No filename provided')
        with open(self.filename, 'w') as f:
            f.write(self.toString())
        

    def toDataFrame(self):
        """ Returns object into one DataFrame, or a dictionary of DataFrames"""
        dfs={}
        low_keys = [s.lower() for s in self.keys()]
        if 'pc_gs_n' in low_keys:
            M = np.column_stack([self['PC_GS_angles']*180/np.pi, self['PC_GS_KP'], self['PC_GS_KI'], self['PC_GS_KD'], self['PC_GS_TF']] )
            cols = ['Pitch_[deg]', 'KP_[-]', 'KI_[s]', 'KD_[1/s]',  'TF_[-]']
            dfs['PitchSchedule'] = pd.DataFrame(data=M, columns=cols)
        if 'ps_bldpitchmin_n' in low_keys:
            M = np.column_stack([self['PS_WindSpeeds'], self['PS_BldPitchMin']])
            cols = ['WindSpeed_[m/s]', 'Pitch_[deg]']
            dfs['PitchSaturation'] = pd.DataFrame(data=M, columns=cols)
        if 'prc_n' in low_keys:
            M = np.column_stack([self['PRC_WindSpeeds'], self['PRC_RotorSpeeds']*30/np.pi])
            cols = ['WindSpeed_[m/s]', 'RotorSpeed_[rpm]']
            dfs['PowerTracking'] = pd.DataFrame(data=M, columns=cols)

        return dfs

    # --- Optional functions
    def __repr__(self):
        """ String that is written to screen when the user calls `print()` on the object. 
        Provide short and relevant information to save time for the user. 
        """
        s='<{} object>:\n'.format(type(self).__name__)
        s+='|Main attributes:\n'
        s+='| - filename: {}\n'.format(self.filename)
        # --- Example printing some relevant information for user
        #s+='|Main keys:\n'
        #s+='| - ID: {}\n'.format(self['ID'])
        #s+='| - data : shape {}\n'.format(self['data'].shape)
        s+='|Main methods:\n'
        s+='| - read, write, toDataFrame, keys'
        return s
    
    def toString(self):
        """ """
        maxKeyLengh = np.max([len(k) for k in self.keys()])
        maxKeyLengh = max(maxKeyLengh, 18)
        fmtKey = '{:' +str(maxKeyLengh)+'s}'
        s=''
        for l in self.lineKeys:
            if len(l)==0:
                s+='\n'
            elif l.startswith('!'):
                s+=l+'\n'
            else:
                param = l
                comment = self.comments[param]
                v = self[param]
                sparam = '! '+fmtKey.format(param)
                # NOTE: could to "param" specific outputs here
                FMTs = {}
                FMTs['{:<4.6f}']=['F_NotchBetaNumDen', 'F_FlCornerFreq', 'F_FlpCornerFreq', 'PC_GS_angles', 'PC_GS_KP', 'PC_GS_KI', 'PC_GS_KD', 'PC_GS_TF', 'IPC_Vramp', 'IPC_aziOffset']
                FMTs['{:<4.3e}']=['IPC_KP','IPC_KI']
                FMTs['{:<4.3f}']=['PRC_WindSpeeds', 'PRC_RotorSpeed','PS_WindSpeeds']
                FMTs['{:<4.4f}']=['WE_FOPoles_v']
                FMTs['{:<10.8f}']=['WE_FOPoles']
                FMTs['{:<10.3f}']=['PS_BldPitchMin']
                fmtFloat='{:<014.5f}'
                for fmt,keys in FMTs.items():
                    if param in keys:
                        fmtFloat=fmt
                        break
                if type(v) is str:
                    sval='"{:15s}"    '.format(v)
                elif hasattr(v, '__len__'):
                    if isinstance(v[0], (np.floating, float)):
                        sval=' '.join([fmtFloat.format(vi) for vi in v]  )+'    '
                    else:
                        sval=' '.join(['{}'.format(vi) for vi in v]  )+'    '
                elif type(v) is int:
                    sval='{:<14d}      '.format(v)
                elif isinstance(v, (np.floating, float)):
                    sval=fmtFloat.format(v) + '     '
                else:
                    sval='{} '.format(v)
                s+='{}{}{}\n'.format(sval, sparam, comment)
        return s






# Some useful constants
pi = np.pi
rad2deg = np.rad2deg(1)
deg2rad = np.deg2rad(1)
rpm2RadSec = 2.0*(np.pi)/60.0
RadSec2rpm = 60/(2.0 * np.pi)

def write_DISCON(turbine, controller, param_file='DISCON.IN', txt_filename='Cp_Ct_Cq.txt', rosco_vt = {}):
    """
    Print the controller parameters to the DISCON.IN input file for the generic controller

    Parameters:
    -----------
    turbine: class
                Turbine class containing turbine operation information (ref speeds, etc...)
    controller: class
                Controller class containing controller operation information (gains, etc...)
    param_file: str, optional
        filename for parameter input file, should be DISCON.IN
    txt_filename: str, optional
                    filename of rotor performance file
    """

    # Get ROSCO var tree if not provided
    if not rosco_vt:
        rosco_vt = DISCON_dict(turbine, controller, txt_filename)

    print('Writing new controller parameter file parameter file: %s.' % param_file)
    # Should be obvious what's going on here...
    file = open(param_file,'w')

    # Write Open loop input
    if rosco_vt['OL_Mode'] and hasattr(controller, 'OpenLoop'):
        write_ol_control(controller)

def read_DISCON(DISCON_filename, DISCON_in = None):
    '''
    Read the DISCON input file.
    Adapted from ROSCO_Toolbox, https:github.com/NREL/ROSCO

    Parameters:
    ----------
    DISCON_filename: string
        Name of DISCON input file to read
    
    Returns:
    --------
    DISCON_in: Dict
        Dictionary containing input parameters from DISCON_in, organized by parameter name
    '''
    
    if DISCON_in is None:
        DISCON_in = OrderedDict()
    comments={}
    lineKeys=[]
    with open(DISCON_filename) as discon:
        for line in discon:
            line=line.strip()
            # empty lines
            if len(line)==0:
                lineKeys.append('')
                continue
            # Pure comments
            if line[0] == '!':
                lineKeys.append(line)
                continue

            if (line.split()[1] != '!'):    # Array valued entries
                sps = line.split()
                array_length = sps.index('!')
                param = sps[array_length+1]
                values = np.array( [float(x) for x in sps[:array_length]] )
            else:                           # All other entries
                param = line.split()[2]
                value = line.split()[0]
                # Remove printed quotations if string is in quotes
                if (value[0] == '"') or (value[0] == "'"):
                    values = value[1:-1]
                else:
                    if value.find('.')>0:
                        values = float(value)
                    else:
                        values = int(value)
            DISCON_in[param] = values
            lineKeys.append(param)

            sp = line.split('!')
            comment = sp[1].strip()
            comment = comment[len(param):].strip()
            comments [param] = comment

    return DISCON_in, comments, lineKeys


if __name__ == '__main__':
    filename = 'DISCON.in'
    rd = ROSCODISCONFile(filename)
    #print(rd.keys())
#     print(rd.toString())
    rd.write(filename+'_WEIO')
    print(rd.toDataFrame())
