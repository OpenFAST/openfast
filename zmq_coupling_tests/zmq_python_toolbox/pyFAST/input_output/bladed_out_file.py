import os
import numpy as np
import re
import pandas as pd
import glob
import shlex
try:
    from .file import File, WrongFormatError, BrokenFormatError
except:
    File = dict
    class WrongFormatError(Exception): pass
    class BrokenFormatError(Exception): pass

        
# --------------------------------------------------------------------------------}
# --- Helper functions 
# --------------------------------------------------------------------------------{
def read_bladed_sensor_file(sensorfile):        
    """ 
    Extract relevant informations from a bladed sensor file
    """
    with open(sensorfile, 'r') as fid:
        sensorLines = fid.readlines()

    dat=dict() # relevant info in sensor file
    
    ## read sensor file line by line (just read up to line 20)
    #while i < 17: 
    for i, t_line in enumerate(sensorLines):
        if i>30:
            break
        t_line = t_line.replace('\t',' ')
        
        if t_line.startswith('NDIMENS'):
            # check what is matrix dimension of the file. For blade & tower,  
            # the matrix is 3-dimensional. 
            temp =  t_line[7:].strip().split()
            dat['NDIMENS'] = int(temp[-1]);

        elif t_line.startswith('DIMENS'):
            # check what is the size of the matrix
            # for example, it can be 11x52500 or 12x4x52500            
            temp =  t_line[6:].strip().split()
            dat['nSensors'] = int(temp[0])
            dat['nMajor']  = int(temp[dat['NDIMENS']-1])
            if dat['NDIMENS'] == 2:
                dat['nSections'] = 1
                dat['SectionList'] = []

        elif t_line.startswith('FORMAT'):
            # precision: n/a, R*4, R*8, I*4 
            temp = t_line[7:].strip()
            dat['Precision'] = np.float32
            if temp[-1] == '8':
                dat['Precision'] = np.float64

        elif t_line.startswith('GENLAB'):
            # category of the file you are reading:
            dat['category'] = t_line[6:].strip().replace('\'','')

        elif t_line.startswith('AXIVAL'):
            # Section on the 3rd dimension you are reading
            # sometimes, the info is written on "AXITICK"
            temp = t_line[7:].split()
            dat['SectionList'] = np.array(temp, dtype=float)
            dat['nSections'] = len(dat['SectionList'])

        elif t_line.startswith('AXITICK'):
            # Section on the 3rd dimension you are reading
            # sometimes, the info is written on "AXIVAL"
            # Check next line, we concatenate if doesnt start with AXISLAB (Might need more cases)
            try:
                # Combine the strings into one string
                combined_string = ''.join(sensorLines)
                
                # Search everything betwee AXITICK and AXISLAB with a regex pattern
                t_line = re.search(r'(?<=AXITICK).+?(?=AXISLAB)', combined_string, flags=re.DOTALL)
                t_line=t_line.group(0)
                # Replace consecutive whitespace characters with a single space
                t_line = re.sub(r'\s+', ' ', t_line)
            except:
                pass

            temp = t_line.strip()
            temp = temp.strip('\'').split('\' \'')
            dat['SectionList'] = np.array(temp, dtype=str)
            dat['nSections'] = len(dat['SectionList'])

        elif t_line.startswith('VARIAB'):
            # channel names, NOTE: either quoted, non-quoted, and a mix of both
            # Check next line, we concatenate if doesnt start with AXISLAB (Might need more cases)
            try:
                nextLine=sensorLines[i+1].strip()
                if not nextLine.startswith('VARUNIT'):
                    t_line = t_line.strip()+' '+nextLine
            except:
                pass
            dat['ChannelName'] = shlex.split(t_line[6:])

        elif t_line.startswith('VARUNIT'):
            # channel units:         
            # Check next line, we concatenate if doesnt start with AXISLAB (Might need more cases)
            try:
                nextLine=sensorLines[i+1].strip()
                if not nextLine.startswith('AXISLAB'):
                    t_line = t_line.strip()+' '+nextLine
            except:
                pass
            def repUnits(s):
                s = s.replace('TT','s^2').replace('T','s').replace('A','rad')
                s = s.replace('P','W').replace('L','m').replace('F','N').replace('M','kg')
                return s
            dat['ChannelUnit']=[repUnits(s) for s in shlex.split(t_line[7:].strip())]

        elif t_line.startswith('MIN '):
            dat['MIN'] = float(t_line[3:].strip()) # Start time?

        elif t_line.startswith('STEP'):
            dat['STEP'] = float(t_line[4:].strip()) # DT?


    NeededKeys=['ChannelName','nSensors','nMajor','nSections']
    if not all(key in dat.keys() for key in NeededKeys):
        raise BrokenFormatError('Broken or unsupported format. Some necessary keys where not found in the bladed sensor file: {}'.format(sensorfile))

    if len(dat['ChannelName']) != dat['nSensors']:
        raise BrokenFormatError('Broken or unsupported format. Wrong number of channels while reading bladed sensor file: {}'.format(sensorfile))
        # if number of channel names are not matching with Sensor number then create dummy ones:
        #dat['ChannelName'] = ['Channel' + str(ss) for ss in range(dat['nSensors'])]

            
    return dat
            
def OrgData(data, **info):
    """ Flatten 3D field into 2D table"""
    # since some of the matrices are 3 dimensional, we want to make all 
    # to 2d matrix, so I am organizing them here:
    if info['NDIMENS'] == 3:
        SName = []
        SUnit = []
        dataOut = np.zeros( (info['nMajor'],len(info['SectionList'])*len(info['ChannelName'])) ) 
        
        col_vec = -1
        for isec,sec in enumerate(info['SectionList']):
            for ichan,(chan,unit) in enumerate(zip(info['ChannelName'], info['ChannelUnit'])):
                try:
                    SName.append(str(np.around(float(sec),2)) + 'm-' + chan)
                except ValueError:
                    SName.append(str(sec) + '-' + chan)
                SUnit.append(unit)
                col_vec +=1
                dataOut[:,col_vec] = data[:,isec,ichan]

        data = dataOut
        info['ChannelName'] = SName
        info['ChannelUnit'] = SUnit
    else:
        pass # Nothing to do for 2D

    return data, info



def read_bladed_output(sensorFilename, readTimeFilesOnly=False):
    """
    read a bladed sensor file and data file, reorganize a 3D file into 2D table
    """
    # --- Read sensor file and extract relevant informations
    sensorInfo = read_bladed_sensor_file(sensorFilename)
    nSensors   = sensorInfo['nSensors']
    nMajor     = sensorInfo['nMajor']
    nSections  = sensorInfo['nSections']
    hasTime = 'MIN' and 'STEP' in sensorInfo.keys()
    # --- Return if caller only wants time series
    if readTimeFilesOnly and not hasTime:
        return [], {}
    
    # --- Read data file
    dataFilename = sensorFilename.replace('%','$')

    if isBinary(dataFilename):            # it is binary            

        with open(os.path.join(dataFilename), 'rb') as fid_2:
            data = np.fromfile(fid_2, sensorInfo['Precision'])

        try:
            if sensorInfo['NDIMENS'] == 3:
                data = np.reshape(data,(nMajor, nSections, nSensors), order='C')

            elif sensorInfo['NDIMENS'] == 2:
                data = np.reshape(data,(nMajor,nSensors), order='C')
        except:
            print('>>> Failed to reshape binary file {}'.format(dataFilename))
            raise

            
    else:
        #print('it is ascii', NDIMENS)
        if sensorInfo['NDIMENS'] == 2:
            try:
                # Data is stored as time, signal, we reshape to signal, time
                data = np.loadtxt(dataFilename)
            except ValueError as e:
                # Most likely this was a binary file...
                data = np.empty((nMajor, nSensors)) * np.nan
                print('>>> Value error while reading 2d ascii file: {}'.format(dataFilename))
                raise e
            except:
                data = np.empty((nMajor, nSensors)) * np.nan
                print('>>> Failed to read 2d ascii file: {}'.format(dataFilename))
                raise

       
        elif sensorInfo['NDIMENS'] == 3:
            try:
                # Data is stored as sections, time, signal, we reshape to signal, section, time
                data = np.loadtxt(dataFilename).reshape((nMajor, nSections, nSensors),order='C')
            except:
                data = np.empty((nMajor, nSections, nSensors)) * np.nan
                print('>>> Failed to read 3d ascii file: {}'.format(dataFilename))

    return OrgData(data, **sensorInfo)


class BladedFile(File):
    r"""
    Read a Bladed out put file (current version is only binary files)
    
    Main methods:
        read: it finds all % and $ files based on selected .$PJ file and calls "DataValue" to read data from all those files
        toDataFrame: create Pandas dataframe output
        
    Main data stored:
         self.dataSets: dictionary of datasets, for each "length" of data

    example: 
        filename = r'h:\004_Loads\Sim\Bladed\003\Ramp_up\Bladed_out_ascii.$04'        
        f = BladedFile(filename)
        print(f.dataSets.keys())
        df = f.toDataFrame()
        
    """ 
    @staticmethod
    def defaultExtensions():
        return ['.%*', '.$*'] 

    @staticmethod
    def formatName():
        return 'Bladed output file'
    
    def __init__(self, filename=None, **kwargs):
        self.filename = filename
        if filename:
            self.read(**kwargs)
        
    def read(self, filename=None, **kwargs):
        """ read self, or read filename if provided """
        if filename:
            self.filename = filename
        if not self.filename:
            raise Exception('No filename provided')
        if not os.path.isfile(self.filename):
            raise OSError(2,'File not found:',self.filename)
        if os.stat(self.filename).st_size == 0:
            raise EmptyFileError('File is empty:',self.filename)
        # Calling children function
        self._read(**kwargs)
    
    def _read(self):
        """ 
        Read a bladed output file, data are in *.$II and sensors in *%II. 
         - If the file is a *$PJ file, all output files are read
         - Otherwise only the current file is read 
        """

        basename, ext = os.path.splitext(self.filename)
        if ext.lower()=='.$pj':
            readTimeFilesOnly=True
            searchPattern = basename + '.%[0-9][0-9]*' # find all files in the folder
        else:
            readTimeFilesOnly=False
            searchPattern = basename + ext.replace('$','%') # sensor file name
        
        # Look for files matching pattern
        files = glob.glob(searchPattern)

        # We'll store the data in "dataSets",dictionaries
        dataSets={}

        if len(files)==0:
            e= FileNotFoundError(searchPattern)
            e.filename=(searchPattern)
            raise e
        elif len(files)==1:
            readTimeFilesOnly=False

        files.sort()

        for i,filename in enumerate(files):

            dataFilename = filename.replace('%','$')
            try:
                # Call "Read_bladed_file" function to Read and store data:
                data, info = read_bladed_output(filename, readTimeFilesOnly=readTimeFilesOnly)    
            except FileNotFoundError as e:
                print('>>> Missing datafile: {}'.format(e.filename))
                if len(files)==1:
                    raise e
                continue
            except ValueError as e:
                print('>>> ValueError while reading: {}'.format(dataFilename))
                if len(files)==1:
                    raise e
                continue
            except:
                raise 
                print('>>> Misc error while reading: {}'.format(dataFilename))
                if len(files)==1:
                    raise 
                continue
            if len(data)==0:
                print('>>> Skipping file since no time present {}'.format(filename))
                continue
            
            # we use number of data as key, but we'll use "name" later
            key = info['nMajor']
            
            if key in dataSets.keys():
                # dataset with this length are already present, we concatenate
                dset = dataSets[key]
                dset['data'] =  np.column_stack((dset['data'], data))
                dset['sensors'] += info['ChannelName']
                dset['units']   += info['ChannelUnit']
                dset['name']  = 'Misc_'+str(key)

            else:
                # We add a new dataset for this length
                dataSets[key] = {}
                dset = dataSets[key]
                # We force a time vector when possible
                if 'MIN' and 'STEP' in info.keys():
                    time = np.arange(info['nMajor'])*info['STEP'] + info['MIN']
                    data = np.column_stack((time, data))
                    info['ChannelName'].insert(0, 'Time')
                    info['ChannelUnit'].insert(0, 's')

                dset['data']    = data
                dset['sensors'] = info['ChannelName']
                dset['units']   = info['ChannelUnit']
                dset['name']    = info['category']

        # Check if we have "many" misc, if only one, replace by "Misc"
        keyMisc = [k for k,v in dataSets.items() if v['name'].startswith('Misc_')]
        if len(keyMisc)==1:
            #dataSets[keyMisc[0]]['name']='Misc'
            # We keep only one dataset for simplicity
            self.dataSets= {'Misc': dataSets[keyMisc[0]]}
        else:
            # Instead of using nMajor as key, we use the "name"
            self.dataSets= {v['name']: v for (k, v) in dataSets.items()}

                
    def toDataFrame(self):        
        dfs={}
        for k,dset in self.dataSets.items():
            BL_ChannelUnit = [ name+' ['+unit+']' for name,unit in zip(dset['sensors'],dset['units'])]
            df = pd.DataFrame(data=dset['data'], columns=BL_ChannelUnit)
            # remove duplicate columns
            df = df.loc[:,~df.columns.duplicated()]
            df.columns.name = k # hack for pyDatView when one dataframe is returned
            dfs[k] = df
        if len(dfs)==1:
            return dfs[next(iter(dfs))]
        else:
            return dfs
    

def isBinary(filename):
    with open(filename, 'r') as f:
        try:
            # first try to read as string
            l = f.readline()
            # then look for weird characters
            for c in l:
                code = ord(c)
                if code<10 or (code>14 and code<31):
                    return True
            return False
        except UnicodeDecodeError:
            return True

if __name__ == '__main__':
    pass
    #filename = r'E:\Work_Google Drive\Bladed_Sims\Bladed_out_binary.$41'
    #Output = BladedFile(filename)
    #df = Output.toDataFrame()


