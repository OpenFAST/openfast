import numpy as np
import pandas as pd
import os
try:
    from .file import File, WrongFormatError, BrokenFormatError
except:
    File = dict
    class WrongFormatError(Exception): pass
    class BrokenFormatError(Exception): pass

# --------------------------------------------------------------------------------}
# --- OUT FILE 
# --------------------------------------------------------------------------------{
class FLEXOutFile(File):

    @staticmethod
    def defaultExtensions():
        return ['.res','.int']

    @staticmethod
    def formatName():
        return 'FLEX output file'

    def _read(self):
        # --- First read the binary file
        dtype=np.float32; # Flex internal data is stored in single precision
        try:
            self.data,self.tmin,self.dt,self.Version,self.DateID,self.title=read_flex_res(self.filename, dtype=dtype)
        except WrongFormatError as e:    
            raise WrongFormatError('FLEX File {}: '.format(self.filename)+'\n'+e.args[0])
        self.nt       = np.size(self.data,0)
        self.nSensors = np.size(self.data,1)
        self.time = np.arange(self.tmin, self.tmin +  self.nt * self.dt, self.dt).reshape(self.nt,1).astype(dtype)

        # --- Then the sensor file
        parentdir = os.path.dirname(self.filename)
        basename = os.path.splitext(os.path.basename(self.filename))[0]
        #print(parentdir)
        #print(basename)
        PossibleFiles=[]
        PossibleFiles+=[os.path.join(parentdir, basename+'.Sensor')]
        PossibleFiles+=[os.path.join(parentdir, 'Sensor_'+basename)]
        PossibleFiles+=[os.path.join(parentdir, 'sensor')]
        # We try allow for other files
        Found =False
        for sf in PossibleFiles:
            if os.path.isfile(sf):
                self.sensors=read_flex_sensor(sf)
                if len(self.sensors['ID'])!=self.nSensors:
                    Found = False
                else:
                    Found = True
                    break
        if not Found:
            # we are being nice and create some fake sensors info
            self.sensors=read_flex_sensor_fake(self.nSensors)

        if len(self.sensors['ID'])!=self.nSensors:
            raise BrokenFormatError('Inconsistent number of sensors: {} (sensor file) {} (out file), for file: {}'.format(len(self.sensors['ID']),self.nSensors,self.filename))

    #def _write(self): # TODO
    #    pass

    def __repr__(self):
        return 'Flex Out File: {}\nVersion:{} - DateID:{} - Title:{}\nSize:{}x{} - tmin:{} - dt:{}]\nSensors:{}'.format(self.filename,self.Version,self.DateID,self.title,self.nt,self.nSensors,self.tmin,self.dt,self.sensors['Name'])

    def _toDataFrame(self):
        # Appending time to form the dataframe
        names = ['Time'] + self.sensors['Name']
        units = ['s']    + self.sensors['Unit']
        units = [u.replace('(','').replace(')','').replace('[','').replace(']','') for u in units]
        data  = np.concatenate((self.time, self.data), axis=1)
        cols=[n+'_['+u+']' for n,u in zip(names,units)]
        return pd.DataFrame(data=data,columns=cols)

# --------------------------------------------------------------------------------}
# --- Helper Functions 
# --------------------------------------------------------------------------------{
def read_flex_res(filename, dtype=np.float32):
    # Read flex file
    with open(filename,'rb') as fid:
        #_ = struct.unpack('i', fid.read(4)) # Dummy
        _ = np.fromfile(fid, 'int32', 1) # Dummy
        # --- Trying to get DateID
        fid.seek(4) # 
        DateID=np.fromfile(fid, 'int32', 6)
        if DateID[0]<32 and DateID[1]<13 and DateID[3]<25 and DateID[4]<61:
            # OK, DateID was present
            title  = fid.read(40).strip()
        else:
            fid.seek(4) # 
            DateID = np.fromfile(fid, 'int32', 1)
            title  = fid.read(60).strip()
        _ = np.fromfile(fid, 'int32', 2) # Dummy
        # FILE POSITION <<< fid.seek(4 * 19) 
        nSensors = np.fromfile(fid, 'int32', 1)[0] 
        IDs = np.fromfile(fid, 'int32', nSensors)
        _ = np.fromfile(fid, 'int32', 1) # Dummy
        # FILE POSITION <<< fid.seek(4*nSensors+4*21)
        Version = np.fromfile(fid, 'int32', 1)[0] 
        # FILE POSITION <<< fid.seek(4*(nSensors)+4*22)
        if Version == 12:
            raise NotImplementedError('Flex out file with version 12, TODO. Implement it!')
            # TODO
            #fseek(o.fid,4*(21+o.nSensors),-1);% seek to the data from beginning of file
            #RL=o.nSensors+5; % calculate the length of each row
            #A = fread(o.fid,[RL,inf],'single'); % read whole file
            #t=A(2,:);% time vector contained in row 2
            #o.SensorData=A(5:end,:);
            # save relevant information 
            #o.tmin = t(1)     ;
            #o.dt   = t(2)-t(1);
            #o.t    = t        ;
            #o.nt   = length(t);
        elif Version in [0,2,3]:
            tmin = np.fromfile(fid, 'f', 1)[0] # Dummy
            dt = np.fromfile(fid, 'f', 1)[0] # Dummy
            scale_factors = np.fromfile(fid, 'f', nSensors).astype(dtype)
        # --- Reading Time series
        # FILE POSITION <<< fid.seek(8*nSensors + 48*2)
        data = np.fromfile(fid, 'int16').astype(dtype) #data = np.fromstring(fid.read(), 'int16').astype(dtype)
        nt   = int(len(data) / nSensors)
        try:
            if Version ==3:
                data = data.reshape(nSensors, nt).transpose()
            else:
                data = data.reshape(nt, nSensors)
        except ValueError:
            raise WrongFormatError("Flat data length {} is not compatible with {}x{} (nt x nSensors)".format(len(data),nt,nSensors))
        for i in range(nSensors):
            data[:, i] *= scale_factors[i]

        return (data,tmin,dt,Version,DateID,title)


def read_flex_sensor(sensor_file):
    with open(sensor_file, 'r') as fid:
        sensor_info_lines = fid.readlines()[2:]
    sensor_info = []
    d=dict({ 'ID':[],'Gain':[],'Offset':[],'Unit':[],'Name':[],'Description':[]});
    for line in sensor_info_lines:
        line   = line.strip().split()
        d['ID']          .append(int(line[0]))
        d['Gain']        .append(float(line[1]))
        d['Offset']      .append(float(line[2]))
        d['Unit']        .append(line[5])
        d['Name']        .append(line[6])
        d['Description'] .append(' '.join(line[7:]))
    return d
 
def read_flex_sensor_fake(nSensors):
    d=dict({ 'ID':[],'Gain':[],'Offset':[],'Unit':[],'Name':[],'Description':[]});
    for i in range(nSensors):
        d['ID']          .append(i+1)
        d['Gain']        .append(1.0)
        d['Offset']      .append(0.0)
        d['Unit']        .append('(NA)')
        d['Name']        .append('S{:04d}'.format(i+1))
        d['Description'] .append('NA')
    return d






# def write_flex_file(filename,data,tmin,dt):
#     ds = dataset
#     # Write int data file
#     f = open(filename, 'wb')
#     f.write(struct.pack('ii', 0, 0))  # 2x empty int
#     title = ("%-60s" % str(ds.name)).encode()
#     f.write(struct.pack('60s', title))  # title
#     f.write(struct.pack('ii', 0, 0))  # 2x empty int
#     ns = len(sensors)
#     f.write(struct.pack('i', ns))
#     f.write(struct.pack('i' * ns, *range(1, ns + 1)))  # sensor number
#     f.write(struct.pack('ii', 0, 0))  # 2x empty int
#     time = ds.basis_attribute()
#     f.write(struct.pack('ff', time[0], time[1] - time[0]))  # start time and time step
# 
#     scale_factors = np.max(np.abs(data), 0) / 32000
#     f.write(struct.pack('f' * len(scale_factors), *scale_factors))
#     # avoid dividing by zero
#     not0 = np.where(scale_factors != 0)
#     data[:, not0] /= scale_factors[not0]
#     #flatten and round
#     data = np.round(data.flatten()).astype(np.int16)
#     f.write(struct.pack('h' * len(data), *data.tolist()))
#     f.close()
# 
#     # write sensor file
#     f = open(os.path.join(os.path.dirname(filename), 'sensor'), 'w')
#     f.write("Sensor list for %s\n" % filename)
#     f.write(" No   forst  offset  korr. c  Volt    Unit   Navn    Beskrivelse------------\n")
#     sensorlineformat = "%3s  %.3f   %.3f      1.00     0.00 %7s %-8s %s\n"
# 
#     if isinstance(ds, FLEX4Dataset):
#         gains = np.r_[ds.gains[1:], np.ones(ds.shape[1] - len(ds.gains))]
#         offsets = np.r_[ds.offsets[1:], np.zeros(ds.shape[1] - len(ds.offsets))]
#         sensorlines = [sensorlineformat % ((nr + 1), gain, offset, att.unit[:7], att.name.replace(" ", "_")[:8], att.description[:512]) for nr, att, gain, offset in zip(range(ns), sensors, gains, offsets)]
#     else:
#         sensorlines = [sensorlineformat % ((nr + 1), 1, 0, att.unit[:7], att.name.replace(" ", "_")[:8], att.description[:512]) for nr, att in enumerate(sensors)]
#     f.writelines(sensorlines)
#     f.close()
