import numpy as np
import pandas as pd
import os
try:
    from .file import File, WrongFormatError, BrokenFormatError
except:
    File = dict
    class WrongFormatError(Exception): pass
    class BrokenFormatError(Exception): pass

class ProfileSet():
    def __init__(self,header,thickness,polars,polar_headers):
        self.header      = header
        self.polars      = polars
        self.thickness   = thickness
        self.polar_headers = polar_headers

    def toString(self,PREFIX=''):
        s =PREFIX+self.header+'\n'
        s+=' '.join([str(t) for t in self.thickness])+'\n'
        s+=str(self.polars[0].shape[0])+'\n'
        for ph,t,polar in zip(self.polar_headers,self.thickness,self.polars):
            s+=ph+'\n'
            s+='\n'.join([' '.join(['{:15.7e}'.format(v) for v in line]) for line in polar])
#             s+=ph+'\n'
        return s

    def __repr__(self):
        s ='Class ProfileSet (attributes: header, polars, thickness, polar_headers)\n'
        s+='   header         : '+self.header+'\n'
        s+='   thickness       : '+str(self.thickness)+'\n'
        s+='   Number of polars: '+str(len(self.thickness))+'\n'
        s+='   Alpha values    : '+str(self.polars[0].shape[0])+'\n'
        for ip,(ph,t) in enumerate(zip(self.polar_headers,self.thickness)):
            s+= '   Polar: {}, Thickness: {}, Header: {}\n'.format(ip+1,t,ph)
        return s

class FLEXProfileFile(File):

    @staticmethod
    def defaultExtensions():
        return ['.pro','.00X'] #'.001 etc..'

    @staticmethod
    def formatName():
        return 'FLEX profile file'

    def _read(self):
        self.ProfileSets=[]
        setNumber=1
        with open(self.filename, 'r', errors="surrogateescape") as f:
            def read_header(allow_empty=False):
                """ Reads the header of a profile set (4 first lines)
                     - The first line may start with "Profile set I:" to indicate a set number
                     - Second line is number of thicknesses
                     - Third is thicnkesses
                     - Fourth is number of alpha values
                """
                header=[]
                for i, line in enumerate(f):
                    header.append(line.strip())
                    if i==3:
                        break
                if len(header)<4:
                    if allow_empty:
                        return [],[],'',False
                    else:
                        raise WrongFormatError('A Flex profile file needs at leats 4 lines of headers')
                try:
                    nThickness=int(header[1])
                except:
                    raise WrongFormatError('Number of thicknesses (integer) should be on line 2')
                try:
                    thickness=np.array(header[2].split()).astype(float)
                except:
                    raise WrongFormatError('Number of thicknesses (integer) should be on line 2')
                if len(thickness)!=nThickness:
                    raise WrongFormatError('Number of thicknesses read ({}) different from the number reported ({})'.format(len(thickness),nThickness))
                try:
                    nAlpha=int(header[3])
                except:
                    raise WrongFormatError('Number of alpha values (integer) should be on line 4')
                if header[0].lower().find('profile set')==0:
                    header[0]=header[0][11:]
                    bHasSets=True
                else:
                    bHasSets=False
                return nAlpha,thickness,header[0],bHasSets

            def read_polars(nAlpha,thickness):
                polars=[]
                polar_headers=[]
                for it,t in enumerate(thickness):
                    polar_headers.append(f.readline().strip())
                    polars.append(np.zeros((nAlpha,4)))
                    try:
                        for ia in range(nAlpha):
                            polars[it][ia,:]=np.array([f.readline().split()]).astype(float)
                    except:
                        raise BrokenFormatError('An error occured while reading set number {}, polar number {}, (thickness {}), value number {}.'.format(setNumber,it+1,t,ia+1))

                return polars,polar_headers

            # Reading headers and polars
            while True:
                nAlpha,thickness,Header,bHasSets = read_header(allow_empty=setNumber>1)
                if len(thickness)==0: 
                    break
                polars,polar_headers             = read_polars(nAlpha,thickness)
                PSet= ProfileSet(Header,thickness,polars,polar_headers)
                self.ProfileSets.append(PSet)
                setNumber=setNumber+1

    def toString(self):
        s=''
        if len(self.ProfileSets)>0:
            prefix='PROFILE SET '
        else:
            prefix=''
        for pset in self.ProfileSets:
            s+=pset.toString(prefix)
        return s

    def _write(self):
        with open(self.filename,'w') as f:
            f.write(self.toString)

    def __repr__(self):
        s ='Class FlexProfileFile (attributes: ProfileSets)\n'
        s+='   Number of profiles sets: '+str(len(self.ProfileSets))+'\n'
        for ps in self.ProfileSets:
            s+=ps.__repr__()
        return s


    def _toDataFrame(self):
        cols=['Alpha_[deg]','Cl_[-]','Cd_[-]','Cm_[-]']
        dfs = {}
        for iset,pset in enumerate(self.ProfileSets):
            for ipol,(thickness,polar) in enumerate(zip(pset.thickness,pset.polars)):
                name='pc_set_{}_t_{}'.format(iset+1,thickness)
                dfs[name] = pd.DataFrame(data=polar, columns=cols)
        return dfs

