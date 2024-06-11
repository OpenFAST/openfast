import numpy as np
import pandas as pd
import os

try:
    from .file import File, WrongFormatError, BrokenFormatError
except:
    File = dict
    class WrongFormatError(Exception): pass
    class BrokenFormatError(Exception): pass

class FLEXBladeFile(File):

    @staticmethod
    def defaultExtensions():
        return ['.bld','.bla','.00X'] #'.001 etc..'

    @staticmethod
    def formatName():
        return 'FLEX blade file'

    def _read(self):
        headers_all = ['r_[m]','EIFlap_[Nm2]','EIEdge_[Nm2]','GKt_[Nm2]','Mass_[kg/m]','Jxx_[kg.m]','PreBendFlap_[m]','PreBendEdge_[m]'\
                    ,'Str.Twist_[deg]','PhiOut_[deg]','Ycog_[m]','Yshc_[m]','CalcOutput_[0/1]'\
                    ,'Chord_[m]','AeroTwist_[deg]','RelThickness_[%]','AeroCenter_[m]','AeroTorsion_[0/1]','ProfileSet_[#]']
        with open(self.filename, 'r', errors="surrogateescape") as f:
            try:
                firstline  = f.readline().strip()
                nSections = int(f.readline().strip().split()[0])
            except:
                raise WrongFormatError('Unable to read first two lines of blade file')
            try:
                self.version=int(firstline[1:4])
            except:
                self.version=0
            # --- Different handling depending on version
            if self.version==0:
                # Version 0 struct has no GKt
                # Version 0 aero has no profile set, no TorsionAero
                nColsStruct = 8
                nColsAero   = 5
                struct_headers = ['r_[m]','EIFlap_[Nm2]','EIEdge_[Nm2]','Mass_[kg/m]','Str.Twist_[deg]','CalcOutput_[0/1]','PreBendFlap_[m]','PreBendEdge_[m]']
                aero_headers   = ['X_BladeRoot_[m]','Chord_[m]','AeroTwist_[deg]','RelThickness_[%]','AeroCenter_[m]']
            elif self.version==1:
                # Version 1 struct has GKt
                # Version 1 aero has no profile set
                nColsStruct = 8 
                nColsAero   = 6
                struct_headers = ['r_[m]','EIFlap_[Nm2]','EIEdge_[Nm2]','Mass_[kg/m]','Str.Twist_[deg]','CalcOutput_[0/1]','PreBendFlap_[m]','PreBendEdge_[m]']
                aero_headers   = ['X_BladeRoot_[m]','Chord_[m]','AeroTwist_[deg]','RelThickness_[%]','AeroCenter_[m]','AeroTorsion_[0/1]']
            elif self.version==2:
                nColsStruct = 9
                nColsAero   = 7
                struct_headers = ['r_[m]','EIFlap_[Nm2]','EIEdge_[Nm2]','Mass_[kg/m]','Str.Twist_[deg]','CalcOutput_[0/1]','PreBendFlap_[m]','PreBendEdge_[m]','GKt_[Nm2]']
                aero_headers   = ['X_BladeRoot_[m]','Chord_[m]','AeroTwist_[deg]','RelThickness_[%]','AeroCenter_[m]','AeroTorsion_[0/1]','ProfileSet_[#]']
            elif self.version==3:
                nColsStruct = 13
                nColsAero   = 7
                struct_headers = ['r_[m]','EIFlap_[Nm2]','EIEdge_[Nm2]','GKt_[Nm2]','Mass_[kg/m]','Jxx_[kg.m]','PreBendFlap_[m]','PreBendEdge_[m]','Str.Twist_[deg]','PhiOut_[deg]','Ycog_[m]','Yshc_[m]','CalcOutput_[0/1]']
                aero_headers   = ['X_BladeRoot_[m]','Chord_[m]','AeroTwist_[deg]','RelThickness_[%]','AeroCenter_[m]','AeroTorsion_[0/1]','ProfileSet_[#]']
            else:
                raise BrokenFormatError('Blade format not implemented')

            struct        = np.zeros((nSections,nColsStruct))
            aero          = np.zeros((nSections,nColsAero))

            #  --- Structural data 
            try:
                for iSec in range(nSections):
                    vals=f.readline().split()
                    #if len(vals)>=nColsStruct:
                    struct[iSec,:]=np.array(vals[0:nColsStruct]).astype(float)
                    #elif self.version==1:
                    #    # version 1 has either 8 or 9 columns
                    #    nColsStruct=nColsStruct-1
                    #    struct_headers=struct_headers[0:-1]
                    #    struct =struct[:,:-1]
                    #    struct[iSec,:]=np.array(vals[0:nColsStruct]).astype(float)
            except:
                raise WrongFormatError('Unable to read structural data')
            try:
                self.BetaC  = float(f.readline().strip().split()[0])
                if self.version==3:
                    f.readline()
                    self.FlapDamping  = [float(v) for v in f.readline().strip().split(';')[0].split()]
                    self.EdgeDamping  = [float(v) for v in f.readline().strip().split(';')[0].split()]
                    self.TorsDamping  = [float(v) for v in f.readline().strip().split(';')[0].split()]
                    f.readline()
                    f.readline()
                else:
                    Damping  = [float(v) for v in f.readline().strip().split()[0:4]]
                    self.FlapDamping  = Damping[0:2]
                    self.EdgeDamping  = Damping[2:4]
                    self.TorsDamping  = []
            except:
                raise
                raise WrongFormatError('Unable to read damping data')

            # --- Aero
            try:
                for iSec in range(nSections):
                    vals=f.readline().split()[0:nColsAero]
                    aero[iSec,:]=np.array(vals).astype(float)
            except:
                raise WrongFormatError('Unable to read aerodynamic data')

            self.ProfileFile=f.readline().strip()

        # --- Concatenating aero and structural data
        self._cols = struct_headers+aero_headers[1:]
        data = np.column_stack((struct,aero[:,1:]))
        dataMiss=pd.DataFrame(data=data, columns=self._cols)
        self._nColsStruct=nColsStruct # to remember where to split
        # --- Making sure all columns are present, irrespectively of version
        self.data=pd.DataFrame(data=[], columns=headers_all)
        for c in self._cols:
            self.data[c]=dataMiss[c]

#     def toString(self):
#         s=''
#         if len(self.ProfileSets)>0:
#             prefix='PROFILE SET '
#         else:
#             prefix=''
#         for pset in self.ProfileSets:
#             s+=pset.toString(prefix)
#         return s
# 
#     def _write(self):
#         with open(self.filename,'w') as f:
#             f.write(self.toString)
# 
    def __repr__(self):
        s ='Class FLEXBladeFile (attributes: data, BetaC, FlapDamping, EdgeDamping, ProfileFile)\n'
        return s

    def _toDataFrame(self):
        return self.data

