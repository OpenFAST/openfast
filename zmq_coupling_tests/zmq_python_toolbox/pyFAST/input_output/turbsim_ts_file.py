from itertools import takewhile
import pandas as pd
import numpy as np
try:
    from .file import File, WrongFormatError, BrokenFormatError
except:
    File = dict
    class WrongFormatError(Exception): pass
    class BrokenFormatError(Exception): pass


class TurbSimTSFile(File):

    @staticmethod
    def defaultExtensions():
        return ['.txt']

    @staticmethod
    def formatName():
        return 'TurbSim time series'

    def _read(self, *args, **kwargs):
        self['header']=[]
        nHeaderMax=10
        # Reading 
        iFirstData=-1
        with open(self.filename, 'r', errors="surrogateescape") as f:
            for i, line in enumerate(f):
                if i>nHeaderMax:
                    raise BrokenFormatError('`nComp` not found in file')
                if line.lower().find('ncomp')>=0:
                    iFirstData=i
                    break
                self['header'].append(line.strip())
            self['nComp'] = int(line.split()[0])
            line          = f.readline().strip()
            nPoints       = int(line.split()[0])
            line          = f.readline().strip()
            self['ID']    = int(line.split()[0])
            f.readline()
            f.readline()
            self['Points']=np.zeros((nPoints,2))
            for i in np.arange(nPoints):
                line = f.readline().strip()
                self['Points'][i,:]= np.array(line.split()).astype(float)
            f.readline()
            f.readline()
            f.readline()
            lines=[]
            # reading full data
            self['data'] = np.array([l.strip().split() for l in takewhile(lambda x: len(x.strip())>0, f.readlines())]).astype(float)

    def columns(self):
        Comp=['u','v','w']
        return ['Time']+['Point{}{}'.format(ip+1,Comp[ic]) for ic in np.arange(self['nComp']) for ip in np.arange(len(self['Points']))]

    def units(self):
        nPoints = self['Points'].shape[0]
        return ['(s)'] +  ['(m/s)']*nPoints*self['nComp']

    def toString(self):

        def toStringVLD(val,lab,descr):
            val='{}'.format(val)
            lab='{}'.format(lab)
            if len(val)<13:
                val='{:13s}'.format(val)
            if len(lab)<13:
                lab='{:13s}'.format(lab)
            return val+' '+lab+' - '+descr.strip().strip('-')+'\n'

        s='\n'.join(self['header'])+'\n'
        nPoints = self['Points'].shape[0]
        s+=toStringVLD(self['nComp'],'nComp'  ,'Number of velocity components in the file'         )
        s+=toStringVLD(nPoints      ,'nPoints','Number of time series points contained in this file(-)')
        s+=toStringVLD(self['ID']   ,'RefPtID','Index of the reference point (1-nPoints)')
        s+='{:^16s}{:^16s} {}\n'.format('Pointyi','Pointzi','! nPoints listed in order of increasing height')
        s+='{:^16s}{:^16s}\n'.format('(m)','(m)')
        for row in self['Points']:
            s+=''.join(['{:16.8e}'.format(v) for v in row])+'\n'

        s+='--------Time Series-------------------------------------------------------------\n'
        s+=''.join(['{:^16s}'.format(c) for c in self.columns()])+'\n'
        s+=''.join(['{:^16s}'.format(c) for c in self.units()])+'\n'
        s+='\n'.join(''.join('{:16.8e}'.format(x) for x in y) for y in self['data'])
        return s

    def _write(self):
        with open(self.filename,'w') as f:
            f.write(self.toString())

        

    def _toDataFrame(self):
        Cols = ['{}_{}'.format(c.replace(' ','_'), u.replace('(','[').replace(')',']')) for c,u in zip(self.columns(),self.units())]
        dfs={}
        dfs['Points']     = pd.DataFrame(data = self['Points'],columns = ['PointYi','PointZi'])
        dfs['TimeSeries'] = pd.DataFrame(data = self['data'] ,columns = Cols)

        return dfs


