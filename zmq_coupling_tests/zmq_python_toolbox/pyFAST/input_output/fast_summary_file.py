import numpy as np
import pandas as pd
import os
# Local 
from .mini_yaml import yaml_read

try:
    from .file import File, EmptyFileError
except:
    EmptyFileError = type('EmptyFileError', (Exception,),{})
    File=dict

# --------------------------------------------------------------------------------}
# --- Main Class 
# --------------------------------------------------------------------------------{
class FASTSummaryFile(File):
    """ 
    Read an OpenFAST summary file (.sum, .yaml). The object behaves as a dictionary.
    NOTE: open new subdyn format supported.

    Main methods
    ------------
    - read, toDataFrame

    Examples
    --------

        # read a subdyn summary file
        sum = FASTSummaryFile('5MW.SD.sum.yaml')
        print(sum['module']) # SubDyn
        M = sum['M'] # Mass matrix
        K = sum['K'] # stiffness matrix

    """

    @staticmethod
    def defaultExtensions():
        return ['.sum','.yaml']

    @staticmethod
    def formatName():
        return 'FAST summary file'

    def __init__(self,filename=None, **kwargs):
        self.filename = None
        if filename:
            self.read(filename, **kwargs)

    def read(self, filename=None, header_only=False):
        """ """
        if filename:
            self.filename = filename
        if not self.filename:
            raise Exception('No filename provided')
        if not os.path.isfile(self.filename):
            raise OSError(2,'File not found:',self.filename)
        if os.stat(self.filename).st_size == 0:
            raise EmptyFileError('File is empty:',self.filename)

        with open(self.filename, 'r', errors="surrogateescape") as fid:
            header= readFirstLines(fid, 4)
        if any(['subdyn' in s.lower() for s in header]):
            self['module']='SubDyn'
            readSubDynSum(self)
        else:
            raise NotImplementedError('This summary file format is not yet supported')

    def toDataFrame(self):
        if 'module' not in self.keys():
            raise Exception('');

        if self['module']=='SubDyn':
            raise Exception('This should not happen since class was added to subdyn object')
        #    dfs=subDynToDataFrame(self)

        return dfs

    def toGraph(self):
        from .fast_input_file_graph import fastToGraph
        return fastToGraph(self)



# --------------------------------------------------------------------------------}
# --- Helper functions 
# --------------------------------------------------------------------------------{
def readFirstLines(fid, nLines):
    lines=[]
    for i, line in enumerate(fid):
        lines.append(line.strip())
        if i==nLines:
            break
    return lines

# --------------------------------------------------------------------------------}
# --- Sub-reader/class for SubDyn summary files
# --------------------------------------------------------------------------------{
def readSubDynSum(self):

    # Read data
    #T=yaml.load(fid, Loader=yaml.SafeLoader)
    yaml_read(self.filename, self)

    # --- Treatement of useful data
    if self['DOF2Nodes'].shape[1]==3:
        self['DOF2Nodes']=np.column_stack((np.arange(self['DOF2Nodes'].shape[0])+1,self['DOF2Nodes']))
    # NOTE: DOFs are reindexed to start at 0
    self['DOF2Nodes'][:,0]-=1
    self['DOF___L'] -=1 # internal DOFs
    self['DOF___B'] -=1 # internal
    self['DOF___F'] -=1 # fixed DOFs

    self['CB_frequencies']=self['CB_frequencies'].ravel()
    self['X'] = self['Nodes'][:,1].astype(float)
    self['Y'] = self['Nodes'][:,2].astype(float)
    self['Z'] = self['Nodes'][:,3].astype(float)

    # --- Useful methods that will be added to the class
    def NodesDisp(self, IDOF, UDOF, maxDisp=None, sortDim=None):
        DOF2Nodes = self['DOF2Nodes']
        # NOTE: SubDyn nodes in the summary files are sorted
        # so the position we give are for all Nodes
        INodes = list(np.sort(np.unique(DOF2Nodes[IDOF,1]))) # Sort
        nShapes = UDOF.shape[1]
        disp=np.empty((len(INodes),3,nShapes)); disp.fill(np.nan)
        pos=np.empty((len(INodes),3))         ; pos.fill(np.nan)
        # TODO
        #   handle T_red for rigid and joints
        for i,iDOF in enumerate(IDOF):
            iNode       = DOF2Nodes[iDOF,1]
            nDOFPerNode = DOF2Nodes[iDOF,2]
            nodeDOF     = DOF2Nodes[iDOF,3]
            iiNode      = INodes.index(iNode)
            if nodeDOF<=3:
                pos[iiNode, 0]=self['X'][iNode-1]
                pos[iiNode, 1]=self['Y'][iNode-1]
                pos[iiNode, 2]=self['Z'][iNode-1]
                for iShape in np.arange(nShapes):
                    disp[iiNode, nodeDOF-1, iShape] = UDOF[i, iShape]
        # Scaling 
        if maxDisp is not None:
            for iShape in np.arange(nShapes):
                mD=np.nanmax(np.abs(disp[:, :, iShape]))
                if mD>1e-5:
                    disp[:, :, iShape] *= maxDisp/mD
        # Sorting according to a dimension
        if sortDim is not None: 
            I=np.argsort(pos[:,sortDim])
            INodes = np.array(INodes)[I]
            disp   = disp[I,:,:]
            pos    = pos[I,:]
        return disp, pos, INodes

    def getModes(data, maxDisp=None, sortDim=None):
        """ return Guyan and CB modes"""
        if maxDisp is None:
            #compute max disp such as it's 10% of maxdimension
            dx = np.max(self['X'])-np.min(self['X'])
            dy = np.max(self['Y'])-np.min(self['Y'])
            dz = np.max(self['Z'])-np.min(self['Z'])
            maxDisp = np.max([dx,dy,dz])*0.1

        # NOTE: DOF have been reindexed -1
        DOF_B = data['DOF___B'].ravel()
        DOF_F = data['DOF___F'].ravel()
        DOF_K = (np.concatenate((DOF_B,data['DOF___L'].ravel(), DOF_F))).astype(int)

        # CB modes
        PhiM      = data['PhiM']
        Phi_CB = np.vstack((np.zeros((len(DOF_B),PhiM.shape[1])),PhiM, np.zeros((len(DOF_F),PhiM.shape[1]))))
        dispCB, posCB, INodesCB = data.NodesDisp(DOF_K, Phi_CB, maxDisp=maxDisp, sortDim=sortDim)
        # Guyan modes
        PhiR      = data['PhiR']
        Phi_Guyan = np.vstack((np.eye(len(DOF_B)),PhiR, np.zeros((len(DOF_F),PhiR.shape[1]))))
        dispGy, posGy, INodesGy = data.NodesDisp(DOF_K, Phi_Guyan, maxDisp=maxDisp, sortDim=sortDim)

        return dispGy, posGy, INodesGy, dispCB, posCB, INodesCB


    def subDynToJson(data, outfile=None):
        """ Convert to a "JSON" format

        TODO: convert to graph and use graph.toJSON

        """
        #return data.toGraph().toJSON(outfile)

        dispGy, posGy, _, dispCB, posCB, _ = data.getModes(sortDim=None) # Sorting mess things up

        Nodes    = self['Nodes'].copy()
        Elements = self['Elements'].copy()
        Elements[:,0]-=1
        Elements[:,1]-=1
        Elements[:,2]-=1
        CB_freq   = data['CB_frequencies'].ravel()

        d=dict();
        d['Connectivity']=Elements[:,[1,2]].astype(int).tolist();
        d['Nodes']=Nodes[:,[1,2,3]].tolist()
        d['ElemProps']=[{'shape':'cylinder','type':int(Elements[iElem,5]),'Diam':np.sqrt(Elements[iElem,7]/np.pi)*4} for iElem in range(len(Elements))] # NOTE: diameter is cranked up
    #  disp[iiNode, nodeDOF-1, iShape] = UDOF[i, iShape]

        d['Modes']=[
                {
                    'name':'GY{:d}'.format(iMode+1),
                    'omega':1,
                    'Displ':dispGy[:,:,iMode].tolist()
                }  for iMode in range(dispGy.shape[2]) ]
        d['Modes']+=[
                {
                    'name':'CB{:d}'.format(iMode+1),
                    'omega':CB_freq[iMode]*2*np.pi, #in [rad/s]
                    'Displ':dispCB[:,:,iMode].tolist()
                }  for iMode in range(dispCB.shape[2]) ]
        d['groundLevel']=np.min(data['Z']) # TODO

        if outfile is not None:
            import json
            with open(outfile, 'w', encoding='utf-8') as f:
                try:
                    f.write(unicode(json.dumps(d, ensure_ascii=False))) #, indent=2)
                except:
                    json.dump(d, f, indent=2) 
        return d


    def subDynToDataFrame(data, sortDim=2, removeZero=True):
        """ Convert to DataFrame containing nodal displacements """
        def toDF(pos,disp,preffix=''):
            disp[np.isnan(disp)]=0
            disptot=disp.copy()
            columns=[]
            for ishape in np.arange(disp.shape[2]):
                disptot[:,:,ishape]= pos + disp[:,:,ishape]
                sMode=preffix+'Mode{:d}'.format(ishape+1)
                columns+=[sMode+'x_[m]',sMode+'y_[m]',sMode+'z_[m]']
            disptot= np.moveaxis(disptot,2,1).reshape(disptot.shape[0],disptot.shape[1]*disptot.shape[2])
            disp   = np.moveaxis(disp,2,1).reshape(disp.shape[0],disp.shape[1]*disp.shape[2])
            df= pd.DataFrame(data = disptot ,columns = columns)
            dfDisp= pd.DataFrame(data = disp ,columns = columns)
            # remove mode components that are fully zero 
            if removeZero:
                df     = df.loc[:,     (dfDisp != 0).any(axis=0)]
                dfDisp = dfDisp.loc[:, (dfDisp != 0).any(axis=0)]
            dfDisp.columns = [c.replace('Mode','Disp') for c in dfDisp.columns.values]
            return df, dfDisp

        dispGy, posGy, _, dispCB, posCB, _ = data.getModes(sortDim=sortDim)

        columns = ['z_[m]','x_[m]','y_[m]']
        dataZXY = np.column_stack((posGy[:,2],posGy[:,0],posGy[:,1]))
        dfZXY   = pd.DataFrame(data = dataZXY, columns=columns)
        df1, df1d = toDF(posGy, dispGy,'Guyan')
        df2, df2d = toDF(posCB, dispCB,'CB')
        df = pd.concat((dfZXY, df1, df2, df1d, df2d), axis=1)
        return df 
    
    # adding method to class dynamically to give it a "SubDyn Summary flavor"
    setattr(FASTSummaryFile, 'NodesDisp'  , NodesDisp) 
    setattr(FASTSummaryFile, 'toDataFrame', subDynToDataFrame)
    setattr(FASTSummaryFile, 'toJSON'     , subDynToJson)
    setattr(FASTSummaryFile, 'getModes'   , getModes)

    return self


if __name__=='__main__':
    T=FASTSummaryFile('../Pendulum.SD.sum.yaml')
    df=T.toDataFrame()
    print(df)
