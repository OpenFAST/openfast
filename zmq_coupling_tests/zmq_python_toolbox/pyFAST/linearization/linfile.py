#############################################################
#	Extract OpenFAST Matrices from linearization files  #
#	Authors: Srinivasa B. Ramisetti    	            #
#	Created:   01-July-2020		  	            #
#	E-mail: ramisettisrinivas@yahoo.com		    #
#	Web:	http://ramisetti.github.io		    #
#############################################################
#!/usr/bin/env python

import os, distutils
import collections
from itertools import islice
import numpy as np
import re
from pyFAST.input_output.fast_linearization_file import FASTLinearizationFile

def _isbool(str):
    flag=0
    if str=='T' or str=='True':
        flag=1
    return flag


def findBladeTriplets(rotFrame, Desc, verbose=True):

    # Find the number of, and indices for, triplets in the rotating frame:
    chkStr = [r'[Bb]lade \d', r'[Bb]lade [Rr]oot \d', r'BD_\d', r'[Bb]\d', r'[Bb]lade\d', r'PitchBearing\d', r'\d']

    origDesc = Desc;

    # hack for ElastoDyn state names (remove unnecessary text in parenthesis)
    for i in range(len(rotFrame)):
        if Desc[i] is None:
            raise Exception('Description not defined, likely a bug.')
        ix = Desc[i].find('(internal DOF index = ')
        if ix>0:
            ix2 = Desc[i].find(')')
            Desc[i] = Desc[i][:ix]+Desc[i][ix2+1:]


    NTriplets = 0;              # first initialize to zero
    Triplets = [];
    for i in range(len(rotFrame)):  # loop through inputs/outputs/states
        if rotFrame[i] == 'T':          # this is in the rotating frame
            Tmp = -1*np.ones(3)
            foundTriplet = False
            foundBladeNumber = False
            for chk in chkStr:
                BldNoCol = re.search(chk,Desc[i])
                if BldNoCol!=None:
                    foundBladeNumber = True

                    Bldstr=BldNoCol.group()
                    # create another regular expression to find the
                    # exact match on a different blade:
                    strng = re.split(Bldstr,Desc[i],1)  #this should return the strings before and after the match

                    
                    FirstStr = strng[0] + Bldstr[:len(Bldstr)-1] + '.'
                    checkThisStr = FirstStr + strng[1]
                    
                    #we need to get rid of the special characters that
                    #may exist in Desc{}:
                    checkThisStr=checkThisStr.replace(')',r'\)').replace('(', r'\(').replace('^',r'\^')
                    FirstStr   = FirstStr.    replace(')',r'\)').replace('(', r'\(').replace('^',r'\^')

                    k = int(Bldstr[len(Bldstr)-1])
                    Tmp[k-1] = int(i)
                    break

                #print(Tmp,j)
            
            # find the other match values
            if foundBladeNumber:
                for j in range((i+1),len(rotFrame)):           # loop through all remaining control inputs
                    #print(i, j)
                    if rotFrame[j]:                       # this is in the rotating frame
                        BldNoCol = re.search(checkThisStr, Desc[j])
                        if BldNoCol!=None:
                            Num = re.search(FirstStr,Desc[j]).group();
                            k = int(Num[len(Num)-1]);
                            Tmp[k-1] = int(j);                             # save the indices for the remaining blades
                            #TmpTmp=Tmp+1
                            #print(BldNoCol.group(),i,j,k)
                            if ( (Tmp>-1).all() ):                         # true if all the elements of Tmp are nonzero; thus, we found a triplet of rotating indices
                                foundTriplet = True;
                                
                                Triplets.append(Tmp);        # these are the indices for control input triplets in the rotating frame

                                NTriplets = NTriplets + 1;          # this  is  the number  of  control input triplets in the rotating frame

                                # we'll set rotFrame to false so that we don't have to check the found channels again; also allows us to throw error if we have a rotating channel that doesn't have a unique match
                                for idx in Tmp:
                                    id=int(idx)
                                    rotFrame[id] = 0;
                                
                                break;
                
                if foundTriplet==False:
                    if verbose:
                        print('Rotating channel "', i, Desc[i], '" does not form a unique blade triplet. Blade(s) not found: ', np.array(np.where(Tmp == -1))+1 )
            else:
                if verbose:
                    print( 'Could not find blade number in rotating channel "', Desc[i], '".')
    #print(NTriplets)
    #print(Triplets)

    return Triplets, NTriplets

def reOrderByIdx_1D(arry,Indi):
    tmp=[None]*len(arry)
    j=0
    for i in Indi:
        tmp[i]=arry[j]
        j=j+1
    return tmp

def reOrderByIdx_2D(arry,Indi,Indj):
    #tmp=[[None]*arry.shape[0]]*arry.shape[1]
    tmp=np.empty((arry.shape[0], arry.shape[1])) 
    #print(arry.shape, len(Indi), len(Indj))
    if arry.shape[0]!= len(Indi):
        Indi=np.arange(0,arry.shape[0])
    if arry.shape[1]!= len(Indj):
        Indj=np.arange(0,arry.shape[1])

    p=0
    q=0
    for i in Indi:
        q=0
        for j in Indj:
            tmp[i][j]=arry[p][q]
            q=q+1
        p=p+1
    return tmp

def getStateOrderingIndx(matData):

    #StateOrderingIndx={}
    StateOrderingIndx = np.arange(0,matData['NumStates'])
    lastModName = '';
    lastModOrd  = 0;
    mod_nDOFs   = 0;    # number of DOFs in each module
    sum_nDOFs2  = 0;    # running total of second-order DOFs
    sum_nDOFs1  = 0;    # running total of first-order DOFs
    indx_start  = 0;    # starting index of the modules
    
    for i in range(0,matData['NumStates']):

        tmp=(matData['DescStates'][i]); # name of the module whose states we are looking at
        modName = tmp.split(' ')[0]
        ModOrd  = matData['StateDerivOrder'][i]

        if lastModName!=modName or lastModOrd!= ModOrd:
            # this is the start of a new set of DOFs, so we'll set the
            # indices for the last matrix
            if lastModOrd == 2:
                mod_nDOFs = int(mod_nDOFs/2);
                #print('mmmm ', mod_nDOFs, indx_start, matData['ndof2'], i)
                StateOrderingIndx[  indx_start           :(indx_start+mod_nDOFs)] = sum_nDOFs2 +                 np.arange(0,mod_nDOFs); # q2 starts at 1
                StateOrderingIndx[ (indx_start+mod_nDOFs):(i)]                  = sum_nDOFs2 + matData['ndof2'] + np.arange(0,mod_nDOFs); # q2_dot starts at matData.ndof2 + 1

                sum_nDOFs2 = sum_nDOFs2 + mod_nDOFs;
            else:
                if indx_start < mod_nDOFs:
                    StateOrderingIndx[indx_start:(indx_start+mod_nDOFs)] = sum_nDOFs1 + matData['NumStates2'] + np.arange(0,mod_nDOFs); # q1 starts at matData.NumStates2 + 1
                
                sum_nDOFs1 = sum_nDOFs1 + mod_nDOFs;
            
            # reset for a new module (or new 1st-order states in the same module)
            mod_nDOFs = 0;
            
            indx_start = i; #start of this module
            lastModName = modName;
            lastModOrd  = matData['StateDerivOrder'][i];
            
        mod_nDOFs = mod_nDOFs+1;
        

    # repeat for the last module found:
    if lastModOrd == 2:
        mod_nDOFs = int(mod_nDOFs/2);
        #print(mod_nDOFs,indx_start)
        StateOrderingIndx[indx_start:(indx_start+mod_nDOFs)] = sum_nDOFs2 + np.arange(0,mod_nDOFs); # q2 starts at 1
        StateOrderingIndx[(indx_start+mod_nDOFs):matData['NumStates']]        = sum_nDOFs2 + matData['ndof2'] + np.arange(0,mod_nDOFs); # q2_dot starts at matData.ndof2 + 1
    else:
        StateOrderingIndx[indx_start:(indx_start+mod_nDOFs)] = sum_nDOFs1 + matData['NumStates2'] + np.arange(0,mod_nDOFs); # q1 starts at matData.NumStates2 + 1

    #print(StateOrderingIndx)
    return StateOrderingIndx

def readOP(fid, n, defaultDerivOrder=2):
    OP=[]
    Var = {'RotatingFrame': [], 'DerivativeOrder': [], 'Description': []}
    colNames=fid.readline().strip()
    dummy=   fid.readline().strip()
    bHasDeriv= colNames.find('Derivative Order')>=0
    for i, line in enumerate(fid):
        sp=line.strip().split()
        if sp[1].find(',')>=0:
            #  Most likely this OP has three values (e.g. orientation angles)
            # For now we discard the two other values
            OP.append(float(sp[1][:-1]))
            iRot=4
        else:
            OP.append(float(sp[1]))
            iRot=2
        Var['RotatingFrame'].append(sp[iRot])
        if bHasDeriv:
            Var['DerivativeOrder'].append(int(sp[iRot+1]))
            Var['Description'].append(' '.join(sp[iRot+2:]).strip())
        else:
            Var['DerivativeOrder'].append(defaultDerivOrder)
            Var['Description'].append(' '.join(sp[iRot+1:]).strip())
        if i>=n-1:
            break

    tmp = dict()
    tmp['x_op']         = OP
    tmp['x_rotFrame']   = Var['RotatingFrame']
    tmp['x_DerivOrder'] = Var['DerivativeOrder']
    tmp['x_desc']       = Var['Description']

    return tmp



def readFASTMatrix(f):
    name=""
    m=0
    tmp=[]
    for line in islice(f,1):
        # copy matrix name to tmp_name
        tmp_name=line.strip().split(':')[0]
        # get matrix dimensions if matrix exists
        if tmp_name != "":
            m=int(line.strip().split(':')[1].split('x')[0])
            n=int(line.strip().split(':')[1].split('x')[1])

    # copy matrix into tmp list
    if m!=0:
        name=tmp_name
        for line in islice(f,m):
            tmp.append([float(num) for num in line.split()])
        tmp=np.array(tmp)

    return name,tmp
    
def ReadFASTLinear(filename, starSub=None, removeStatesPattern=None):
    """ 
    Read one lin file.

    INPUTS:
     - filename: linfile name, string.
     - starSub: if None, raise an error if `****` are present
                otherwise replace *** with `starSub` (e.g. 0)
                see FASTLinearizationFile. 
     - removeStatesPattern: remove states matching a giving description pattern.
               e.g:  'tower|Drivetrain'  or '^AD'
               see FASTLinearizationFile. 
    OUTPUTS:
     - data: a dictionary with fields matching the MATLAB implementation. 
    """

    # --- Read lin file
    f = FASTLinearizationFile(filename, starSub=starSub, removeStatesPattern=removeStatesPattern)

    # --- Legacy structure. TODO, just use "f"
    data={}
    data['n_x'] = f.nx
    data['n_u'] = f.nu
    data['n_y'] = f.ny
    data['n_z'] = f.nz
    KeyMap={'t':'t'}
    KeyMap.update({'Azimuth':'Azimuth', 'WindSpeed':'WindSpeed', 'RotSpeed':'RotSpeed'})
    KeyMap.update({'x_op':'x', 'xdot_op':'xdot', 'z_op':'z', 'y_op':'y'})
    for knew,kold in KeyMap.items():
        if kold in f.keys():
            data[knew] = f[kold]
    data['WindSpeed'] = data['WindSpeed'] if data['WindSpeed'] is not None else np.nan
    data['RotSpeed']  = data['RotSpeed'] if data['RotSpeed'] is not None else np.nan
    data['Azimuth']   = np.mod(data['Azimuth'],2.0*np.pi)
    if data['n_x']>0:
        data['x_rotFrame']   = f['x_info']['RotatingFrame']
        data['x_DerivOrder'] = f['x_info']['DerivativeOrder']
        data['x_desc']       = f['x_info']['Description']
        data['xdot_desc']    = f['xdot_info']['Description']
        data['n_x2'] = np.sum(data['x_DerivOrder']== 2)
        data['x_op'] = f['x']
    else:
        data['n_x2'] = 0;
    if data['n_z']>0:
        data['z_desc']     = f['z_info']['Description']
        data['z_op']       = f['z']
    if data['n_u']>0:
        data['u_rotFrame'] = f['u_info']['RotatingFrame']
        data['u_desc']     = f['u_info']['Description']
        data['u_op']       = f['u']
    if data['n_y']>0:
        data['y_rotFrame'] = f['y_info']['RotatingFrame']
        data['y_desc']     = f['y_info']['Description']
        data['y_op']       = f['y']
    for k in ['A','B','C','D']:
        if k in f.keys():
            data[k] = f[k]
    return data, None

def ReadFASTLinearLegacy(filename):
    """ 
    Legacy reader, TODO remove in future release
    """

    def extractVal(lines, key):
        for l in lines:
            if l.find(key)>=0:
                return l.split(key)[1].split()[0]
        return None

    def readToMarker(fid, marker, nMax):
        lines=[]
        for i, line in enumerate(fid):
            if i>nMax:
                raise BrokenFormatError('`{}` not found in file'.format(marker))
            if line.find(marker)>=0:
                break
            lines.append(line.strip())
        return lines, line

    with open(filename) as f:
        info = {}
        data = {}
        SetOfMatrices = 1
        info['name'] = os.path.splitext(os.path.basename(filename))[0]

        # --- 
        header, lastLine=readToMarker(f, 'Jacobians included', 30)
        header.append(lastLine)
        data['t']   = float(extractVal(header,'Simulation time:'          ))
        data['n_x'] = int(extractVal(header,'Number of continuous states:'))
        data['n_xd'] = int(extractVal(header,'Number of discrete states:'  ))
        data['n_z'] = int(extractVal(header,'Number of constraint states:'))
        data['n_u'] = int(extractVal(header,'Number of inputs:'           ))
        data['n_y'] = int(extractVal(header,'Number of outputs:'          ))
        bJac =    extractVal(header,'Jacobians included in this file?')
        if bJac: 
            SetOfMatrices = 2
        try:
            data['Azimuth'] = float(extractVal(header,'Azimuth:'))
        except:
            data['Azimuth'] = None
        try:
            data['RotSpeed'] = float(extractVal(header,'Rotor Speed:')) # rad/s
        except:
            data['RotSpeed'] = np.nan
        try:
            data['WindSpeed'] = float(extractVal(header,'Wind Speed:'))
        except:
            data['WindSpeed'] = np.nan

        data['Azimuth']=np.mod(data['Azimuth'],2.0*np.pi)
        try:
            # skip next three lines
            for line in islice(f,2):
                pass
            if data['n_x'] > 0:
                temp = readOP(f, data['n_x'])
                data['x_op']=temp['x_op']
                data['x_rotFrame']=temp['x_rotFrame']
                data['x_DerivOrder']=temp['x_DerivOrder']
                data['x_desc']=temp['x_desc']

                # skip next three lines
                for line in islice(f,2):
                    pass

                temp = readOP(f, data['n_x'], defaultDerivOrder=2)
                data['xdot_op']=temp['x_op']
                data['xdot_desc']=temp['x_desc']

                #(number of second-order states)
                data['n_x2'] = sum(1 for i in data['x_DerivOrder'] if i == 2)
            else:
                data['n_x2'] = 0;
            
            
            if data['n_xd'] > 0:
                # skip next three lines
                for line in islice(f,2):
                    pass
                temp = readOP(f, data['n_xd'], defaultDerivOrder=2)
                data['xd_op']=temp['x_op']
                data['xd_desc']=temp['x_desc']
            if data['n_z'] > 0:
                # skip next three lines
                for line in islice(f,2):
                    pass
                temp = readOP(f, data['n_z'], defaultDerivOrder=0)
                data['z_op']=temp['x_op']
                data['z_desc']=temp['x_desc']
            if data['n_u'] > 0:
                # skip next three lines
                for line in islice(f,2):
                    pass
                temp = readOP(f, data['n_u'], defaultDerivOrder=0)
                data['u_op']=temp['x_op']
                data['u_desc']=temp['x_desc']
                data['u_rotFrame']=temp['x_rotFrame']
            if data['n_y'] > 0:
                # skip next three lines
                for line in islice(f,2):
                    pass
                temp = readOP(f, data['n_y'], defaultDerivOrder=0)
                data['y_op']=temp['x_op']
                data['y_desc']=temp['x_desc']
                data['y_rotFrame']=temp['x_rotFrame']

            # skip next one line
            for line in islice(f,4):
                pass

            mat_names=[]
            while True:
                name,mat=readFASTMatrix(f)
                if not name:
                    break;
                mat_names.append(name)
                data[name]=mat

            return data, info

        except (ValueError, AssertionError):
            raise


def get_Mats(FileNames, verbose=True, removeTwrAzimuth=False, starSub=None, removeStatesPattern=None):
    """ 
    Extra main data from a list of lin files.
    INPUTS:
     - FileNames : list of lin files
     - starSub: if None, raise an error if `****` are present
                otherwise replace *** with `starSub` (e.g. 0)
                see FASTLinearizationFile. 
     - removeStatesPattern: remove states matching a giving description pattern.
                e.g:  'tower|Drivetrain'  or '^AD'
                see FASTLinearizationFile. 
     - removeTwrAzimuth: if False do nothing
                otherwise discard lin files where azimuth in [60, 180, 300]+/-4deg (close to tower). 

    """
    NAzimStep = len(FileNames)
    data     = [None]*NAzimStep
    # --- Read all files
    for iFile, filename in enumerate(FileNames):
        #data[iFile], _ = ReadFASTLinearLegacy(FileNames[iFile]);
        data[iFile],_= ReadFASTLinear(FileNames[iFile], starSub=starSub, removeStatesPattern=removeStatesPattern);
    Azimuth  = np.array([d['Azimuth'] for d in data])*180/np.pi
    # --- Sort by azimuth, not required, but filenames are not sorted, nor are the lin times azimuth
    ISort = np.argsort(Azimuth)
    Azimuth   = Azimuth[ISort]
    data      = [data[i]      for i in ISort]
    FileNames = [FileNames[i] for i in ISort]
    # --- Remove some azimuth
    if removeTwrAzimuth:
        IDiscard = [i  for i,a in enumerate(Azimuth) if np.any(np.abs(np.array([60,180,300])-a)<4)  ]
        n=len(FileNames[0]); sfmt='{:'+str(n+2)+'s}'
        for i in IDiscard:
            print('          discard: '+sfmt.format(FileNames[i]) + '  (psi={:5.1f})'.format(Azimuth[i]))
        data      = [d for i,d in enumerate(data)      if i not in IDiscard]
        NAzimStep= len(data)


    matData={}
    matData['NAzimStep']= NAzimStep;

    # --- Using last azimuth to initial some of the "matData" variables
    # Input data from linearization files
    matData['NumStates']       = int(data[NAzimStep-1]['n_x']);
    matData['NumStates2']      = int(data[NAzimStep-1]['n_x2']);

    #return data, matData
    matData['ndof1']           = int(matData['NumStates'] - matData['NumStates2']); # number of first-order states = number of first-order DOFs
    matData['ndof2']           = int(data[NAzimStep-1]['n_x2'] / 2); # half the number of second-order states = number of second-order DOFs
    #% matData['ndof']            = matData['ndof2'] + matData['ndof1']; #half the number of second-order states plus the number of first-order states (i.e., states that aren't derivatives)

    matData['NumInputs']       = int(data[NAzimStep-1]['n_u']);
    matData['NumOutputs']      = int(data[NAzimStep-1]['n_y']);

    # allocate space for these variables
    matData['Azimuth']   = np.zeros(NAzimStep);
    matData['Omega']     = np.zeros(NAzimStep);
    matData['OmegaDot']  = np.zeros(NAzimStep);
    matData['WindSpeed'] = np.zeros(NAzimStep)*np.nan;

    if matData['NumStates'] > 0:
        matData['DescStates']      = data[NAzimStep-1]['x_desc'];
        matData['StateDerivOrder'] = data[NAzimStep-1]['x_DerivOrder'];
        matData['xdop']            = np.zeros((matData['NumStates'], NAzimStep))
        matData['xop']             = np.zeros((matData['NumStates'], NAzimStep))
        matData['A']               = np.zeros((matData['NumStates'], matData['NumStates'], NAzimStep))

    if matData['NumInputs'] > 0:
        matData['DescCntrlInpt'] = data[NAzimStep-1]['u_desc'];    
        matData['u_op']          = np.zeros((matData['NumInputs'],NAzimStep))
        if matData['NumStates']>0:
            matData['B'] = np.zeros((matData['NumStates'], matData['NumInputs'],NAzimStep))

    if matData['NumOutputs'] > 0:
        matData['DescOutput']    = data[NAzimStep-1]['y_desc']
        matData['y_op']          = np.zeros((matData['NumOutputs'],NAzimStep))

        if matData['NumStates'] > 0:
            matData['C']         = np.zeros((matData['NumOutputs'], matData['NumStates'], NAzimStep))
        if matData['NumInputs'] > 0:
            matData['D']         = np.zeros((matData['NumOutputs'], matData['NumInputs'], NAzimStep))

    # Reorder state matrices so that they follow the {q2, q2_dot, q1}
    # format that is assumed in the MBC3 equations.
    if matData['NumStates'] > 0:
        # keep StateOrderingIndx for applying inverse of MBC3 later
        # (to visualize mode shapes)
        matData['StateOrderingIndx'] = getStateOrderingIndx(matData)

        sortedIndx=matData['StateOrderingIndx']
        #print(sortedIndx)
        x_rotFrame= reOrderByIdx_1D(data[NAzimStep-1]['x_rotFrame'],sortedIndx)
        matData['DescStates'] = reOrderByIdx_1D(data[NAzimStep-1]['x_desc'],sortedIndx)
        matData['StateDerivOrder'] = reOrderByIdx_1D(data[NAzimStep-1]['x_DerivOrder'],sortedIndx)

    # --- Store file data into matData
    for iFile in np.arange(0,NAzimStep):
        matData['Omega'][iFile]     = data[iFile]['RotSpeed']         ; 
        matData['Azimuth'][iFile]   = data[iFile]['Azimuth']*180/np.pi; 
        matData['WindSpeed'][iFile] = data[iFile]['WindSpeed']
    
        if 'A' in data[iFile]:
            matData['A'][:,:,iFile]=reOrderByIdx_2D(data[iFile]['A'],sortedIndx,sortedIndx)
            #print('size of matData[A] for file ', iFile, ' is :', matData['A'].shape)
        if 'B' in data[iFile]:
            matData['B'][:,:,iFile]=reOrderByIdx_2D(data[iFile]['B'],sortedIndx,np.arange(0,matData['NumStates']))
            #print('size of matData[B] for file ', iFile, ' is :', matData['B'].shape)
        if 'C' in data[iFile]:
            matData['C'][:,:,iFile]=reOrderByIdx_2D(data[iFile]['C'],np.arange(0,matData['NumStates']),sortedIndx)
            #print('size of matData[C] for file ', iFile, ' is :', matData['C'].shape)
        if 'D' in data[iFile]:
            matData['D'][:,:,iFile]=reOrderByIdx_2D(data[iFile]['D'],sortedIndx,sortedIndx)
            #print('size of matData[D] for file ', iFile, ' is :', matData['D'].shape)

        if 'x_op' in data[iFile]:
            matData['xop'][:,iFile] = reOrderByIdx_1D(data[iFile]['x_op'],sortedIndx)
            #print('matData[xop] for file ', iFile,' is :',(matData['xop'][:,iFile]))
        if 'xdot_op' in data[iFile]:
            matData['xdop'][:,iFile] = reOrderByIdx_1D(data[iFile]['xdot_op'],sortedIndx)
            #print('matData[xdop] for file ', iFile,' is :',(matData['xdop'][:,iFile]))

        if 'u_op' in data[iFile]:
            matData['u_op'][:,iFile] = data[iFile]['u_op']

        if 'y_op' in data[iFile]:
            matData['y_op'][:,iFile] = data[iFile]['y_op']


    # Find the azimuth-averaged linearized 1st order state matrices:
    if 'A' in matData:
        matData['Avgxdop'] = np.mean(matData['xdop'],axis=1)
        matData['Avgxop']  = np.mean(matData['xop'], axis=1)
        matData['AvgA']    = np.mean(matData['A'],axis=2)

    #print(matData['AvgA'])
    if 'B' in matData:
        matData['AvgB']    = np.mean(matData['B'],axis=2)

    if 'C' in matData:
        matData['AvgC']    = np.mean(matData['C'],axis=2)

    if 'D' in matData:
        matData['AvgD']    = np.mean(matData['D'],axis=2)


    foundED  = True;
    for i in range(matData['ndof2']):
        # find the starting index of the string 'DOF_GeAz'
        if (matData['DescStates'][i].find('DOF_GeAz') != -1):
            matData['Omega']    = matData['xdop'][i,:]
            matData['OmegaDot'] = matData['xdop'][i+matData['ndof2'],:]
            foundED = True
            break

    for i in range(matData['ndof2']):
        # find the starting index of the string 'DOF_DrTr'
        if (matData['DescStates'][i].find('DOF_DrTr') != -1):
            matData['Omega']    = matData['Omega']    + matData['xdop'][i,:] #This always comes after DOF_GeAz so let's just add it here (it won't get written over later).
            matData['OmegaDot'] = matData['OmegaDot'] + matData['xdop'][i+matData['ndof2'],:]
            foundED = True
            break

    if not foundED:
        for i in range(matData['ndof2']):
            # find the starting index of the string 'Gearbox_Rot'
            if (matData['DescStates'][i].find('MBD Gearbox_Rot') != -1):
                matData['Omega']    = matData['xdop'][i,:]
                matData['OmegaDot'] = matData['xdop'][i+matData['ndof2'],:]
                break

    #print("\n".join(matData['DescStates']))
    #exit()
    # ----------- Find multi-blade coordinate (MBC) transformation indices ----

    # Find the indices for, state triplets in the rotating frame
    #  (note that we avoid the "first time derivative" states)
    if matData['ndof2'] > 0:
        matData['RotTripletIndicesStates2'], matData['n_RotTripletStates2'] = findBladeTriplets(x_rotFrame[0:matData['ndof2']],matData['DescStates'][0:matData['ndof2']], verbose=verbose)
    else:
        matData['RotTripletIndicesStates2'] = [];
        matData['n_RotTripletStates2'] = 0;

    if matData['ndof1'] > 0:
        matData['RotTripletIndicesStates1'], matData['n_RotTripletStates1'] = findBladeTriplets( x_rotFrame[matData['NumStates2']:] ,matData['DescStates'][matData['NumStates2']:] , verbose=verbose);
    else:
        matData['RotTripletIndicesStates1'] = [];
        matData['n_RotTripletStates1'] = 0;

    # Find the indices for control input triplets in the rotating frame:
    if matData['NumInputs'] > 0:
        matData['RotTripletIndicesCntrlInpt'], matData['n_RotTripletInputs'] = findBladeTriplets(data[0]['u_rotFrame'],matData['DescCntrlInpt'], verbose=verbose );
    else:
        matData['RotTripletIndicesCntrlInpt'] = [];
        matData['n_RotTripletInputs'] = 0;

    # Find the indices for output measurement triplets in the rotating frame:
    if (matData['NumOutputs'] > 0 ):
        matData['RotTripletIndicesOutput'], matData['n_RotTripletOutputs'] = findBladeTriplets(data[0]['y_rotFrame'],matData['DescOutput'], verbose=verbose );
    else:
        matData['RotTripletIndicesOutput'] = [];
        matData['n_RotTripletOutputs'] = 0;
        
    return matData, data

if __name__ == '__main__':
    madData,dat = get_Mats(['file.lin'])
