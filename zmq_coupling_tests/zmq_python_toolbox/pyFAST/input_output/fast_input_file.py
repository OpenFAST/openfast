import numpy as np
import os
import pandas as pd
import re
try:
    from .file import File, WrongFormatError, BrokenFormatError
except:
    File = dict
    class WrongFormatError(Exception): pass
    class BrokenFormatError(Exception): pass

__all__  = ['FASTInputFile']

TABTYPE_NOT_A_TAB          = 0
TABTYPE_NUM_WITH_HEADER    = 1
TABTYPE_NUM_WITH_HEADERCOM = 2
TABTYPE_NUM_NO_HEADER      = 4
TABTYPE_NUM_BEAMDYN        = 5
TABTYPE_NUM_SUBDYNOUT      = 7
TABTYPE_MIX_WITH_HEADER    = 6
TABTYPE_FIL                = 3
TABTYPE_FMT                = 9999 # TODO



class FASTInputFile(File):
    """ 
    Read/write an OpenFAST input file. The object behaves like a dictionary.
    A generic reader/writer is used at first.
    If a dedicated OpenFAST input file is detected, additional functionalities are added.
    See at the end of this file for dedicated class that can be used instead of this generic reader.

    Main methods
    ------------
    - read, write, toDataFrame, keys, toGraph


    Return an object which inherits from FASTInputFileBase
     - The generic file reader is run first
     - If a specific file format/module is detected, a fixed file format object is returned
       The fixed file format have additional outputs, sanity checks and methods
    """

    @staticmethod
    def defaultExtensions():
        return ['.dat','.fst','.txt','.fstf','.dvr']

    @staticmethod
    def formatName():
        return 'FAST input file'

    def __init__(self, filename=None, **kwargs):
        self._fixedfile = None
        self.basefile = FASTInputFileBase(filename, **kwargs) # Generic fileformat

    @property
    def fixedfile(self):
        if self._fixedfile is not None:
            return self._fixedfile
        elif len(self.basefile.data)>0:
            self._fixedfile=self.fixedFormat()
            return self._fixedfile
        else:
            return self.basefile

    @property
    def module(self):
        if self._fixedfile is None:
            return self.basefile.module
        else:
            return self._fixedfile.module

    @property
    def hasNodal(self):
        if self._fixedfile is None:
            return self.basefile.hasNodal
        else:
            return self._fixedfile.hasNodal

    def getID(self, label):
        return self.basefile.getID(label)

    @property
    def data(self):
        return self.basefile.data

    def fixedFormat(self):
        # --- Creating a dedicated Child
        KEYS = list(self.basefile.keys())
        if 'NumBlNds' in KEYS:
            return ADBladeFile.from_fast_input_file(self.basefile)
        elif 'rhoinf' in KEYS:
            return BDFile.from_fast_input_file(self.basefile)
        elif 'NBlInpSt' in KEYS:
            return EDBladeFile.from_fast_input_file(self.basefile)
        elif 'NTwInpSt' in KEYS:
            return EDTowerFile.from_fast_input_file(self.basefile)
        elif 'MassMatrix' in KEYS and self.module == 'ExtPtfm':
            return ExtPtfmFile.from_fast_input_file(self.basefile)
        elif 'NumCoords' in KEYS and 'InterpOrd' in KEYS:
            return ADPolarFile.from_fast_input_file(self.basefile)
        else:
            # TODO: HD, SD, SvD, ED, AD, EDbld, BD,
            #print('>>>>>>>>>>>> NO FILEFORMAT', KEYS)
            return self.basefile

    def read(self, filename=None):
        return self.fixedfile.read(filename)

    def write(self, filename=None):
        return self.fixedfile.write(filename)

    def toDataFrame(self):
        return self.fixedfile.toDataFrame()

    def toString(self):
        return self.fixedfile.toString()

    def keys(self):
        return self.fixedfile.keys()

    def toGraph(self, **kwargs):
        return self.fixedfile.toGraph(**kwargs)

    @property
    def filename(self):
        return self.fixedfile.filename

    @property
    def comment(self):
        return self.fixedfile.comment

    @comment.setter
    def comment(self,comment):
        self.fixedfile.comment = comment

    def __iter__(self):
        return self.fixedfile.__iter__()

    def __next__(self): 
        return self.fixedfile.__next__()

    def __setitem__(self,key,item):
        return self.fixedfile.__setitem__(key,item)

    def __getitem__(self,key):
        return self.fixedfile.__getitem__(key)

    def __repr__(self):
        return self.fixedfile.__repr__()
        #s ='Fast input file: {}\n'.format(self.filename)
        #return s+'\n'.join(['{:15s}: {}'.format(d['label'],d['value']) for i,d in enumerate(self.data)])


# --------------------------------------------------------------------------------}
# --- BASE INPUT FILE 
# --------------------------------------------------------------------------------{
class FASTInputFileBase(File):
    """ 
    Read/write an OpenFAST input file. The object behaves like a dictionary.

    Main methods
    ------------
    - read, write, toDataFrame, keys

    Main keys
    ---------
    The keys correspond to the keys used in the file. For instance for a .fst file: 'DT','TMax'

    Examples
    --------

        filename = 'AeroDyn.dat'
        f = FASTInputFile(filename)
        f['TwrAero'] = True
        f['AirDens'] = 1.225
        f.write('AeroDyn_Changed.dat')

    """
    @staticmethod
    def defaultExtensions():
        return ['.dat','.fst','.txt','.fstf','.dvr']

    @staticmethod
    def formatName():
        return 'FAST input file Base'

    def __init__(self, filename=None, **kwargs):
        self._size=None
        self.setData() # Init data
        if filename:
            self.filename = filename
            self.read()

    def setData(self, filename=None, data=None, hasNodal=False, module=None):
        """ Set the data of this object. This object shouldn't store anything else. """
        if data is None:
            self.data     = []
        else:
            self.data     = data
        self.hasNodal = hasNodal
        self.module   = module
        self.filename = filename

    def keys(self):
        self.labels = [ d['label'] for i,d in enumerate(self.data) if (not d['isComment']) and (i not in self._IComment)]
        return self.labels

    def getID(self,label):
        i=self.getIDSafe(label)
        if i<0:
            raise KeyError('Variable `'+ label+'` not found in FAST file:'+self.filename)
        else:
            return i
    def getIDs(self,label):
        I=[]
        # brute force search
        for i in range(len(self.data)):
            d = self.data[i]
            if d['label'].lower()==label.lower():
                I.append(i)
        if len(I)<0:
            raise KeyError('Variable `'+ label+'` not found in FAST file:'+self.filename)
        else:
            return I

    def getIDSafe(self,label):
        # brute force search
        for i in range(len(self.data)):
            d = self.data[i]
            if d['label'].lower()==label.lower():
                return i
        return -1

    # Making object an iterator
    def __iter__(self):
        self.iCurrent=-1
        self.iMax=len(self.data)-1
        return self

    def __next__(self): # Python 2: def next(self)
        if self.iCurrent > self.iMax:
            raise StopIteration
        else:
            self.iCurrent += 1
            return self.data[self.iCurrent]

    # Making it behave like a dictionary
    def __setitem__(self, key, item):
        I = self.getIDs(key)
        for i in I: 
            if self.data[i]['tabType'] != TABTYPE_NOT_A_TAB:
                # For tables, we automatically update variable that stores the dimension 
                nRows   = len(item)
                if 'tabDimVar' in self.data[i].keys():
                    dimVar  = self.data[i]['tabDimVar']
                    iDimVar = self.getID(dimVar)
                    self.data[iDimVar]['value'] = nRows # Avoiding a recursive call to __setitem__ here
                else:
                    pass
            self.data[i]['value'] = item

    def __getitem__(self,key):
        i = self.getID(key)
        return self.data[i]['value']

    def __repr__(self):
        s ='Fast input file base: {}\n'.format(self.filename)
        return s+'\n'.join(['{:15s}: {}'.format(d['label'],d['value']) for i,d in enumerate(self.data)])

    def addKeyVal(self, key, val, descr=None):
        i=self.getIDSafe(key)
        if i<0:
            d = getDict()
        else:
            d = self.data[i]
        d['label']=key
        d['value']=val
        if descr is not None:
            d['descr']=descr
        if i<0:
            self.data.append(d)

    def addValKey(self,val,key,descr=None):
        self.addKeyVal(key, val, descr)

    def addComment(self, comment='!'):
        d=getDict()
        d['isComment'] = True
        d['value']     = comment
        self.data.append(d)

    def addTable(self, label, tab, cols=None, units=None, tabType=1, tabDimVar=None):
        d=getDict()
        d['label']          = label
        d['value']          = tab
        d['tabType']        = tabType
        d['tabDimVar']      = tabDimVar
        d['tabColumnNames'] = cols
        d['tabUnits']       = units
        self.data.append(d)

    @property
    def comment(self):
        return '\n'.join([self.data[i]['value'] for i in self._IComment])

    @comment.setter
    def comment(self, comment):
        splits = comment.split('\n')
        for i,com in zip(self._IComment, splits):
            self.data[i]['value'] = com
            self.data[i]['label'] = ''
            self.data[i]['descr'] = ''
            self.data[i]['isComment'] = True

    @property
    def _IComment(self):
        """ return indices of comment line"""
        return [1] # Typical OpenFAST files have comment on second line [1]


    def read(self, filename=None):
        if filename:
            self.filename = filename
        if self.filename:
            if not os.path.isfile(self.filename):
                raise OSError(2,'File not found:',self.filename)
            if os.stat(self.filename).st_size == 0:
                raise EmptyFileError('File is empty:',self.filename)
            self._read()
        else:  
            raise Exception('No filename provided')

    def _read(self):

        # --- Tables that can be detected based on the "Value" (first entry on line)
        # TODO members for  BeamDyn with mutliple key point                                                                                                                                                                                                                                                                                                        ####### TODO PropSetID is Duplicate SubDyn and used in HydroDyn
        NUMTAB_FROM_VAL_DETECT  = ['HtFract'  , 'TwrElev'   , 'BlFract'  , 'Genspd_TLU' , 'BlSpn'        , 'HvCoefID' , 'AxCoefID' , 'JointID'  , 'Dpth'      , 'FillNumM'    , 'MGDpth'    , 'SimplCd'  , 'RNodes'       , 'kp_xr'      , 'mu1'           , 'TwrHtFr'   , 'TwrRe'  , 'WT_X']
        NUMTAB_FROM_VAL_DIM_VAR = ['NTwInpSt' , 'NumTwrNds' , 'NBlInpSt' , 'DLL_NumTrq' , 'NumBlNds'     , 'NHvCoef'  , 'NAxCoef'  , 'NJoints'  , 'NCoefDpth' , 'NFillGroups' , 'NMGDepths' , 1          , 'BldNodes'     , 'kp_total'   , 1               , 'NTwrHt'    , 'NTwrRe' , 'NumTurbines']
        NUMTAB_FROM_VAL_VARNAME = ['TowProp'  , 'TowProp'   , 'BldProp'  , 'DLLProp'    , 'BldAeroNodes' , 'HvCoefs'  , 'AxCoefs'  , 'Joints'   , 'DpthProp'  , 'FillGroups'  , 'MGProp'    , 'SmplProp' , 'BldAeroNodes' , 'MemberGeom' , 'DampingCoeffs' , 'TowerProp' , 'TowerRe', 'WindTurbines']
        NUMTAB_FROM_VAL_NHEADER = [2          , 2           , 2          , 2            , 2              , 2          , 2          , 2          , 2           , 2             , 2           , 2          , 1              , 2            , 2               , 1           , 1        , 2 ]
        NUMTAB_FROM_VAL_TYPE    = ['num'      , 'num'       , 'num'      , 'num'        , 'num'          , 'num'      , 'num'      , 'num'      , 'num'       , 'num'         , 'num'       , 'num'      , 'mix'          , 'num'        , 'num'           , 'num'       , 'num'    , 'mix']
        # SubDyn
        NUMTAB_FROM_VAL_DETECT  += [ 'RJointID'        , 'IJointID'        , 'COSMID'             , 'CMJointID'         ]
        NUMTAB_FROM_VAL_DIM_VAR += [ 'NReact'          , 'NInterf'         , 'NCOSMs'             , 'NCmass'            ]
        NUMTAB_FROM_VAL_VARNAME += [ 'BaseJoints'      , 'InterfaceJoints' , 'MemberCosineMatrix' , 'ConcentratedMasses']
        NUMTAB_FROM_VAL_NHEADER += [ 2                 , 2                 , 2                    , 2                   ]
        NUMTAB_FROM_VAL_TYPE    += [ 'mix'             , 'num'             , 'num'                , 'num'               ]
        # AD Driver old and new
        NUMTAB_FROM_VAL_DETECT  += [ 'WndSpeed' , 'HWndSpeed' ]
        NUMTAB_FROM_VAL_DIM_VAR += [ 'NumCases' , 'NumCases'  ]
        NUMTAB_FROM_VAL_VARNAME += [ 'Cases'    , 'Cases'     ]
        NUMTAB_FROM_VAL_NHEADER += [ 2          , 2           ]
        NUMTAB_FROM_VAL_TYPE    += [ 'num'      , 'num'       ]

        # --- Tables that can be detected based on the "Label" (second entry on line)
        # NOTE: MJointID1, used by SubDyn and HydroDyn
        NUMTAB_FROM_LAB_DETECT   = ['NumAlf'  , 'F_X'       , 'MemberCd1'    , 'MJointID1' , 'NOutLoc'    , 'NOutCnt'    , 'PropD'       ]
        NUMTAB_FROM_LAB_DIM_VAR  = ['NumAlf'  , 'NKInpSt'   , 'NCoefMembers' , 'NMembers'  , 'NMOutputs'  , 'NMOutputs'  , 'NPropSets'   ]
        NUMTAB_FROM_LAB_VARNAME  = ['AFCoeff' , 'TMDspProp' , 'MemberProp'   , 'Members'   , 'MemberOuts' , 'MemberOuts' , 'SectionProp' ]
        NUMTAB_FROM_LAB_NHEADER  = [2         , 2           , 2              , 2           , 2            , 2            , 2             ]
        NUMTAB_FROM_LAB_NOFFSET  = [0         , 0           , 0              , 0           , 0            , 0            , 0             ]
        NUMTAB_FROM_LAB_TYPE     = ['num'     , 'num'       , 'num'          , 'mix'       , 'num'        , 'sdout'      , 'num'         ]
        # MoorDyn Version 1 and 2 (with AUTO for LAB_DIM_VAR)
        NUMTAB_FROM_LAB_DETECT   += ['Diam'       ,'Type'           ,'LineType'    , 'Attachment']
        NUMTAB_FROM_LAB_DIM_VAR  += ['NTypes:AUTO','NConnects'      ,'NLines:AUTO' , 'AUTO']
        NUMTAB_FROM_LAB_VARNAME  += ['LineTypes'  ,'ConnectionProp' ,'LineProp'    , 'Points']
        NUMTAB_FROM_LAB_NHEADER  += [ 2           , 2               , 2            , 2     ]
        NUMTAB_FROM_LAB_NOFFSET  += [ 0           , 0               , 0            , 0     ]
        NUMTAB_FROM_LAB_TYPE     += ['mix'        ,'mix'            ,'mix'         , 'mix']
        # SubDyn
        NUMTAB_FROM_LAB_DETECT   += ['GuyanDampSize'     , 'YoungE'   , 'YoungE'    , 'EA'             , 'MatDens'       ]
        NUMTAB_FROM_LAB_DIM_VAR  += [6                   , 'NPropSets', 'NXPropSets', 'NCablePropSets' , 'NRigidPropSets']
        NUMTAB_FROM_LAB_VARNAME  += ['GuyanDampMatrix'   , 'BeamProp' , 'BeamPropX' , 'CableProp'      , 'RigidProp'     ]
        NUMTAB_FROM_LAB_NHEADER  += [0                   , 2          , 2           , 2                , 2               ]
        NUMTAB_FROM_LAB_NOFFSET  += [1                   , 0          , 0           , 0                , 0               ]
        NUMTAB_FROM_LAB_TYPE     += ['num'               , 'num'      , 'num'       , 'num'            , 'num'           ]
        # OLAF
        NUMTAB_FROM_LAB_DETECT   += ['GridName'   ]
        NUMTAB_FROM_LAB_DIM_VAR  += ['nGridOut'   ]
        NUMTAB_FROM_LAB_VARNAME  += ['GridOutputs']
        NUMTAB_FROM_LAB_NHEADER  += [0            ]
        NUMTAB_FROM_LAB_NOFFSET  += [2            ]
        NUMTAB_FROM_LAB_TYPE     += ['mix'        ]

        FILTAB_FROM_LAB_DETECT   = ['FoilNm' ,'AFNames']
        FILTAB_FROM_LAB_DIM_VAR  = ['NumFoil','NumAFfiles']
        FILTAB_FROM_LAB_VARNAME  = ['FoilNm' ,'AFNames']

        # Using lower case to be more tolerant..
        NUMTAB_FROM_VAL_DETECT_L = [s.lower() for s in NUMTAB_FROM_VAL_DETECT]
        NUMTAB_FROM_LAB_DETECT_L = [s.lower() for s in NUMTAB_FROM_LAB_DETECT]                                         
        FILTAB_FROM_LAB_DETECT_L = [s.lower() for s in FILTAB_FROM_LAB_DETECT]

        # Reset data
        self.data   = []
        self.hasNodal=False
        self.module = None
        #with open(self.filename, 'r', errors="surrogateescape") as f:
        with open(self.filename, 'r', errors="surrogateescape") as f:
            lines=f.read().splitlines()
        # IF NEEDED> DO THE FOLLOWING FORMATTING:
            #lines = [str(l).encode('utf-8').decode('ascii','ignore') for l in f.read().splitlines()]

        # Fast files start with ! or -
        #if lines[0][0]!='!' and lines[0][0]!='-':
        #    raise Exception('Fast file do not start with ! or -, is it the right format?')

        # Special filetypes
        if detectAndReadExtPtfmSE(self, lines):
            return
        if self.detectAndReadAirfoilAD14(lines):
            return

        # Parsing line by line, storing each line into a dictionary
        i=0    
        nComments  = 0
        nWrongLabels = 0
        allowSpaceSeparatedList=False
        iTab = 0

        labOffset=''
        while i<len(lines):
            line = lines[i]

            # --- Read special sections
            if line.upper().find('ADDITIONAL OUTPUTS')>0 \
            or line.upper().find('MESH-BASED OUTPUTS')>0 \
            or line.upper().find('OUTPUT CHANNELS'   )>0: # "OutList - The next line(s) contains a list of output parameters. See OutListParameters.xlsx for a listing of available output channels, (-)'"
                # TODO, lazy implementation so far, MAKE SUB FUNCTION
                parts = re.match(r'^\W*\w+', line)
                if parts:
                    firstword = parts.group(0).strip()
                else:
                    raise NotImplementedError
                remainer  = re.sub(r'^\W*\w+\W*', '', line)
                # Parsing outlist, and then we continue at a new "i" (to read END etc.)
                OutList,i = parseFASTOutList(lines,i+1) 
                d = getDict()
                if self.hasNodal and not firstword.endswith('_Nodal'):
                    d['label']   = firstword+'_Nodal'
                else:
                    d['label']   = firstword
                d['descr']   = remainer
                d['tabType'] = TABTYPE_FIL # TODO
                d['value']   = ['']+OutList
                self.data.append(d)
                if i>=len(lines):
                    break

                # --- Here we cheat and force an exit of the input file
                # The reason for this is that some files have a lot of things after the END, which will result in the file being intepreted as a wrong format due to too many comments
                if i+2<len(lines) and (lines[i+2].lower().find('bldnd_bladesout')>0 or lines[i+2].lower().find('bldnd_bloutnd')>0):
                    self.hasNodal=True
                else:
                    self.data.append(parseFASTInputLine('END of input file (the word "END" must appear in the first 3 columns of this last OutList line)',i+1))
                    self.data.append(parseFASTInputLine('---------------------------------------------------------------------------------------',i+2))
                    break
            elif line.upper().find('SSOUTLIST'   )>0 or line.upper().find('SDOUTLIST'   )>0:
                # SUBDYN Outlist doesn not follow regular format
                self.data.append(parseFASTInputLine(line,i))
                # OUTLIST Exception for BeamDyn
                OutList,i = parseFASTOutList(lines,i+1) 
                # TODO
                for o in OutList:
                    d = getDict()
                    d['isComment'] = True
                    d['value']=o
                    self.data.append(d)
                # --- Here we cheat and force an exit of the input file
                self.data.append(parseFASTInputLine('END of input file (the word "END" must appear in the first 3 columns of this last OutList line)',i+1))
                self.data.append(parseFASTInputLine('---------------------------------------------------------------------------------------',i+2))
                break
            elif line.upper().find('ADDITIONAL STIFFNESS')>0:
                # TODO, lazy implementation so far, MAKE SUB FUNCTION
                self.data.append(parseFASTInputLine(line,i))
                i +=1
                KDAdd = []
                for _ in range(19):
                    KDAdd.append(lines[i])
                    i +=1
                d = getDict()
                d['label']   = 'KDAdd'   # TODO
                d['tabType'] = TABTYPE_FIL # TODO
                d['value']   = KDAdd
                self.data.append(d)
                if i>=len(lines):
                    break
            elif line.upper().find('DISTRIBUTED PROPERTIES')>0:
                self.data.append(parseFASTInputLine(line,i));
                i+=1;
                self.readBeamDynProps(lines,i)
                return
            elif line.upper().find('OUTPUTS')>0:
                if 'Points' in self.keys() and 'dtM' in self.keys():
                    OutList,i = parseFASTOutList(lines,i+1) 
                    d = getDict()
                    d['label']   = 'Outlist'
                    d['descr']   = ''
                    d['tabType'] = TABTYPE_FIL # TODO
                    d['value']   = OutList
                    self.addComment('------------------------ OUTPUTS --------------------------------------------')
                    self.data.append(d)
                    self.addComment('END')
                    self.addComment('------------------------- need this line --------------------------------------')
                    return

            # --- Parsing of standard lines: value(s) key comment
            line = lines[i]
            d = parseFASTInputLine(line,i,allowSpaceSeparatedList)
            labelRaw =d['label'].lower()
            d['label']+=labOffset

            # --- Handling of special files
            if labelRaw=='kp_total':
                # BeamDyn has weird space speparated list around keypoint definition
                allowSpaceSeparatedList=True
            elif labelRaw=='numcoords':
                # TODO, lazy implementation so far, MAKE SUB FUNCTION
                if isStr(d['value']):
                    if d['value'][0]=='@':
                        # it's a ref to the airfoil coord file
                        pass
                else:
                    if not strIsInt(d['value']): 
                        raise WrongFormatError('Wrong value of NumCoords')
                    if int(d['value'])<=0:
                        pass
                    else:
                        self.data.append(d); i+=1;
                        # 3 comment lines
                        self.data.append(parseFASTInputLine(lines[i],i)); i+=1;
                        self.data.append(parseFASTInputLine(lines[i],i)); i+=1;
                        self.data.append(parseFASTInputLine(lines[i],i)); i+=1;
                        splits=cleanAfterChar(cleanLine(lines[i]),'!').split()
                        # Airfoil ref point
                        try:
                            pos=[float(splits[0]), float(splits[1])]
                        except:
                            raise WrongFormatError('Wrong format while reading coordinates of airfoil reference')
                        i+=1
                        d = getDict()
                        d['label'] = 'AirfoilRefPoint'
                        d['value'] = pos
                        self.data.append(d)
                        # 2 comment lines
                        self.data.append(parseFASTInputLine(lines[i],i)); i+=1;
                        self.data.append(parseFASTInputLine(lines[i],i)); i+=1;
                        # Table of coordinats itself
                        d = getDict()
                        d['label']     = 'AirfoilCoord'
                        d['tabDimVar'] = 'NumCoords'
                        d['tabType']   = TABTYPE_NUM_WITH_HEADERCOM
                        nTabLines = self[d['tabDimVar']]-1  # SOMEHOW ONE DATA POINT LESS
                        d['value'], d['tabColumnNames'],_  = parseFASTNumTable(self.filename,lines[i:i+nTabLines+1],nTabLines,i,1)
                        d['tabUnits'] = ['(-)','(-)']
                        self.data.append(d)
                        break

            elif labelRaw=='re':
                try:
                    nAirfoilTab = self['NumTabs']
                    iTab +=1
                    if nAirfoilTab>1:
                        labOffset ='_'+str(iTab)
                    d['label']=labelRaw+labOffset
                except:
                    # Unsteady driver input file...
                    pass


            #print('label>',d['label'],'<',type(d['label']));
            #print('value>',d['value'],'<',type(d['value']));
            #print(isStr(d['value']))
            #if isStr(d['value']):
            #    print(d['value'].lower() in NUMTAB_FROM_VAL_DETECT_L)

                
            # --- Handling of tables
            if isStr(d['value']) and d['value'].lower() in NUMTAB_FROM_VAL_DETECT_L:
                # Table with numerical values, 
                ii             = NUMTAB_FROM_VAL_DETECT_L.index(d['value'].lower())
                tab_type       = NUMTAB_FROM_VAL_TYPE[ii]
                if tab_type=='num':
                    d['tabType']   = TABTYPE_NUM_WITH_HEADER
                else:
                    d['tabType']   = TABTYPE_MIX_WITH_HEADER
                d['label']     = NUMTAB_FROM_VAL_VARNAME[ii]+labOffset
                d['tabDimVar'] = NUMTAB_FROM_VAL_DIM_VAR[ii]
                nHeaders       = NUMTAB_FROM_VAL_NHEADER[ii]
                nTabLines=0
                if isinstance(d['tabDimVar'],int):
                    nTabLines = d['tabDimVar']
                else:
                    nTabLines = self[d['tabDimVar']]
                #print('Reading table {} Dimension {} (based on {})'.format(d['label'],nTabLines,d['tabDimVar']));
                d['value'], d['tabColumnNames'], d['tabUnits'] = parseFASTNumTable(self.filename,lines[i:i+nTabLines+nHeaders], nTabLines, i, nHeaders, tableType=tab_type, varNumLines=d['tabDimVar'])
                i += nTabLines+nHeaders-1

                # --- Temporary hack for e.g. SubDyn, that has duplicate table, impossible to detect in the current way...
                # So we remove the element form the list one read
                del NUMTAB_FROM_VAL_DETECT[ii]  
                del NUMTAB_FROM_VAL_DIM_VAR[ii] 
                del NUMTAB_FROM_VAL_VARNAME[ii] 
                del NUMTAB_FROM_VAL_NHEADER[ii] 
                del NUMTAB_FROM_VAL_TYPE   [ii] 
                del NUMTAB_FROM_VAL_DETECT_L[ii]  

            elif isStr(labelRaw) and labelRaw in NUMTAB_FROM_LAB_DETECT_L:
                ii      = NUMTAB_FROM_LAB_DETECT_L.index(labelRaw)
                tab_type       = NUMTAB_FROM_LAB_TYPE[ii]
                # Special case for airfoil data, the table follows NumAlf, so we add d first
                doDelete =True
                if labelRaw=='numalf':
                    doDelete =False
                    d['tabType']=TABTYPE_NOT_A_TAB
                    self.data.append(d)
                    # Creating a new dictionary for the table
                    d = {'value':None, 'label':'NumAlf'+labOffset, 'isComment':False, 'descr':'', 'tabType':None}
                    i += 1
                nHeaders       = NUMTAB_FROM_LAB_NHEADER[ii]
                nOffset        = NUMTAB_FROM_LAB_NOFFSET[ii]
                if nOffset>0:
                    # Creating a dictionary for that entry
                    dd = {'value':d['value'], 'label':d['label']+labOffset, 'isComment':False, 'descr':d['descr'], 'tabType':TABTYPE_NOT_A_TAB}
                    self.data.append(dd)

                d['label']     = NUMTAB_FROM_LAB_VARNAME[ii]
                if d['label'].lower()=='afcoeff' :
                    d['tabType']        = TABTYPE_NUM_WITH_HEADERCOM
                else:
                    if tab_type=='num':
                        d['tabType']   = TABTYPE_NUM_WITH_HEADER
                    elif tab_type=='sdout':
                        d['tabType']   = TABTYPE_NUM_SUBDYNOUT
                    else:
                        d['tabType']   = TABTYPE_MIX_WITH_HEADER
                # Finding table dimension (number of lines)
                tabDimVar = NUMTAB_FROM_LAB_DIM_VAR[ii]
                if isinstance(tabDimVar, int): # dimension hardcoded
                    d['tabDimVar'] = tabDimVar
                    nTabLines = d['tabDimVar']
                else:
                    # We either use a variable name or "AUTO" to find the number of rows
                    tabDimVars = tabDimVar.split(':')
                    for tabDimVar in tabDimVars:
                        d['tabDimVar'] = tabDimVar
                        if tabDimVar=='AUTO':
                            # Determine table dimension automatically
                            nTabLines = findNumberOfTableLines(lines[i+nHeaders:], break_chars=['---','!','#'])
                            break
                        else:
                            try:
                                nTabLines = self[tabDimVar+labOffset]
                                break
                            except KeyError:
                                #print('Cannot determine table dimension using {}'.format(tabDimVar))
                                # Hopefully this table has AUTO as well
                                pass

                d['label']  += labOffset
                #print('Reading table {} Dimension {} (based on {})'.format(d['label'],nTabLines,d['tabDimVar']));
                d['value'], d['tabColumnNames'], d['tabUnits'] = parseFASTNumTable(self.filename,lines[i:i+nTabLines+nHeaders+nOffset],nTabLines,i, nHeaders, tableType=tab_type, nOffset=nOffset, varNumLines=d['tabDimVar'])
                i += nTabLines+1-nOffset

                # --- Temporary hack for e.g. SubDyn, that has duplicate table, impossible to detect in the current way...
                # So we remove the element form the list one read
                if doDelete:
                    del NUMTAB_FROM_LAB_DETECT[ii]  
                    del NUMTAB_FROM_LAB_DIM_VAR[ii] 
                    del NUMTAB_FROM_LAB_VARNAME[ii] 
                    del NUMTAB_FROM_LAB_NHEADER[ii] 
                    del NUMTAB_FROM_LAB_NOFFSET[ii] 
                    del NUMTAB_FROM_LAB_TYPE   [ii] 
                    del NUMTAB_FROM_LAB_DETECT_L[ii]  

            elif isStr(d['label']) and d['label'].lower() in FILTAB_FROM_LAB_DETECT_L:
                ii             = FILTAB_FROM_LAB_DETECT_L.index(d['label'].lower())
                d['label']     = FILTAB_FROM_LAB_VARNAME[ii]+labOffset
                d['tabDimVar'] = FILTAB_FROM_LAB_DIM_VAR[ii]
                d['tabType']   = TABTYPE_FIL
                nTabLines = self[d['tabDimVar']]
                #print('Reading table {} Dimension {} (based on {})'.format(d['label'],nTabLines,d['tabDimVar']));
                d['value'] = parseFASTFilTable(lines[i:i+nTabLines],nTabLines,i)
                i += nTabLines-1



            self.data.append(d)
            i += 1
            # --- Safety checks
            if d['isComment']:
                #print(line)
                nComments +=1
            else:
                if hasSpecialChars(d['label']):
                    nWrongLabels +=1
                    #print('label>',d['label'],'<',type(d['label']),line);
                    if i>3: # first few lines may be comments, we allow it
                        #print('Line',i,'Label:',d['label'])
                        raise WrongFormatError('Special Character found in Label: `{}`, for line: `{}`'.format(d['label'],line))
                if len(d['label'])==0:
                    nWrongLabels +=1
            if nComments>len(lines)*0.35:
                #print('Comment fail',nComments,len(lines),self.filename)
                raise WrongFormatError('Most lines were read as comments, probably not a FAST Input File: {}'.format(self.filename))
            if nWrongLabels>len(lines)*0.10:
                #print('Label fail',nWrongLabels,len(lines),self.filename)
                raise WrongFormatError('Too many lines with wrong labels, probably not a FAST Input File {}:'.format(self.filename))

            # --- END OF FOR LOOP ON LINES

        # --- PostReading checks
        labels = self.keys()
        duplicates = set([x for x in labels if (labels.count(x) > 1) and x!='OutList' and x.strip()!='-'])
        if len(duplicates)>0:
            print('[WARN] Duplicate labels found in file: '+self.filename)
            print('       Duplicates: '+', '.join(duplicates))
            print('       It\'s strongly recommended to make them unique! ')
#         except WrongFormatError as e:    
#             raise WrongFormatError('Fast File {}: '.format(self.filename)+'\n'+e.args[0])
#         except Exception as e:    
#             raise e
# #             print(e)
#             raise Exception('Fast File {}: '.format(self.filename)+'\n'+e.args[0])
        self._lines = lines 

            
    def toString(self):
        s=''
        # Special file formats, TODO subclass
        def toStringVLD(val,lab,descr):
            val='{}'.format(val)
            lab='{}'.format(lab)
            if len(val)<13:
                val='{:13s}'.format(val)
            if len(lab)<13:
                lab='{:13s}'.format(lab)
            return val+' '+lab+' - '+descr.strip().lstrip('-').lstrip()

        def toStringIntFloatStr(x):
            try:
                if int(x)==x:
                    s='{:15.0f}'.format(x)
                else:
                    s='{:15.8e}'.format(x)
            except:
                s=x
            return s

        def beamdyn_section_mat_tostring(x,K,M):
            def mat_tostring(M,fmt='24.16e'):
                return '\n'.join(['   '+' '.join(['{:24.16E}'.format(m) for m in M[i,:]]) for i in range(np.size(M,1))])
            s=''
            s+='{:.6f}\n'.format(x)
            s+=mat_tostring(K)
            #s+=np.array2string(K)
            s+='\n'
            s+='\n'
            s+=mat_tostring(M)
            #s+=np.array2string(M)
            s+='\n'
            s+='\n'
            return s

        for i in range(len(self.data)):
            d=self.data[i]
            if d['isComment']:
                s+='{}'.format(d['value'])
            elif d['tabType']==TABTYPE_NOT_A_TAB:
                if isinstance(d['value'], list):
                    sList=', '.join([str(x) for x in d['value']])
                    s+=toStringVLD(sList, d['label'], d['descr'])
                else:
                    s+=toStringVLD(d['value'],d['label'],d['descr'])
            elif d['tabType']==TABTYPE_NUM_WITH_HEADER:
                if d['tabColumnNames'] is not None:
                    s+='{}'.format(' '.join(['{:15s}'.format(s) for s in d['tabColumnNames']]))
                #s+=d['descr'] # Not ready for that
                    if d['tabUnits'] is not None:
                        s+='\n'
                        s+='{}'.format(' '.join(['{:15s}'.format(s) for s in d['tabUnits']]))
                    newline='\n'
                else:
                    newline=''
                if np.size(d['value'],0) > 0 :
                    s+=newline
                    s+='\n'.join('\t'.join( ('{:15.0f}'.format(x) if int(x)==x else '{:15.8e}'.format(x) )  for x in y) for y in d['value'])
            elif d['tabType']==TABTYPE_MIX_WITH_HEADER:
                s+='{}'.format(' '.join(['{:15s}'.format(s) for s in d['tabColumnNames']]))
                if d['tabUnits'] is not None:
                    s+='\n'
                    s+='{}'.format(' '.join(['{:15s}'.format(s) for s in d['tabUnits']]))
                if np.size(d['value'],0) > 0 :
                    s+='\n'
                    s+='\n'.join('\t'.join(toStringIntFloatStr(x) for x in y) for y in d['value'])
            elif d['tabType']==TABTYPE_NUM_WITH_HEADERCOM:
                s+='! {}\n'.format(' '.join(['{:15s}'.format(s) for s in d['tabColumnNames']]))
                s+='! {}\n'.format(' '.join(['{:15s}'.format(s) for s in d['tabUnits']]))
                s+='\n'.join('\t'.join('{:15.8e}'.format(x) for x in y) for y in d['value'])
            elif d['tabType']==TABTYPE_FIL:
                #f.write('{} {} {}\n'.format(d['value'][0],d['tabDetect'],d['descr']))
                label = d['label']
                if 'kbot' in self.keys(): # Moordyn has no 'OutList' label..
                    label=''
                if len(d['value'])==1:
                    s+='{} {} {}'.format(d['value'][0], label, d['descr']) # TODO?
                else:
                    s+='{} {} {}\n'.format(d['value'][0], label, d['descr']) # TODO?
                    s+='\n'.join(fil for fil in d['value'][1:])
            elif d['tabType']==TABTYPE_NUM_BEAMDYN:
                # TODO use dedicated sub-class
                data = d['value']
                Cols =['Span'] 
                Cols+=['K{}{}'.format(i+1,j+1) for i in range(6) for j in range(6)] 
                Cols+=['M{}{}'.format(i+1,j+1) for i in range(6) for j in range(6)] 
                for i in np.arange(len(data['span'])):
                    x = data['span'][i]
                    K = data['K'][i]
                    M = data['M'][i]
                    s += beamdyn_section_mat_tostring(x,K,M)
            elif d['tabType']==TABTYPE_NUM_SUBDYNOUT:
                data = d['value']
                s+='{}\n'.format(' '.join(['{:15s}'.format(s) for s in d['tabColumnNames']]))
                s+='{}'.format(' '.join(['{:15s}'.format(s) for s in d['tabUnits']]))
                if np.size(d['value'],0) > 0 :
                    s+='\n'
                    s+='\n'.join('\t'.join('{:15.0f}'.format(x) for x in y) for y in data)
            else:
                raise Exception('Unknown table type for variable {}'.format(d))
            if i<len(self.data)-1:
                s+='\n'
        return s

    def write(self, filename=None):
        if filename:
            self.filename = filename
        if self.filename:
            dirname = os.path.dirname(self.filename)
            if not os.path.exists(dirname) and len(dirname)>0:
                print('[WARN] Creating directory: ',dirname)
                os.makedirs(dirname)

            self._write()
        else:
            raise Exception('No filename provided')

    def _writeSanityChecks(self):
        """ Sanity checks before write"""
        pass

    def _write(self):
        self._writeSanityChecks()
        with open(self.filename,'w') as f:
            f.write(self.toString())

    def toDataFrame(self):
        return self._toDataFrame()

    def _toDataFrame(self):
        dfs={}

        for i in range(len(self.data)): 
            d=self.data[i]
            if d['tabType'] in [TABTYPE_NUM_WITH_HEADER, TABTYPE_NUM_WITH_HEADERCOM, TABTYPE_NUM_NO_HEADER, TABTYPE_MIX_WITH_HEADER]:
                Val= d['value']
                if d['tabUnits'] is None:
                    Cols=d['tabColumnNames']
                else:
                    Cols=['{}_{}'.format(c,u.replace('(','[').replace(')',']')) for c,u in zip(d['tabColumnNames'],d['tabUnits'])]
                #print(Val)
                #print(Cols)

                # --- Adding some useful tabulated data for some files (Shapefunctions, polar)
                  
                if self.getIDSafe('TwFAM1Sh(2)')>0:
                    # Hack for tower files, we add the modes
                    # NOTE: we provide interpolated shape function just in case the resolution of the input file is low..
                    x=Val[:,0]
                    Modes=np.zeros((x.shape[0],4))
                    Modes[:,0] = x**2 * self['TwFAM1Sh(2)'] + x**3 * self['TwFAM1Sh(3)'] + x**4 * self['TwFAM1Sh(4)'] + x**5 * self['TwFAM1Sh(5)'] + x**6 * self['TwFAM1Sh(6)']
                    Modes[:,1] = x**2 * self['TwFAM2Sh(2)'] + x**3 * self['TwFAM2Sh(3)'] + x**4 * self['TwFAM2Sh(4)'] + x**5 * self['TwFAM2Sh(5)'] + x**6 * self['TwFAM2Sh(6)']
                    Modes[:,2] = x**2 * self['TwSSM1Sh(2)'] + x**3 * self['TwSSM1Sh(3)'] + x**4 * self['TwSSM1Sh(4)'] + x**5 * self['TwSSM1Sh(5)'] + x**6 * self['TwSSM1Sh(6)']
                    Modes[:,3] = x**2 * self['TwSSM2Sh(2)'] + x**3 * self['TwSSM2Sh(3)'] + x**4 * self['TwSSM2Sh(4)'] + x**5 * self['TwSSM2Sh(5)'] + x**6 * self['TwSSM2Sh(6)']
                    Val = np.hstack((Val,Modes))
                    ShapeCols = [c+'_[-]' for c in ['ShapeForeAft1','ShapeForeAft2','ShapeSideSide1','ShapeSideSide2']]
                    Cols = Cols + ShapeCols

                name=d['label']

                if name=='DampingCoeffs':
                    pass
                else:
                    dfs[name]=pd.DataFrame(data=Val,columns=Cols)
            elif d['tabType'] in [TABTYPE_NUM_BEAMDYN]:
                span = d['value']['span']
                M    = d['value']['M']
                K    = d['value']['K']
                nSpan=len(span)
                MM=np.zeros((nSpan,1+36+36))
                MM[:,0]    = span
                MM[:,1:37] = K.reshape(nSpan,36)
                MM[:,37:]  = M.reshape(nSpan,36)
                Cols =['Span'] 
                Cols+=['K{}{}'.format(i+1,j+1) for i in range(6) for j in range(6)] 
                Cols+=['M{}{}'.format(i+1,j+1) for i in range(6) for j in range(6)] 
                # Putting the main terms first
                IAll = range(1+36+36)
                IMain= [0] + [i*6+i+1 for i in range(6)] + [i*6+i+37 for i in range(6)]
                IOrg = IMain + [i for i in range(1+36+36) if i not in IMain]
                Cols = [Cols[i] for i in IOrg]
                data = MM[:,IOrg]
                name=d['label']
                dfs[name]=pd.DataFrame(data=data,columns=Cols)
        if len(dfs)==1:
            dfs=dfs[list(dfs.keys())[0]]
        return dfs

    def toGraph(self, **kwargs):
        from .fast_input_file_graph import fastToGraph
        return fastToGraph(self, **kwargs)
        


# --------------------------------------------------------------------------------}
# --- SubReaders /detectors
# --------------------------------------------------------------------------------{


    def detectAndReadAirfoilAD14(self,lines):
        if len(lines)<14:
            return False
        # Reading number of tables
        L3 = lines[2].strip().split()
        if len(L3)<=0:
            return False
        if not strIsInt(L3[0]):
            return False
        nTables=int(L3[0])
        # Reading table ID
        L4 = lines[3].strip().split()
        if len(L4)<=nTables:
            return False
        TableID=L4[:nTables]
        if nTables==1:
            TableID=['']
        # Keywords for file format
        KW1=lines[12].strip().split()
        KW2=lines[13].strip().split()
        if len(KW1)>nTables and len(KW2)>nTables:
            if KW1[nTables].lower()=='angle' and KW2[nTables].lower()=='minimum':
                d = getDict(); d['isComment'] = True; d['value'] = lines[0]; self.data.append(d);
                d = getDict(); d['isComment'] = True; d['value'] = lines[1]; self.data.append(d);
                for i in range(2,14):
                    splits = lines[i].split()
                    #print(splits)
                    d = getDict()
                    d['label'] = ' '.join(splits[1:]) # TODO
                    d['descr'] = ' '.join(splits[1:]) # TODO
                    d['value'] = float(splits[0])
                    self.data.append(d)
                #pass
                #for i in range(2,14):
                nTabLines=0
                while 14+nTabLines<len(lines) and  len(lines[14+nTabLines].strip())>0 :
                    nTabLines +=1
                #data = np.array([lines[i].strip().split() for i in range(14,len(lines)) if len(lines[i])>0]).astype(float)
                #data = np.array([lines[i].strip().split() for i in takewhile(lambda x: len(lines[i].strip())>0, range(14,len(lines)-1))]).astype(float)
                data = np.array([lines[i].strip().split() for i in range(14,nTabLines+14)]).astype(float)
                #print(data)
                d = getDict()
                d['label']     = 'Polar'
                d['tabDimVar'] = nTabLines
                d['tabType']   = TABTYPE_NUM_NO_HEADER
                d['value']     = data
                if np.size(data,1)==1+nTables*3:
                    d['tabColumnNames'] = ['Alpha']+[n+l for l in TableID for n in ['Cl','Cd','Cm']]
                    d['tabUnits']       = ['(deg)']+['(-)' , '(-)' , '(-)']*nTables
                elif np.size(data,1)==1+nTables*2:
                    d['tabColumnNames'] = ['Alpha']+[n+l for l in TableID for n in ['Cl','Cd']]
                    d['tabUnits']       = ['(deg)']+['(-)' , '(-)']*nTables
                else:
                    d['tabColumnNames'] = ['col{}'.format(j) for j in range(np.size(data,1))]
                self.data.append(d)
                return True

    def readBeamDynProps(self,lines,iStart):
        nStations=self['station_total']
        #M=np.zeros((nStations,1+36+36))
        M    = np.zeros((nStations,6,6))
        K    = np.zeros((nStations,6,6))
        span = np.zeros(nStations)
        i=iStart;
        try:
            for j in range(nStations):
                # Read span location
                span[j]=float(lines[i]); i+=1;
                # Read stiffness matrix
                K[j,:,:]=np.array((' '.join(lines[i:i+6])).split()).astype(float).reshape(6,6)
                i+=7
                # Read mass matrix
                M[j,:,:]=np.array((' '.join(lines[i:i+6])).split()).astype(float).reshape(6,6)
                i+=7
        except: 
            raise WrongFormatError('An error occured while reading section {}/{}'.format(j+1,nStations))
        d = getDict()
        d['label']   = 'BeamProperties'
        d['descr']   = ''
        d['tabType'] = TABTYPE_NUM_BEAMDYN
        d['value']   = {'span':span, 'K':K, 'M':M}
        self.data.append(d)


# --------------------------------------------------------------------------------}
# --- Helper functions 
# --------------------------------------------------------------------------------{
def isStr(s):
   return isinstance(s, str)

def strIsFloat(s):
    #return s.replace('.',',1').isdigit()
    try:
        float(s)
        return True
    except:
        return False

def strIsBool(s):
    return s.lower() in ['true','false','t','f']

def strIsInt(s):
    s = str(s)
    if s[0] in ('-', '+'):
        return s[1:].isdigit()
    return s.isdigit()    

def strToBool(s):
    return s.lower() in ['true','t']

def hasSpecialChars(s):
    # fast allows for parenthesis
    # For now we allow for - but that's because of BeamDyn geometry members 
    return not re.match("^[\"\'a-zA-Z0-9_()-]*$", s)

def cleanLine(l):
    # makes a string single space separated
    l = l.replace('\t',' ')
    l = ' '.join(l.split())
    l = l.strip()
    return l

def cleanAfterChar(l,c):
    # remove whats after a character
    n = l.find(c);
    if n>0:
        return l[:n]
    else:
        return l

def getDict():
    return {'value':None, 'label':'', 'isComment':False, 'descr':'', 'tabType':TABTYPE_NOT_A_TAB}

def _merge_value(splits):

    merged = splits.pop(0)
    if merged[0] == '"':
        while merged[-1] != '"':
            merged += " "+splits.pop(0)
    splits.insert(0, merged)




def parseFASTInputLine(line_raw,i,allowSpaceSeparatedList=False):
    d = getDict()
    #print(line_raw)
    try:
        # preliminary cleaning (Note: loss of formatting)
        line = cleanLine(line_raw)
        # Comment
        if any(line.startswith(c) for c in ['#','!','--','==']) or len(line)==0:
            d['isComment']=True
            d['value']=line_raw
            return d
        if line.lower().startswith('end'):
            sp =line.split()
            if len(sp)>2 and sp[1]=='of':
                d['isComment']=True
                d['value']=line_raw

        # Detecting lists
        List=[];
        iComma=line.find(',')
        if iComma>0 and iComma<30:
            fakeline=line.replace(' ',',')
            fakeline=re.sub(',+',',',fakeline)
            csplits=fakeline.split(',')
            # Splitting based on comma and looping while it's numbers of booleans
            ii=0
            s=csplits[ii]
            #print(csplits)
            while strIsFloat(s) or strIsBool(s) and ii<len(csplits):
                if strIsInt(s):
                    List.append(int(s))
                elif strIsFloat(s):
                    List.append(float(s))
                elif strIsBool(s):
                    List.append(strToBool(s))
                else:
                    raise WrongFormatError('Lists of strings not supported.')
                ii =ii+1
                if ii>=len(csplits):
                    raise WrongFormatError('Wrong number of list values')
                s = csplits[ii]
            #print('[INFO] Line {}: Found list: '.format(i),List)
        # Defining value and remaining splits
        if len(List)>=2:
            d['value']=List
            line_remaining=line
            # eating line, removing each values
            for iii in range(ii):
                sValue=csplits[iii]
                ipos=line_remaining.find(sValue)
                line_remaining = line_remaining[ipos+len(sValue):]
            splits=line_remaining.split()
            iNext=0
        else:
            # It's not a list, we just use space as separators
            splits=line.split(' ')
            _merge_value(splits)
            s=splits[0]

            if strIsInt(s):
                d['value']=int(s)
                if allowSpaceSeparatedList and len(splits)>1:
                    if strIsInt(splits[1]):
                        d['value']=splits[0]+ ' '+splits[1]
            elif strIsFloat(s):
                d['value']=float(s)
            elif strIsBool(s):
                d['value']=strToBool(s)
            else:
                d['value']=s
            iNext=1

        # Extracting label (TODO, for now only second split)
        bOK=False
        while (not bOK) and iNext<len(splits):
            # Nasty handling of !XXX: comments
            if splits[iNext][0]=='!' and splits[iNext][-1]==':': 
                iNext=iNext+2
                continue
            # Nasty handling of the fact that sometimes old values are repeated before the label
            if strIsFloat(splits[iNext]):
                iNext=iNext+1
                continue
            else:
                bOK=True
        if bOK:
            d['label']= splits[iNext].strip()
            iNext = iNext+1
        else:
            #print('[WARN] Line {}: No label found -> comment assumed'.format(i+1))
            d['isComment']=True
            d['value']=line_raw
            iNext = len(splits)+1
        
        # Recombining description
        if len(splits)>=iNext+1:
            d['descr']=' '.join(splits[iNext:])
    except WrongFormatError as e:
        raise WrongFormatError('Line {}: '.format(i+1)+e.args[0])
    except Exception as e:
        raise Exception('Line {}: '.format(i+1)+e.args[0])

    return d

def parseFASTOutList(lines,iStart):
    OutList=[]
    i = iStart
    while i<len(lines) and lines[i].upper().find('END')!=0 and lines[i].upper().find('---')!=0 and lines[i].upper().find('===')!=0:
        OutList.append(lines[i]) #TODO better parsing
        #print('OutList',lines[i])
        i += 1
        if i>=len(lines):
            print('[WARN] End of file reached while reading Outlist')
    #i=min(i+1,len(lines))
    return OutList,iStart+len(OutList)


def extractWithinParenthesis(s):
    mo = re.search(r'\((.*)\)', s)
    if mo:
        return mo.group(1)
    return ''

def extractWithinBrackets(s):
    mo = re.search(r'\((.*)\)', s)
    if mo:
        return mo.group(1)
    return ''

def detectUnits(s,nRef):
    nPOpen=s.count('(')
    nPClos=s.count(')')
    nBOpen=s.count('[')
    nBClos=s.count(']')

    sep='!#@#!'
    if (nPOpen == nPClos) and (nPOpen>=nRef):
        #that's pretty good
        Units=s.replace('(','').replace(')',sep).split(sep)[:-1]
    elif (nBOpen == nBClos) and (nBOpen>=nRef):
        Units=s.replace('[','').replace(']',sep).split(sep)[:-1]
    else:
        Units=s.split()
    return Units


def findNumberOfTableLines(lines, break_chars):
    """ Loop through lines until a one of the "break character is found"""
    for i, l in enumerate(lines):
        for bc in break_chars:
            if l.startswith(bc):
                return i
    # Not found
    print('[FAIL] end of table not found')
    return len(lines)


def parseFASTNumTable(filename,lines,n,iStart,nHeaders=2,tableType='num',nOffset=0, varNumLines=''):
    """ 
    First lines of data starts at: nHeaders+nOffset
    
    """
    Tab = None
    ColNames = None
    Units = None
    

    if len(lines)!=n+nHeaders+nOffset:
        raise BrokenFormatError('Not enough lines in table: {} lines instead of {}\nFile:{}'.format(len(lines)-nHeaders,n,filename))
    try:
        if nHeaders==0:
            # Extract number of values from number of numerical values on first line
            numeric_const_pattern = r'[-+]? (?: (?: \d* \. \d+ ) | (?: \d+ \.? ) )(?: [Ee] [+-]? \d+ ) ?'
            rx = re.compile(numeric_const_pattern, re.VERBOSE)
            header = cleanAfterChar(lines[nOffset], '!')
            if tableType=='num':
                dat= np.array(rx.findall(header)).astype(float)
                ColNames=['C{}'.format(j) for j in range(len(dat))]
            else:
                raise NotImplementedError('Reading FAST tables with no headers for type different than num not implemented yet')

        elif nHeaders>=1:
            # Extract column names
            i = 0
            sTmp = cleanLine(lines[i])
            sTmp = cleanAfterChar(sTmp,'[')
            sTmp = cleanAfterChar(sTmp,'(')
            sTmp = cleanAfterChar(sTmp,'!')
            sTmp = cleanAfterChar(sTmp,'#')
            if sTmp.startswith('!'):
                sTmp=sTmp[1:].strip()
            ColNames=sTmp.split()
        if nHeaders>=2:
            # Extract units
            i = 1
            sTmp = cleanLine(lines[i])
            sTmp = cleanAfterChar(sTmp,'!')
            sTmp = cleanAfterChar(sTmp,'#')
            if sTmp.startswith('!'):
                sTmp=sTmp[1:].strip()

            Units = detectUnits(sTmp,len(ColNames))
            Units = ['({})'.format(u.strip()) for u in Units]
            # Forcing user to match number of units and column names
            if len(ColNames) != len(Units):
                print(ColNames)
                print(Units)
                print('[WARN] {}: Line {}: Number of column names different from number of units in table'.format(filename, iStart+i+1))

        nCols=len(ColNames)

        if tableType=='num':
            if n==0:
                Tab = np.zeros((n, nCols))
            for i in range(nHeaders+nOffset,n+nHeaders+nOffset):
                l = cleanAfterChar(lines[i].lower(),'!')
                l = cleanAfterChar(l,'#')
                v = l.split()
                if len(v) != nCols:
                    # Discarding SubDyn special cases
                    if ColNames[-1].lower() not in ['nodecnt']:
                        print('[WARN] {}: Line {}: number of data different from number of column names. ColumnNames: {}'.format(filename, iStart+i+1, ColNames))
                if i==nHeaders+nOffset:
                    # Node Cnt
                    if len(v) != nCols:
                        if ColNames[-1].lower()== 'nodecnt':
                            ColNames = ColNames+['Col']*(len(v)-nCols)
                            Units    = Units+['Col']*(len(v)-nCols)

                    nCols=len(v)
                    Tab = np.zeros((n, nCols))
                # Accounting for TRUE FALSE and converting to float
                v = [s.replace('true','1').replace('false','0').replace('noprint','0').replace('print','1') for s in v]
                v = [float(s) for s in v[0:nCols]]
                if len(v) < nCols:
                    raise Exception('Number of data is lower than number of column names')
                Tab[i-nHeaders-nOffset,:] = v
        elif tableType=='mix':
            # a mix table contains a mixed of strings and floats
            # For now, we are being a bit more relaxed about the number of columns
            if n==0:
                Tab = np.zeros((n, nCols)).astype(object)
            for i in range(nHeaders+nOffset,n+nHeaders+nOffset):
                l = lines[i]
                l = cleanAfterChar(l,'!')
                l = cleanAfterChar(l,'#')
                v = l.split()
                if l.startswith('---'):
                    raise BrokenFormatError('Error reading line {} while reading table. Is the variable `{}` set correctly?'.format(iStart+i+1, varNumLines))
                if len(v) != nCols:
                    # Discarding SubDyn special cases
                    if ColNames[-1].lower() not in ['cosmid', 'ssifile']:
                        print('[WARN] {}: Line {}: Number of data is different than number of column names. Column Names: {}'.format(filename,iStart+1+i, ColNames))
                if i==nHeaders+nOffset:
                    if len(v)>nCols:
                        ColNames = ColNames+['Col']*(len(v)-nCols)
                        Units    = Units+['Col']*(len(v)-nCols)
                    nCols=len(v)
                    Tab = np.zeros((n, nCols)).astype(object)
                v=v[0:min(len(v),nCols)]
                Tab[i-nHeaders-nOffset,0:len(v)] = v
            # If all values are float, we convert to float
            if all([strIsFloat(x) for x in Tab.ravel()]):
                Tab=Tab.astype(float)
        elif tableType=='sdout':
            header = lines[0]
            units  = lines[1]
            Tab=[]
            for i in range(nHeaders+nOffset,n+nHeaders+nOffset):
                l = cleanAfterChar(lines[i].lower(),'!')
                Tab.append( np.array(l.split()).astype(int))
        else:
            raise Exception('Unknown table type')

        ColNames = ColNames[0:nCols]
        if Units is not None:
            Units    = Units[0:nCols]
            Units    = ['('+u.replace('(','').replace(')','')+')' for u in Units]
        if nHeaders==0:
            ColNames=None
            
    except Exception as e:    
        raise BrokenFormatError('Line {}: {}'.format(iStart+i+1,e.args[0]))
    return Tab, ColNames, Units


def parseFASTFilTable(lines,n,iStart):
    Tab = []
    try:
        i=0
        if len(lines)!=n:
            raise WrongFormatError('Not enough lines in table: {} lines instead of {}'.format(len(lines),n))
        for i in range(n):
            l = lines[i].split()
            #print(l[0].strip())
            Tab.append(l[0].strip())
            
    except Exception as e:    
        raise Exception('Line {}: '.format(iStart+i+1)+e.args[0])
    return Tab



# --------------------------------------------------------------------------------}
# --------------------------------------------------------------------------------}
# --------------------------------------------------------------------------------}
# --- Predefined types (may change with OpenFAST version..)
# --------------------------------------------------------------------------------{
# --------------------------------------------------------------------------------{
# --------------------------------------------------------------------------------{


# --------------------------------------------------------------------------------}
# --- BeamDyn 
# --------------------------------------------------------------------------------{
class BDFile(FASTInputFileBase):
    @classmethod
    def from_fast_input_file(cls, parent):
        self = cls()
        self.setData(filename=parent.filename, data=parent.data, hasNodal=parent.hasNodal, module='BD')
        return self

    def __init__(self, filename=None, **kwargs):
        FASTInputFileBase.__init__(self, filename, **kwargs)
        if filename is None:
            # Define a prototype for this file format
            self.addComment('--------- BEAMDYN with OpenFAST INPUT FILE -------------------------------------------')
            self.addComment('BeamDyn input file, written by BDFile')
            self.addComment('---------------------- SIMULATION CONTROL --------------------------------------')
            self.addValKey(False        , 'Echo'            , 'Echo input data to "<RootName>.ech"? (flag)')
            self.addValKey(True         , 'QuasiStaticInit' , 'Use quasi-static pre-conditioning with centripetal accelerations in initialization? (flag) [dynamic solve only]')
            self.addValKey(          0  , 'rhoinf'          , 'Numerical damping parameter for generalized-alpha integrator')
            self.addValKey(          2  , 'quadrature'      , 'Quadrature method: 1=Gaussian; 2=Trapezoidal (switch)')
            self.addValKey("DEFAULT"    , 'refine'          , 'Refinement factor for trapezoidal quadrature (-) [DEFAULT = 1; used only when quadrature=2]')
            self.addValKey("DEFAULT"    , 'n_fact'          , 'Factorization frequency for the Jacobian in N-R iteration(-) [DEFAULT = 5]')
            self.addValKey("DEFAULT"    , 'DTBeam'          , 'Time step size (s)')
            self.addValKey("DEFAULT"    , 'load_retries'    , 'Number of factored load retries before quitting the simulation [DEFAULT = 20]')
            self.addValKey("DEFAULT"    , 'NRMax'           , 'Max number of iterations in Newton-Raphson algorithm (-) [DEFAULT = 10]')
            self.addValKey("DEFAULT"    , 'stop_tol'        , 'Tolerance for stopping criterion (-) [DEFAULT = 1E-5]')
            self.addValKey("DEFAULT"    , 'tngt_stf_fd'     , 'Use finite differenced tangent stiffness matrix? (flag)')
            self.addValKey("DEFAULT"    , 'tngt_stf_comp'   , 'Compare analytical finite differenced tangent stiffness matrix? (flag)')
            self.addValKey("DEFAULT"    , 'tngt_stf_pert'   , 'Perturbation size for finite differencing (-) [DEFAULT = 1E-6]')
            self.addValKey("DEFAULT"    , 'tngt_stf_difftol', 'Maximum allowable relative difference between analytical and fd tangent stiffness (-); [DEFAULT = 0.1]')
            self.addValKey(True         , 'RotStates'       , 'Orient states in the rotating frame during linearization? (flag) [used only when linearizing] ')
            self.addComment('---------------------- GEOMETRY PARAMETER --------------------------------------')
            self.addValKey(          1  , 'member_total'    , 'Total number of members (-)')
            self.addValKey(          0  , 'kp_total'        , 'Total number of key points (-) [must be at least 3]')
            self.addValKey(      [1, 0] , 'kp_per_member'   , 'Member number; Number of key points in this member')
            self.addTable('MemberGeom', np.zeros((0,4)), tabType=1, tabDimVar='kp_total', 
                    cols=['kp_xr', 'kp_yr', 'kp_zr', 'initial_twist'], 
                    units=['(m)', '(m)', '(m)', '(deg)'])
            self.addComment('---------------------- MESH PARAMETER ------------------------------------------')
            self.addValKey(          5  , 'order_elem'     , 'Order of interpolation (basis) function (-)')
            self.addComment('---------------------- MATERIAL PARAMETER --------------------------------------')
            self.addValKey('"undefined"', 'BldFile'        ,  'Name of file containing properties for blade (quoted string)')
            self.addComment('---------------------- PITCH ACTUATOR PARAMETERS -------------------------------')
            self.addValKey(False        , 'UsePitchAct'    , 'Whether a pitch actuator should be used (flag)')
            self.addValKey(          1  , 'PitchJ'         , 'Pitch actuator inertia (kg-m^2) [used only when UsePitchAct is true]')
            self.addValKey(          0  , 'PitchK'         , 'Pitch actuator stiffness (kg-m^2/s^2) [used only when UsePitchAct is true]')
            self.addValKey(          0  , 'PitchC'         , 'Pitch actuator damping (kg-m^2/s) [used only when UsePitchAct is true]')
            self.addComment('---------------------- OUTPUTS -------------------------------------------------')
            self.addValKey(False        , 'SumPrint'      , 'Print summary data to "<RootName>.sum" (flag)')
            self.addValKey('"ES10.3E2"' , 'OutFmt'        , 'Format used for text tabular output, excluding the time channel.')
            self.addValKey(          0  , 'NNodeOuts'     , 'Number of nodes to output to file [0 - 9] (-)')
            self.addValKey(         [1] , 'OutNd'         , 'Nodes whose values will be output  (-)')
            self.addValKey(        [''] , 'OutList'       , 'The next line(s) contains a list of output parameters. See OutListParameters.xlsx, BeamDyn tab for a listing of available output channels, (-)')
            self.addComment('END of OutList (the word "END" must appear in the first 3 columns of this last OutList line)')
            self.addComment('---------------------- NODE OUTPUTS --------------------------------------------')
            self.addValKey(          99 , 'BldNd_BlOutNd' , 'Blade nodes on each blade (currently unused)')
            self.addValKey(        [''] , 'OutList_Nodal' , 'The next line(s) contains a list of output parameters.  See OutListParameters.xlsx, BeamDyn_Nodes tab for a listing of available output channels, (-)')
            self.addComment('END of input file (the word "END" must appear in the first 3 columns of this last OutList line)')
            self.addComment('--------------------------------------------------------------------------------')
            self.hasNodal=True
            #"RootFxr, RootFyr, RootFzr"  
            #"RootMxr, RootMyr, RootMzr"  
            #"TipTDxr, TipTDyr, TipTDzr"  
            #"TipRDxr, TipRDyr, TipRDzr"  

        else:
            # fix some stuff that generic reader fail at
            self.data[1] =  {'value':self._lines[1], 'label':'', 'isComment':True, 'descr':'', 'tabType':0}
            i  = self.getID('kp_total')
            listval = [int(v) for v in str(self.data[i+1]['value']).split()]
            self.data[i+1]['value']=listval
            self.data[i+1]['label']='kp_per_member'
            self.data[i+1]['isComment']=False
        self.module='BD'

    def _writeSanityChecks(self):
        """ Sanity checks before write """
        self['kp_total']=self['MemberGeom'].shape[0]
        i  = self.getID('kp_total')
        self.data[i+1]['value']=[1, self['MemberGeom'].shape[0]] # kp_per_member
        self.data[i+1]['label']='kp_per_member'
        # Could check length of OutNd

    def _toDataFrame(self):
        df = FASTInputFileBase._toDataFrame(self)
        # TODO add quadrature points based on trapz/gauss
        return df

    @property
    def _IComment(self): return [1]

# --------------------------------------------------------------------------------}
# --- ElastoDyn Blade 
# --------------------------------------------------------------------------------{
class EDBladeFile(FASTInputFileBase):
    @classmethod
    def from_fast_input_file(cls, parent):
        self = cls()
        self.setData(filename=parent.filename, data=parent.data, hasNodal=parent.hasNodal, module='EDBlade')
        return self

    def __init__(self, filename=None, **kwargs):
        FASTInputFileBase.__init__(self, filename, **kwargs)
        if filename is None:
            # Define a prototype for this file format
            self.addComment('------- ELASTODYN V1.00.* INDIVIDUAL BLADE INPUT FILE --------------------------')
            self.addComment('ElastoDyn blade definition, written by EDBladeFile.')
            self.addComment('---------------------- BLADE PARAMETERS ----------------------------------------')
            self.addValKey(         0  , 'NBlInpSt'   , 'Number of blade input stations (-)')
            self.addValKey(         1. , 'BldFlDmp(1)', 'Blade flap mode #1 structural damping in percent of critical (%)')
            self.addValKey(         1. , 'BldFlDmp(2)', 'Blade flap mode #2 structural damping in percent of critical (%)')
            self.addValKey(         1. , 'BldEdDmp(1)', 'Blade edge mode #1 structural damping in percent of critical (%)')
            self.addComment('---------------------- BLADE ADJUSTMENT FACTORS --------------------------------')
            self.addValKey(         1. , 'FlStTunr(1)', 'Blade flapwise modal stiffness tuner, 1st mode (-)')
            self.addValKey(         1. , 'FlStTunr(2)', 'Blade flapwise modal stiffness tuner, 2nd mode (-)')
            self.addValKey(         1. , 'AdjBlMs'    , 'Factor to adjust blade mass density (-)')
            self.addValKey(         1. , 'AdjFlSt'    , 'Factor to adjust blade flap stiffness (-)')
            self.addValKey(         1. , 'AdjEdSt'    , 'Factor to adjust blade edge stiffness (-)')
            self.addComment('---------------------- DISTRIBUTED BLADE PROPERTIES ----------------------------')
            self.addTable('BldProp', np.zeros((0,6)), tabType=1, tabDimVar='NBlInpSt', cols=['BlFract', 'PitchAxis', 'StrcTwst', 'BMassDen', 'FlpStff', 'EdgStff'], units=['(-)', '(-)', '(deg)', '(kg/m)', '(Nm^2)', '(Nm^2)'])
            self.addComment('---------------------- BLADE MODE SHAPES ---------------------------------------')
            self.addValKey(     1.0   , 'BldFl1Sh(2)', 'Flap mode 1, coeff of x^2')
            self.addValKey(     0.0   , 'BldFl1Sh(3)', '           , coeff of x^3')
            self.addValKey(     0.0   , 'BldFl1Sh(4)', '           , coeff of x^4')
            self.addValKey(     0.0   , 'BldFl1Sh(5)', '           , coeff of x^5')
            self.addValKey(     0.0   , 'BldFl1Sh(6)', '           , coeff of x^6')
            self.addValKey(     0.0   , 'BldFl2Sh(2)', 'Flap mode 2, coeff of x^2') # NOTE: using something not too bad just incase user uses these as is..
            self.addValKey(     0.0   , 'BldFl2Sh(3)', '           , coeff of x^3')
            self.addValKey(   -13.0   , 'BldFl2Sh(4)', '           , coeff of x^4')
            self.addValKey(    27.0   , 'BldFl2Sh(5)', '           , coeff of x^5')
            self.addValKey(   -13.0   , 'BldFl2Sh(6)', '           , coeff of x^6')
            self.addValKey(     1.0   , 'BldEdgSh(2)', 'Edge mode 1, coeff of x^2')
            self.addValKey(     0.0   , 'BldEdgSh(3)', '           , coeff of x^3')
            self.addValKey(     0.0   , 'BldEdgSh(4)', '           , coeff of x^4')
            self.addValKey(     0.0   , 'BldEdgSh(5)', '           , coeff of x^5')
            self.addValKey(     0.0   , 'BldEdgSh(6)', '           , coeff of x^6')
        else:
            # fix some stuff that generic reader fail at
            self.data[1] =  {'value':self._lines[1], 'label':'', 'isComment':True, 'descr':'', 'tabType':0}
        self.module='EDBlade'

    def _writeSanityChecks(self):
        """ Sanity checks before write """
        self['NBlInpSt']=self['BldProp'].shape[0]
        # Sum of Coeffs should be 1
        for s in ['BldFl1Sh','BldFl2Sh','BldEdgSh']:
            sumcoeff=np.sum([self[s+'('+str(i)+')'] for i in [2,3,4,5,6] ])
            if np.abs(sumcoeff-1)>1e-4:
                print('[WARN] Sum of coefficients for polynomial {} not equal to 1 ({}). File: {}'.format(s, sumcoeff, self.filename))

    def _toDataFrame(self):
        df = FASTInputFileBase._toDataFrame(self)
        # We add the shape functions for EDBladeFile
        x=df['BlFract_[-]'].values
        Modes=np.zeros((x.shape[0],3))
        Modes[:,0] = x**2 * self['BldFl1Sh(2)'] + x**3 * self['BldFl1Sh(3)'] + x**4 * self['BldFl1Sh(4)'] + x**5 * self['BldFl1Sh(5)'] + x**6 * self['BldFl1Sh(6)']
        Modes[:,1] = x**2 * self['BldFl2Sh(2)'] + x**3 * self['BldFl2Sh(3)'] + x**4 * self['BldFl2Sh(4)'] + x**5 * self['BldFl2Sh(5)'] + x**6 * self['BldFl2Sh(6)']
        Modes[:,2] = x**2 * self['BldEdgSh(2)'] + x**3 * self['BldEdgSh(3)'] + x**4 * self['BldEdgSh(4)'] + x**5 * self['BldEdgSh(5)'] + x**6 * self['BldEdgSh(6)']
        df[['ShapeFlap1_[-]','ShapeFlap2_[-]','ShapeEdge1_[-]']]=Modes
        return df

    @property
    def _IComment(self): return [1]

# --------------------------------------------------------------------------------}
# --- ElastoDyn Tower 
# --------------------------------------------------------------------------------{
class EDTowerFile(FASTInputFileBase):
    @classmethod
    def from_fast_input_file(cls, parent):
        self = cls()
        self.setData(filename=parent.filename, data=parent.data, hasNodal=parent.hasNodal, module='EDTower')
        return self

    def __init__(self, filename=None, **kwargs):
        FASTInputFileBase.__init__(self, filename, **kwargs)
        if filename is None:
            # Define a prototype for this file format
            self.addComment('------- ELASTODYN V1.00.* TOWER INPUT FILE -------------------------------------')
            self.addComment('ElastoDyn tower definition, written by EDTowerFile.')
            self.addComment('---------------------- TOWER PARAMETERS ----------------------------------------')
            self.addValKey(         0  , 'NTwInpSt'      , 'Number of blade input stations (-)')
            self.addValKey(         1. , 'TwrFADmp(1)'   , 'Tower 1st fore-aft mode structural damping ratio (%)')
            self.addValKey(         1. , 'TwrFADmp(2)'   , 'Tower 2nd fore-aft mode structural damping ratio (%)')
            self.addValKey(         1. , 'TwrSSDmp(1)'   , 'Tower 1st side-to-side mode structural damping ratio (%)')
            self.addValKey(         1. , 'TwrSSDmp(2)'   , 'Tower 2nd side-to-side mode structural damping ratio (%)')
            self.addComment('---------------------- TOWER ADJUSTMENT FACTORS --------------------------------')
            self.addValKey(         1. , 'FAStTunr(1)'   , 'Tower fore-aft modal stiffness tuner, 1st mode (-)')
            self.addValKey(         1. , 'FAStTunr(2)'   , 'Tower fore-aft modal stiffness tuner, 2nd mode (-)')
            self.addValKey(         1. , 'SSStTunr(1)'   , 'Tower side-to-side stiffness tuner, 1st mode (-)')
            self.addValKey(         1. , 'SSStTunr(2)'   , 'Tower side-to-side stiffness tuner, 2nd mode (-)')
            self.addValKey(         1. , 'AdjTwMa'       , 'Factor to adjust tower mass density (-)')
            self.addValKey(         1. , 'AdjFASt'       , 'Factor to adjust tower fore-aft stiffness (-)')
            self.addValKey(         1. , 'AdjSSSt'       , 'Factor to adjust tower side-to-side stiffness (-)')
            self.addComment('---------------------- DISTRIBUTED TOWER PROPERTIES ----------------------------')
            self.addTable('TowProp', np.zeros((0,6)), tabType=1, tabDimVar='NTwInpSt', 
                    cols=['HtFract','TMassDen','TwFAStif','TwSSStif'], 
                    units=['(-)', '(kg/m)', '(Nm^2)', '(Nm^2)'])
            self.addComment('---------------------- TOWER FORE-AFT MODE SHAPES ------------------------------')
            self.addValKey(     1.0   , 'TwFAM1Sh(2)', 'Mode 1, coefficient of x^2 term')
            self.addValKey(     0.0   , 'TwFAM1Sh(3)', '      , coefficient of x^3 term')
            self.addValKey(     0.0   , 'TwFAM1Sh(4)', '      , coefficient of x^4 term')
            self.addValKey(     0.0   , 'TwFAM1Sh(5)', '      , coefficient of x^5 term')
            self.addValKey(     0.0   , 'TwFAM1Sh(6)', '      , coefficient of x^6 term')
            self.addValKey(    -26.   , 'TwFAM2Sh(2)', 'Mode 2, coefficient of x^2 term') # NOTE: using something not too bad just incase user uses these as is..
            self.addValKey(     0.0   , 'TwFAM2Sh(3)', '      , coefficient of x^3 term')
            self.addValKey(     27.   , 'TwFAM2Sh(4)', '      , coefficient of x^4 term')
            self.addValKey(     0.0   , 'TwFAM2Sh(5)', '      , coefficient of x^5 term')
            self.addValKey(     0.0   , 'TwFAM2Sh(6)', '      , coefficient of x^6 term')
            self.addComment('---------------------- TOWER SIDE-TO-SIDE MODE SHAPES --------------------------')
            self.addValKey(     1.0   , 'TwSSM1Sh(2)', 'Mode 1, coefficient of x^2 term')
            self.addValKey(     0.0   , 'TwSSM1Sh(3)', '      , coefficient of x^3 term')
            self.addValKey(     0.0   , 'TwSSM1Sh(4)', '      , coefficient of x^4 term')
            self.addValKey(     0.0   , 'TwSSM1Sh(5)', '      , coefficient of x^5 term')
            self.addValKey(     0.0   , 'TwSSM1Sh(6)', '      , coefficient of x^6 term')
            self.addValKey(    -26.   , 'TwSSM2Sh(2)', 'Mode 2, coefficient of x^2 term') # NOTE: using something not too bad just incase user uses these as is..
            self.addValKey(     0.0   , 'TwSSM2Sh(3)', '      , coefficient of x^3 term')
            self.addValKey(     27.   , 'TwSSM2Sh(4)', '      , coefficient of x^4 term')
            self.addValKey(     0.0   , 'TwSSM2Sh(5)', '      , coefficient of x^5 term')
            self.addValKey(     0.0   , 'TwSSM2Sh(6)', '      , coefficient of x^6 term')
        else:
            # fix some stuff that generic reader fail at
            self.data[1] =  {'value':self._lines[1], 'label':'', 'isComment':True, 'descr':'', 'tabType':0}
        self.module='EDTower'

    def _writeSanityChecks(self):
        """ Sanity checks before write """
        self['NTwInpSt']=self['TowProp'].shape[0]
        # Sum of Coeffs should be 1
        for s in ['TwFAM1Sh','TwFAM2Sh','TwSSM1Sh','TwSSM2Sh']:
            sumcoeff=np.sum([self[s+'('+str(i)+')'] for i in [2,3,4,5,6] ])
            if np.abs(sumcoeff-1)>1e-4:
                print('[WARN] Sum of coefficients for polynomial {} not equal to 1 ({}). File: {}'.format(s, sumcoeff, self.filename))

    def _toDataFrame(self):
        df = FASTInputFileBase._toDataFrame(self)
        # We add the shape functions for EDBladeFile
        # NOTE: we provide interpolated shape function just in case the resolution of the input file is low..
        x = df['HtFract_[-]'].values
        Modes=np.zeros((x.shape[0],4))
        Modes[:,0] = x**2 * self['TwFAM1Sh(2)'] + x**3 * self['TwFAM1Sh(3)'] + x**4 * self['TwFAM1Sh(4)'] + x**5 * self['TwFAM1Sh(5)'] + x**6 * self['TwFAM1Sh(6)']
        Modes[:,1] = x**2 * self['TwFAM2Sh(2)'] + x**3 * self['TwFAM2Sh(3)'] + x**4 * self['TwFAM2Sh(4)'] + x**5 * self['TwFAM2Sh(5)'] + x**6 * self['TwFAM2Sh(6)']
        Modes[:,2] = x**2 * self['TwSSM1Sh(2)'] + x**3 * self['TwSSM1Sh(3)'] + x**4 * self['TwSSM1Sh(4)'] + x**5 * self['TwSSM1Sh(5)'] + x**6 * self['TwSSM1Sh(6)']
        Modes[:,3] = x**2 * self['TwSSM2Sh(2)'] + x**3 * self['TwSSM2Sh(3)'] + x**4 * self['TwSSM2Sh(4)'] + x**5 * self['TwSSM2Sh(5)'] + x**6 * self['TwSSM2Sh(6)']
        ShapeCols = [c+'_[-]' for c in ['ShapeForeAft1','ShapeForeAft2','ShapeSideSide1','ShapeSideSide2']]
        df[ShapeCols]=Modes
        return df

    @property
    def _IComment(self): return [1]

# --------------------------------------------------------------------------------}
# --- AeroDyn Blade 
# --------------------------------------------------------------------------------{
class ADBladeFile(FASTInputFileBase):
    @classmethod
    def from_fast_input_file(cls, parent):
        self = cls()
        self.setData(filename=parent.filename, data=parent.data, hasNodal=parent.hasNodal, module='ADBlade')
        return self

    def __init__(self, filename=None, **kwargs):
        FASTInputFileBase.__init__(self, filename, **kwargs)
        if filename is None:
            # Define a prototype for this file format
            self.addComment('------- AERODYN BLADE DEFINITION INPUT FILE ----------------------------------------------')
            self.addComment('Aerodynamic blade definition, written by ADBladeFile')
            self.addComment('======  Blade Properties =================================================================')
            self.addKeyVal('NumBlNds', 0, 'Number of blade nodes used in the analysis (-)')
            self.addTable('BldAeroNodes', np.zeros((0,7)), tabType=1, tabDimVar='NumBlNds', cols=['BlSpn', 'BlCrvAC', 'BlSwpAC', 'BlCrvAng', 'BlTwist', 'BlChord', 'BlAFID'], units=['(m)', '(m)', '(m)', '(deg)', '(deg)', '(m)', '(-)'])
        self.module='ADBlade'

    def _writeSanityChecks(self):
        """ Sanity checks before write"""
        self['NumBlNds']=self['BldAeroNodes'].shape[0]
        aeroNodes = self['BldAeroNodes']
        # TODO double check this calculation with gradient
        dr = np.gradient(aeroNodes[:,0])
        dx = np.gradient(aeroNodes[:,1])
        crvAng = np.degrees(np.arctan2(dx,dr))
        if np.mean(np.abs(crvAng-aeroNodes[:,3]))>0.1:
            print('[WARN] BlCrvAng might not be computed correctly')

    def _toDataFrame(self):
        df = FASTInputFileBase._toDataFrame(self)
        aeroNodes = self['BldAeroNodes']
        r     = aeroNodes[:,0]
        chord = aeroNodes[:,5]
        twist = aeroNodes[:,4]*np.pi/180
        prebendAC = aeroNodes[:,1]
        sweepAC   = aeroNodes[:,2]

        # --- IEA 15
        ##'le_location: 'Leading-edge positions from a reference blade axis (usually blade pitch axis). Locations are normalized by the local chord length. Positive in -x direction for airfoil-aligned coordinate system')
        ## pitch_axis
        ##'1D array of the chordwise position of the pitch axis (0-LE, 1-TE), defined along blade span.')
        #grid     = [0.0, 0.02040816326530612, 0.04081632653061224, 0.061224489795918366, 0.08163265306122448, 0.1020408163265306, 0.12244897959183673, 0.14285714285714285, 0.16326530612244897, 0.18367346938775508, 0.2040816326530612, 0.22448979591836732, 0.24489795918367346, 0.26530612244897955, 0.2857142857142857, 0.3061224489795918, 0.32653061224489793, 0.3469387755102041, 0.36734693877551017, 0.3877551020408163, 0.4081632653061224, 0.42857142857142855, 0.44897959183673464, 0.4693877551020408, 0.4897959183673469, 0.5102040816326531, 0.5306122448979591, 0.5510204081632653, 0.5714285714285714, 0.5918367346938775, 0.6122448979591836, 0.6326530612244897, 0.6530612244897959, 0.673469387755102, 0.6938775510204082, 0.7142857142857142, 0.7346938775510203, 0.7551020408163265, 0.7755102040816326, 0.7959183673469387, 0.8163265306122448, 0.836734693877551, 0.8571428571428571, 0.8775510204081632, 0.8979591836734693, 0.9183673469387754, 0.9387755102040816, 0.9591836734693877, 0.9795918367346939, 1.0]
        #values   = [0.5045454545454545, 0.4900186808012221, 0.47270018284548393, 0.4540147730610375, 0.434647782591965, 0.4156278851950606, 0.3979378721273935, 0.38129960745617403, 0.3654920515699109, 0.35160780834472827, 0.34008443128769117, 0.3310670675965599, 0.3241031342163746, 0.3188472934612394, 0.3146895762675238, 0.311488897995355, 0.3088429219529899, 0.3066054031112312, 0.3043613335231313, 0.3018756624023877, 0.2992017656131912, 0.29648581499532917, 0.29397119399704474, 0.2918571873240831, 0.2901098902886204, 0.28880659979944606, 0.28802634398115073, 0.28784151044623507, 0.28794253614539367, 0.28852264941156663, 0.28957685074559625, 0.2911108045758606, 0.2930139151081327, 0.2952412111444283, 0.2977841397364215, 0.300565286724993, 0.3035753776130124, 0.30670446458784534, 0.30988253764299156, 0.3130107259708016, 0.31639042766652853, 0.32021109189825026, 0.32462311714967124, 0.329454188784972, 0.33463306413024474, 0.3401190402144396, 0.3460555975714659, 0.3527211856428439, 0.3600890296396286, 0.36818181818181805]
        ##'ref_axis_blade' desc='2D array of the coordinates (x,y,z) of the blade reference axis, defined along blade span. The coordinate system is the one of BeamDyn: it is placed at blade root with x pointing the suction side of the blade, y pointing the trailing edge and z along the blade span. A standard configuration will have negative x values (prebend), if swept positive y values, and positive z values.')
        #x_grid   = [0.0, 0.02040816326530612, 0.04081632653061224, 0.061224489795918366, 0.08163265306122448, 0.1020408163265306, 0.12244897959183673, 0.14285714285714285, 0.16326530612244897, 0.18367346938775508, 0.2040816326530612, 0.22448979591836732, 0.24489795918367346, 0.26530612244897955, 0.2857142857142857, 0.3061224489795918, 0.32653061224489793, 0.3469387755102041, 0.36734693877551017, 0.3877551020408163, 0.4081632653061224, 0.42857142857142855, 0.44897959183673464, 0.4693877551020408, 0.4897959183673469, 0.5102040816326531, 0.5306122448979591, 0.5510204081632653, 0.5714285714285714, 0.5918367346938775, 0.6122448979591836, 0.6326530612244897, 0.6530612244897959, 0.673469387755102, 0.6938775510204082, 0.7142857142857142, 0.7346938775510203, 0.7551020408163265, 0.7755102040816326, 0.7959183673469387, 0.8163265306122448, 0.836734693877551, 0.8571428571428571, 0.8775510204081632, 0.8979591836734693, 0.9183673469387754, 0.9387755102040816, 0.9591836734693877, 0.9795918367346939, 1.0]
        #x_values = [0.0, 0.018400065266506227, 0.04225083661157623, 0.0713435070518306, 0.1036164118664373, 0.13698065932882636, 0.16947761902506267, 0.19850810716711273, 0.22314347791028566, 0.24053558565655847, 0.24886598803245524, 0.2502470372487695, 0.24941257744761433, 0.24756615214432298, 0.24481686563607896, 0.24130290560673967, 0.23698965095246982, 0.23242285078249267, 0.22531163517427788, 0.2110134548882222, 0.18623119147117725, 0.1479307251853749, 0.09847131457569316, 0.04111540547132665, -0.02233952894219675, -0.08884150619038655, -0.15891966620096387, -0.2407441175807782, -0.3366430472730907, -0.44693576549987823, -0.5680658106768092, -0.6975208703059096, -0.8321262196998409, -0.9699653368698024, -1.1090930486685822, -1.255144506570033, -1.4103667735456449, -1.5733007007462756, -1.7434963771088456, -1.9194542609028804, -2.1000907378795275, -2.285501961499942, -2.4756894577736315, -2.6734165188032692, -2.8782701025304545, -3.090085737186208, -3.308459127246535, -3.533712868740941, -3.7641269864926348, -4.0]
        #y_grid   = [0.0, 1.0]
        #y_values = [0.0, 0.0]
        #z_grid   = [0.0, 0.02040816326530612, 0.04081632653061224, 0.061224489795918366, 0.08163265306122448, 0.1020408163265306, 0.12244897959183673, 0.14285714285714285, 0.16326530612244897, 0.18367346938775508, 0.2040816326530612, 0.22448979591836732, 0.24489795918367346, 0.26530612244897955, 0.2857142857142857, 0.3061224489795918, 0.32653061224489793, 0.3469387755102041, 0.36734693877551017, 0.3877551020408163, 0.4081632653061224, 0.42857142857142855, 0.44897959183673464, 0.4693877551020408, 0.4897959183673469, 0.5102040816326531, 0.5306122448979591, 0.5510204081632653, 0.5714285714285714, 0.5918367346938775, 0.6122448979591836, 0.6326530612244897, 0.6530612244897959, 0.673469387755102, 0.6938775510204082, 0.7142857142857142, 0.7346938775510203, 0.7551020408163265, 0.7755102040816326, 0.7959183673469387, 0.8163265306122448, 0.836734693877551, 0.8571428571428571, 0.8775510204081632, 0.8979591836734693, 0.9183673469387754, 0.9387755102040816, 0.9591836734693877, 0.9795918367346939, 1.0]
        #z_values = [0.0, 2.387755102040816, 4.775510204081632, 7.163265306122448, 9.551020408163264, 11.938775510204081, 14.326530612244898, 16.714285714285715, 19.10204081632653, 21.489795918367346, 23.877551020408163, 26.265306122448976, 28.653061224489797, 31.04081632653061, 33.42857142857143, 35.81632653061224, 38.20408163265306, 40.59183673469388, 42.979591836734684, 45.36734693877551, 47.75510204081632, 50.14285714285714, 52.53061224489795, 54.91836734693877, 57.30612244897959, 59.69387755102041, 62.08163265306122, 64.46938775510203, 66.85714285714285, 69.24489795918367, 71.63265306122447, 74.0204081632653, 76.40816326530611, 78.79591836734693, 81.18367346938776, 83.57142857142857, 85.95918367346938, 88.3469387755102, 90.73469387755102, 93.12244897959182, 95.51020408163265, 97.89795918367345, 100.28571428571428, 102.6734693877551, 105.0612244897959, 107.44897959183673, 109.83673469387753, 112.22448979591836, 114.61224489795919, 117.0]
        #r_ = [0.0, 0.02, 0.15, 0.245170, 1.0]
        #ac = [0.5, 0.5, 0.316, 0.25, 0.25]
        #r0 = r/r[-1]
        #z = np.interp(r0, z_grid, z_values)
        #x = np.interp(r0, x_grid, x_values)
        #y = np.interp(r0, y_grid, y_values)
        #xp = np.interp(r0, grid, values)
        #df['z'] = z
        #df['x'] = x
        #df['y'] = y
        #df['xp'] = xp
        #ACloc = np.interp(r0, r_,ac)

        ## Get the absolute offset between pitch axis (rotation center) and aerodynamic center
        #ch_offset = inputs['chord'] * (inputs['ac'] - inputs['le_location'])
        ## Rotate it by the twist using the AD15 coordinate system
        #x , y = util.rotate(0., 0., 0., ch_offset, -np.deg2rad(inputs['theta']))
        ## Apply offset to determine the AC axis
        #BlCrvAC = inputs['ref_axis_blade'][:,0] + x
        #BlSwpAC = inputs['ref_axis_blade'][:,1] + y

        # --- Adding C2 axis
        ACloc = r*0 + 0.25 # distance (in chord) from leading edge to aero center
        n=int(len(r)*0.15) # 15% span
        ACloc[:n]=np.linspace(0.5,0.25, n) # Root is at 0

        dx = chord*(0.5-ACloc) * np.sin(twist)  # Should be mostly >0
        dy = chord*(0.5-ACloc) * np.cos(twist)  # Should be mostly >0
        prebend = prebendAC + dx
        sweep   = sweepAC   + dy
        df['c2_Crv_Approx_[m]'] = prebend
        df['c2_Swp_Approx_[m]'] = sweep
        df['AC_Approx_[-]'] = ACloc
        # --- Calc CvrAng
        dr = np.gradient(aeroNodes[:,0])
        dx = np.gradient(aeroNodes[:,1])
        df['CrvAng_Calc_[-]'] =  np.degrees(np.arctan2(dx,dr))
        return df

    @property
    def _IComment(self): return [1]


# --------------------------------------------------------------------------------}
# --- AeroDyn Polar 
# --------------------------------------------------------------------------------{
class ADPolarFile(FASTInputFileBase):
    @staticmethod
    def formatName():
        return 'FAST AeroDyn polar file'

    @classmethod
    def from_fast_input_file(cls, parent):
        self = cls()
        self.setData(filename=parent.filename, data=parent.data, hasNodal=parent.hasNodal, module='ADPolar')
        return self

    def __init__(self, filename=None, hasUA=True, numTabs=1, **kwargs):
        FASTInputFileBase.__init__(self, filename, **kwargs)
        if filename is None:
            # Define a prototype for this file format
            self.addComment('! ------------ AirfoilInfo Input File ------------------------------------------')
            self.addComment('! Airfoil definition, written by ADPolarFile')
            self.addComment('! ')
            self.addComment('! ')
            self.addComment('! ------------------------------------------------------------------------------')
            self.addValKey("DEFAULT", 'InterpOrd' , 'Interpolation order to use for quasi-steady table lookup {1=linear; 3=cubic spline; "default"} [default=3]')
            self.addValKey(        1, 'NonDimArea', 'The non-dimensional area of the airfoil (area/chord^2) (set to 1.0 if unsure or unneeded)')
            self.addValKey(        0, 'NumCoords' , 'The number of coordinates in the airfoil shape file.  Set to zero if coordinates not included.')
            self.addValKey( numTabs , 'NumTabs'   , 'Number of airfoil tables in this file.  Each table must have lines for Re and Ctrl.')
            # TODO multiple tables
            for iTab in range(numTabs):
                if numTabs==1:
                    labOffset =''
                else:
                    labOffset ='_'+str(iTab+1)
                self.addComment('! ------------------------------------------------------------------------------')
                self.addComment('! data for table {}'.format(iTab+1))
                self.addComment('! ------------------------------------------------------------------------------')
                self.addValKey(     1.0 ,    'Re'  +labOffset , 'Reynolds number in millions')
                self.addValKey(       0 ,    'Ctrl'+labOffset , 'Control setting')
                if hasUA:
                    self.addValKey(True  ,    'InclUAdata', 'Is unsteady aerodynamics data included in this table? If TRUE, then include 30 UA coefficients below this line')
                    self.addComment('!........................................')
                    self.addValKey(     np.nan   , 'alpha0'    + labOffset, r"0-lift angle of attack, depends on airfoil.")
                    self.addValKey(     np.nan   , 'alpha1'    + labOffset, r"Angle of attack at f=0.7, (approximately the stall angle) for AOA>alpha0. (deg)")
                    self.addValKey(     np.nan   , 'alpha2'    + labOffset, r"Angle of attack at f=0.7, (approximately the stall angle) for AOA<alpha0. (deg)")
                    self.addValKey(          1   , 'eta_e'     + labOffset, r"Recovery factor in the range [0.85 - 0.95] used only for UAMOD=1, it is set to 1 in the code when flookup=True. (-)")
                    self.addValKey(     np.nan   , 'C_nalpha'  + labOffset, r"Slope of the 2D normal force coefficient curve. (1/rad)")
                    self.addValKey(   "DEFAULT"  , 'T_f0'      + labOffset, r"Initial value of the time constant associated with Df in the expression of Df and f''. [default = 3]")
                    self.addValKey(   "DEFAULT"  , 'T_V0'      + labOffset, r"Initial value of the time constant associated with the vortex lift decay process; it is used in the expression of Cvn. It depends on Re,M, and airfoil class. [default = 6]")
                    self.addValKey(   "DEFAULT"  , 'T_p'       + labOffset, r"Boundary-layer,leading edge pressure gradient time constant in the expression of Dp. It should be tuned based on airfoil experimental data. [default = 1.7]")
                    self.addValKey(   "DEFAULT"  , 'T_VL'      + labOffset, r"Initial value of the time constant associated with the vortex advection process; it represents the non-dimensional time in semi-chords, needed for a vortex to travel from LE to trailing edge (TE); it is used in the expression of Cvn. It depends on Re, M (weakly), and airfoil. [valid range = 6 - 13, default = 11]")
                    self.addValKey(   "DEFAULT"  , 'b1'        + labOffset, r"Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils, but may be different for turbine airfoils. [from experimental results, defaults to 0.14]")
                    self.addValKey(   "DEFAULT"  , 'b2'        + labOffset, r"Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils, but may be different for turbine airfoils. [from experimental results, defaults to 0.53]")
                    self.addValKey(   "DEFAULT"  , 'b5'        + labOffset, r"Constant in the expression of K'''_q,Cm_q^nc, and k_m,q.  [from  experimental results, defaults to 5]")
                    self.addValKey(   "DEFAULT"  , 'A1'        + labOffset, r"Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils, but may be different for turbine airfoils. [from experimental results, defaults to 0.3]")
                    self.addValKey(   "DEFAULT"  , 'A2'        + labOffset, r"Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils, but may be different for turbine airfoils. [from experimental results, defaults to 0.7]")
                    self.addValKey(   "DEFAULT"  , 'A5'        + labOffset, r"Constant in the expression of K'''_q,Cm_q^nc, and k_m,q. [from experimental results, defaults to 1]")
                    self.addValKey(          0   , 'S1'        + labOffset, r"Constant in the f curve best-fit for alpha0<=AOA<=alpha1; by definition it depends on the airfoil. [ignored if UAMod<>1]")
                    self.addValKey(          0   , 'S2'        + labOffset, r"Constant in the f curve best-fit for         AOA> alpha1; by definition it depends on the airfoil. [ignored if UAMod<>1]")
                    self.addValKey(          0   , 'S3'        + labOffset, r"Constant in the f curve best-fit for alpha2<=AOA< alpha0; by definition it depends on the airfoil. [ignored if UAMod<>1]")
                    self.addValKey(          0   , 'S4'        + labOffset, r"Constant in the f curve best-fit for         AOA< alpha2; by definition it depends on the airfoil. [ignored if UAMod<>1]")
                    self.addValKey(     np.nan   , 'Cn1'       + labOffset, r"Critical value of C0n at leading edge separation. It should be extracted from airfoil data at a given Mach and Reynolds number. It can be calculated from the static value of Cn at either the break in the pitching moment or the loss of chord force at the onset of stall. It is close to the condition of maximum lift of the airfoil at low Mach numbers.")
                    self.addValKey(     np.nan   , 'Cn2'       + labOffset, r"As Cn1 for negative AOAs.")
                    self.addValKey(   "DEFAULT"  , 'St_sh'     + labOffset, r"Strouhal's shedding frequency constant.  [default = 0.19]")
                    self.addValKey(     np.nan   , 'Cd0'       + labOffset, r"2D drag coefficient value at 0-lift.")
                    self.addValKey(     np.nan   , 'Cm0'       + labOffset, r"2D pitching moment coefficient about 1/4-chord location, at 0-lift, positive if nose up. [If the aerodynamics coefficients table does not include a column for Cm, this needs to be set to 0.0]")
                    self.addValKey(          0   , 'k0'        + labOffset, r"Constant in the \hat(x)_cp curve best-fit; = (\hat(x)_AC-0.25).  [ignored if UAMod<>1]")
                    self.addValKey(          0   , 'k1'        + labOffset, r"Constant in the \hat(x)_cp curve best-fit.  [ignored if UAMod<>1]")
                    self.addValKey(          0   , 'k2'        + labOffset, r"Constant in the \hat(x)_cp curve best-fit.  [ignored if UAMod<>1]")
                    self.addValKey(          0   , 'k3'        + labOffset, r"Constant in the \hat(x)_cp curve best-fit.  [ignored if UAMod<>1]")
                    self.addValKey(          0   , 'k1_hat'    + labOffset, r"Constant in the expression of Cc due to leading edge vortex effects.  [ignored if UAMod<>1]")
                    self.addValKey(   "DEFAULT"  , 'x_cp_bar'  + labOffset, r"Constant in the expression of \hat(x)_cp^v. [ignored if UAMod<>1, default = 0.2]")
                    self.addValKey(   "DEFAULT"  , 'UACutout'  + labOffset, r"Angle of attack above which unsteady aerodynamics are disabled (deg). [Specifying the string 'Default' sets UACutout to 45 degrees]")
                    self.addValKey(   "DEFAULT"  , 'filtCutOff'+ labOffset, r"Reduced frequency cut-off for low-pass filtering the AoA input to UA, as well as the 1st and 2nd derivatives (-) [default = 0.5]")
                    self.addComment('!........................................')
                else:             
                    self.addValKey(False ,    'InclUAdata'+labOffset, 'Is unsteady aerodynamics data included in this table? If TRUE, then include 30 UA coefficients below this line')
                self.addComment('! Table of aerodynamics coefficients')
                self.addValKey(0 ,    'NumAlf'+labOffset, '! Number of data lines in the following table')
                self.addTable('AFCoeff'+labOffset, np.zeros((0,4)), tabType=2, tabDimVar='NumAlf', cols=['Alpha', 'Cl', 'Cd', 'Cm'], units=['(deg)', '(-)', '(-)', '(-)'])
        self.module='ADPolar'

    def _writeSanityChecks(self):
        """ Sanity checks before write"""
        nTabs = self['NumTabs']
        if nTabs==1:
            self['NumAlf']=self['AFCoeff'].shape[0]
        else:
            for iTab in range(nTabs):
                labOffset='_{}'.format(iTab+1)
                self['NumAlf'+labOffset] = self['AFCoeff'+labOffset].shape[0]
        # Potentially compute unsteady params here

    def _write(self):
        nTabs = self['NumTabs']
        if nTabs==1:
            FASTInputFileBase._write(self)
        else:
            self._writeSanityChecks()
            Labs=['Re','Ctrl','UserProp','alpha0','alpha1','alpha2','eta_e','C_nalpha','T_f0','T_V0','T_p','T_VL','b1','b2','b5','A1','A2','A5','S1','S2','S3','S4','Cn1','Cn2','St_sh','Cd0','Cm0','k0','k1','k2','k3','k1_hat','x_cp_bar','UACutout','filtCutOff','InclUAdata','NumAlf','AFCoeff']
            # Store all labels
            AllLabels=[self.data[i]['label'] for i in range(len(self.data))]
            # Removing lab Offset - TODO TEMPORARY HACK
            for iTab in range(nTabs):
                labOffset='_{}'.format(iTab+1)
                for labRaw in Labs:
                    i = self.getIDSafe(labRaw+labOffset)
                    if i>0:
                        self.data[i]['label'] = labRaw
            # Write
            with open(self.filename,'w') as f:
                f.write(self.toString())
            # Restore labels 
            for i,labFull in enumerate(AllLabels):
                self.data[i]['label'] = labFull

    def _toDataFrame(self):
        dfs = FASTInputFileBase._toDataFrame(self)
        if not isinstance(dfs, dict):
            dfs={'AFCoeff':dfs}

        for k,df in dfs.items():
            sp = k.split('_')
            if len(sp)==2:
                labOffset='_'+sp[1]
            else:
                labOffset=''
            alpha = df['Alpha_[deg]'].values*np.pi/180.
            Cl    = df['Cl_[-]'].values
            Cd    = df['Cd_[-]'].values

            # Cn with Cd0
            try:
                Cd0   = self['Cd0'+labOffset]
                # Cn (with or without Cd0)
                Cn1 = Cl*np.cos(alpha)+ (Cd-Cd0)*np.sin(alpha) 
                df['Cn_Cd0off_[-]'] = Cn1
            except:
                pass

            # Regular Cn
            Cn  = Cl*np.cos(alpha)+ Cd*np.sin(alpha) 
            df['Cn_[-]'] = Cn

            # Linear Cn
            try:
                CnLin_ = self['C_nalpha'+labOffset]*(alpha-self['alpha0'+labOffset]*np.pi/180.)
                CnLin = CnLin_.copy()
                CnLin[alpha<-20*np.pi/180]=np.nan
                CnLin[alpha> 30*np.pi/180]=np.nan
                df['Cn_pot_[-]'] = CnLin
            except:
                pass

            # Highlighting points surrounding 0 1 2 Cn points
            CnPoints = Cn*np.nan
            try:
                iBef2 = np.where(alpha<self['alpha2'+labOffset]*np.pi/180.)[0][-1]
                iBef1 = np.where(alpha<self['alpha1'+labOffset]*np.pi/180.)[0][-1]
                iBef0 = np.where(alpha<self['alpha0'+labOffset]*np.pi/180.)[0][-1]
                CnPoints[iBef2:iBef2+2] = Cn[iBef2:iBef2+2]
                CnPoints[iBef1:iBef1+2] = Cn[iBef1:iBef1+2]
                CnPoints[iBef0:iBef0+2] = Cn[iBef0:iBef0+2]
                df['Cn_012_[-]'] = CnPoints
            except:
                pass

            # Cnf
            try:
                df['Cn_f_[-]'] = CnLin_ * ((1 + np.sqrt(0.7)) / 2) ** 2
            except:
                pass

        if len(dfs)==1:
            dfs=dfs[list(dfs.keys())[0]]
        return dfs

    @property
    def comment(self):
        return '\n'.join([self.data[i]['value'] for i in self._IComment])

    @comment.setter
    def comment(self, comment):
        splits = comment.split('\n')
        for i,com in zip(self._IComment, splits):
            self.data[i]['value'] = '! ' +com

    @property
    def _IComment(self):
        """ return indices of comment line"""
        I=[]
        for i in [1,2,3]:
            if self.data[i]['value'].startswith('!'):
                if not self.data[i]['value'].startswith('! ---'):
                    I.append(i)
        return I


# --------------------------------------------------------------------------------}
# --- ExtPtfm 
# --------------------------------------------------------------------------------{
def detectAndReadExtPtfmSE(self, lines):
    # TODO place this in ExtPtfmFile
    def readmat(n,m,lines,iStart):
        M=np.zeros((n,m))
        for j in np.arange(n):
            i=iStart+j
            M[j,:]=np.array(lines[i].split()).astype(float)
        return M
    if len(lines)<10:
        return False
    if not (lines[0][0]=='!' and lines[1][0]=='!'):
        return False
    if lines[1].lower().find('flex')<0:
        return
    if  lines[2].lower().find('!dimension')<0:
        return
    # --- At this stage we assume it's in the proper format
    self.module='ExtPtfm'
    nDOFCommon = -1
    i=2;
    try:
        while i<len(lines):
            l=lines[i].lower()
            if l.find('!mass')==0:
                l=lines[i+1]
                nDOF=int(l.split(':')[1])
                if nDOF<-1 or nDOF!=nDOFCommon:
                    raise BrokenFormatError('ExtPtfm stiffness matrix nDOF issue. nDOF common: {}, nDOF provided: {}'.format(nDOFCommon,nDOF))
                self.addKeyVal('MassMatrix',readmat(nDOF,nDOF,lines,i+2))
                i=i+1+nDOF
            elif l.find('!stiffness')==0:
                l=lines[i+1]
                nDOF=int(l.split(':')[1])
                if nDOF<-1 or nDOF!=nDOFCommon:
                    raise BrokenFormatError('ExtPtfm stiffness matrix nDOF issue nDOF common: {}, nDOF provided: {}'.format(nDOFCommon,nDOF))
                self.addKeyVal('StiffnessMatrix',readmat(nDOF,nDOF,lines,i+2))
                i=i+1+nDOF
            elif l.find('!damping')==0:
                l=lines[i+1]
                nDOF=int(l.split(':')[1])
                if nDOF<-1 or nDOF!=nDOFCommon:
                    raise BrokenFormatError('ExtPtfm damping matrix nDOF issue nDOF common: {}, nDOF provided: {}'.format(nDOFCommon,nDOF))
                self.addKeyVal('DampingMatrix',readmat(nDOF,nDOF,lines,i+2))
                i=i+1+nDOF
            elif l.find('!loading')==0:
                try: 
                    nt=int(self['T']/self['dt'])+1
                except:
                    raise BrokenFormatError('Cannot read loading since time step and simulation time not properly set.')
                self.addKeyVal('Loading',readmat(nt,nDOFCommon+1,lines,i+2))
                i=i+nt+1
            elif len(l)>0:
                if l[0]=='!':
                    if l.find('!dimension')==0:
                        self.addKeyVal('nDOF',int(l.split(':')[1]))
                        nDOFCommon=self['nDOF']
                    elif l.find('!time increment')==0:
                        self.addKeyVal('dt',float(l.split(':')[1]))
                    elif l.find('!total simulation time')==0:
                        self.addKeyVal('T',float(l.split(':')[1]))
                elif len(l.strip())==0:
                    pass
                else:
                    raise BrokenFormatError('Unexcepted content found on line {}'.format(i))
            i+=1
    except BrokenFormatError as e:
        raise e
    except: 
        raise

    return True

class ExtPtfmFile(FASTInputFileBase):
    @classmethod
    def from_fast_input_file(cls, parent):
        self = cls()
        self.setData(filename=parent.filename, data=parent.data, hasNodal=parent.hasNodal, module='ExtPtfm')
        return self

    def __init__(self, filename=None, **kwargs):
        FASTInputFileBase.__init__(self, filename, **kwargs)
        if filename is None:
            # Define a prototype for this file format
            self.addValKey(0 ,    'nDOF', '')
            self.addValKey(1 ,    'dt'  , '')
            self.addValKey(0 ,    'T'   , '')
            self.addTable('MassMatrix'     , np.zeros((0,0)), tabType=0) 
            self.addTable('StiffnessMatrix', np.zeros((0,0)), tabType=0) 
            self.addTable('DampingMatrix'  , np.zeros((0,0)), tabType=0) 
            self.addTable('Loading'        , np.zeros((0,0)), tabType=0) 
            self.comment=''
        self.module='ExtPtfm'


    def _read(self):
        with open(self.filename, 'r', errors="surrogateescape") as f:
            lines=f.read().splitlines()
        detectAndReadExtPtfmSE(self, lines)

    def toString(self):
        s=''
        s+='!Comment\n'
        s+='!Comment Flex 5 Format\n'
        s+='!Dimension: {}\n'.format(self['nDOF'])
        s+='!Time increment in simulation: {}\n'.format(self['dt'])
        s+='!Total simulation time in file: {}\n'.format(self['T'])

        s+='\n!Mass Matrix\n'
        s+='!Dimension: {}\n'.format(self['nDOF'])
        s+='\n'.join(''.join('{:16.8e}'.format(x) for x in y) for y in self['MassMatrix'])

        s+='\n\n!Stiffness Matrix\n'
        s+='!Dimension: {}\n'.format(self['nDOF'])
        s+='\n'.join(''.join('{:16.8e}'.format(x) for x in y) for y in self['StiffnessMatrix'])

        s+='\n\n!Damping Matrix\n'
        s+='!Dimension: {}\n'.format(self['nDOF'])
        s+='\n'.join(''.join('{:16.8e}'.format(x) for x in y) for y in self['DampingMatrix'])

        s+='\n\n!Loading and Wave Elevation\n'
        s+='!Dimension: 1 time column -  {} force columns\n'.format(self['nDOF'])
        s+='\n'.join(''.join('{:16.8e}'.format(x) for x in y) for y in self['Loading'])
        return s

    def _writeSanityChecks(self):
        """ Sanity checks before write"""
        assert self['MassMatrix'].shape[0]      == self['nDOF']
        assert self['StiffnessMatrix'].shape[0] == self['nDOF']
        assert self['DampingMatrix'].shape[0]   == self['nDOF']
        assert self['MassMatrix'].shape[0]      == self['MassMatrix'].shape[1]
        assert self['StiffnessMatrix'].shape[0] == self['StiffnessMatrix'].shape[1]
        assert self['DampingMatrix'].shape[0]   == self['DampingMatrix'].shape[1]
        #         if self['T']>0:
        #             assert self['Loading'].shape[0]     == (int(self['T']/self['dT'])+1

    def _toDataFrame(self):
        # Special types, TODO Subclass
        nDOF=self['nDOF']
        Cols=['Time_[s]','InpF_Fx_[N]', 'InpF_Fy_[N]', 'InpF_Fz_[N]', 'InpF_Mx_[Nm]', 'InpF_My_[Nm]', 'InpF_Mz_[Nm]']
        Cols+=['CBF_{:03d}_[-]'.format(iDOF+1) for iDOF in np.arange(nDOF)]
        Cols=Cols[:nDOF+1]
        #dfs['Loading']         = pd.DataFrame(data = self['Loading'],columns  = Cols)
        dfs = pd.DataFrame(data = self['Loading'],columns  = Cols)

        #Cols=['SurgeAcc_[m/s]', 'SwayAcc_[m/s]', 'HeaveAcc_[m/s]', 'RollAcc_[rad/s]', 'PitchAcc_[rad/s]', 'YawAcc_[rad/s]']
        #Cols+=['CBQD_{:03d}_[-]'.format(iDOF+1) for iDOF in np.arange(nDOF)]
        #Cols=Cols[:nDOF]
        #dfs['MassMatrix']      = pd.DataFrame(data = self['MassMatrix'], columns=Cols)

        #Cols=['SurgeVel_[m/s]', 'SwayVel_[m/s]', 'HeaveVel_[m/s]', 'RollVel_[rad/s]', 'PitchVel_[rad/s]', 'YawVel_[rad/s]']
        #Cols+=['CBQD_{:03d}_[-]'.format(iDOF+1) for iDOF in np.arange(nDOF)]
        #Cols=Cols[:nDOF]
        #dfs['DampingMatrix']   = pd.DataFrame(data = self['DampingMatrix'], columns=Cols)

        #Cols=['Surge_[m]', 'Sway_[m]', 'Heave_[m]', 'Roll_[rad]', 'Pitch_[rad]', 'Yaw_[rad]']
        #Cols+=['CBQ_{:03d}_[-]'.format(iDOF+1) for iDOF in np.arange(nDOF)]
        #Cols=Cols[:nDOF]
        #dfs['StiffnessMatrix'] = pd.DataFrame(data = self['StiffnessMatrix'], columns=Cols)
        return dfs



if __name__ == "__main__":
    f = FASTInputFile('tests/example_files/FASTIn_HD_SeaState.dat')
    print(f)
    pass
    #B=FASTIn('Turbine.outb')



