""" 
Tools to manipulate:
 - "CDDOP":  campbell data  at one operating point 
 - "CampbellData" list of CDDOP (for multiple operating points)

"""

import re
import numpy as np

def campbell_diagram_data_oneOP(mbc_data, BladeLen=None, TowerLen=None):
    """ 
    Return Campbell Data for one Operating Point "CDDOP"
    
    """
    CDDOP={}
    usePercent = False;
    #
    # mbc_data.eigSol = eiganalysis(mbc_data.AvgA);
    ndof = mbc_data['ndof2'] + mbc_data['ndof1']; #size(mbc_data.AvgA,1)/2;          # number of translational states
    nModes = len(mbc_data['eigSol']['Evals'])
    #print(nModes)
    DescStates = PrettyStateDescriptions(mbc_data['DescStates'], mbc_data['ndof2'], mbc_data['performedTransformation']);

    ## store indices of max mode for state and to order natural frequencies
    #StatesMaxMode_vals = np.amax(mbc_data['eigSol']['MagnitudeModes'],axis=1); # find which mode has the maximum value for each state (max of each row before scaling)
    StatesMaxMode = np.argmax(mbc_data['eigSol']['MagnitudeModes'],axis=1); # find which mode has the maximum value for each state (max of each row before scaling)
    SortedFreqIndx = np.argsort((mbc_data['eigSol']['NaturalFreqs_Hz']).flatten(),kind="heapsort");
    #print(SortedFreqIndx)


    if BladeLen is not None and TowerLen is not None:
        ## get the scaling factors for the mode rows
        ScalingFactor = getScaleFactors(DescStates, TowerLen, BladeLen);

        ## scale the magnitude of the modes by ScalingFactor (for consistent units)
        #  and then scale the columns so that their maximum is 1

        ModesMagnitude = np.matmul(np.diag(ScalingFactor), mbc_data['eigSol']['MagnitudeModes']); # scale the rows
        #print(ModesMagnitude)
        
        CDDOP['ScalingFactor'] = ScalingFactor;
    else:
        ModesMagnitude = mbc_data['eigSol']['MagnitudeModes'];

    if usePercent:
        scaleCol = np.sum( ModesMagnitude )/100; # find the sum of the column, and multiply by 100 (divide here) to get a percentage
    else:
        scaleCol = np.amax(ModesMagnitude,axis=0); #find the maximum value in the column, so the first element has value of 1

    ModesMagnitude = np.matmul(ModesMagnitude,np.diag(1./scaleCol)) # scale the columns

    # --- Summary data (array for all modes)
    CDDOP['NaturalFreq_Hz'] = mbc_data['eigSol']['NaturalFreqs_Hz'][SortedFreqIndx]
    CDDOP['DampingRatio']   = mbc_data['eigSol']['DampRatios'][SortedFreqIndx]
    CDDOP['RotSpeed_rpm']   = mbc_data['RotSpeed_rpm']
    if 'WindSpeed' in mbc_data:
        CDDOP['WindSpeed']  = mbc_data['WindSpeed']
        
    #print(ModesMagnitude)

    # --- Storing data per mode
    CDDOP['Modes']=[]
    for i in range(nModes):
        CData={}
        CData['NaturalFreq_Hz'] = mbc_data['eigSol']['NaturalFreqs_Hz'][SortedFreqIndx[i]][0]
        CData['DampedFreq_Hz']  = mbc_data['eigSol']['DampedFreqs_Hz'][SortedFreqIndx[i]][0];
        CData['DampingRatio']   = mbc_data['eigSol']['DampRatios'][SortedFreqIndx[i]][0];

        
        #print(np.argsort(ModesMagnitude[:,SortedFreqIndx[0]])[::-1])
        # sort indices in descending order
        sort_state = np.argsort( ModesMagnitude[:,SortedFreqIndx[i]])[::-1];

        #print(type(sort_state))
        CData['DescStates']=[DescStates[i] for i in sort_state]
        CData['MagnitudePhase']=ModesMagnitude[sort_state,SortedFreqIndx[i]];
        Phase =                mbc_data['eigSol']['PhaseModes_deg'][sort_state,SortedFreqIndx[i]];
        # if the phase is more than +/- 90 degrees different than the first
        # one (whose value == 1 or is the largest %), we'll stick a negative value on the magnitude:
        Phase = np.mod(Phase, 360);
    
        CData['PhaseDiff'] = np.mod( Phase - Phase[0], 360); # difference in range [0, 360)
        PhaseIndx = CData['PhaseDiff'] > 180;
        CData['PhaseDiff'][PhaseIndx] = CData['PhaseDiff'][PhaseIndx] - 360;   # move to range (-180, 180]
    
        if ~usePercent:
            PhaseIndx = CData['PhaseDiff'] > 90;
            CData['MagnitudePhase'][PhaseIndx] = -1*CData['MagnitudePhase'][PhaseIndx];
            CData['PhaseDiff'][PhaseIndx] = CData['PhaseDiff'][PhaseIndx] - 180;

            PhaseIndx = CData['PhaseDiff'] <= -90;
            CData['MagnitudePhase'][PhaseIndx] = -1*CData['MagnitudePhase'][PhaseIndx];
            CData['PhaseDiff'][PhaseIndx] = CData['PhaseDiff'][PhaseIndx] + 180;

        #print(CData['MagnitudePhase'])
        #print(CData['PhaseDiff'])

        CData['StateHasMaxAtThisMode'] = np.ones(ndof, dtype=bool);
        ix = (StatesMaxMode == SortedFreqIndx[i]);
        tmp=ix[sort_state]
        CData['StateHasMaxAtThisMode']=tmp
                    
        #print(CData['StateHasMaxAtThisMode'])
        #print(CData['NaturalFreq_Hz'])
        CDDOP['Modes'].append(CData)

    #print(CDDOP[0]['MagnitudePhase'])
    # Adding short description to summary
    CDDOP['ShortModeDescr'] =[extractShortModeDescription(CDDOP['Modes'][i]) for i in range(nModes)]


    CDDOP['nColsPerMode'] = 5;
    CDDOP['ModesTable'] = {}

    for i in range(nModes):
        colStart = i*CDDOP['nColsPerMode'];
        CDDOP['ModesTable'][1, colStart+1 ] = 'Mode number:';
        CDDOP['ModesTable'][1, colStart+2 ] = i;

        CDDOP['ModesTable'][2, colStart+1 ] = 'Natural (undamped) frequency (Hz):';
        CDDOP['ModesTable'][2, colStart+2 ] = CDDOP['Modes'][i]['NaturalFreq_Hz']

        CDDOP['ModesTable'][3, colStart+1 ] = 'Damped frequency (Hz):';
        CDDOP['ModesTable'][3, colStart+2 ] = CDDOP['Modes'][i]['DampedFreq_Hz']

        CDDOP['ModesTable'][4, colStart+1 ] = 'Damping ratio (-):';
        CDDOP['ModesTable'][4, colStart+2 ] = CDDOP['Modes'][i]['DampingRatio']
        
        CDDOP['ModesTable'][5, colStart+1 ] = 'Mode ' + str(i) + ' state description';
        CDDOP['ModesTable'][5, colStart+2 ] = 'State has max at mode ' + str(i);
        if usePercent:
            CDDOP['ModesTable'][5, colStart+3 ] = 'Mode ' + str(i) + ' contribution (%)';
        else:
            CDDOP['ModesTable'][5, colStart+3 ] = 'Mode ' + str(i) + ' signed magnitude';

        CDDOP['ModesTable'][5, colStart+4 ] = 'Mode ' + str(i) + ' phase (deg)';

        # need to cross check these 4 lines
        CDDOP['ModesTable'][6,colStart+1] = CDDOP['Modes'][i]['DescStates'];
        CDDOP['ModesTable'][6,colStart+2] = CDDOP['Modes'][i]['StateHasMaxAtThisMode'];
        CDDOP['ModesTable'][6,colStart+3] = CDDOP['Modes'][i]['MagnitudePhase'];
        CDDOP['ModesTable'][6,colStart+4] = CDDOP['Modes'][i]['PhaseDiff'];

    #print(CDDOP['ModesTable'])
    return CDDOP


def getScaleFactors(DescStates, TowerLen, BladeLen):
    
    ScalingFactor = np.ones(len(DescStates))    
    
    # look at the state description strings for tower and blade
    # translational dofs:
    for i in range(len(ScalingFactor)):
        
        # look for blade dofs:
        if DescStates[i].find('blade')!=-1 or DescStates[i].find('Blade')!=-1:
           if DescStates[i].find('rotational')==-1: # make sure this isn't a rotational dof from BeamDyn
               ScalingFactor[i] = 1.0/BladeLen;
               #print(DescStates[i])
               
            # look for tower translational dofs:
        if DescStates[i].find('tower')!=-1 or DescStates[i].find('Tower')!=-1:
            ScalingFactor[i] = 1.0/TowerLen;

            # look for blade dofs:
        elif DescStates[i].find('blade')!=-1 or DescStates[i].find('Blade')!=-1:
            if DescStates[i].find('rotational')==-1: # make sure this isn't a rotational dof from BeamDyn
                ScalingFactor[i] =  1.0/BladeLen;

    return ScalingFactor


def printCampbellDataOP(CDDOP, nModesMax=15, nCharMaxDesc=50) :
    """ 
    Print frequencies, damping and mode description for Campbell diagram data at one operating points (CDDOP)
    INPUTS:
    - CDDOP: dictionary as returned by campbell_diagram_data_oneOP
    - nModesMax   : Maximum number of modes to be shown
    - nCharMaxDesc: Maximum number of characters for description written to screen
    OUTPUTS
     - Freq: frequencies for each modes
     - Damp: Damping for each modes

    See matlab-toolbox/Campbell/printCDDDOP.m
    """
    nModesMax = np.min([len(CDDOP['Modes']),nModesMax])
    Freq = np.array([CDDOP['Modes'][i]['NaturalFreq_Hz'] for i in np.arange(nModesMax)])
    Damp = np.array([CDDOP['Modes'][i]['DampingRatio']   for i in np.arange(nModesMax)])
    print('Mode, NatFreq_[Hz], Damp_Ratio_[-], LogDec._[%], Mode_content_[-]')
    for i in np.arange(nModesMax):
        Mode = CDDOP['Modes'][i]
        # Extracting description the best we can
        Desc = extractShortModeDescription(Mode)
        print('{:3d} ,{:12.3f}, {:8.5f}       , {:7.4f},  {:s}'.format(i+1,Mode['NaturalFreq_Hz'],Mode['DampingRatio'],Mode['DampingRatio']*100*2*np.pi, Desc[:min(nCharMaxDesc,len(Desc))]))
    return Freq, Damp


# --------------------------------------------------------------------------------}
# --- Manipulation of multiple operating points
# --------------------------------------------------------------------------------{

def campbellData2TXT(CD, txtFileName=None, nFreqOut=500, freqRange=None, posDampRange=None, skipNonEDBD=True):
    """ Write frequencies, damping, and mode contents for each operating points to a string
    Write to file if filename provided
    INPUTS:
     - nFreqOut: maximum number of frequencies to write to Campbell_Summary.txt file
     - freqRange:    range in which frequencies are "accepted",  if None: [-np.inf, np.inf]
     - posDampRange: range in which damping are "accepted'    ,  if None: [1e-5, 0.96]
    """
    if not isinstance(CD, list):
        CD=[CD]
    if freqRange is None:
        freqRange =[-np.inf, np.inf] 
    if posDampRange is None:
        posDampRange =[1e-5, 0.96] 

    def mode2txtline(cd, im):
        m    = cd['Modes'][im]
        Desc = cd['ShortModeDescr'][im]
        zeta = m['DampingRatio']
        return '{:03d} ; {:8.3f} ; {:7.4f} ; {:s}\n'.format(im+1,m['NaturalFreq_Hz'],m['DampingRatio'],Desc)

    txt=''
    for iOP,cd in enumerate(CD):
        WS  = cd['WindSpeed']
        RPM = cd['RotSpeed_rpm']
        nFreqOut_loc = np.min([len(cd['Modes']),nFreqOut])
        txt+='# -----------------------------------------------------------------------\n'
        txt+='# --- OP {:d} - WS {:.1f} - RPM {:.2f} \n'.format(iOP+1, WS, RPM)
        txt+='# -----------------------------------------------------------------------\n'
        txt+='# --- "Selected" modes\n'
        txt+='# ID; Freq [Hz]; Zeta [-]; Mode content\n'
        skippedDamp=[]
        skippedFreq=[]
        skippedEDBD=[]
        for im in np.arange(nFreqOut_loc):
            m    = cd['Modes'][im]
            Desc = cd['ShortModeDescr'][im]
            zeta = m['DampingRatio']
            freq = m['NaturalFreq_Hz']
            hasED = Desc.find('ED')>=0
            hasBD = Desc.find('BD')>=0
            hasAD = Desc.find('AD')>=0
            if (freq>freqRange[1] or freq<freqRange[0]):
                skippedFreq.append(im)
            elif (zeta>posDampRange[1] or abs(zeta)<posDampRange[0]):
                skippedDamp.append(im)
            elif skipNonEDBD and (not (hasBD or hasED)):
                skippedEDBD.append(im)
            else:
                txt+=mode2txtline(cd, im)
        if len(skippedEDBD)>0:
            txt+='# --- Skipped (No ED/BD)\n'
            for im in skippedEDBD:
                txt+=mode2txtline(cd, im)
        if len(skippedFreq)>0:
            txt+='# --- Skipped (Frequency outside of `freqRange`={})\n'.format(freqRange)
            for im in skippedFreq:
                txt+=mode2txtline(cd, im)
        if len(skippedDamp)>0:
            txt+='# --- Skipped (Damping outside of `posDampRange`={})\n'.format(posDampRange)
            for im in skippedDamp:
                txt+=mode2txtline(cd, im)
    if txtFileName is not None:
        with open(txtFileName, 'w') as f:
            f.write(txt)
    return txt

def campbellData2CSV(baseName, CD, modeID_table, modesDesc):
    # Write summary of Campbell data and modes identification to several CSV files
    # The files generated will be:
    #    - [BaseName '_ModesID.csv ']
    #    - [BaseName '_OP.csv      ']
    #    - [BaseName '_PointsXX.csv'] for each operating point where XX is the index.
    #
    # INPUTS:
    #   - BaseName    :  basename that will be used to create the different CSV Files
    #   - ModesData: structure as returned by IdentifyModes

    nModes = modeID_table.shape[0]
    nOP    = modeID_table.shape[1]
    filenames=[]

    # --- Write ModeID using Matlab format
    filename='{:s}_ModesID.csv'.format(baseName)
    filenames.append(filename)
    with open(filename, 'w') as f:
        f.write('Mode Number Table,' +','.join(['']*nOP) +'\n')
        if np.isnan(CD[0]['WindSpeed']):
            f.write('Rotor Speed (rpm),' +','.join([str(cd['RotSpeed_rpm']) for cd in CD]) +'\n')
        else:
            f.write('Wind Speed (mps),' +','.join([str(cd['WindSpeed']) for cd in CD]) +'\n')
        for im, v in enumerate(modesDesc):
            f.write(v[0]+',' +','.join([str(ID) for ID in modeID_table[im,:]]) +'\n')

    # --- Write OP using Matlab format
    filename='{:s}_OP.csv'.format(baseName)
    filenames.append(filename)
    with open(filename, 'w') as f:
        f.write('Operating Points,'  +','.join(['']*nOP) +'\n')
        f.write('Wind Speed (mps),'  +','.join([str(cd['WindSpeed']) for cd in CD]) +'\n')
        f.write('Rotor Speed (rpm),' +','.join([str(cd['RotSpeed_rpm']) for cd in CD]) +'\n')

    # --- Write Modes for each OP
    for iOP,cd in enumerate(CD):
        filename='{:s}_Point{:02d}.csv'.format(baseName,iOP+1)
        filenames.append(filename)
        with open(filename, 'w') as f:
            Modes = cd['Modes']
            f.write(','.join(['Mode Number:,{:d},,,'.format(im+1) for im,m in enumerate(Modes)]) +'\n')
            f.write(','.join(['Natural (undamped) frequency (Hz):, {:f},,,'.format(m['NaturalFreq_Hz']) for im,m in enumerate(Modes)]) +'\n')
            f.write(','.join(['Damped frequency (Hz):, {:f},,,'.format(m['DampedFreq_Hz']) for im,m in enumerate(Modes)]) +'\n')
            f.write(','.join(['Damping Ratio (-):, {:f},,,'.format(m['DampingRatio']) for im,m in enumerate(Modes)]) +'\n')
            f.write(','.join(['Mode {:d} state description,State has max at mode {:d},Mode {:d} signed magnitude,Mode {:d} phase (deg),'.format(im+1,im+1,im+1,im+1) for im,m in enumerate(Modes)]) +'\n')
            nComp = len(Modes[0]['DescStates'])
            for iC in np.arange(nComp):
                f.write(','.join(['{:s},{:d},{:f},{:f},'.format(m['DescStates'][iC].replace(',',' '),m['StateHasMaxAtThisMode'][iC],m['MagnitudePhase'][iC],m['PhaseDiff'][iC]) for im,m in enumerate(Modes)]) +'\n')
                   
    return filenames[0]


def IdentifyModes(CampbellData):
    """ 
    Attempts to perform an identification of the modes.
    For now, the method is based on the energy content of the modes, and the state descriptions where the energy is maximum

    Original contribution by: Srinivasa B. Ramisett, ramisettisrinivas@yahoo.com, http://ramisetti.github.io
    """
    #import pickle
    #pickle.dump(CampbellData, open('C:/Work/_libs/python-toolbox/data/_CampbellData_UA4_DB2.pkl','wb'))


    # --- Looking at states descriptions (of first run, first mode), to see if we are offshore. 
    # NOTE: DescStates is likely the same for all modes
    DescStates  = CampbellData[0]['Modes'][0]['DescStates']
    hasHeave    = any(['heave' in s.lower() for s in DescStates])
    hasSurge    = any(['surge' in s.lower() for s in DescStates])
    hasSway     = any(['sway' in s.lower() for s in DescStates])
    hasYaw      = any(['platform yaw' in s.lower() for s in DescStates])
    hasRoll     = any(['platform roll tilt' in s.lower() for s in DescStates])
    hasPitch    = any(['platform pitch tilt' in s.lower() for s in DescStates])
    hasEdge1Col = any(['1st edgewise bending-mode dof of blade collective' in s.lower() for s in DescStates])

    # --- Setting up a list of modes with 
    modesDesc = []
    modesDesc.append( ['Generator DOF (not shown)'     , 'ED Variable speed generator DOF, rad'] )
    # Platform DOFs
    if hasSurge:
        modesDesc.append( ['Platform surge', 'ED Platform horizontal surge translation DOF, m'] )
    if hasSway:
        modesDesc.append( ['Platform sway', 'ED Platform horizontal sway translation DOF, m'] )
    if hasHeave:
        modesDesc.append( ['Platform heave', 'ED Platform vertical heave translation DOF, m'] )
    if hasRoll:
        modesDesc.append( ['Platform roll', 'ED Platform roll tilt rotation DOF, rad'] )
    if hasPitch:
        modesDesc.append( ['Platform pitch', 'ED Platform pitch tilt rotation DOF, rad'] )
    if hasYaw:
        modesDesc.append( ['Platform yaw', 'ED Platform yaw rotation DOF, rad'] )

    modesDesc.append( ['1st Tower FA'                  , 'ED 1st tower fore-aft bending mode DOF, m'] )
    modesDesc.append( ['1st Tower SS'                  , 'ED 1st tower side-to-side bending mode DOF, m'] )
    modesDesc.append( ['1st Blade Flap (Regressive)'   , 'ED 1st flapwise bending-mode DOF of blade (sine|cosine), m', r'Blade (sine|cosine) finite element node \d rotational displacement in Y, rad'] )
    modesDesc.append( ['1st Blade Flap (Collective)'   , 'ED 1st flapwise bending-mode DOF of blade collective, m', r'Blade collective finite element node \d rotational displacement in Y, rad'] )
    modesDesc.append( ['1st Blade Flap (Progressive)'  , 'ED 1st flapwise bending-mode DOF of blade (sine|cosine), m']                                                                          )      # , ... # 'Blade (sine|cosine) finite element node \d rotational displacement in Y, rad']
    modesDesc.append( ['1st Blade Edge (Regressive)'   , 'ED 1st edgewise bending-mode DOF of blade (sine|cosine), m', r'Blade (sine|cosine) finite element node \d rotational displacement in X, rad'] )
    if hasEdge1Col:
        modesDesc.append(['1st Blade Edge (Collective)', 'ED 1st edgewise bending-mode DOF of blade collective, m']                                                                          )      # , ... # 'Blade (sine|cosine) finite element node \d rotational displacement in Y, rad']
    modesDesc.append( ['1st Blade Edge (Progressive)'  , 'ED 1st edgewise bending-mode DOF of blade (sine|cosine), m'] )
    modesDesc.append( ['1st Drivetrain Torsion'        , 'ED Drivetrain rotational-flexibility DOF, rad'] )
    modesDesc.append( ['2nd Tower FA'                  , 'ED 2nd tower fore-aft bending mode DOF, m'] )
    modesDesc.append( ['2nd Tower SS'                  , 'ED 2nd tower side-to-side bending mode DOF, m'] )
    modesDesc.append( ['2nd Blade Flap (Regressive)'   , 'ED 2nd flapwise bending-mode DOF of blade (sine|cosine), m'] )
    modesDesc.append( ['2nd Blade Flap (Collective)'   , 'ED 2nd flapwise bending-mode DOF of blade collective, m', r'Blade collective finite element node \d rotational displacement in Y, rad'] )
    modesDesc.append( ['2nd Blade Flap (Progressive)'  , 'ED 2nd flapwise bending-mode DOF of blade (sine|cosine), m'] )
    modesDesc.append( ['Nacelle Yaw (not shown)'  , 'ED Nacelle yaw DOF, rad'] )


    nModes = int(len(modesDesc))
    nRuns = int(len(CampbellData))
    modeID_table=np.zeros((nModes,nRuns)).astype(int)



    def doesDescriptionMatch(description, listOfModePatterns):
        """ loop through all mode desrption  """
        for iModePattern, modeIDdescList in enumerate(listOfModePatterns):
            modeIDName     = modeIDdescList[0]
            patternList    = modeIDdescList[1:] # list of patterns for a given mode
            found, pattern = doesDescriptionMatchPatternList(description, patternList)
            if found:
                return True, iModePattern, pattern
        return False, -1, ''

    def doesDescriptionMatchPatternList(description, patternList):
        """ loop through all patterns to find a match  """
        for pattern in patternList:
            # Looking for targetDesc into description
            if re.search(pattern ,description, re.IGNORECASE)!=None:
                return True, pattern
        return False, ''


    verbose=False
    Levels=[1,2,3]


    # --- Loop on operating points/Runs
    for i in range(nRuns):
        Modes  = CampbellData[i]['Modes']
        nModes = len(Modes)
        # Array of logical, False for Modes that are not identified
        modesIdentified = [False] * nModes
        modesSkipped    = [False] * nModes
        #verbose=verbose and i==0 # only display for first mode for now..
        #verbose=i==1

        # --- Give an index to each mode so that we can easily identify them
        for im, mode in enumerate(Modes):
            mode['index']=im

        # --- Skip modes based on simple criteria
        for im, mode in enumerate(Modes):
            if mode['NaturalFreq_Hz'] < 1e-5 or mode['DampingRatio'] > 0.98: 
                modesSkipped[im]=True



        modesDescLocal = modesDesc.copy()

        # --- Level 1 - Find well-defined modes (modes which have only one Max)
        if 1 in Levels:
            for im, mode in enumerate(Modes):
                if modesIdentified[im] or modesSkipped[im]:
                    continue # skip this mode since it has already been identified
                stateMax=np.argwhere((mode['StateHasMaxAtThisMode']==1)).flatten()
                if len(stateMax)==1:
                    description = mode['DescStates'][stateMax[0]]
                    if description.startswith('AD'):
                        # We skipp the pure "AD" modes
                        modesSkipped[im] = True
                        continue
                    found, modeID, patternMatched = doesDescriptionMatch(description, modesDescLocal)
                    if found and modeID_table[modeID,i]==0:
                        modesDescLocal[modeID] = [None,]   # we empty this mode patternlist so that it cannot be matched again
                        modesIdentified[im]    = True
                        modeID_table[modeID,i] = im+1 # Using Matlab Indexing
                        if verbose:
                            print('L1 Mode {} identified using pattern: {}'.format(im+1,patternMatched))
                            print('     String was: ', description)
                    #else:
                    #    print('>>> Cannot identify mode with description {}. Update MBC3 script.'.format(description))

        # --- Level 2 - Find modes with several max - Looping on mode pattern to respect expected frequency order
        if 2 in Levels:
            for modeID, modeIDdescList in enumerate(modesDescLocal):
                modeIDName     = modeIDdescList[0]
                patternList   = modeIDdescList[1:]
                # Skip modes already identified above
                if modeID_table[modeID,i]>0:
                    continue
                if verbose:
                    print('------------------------- LEVEL 2 - LOOKING FOR MODE ',modeIDName)

                found = False;
                # --- Loop on all non-identified modes in increasing order
                im = 0;
                for im, mode in enumerate(Modes):
                    if modesIdentified[im] or modesSkipped[im]:
                        continue # move to next mode
                    # List of component descriptions where this mode has maximum values
                    stateMax=np.argwhere((mode['StateHasMaxAtThisMode']==1)).flatten()
                    descriptions     = np.array(mode['DescStates'])[stateMax]
                    descriptions     = descriptions[:7] # we keep only the first 7 descriptions
                    descriptionsED   = [d for d in descriptions  if d.startswith('ED')]
                    descriptionsBD   = [d for d in descriptions  if d.startswith('BD')]
                    descriptionsAD   = [d for d in descriptions  if d.startswith('AD')]
                    descriptionsMisc = [d for d in descriptions  if d not in descriptionsED+descriptionsBD+descriptionsAD]
                    descriptions     = descriptionsED+descriptionsBD+descriptionsMisc # NOTE: we skipp AD descriptions
                    j = 0;
                    for description in descriptions:
                        found, pattern = doesDescriptionMatchPatternList(description, patternList)
                        if found:
                            if verbose:
                                print('L2 Mode {} identified using pattern {}'.format(im+1,pattern))
                            modeID_table[modeID,i] = im+1 # Using Matlab Indexing
                            modesDescLocal[modeID] = [None,]   # we empty this mode patternlist so that it cannot be matched again
                            modesIdentified[im] = True;
                            break
                    if found:
                        break
            if verbose:
                print('>> modeIDTable',modeID_table[:,i])
            # We disqualify modes that had max and that didn't match anything:
            for im, mode in enumerate(Modes):
                if modesIdentified[im] or modesSkipped[im]:
                    continue # move to next mode
                stateMax=np.argwhere((mode['StateHasMaxAtThisMode']==1)).flatten()
                if len(stateMax)>=1:
                    modesSkipped[im]=True
                    shortdescr = CampbellData[i]['ShortModeDescr'][im]
                    if verbose:
                        if shortdescr.find('ED')>=0:
                            print('>>>> short', CampbellData[i]['ShortModeDescr'][im])
                            print('>>>> Problem in IdentifyModes. ED DOF found in level 2')

        if 3 in Levels:
            # --- Level 3 - Try our best for modes with no max
            # Loop on modes to be identified
            for modeID, modeIDdescList in enumerate(modesDescLocal):
                modeIDName  = modeIDdescList[0]
                patternList = modeIDdescList[1:]

                # Skip modes already identified above
                if modeID_table[modeID,i]>0:
                    continue
                if verbose:
                    print('------------------------- LEVEL 3 - LOOKING FOR MODE ',modeIDName)

                found = False;
                # --- Loop on all non-identified modes in increasing order
                im = 0;
                while not found and im < nModes: # Loop on modes
                    mode = Modes[im]
                    if modesIdentified[im] or modesSkipped[im]:
                        pass
                    else:
                        # --- Otherwise, use as mode descriptions the other ones. Seems weird
                        stateMax=np.argwhere((mode['StateHasMaxAtThisMode']==0)).flatten()
                        descriptions     = np.array(mode['DescStates'])[stateMax]
                        ADcounts = np.sum([s.startswith('AD') for s in descriptions[:5]])

                        descriptions2= np.array(mode['DescStates'])[mode['StateHasMaxAtThisMode']]
                        if len(descriptions2) == 0:
                            noMax=True
                            descriptions3 = mode['DescStates'][:5]
                        else:
                            noMax=False
                        if ADcounts<5:
                            descriptions=[d for d in descriptions if not d.startswith('AD')]
#                                 descriptionsED   = [d for d in descriptions  if d.startswith('ED')]
#                                 descriptionsBD   = [d for d in descriptions  if d.startswith('BD')]
#                                 descriptionsAD   = [d for d in descriptions  if d.startswith('AD')]
#                                 descriptionsMisc = [d for d in descriptions  if d not in descriptionsED+descriptionsBD+descriptionsAD]
#                                 descriptions     = descriptionsED+descriptionsBD+descriptionsMisc # NOTE: we skipp AD descriptions
                            descriptions     = descriptions[:5] # we keep only the first 7 descriptions
                            if verbose:
                                print('>>> Mode',mode['index'], modesIdentified[im], modesSkipped[im])
                                print('>>>> descr', [replaceModeDescription(s) for s in descriptions])
                                print('>>>> short', CampbellData[i]['ShortModeDescr'][im])
                        else:
                            descriptions=[]
#                         #descriptions     = descriptions[:7] # we keep only the first 7 descriptions

                        j = 0;
                        while not found and j < len(descriptions):
                            j = j + 1;
                            if not found:
                                for targetDesc in patternList:
                                    # Looking for targetDesc into list of descriptions
                                    if re.search(targetDesc ,descriptions[j-1],re.IGNORECASE)!=None:
                                        modeID_table[modeID,i] = im+1 # Using Matlab Indexing
                                        if verbose:
                                            print('L3 Mode {} identified as {}'.format(im+1,targetDesc))
                                            print('     String was: ', descriptions[j-1])
                                        modesIdentified[im] = True;
                                        found = True;
                                        break;
                    im=im+1; # increment counter
            if verbose:
                print('>> modeIDTable',modeID_table[:,i])

        if verbose:
            print('---------- Summary')
            for j in np.arange(len(modeID_table)):
                print('{:32s}  {:d}'.format(modesDesc[j][0],modeID_table[j,i]))
            print('---------- ')


    return modeID_table,modesDesc

def IdentifiedModesDict(CampbellData,modeID_table,modesDesc):
    """
    To be called with the results of IdentifyModes.
    Create a list of dictionaries to more easily interprete the result
    """
    nOP = modeID_table.shape[1]
    modesInfoPerOP=[]
    for iOP in np.arange(nOP):
        modesInfo={}
        for i in np.arange(len(modesDesc)):
            desc = modesDesc[i][0]
            ID   = int(modeID_table[i,iOP])
            if ID==0:
                freq=np.nan
                damp=np.nan
                cont=''
            else:
                freq = np.around(CampbellData[iOP]['Modes'][ID-1]['NaturalFreq_Hz'],5)
                damp = np.around(CampbellData[iOP]['Modes'][ID-1]['DampingRatio'],5)
                cont = CampbellData[iOP]['ShortModeDescr'][ID-1]
            modesInfo[desc]={'ID':ID,'f0':freq,'zeta':damp,'cont':cont}
        modesInfoPerOP.append(modesInfo)
    return modesInfoPerOP

def extractShortModeDescription(mode):
    """ 
    Look at which modes have max, append these description, perform shortening substitution
    The description starts with noMax if no maximum exsits in this mode
    The ElastoDyn modes are placed first
    """
    descriptions = np.array(mode['DescStates'])[mode['StateHasMaxAtThisMode']]
    if len(descriptions) == 0:
        noMax=True
        descriptions = mode['DescStates'][:5]
    else:
        noMax=False
    sED   = [s for s in descriptions if s.startswith('ED')]
    sBD   = [s for s in descriptions if s.startswith('BD')]
    sMisc = [s for s in descriptions if s not in sED+sBD]
    sAll  = [replaceModeDescription(s) for s in sED+sBD+sMisc]
    shortdescr = ' - '.join(sAll)
    if noMax:
        shortdescr = 'NoMax - ' + shortdescr
    return shortdescr


def replaceModeDescription(s):
    """ Perform substitutions to make the mode description shorter"""
    s = s.replace('Blade','Bld')
    s = s.replace('blade','Bld')
    s = s.replace('First time derivative of','d/dt of')
    s = s.replace('fore-aft bending mode DOF, m','FA')
    s = s.replace('side-to-side bending mode DOF, m','SS')
    s = s.replace('bending-mode DOF of Bld ','')
    s = s.replace(' rotational-flexibility DOF, rad','-rot')
    s = s.replace('rotational displacement in ','rot')
    s = s.replace('translational displacement in ','trans')
    s = s.replace('Platform horizontal surge translation DOF','Platform surge')
    s = s.replace('Platform vertical heave translation DOF','Platform heave')
    s = s.replace('Platform pitch tilt rotation DOF','Platform pitch')
    s = s.replace(', rad','')
    s = s.replace(', m','')
    s = s.replace('finite element node ','N')
    s = s.replace('cosine','cos')
    s = s.replace('sine','sin')
    s = s.replace('flapwise','FLAP')
    s = s.replace('edgewise','EDGE')
    s = s.replace('collective','coll.')
    s = s.replace('rotZ','TORS-ROT')
    s = s.replace('transX','FLAP-DISP')
    s = s.replace('transY','EDGE-DISP')
    s = s.replace('rotX','EDGE')
    s = s.replace('rotY','FLAP')
    return s


def PrettyStateDescriptions(DescStates, ndof2, performedTransformation):
    idx=np.array(list(range(0,ndof2))+list(range(ndof2*2+1,len(DescStates))))    
    tmpDS = [DescStates[i] for i in idx]
    
    if performedTransformation:
        key_vals=[['BD_1','Blade collective'],['BD_2','Blade cosine'],['BD_3','Blade sine'],['blade 1','blade collective'], ['blade 2','blade cosine'], ['blade 3','Blade sine '],
                ['PitchBearing1','Pitch bearing collective '], ['PitchBearing2','Pitch bearing cosine '], ['PitchBearing3','Pitch bearing sine '],
                ['at Blade', 'at blade'],
                ['of Blade', 'of blade'],
                ]
        # Replace Substrings from String List 
        sub = dict(key_vals)
        for key, val in sub.items(): 
            for idx, ele in enumerate(tmpDS):
                if key in ele: 
                    tmpDS[idx] = ele.replace(key, val)

        StateDesc=tmpDS
    else:
        StateDesc = tmpDS
    
    for i in range(len( StateDesc )):
        First = re.split(r'\(',StateDesc[i],2)
        Last  = re.split(r'\)',StateDesc[i],2)

        if len(First)>0 and len(Last)>0 and len(First[0]) != len(StateDesc[i]) and len( Last[-1] ) != len(StateDesc[i]):
            StateDesc[i] = (First[0]).strip() + Last[-1];
            #print(StateDesc[i])
        
    return StateDesc



