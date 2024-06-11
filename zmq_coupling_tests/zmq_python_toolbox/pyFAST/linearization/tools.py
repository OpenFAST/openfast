""" 
Simple tools to assist in doing linearization analyses with OpenFAST
"""
import numpy as np
import re
import os
import glob
import struct

from pyFAST.linearization.mbc import fx_mbc3, formatModesForViz
from pyFAST.linearization.campbell_data import campbell_diagram_data_oneOP # 


def getCampbellDataOP(fstFile_or_linFiles, writeModes=None, BladeLen=None, TowerLen=None, 
        removeTwrAzimuth=False, starSub=None, removeStatesPattern=None, verbose=False, 
        writeViz=False, **kwargs):
    """ 
    Return Campbell Data at one operating point from a .fst file or a list of lin files
    INPUTS:
     - fstFile_or_linFiles: filename of one .fst file, or, list of .lin files
                examples:   'main.fst'  or ['main.1.lin', 'main.4.lin']

     - writeModes: if True, a binary file and a .viz file is written to disk for OpenFAST VTK visualization.
                   if None, the binary file is written only if a checkpoint file is present.
                          For instance, if the main file is     : 'main.fst', 
                          the binary file will be               : 'main.ModeShapeVTK.pyPostMBC'
                          the viz file will be                  : 'main.ModeShapeVTK.viz'
                          the checkpoint file is expected to be : 'main.ModeShapeVTK.chkp'

     - BladeLen: blade length needed to scale the Campbell diagram data. 
                 if None: the length is inferred by reading the .fst file (and ED file)
     - TowerLen: tower length needed to scale the Campbell diagram data. 
                 if None: the length is inferred by reading the .fst file (and ED file)

     - starSub: if None, raise an error if `****` are present
                otherwise replace *** with `starSub` (e.g. 0)
                see FASTLinearizationFile. 
     - removeStatesPattern: remove states matching a giving description pattern.
               e.g:  'tower|Drivetrain'  or '^AD'
               see FASTLinearizationFile. 
     - removeTwrAzimuth: if False do nothing
                otherwise discard lin files where azimuth in [60, 180, 300]+/-4deg (close to tower). 
     - verbose: if True, more info is written to stdout

     - **kwargs: list of key/values to be passed to writeVizFile (see function below)
           VTKLinModes=15, VTKLinScale=10, VTKLinTim=1, VTKLinTimes1=True, VTKLinPhase=0, VTKModes=None

    OUTPUTS:
     - CDDOP: campbell diagram data for operating point (dictionary with many keys)
     - MBCOP: MBC data for operating point (dictionary)
    """
    # --- Figure out if the user provided a .fst file or a list of .lin files
    fstFile, linFiles =  getFST_and_LinFiles(fstFile_or_linFiles, verbose=verbose)

    # --- Open lin files for given OP/fst file, perform MBC
    MBCOP, matData = getMBCOP(fstFile=fstFile, linFiles=linFiles, verbose=verbose, removeTwrAzimuth=removeTwrAzimuth, starSub=starSub, removeStatesPattern=removeStatesPattern)
    if MBCOP is None:
        return None, None

    # --- Check if checkpoint file exists. If it does, we write the Modes
    fullpathbase, ext = os.path.splitext(fstFile)
    fullpath_chkp  = fullpathbase + '.ModeShapeVTK.chkp'
    fullpath_modes = fullpathbase + '.ModeShapeVTK.pyPostMBC'
    if writeModes is None:
        writeModes = os.path.exists(fullpath_chkp)

    # --- Write Modes for OpenFAST VTK visualization
    if writeModes:
        writeMBCOPForViz(MBCOP, matData, fullpath_modes, verbose=verbose)
        MBCOP['modeFile'] = fullpath_modes
    else:
        MBCOP['modeFile'] = None

    # --- Write Viz file if requested (NOTE: better to do this outside)
    if writeViz:
        vizfile = writeVizFile(fstFile, verbose=verbose, **kwargs)
        MBCOP['vizFile'] = vizfile
    else:
        MBCOP['vizFile'] = None

    # --- Estimate blade and tower length for scaling
    if BladeLen is None and TowerLen is None:
        BladeLen, TowerLen = estimateLengths(fstFile, verbose=verbose)

    # --- put data into "CampbellData" format
    CDDOP = campbell_diagram_data_oneOP(MBCOP, BladeLen, TowerLen)

    return CDDOP, MBCOP


def getCampbellDataOPs(fstFiles, BladeLen=None, TowerLen=None, verbose=False, **kwargs):
    """ 
    Return Campbell Data at several operating points from a list of .fst files
        see getCampbellDataOP for input arguments
    """
    # --- Estimate blade length and tower length for scaling
    if BladeLen is None and TowerLen is None:
        BladeLen, TowerLen = estimateLengths(fstFiles[0], verbose=verbose)

    # --- Run MBC for all operating points
    MBC = []
    CDD = []
    for i_lin, fstFile in enumerate(fstFiles):
        CDDOP, MBCOP = getCampbellDataOP(fstFile, BladeLen=BladeLen, TowerLen=TowerLen, verbose=verbose, **kwargs)
        if MBCOP is not None:
            CDD.append(CDDOP)
            MBC.append(MBCOP)
    # Remove missing data
    if len(CDD)==0:
        raise Exception('No linearization file found')
    return CDD, MBC

def getMBCOP(fstFile, linFiles=None, verbose=False, removeTwrAzimuth=False, starSub=None, removeStatesPattern=None):
    """ 
    Run necessary MBC for an OpenFAST file (one operating point)

    INPUTS:
     - fstFile: main openfast `.fst` filename 
     - linFiles: list of linfiles, inferred from fstfile if None provided
     - starSub: if None, raise an error if `****` are present
                otherwise replace *** with `starSub` (e.g. 0)
                see FASTLinearizationFile. 
     - removeStatesPattern: remove states matching a giving description pattern.
               e.g. r'^AD' : remove the AeroDyn states
               see FASTLinearizationFile. 
     - removeTwrAzimuth: if False do nothing
                otherwise discard lin files where azimuth in [60, 180, 300]+/-4deg (close to tower). 

    """

    # --- Find available lin files
    if linFiles is None:
        linFiles = findLinFiles(fstFile, verbose=verbose)

    # --- run MBC3 and campbell post_pro on lin files, generate postMBC file if needed 
    if len(linFiles)>0:
        MBC, matData = fx_mbc3(linFiles, verbose=False, removeTwrAzimuth=removeTwrAzimuth, starSub=starSub, removeStatesPattern=removeStatesPattern)
    else:
        return None, None

    return MBC, matData

def getMBCOPs(fstfiles, verbose=True, removeTwrAzimuth=False):
    """
    Run MBC transform on set of openfast linear outputs (multiple operating points)

    INPUTS:
      - fstfiles: list of .fst files
    """
    MBC = [None]*len(fstfiles)
    for i_lin, fstfile in enumerate(fstfiles):
        # MBC for a given operating point (OP)
        MBC[i_lin], matData = getMBCOP(fstfile, verbose=verbose, removeTwrAzimuth=removeTwrAzimuth)
    return MBC



def getFST_and_LinFiles(fstFile_or_linFiles, verbose=False):
    """ 
    Given a .fst or a list of .lin files, return both:
    - if a fst file is provided, .lin file next to it are sought for
    - if a list of line files are provided, return the .fst file from which they originated
    """
    if isinstance(fstFile_or_linFiles,str):
        # The user provided a string, we expect it's a .fst file
        fullpathbase, ext = os.path.splitext(fstFile_or_linFiles)
        if ext.lower()!='.fst':
            raise Exception('Provide either one fst file, or a list of .lin files')
        fstFile = fstFile_or_linFiles
        # --- Find available lin files
        linFiles = findLinFiles(fstFile, verbose=verbose)
    else:
        # the user provided a list (hopefully)
        fullpathbase, ext = os.path.splitext(fstFile_or_linFiles[0])
        if ext.lower()!='.lin':
            print(fstFile_or_linFiles)
            raise Exception('Provide either one fst file, or a list of .lin files')
        linFiles = fstFile_or_linFiles
        fullpathbase, ext = os.path.splitext(fullpathbase)
        fstFile = fullpathbase+'.fst'

    return fstFile, linFiles

def writeVizFile(fstFile, VTKLinModes=15, VTKLinScale=10, VTKLinTim=1, VTKLinTimes1=True, VTKLinPhase=0, VTKModes=None, verbose=False):
    fullpathbase, ext = os.path.splitext(fstFile)
    filebase          = os.path.basename(fullpathbase)
    fullpath_viz      = fullpathbase + '.ModeShapeVTK.viz'
    base_chkp         = filebase     + '.ModeShapeVTK'
    base_modes        = filebase     + '.ModeShapeVTK.pyPostMBC'

    if VTKModes is None:
        VTKModes=','.join((np.arange(VTKLinModes)+1).astype(str))

    with open(fullpath_viz, 'w') as f:
        f.write('------- OpenFAST MODE-SHAPE INPUT FILE -------------------------------------------\n');
        f.write('# Options for visualizing mode shapes\n');
        f.write('---------------------- FILE NAMES ----------------------------------------------\n');
        f.write('"{:s}"   CheckpointRoot - Rootname of the checkpoint file written when OpenFAST generated the linearization files (without the ".chkp" extension)\n'.format(base_chkp))
        f.write('"{:s}"   ModesFileName - Name of the mode-shape file (with eigenvectors)\n'.format(base_modes))
        f.write('---------------------- VISUALIZATION OPTIONS -----------------------------------\n')
        f.write('{:d}        VTKLinModes   - Number of modes to visualize (0 <= VTKLinModes <= NumModes)\n'.format(VTKLinModes))
        f.write('{:s}        VTKModes      - List of which VTKLinModes modes will be visualized (modes will be added sequentially from the last value entered)\n'.format(VTKModes))
        f.write('{:f}        VTKLinScale   - Mode shape visualization scaling factor (exaggerates mode shapes: try 10 for ElastoDyn; 0.1 for BeamDyn)\n'.format(VTKLinScale)) 
        f.write('{:d}        VTKLinTim     - Switch to make one animation for all LinTimes together (VTKLinTim=1) or separate animations for each LinTimes (VTKLinTim=2)\n'.format(VTKLinTim))
        f.write('{}      VTKLinTimes1  - If VTKLinTim=2, visualize modes at LinTimes(1) only? (if false, files will be generated at all LinTimes)\n'.format(VTKLinTimes1))
        f.write('{:f}       VTKLinPhase   - Phase used when making one animation for all LinTimes together (used only when VTKLinTim=1)\n'.format(VTKLinPhase))
#             if isnan(opts.VTKLinScale)
#                 % Then user didn't specify it, we use some logic
#                 if CompElast==1 % ElastoDyn - VTKLinScale=10
#                     fprintf(fid,'10    VTKLinScale   - Mode shape visualization scaling factor (exaggerates mode shapes: try 10 for ElastoDyn; 0.1 for BeamDyn)\n');
#                 elseif CompElast==2 % BeamDyn - VTKLinScale=0.1
#                     fprintf(fid,'0.1   VTKLinScale   - Mode shape visualization scaling factor (exaggerates mode shapes: try 10 for ElastoDyn; 0.1 for BeamDyn)\n'); 
#                     fprintf(fid,'%f    VTKLinScale   - Mode shape visualization scaling factor (exaggerates mode shapes: try 10 for ElastoDyn; 0.1 for BeamDyn)\n',opts.VTKLinScale); 
#             end
#             if isnan(opts.VTKLinTim)
#                 % The user didn't specify this, we use some logic
#                 fprintf(fid,'2         VTKLinTim     - Switch to make one animation for all LinTimes together (VTKLinTim=1) or separate animations for each LinTimes (VTKLinTim=2)\n');
#             else
#                 if (RPM<1e-3) % When RPM =0, VTKLinTim=1 would only produce one VTK
#                     fprintf(fid,'2         VTKLinTim     - Switch to make one animation for all LinTimes together (VTKLinTim=1) or separate animations for each LinTimes (VTKLinTim=2)\n');
#                 else
#                     fprintf(fid,'%d        VTKLinTim     - Switch to make one animation for all LinTimes together (VTKLinTim=1) or separate animations for each LinTimes (VTKLinTim=2)\n',opts.VTKLinTim);
#                 end
#             end
#             if length(opts.VTKLinTimes1)==0
#                 % The user didn't specify this
#                 fprintf(fid,'true      VTKLinTimes1  - If VTKLinTim=2, visualize modes at LinTimes(1) only? (if false, files will be generated at all LinTimes)\n');
#             else
#                 fprintf(fid,'%s        VTKLinTimes1  - If VTKLinTim=2, visualize modes at LinTimes(1) only? (if false, files will be generated at all LinTimes)\n',opts.VTKLinTimes1);
#             end
#             fclose(fid);
        if verbose:
            print(' Written Viz File:', fullpath_viz)
    return fullpath_viz

def writeVizFiles(fstFiles, **kwargs):
    """ write viz file for a set of fst files
    see writeVizFile
    """
    return [writeVizFile(fst, **kwargs) for fst in fstFiles]

# ------------------------------------------------------------------------
def fread(fid, n, type):
    """ Mimic the matlab function fread"""
    fmt, nbytes = {'uint8': ('B', 1), 'int16':('h', 2), 'int32':('i', 4), 'float32':('f', 4), 'float64':('d', 8)}[type]
    v = struct.unpack(fmt * n, fid.read(nbytes * n))
    if n==1:
        return v[0]
    else:
        return np.asarray(v)

def fwrite(fid, data, type):
    """ Mimic the matlab function fwrite"""
    # @ is used for packing in native byte order
    #  B - unsigned integer 8 bits
    #  h - integer 16 bits
    #  i - integer 32 bits
    #  f - float 32 bits
    #  d - float 64 bits
    fmt, _ = {'uint8': ('B', 1), 'int16':('h', 2), 'int32':('i', 4), 'float32':('f', 4), 'float64':('d', 8)}[type]
    if hasattr(data, '__len__'):
        data = data.flatten(order='F')
        n=len(data)
        fid.write(struct.pack('@'+str(n)+fmt, *data))
    else:
        fid.write(struct.pack('@'+fmt, data))


def writeMBCOPForViz(MBCOP, matData, modesFilename, nModesOut=None, nDigits=None, verbose=False, hack=False):
    VTK = formatModesForViz(MBCOP, matData, MBCOP['nb'], MBCOP['EigenVects_save'])
    return writeModesForViz(VTK, modesFilename, nModesOut=nModesOut, nDigits=nDigits, verbose=verbose, hack=hack)

# def writeModesForViz(f0, fd, zeta, Q, modesFilename, format='OFBinary'):
def writeModesForViz(VTK, modesFilename, nModesOut=None, nDigits=None, verbose=False, hack=False):
    """ 
    write binary file that will be read by OpenFAST to export modes to VTK
    """

    reFmt = 'float64' #8-byte real numbers
    nStates, nModes, nLinTimes = VTK['x_eig_magnitude'].shape
    if nModesOut is None:
        nModesOut=nModes

    #------- HACK to compare with matlab
    if hack:
        VTK['NaturalFreq_Hz'] =VTK['NaturalFreq_Hz']*0 + 1
        VTK['DampingRatio']   =VTK['DampingRatio']  *0 + 2
        VTK['DampedFreq_Hz']  =VTK['DampedFreq_Hz'] *0 + 3
        for iMode in range(nModes):
            VTK['x_eig_magnitude'][:,iMode,:] = np.zeros((nStates,nLinTimes)) + iMode+1
            VTK['x_eig_phase']    [:,iMode,:] = np.zeros((nStates,nLinTimes)) + iMode+1
            VTK['x_eig_magnitude'][2,iMode,:] = 12
            VTK['x_eig_phase']    [4,iMode,:] = 11
        nModesOut=1
    # ------END HACK
    # --- Reduce differences python/Matlab by rounding
    if nDigits is not None:
        res = 10**nDigits
        VTK['NaturalFreq_Hz']  = np.around(VTK['NaturalFreq_Hz'] , nDigits)
        VTK['DampingRatio']    = np.around(VTK['DampingRatio']   , nDigits)
        VTK['DampedFreq_Hz']   = np.around(VTK['DampedFreq_Hz']  , nDigits)
        VTK['x_eig_magnitude'] = np.around(VTK['x_eig_magnitude'], nDigits)
        VTK['x_eig_phase'    ] = np.around(VTK['x_eig_phase']    , nDigits)

    # --- Write to disk
    with open(modesFilename, 'wb') as fid:
        fwrite(fid, 1,        'int32' )# write a file identifier in case we ever change this format
        fwrite(fid, nModesOut, 'int32' )# number of modes (for easier file reading)
        fwrite(fid, nStates,  'int32' )# number of states (for easier file reading)
        fwrite(fid, nLinTimes,'int32' )# number of azimuths (i.e., LinTimes) (for easier file reading)
        # Freq and damping (not used in the FAST visualization algorithm)
        fwrite(fid, VTK['NaturalFreq_Hz'], reFmt)
        fwrite(fid, VTK['DampingRatio'],   reFmt)
        fwrite(fid, VTK['DampedFreq_Hz'],  reFmt)
        # Writing data mode by mode
        for iMode in range(nModesOut):
            fwrite(fid, VTK['x_eig_magnitude'][:,iMode,:], reFmt)
            fwrite(fid, VTK['x_eig_phase']    [:,iMode,:], reFmt)
    if verbose:
        print(' Written ModeFile:', modesFilename)

def readModesForViz(modesFilename):
    """ """
    reFmt = 'float64' #8-byte real numbers
    with open(modesFilename, 'rb') as fid:
        fformat  = fread(fid, 1, 'int32' ) # format identifier 
        nModes   = fread(fid, 1, 'int32' ) # number of modes (for easier file reading)
        nStates  = fread(fid, 1, 'int32' ) # number of states (for easier file reading)
        nLinTimes= fread(fid, 1, 'int32' ) # number of azimuths (i.e., LinTimes) (for easier file reading)
        # Freq and damping (not used in the FAST visualization algorithm)
        f0   = fread(fid, nModes, reFmt ) # number of modes (for easier file reading)
        zeta = fread(fid, nModes, reFmt ) # number of states (for easier file reading)
        fd   = fread(fid, nModes, reFmt ) # number of azimuths (i.e., LinTimes) (for easier file reading)
        # Reading data mode by mode
        Qmag = np.zeros((nStates, nModes, nLinTimes))
        Qphi = np.zeros((nStates, nModes, nLinTimes))
        for iMode in range(nModes):
            mag = fread(fid, nStates*nLinTimes, reFmt).reshape((nStates,nLinTimes), order='F')
            phi = fread(fid, nStates*nLinTimes, reFmt).reshape((nStates,nLinTimes), order='F')
            Qmag[:,iMode,:] = mag
            Qphi[:,iMode,:] = phi

    VTK = {}
    VTK['NaturalFreq_Hz']  = f0
    VTK['DampingRatio']    = zeta
    VTK['DampedFreq_Hz']   = fd
    VTK['x_eig_magnitude'] = Qmag
    VTK['x_eig_phase']     = Qphi

    return VTK


def findLinFiles(fstFile, verbose=False):
    """ 
    Find .lin files given a .fst file
    """
    fullpathbase, ext = os.path.splitext(fstFile)
    # NOTE: the code below is problematic for module lin files ED.1.lin 
    # So we do a re search to filter these out
    # First use glob
    lin_file_fmt    = '{}.*.lin'.format(fullpathbase)
    lin_files       = glob.glob(lin_file_fmt)
    # Then use re for stricter search
    lin_file_fmt_re    = r'.*\.[0-9]+\.lin'
    lin_files = glob_re(lin_file_fmt_re, lin_files)
    if len(lin_files)==0:
        if verbose:
            print('[WARN] Lin. files: {} ({})'.format(lin_file_fmt, len(lin_files)))
    else:
        if verbose:
            print('       Lin. files: {} ({})'.format(lin_file_fmt, len(lin_files)))
    return lin_files

def estimateLengths(fstFile, verbose=False):
    if os.path.exists(fstFile):
        # try to read BladeLen and TowerLen from fst file
        # TODO: can be done with pyFAST.io or AeroElasticSE
        # The interface is very similar
        from pyFAST.input_output.fast_input_deck import FASTInputDeck
        fst = FASTInputDeck(fstFile, 'ED')
        ED = fst.fst_vt['ElastoDyn']
        if ED is None:
            raise Exception('Unable to infer BladeLen and TowerLen, ElastoDyn file not found. FST file is: ',fstFile)
        BladeLen = ED['TipRad'] - ED['HubRad']
        TowerLen = ED['TowerHt']
    else:
        raise Exception('Provide `BladeLen` and `TowerLen`, or, an existing fst and ED file')
    if verbose:
        print('BladeLength (for scaling): ',BladeLen)
        print('TowerLength (for scaling): ',TowerLen)

    return BladeLen, TowerLen

def glob_re(pattern, strings):
    """ Apply a pattern to a list of strings 
    Typically used as "glob" and "re" to select a list of files matting a given mattern"""
    return list(filter(re.compile(pattern).match, strings))
