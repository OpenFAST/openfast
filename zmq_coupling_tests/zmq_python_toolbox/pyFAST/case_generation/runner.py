# --- For cmd.py
import os
import subprocess
import multiprocessing

import collections
import glob
import pandas as pd
import numpy as np
import shutil 
import stat
import re

# --- Fast libraries
from pyFAST.input_output.fast_input_file import FASTInputFile
from pyFAST.input_output.fast_output_file import FASTOutputFile

FAST_EXE='openfast'

# --------------------------------------------------------------------------------}
# --- Tools for executing FAST
# --------------------------------------------------------------------------------{
# --- START cmd.py
def run_cmds(inputfiles, exe, parallel=True, showOutputs=True, nCores=None, showCommand=True, flags=[], verbose=True): 
    """ Run a set of simple commands of the form `exe input_file`
    By default, the commands are run in "parallel" (though the method needs to be improved)
    The stdout and stderr may be displayed on screen (`showOutputs`) or hidden. 
    A better handling is yet required.
    """
    Failed=[]
    def _report(p):
        if p.returncode==0:
            if verbose:
                print('[ OK ] Input    : ',p.input_file)
        else:
            Failed.append(p)
            if verbose:
                print('[FAIL] Input    : ',p.input_file)
                print('       Directory: '+os.getcwd())
                print('       Command  : '+p.cmd)
                print('       Use `showOutputs=True` to debug, or run the command above.')
            #out, err = p.communicate()
            #print('StdOut:\n'+out)
            #print('StdErr:\n'+err)
    ps=[]
    iProcess=0
    if nCores is None:
        nCores=multiprocessing.cpu_count()
    if nCores<0:
        nCores=len(inputfiles)+1
    for i,f in enumerate(inputfiles):
        if len(flags)>0:
            f=flags + [f]
        #print('Process {}/{}: {}'.format(i+1,len(inputfiles),f))
        ps.append(run_cmd(f, exe, wait=(not parallel), showOutputs=showOutputs, showCommand=showCommand))
        iProcess += 1
        # waiting once we've filled the number of cores
        # TODO: smarter method with proper queue, here processes are run by chunks
        if parallel:
            if iProcess==nCores:
                for p in ps:
                    p.wait()
                for p in ps:
                    _report(p)
                ps=[]
                iProcess=0
    # Extra process if not multiptle of nCores (TODO, smarter method)
    for p in ps:
        p.wait()
    for p in ps:
        _report(p)
    # --- Giving a summary
    if len(Failed)==0:
        if verbose:
            print('[ OK ] All simulations run successfully.')
        return True, Failed
    else:
        print('[FAIL] {}/{} simulations failed:'.format(len(Failed),len(inputfiles)))
        for p in Failed:
            print('      ',p.input_file)
        return False, Failed

def run_cmd(input_file_or_arglist, exe, wait=True, showOutputs=False, showCommand=True):
    """ Run a simple command of the form `exe input_file` or `exe arg1 arg2`  """
    # TODO Better capture STDOUT
    if not os.path.exists(exe):
        raise Exception('Executable not found: {}'.format(exe))
    if isinstance(input_file_or_arglist, list):
        input_file     = ' '.join(input_file_or_arglist)
        args= [exe] + input_file_or_arglist
    else:
        input_file=input_file_or_arglist
        args= [exe,input_file]
    args = [a.strip() for a in args] # No surounding spaces, could cause issue
    shell=False
    if showOutputs:
        STDOut= None
    else:
        STDOut= open(os.devnull, 'w') 
    if showCommand:
        print('Running: '+' '.join(args))
    if wait:
        class Dummy():
            pass
        p=Dummy()
        p.returncode=subprocess.call(args , stdout=STDOut, stderr=subprocess.STDOUT, shell=shell)
    else:
        p=subprocess.Popen(args, stdout=STDOut, stderr=subprocess.STDOUT, shell=shell)
    # Storing some info into the process
    p.cmd            = ' '.join(args)
    p.args           = args
    p.input_file     = input_file
    p.exe            = exe
    return p

def runBatch(batchfiles, showOutputs=True, showCommand=True, verbose=True):
    """ 
    Run one or several batch files
    TODO: error handling, status, parallel
    """

    if showOutputs:
        STDOut= None
    else:
        STDOut= open(os.devnull, 'w') 

    curDir = os.getcwd()
    print('Current directory', curDir)

    def runOneBatch(batchfile):
        batchfile = batchfile.strip()
        batchfile = batchfile.replace('\\','/')
        batchDir  = os.path.dirname(batchfile)
        if showCommand:
            print('>>>> Running batch file:', batchfile)
            print('           in directory:', batchDir)
        try:
            os.chdir(batchDir)
            returncode=subprocess.call([batchfile], stdout=STDOut, stderr=subprocess.STDOUT, shell=shell)
        except:
            os.chdir(curDir)
            returncode=-10
        return returncode

    shell=False
    if isinstance(batchfiles,list):
        returncodes = []
        Failed=[]
        for batchfile in batchfiles:
            returncode = runOneBatch(batchfile)
            if returncode!=0:
                Failed.append(batchfile)
        if len(Failed)>0:
            returncode=1
            print('[FAIL] {}/{} Batch files failed.'.format(len(Failed),len(batchfiles)))
            print(Failed)
        else:
            returncode=0
            if verbose:
                print('[ OK ] {} batch filse ran successfully.'.format(len(batchfiles)))
        # TODO
    else:
        returncode = runOneBatch(batchfiles)
        if returncode==0:
            if verbose:
                print('[ OK ] Batch file ran successfully.')
        else:
            print('[FAIL] Batch file failed:',batchfiles)

    return returncode


# --- END cmd.py

def run_fastfiles(fastfiles, fastExe=None, parallel=True, showOutputs=True, nCores=None, showCommand=True, reRun=True, verbose=True):
    if fastExe is None:
        fastExe=FAST_EXE
    if not reRun:
        # Figure out which files exist
        newfiles=[]
        for f in fastfiles:
            base=os.path.splitext(f)[0]
            if os.path.exists(base+'.outb') or os.path.exists(base+'.out'):
                print('>>> Skipping existing simulation for: ',f)
                pass
            else:
                newfiles.append(f)
        fastfiles=newfiles

    return run_cmds(fastfiles, fastExe, parallel=parallel, showOutputs=showOutputs, nCores=nCores, showCommand=showCommand, verbose=verbose)

def run_fast(input_file, fastExe=None, wait=True, showOutputs=False, showCommand=True):
    if fastExe is None:
        fastExe=FAST_EXE
    return run_cmd(input_file, fastExe, wait=wait, showOutputs=showOutputs, showCommand=showCommand)


def writeBatch(batchfile, fastfiles, fastExe=None, nBatches=1, pause=False, flags='', flags_after='',
        run_if_ext_missing=None,
        discard_if_ext_present=None,
        dispatch=False,
        stdOutToFile=False,
        echo=True):
    """ Write one or several batch file, all paths are written relative to the batch file directory.
    The batch file will consist of lines of the form:
         [CONDITION] EXE [FLAGS] FILENAME [FLAGS_AFTER]

    INPUTS:
    - batchfile: path of the batch file to be written. 
                 If several files are requested (using nBatches) _i is inserted before the extension
    - nBatches: split into nBatches files.
    - pause: insert a pause statement at the end so that batch file is not closed after execution
    - flags: flags (string) to be placed between the executable and the filename
    - flags_after: flags (string) to be placed after the filename
    - run_if_ext_missing: add a line in the batch file so that the command is only run if
                          the file `f.EXT` is missing, where .EXT is specified in run_if_ext_missing
                          If None, the command is always run
    - discard_if_ext_present: similar to run_if_ext_missing, but this time, the lines are not written to the batch file
                          The test for existing outputs is done before writing the batch file
    - dispatch: if True, the input files are dispatched (the first nBatches files are dispathced on the nBatches)
    - stdOutToFile: if True, the output of the command is redirected to filename.stdout

    example:
       writeBatch('dir/MyBatch.bat', ['dir/c1.fst','dir/c2.fst'], 'op.exe', flags='-v', run_if_ext_missing='.outb')

       will generate a file with the following content:
         if not exist c1.outb (../of.exe c1.fst) else (echo "Skipping c1.fst")
         if not exist c2.outb (../of.exe c2.fst) else (echo "Skipping c2.fst")


    """
    if fastExe is None:
        fastExe=FAST_EXE
    fastExe_abs   = os.path.abspath(fastExe)
    batchfile_abs = os.path.abspath(batchfile)
    batchdir      = os.path.dirname(batchfile_abs)
    fastExe_rel   = os.path.relpath(fastExe_abs, batchdir)
    if len(flags)>0:
        flags=' '+flags
    if len(flags_after)>0:
        flags_after=' '+flags_after

    # Remove commandlines if outputs are already present
    if discard_if_ext_present:
        outfiles = [os.path.splitext(f)[0] + discard_if_ext_present for f in fastfiles]
        nIn=len(fastfiles)
        fastfiles =[f for f,o in zip(fastfiles,outfiles) if not os.path.exists(o)]
        nMiss=len(fastfiles)
        if nIn>nMiss:
            print('[INFO] WriteBatch: discarding simulations, only {}/{} needed'.format(nMiss, nIn))


    def writeb(batchfile, fastfiles):
        with open(batchfile,'w') as f:
            if not echo:
                if os.name == 'nt':
                    f.write('@echo off\n')
            for ff in fastfiles:
                ff_abs = os.path.abspath(ff)
                ff_rel = os.path.relpath(ff_abs, batchdir)
                cmd = fastExe_rel + flags + ' '+ ff_rel + flags_after
                if stdOutToFile:
                    stdout = os.path.splitext(ff_rel)[0]+'.stdout'
                    cmd += ' > ' +stdout
                if run_if_ext_missing is not None:
                    # TODO might be windows only
                    ff_out = os.path.splitext(ff_rel)[0] + run_if_ext_missing
                    if os.name == 'nt':
                        cmd = 'if not exist {} ({}) else (echo Skipping {})'.format(ff_out, cmd, ff_rel)
                    else:
                        cmd = 'if [[ ! -f {} ]] ; then {}; else echo Skipping {} ; fi'.format(ff_out, cmd, ff_rel)
                f.write("{:s}\n".format(cmd))
            if pause:
                f.write("pause\n") # might be windows only..



    if nBatches==1:
        writeb(batchfile, fastfiles)
        return batchfile
    else:

        if dispatch:
            # TODO this can probably be done with a one liner
            fastfiles2=[]
            for i in range(nBatches):
                fastfiles2+=fastfiles[i::nBatches]
            fastfiles = fastfiles2

        splits = np.array_split(fastfiles,nBatches)
        base, ext = os.path.splitext(batchfile)
        batchfiles=[]
        for i in np.arange(nBatches):
            batchfile = base+'_{:d}'.format(i+1) + ext
            batchfiles.append(batchfile)
            writeb(batchfile, splits[i])
        return batchfiles






def removeFASTOuputs(workDir):
    # Cleaning folder
    for f in glob.glob(os.path.join(workDir,'*.out')):
        os.remove(f)
    for f in glob.glob(os.path.join(workDir,'*.outb')):
        os.remove(f)
    for f in glob.glob(os.path.join(workDir,'*.ech')):
        os.remove(f)
    for f in glob.glob(os.path.join(workDir,'*.sum')):
        os.remove(f)

if __name__=='__main__':
    run_cmds(['main1.fst','main2.fst'], './Openfast.exe', parallel=True, showOutputs=False, nCores=4, showCommand=True)
    pass
    # --- Test of templateReplace

