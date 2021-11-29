#
# Copyright 2017 National Renewable Energy Laboratory
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

"""
    This library contains utility functions for the custom python programs making
    up the regression test system.
"""

import sys
import os
from stat import ST_MODE
import shutil 

def exitWithError(error, code=1):
    print(error)
    sys.exit(code)

def validInput(argv, nArgsExpected):
    valid = True if len(argv) == nArgsExpected else False
    return valid

def validateInputOrExit(argv, nArgsExpected, usage):
    if len(argv) != nArgsExpected:
        exitWithError(
            "Error: {} arguments given, expected {}\n".format(len(argv), nArgsExpected) +
            "Usage: {}".format(usage)
        )

def validateFileOrExit(path):
    if not os.path.isfile(path):
        exitWithError("Error: file does not exist at {}".format(path))

def validateDirOrExit(path):
    if not os.path.isdir(path):
        exitWithError("Error: directory does not exist at {}".format(path))

def validateDirOrMkdir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def validateExeOrExit(path):
    # On windows, check if maybe user forgot the .exe
    if os.path.isfile(path+'.exe'):
        path=path+'.exe'
    validateFileOrExit(path)
    permissionsMask = oct(os.stat(path)[ST_MODE])[-1:]
    if not int(permissionsMask)%2 == 1:
        exitWithError("Error: executable at {} does not have proper permissions.".format(path))

def copyTree(src, dst, excludeExt=[], renameDict={}, renameExtDict={}, includeExt=None):
    """ 
    Copy a directory to another one, overwritting files if necessary.
    copy_tree from distutils and copytree from shutil fail on Windows (in particular on git files)
    INPUTS:
     - src: source directory
     - dst: destination directory where files will be written/overwritten
     - includeExt: if provided, list of file extensions used for the copy
     - excludeExt: if provided, list of file extensions to be excluded from the copy
     - renameDict: dictionary used to rename files (the key is replaced by the value)
     - renameExt: dictionary used to rename extensions (the key is replaced by the value)
    """
    def forceMergeFlatDir(srcDir, dstDir):
        if not os.path.exists(dstDir):
            os.makedirs(dstDir)
        for item in os.listdir(srcDir):
            srcFile = os.path.join(srcDir, item)
            dstFile = os.path.join(dstDir, item)
            forceCopyFile(srcFile, dstFile)

    def forceCopyFile (sfile, dfile):
        # ---- Handling error due to wrong mod
        if os.path.isfile(dfile):
            if not os.access(dfile, os.W_OK):
                os.chmod(dfile, stat.S_IWUSR)
        #print(sfile, ' > ', dfile)
        shutil.copy2(sfile, dfile)

    def isAFlatDir(sDir):
        for item in os.listdir(sDir):
            sItem = os.path.join(sDir, item)
            if os.path.isdir(sItem):
                return False
        return True

    if includeExt is not None and len(excludeExt)>0:
        raise Exception('Provide includeExt or excludeExt, not both')

    for item in os.listdir(src):
        filebase, ext = os.path.splitext(item)
        if ext in excludeExt:
            continue
        if includeExt is not None:
            if ext not in includeExt:
                continue
        s = os.path.join(src, item)
        if item in renameDict.keys():
            item = renameDict[item] # renaming filename base on rename dict
        if ext in renameExtDict.keys():
            item = filebase + renameExtDict[ext] # changing extension based on rename ext dict
        d = os.path.join(dst, item)
        if os.path.isfile(s):
            if not os.path.exists(dst):
                os.makedirs(dst)
            forceCopyFile(s,d)
        if os.path.isdir(s):
            isRecursive = not isAFlatDir(s)
            if isRecursive:
                copyTree(s, d)
            else:
                forceMergeFlatDir(s, d)


def deleteOutputs(inputFile, extensions=['.out','.ech','.yaml','.sum']):
    """ 
    Delete output files from a given OpenFAST/driver input file 
    assuming the outputs have the same filebase as the inputfile
    """
    filebase, ext = os.path.splitext(inputFile)
    for e in extensions:
        outputFilename = filebase+e
        if os.path.exists(outputFilename):
            try:
                os.remove(outputFilename)
            except:
                print('[FAIL] deleting ',outputFilename)


