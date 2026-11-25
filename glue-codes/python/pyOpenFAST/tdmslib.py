# -*- coding: utf-8 -*-
#*********************************************************************************************************************************
# LICENSING
# Copyright (C) 2025  National Renewable Energy Laboratory
# Copyright (C) 2025-2026  National Laboratory of the Rockies
#
#    This file is part of OpenFAST.
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
#*********************************************************************************************************************************
"""
Created on Wed Jun  5 11:26:52 2024

@author: schamot

Description:
    Code takes a .tdms file path as input and outputs the data as a dictionary.
    
    Parameters
    ----------
    path : str
        path to the .tdms file.

    Returns
    -------
    pyDict : dict
        Python dictionary of the tdms file.

"""

import nptdms

def main():
    fileName = "Oscilloscope Data/UF HASEL tests/UF2down_3kVcycleForce.tdms"
    output = TdmsToDict(fileName)
    print(output)

def TdmsToDict(path):
    '''
    Code takes a .tdms file path as input and outputs the data as a dictionary.

    Parameters
    ----------
    path : str
        path to the .tdms file.

    Returns
    -------
    pyDict : dict
        Python dictionary of the tdms file.

    '''
    # Get the files requested
    tdmsFile = path
    
    # open/read the tdms file
    tdms_file = nptdms.TdmsFile.read(tdmsFile)

    # Set up python dictionary
    pyDict = {}    
    # Get group names
    for group in tdms_file.groups():
        #print(f'''Group: {group.name}''')
        # Create group dictionary key
        pyDict[group.name] = {}
        # Get channel names per group
        for channel in group.channels():
            #print(f'''\tChannel: {channel.name}''')
            pyDict[group.name][channel.name] = channel[:]
            # if the data has properties add them
            if channel.properties != {}:
                pyDict[group.name][channel.name + "_properties"] = channel.properties

    return pyDict

if __name__ == "__main__":
    main()