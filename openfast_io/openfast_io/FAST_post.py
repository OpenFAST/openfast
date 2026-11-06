#*********************************************************************************************************************************
# LICENSING
# Copyright (C) 2024-2025  National Renewable Energy Laboratory
# Copyright (C) 2025-2026  National Laboratory of the Rockies
#
#    This file is part of openfast_io.
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
from __future__ import print_function
from rosco.toolbox.ofTools.fast_io.output_processing import output_processing
import rosco.toolbox

def FAST_IO_timeseries(fname):
    # interface to FAST_IO data load
    try:
        test = rosco.toolbox.__file__
    except:
        print('WARNING: rosco.toolbox required for wisdem.aeroelasticse.FAST_post.FAST_IO_timeseries')
    
    fast_out = output_processing()
    fast_data = fast_out.load_fast_out(fname, verbose=True)[0]
    return fast_data
