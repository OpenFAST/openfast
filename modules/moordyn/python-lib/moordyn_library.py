#**********************************************************************************************************************************
# LICENSING
# Copyright (C) 2021 National Renewable Energy Laboratory
# Author: Nicole Mendoza
#
# This file is part of MoorDyn.
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
#**********************************************************************************************************************************
#
# This is the Python-C interface library for MoorDyn
# Usage: THIS LIBRARY IS NOT TO BE CHANGED OR EDITED BY THE USER
from ctypes import (
    CDLL,
    POINTER,
    create_string_buffer,
    byref,
    c_byte,
    c_int,
    c_double,
    c_float, 
    c_char,
    c_char_p, 
    c_bool
)
import numpy as np
import datetime

class MoorDynLib(CDLL):

    # Human readable error levels
    error_levels = {
        0: "None",
        1: "Info",
        2: "Warning",
        3: "Severe Error",
        4: "Fatal Error"
    }

    #   NOTE:   the error message length in Fortran is controlled by the
    #           ErrMsgLen variable in the NWTC_Base.f90 file.  If that ever
    #           changes, it may be necessary to update the corresponding size
    #           here.
    error_msg_c_len = 1025

    #   NOTE:   the length of the name used for any output file written by the
    #           HD Fortran code is 1025.
    default_str_c_len = 1025

    def __init__(self, library_path):
        super().__init__(library_path)
        self.library_path = library_path

        self._initialize_routines()

        # Create buffers for class data
        self.abort_error_level = 4
        self.error_status_c = c_int(0)
        self.error_message_c = create_string_buffer(self.error_msg_c_len)
        self.error_message     = create_string_buffer(1025)
        self.ended             = False   # For error handling at end

        self._channel_names    = create_string_buffer(256*1000)
        self._channel_units    = create_string_buffer(256*1000)

        self.dt                = c_double(0)
        self.total_time        = c_double(0)
        self.numTimeSteps      = c_int(0)

    # Initialize routines ------------------------------------------------------------------------------------------------------------
    def _initialize_routines(self):
        self.MD_C_Init.argtypes = [
            POINTER(c_char_p),                    # IN: input file string
            POINTER(c_int),                       # IN: input file string length
            POINTER(c_double),                    # IN: dt
            POINTER(c_float),                     # IN: g
            POINTER(c_float),                     # IN: rho_water
            POINTER(c_float),                     # IN: depth_water
            POINTER(c_float),                     # IN: platform initial position
            POINTER(c_int),                       # IN: interpolation order
            POINTER(c_int),                       # OUT: number of channels
            POINTER(c_char),                      # OUT: output channel names
            POINTER(c_char),                      # OUT: output channel units
            POINTER(c_int),                       # OUT: ErrStat_C
            POINTER(c_char)                       # OUT: ErrMsg_C
        ]
        self.MD_C_Init.restype = c_int

        self.MD_C_CalcOutput.argtypes = [
            POINTER(c_double),                    # IN: Time @ n
            POINTER(c_float),                     # IN: Positions -- node positions    (1 x 6 array)  
            POINTER(c_float),                     # IN: Velocities -- node velocities  (1 x 6 array)
            POINTER(c_float),                     # IN: Accelerations -- node accelerations  (1 x 6 array)
            POINTER(c_float),                     # OUT: Forces (3 forces and 3 moments)
            POINTER(c_float),                     # OUT: Output Channel Values
            POINTER(c_int),                       # OUT: ErrStat_C
            POINTER(c_char)                       # OUT: ErrMsg_C
        ]
        self.MD_C_CalcOutput.restype = c_int
        
        self.MD_C_UpdateStates.argtypes = [
            POINTER(c_double),                    # IN: time @ n
            POINTER(c_double),                    # IN: time @ n+1
            POINTER(c_float),                     # IN: Positions -- node positions    (1 x 6 array)
            POINTER(c_float),                     # IN: Velocities -- node velocities  (1 x 6 array)
            POINTER(c_float),                     # IN: Accelerations -- node accelerations  (1 x 6 array)
            POINTER(c_int),                       # OUT: ErrStat_C
            POINTER(c_char)                       # OUT: ErrMsg_C
        ]
        self.MD_C_UpdateStates.restype = c_int

        self.MD_C_End.argtypes = [
            POINTER(c_int),                       # OUT: ErrStat_C
            POINTER(c_char)                       # OUT: ErrMsg_C
        ]
        self.MD_C_End.restype = c_int

    # md_init ------------------------------------------------------------------------------------------------------------
    def md_init(self, input_string_array, g, rho_water, depth_water, platform_init_pos, interpOrder):

        # Convert the string into a c_char byte array
        input_string = '\x00'.join(input_string_array)
        input_string = input_string.encode('utf-8')
        input_string_length = len(input_string)

        # Convert the initial positions array into c_float array
        init_positions_c = (c_float * 6)(0.0, )
        for i, p in enumerate(platform_init_pos):
            init_positions_c[i] = c_float(p)

        self._numChannels = c_int(0)

        self.MD_C_Init(
            c_char_p(input_string),                # IN: input file string
            byref(c_int(input_string_length)),     # IN: input file string length
            byref(c_double(self.dt)),              # IN: time step (dt)
            byref(c_float(g)),                     # IN: g
            byref(c_float(rho_water)),             # IN: rho_water
            byref(c_float(depth_water)),           # IN: depth_water
            init_positions_c,                      # IN: platform initial position
            byref(c_int(interpOrder)),             # IN: interpolation order
            byref(self._numChannels),              # OUT: number of channels
            self._channel_names,                   # OUT: output channel names
            self._channel_units,                   # OUT: output channel units
            byref(self.error_status_c),            # OUT: ErrStat_C
            self.error_message_c                   # OUT: ErrMsg_C
        )
        
        self.check_error()

    # md_calcOutput ------------------------------------------------------------------------------------------------------------
    def md_calcOutput(self, t, positions, velocities, accelerations, forces, output_channel_values):

        positions_c = (c_float * 6)(0.0,)
        for i, p in enumerate(positions):
            positions_c[i] = c_float(p)

        velocities_c = (c_float * 6)(0.0,)
        for i, p in enumerate(velocities):
            velocities_c[i] = c_float(p)

        accelerations_c = (c_float * 6)(0.0,)
        for i, p in enumerate(accelerations):
            accelerations_c[i] = c_float(p)

        forces_c = (c_float * 6)(0.0,)
        for i, p in enumerate(forces):
            forces_c[i] = c_float(p)

        outputs_c = (c_float * self._numChannels.value)(0.0,)
        for i, p in enumerate(output_channel_values):
            outputs_c[i] = c_float(p)

        self.MD_C_CalcOutput(
            byref(c_double(t)),                    # IN: time
            positions_c,                           # IN: positions
            velocities_c,                          # IN: velocities
            accelerations_c,                       # IN: accelerations
            forces_c,                              # OUT: forces
            outputs_c,                             # OUT: output channel values
            byref(self.error_status_c),            # OUT: ErrStat_C
            self.error_message_c                   # OUT: ErrMsg_C
        )

        for i in range(0,len(forces_c)):
            forces[i] = c_float(forces_c[i]).value

        for i in range(0,len(outputs_c)):
            output_channel_values[i] = c_float(outputs_c[i]).value
        
        self.check_error()

    # md_updateStates ------------------------------------------------------------------------------------------------------------
    def md_updateStates(self, t1, t2, positions, velocities, accelerations):

        positions_c = (c_float * 6)(0.0,)
        for i, p in enumerate(positions):
            positions_c[i] = c_float(p)

        velocities_c = (c_float * 6)(0.0,)
        for i, p in enumerate(velocities):
            velocities_c[i] = c_float(p)

        accelerations_c = (c_float * 6)(0.0,)
        for i, p in enumerate(accelerations):
            accelerations_c[i] = c_float(p)

        self.MD_C_UpdateStates(
            byref(c_double(t1)),                   # IN: current time (t)
            byref(c_double(t2)),                   # IN: next time step (t+1)
            positions_c,                           # IN: positions
            velocities_c,                          # IN: velocities
            accelerations_c,                       # IN: accelerations
            byref(self.error_status_c),            # OUT: ErrStat_C
            self.error_message_c                   # OUT: ErrMsg_C
        )
        
        self.check_error()

    # md_end ------------------------------------------------------------------------------------------------------------
    def md_end(self):

        if not self.ended:
            self.ended = True
            self.MD_C_End(
                byref(self.error_status_c),            # OUT: ErrStat_C
                self.error_message_c                   # OUT: ErrMsg_C
            )
            self.check_error()

    #===============================================================================
    # OTHER FUNCTIONS --------------------------------------------------------------------------------------------------

    # Error Handling Functions
    def check_error(self):
        if self.error_status_c.value == 0:
            return
        elif self.error_status_c.value < self.abort_error_level:
            print(f"MoorDyn error status: {self.error_levels[self.error_status_c.value]}: {self.error_message_c.value.decode('ascii')}")
        else:
            print(f"MoorDyn error status: {self.error_levels[self.error_status_c.value]}: {self.error_message_c.value.decode('ascii')}")
            self.md_end()
            raise Exception("\nMoorDyn terminated prematurely.")

    # Output Channel Functions
    @property
    def output_channel_names(self):
        if len(self._channel_names.value.split()) == 0:
             return []
        output_channel_names = self._channel_names.value.split()
        output_channel_names = [n.decode('UTF-8') for n in output_channel_names]
        return output_channel_names

    @property
    def output_channel_units(self):
        if len(self._channel_units.value.split()) == 0:
            return []
        output_channel_units = self._channel_units.value.split()
        output_channel_units = [n.decode('UTF-8') for n in output_channel_units]
        return output_channel_units

#===============================================================================
#   Helper class for writing output channels to file.
#   For the regression testing to mirror the output from the InfowWind Fortran
#   driver.  This may also have value for debugging the interfacing to MD.

class WriteOutChans():
    """
    This is only for testing purposes. Since we are not returning the
    output channels to anything, we will write them to file.  When coupled to
    another code, this data would be passed back for inclusion the any output
    file there.
    """
    def __init__(self,filename,chan_names,chan_units):
        chan_names.insert(0,'Time')             # add time index header
        chan_units.insert(0,'(s)')              # add time index unit
        self.OutFile=open(filename,'wt')        # open output file and write header info
        # write file header
        t_string=datetime.datetime.now()
        dt_string=datetime.date.today()
        self.OutFile.write(f"## This file was generated by MoorDyn c-bindings library on {dt_string.strftime('%b-%d-%Y')} at {t_string.strftime('%H:%M:%S')}\n")
        self.OutFile.write(f"## This file contains output channels requested from the OutList section of the input file")
        self.OutFile.write(f"{filename}\n")
        self.OutFile.write("#\n")
        self.OutFile.write("#\n")
        self.OutFile.write("#\n")
        self.OutFile.write("#\n")
        l = len(chan_names)
        f_string = "{:^15s}"+"   {:^20s}  "*(l-1)
        self.OutFile.write(f_string.format(*chan_names) + '\n')
        self.OutFile.write(f_string.format(*chan_units) + '\n')
        self.opened = True

    def write(self,chan_data):
        l = chan_data.shape[1]
        f_string = "{:10.4f}"+"{:25.7f}"*(l-1)
        for i in range(0,chan_data.shape[0]):
            self.OutFile.write(f_string.format(*chan_data[i,:]) + '\n')
            #if i==0:
            #    print(f"{chan_data[i,:]}")

    def end(self):
        if self.opened:
            self.OutFile.close()
            self.opened = False

#===============================================================================
#   Helper class for debugging the interface.  This will write out all the
#   input position/orientation, velocities, accelerations, and the resulting
#   forces and moments at each input node.  If all is functioning correctly,
#   this will be identical to the corresponding values in the MoorDyn output
#   channels.

class DriverDbg():
    """
    This is only for debugging purposes only.  The input motions and resulting
    forces can be written to file with this class to verify the data I/O to the
    Fortran library. When coupled to another code, the force/moment array would
    be passed back to the calling code for use in the structural solver.
    """
    def __init__(self,filename):
        self.DbgFile=open(filename,'wt')        # open output file and write header info
        self.numNodePts=1
        # write file header
        t_string=datetime.datetime.now()
        dt_string=datetime.date.today()
        self.DbgFile.write(f"## This file was generated by moordyn_c_lib on {dt_string.strftime('%b-%d-%Y')} at {t_string.strftime('%H:%M:%S')}\n")
        self.DbgFile.write(f"## This file contains the resulting forces/moments on the substructure passed into the moordyn_c_lib\n")
        self.DbgFile.write("#\n")
        f_string = "{:^25s}"
        self.DbgFile.write("#        T     ")
        for i in range(1,2):
            f_num = "N{0:04d}_".format(i)
            self.DbgFile.write(f_string.format(f_num+"x"  ))
            self.DbgFile.write(f_string.format(f_num+"y"  ))
            self.DbgFile.write(f_string.format(f_num+"z"  ))
            self.DbgFile.write(f_string.format(f_num+"Rx" ))
            self.DbgFile.write(f_string.format(f_num+"Ry" ))
            self.DbgFile.write(f_string.format(f_num+"Rz" ))
            self.DbgFile.write(f_string.format(f_num+"Vx" ))
            self.DbgFile.write(f_string.format(f_num+"Vy" ))
            self.DbgFile.write(f_string.format(f_num+"Vz" ))
            self.DbgFile.write(f_string.format(f_num+"RVx"))
            self.DbgFile.write(f_string.format(f_num+"RVy"))
            self.DbgFile.write(f_string.format(f_num+"RVz"))
            self.DbgFile.write(f_string.format(f_num+"Ax" ))
            self.DbgFile.write(f_string.format(f_num+"Ay" ))
            self.DbgFile.write(f_string.format(f_num+"Az" ))
            self.DbgFile.write(f_string.format(f_num+"RAx"))
            self.DbgFile.write(f_string.format(f_num+"RAy"))
            self.DbgFile.write(f_string.format(f_num+"RAz"))
            self.DbgFile.write(f_string.format(f_num+"Fx" ))
            self.DbgFile.write(f_string.format(f_num+"Fy" ))
            self.DbgFile.write(f_string.format(f_num+"Fz" ))
            self.DbgFile.write(f_string.format(f_num+"Mx" ))
            self.DbgFile.write(f_string.format(f_num+"My" ))
            self.DbgFile.write(f_string.format(f_num+"Mz" ))
        self.DbgFile.write("\n")
        self.DbgFile.write("#       (s)    ")
        for i in range(1,2):
            self.DbgFile.write(f_string.format("(m)"      ))
            self.DbgFile.write(f_string.format("(m)"      ))
            self.DbgFile.write(f_string.format("(m)"      ))
            self.DbgFile.write(f_string.format("(rad)"    ))
            self.DbgFile.write(f_string.format("(rad)"    ))
            self.DbgFile.write(f_string.format("(rad)"    ))
            self.DbgFile.write(f_string.format("(m/s)"    ))
            self.DbgFile.write(f_string.format("(m/s)"    ))
            self.DbgFile.write(f_string.format("(m/s)"    ))
            self.DbgFile.write(f_string.format("(rad/s)"  ))
            self.DbgFile.write(f_string.format("(rad/s)"  ))
            self.DbgFile.write(f_string.format("(rad/s)"  ))
            self.DbgFile.write(f_string.format("(m/s^2)"  ))
            self.DbgFile.write(f_string.format("(m/s^2)"  ))
            self.DbgFile.write(f_string.format("(m/s^2)"  ))
            self.DbgFile.write(f_string.format("(rad/s^2)"))
            self.DbgFile.write(f_string.format("(rad/s^2)"))
            self.DbgFile.write(f_string.format("(rad/s^2)"))
            self.DbgFile.write(f_string.format("(N)"      ))
            self.DbgFile.write(f_string.format("(N)"      ))
            self.DbgFile.write(f_string.format("(N)"      ))
            self.DbgFile.write(f_string.format("(N-m)"    ))
            self.DbgFile.write(f_string.format("(N-m)"    ))
            self.DbgFile.write(f_string.format("(N-m)"    ))
        self.DbgFile.write("\n")
        self.opened = True

    def write(self,t,Positions,Velocities,Accelerations,Forces):
        t_string = "{:10.4f}"
        f_string = "{:25.7f}"*6
        self.DbgFile.write(t_string.format(t))
        self.DbgFile.write(f_string.format(*Positions[:]))
        self.DbgFile.write(f_string.format(*Velocities[:]))
        self.DbgFile.write(f_string.format(*Accelerations[:]))
        self.DbgFile.write(f_string.format(*Forces[:]))
        self.DbgFile.write("\n")

    def end(self):
        if self.opened:
            self.DbgFile.close()
            self.opened = False
