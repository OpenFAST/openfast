#**********************************************************************************************************************************
# LICENSING
# Copyright (C) 2021 National Renewable Energy Laboratory
#
# This file is part of HydroDyn.
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
# This is the Python-C interface library for HydroDyn.  This may be used
# directly with Python based codes to call and run HydroDyn.  An example of
# using this library from Python is given in the accompanying Python driver
# program.  Additional notes and information on the interfacing is included
# there.
#
#   Note on angles:
#       All angles passed in are assumed to be given in radians as an Euler
#       angle sequence R(z)*R(y)*R(x) (notice the order as this important when
#       passing angles in). Written in matrix form as (doxygen formatted):
#
#          \f{eqnarray*}{   
#          M & = & R(\theta_z) R(\theta_y) R(\theta_x) \\
#            & = & \begin{bmatrix}  \cos(\theta_z) & \sin(\theta_z) & 0 \\
#                                  -\sin(\theta_z) & \cos(\theta_z) & 0 \\
#                                    0      &  0      & 1 \end{bmatrix}
#                  \begin{bmatrix}  \cos(\theta_y) & 0 & -\sin(\theta_y) \\
#                                         0 & 1 & 0        \\
#                                   \sin(\theta_y) & 0 & \cos(\theta_y)  \end{bmatrix}
#                  \begin{bmatrix}   1 &  0       & 0       \\
#                                    0 &  \cos(\theta_x) & \sin(\theta_x) \\
#                                    0 & -\sin(\theta_x) & \cos(\theta_x) \end{bmatrix} \\
#            & = & \begin{bmatrix}  
#             \cos(\theta_y)\cos(\theta_z) &   \cos(\theta_x)\sin(\theta_z)+\sin(\theta_x)\sin(\theta_y)\cos(\theta_z) &
#                                              \sin(\theta_x)\sin(\theta_z)-\cos(\theta_x)\sin(\theta_y)\cos(\theta_z) \\
#             -\cos(\theta_y)\sin(\theta_z)  & \cos(\theta_x)\cos(\theta_z)-\sin(\theta_x)\sin(\theta_y)\sin(\theta_z) & 
#                                              \sin(\theta_x)\cos(\theta_z)+\cos(\theta_x)\sin(\theta_y)\sin(\theta_z) \\
#             \sin(\theta_y)                & -\sin(\theta_x)\cos(\theta_y) & \cos(\theta_x)\cos(\theta_y) \\
#                  \end{bmatrix}   
#          \f}
#
#       When passed into the Fortran library, this Euler angle set is converted
#       into a DCM (direction cosine matrix) and stored on the input mesh.  All
#       calculations internally in HD are then performed using the DCM form of
#       the input angles.  
#
#       It should be noted that a small angle assumption when returning the
#       outputs for the platform roll, pitch, and yaw.  These angles are
#       assumed to be small enough that treating them as independent angles
#       does not introduce significant error in the output channels.  This may
#       yield a small discrepency between the values passed in and the returned
#       output channel values.  These output channels should not be directly
#       used -- they are only for reporting purposes.
#
from ctypes import (
    CDLL,
    POINTER,
    create_string_buffer,
    byref,
    c_int,
    c_double,
    c_float, 
    c_char,
    c_char_p,
)
import numpy as np
import datetime

class HydroDynLib(CDLL):
    # Human readable error levels from IfW.
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
        self.ended = False                  # For error handling at end


        # Create buffers for class data
        self.abort_error_level = 4
        self.error_status_c = c_int(0)
        self.error_message_c = create_string_buffer(self.error_msg_c_len)

        # This is not sufficient for HD
        #FIXME: ChanLen may not always be 20 -- could be as much as 256
        #       Possible fix is to pass this length over to Fortran side.
        #       Also may want to convert this at some point to C_NULL_CHAR
        #       delimeter instead of fixed width.  Future problem though.
        # Number of channel names may exceeed 5000
        self._channel_names_c = create_string_buffer(20 * 4000)
        self._channel_units_c = create_string_buffer(20 * 4000)

        # Initial environmental conditions
        self.gravity     = 9.80665  # Gravity (m/s^2)
        self.defWtrDens  = 1025.0   # Water density (kg/m^3)
        self.defWtrDpth  = 200.0    # Water depth (m)
        self.defMSL2SWL  = 0.0      # Offset between still-water level and mean sea level (m) [positive upward]

        # Interpolation order (must be 1: linear, or 2: quadratic)
        self.InterpOrder    =   1   # default of linear interpolation

        # Initial time related variables
        self.dt = 0.1                   # typical default for HD
        self.t_start = 0.0              # initial time 
        self.tmax = 600.0               # typical default for HD waves FFT
        #FIXME: check tmax/total_time and note exactly what is different between them.
        self.total_time = 0.0           # may be longer than tmax
        self.numTimeSteps = 0

        self.numChannels = 0                # Number of channels returned

        # Number of bodies and initial reference point
        #   The initial position is only set as (X,Y).  The Z value and
        #   orientation is set by HD and will be returned along with the full
        #   set of numBodies where it is expecting loads inputs.
        self.ptfmRefPt_x    = 0.0
        self.ptfmRefPt_y    = 0.0

        # Nodes
        #   The number of nodes must be constant throughout simulation.  The
        #   initial position is given in the initNodePos array (resize as
        #   needed, should be Nx6).
        #   Rotations are given in radians assuming small angles.  See note at
        #   top of this file.
        self.numNodePts     = 1     # Single ptfm attachment point for floating rigid
        self.initNodePos = np.zeros((self.numNodePts,6))      # N x 6 array [x,y,z,Rx,Ry,Rz]

        # OutRootName
        #   If HD writes a file (echo, summary, or other), use this for the
        #   root of the file name.
        self.outRootName = "Output_HDlib_default"

    # _initialize_routines() ------------------------------------------------------------------------------------------------------------
    def _initialize_routines(self):
        self.HydroDyn_C_Init.argtypes = [
            POINTER(c_char),                    # OutRootName 
            POINTER(c_char_p),                  # SeaState input file string
            POINTER(c_int),                     # SeaState input file string length
            POINTER(c_char_p),                  # HydroDyn input file string
            POINTER(c_int),                     # HydroDyn input file string length
            POINTER(c_float),                   # gravity
            POINTER(c_float),                   # defWtrDens
            POINTER(c_float),                   # defWtrDpth
            POINTER(c_float),                   # defMSL2SWL
            POINTER(c_float),                   # PtfmRefPt_X
            POINTER(c_float),                   # PtfmRefPt_Y
            POINTER(c_int),                     # numNodePts -- number of points expecting motions/loads
            POINTER(c_float),                   # initNodePos -- initial node positions in flat array of 6*numNodePts
            POINTER(c_int),                     # InterpOrder
            POINTER(c_double),                  # t_initial 
            POINTER(c_double),                  # dt
            POINTER(c_double),                  # tmax 
            POINTER(c_int),                     # number of channels
            POINTER(c_char),                    # output channel names
            POINTER(c_char),                    # output channel units
            POINTER(c_int),                     # ErrStat_C
            POINTER(c_char)                     # ErrMsg_C
        ]
        self.HydroDyn_C_Init.restype = c_int 

        self.HydroDyn_C_CalcOutput.argtypes = [
            POINTER(c_double),                  # Time_C
            POINTER(c_int),                     # numNodePts -- number of points expecting motions/loads
            POINTER(c_float),                   # nodePos -- node positions      in flat array of 6*numNodePts
            POINTER(c_float),                   # nodeVel -- node velocities     in flat array of 6*numNodePts
            POINTER(c_float),                   # nodeAcc -- node accelerations  in flat array of 6*numNodePts
            POINTER(c_float),                   # nodeFrc -- node forces/moments in flat array of 6*numNodePts
            POINTER(c_float),                   # Output Channel Values
            POINTER(c_int),                     # ErrStat_C
            POINTER(c_char)                     # ErrMsg_C
        ]
        self.HydroDyn_C_CalcOutput.restype = c_int

        self.HydroDyn_C_UpdateStates.argtypes = [
            POINTER(c_double),                  # Time_C
            POINTER(c_double),                  # TimeNext_C
            POINTER(c_int),                     # numNodePts -- number of points expecting motions/loads
            POINTER(c_float),                   # nodePos -- node positions      in flat array of 6*numNodePts
            POINTER(c_float),                   # nodeVel -- node velocities     in flat array of 6*numNodePts
            POINTER(c_float),                   # nodeAcc -- node accelerations  in flat array of 6*numNodePts
            POINTER(c_int),                     # ErrStat_C
            POINTER(c_char)                     # ErrMsg_C
        ]
        self.HydroDyn_C_UpdateStates.restype = c_int

        self.HydroDyn_C_End.argtypes = [
            POINTER(c_int),                     # ErrStat_C
            POINTER(c_char)                     # ErrMsg_C
        ]
        self.HydroDyn_C_End.restype = c_int

    # hydrodyn_init ------------------------------------------------------------------------------------------------------------
    def hydrodyn_init(self, seast_input_string_array, hd_input_string_array):
        # nodePositions -- N x 6 array  -- position info as [x1,y1,z1,Rx1,Ry1,Rz1]

        # Primary input file will be passed as a single string joined by
        # C_NULL_CHAR.
        seast_input_string = '\x00'.join(seast_input_string_array)
        seast_input_string = seast_input_string.encode('utf-8')
        seast_input_string_length = len(seast_input_string)
        hd_input_string = '\x00'.join(hd_input_string_array)
        hd_input_string = hd_input_string.encode('utf-8')
        hd_input_string_length = len(hd_input_string)
        
        self._numChannels_c = c_int(0)

        # Rootname for HD output files (echo etc).
        _outRootName_c = create_string_buffer((self.outRootName.ljust(self.default_str_c_len)).encode('utf-8'))

        # store initial number of node points for error handling at calls
        self._initNumNodePts = self.numNodePts

        # initNodePos
        #   Verify that the shape of initNodePos is correct
        if self.initNodePos.shape[1] != 6:
            print("Expecting a Nx6 array of initial node locations (initNodePos) with second index for [x,y,z,Rx,Ry,Rz]")
            self.hydrodyn_end()
            raise Exception("\nHydroDyn terminated prematurely.")
        if self.initNodePos.shape[0] != self.numNodePts:
            print("Expecting a Nx6 array of initial node locations (initNodePos) with first index for node number.")
            self.hydrodyn_end()
            raise Exception("\nHydroDyn terminated prematurely.")
        #FIXME: for now we only allow for a single rigid motion platform.  If
        #       multiple nodes are allowed for substructure coupling or
        #       flexible floating platforms, then remove this check and verify
        #       the fortran side works as expected (theoretically should, but
        #       it is untested).
        if self.numNodePts != 1:
            print(f"HydroDyn C interface does not currently support flexible structures ({self.numNodePts} input motion points were requested, but only 1 is currently supported).")
            self.hydrodyn_end()
            raise Exception("\nHydroDyn terminated prematurely.")

        #   Make a flat 1D array of position info:
        #       [x2,y1,z1,Rx1,Ry1,Rz1, x2,y2,z2,Rx2,Ry2,Rz2 ...]
        nodeInitLoc_flat = [pp for p in self.initNodePos for pp in p]
        nodeInitLoc_flat_c = (c_float * (6 * self.numNodePts))(0.0,)
        for i, p in enumerate(nodeInitLoc_flat):
            nodeInitLoc_flat_c[i] = c_float(p)


        # call HydroDyn_C_Init
        self.HydroDyn_C_Init(
            _outRootName_c,                         # IN: rootname for HD file writing
            c_char_p(seast_input_string),           # IN: SeaState input file string
            byref(c_int(seast_input_string_length)),# IN: SeaState input file string length
            c_char_p(hd_input_string),              # IN: HydroDyn input file string
            byref(c_int(hd_input_string_length)),   # IN: HydroDyn input file string length
            byref(c_float(self.gravity)),           # IN: gravity
            byref(c_float(self.defWtrDens)),        # IN: default water density
            byref(c_float(self.defWtrDpth)),        # IN: default water depth
            byref(c_float(self.defMSL2SWL)),        # IN: default offset between still-water level and mean sea level
            byref(c_float(self.ptfmRefPt_x)),       # IN: Platform initial position (X)
            byref(c_float(self.ptfmRefPt_y)),       # IN: Platform initial position (Y)
            byref(c_int(self.numNodePts)),          # IN: number of attachment points expected (where motions are transferred into HD)
            nodeInitLoc_flat_c,                     # IN: initNodePos -- initial node positions in flat array of 6*numNodePts
            byref(c_int(self.InterpOrder)),         # IN: InterpOrder (1: linear, 2: quadratic)
            byref(c_double(self.t_start)),          # IN: time initial 
            byref(c_double(self.dt)),               # IN: time step (dt)
            byref(c_double(self.tmax)),             # IN: tmax
            byref(self._numChannels_c),             # OUT: number of channels
            self._channel_names_c,                  # OUT: output channel names
            self._channel_units_c,                  # OUT: output channel units
            byref(self.error_status_c),             # OUT: ErrStat_C
            self.error_message_c                    # OUT: ErrMsg_C
        )

        self.check_error()
        
        # Initialize output channels
        self.numChannels = self._numChannels_c.value


    # hydrodyn_calcOutput ------------------------------------------------------------------------------------------------------------
    def hydrodyn_calcOutput(self, time, nodePos, nodeVel, nodeAcc, nodeFrcMom, outputChannelValues):

        # Check input motion info
        self.check_input_motions(nodePos,nodeVel,nodeAcc)

        # set flat arrays for inputs of motion
        #   Position -- [x2,y1,z1,Rx1,Ry1,Rz1, x2,y2,z2,Rx2,Ry2,Rz2 ...]
        nodePos_flat = [pp for p in nodePos for pp in p]
        nodePos_flat_c = (c_float * (6 * self.numNodePts))(0.0,)
        for i, p in enumerate(nodePos_flat):
            nodePos_flat_c[i] = c_float(p)

        #   Velocity -- [Vx2,Vy1,Vz1,RVx1,RVy1,RVz1, Vx2,Vy2,Vz2,RVx2,RVy2,RVz2 ...]
        nodeVel_flat = [pp for p in nodeVel for pp in p]
        nodeVel_flat_c = (c_float * (6 * self.numNodePts))(0.0,)
        for i, p in enumerate(nodeVel_flat):
            nodeVel_flat_c[i] = c_float(p)

        #   Acceleration -- [Ax1,Ay1,Az1,RAx1,RAy1,RAz1, Ax2,Ay2,Az2,RAx2,RAy2,RAz2 ...]
        nodeAcc_flat = [pp for p in nodeAcc for pp in p]
        nodeAcc_flat_c = (c_float * (6 * self.numNodePts))(0.0,)
        for i, p in enumerate(nodeAcc_flat):
            nodeAcc_flat_c[i] = c_float(p)

        # Resulting Forces/moments --  [Fx1,Fy1,Fz1,Mx1,My1,Mz1, Fx2,Fy2,Fz2,Mx2,My2,Mz2 ...]
        nodeFrc_flat_c = (c_float * (6 * self.numNodePts))(0.0,)

        # Set up output channels
        outputChannelValues_c = (c_float * self.numChannels)(0.0,)

        # Run HydroDyn_C_CalcOutput
        self.HydroDyn_C_CalcOutput(
            byref(c_double(time)),                  # IN: time at which to calculate output forces 
            byref(c_int(self.numNodePts)),          # IN: number of attachment points expected (where motions are transferred into HD)
            nodePos_flat_c,                         # IN: positions - specified by user
            nodeVel_flat_c,                         # IN: velocities at desired positions
            nodeAcc_flat_c,                         # IN: accelerations at desired positions
            nodeFrc_flat_c,                         # OUT: resulting forces/moments array
            outputChannelValues_c,                  # OUT: output channel values as described in input file
            byref(self.error_status_c),             # OUT: ErrStat_C
            self.error_message_c                    # OUT: ErrMsg_C
        )

        self.check_error()

        ## Reshape Force/Moment into [N,6]
        count = 0
        for j in range(0,self.numNodePts):
            nodeFrcMom[j,0] = nodeFrc_flat_c[count]
            nodeFrcMom[j,1] = nodeFrc_flat_c[count+1]
            nodeFrcMom[j,2] = nodeFrc_flat_c[count+2]
            nodeFrcMom[j,3] = nodeFrc_flat_c[count+3]
            nodeFrcMom[j,4] = nodeFrc_flat_c[count+4]
            nodeFrcMom[j,5] = nodeFrc_flat_c[count+5]
            count = count + 6
        
        # Convert output channel values back into python
        for k in range(0,self.numChannels):
            outputChannelValues[k] = float(outputChannelValues_c[k])

    # hydrodyn_updateStates ------------------------------------------------------------------------------------------------------------
    def hydrodyn_updateStates(self, time, timeNext, nodePos, nodeVel, nodeAcc, nodeFrcMom):

        # Check input motion info
        self.check_input_motions(nodePos,nodeVel,nodeAcc)

        # set flat arrays for inputs of motion
        #   Position -- [x2,y1,z1,Rx1,Ry1,Rz1, x2,y2,z2,Rx2,Ry2,Rz2 ...]
        nodePos_flat = [pp for p in nodePos for pp in p]
        nodePos_flat_c = (c_float * (6 * self.numNodePts))(0.0,)
        for i, p in enumerate(nodePos_flat):
            nodePos_flat_c[i] = c_float(p)

        #   Velocity -- [Vx2,Vy1,Vz1,RVx1,RVy1,RVz1, Vx2,Vy2,Vz2,RVx2,RVy2,RVz2 ...]
        nodeVel_flat = [pp for p in nodeVel for pp in p]
        nodeVel_flat_c = (c_float * (6 * self.numNodePts))(0.0,)
        for i, p in enumerate(nodeVel_flat):
            nodeVel_flat_c[i] = c_float(p)

        #   Acceleration -- [Ax1,Ay1,Az1,RAx1,RAy1,RAz1, Ax2,Ay2,Az2,RAx2,RAy2,RAz2 ...]
        nodeAcc_flat = [pp for p in nodeAcc for pp in p]
        nodeAcc_flat_c = (c_float * (6 * self.numNodePts))(0.0,)
        for i, p in enumerate(nodeAcc_flat):
            nodeAcc_flat_c[i] = c_float(p)

        # Resulting Forces/moments --  [Fx1,Fy1,Fz1,Mx1,My1,Mz1, Fx2,Fy2,Fz2,Mx2,My2,Mz2 ...]
        nodeFrc_flat_c = (c_float * (6 * self.numNodePts))(0.0,)

        # Run HydroDyn_UpdateStates_c
        self.HydroDyn_C_UpdateStates(
            byref(c_double(time)),                  # IN: time at which to calculate output forces 
            byref(c_double(timeNext)),              # IN: time T+dt we are stepping to 
            byref(c_int(self.numNodePts)),          # IN: number of attachment points expected (where motions are transferred into HD)
            nodePos_flat_c,                         # IN: positions - specified by user
            nodeVel_flat_c,                         # IN: velocities at desired positions
            nodeAcc_flat_c,                         # IN: accelerations at desired positions
            byref(self.error_status_c),             # OUT: ErrStat_C
            self.error_message_c                    # OUT: ErrMsg_C
        )

        self.check_error()

    # hydrodyn_end ------------------------------------------------------------------------------------------------------------
    def hydrodyn_end(self):
        if not self.ended:
            self.ended = True
            # Run HydroDyn_C_End
            self.HydroDyn_C_End(
                byref(self.error_status_c),
                self.error_message_c
            )

            self.check_error()

    # other functions ----------------------------------------------------------------------------------------------------------
    def check_error(self):
        if self.error_status_c.value == 0:
            return
        elif self.error_status_c.value < self.abort_error_level:
            print(f"HydroDyn error status: {self.error_levels[self.error_status_c.value]}: {self.error_message_c.value.decode('ascii')}")
        else:
            print(f"HydroDyn error status: {self.error_levels[self.error_status_c.value]}: {self.error_message_c.value.decode('ascii')}")
            self.hydrodyn_end()
            raise Exception("\nHydroDyn terminated prematurely.")


    def check_input_motions(self,nodePos,nodeVel,nodeAcc):
        # make sure number of nodes didn't change for some reason
        if self._initNumNodePts != self.numNodePts:
            # @ANDY TODO: `time` is not available here so this would be a runtime error
            print(f"At time {time}, the number of node points changed from initial value of {self._initNumNodePts}.  This is not permitted during the simulation.")
            self.hydrodyn_end()
            raise Exception("\nError in calling HydroDyn library.")

        #   Verify that the shape of positions array is correct
        if nodePos.shape[1] != 6:
            print("Expecting a Nx6 array of node positions (nodePos) with second index for [x,y,z,Rx,Ry,Rz]")
            self.hydrodyn_end()
            raise Exception("\nHydroDyn terminated prematurely.")
        if nodePos.shape[0] != self.numNodePts:
            print("Expecting a Nx6 array of node positions (nodePos) with first index for node number.")
            self.hydrodyn_end()
            raise Exception("\nHydroDyn terminated prematurely.")


        #   Verify that the shape of velocities array is correct
        if nodeVel.shape[1] != 6:
            print("Expecting a Nx6 array of node velocities (nodeVel) with second index for [x,y,z,Rx,Ry,Rz]")
            self.hydrodyn_end()
            raise Exception("\nHydroDyn terminated prematurely.")
        if nodeVel.shape[0] != self.numNodePts:
            print("Expecting a Nx6 array of node velocities (nodeVel) with first index for node number.")
            self.hydrodyn_end()
            raise Exception("\nHydroDyn terminated prematurely.")


        #   Verify that the shape of accelerations array is correct
        if nodeAcc.shape[1] != 6:
            print("Expecting a Nx6 array of node accelerations (nodeAcc) with second index for [x,y,z,Rx,Ry,Rz]")
            self.hydrodyn_end()
            raise Exception("\nHydroDyn terminated prematurely.")
        if nodeAcc.shape[0] != self.numNodePts:
            print("Expecting a Nx6 array of node accelerations (nodeAcc) with first index for node number.")
            self.hydrodyn_end()
            raise Exception("\nHydroDyn terminated prematurely.")



    @property
    def output_channel_names(self):
        if len(self._channel_names_c.value.split()) == 0:
            return []
        output_channel_names = self._channel_names_c.value.split()
        output_channel_names = [n.decode('UTF-8') for n in output_channel_names]
        return output_channel_names

    @property
    def output_channel_units(self):
        if len(self._channel_units_c.value.split()) == 0:
            return []
        output_channel_units = self._channel_units_c.value.split()
        output_channel_units = [n.decode('UTF-8') for n in output_channel_units]
        return output_channel_units


#===============================================================================
#   Helper class for debugging the interface.  This will write out all the
#   input position/orientation, velocities, accelerations, and the resulting
#   forces and moments at each input node.  If all is functioning correctly,
#   this will be identical to the corresponding values in the HydroDyn output
#   channels.

class DriverDbg():
    """
    This is only for debugging purposes only.  The input motions and resulting
    forces can be written to file with this class to verify the data I/O to the
    Fortran library.
    When coupled to another code, the force/moment array would be passed back
    to the calling code for use in the structural solver.
    """
    def __init__(self,filename,numNodePts):
        self.DbgFile=open(filename,'wt')        # open output file and write header info
        self.numNodePts=numNodePts
        # write file header
        t_string=datetime.datetime.now()
        dt_string=datetime.date.today()
        self.DbgFile.write(f"## This file was generated by hydrodyn_c_lib on {dt_string.strftime('%b-%d-%Y')} at {t_string.strftime('%H:%M:%S')}\n")
        self.DbgFile.write(f"## This file contains the resulting forces/moments at each of {self.numNodePts} node(s) passed into the hydrodyn_c_lib\n")
        self.DbgFile.write("#\n")
        f_string = "{:^25s}"
        self.DbgFile.write("#        T     ")
        for i in range(1,self.numNodePts+1):
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
        for i in range(1,self.numNodePts+1):
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

    def write(self,t,nodePos,nodeVel,nodeAcc,nodeFrc):
        t_string = "{:10.4f}"
        f_string = "{:25.7f}"*6
        self.DbgFile.write(t_string.format(t))
        for i in range(0,self.numNodePts):
            self.DbgFile.write(f_string.format(*nodePos[i,:]))
            self.DbgFile.write(f_string.format(*nodeVel[i,:]))
            self.DbgFile.write(f_string.format(*nodeAcc[i,:]))
            self.DbgFile.write(f_string.format(*nodeFrc[i,:]))
        self.DbgFile.write("\n")

    def end(self):
        if self.opened:
            self.DbgFile.close()
            self.opened = False


#===============================================================================
#   Helper class for writing channels to file.
#   for the regression testing to mirror the output from the InfowWind Fortran
#   driver.  This may also have value for debugging the interfacing to IfW.

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
        self.OutFile.write(f"## This file was generated by HydroDyn_Driver on {dt_string.strftime('%b-%d-%Y')} at {t_string.strftime('%H:%M:%S')}\n")
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
