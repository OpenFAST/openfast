#**********************************************************************************************************************************
# LICENSING
# Copyright (C) 2021 National Renewable Energy Laboratory
#
# This file is part of InflowWind.
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
# This is the Python-C interface library for AeroDyn with InflowWind.  This may
# be used directly with Python based codes to call and run AeroDyn and
# InflowWind together.  An example of using this library from Python is given
# in the accompanying Python driver program.  Additional notes and information
# on the interfacing is included there.
#
#
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
    c_wchar,
    c_wchar_p,
    c_bool
)
import numpy as np
import datetime

class AeroDynInflowLib(CDLL):
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

        # Input file handling
        self.ADinputPass  = True            # Assume passing of input file as a string
        self.IfWinputPass = True            # Assume passing of input file as a string

        # Create buffers for class data
        self.abort_error_level = 4
        self.error_status_c = c_int(0)
        self.error_message_c = create_string_buffer(self.error_msg_c_len)

        # This is not sufficient for AD
        #FIXME: ChanLen may not always be 20 -- could be as much as 256
        #       Possible fix is to pass this length over to Fortran side.
        #       Also may want to convert this at some point to C_NULL_CHAR
        #       delimeter instead of fixed width.  Future problem though.
        # Number of channel names may exceeed 5000
        self._channel_names_c = create_string_buffer(20 * 8000)
        self._channel_units_c = create_string_buffer(20 * 8000)

        # Initial environmental conditions
        #self.MHK = false    #  MHK turbine type switch -- disabled for now
        self.gravity     =   9.80665  # Gravitational acceleration (m/s^2)
        self.defFldDens  =     1.225  # Air density (kg/m^3)
        self.defKinVisc  = 1.464E-05  # Kinematic viscosity of working fluid (m^2/s)
        self.defSpdSound =     335.0  # Speed of sound in working fluid (m/s)
        self.defPatm     =  103500.0  # Atmospheric pressure (Pa) [used only for an MHK turbine cavitation check]
        self.defPvap     =    1700.0  # Vapour pressure of working fluid (Pa) [used only for an MHK turbine cavitation check]
        self.WtrDpth     =       0.0  # Water depth (m)
        self.MSL2SWL     =       0.0  # Offset between still-water level and mean sea level (m) [positive upward]

        # flags 
        self.storeHHVel  = False
        self.transposeDCM= False

        # VTK
        self.WrVTK       = 0          # default of no vtk output
        self.WrVTK_Type  = 1          # default of surface meshes
        self.VTKNacDim   = np.array([-2.5,-2.5,0,10,5,5], dtype="float32")        # default nacelle dimension for VTK surface rendering [x0,y0,z0,Lx,Ly,Lz] (m)
        self.VTKHubRad   = 1.5        # default hub radius for VTK surface rendering

        # Output file
        self.wrOuts      = 0          # wrOuts -- file format for writing outputs
        self.DT_Outs     = 0.0        # DT_Outs -- timestep for outputs to file

        # Interpolation order (must be 1: linear, or 2: quadratic)
        self.InterpOrder = 1          # default of linear interpolation

        # Initial time related variables
        self.dt          = 0.1        # typical default for HD
        self.tmax        = 600.0      # typical default for HD waves FFT
        #FIXME: check tmax/total_time and note exactly what is different between them.
        self.total_time  = 0.0        # may be longer than tmax
        self.numTimeSteps= 0

        # number of output channels
        self.numChannels = 0          # Number of channels returned

        # Aero calculation method -- AeroProjMod
        #     APM_BEM_NoSweepPitchTwist - 1 -  "Original AeroDyn model where momentum balance is done in the WithoutSweepPitchTwist system"
        #     APM_BEM_Polar             - 2 -  "Use staggered polar grid for momentum balance in each annulus"
        #     APM_LiftingLine           - 3 -  "Use the blade lifting line (i.e. the structural) orientation (currently for OLAF with VAWT)"
        self.AeroProjMod = 1

        # Initial position of hub and blades
        #   used for setup of AD, not used after init.
        self.initHubPos         = np.zeros(shape=(3),dtype=c_float)
        self.initHubOrient      = np.zeros(shape=(9),dtype=c_double)
        self.initNacellePos     = np.zeros(shape=(3),dtype=c_float)
        self.initNacelleOrient  = np.zeros(shape=(9),dtype=c_double)
        self.numBlades          = 3
        self.initRootPos        = np.zeros(shape=(self.numBlades,3),dtype=c_float)
        self.initRootOrient     = np.zeros(shape=(self.numBlades,9),dtype=c_double)

        # Structural Mesh
        #   The number of nodes must be constant throughout simulation.  The
        #   initial position is given in the initMeshPos array (resize as
        #   needed, should be Nx6).
        #   Rotations are given in radians assuming small angles.  See note at
        #   top of this file.
        self.numMeshPts     = 1
        self.initMeshPos    = np.zeros(shape=(self.numMeshPts,3),dtype=c_float )    # Nx3 array [x,y,z]
        self.initMeshOrient = np.zeros(shape=(self.numMeshPts,9),dtype=c_double)    # Nx9 array [r11,r12,r13,r21,r22,r23,r31,r32,r33]

        # OutRootName
        #   If HD writes a file (echo, summary, or other), use this for the
        #   root of the file name.
        self.outRootName = "Output_ADIlib_default"

    # _initialize_routines() ------------------------------------------------------------------------------------------------------------
    def _initialize_routines(self):
        self.AeroDyn_Inflow_C_Init.argtypes = [
            POINTER(c_bool),                    # AD input file passed as string
            POINTER(c_char_p),                  # AD input file as string
            POINTER(c_int),                     # AD input file string length
            POINTER(c_bool),                    # IfW input file passed as string
            POINTER(c_char_p),                  # IfW input file as string
            POINTER(c_int),                     # IfW input file string length
            POINTER(c_char),                    # OutRootName 
            POINTER(c_float),                   # gravity
            POINTER(c_float),                   # defFldDens
            POINTER(c_float),                   # defKinVisc
            POINTER(c_float),                   # defSpdSound
            POINTER(c_float),                   # defPatm
            POINTER(c_float),                   # defPvap
            POINTER(c_float),                   # WtrDpth
            POINTER(c_float),                   # MSL2SWL
            POINTER(c_int),                     # AeroProjMod
            POINTER(c_int),                     # InterpOrder 
            POINTER(c_double),                  # dt
            POINTER(c_double),                  # tmax 
            POINTER(c_bool),                    # storeHHVel
            POINTER(c_bool),                    # transposeDCM
            POINTER(c_int),                     # WrVTK
            POINTER(c_int),                     # WrVTK_Type
            POINTER(c_float),                   # VTKNacDim
            POINTER(c_float),                   # VTKHubRad
            POINTER(c_int),                     # wrOuts -- file format for writing outputs
            POINTER(c_double),                  # DT_Outs -- timestep for outputs to file
            POINTER(c_float),                   # initHubPos
            POINTER(c_double),                  # initHubOrient_flat
            POINTER(c_float),                   # initNacellePos
            POINTER(c_double),                  # initNacelleOrient_flat
            POINTER(c_int),                     # numBlades
            POINTER(c_float),                   # initRootPos_flat
            POINTER(c_double),                  # initRootOrient_flat
            POINTER(c_int),                     # numMeshPts
            POINTER(c_float),                   # initMeshPos_flat
            POINTER(c_double),                  # initMeshOrient_flat
            POINTER(c_int),                     # number of channels
            POINTER(c_char),                    # output channel names
            POINTER(c_char),                    # output channel units
            POINTER(c_int),                     # ErrStat_C
            POINTER(c_char)                     # ErrMsg_C
        ]
        self.AeroDyn_Inflow_C_Init.restype = c_int 

        #self.AeroDyn_Inflow_C_ReInit.argtypes = [
        #    POINTER(c_double),                  # t_initial 
        #    POINTER(c_double),                  # dt
        #    POINTER(c_double),                  # tmax 
        #    POINTER(c_int),                     # ErrStat_C
        #    POINTER(c_char)                     # ErrMsg_C
        #]
        #self.AeroDyn_Inflow_C_ReInit.restype = c_int 

        self.AeroDyn_Inflow_C_CalcOutput.argtypes = [
            POINTER(c_double),                  # Time_C
            POINTER(c_float),                   # HubPos
            POINTER(c_double),                  # HubOrient_flat
            POINTER(c_float),                   # HubVel
            POINTER(c_float),                   # HubAcc
            POINTER(c_float),                   # NacPos
            POINTER(c_double),                  # NacOrient_flat
            POINTER(c_float),                   # NacVel
            POINTER(c_float),                   # NacAcc
            POINTER(c_float),                   # RootPos
            POINTER(c_double),                  # RootOrient_flat
            POINTER(c_float),                   # RootVel
            POINTER(c_float),                   # RootAcc
            POINTER(c_int),                     # numMeshPts
            POINTER(c_float),                   # MeshPos
            POINTER(c_double),                  # MeshOrient_flat
            POINTER(c_float),                   # MeshVel
            POINTER(c_float),                   # MeshAcc
            POINTER(c_float),                   # meshFrc -- mesh forces/moments in flat array of 6*numMeshPts
            POINTER(c_float),                   # Output Channel Values
            POINTER(c_int),                     # ErrStat_C
            POINTER(c_char)                     # ErrMsg_C
        ]
        self.AeroDyn_Inflow_C_CalcOutput.restype = c_int

        self.AeroDyn_Inflow_C_UpdateStates.argtypes = [
            POINTER(c_double),                  # Time_C
            POINTER(c_double),                  # TimeNext_C
            POINTER(c_float),                   # HubPos
            POINTER(c_double),                  # HubOrient_flat
            POINTER(c_float),                   # HubVel
            POINTER(c_float),                   # HubAcc
            POINTER(c_float),                   # NacPos
            POINTER(c_double),                  # NacOrient_flat
            POINTER(c_float),                   # NacVel
            POINTER(c_float),                   # NacAcc
            POINTER(c_float),                   # RootPos
            POINTER(c_double),                  # RootOrient_flat
            POINTER(c_float),                   # RootVel
            POINTER(c_float),                   # RootAcc
            POINTER(c_int),                     # numMeshPts
            POINTER(c_float),                   # MeshPos
            POINTER(c_double),                  # MeshOrient_flat
            POINTER(c_float),                   # MeshVel
            POINTER(c_float),                   # MeshAcc
            POINTER(c_int),                     # ErrStat_C
            POINTER(c_char)                     # ErrMsg_C
        ]
        self.AeroDyn_Inflow_C_UpdateStates.restype = c_int

        self.AeroDyn_Inflow_C_End.argtypes = [
            POINTER(c_int),                     # ErrStat_C
            POINTER(c_char)                     # ErrMsg_C
        ]
        self.AeroDyn_Inflow_C_End.restype = c_int

    # aerodyn_inflow_init ------------------------------------------------------------------------------------------------------------
    def aerodyn_inflow_init(self, AD_input_string_array, IfW_input_string_array):
        # some bookkeeping initialization
        self._numChannels_c = c_int(0)
        self._initNumMeshPts = self.numMeshPts
        self._initNumBlades  = self.numBlades

        # Primary input file will be passed as a single string joined by
        # C_NULL_CHAR.
        AD_input_string = '\x00'.join(AD_input_string_array)
        AD_input_string = AD_input_string.encode('utf-8')
        AD_input_string_length = len(AD_input_string)

        # Primary input file will be passed as a single string joined by
        # C_NULL_CHAR.
        IfW_input_string = '\x00'.join(IfW_input_string_array)
        IfW_input_string = IfW_input_string.encode('utf-8')
        IfW_input_string_length = len(IfW_input_string)

        # Rootname for ADI output files (echo etc).
        _outRootName_c = create_string_buffer((self.outRootName.ljust(self.default_str_c_len)).encode('utf-8'))

        # check hub and root points for initialization
        self.check_init_hubroot()

        # Check initial mesh positions
        self.check_init_mesh()

        #   Flatten arrays to pass
        #       [x2,y1,z1, x2,y2,z2 ...]
        VTKNacDim_c             = (c_float  * len(self.VTKNacDim        ))(*self.VTKNacDim        )
        initHubPos_c            = (c_float  * len(self.initHubPos       ))(*self.initHubPos       )
        initHubOrient_c         = (c_double * len(self.initHubOrient    ))(*self.initHubOrient    )
        initNacellePos_c        = (c_float  * len(self.initNacellePos   ))(*self.initNacellePos   )
        initNacelleOrient_c     = (c_double * len(self.initNacelleOrient))(*self.initNacelleOrient)
        initRootPos_flat_c      = self.flatPosArr(   self._initNumBlades, self.numBlades,self.initRootPos,    'Init','RootPos')
        initRootOrient_flat_c   = self.flatOrientArr(self._initNumBlades, self.numBlades,self.initRootOrient, 'Init','RootOrient')
        initMeshPos_flat_c      = self.flatPosArr(   self._initNumMeshPts,self.numMeshPts,self.initMeshPos,   'Init','MeshPos')
        initMeshOrient_flat_c   = self.flatOrientArr(self._initNumMeshPts,self.numMeshPts,self.initMeshOrient,'Init','MeshOrient')


        # call AeroDyn_Inflow_C_Init
        self.AeroDyn_Inflow_C_Init(
            byref(c_bool(self.ADinputPass)),        # IN: AD input file is passed
            c_char_p(AD_input_string),              # IN: AD input file as string (or filename if ADinputPass is false)
            byref(c_int(AD_input_string_length)),   # IN: AD input file string length
            byref(c_bool(self.IfWinputPass)),       # IN: IfW input file is passed
            c_char_p(IfW_input_string),             # IN: IfW input file as string (or filename if IfWinputPass is false)
            byref(c_int(IfW_input_string_length)),  # IN: IfW input file string length
            _outRootName_c,                         # IN: rootname for ADI file writing
            byref(c_float(self.gravity)),           # IN: gravity
            byref(c_float(self.defFldDens)),        # IN: defFldDens
            byref(c_float(self.defKinVisc)),        # IN: defKinVisc
            byref(c_float(self.defSpdSound)),       # IN: defSpdSound
            byref(c_float(self.defPatm)),           # IN: defPatm
            byref(c_float(self.defPvap)),           # IN: defPvap
            byref(c_float(self.WtrDpth)),           # IN: WtrDpth
            byref(c_float(self.MSL2SWL)),           # IN: MSL2SWL
            byref(c_int(self.AeroProjMod)),         # IN: AeroProjMod
            byref(c_int(self.InterpOrder)),         # IN: InterpOrder (1: linear, 2: quadratic)
            byref(c_double(self.dt)),               # IN: time step (dt)
            byref(c_double(self.tmax)),             # IN: tmax
            byref(c_bool(self.storeHHVel)),         # IN: storeHHVel
            byref(c_bool(self.transposeDCM)),       # IN: transposeDCM
            byref(c_int(self.WrVTK)),               # IN: WrVTK
            byref(c_int(self.WrVTK_Type)),          # IN: WrVTK_Type
            VTKNacDim_c,                            # IN: VTKNacDim
            byref(c_float(self.VTKHubRad)),         # IN: VTKHubRad
            byref(c_int(self.wrOuts)),              # IN: wrOuts -- file format for writing outputs
            byref(c_double(self.DT_Outs)),          # IN: DT_Outs -- timestep for outputs to file
            initHubPos_c,                           # IN: initHubPos -- initial hub position
            initHubOrient_c,                        # IN: initHubOrient -- initial hub orientation DCM in flat array of 9 elements
            initNacellePos_c,                       # IN: initNacellePos -- initial hub position
            initNacelleOrient_c,                    # IN: initNacelleOrient -- initial hub orientation DCM in flat array of 9 elements
            byref(c_int(self.numBlades)),           # IN: number of blades (matches number of blade root positions)
            initRootPos_flat_c,                     # IN: initBladeRootPos -- initial node positions in flat array of 3*numBlades
            initRootOrient_flat_c,                  # IN: initBladeRootOrient -- initial blade root orientation DCMs in flat array of 9*numBlades
            byref(c_int(self.numMeshPts)),          # IN: number of attachment points expected (where motions are transferred into HD)
            initMeshPos_flat_c,                     # IN: initMeshPos -- initial node positions in flat array of 3*numMeshPts
            initMeshOrient_flat_c,                  # IN: initMeshOrient -- initial node orientation DCMs in flat array of 9*numMeshPts
            byref(self._numChannels_c),             # OUT: number of channels
            self._channel_names_c,                  # OUT: output channel names
            self._channel_units_c,                  # OUT: output channel units
            byref(self.error_status_c),             # OUT: ErrStat_C
            self.error_message_c                    # OUT: ErrMsg_C
        )

        self.check_error()
        
        # Initialize output channels
        self.numChannels = self._numChannels_c.value


    ## aerodyn_inflow_reinit ------------------------------------------------------------------------------------------------------------
    #def aerodyn_inflow_reinit(self):
    #    #FIXME: need to pass something in here I think.  Not sure what.
    #
    #    # call AeroDyn_Inflow_C_ReInit
    #    self.AeroDyn_Inflow_C_ReInit(
    #        byref(c_double(self.dt)),               # IN: time step (dt)
    #        byref(c_double(self.tmax)),             # IN: tmax
    #        byref(self.error_status_c),             # OUT: ErrStat_C
    #        self.error_message_c                    # OUT: ErrMsg_C
    #    )
    #
    #    self.check_error()
    #    #FIXME: anything coming out that needs handling/passing?


    # aerodyn_inflow_calcOutput ------------------------------------------------------------------------------------------------------------
    def aerodyn_inflow_calcOutput(self, time, hubPos, hubOrient, hubVel, hubAcc, \
                            nacPos, nacOrient, nacVel, nacAcc, \
                            rootPos, rootOrient, rootVel, rootAcc, \
                            meshPos, meshOrient, meshVel, meshAcc, \
                            meshFrcMom, outputChannelValues):

        # Check input motion info
        self.check_input_motions_hubNac(hubPos,hubOrient,hubVel,hubAcc,'hub')
        self.check_input_motions_hubNac(nacPos,nacOrient,nacVel,nacAcc,'nacelle')
        self.check_input_motions_root(rootPos,rootOrient,rootVel,rootAcc)
        self.check_input_motions_mesh(meshPos,meshOrient,meshVel,meshAcc)

        _hubPos_c    = (c_float  * len(np.squeeze(hubPos)   ))(*np.squeeze(hubPos)   )
        _hubOrient_c = (c_double * len(np.squeeze(hubOrient)))(*np.squeeze(hubOrient))
        _hubVel_c    = (c_float  * len(np.squeeze(hubVel)   ))(*np.squeeze(hubVel)   )
        _hubAcc_c    = (c_float  * len(np.squeeze(hubAcc)   ))(*np.squeeze(hubAcc)   )
        _nacPos_c    = (c_float  * len(np.squeeze(nacPos)   ))(*np.squeeze(nacPos)   )
        _nacOrient_c = (c_double * len(np.squeeze(nacOrient)))(*np.squeeze(nacOrient))
        _nacVel_c    = (c_float  * len(np.squeeze(nacVel)   ))(*np.squeeze(nacVel)   )
        _nacAcc_c    = (c_float  * len(np.squeeze(nacAcc)   ))(*np.squeeze(nacAcc)   )
        #   Make a flat 1D arrays of motion info:
        #       [x2,y1,z1, x2,y2,z2 ...]
        _rootPos_flat_c    = self.flatPosArr(   self._initNumBlades,self.numBlades,rootPos,   time,'MeshPos')
        _rootOrient_flat_c = self.flatOrientArr(self._initNumBlades,self.numBlades,rootOrient,time,'MeshOrient')
        _rootVel_flat_c    = self.flatVelAccArr(self._initNumBlades,self.numBlades,rootVel,   time,'MeshVel')
        _rootAcc_flat_c    = self.flatVelAccArr(self._initNumBlades,self.numBlades,rootAcc,   time,'MeshAcc')
        #   Make a flat 1D arrays of motion info:
        #       [x2,y1,z1, x2,y2,z2 ...]
        _meshPos_flat_c    = self.flatPosArr(   self._initNumMeshPts,self.numMeshPts,meshPos,   time,'MeshPos')
        _meshOrient_flat_c = self.flatOrientArr(self._initNumMeshPts,self.numMeshPts,meshOrient,time,'MeshOrient')
        _meshVel_flat_c    = self.flatVelAccArr(self._initNumMeshPts,self.numMeshPts,meshVel,   time,'MeshVel')
        _meshAcc_flat_c    = self.flatVelAccArr(self._initNumMeshPts,self.numMeshPts,meshAcc,   time,'MeshAcc')

        # Resulting Forces/moments --  [Fx1,Fy1,Fz1,Mx1,My1,Mz1, Fx2,Fy2,Fz2,Mx2,My2,Mz2 ...]
        _meshFrc_flat_c = (c_float * (6 * self.numMeshPts))(0.0,)

        # Set up output channels
        outputChannelValues_c = (c_float * self.numChannels)(0.0,)

        # Run AeroDyn_Inflow_C_CalcOutput
        self.AeroDyn_Inflow_C_CalcOutput(
            byref(c_double(time)),                  # IN: time at which to calculate output forces 
            _hubPos_c,                              # IN: hub positions
            _hubOrient_c,                           # IN: hub orientations
            _hubVel_c,                              # IN: hub velocity [TVx,TVy,TVz,RVx,RVy,RVz]
            _hubAcc_c,                              # IN: hub acclerations [TAx,TAy,TAz,RAx,RAy,RAz]
            _nacPos_c,                              # IN: nac positions
            _nacOrient_c,                           # IN: nac orientations
            _nacVel_c,                              # IN: nac velocity [TVx,TVy,TVz,RVx,RVy,RVz]
            _nacAcc_c,                              # IN: nac acclerations [TAx,TAy,TAz,RAx,RAy,RAz]
            _rootPos_flat_c,                        # IN: positions
            _rootOrient_flat_c,                     # IN: Orientations (DCM)
            _rootVel_flat_c,                        # IN: velocities at desired positions
            _rootAcc_flat_c,                        # IN: accelerations at desired positions
            byref(c_int(self.numMeshPts)),          # IN: number of attachment points expected (where motions are transferred into HD)
            _meshPos_flat_c,                        # IN: positions
            _meshOrient_flat_c,                     # IN: Orientations (DCM)
            _meshVel_flat_c,                        # IN: velocities at desired positions
            _meshAcc_flat_c,                        # IN: accelerations at desired positions
            _meshFrc_flat_c,                        # OUT: resulting forces/moments array
            outputChannelValues_c,                  # OUT: output channel values as described in input file
            byref(self.error_status_c),             # OUT: ErrStat_C
            self.error_message_c                    # OUT: ErrMsg_C
        )

        self.check_error()

        ## Reshape Force/Moment into [N,6]
        count = 0
        for j in range(0,self.numMeshPts):
            meshFrcMom[j,0] = _meshFrc_flat_c[count]
            meshFrcMom[j,1] = _meshFrc_flat_c[count+1]
            meshFrcMom[j,2] = _meshFrc_flat_c[count+2]
            meshFrcMom[j,3] = _meshFrc_flat_c[count+3]
            meshFrcMom[j,4] = _meshFrc_flat_c[count+4]
            meshFrcMom[j,5] = _meshFrc_flat_c[count+5]
            count = count + 6
        
        # Convert output channel values back into python
        for k in range(0,self.numChannels):
            outputChannelValues[k] = float(outputChannelValues_c[k])

    # aerodyn_inflow_updateStates ------------------------------------------------------------------------------------------------------------
    def aerodyn_inflow_updateStates(self, time, timeNext, \
                            hubPos, hubOrient, hubVel, hubAcc, \
                            nacPos, nacOrient, nacVel, nacAcc, \
                            rootPos, rootOrient, rootVel, rootAcc, \
                            meshPos, meshOrient, meshVel, meshAcc):

        # Check input motion info
        self.check_input_motions_hubNac(hubPos,hubOrient,hubVel,hubAcc,'hub')
        self.check_input_motions_hubNac(nacPos,nacOrient,nacVel,nacAcc,'nacelle')
        self.check_input_motions_root(rootPos,rootOrient,rootVel,rootAcc)
        self.check_input_motions_mesh(meshPos,meshOrient,meshVel,meshAcc)

        _hubPos_c    = (c_float  * len(np.squeeze(hubPos)   ))(*np.squeeze(hubPos)   )
        _hubOrient_c = (c_double * len(np.squeeze(hubOrient)))(*np.squeeze(hubOrient))
        _hubVel_c    = (c_float  * len(np.squeeze(hubVel)   ))(*np.squeeze(hubVel)   )
        _hubAcc_c    = (c_float  * len(np.squeeze(hubAcc)   ))(*np.squeeze(hubAcc)   )
        _nacPos_c    = (c_float  * len(np.squeeze(nacPos)   ))(*np.squeeze(nacPos)   )
        _nacOrient_c = (c_double * len(np.squeeze(nacOrient)))(*np.squeeze(nacOrient))
        _nacVel_c    = (c_float  * len(np.squeeze(nacVel)   ))(*np.squeeze(nacVel)   )
        _nacAcc_c    = (c_float  * len(np.squeeze(nacAcc)   ))(*np.squeeze(nacAcc)   )
        #   Make a flat 1D arrays of motion info:
        #       [x2,y1,z1, x2,y2,z2 ...]
        _rootPos_flat_c    = self.flatPosArr(   self._initNumBlades,self.numBlades,rootPos,   time,'MeshPos')
        _rootOrient_flat_c = self.flatOrientArr(self._initNumBlades,self.numBlades,rootOrient,time,'MeshOrient')
        _rootVel_flat_c    = self.flatVelAccArr(self._initNumBlades,self.numBlades,rootVel,   time,'MeshVel')
        _rootAcc_flat_c    = self.flatVelAccArr(self._initNumBlades,self.numBlades,rootAcc,   time,'MeshAcc')
        #   Make a flat 1D arrays of motion info:
        #       [x2,y1,z1, x2,y2,z2 ...]
        _meshPos_flat_c    = self.flatPosArr(   self._initNumMeshPts,self.numMeshPts,meshPos,   time,'MeshPos')
        _meshOrient_flat_c = self.flatOrientArr(self._initNumMeshPts,self.numMeshPts,meshOrient,time,'MeshOrient')
        _meshVel_flat_c    = self.flatVelAccArr(self._initNumMeshPts,self.numMeshPts,meshVel,   time,'MeshVel')
        _meshAcc_flat_c    = self.flatVelAccArr(self._initNumMeshPts,self.numMeshPts,meshAcc,   time,'MeshAcc')

        # Resulting Forces/moments --  [Fx1,Fy1,Fz1,Mx1,My1,Mz1, Fx2,Fy2,Fz2,Mx2,My2,Mz2 ...]
        _meshFrc_flat_c = (c_float * (6 * self.numMeshPts))(0.0,)

        # Run AeroDyn_Inflow_UpdateStates_c
        self.AeroDyn_Inflow_C_UpdateStates(
            byref(c_double(time)),                  # IN: time at which to calculate output forces 
            byref(c_double(timeNext)),              # IN: time T+dt we are stepping to
            _hubPos_c,                              # IN: hub positions
            _hubOrient_c,                           # IN: hub orientations
            _hubVel_c,                              # IN: hub velocity [TVx,TVy,TVz,RVx,RVy,RVz]
            _hubAcc_c,                              # IN: hub acclerations [TAx,TAy,TAz,RAx,RAy,RAz]
            _nacPos_c,                              # IN: nac positions
            _nacOrient_c,                           # IN: nac orientations
            _nacVel_c,                              # IN: nac velocity [TVx,TVy,TVz,RVx,RVy,RVz]
            _nacAcc_c,                              # IN: nac acclerations [TAx,TAy,TAz,RAx,RAy,RAz]
            _rootPos_flat_c,                        # IN: positions
            _rootOrient_flat_c,                     # IN: Orientations (DCM)
            _rootVel_flat_c,                        # IN: velocities at desired positions
            _rootAcc_flat_c,                        # IN: accelerations at desired positions
            byref(c_int(self.numMeshPts)),          # IN: number of attachment points expected (where motions are transferred into HD)
            _meshPos_flat_c,                        # IN: positions
            _meshOrient_flat_c,                     # IN: Orientations (DCM)
            _meshVel_flat_c,                        # IN: velocities at desired positions
            _meshAcc_flat_c,                        # IN: accelerations at desired positions
            byref(self.error_status_c),             # OUT: ErrStat_C
            self.error_message_c                    # OUT: ErrMsg_C
        )

        self.check_error()

    # aerodyn_inflow_end ------------------------------------------------------------------------------------------------------------
    def aerodyn_inflow_end(self):
        if not self.ended:
            self.ended = True
            # Run AeroDyn_Inflow_C_End
            self.AeroDyn_Inflow_C_End(
                byref(self.error_status_c),
                self.error_message_c
            )

            self.check_error()

    # other functions ----------------------------------------------------------------------------------------------------------
    def check_error(self):
        if self.error_status_c.value == 0:
            return
        elif self.error_status_c.value < self.abort_error_level:
            print(f"AeroDyn/InflowWind error status: {self.error_levels[self.error_status_c.value]}: {self.error_message_c.value.decode('ascii')}")
        else:
            print(f"AeroDyn/InflowWind error status: {self.error_levels[self.error_status_c.value]}: {self.error_message_c.value.decode('ascii')}")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn/InflowWind terminated prematurely.")


    def flatPosArr(self,initNumMeshPts,numPts,MeshPosArr,time,name):
        if initNumMeshPts != numPts:
            print(f"At time {time}, the number of {name} points changed from initial value of {initNumMeshPts}.  This is not permitted during the simulation.")
            self.aerodyn_inflow_end()
            raise Exception("\nError in calling AeroDyn/InflowWind library.")
        meshPos_flat = [pp for p in MeshPosArr for pp in p]
        meshPos_flat_c = (c_float * (3 * numPts))(0.0,)
        for i, p in enumerate(meshPos_flat):
            meshPos_flat_c[i] = c_float(p)
        return meshPos_flat_c


    def flatOrientArr(self,initNumMeshPts,numPts,MeshOrientArr,time,name):
        if initNumMeshPts != numPts:
            print(f"At time {time}, the number of {name} points changed from initial value of {initNumMeshPts}.  This is not permitted during the simulation.")
            self.aerodyn_inflow_end()
            raise Exception("\nError in calling AeroDyn/InflowWind library.")
        meshOrient_flat = [pp for p in MeshOrientArr for pp in p]
        meshOrient_flat_c = (c_double * (9 * numPts))(0.0,)
        for i, p in enumerate(meshOrient_flat):
            meshOrient_flat_c[i] = c_double(p)
        return meshOrient_flat_c


    def flatVelAccArr(self,initNumMeshPts,numPts,MeshArr,time,name):
        if initNumMeshPts != numPts:
            print(f"At time {time}, the number of {name} points changed from initial value of {initNumMeshPts}.  This is not permitted during the simulation.")
            self.aerodyn_inflow_end()
            raise Exception("\nError in calling AeroDyn/InflowWind library.")
        #   Velocity -- [Vx2,Vy1,Vz1,RVx1,RVy1,RVz1, Vx2,Vy2,Vz2,RVx2,RVy2,RVz2 ...]
        meshVel_flat = [pp for p in MeshArr for pp in p]
        meshVel_flat_c = (c_float * (6 * self.numMeshPts))(0.0,)
        for i, p in enumerate(meshVel_flat):
            meshVel_flat_c[i] = c_float(p)
        return meshVel_flat_c


    def check_init_hubroot(self):
        #print("shape of initRootPos       ",   self.initRootPos.shape)
        #print("               ndim        ",   np.squeeze(self.initRootPos.ndim))
        #print("               size 0      ",   self.initRootPos.shape[0])
        #print("               size 1      ",   self.initRootPos.shape[1])
        #print("shape of initRootOrient    ",   self.initRootOrient.shape)
        #print("               ndim        ",   np.squeeze(self.initRootPos.ndim))
        #print("               size 0      ",   self.initRootOrient.shape[0])
        #print("               size 1      ",   self.initRootOrient.shape[1])
        #print("               float       ",   type(self.initRootOrient[0,0]))
        #print("shape of initHubPos        ",   self.initHubPos.shape)
        #print("               ndim        ",   np.squeeze(self.initHubPos.ndim))
        #print("               size 0      ",   self.initHubPos.shape[0])
        #print("shape of initHubOrient     ",   self.initHubOrient.shape)
        #print("               ndim        ",   np.squeeze(self.initHubOrient.ndim))
        #print("               size 0      ",   self.initHubOrient.shape[0])
        #print("shape of initNacellePos    ",   self.initNacellePos.shape)
        #print("               ndim        ",   np.squeeze(self.initNacellePos.ndim))
        #print("               size 0      ",   self.initNacellePos.shape[0])
        #print("shape of initNacelleOrient ",   self.initNacelleOrient.shape)
        #print("               ndim        ",   np.squeeze(self.initNacelleOrient.ndim))
        #print("               size 0      ",   self.initNacelleOrient.shape[0])
        if self.numBlades < 1:
            print("No blades.  Set numBlades to number of AD blades in the model")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn terminated prematurely.")
        if self.initRootPos.shape[1] != 3:
            print("Expecting a Nx3 array of blade root positions (initRootPos) with second index for [x,y,z]")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn terminated prematurely.")
        if self.initRootPos.shape[0] != self.numBlades:
            print("Expecting a Nx3 array of blade root positions (initRootPos) with first index for blade number")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn terminated prematurely.")
        if self.initRootOrient.shape[1] != 9:
            print("Expecting a Nx9 array of blade root orientations as DCMs (initRootOrient) with second index for [r11,r12,r13,r21,r22,r23,r31,r32,r33]")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn terminated prematurely.")
        if self.initRootOrient.shape[0] != self.numBlades:
            print("Expecting a Nx3 array of blade root orientations (initRootOrient) with first index for blade number")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn terminated prematurely.")
        if np.squeeze(self.initHubPos.ndim) > 1 or self.initHubPos.shape[0] != 3:
            print("Expecting a 3 element array for initHubPos [x,y,z]")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn terminated prematurely.")
        if np.squeeze(self.initHubOrient.ndim) > 1 or self.initHubOrient.shape[0] != 9:
            print("Expecting a 9 element array for initHubOrient DCM [r11,r12,r13,r21,r22,r23,r31,r32,r33]")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn terminated prematurely.")
        if np.squeeze(self.initNacellePos.ndim) > 1 or self.initNacellePos.shape[0] != 3:
            print("Expecting a 3 element array for initNacellePos [x,y,z]")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn terminated prematurely.")
        if np.squeeze(self.initNacelleOrient.ndim) > 1 or self.initNacelleOrient.shape[0] != 9:
            print("Expecting a 9 element array for initNacelleOrient DCM [r11,r12,r13,r21,r22,r23,r31,r32,r33]")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn terminated prematurely.")
            

    def check_init_mesh(self):
        #print("shape of initMeshPos       ",   self.initMeshPos.shape)
        #print("               size 0      ",   self.initMeshPos.shape[0])
        #print("               size 1      ",   self.initMeshPos.shape[1])
        #print("shape of initMeshOrient    ",   self.initMeshOrient.shape)
        #print("               size 0      ",   self.initMeshOrient.shape[0])
        #print("               size 1      ",   self.initMeshOrient.shape[1])
        #print("               float       ",   type(self.initMeshOrient[0,0]))
        # initMeshPos
        #   Verify that the shape of initMeshPos is correct
        if self.initMeshPos.shape[0] != self.initMeshOrient.shape[0]:
            print("Different number of meshs in inital position and orientation arrays")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn terminated prematurely.")
        if self.initMeshPos.shape[1] != 3:
            print("Expecting a Nx3 array of initial mesh positions (initMeshPos) with second index for [x,y,z]")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn terminated prematurely.")
        if self.initMeshPos.shape[0] != self.numMeshPts:
            print("Expecting a Nx3 array of initial mesh positions (initMeshPos) with first index for mesh number.")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn terminated prematurely.")
        if self.initMeshOrient.shape[1] != 9:
            print("Expecting a Nx9 array of initial mesh orientations as DCMs (initMeshOrient) with second index for [r11,r12,r13,r21,r22,r23,r31,r32,r33]")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn terminated prematurely.")
        if self.initMeshOrient.shape[0] != self.numMeshPts:
            print("Expecting a Nx3 array of initial mesh orientations (initMeshOrient) with first index for mesh number.")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn terminated prematurely.")


    def check_input_motions_hubNac(self,nodePos,nodeOrient,nodeVel,nodeAcc,_name):
        #   Verify that the shape of positions array is correct
        if nodePos.size != 3:
            print("Expecting a Nx3 array of "+_name+" positions with second index for [x,y,z]")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn/AeroDyn terminated prematurely.")

        #   Verify that the shape of orientations array is correct
        if nodeOrient.size != 9:
            print("Expecting a Nx9 array of "+_name+" orientations with second index for [r11,r12,r13,r21,r22,r23,r31,r32,r33]")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn/AeroDyn terminated prematurely.")

        #   Verify that the shape of velocities array is correct
        if nodeVel.size != 6:
            print("Expecting a Nx6 array of "+_name+" velocities with second index for [x,y,z,Rx,Ry,Rz]")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn/AeroDyn terminated prematurely.")

        #   Verify that the shape of accelerations array is correct
        if nodeAcc.size != 6:
            print("Expecting a Nx6 array of "+_name+" accelerations with second index for [x,y,z,Rx,Ry,Rz]")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn/AeroDyn terminated prematurely.")

    def check_input_motions_root(self,rootPos,rootOrient,rootVel,rootAcc):
        # make sure number of roots didn't change for some reason
        if self._initNumBlades != self.numBlades:
            print(f"At time {time}, the number of root points changed from initial value of {self._initNumBlades}.  This is not permitted during the simulation.")
            self.aerodyn_inflow_end()
            raise Exception("\nError in calling AeroDyn/AeroDyn library.")

        #   Verify that the shape of positions array is correct
        if rootPos.shape[1] != 3:
            print("Expecting a Nx3 array of root positions (rootOrient) with second index for [x,y,z]")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn/AeroDyn terminated prematurely.")
        if rootPos.shape[0] != self.numBlades:
            print("Expecting a Nx3 array of root positions (rootOrient) with first index for root number.")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn/AeroDyn terminated prematurely.")

        #   Verify that the shape of orientations array is correct
        if rootOrient.shape[1] != 9:
            print("Expecting a Nx9 array of root orientations (rootPos) with second index for [r11,r12,r13,r21,r22,r23,r31,r32,r33]")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn/AeroDyn terminated prematurely.")
        if rootOrient.shape[0] != self.numBlades:
            print("Expecting a Nx9 array of root orientations (rootPos) with first index for root number.")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn/AeroDyn terminated prematurely.")

        #   Verify that the shape of velocities array is correct
        if rootVel.shape[1] != 6:
            print("Expecting a Nx6 array of root velocities (rootVel) with second index for [x,y,z,Rx,Ry,Rz]")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn/AeroDyn terminated prematurely.")
        if rootVel.shape[0] != self.numBlades:
            print("Expecting a Nx6 array of root velocities (rootVel) with first index for root number.")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn/AeroDyn terminated prematurely.")

        #   Verify that the shape of accelerations array is correct
        if rootAcc.shape[1] != 6:
            print("Expecting a Nx6 array of root accelerations (rootAcc) with second index for [x,y,z,Rx,Ry,Rz]")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn/AeroDyn terminated prematurely.")
        if rootAcc.shape[0] != self.numBlades:
            print("Expecting a Nx6 array of root accelerations (rootAcc) with first index for root number.")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn/AeroDyn terminated prematurely.")


    def check_input_motions_mesh(self,meshPos,meshOrient,meshVel,meshAcc):
        # make sure number of meshs didn't change for some reason
        if self._initNumMeshPts != self.numMeshPts:
            print(f"At time {time}, the number of mesh points changed from initial value of {self._initNumMeshPts}.  This is not permitted during the simulation.")
            self.aerodyn_inflow_end()
            raise Exception("\nError in calling AeroDyn/AeroDyn library.")

        #   Verify that the shape of positions array is correct
        if meshPos.shape[1] != 3:
            print("Expecting a Nx3 array of mesh positions (meshOrient) with second index for [x,y,z]")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn/AeroDyn terminated prematurely.")
        if meshPos.shape[0] != self.numMeshPts:
            print("Expecting a Nx3 array of mesh positions (meshOrient) with first index for mesh number.")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn/AeroDyn terminated prematurely.")

        #   Verify that the shape of orientations array is correct
        if meshOrient.shape[1] != 9:
            print("Expecting a Nx9 array of mesh orientations (meshPos) with second index for [r11,r12,r13,r21,r22,r23,r31,r32,r33]")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn/AeroDyn terminated prematurely.")
        if meshOrient.shape[0] != self.numMeshPts:
            print("Expecting a Nx9 array of mesh orientations (meshPos) with first index for mesh number.")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn/AeroDyn terminated prematurely.")

        #   Verify that the shape of velocities array is correct
        if meshVel.shape[1] != 6:
            print("Expecting a Nx6 array of mesh velocities (meshVel) with second index for [x,y,z,Rx,Ry,Rz]")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn/AeroDyn terminated prematurely.")
        if meshVel.shape[0] != self.numMeshPts:
            print("Expecting a Nx6 array of mesh velocities (meshVel) with first index for mesh number.")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn/AeroDyn terminated prematurely.")

        #   Verify that the shape of accelerations array is correct
        if meshAcc.shape[1] != 6:
            print("Expecting a Nx6 array of mesh accelerations (meshAcc) with second index for [x,y,z,Rx,Ry,Rz]")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn/AeroDyn terminated prematurely.")
        if meshAcc.shape[0] != self.numMeshPts:
            print("Expecting a Nx6 array of mesh accelerations (meshAcc) with first index for mesh number.")
            self.aerodyn_inflow_end()
            raise Exception("\nAeroDyn/AeroDyn terminated prematurely.")



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
#   forces and moments at each input mesh point.  If all is functioning
#   correctly, this will be identical to the corresponding values in the
#   AeroDyn/InflowWind output channels.

#FIXME: this is incorrect
class DriverDbg():
    """
    This is only for debugging purposes only.  The input motions and resulting
    forces can be written to file with this class to verify the data I/O to the
    Fortran library.
    When coupled to another code, the force/moment array would be passed back
    to the calling code for use in the structural solver.
    """
    def __init__(self,filename,numMeshPts):
        self.DbgFile=open(filename,'wt')        # open output file and write header info
        self.numMeshPts=numMeshPts
        # write file header
        t_string=datetime.datetime.now()
        dt_string=datetime.date.today()
        self.DbgFile.write(f"## This file was generated by aerodyn_inflow_c_lib on {dt_string.strftime('%b-%d-%Y')} at {t_string.strftime('%H:%M:%S')}\n")
        self.DbgFile.write(f"## This file contains the resulting forces/moments at each of {self.numMeshPts} mesh points passed into the aerodyn_inflow_c_lib\n")
        self.DbgFile.write("#\n")
        self.DbgFile.write("#\n")
        self.DbgFile.write("#\n")
        self.DbgFile.write("#\n")
        f_string = "{:^25s}"
        self.DbgFile.write("       Time    ")
        for i in range(1,self.numMeshPts+1):
            f_num = "N{0:04d}_".format(i)
            self.DbgFile.write(f_string.format(f_num+"x"  ))
            self.DbgFile.write(f_string.format(f_num+"y"  ))
            self.DbgFile.write(f_string.format(f_num+"z"  ))
            #self.DbgFile.write(f_string.format(f_num+"Rx" ))
            #self.DbgFile.write(f_string.format(f_num+"Ry" ))
            #self.DbgFile.write(f_string.format(f_num+"Rz" ))
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
        self.DbgFile.write("       (s)     ")
        for i in range(1,self.numMeshPts+1):
            self.DbgFile.write(f_string.format("(m)"      ))
            self.DbgFile.write(f_string.format("(m)"      ))
            self.DbgFile.write(f_string.format("(m)"      ))
            #self.DbgFile.write(f_string.format("(rad)"    ))
            #self.DbgFile.write(f_string.format("(rad)"    ))
            #self.DbgFile.write(f_string.format("(rad)"    ))
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

    def write(self,t,meshPos,meshVel,meshAcc,meshFrc):
        t_string  = "{:10.4f}"
        f_string3 = "{:25.7e}"*3
        f_string6 = "{:25.7e}"*6
        self.DbgFile.write(t_string.format(t))
        for i in range(0,self.numMeshPts):
            self.DbgFile.write(f_string3.format(*meshPos[i,:]))
            self.DbgFile.write(f_string6.format(*meshVel[i,:]))
            self.DbgFile.write(f_string6.format(*meshAcc[i,:]))
            self.DbgFile.write(f_string6.format(*meshFrc[i,:]))
        self.DbgFile.write("\n")

    def end(self):
        if self.opened:
            self.DbgFile.close()
            self.opened = False


#===============================================================================
#   Helper class for writing channels to file.
#   for the regression testing to mirror the output from the AD15 and InfowWind
#   from an OpenFAST simulation.  This may also have value for debugging
#   interfacing to the AeroDyn_Inflow_C_Binding library

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
        self.OutFile.write(f"## This file was generated by AeroDyn_Inflow_Driver on {dt_string.strftime('%b-%d-%Y')} at {t_string.strftime('%H:%M:%S')}\n")
        self.OutFile.write(f"## This file contains output channels requested from the OutList section of the AD15 and IfW input files")
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
