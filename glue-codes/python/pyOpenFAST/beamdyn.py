#**********************************************************************************************************************************
# LICENSING
# Copyright (C) 2026 National Renewable Energy Laboratory
#
# This file is part of BeamDyn.
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
# This is the Python-C interface library for BeamDyn
# Usage: THIS LIBRARY IS NOT TO BE CHANGED OR EDITED BY THE USER
#
# Conventions of the BeamDyn C interface (see BeamDyn_C_Binding.f90 / BeamDyn_C_Binding.h):
#   - Direction cosine matrices (DCMs) map a vector from the global frame to the local frame.  They cross the
#     interface as flat arrays stored row by row; here they are passed and returned as (3,3) or (N,3,3) arrays.
#   - Root kinematics (position, orientation, velocity, acceleration), gravity, and node positions are double precision.
#   - Loads, node velocities and accelerations, reaction loads, and output channel values are single precision.
#   - bd_calcOutput must be called at the initial time before the first bd_updateStates.
from ctypes import (
    POINTER,
    create_string_buffer,
    byref,
    c_int,
    c_double,
    c_float,
    c_char,
    c_char_p,
)
import datetime

import numpy as np

from .interface_abc import OpenFASTInterfaceType


class BeamDynLib(OpenFASTInterfaceType):

    def __init__(self, library_path):
        super().__init__(library_path)

        self._initialize_routines()
        self.ended = False   # For error handling at end

        # Buffers for the output channel names and units (ChanLen characters per channel)
        self._channel_names = create_string_buffer(self.CHANNEL_NAME_LENGTH * self.MAX_CHANNELS + 1)
        self._channel_units = create_string_buffer(self.CHANNEL_NAME_LENGTH * self.MAX_CHANNELS + 1)

        # Sizes returned by BD_C_Init
        self._numOutputNodes    = c_int(0)
        self._numPointLoadNodes = c_int(0)
        self._numDistrLoadNodes = c_int(0)
        self._numChannels       = c_int(0)

        # Time step (set by bd_init)
        self.dt = 0.0

    # Initialize routines ------------------------------------------------------------------------------------------------------------
    def _initialize_routines(self):
        self.BD_C_SetVTK.argtypes = [
            POINTER(c_int),                       # IN: VTK output switch
            POINTER(c_double),                    # IN: time between animation frames
            POINTER(c_int),                       # OUT: ErrStat_C
            POINTER(c_char)                       # OUT: ErrMsg_C
        ]
        self.BD_C_SetVTK.restype = None

        self.BD_C_Init.argtypes = [
            POINTER(c_int),                       # IN: input file passed
            POINTER(c_char_p),                    # IN: input file string
            POINTER(c_int),                       # IN: input file string length
            c_char_p,                             # IN: output root name
            POINTER(c_double),                    # IN: root position (3)
            POINTER(c_double),                    # IN: root orientation DCM (9)
            POINTER(c_double),                    # IN: root velocity (6)
            POINTER(c_double),                    # IN: gravity vector (3)
            POINTER(c_double),                    # IN: dt
            POINTER(c_int),                       # IN: interpolation order
            POINTER(c_int),                       # IN: dynamic solve flag
            POINTER(c_int),                       # IN: write output file flag
            POINTER(c_int),                       # OUT: number of output nodes
            POINTER(c_int),                       # OUT: number of point load nodes
            POINTER(c_int),                       # OUT: number of distributed load nodes
            POINTER(c_int),                       # OUT: number of channels
            POINTER(c_char),                      # OUT: output channel names
            POINTER(c_char),                      # OUT: output channel units
            POINTER(c_int),                       # OUT: ErrStat_C
            POINTER(c_char)                       # OUT: ErrMsg_C
        ]
        self.BD_C_Init.restype = None

        self.BD_C_GetRefPositions.argtypes = [
            POINTER(c_double),                     # OUT: output node reference positions (3 x NumOutputNodes)
            POINTER(c_double),                    # OUT: output node reference orientations (9 x NumOutputNodes)
            POINTER(c_double),                     # OUT: point load node reference positions (3 x NumPointLoadNodes)
            POINTER(c_double),                     # OUT: distributed load node reference positions (3 x NumDistrLoadNodes)
            POINTER(c_int),                       # OUT: ErrStat_C
            POINTER(c_char)                       # OUT: ErrMsg_C
        ]
        self.BD_C_GetRefPositions.restype = None

        self.BD_C_SetRootMotion.argtypes = [
            POINTER(c_double),                    # IN: root displacement (3)
            POINTER(c_double),                    # IN: root orientation DCM (9)
            POINTER(c_double),                    # IN: root velocity (6)
            POINTER(c_double),                    # IN: root acceleration (6)
            POINTER(c_int),                       # OUT: ErrStat_C
            POINTER(c_char)                       # OUT: ErrMsg_C
        ]
        self.BD_C_SetRootMotion.restype = None

        self.BD_C_SetPointLoads.argtypes = [
            POINTER(c_float),                     # IN: point loads (6 x NumPointLoadNodes)
            POINTER(c_int),                       # OUT: ErrStat_C
            POINTER(c_char)                       # OUT: ErrMsg_C
        ]
        self.BD_C_SetPointLoads.restype = None

        self.BD_C_SetDistrLoads.argtypes = [
            POINTER(c_float),                     # IN: distributed loads (6 x NumDistrLoadNodes)
            POINTER(c_int),                       # OUT: ErrStat_C
            POINTER(c_char)                       # OUT: ErrMsg_C
        ]
        self.BD_C_SetDistrLoads.restype = None

        self.BD_C_UpdateStates.argtypes = [
            POINTER(c_double),                    # IN: time @ n
            POINTER(c_double),                    # IN: time @ n+1
            POINTER(c_int),                       # OUT: ErrStat_C
            POINTER(c_char)                       # OUT: ErrMsg_C
        ]
        self.BD_C_UpdateStates.restype = None

        self.BD_C_CalcOutput.argtypes = [
            POINTER(c_double),                    # IN: time
            POINTER(c_double),                     # OUT: node positions (3 x NumOutputNodes)
            POINTER(c_double),                    # OUT: node orientations (9 x NumOutputNodes)
            POINTER(c_float),                     # OUT: node velocities (6 x NumOutputNodes)
            POINTER(c_float),                     # OUT: node accelerations (6 x NumOutputNodes)
            POINTER(c_float),                     # OUT: root reaction force and moment (6)
            POINTER(c_float),                     # OUT: output channel values
            POINTER(c_int),                       # OUT: ErrStat_C
            POINTER(c_char)                       # OUT: ErrMsg_C
        ]
        self.BD_C_CalcOutput.restype = None

        self.BD_C_PackStates.argtypes = [
            c_char_p,                             # IN: checkpoint root name
            POINTER(c_int),                       # OUT: ErrStat_C
            POINTER(c_char)                       # OUT: ErrMsg_C
        ]
        self.BD_C_PackStates.restype = None

        self.BD_C_UnpackStates.argtypes = [
            c_char_p,                             # IN: checkpoint root name
            POINTER(c_int),                       # OUT: ErrStat_C
            POINTER(c_char)                       # OUT: ErrMsg_C
        ]
        self.BD_C_UnpackStates.restype = None

        self.BD_C_End.argtypes = [
            POINTER(c_int),                       # OUT: ErrStat_C
            POINTER(c_char)                       # OUT: ErrMsg_C
        ]
        self.BD_C_End.restype = None

    # bd_setVTK ------------------------------------------------------------------------------------------------------------
    def bd_setVTK(self, write_vtk, vtk_dt=0.0):
        """Request VTK output of the meshes (call before bd_init).

        write_vtk -- 0: none; 1: reference meshes at initialization; 2: reference meshes and animation
        vtk_dt    -- time between animation frames (s); rounded to a multiple of the time step, at least one step
        """
        self.BD_C_SetVTK(
            byref(c_int(write_vtk)),               # IN: VTK output switch
            byref(c_double(vtk_dt)),               # IN: time between animation frames
            byref(self.error_status_c),            # OUT: ErrStat_C
            self.error_message_c                   # OUT: ErrMsg_C
        )

        self.check_error()

    # bd_init ------------------------------------------------------------------------------------------------------------
    def bd_init(self, input_file_passed, input_string_array, out_root_name,
                root_position, root_orientation, root_velocity, gravity, dt,
                interp_order=1, dynamic_solve=True, write_outputs=False):
        """Initialize the beam.

        input_file_passed  -- 1: input_string_array holds the lines of the primary input file;
                              0: input_string_array is the path to the primary input file
        out_root_name      -- root name for the summary, echo, and (optional) output files
        root_position      -- (3) initial root position, global frame (m)
        root_orientation   -- (3,3) initial root orientation DCM (global to root)
        root_velocity      -- (6) initial root translational and rotational velocity, global frame (m/s, rad/s)
        gravity            -- (3) gravity vector, global frame (m/s^2)
        dt                 -- time step (s), > 0
        """
        self.dt = float(dt)
        self.ended = False

        # Convert the string into a c_char byte array (a single string is a path, not a list of lines)
        if isinstance(input_string_array, str):
            input_string_array = [input_string_array]
        input_string = '\x00'.join(input_string_array)
        input_string = input_string.encode('utf-8')
        input_string_length = len(input_string)

        root_pos_c  = self._to_c(root_position, c_double, 3)
        root_ori_c  = self._to_c(root_orientation, c_double, 9)
        root_vel_c  = self._to_c(root_velocity, c_double, 6)
        gravity_c   = self._to_c(gravity, c_double, 3)

        # The root name is read up to the first null character (at most IntfStrLen characters)
        root_name_c = create_string_buffer(str(out_root_name).encode('utf-8'), self.IntfStrLen)

        self.BD_C_Init(
            byref(c_int(input_file_passed)),       # IN: input file passed
            c_char_p(input_string),                # IN: input file string
            byref(c_int(input_string_length)),     # IN: input file string length
            root_name_c,                           # IN: output root name
            root_pos_c,                            # IN: root position
            root_ori_c,                            # IN: root orientation DCM
            root_vel_c,                            # IN: root velocity
            gravity_c,                             # IN: gravity vector
            byref(c_double(self.dt)),              # IN: time step (dt)
            byref(c_int(interp_order)),            # IN: interpolation order
            byref(c_int(1 if dynamic_solve else 0)),   # IN: dynamic solve flag
            byref(c_int(1 if write_outputs else 0)),   # IN: write output file flag
            byref(self._numOutputNodes),           # OUT: number of output nodes
            byref(self._numPointLoadNodes),        # OUT: number of point load nodes
            byref(self._numDistrLoadNodes),        # OUT: number of distributed load nodes
            byref(self._numChannels),              # OUT: number of channels
            self._channel_names,                   # OUT: output channel names
            self._channel_units,                   # OUT: output channel units
            byref(self.error_status_c),            # OUT: ErrStat_C
            self.error_message_c                   # OUT: ErrMsg_C
        )

        self.check_error()

    # bd_getRefPositions ------------------------------------------------------------------------------------------------------------
    def bd_getRefPositions(self):
        """Return the reference positions (N,3) and orientations (N,3,3) of the output nodes and the reference
        positions of the point load nodes (Np,3) and distributed load nodes (Nd,3)."""
        n_out, n_pl, n_dl = self.num_output_nodes, self.num_point_load_nodes, self.num_distr_load_nodes
        out_pos_c = (c_double * (3 * n_out))()
        out_ori_c = (c_double * (9 * n_out))()
        pl_pos_c  = (c_double * (3 * n_pl))()
        dl_pos_c  = (c_double * (3 * n_dl))()

        self.BD_C_GetRefPositions(
            out_pos_c,                             # OUT: output node reference positions
            out_ori_c,                             # OUT: output node reference orientations
            pl_pos_c,                              # OUT: point load node reference positions
            dl_pos_c,                              # OUT: distributed load node reference positions
            byref(self.error_status_c),            # OUT: ErrStat_C
            self.error_message_c                   # OUT: ErrMsg_C
        )

        self.check_error()

        return (
            np.array(out_pos_c, dtype=float).reshape(n_out, 3),
            np.array(out_ori_c, dtype=float).reshape(n_out, 3, 3),
            np.array(pl_pos_c,  dtype=float).reshape(n_pl, 3),
            np.array(dl_pos_c,  dtype=float).reshape(n_dl, 3),
        )

    # bd_setRootMotion ------------------------------------------------------------------------------------------------------------
    def bd_setRootMotion(self, displacement, orientation, velocity, acceleration):
        """Set the root motion used by the next bd_calcOutput (motion at the current time) or bd_updateStates
        (motion at the next time) call."""
        self.BD_C_SetRootMotion(
            self._to_c(displacement, c_double, 3),  # IN: root displacement
            self._to_c(orientation,  c_double, 9),  # IN: root orientation DCM
            self._to_c(velocity,     c_double, 6),  # IN: root velocity
            self._to_c(acceleration, c_double, 6),  # IN: root acceleration
            byref(self.error_status_c),             # OUT: ErrStat_C
            self.error_message_c                    # OUT: ErrMsg_C
        )

        self.check_error()

    # bd_setPointLoads ------------------------------------------------------------------------------------------------------------
    def bd_setPointLoads(self, point_loads):
        """Set the point loads: (NumPointLoadNodes,6) force and moment at each point load node, global frame."""
        self.BD_C_SetPointLoads(
            self._to_c(point_loads, c_float, 6 * self.num_point_load_nodes),   # IN: point loads
            byref(self.error_status_c),            # OUT: ErrStat_C
            self.error_message_c                   # OUT: ErrMsg_C
        )

        self.check_error()

    # bd_setDistrLoads ------------------------------------------------------------------------------------------------------------
    def bd_setDistrLoads(self, distr_loads):
        """Set the distributed loads: (NumDistrLoadNodes,6) force and moment per unit length at each distributed
        load node, global frame."""
        self.BD_C_SetDistrLoads(
            self._to_c(distr_loads, c_float, 6 * self.num_distr_load_nodes),   # IN: distributed loads
            byref(self.error_status_c),            # OUT: ErrStat_C
            self.error_message_c                   # OUT: ErrMsg_C
        )

        self.check_error()

    # bd_updateStates ------------------------------------------------------------------------------------------------------------
    def bd_updateStates(self, t1, t2):
        """Advance the states from t1 to t2 using the inputs set for t2."""
        self.BD_C_UpdateStates(
            byref(c_double(t1)),                   # IN: current time (t)
            byref(c_double(t2)),                   # IN: next time step (t+dt)
            byref(self.error_status_c),            # OUT: ErrStat_C
            self.error_message_c                   # OUT: ErrMsg_C
        )

        self.check_error()

    # bd_calcOutput ------------------------------------------------------------------------------------------------------------
    def bd_calcOutput(self, t, output_channel_values):
        """Calculate the outputs at time t.  The output channel values are written into output_channel_values.

        Returns the node positions (N,3), orientations (N,3,3), velocities (N,6), accelerations (N,6), and the
        root reaction force and moment (6)."""
        n_out = self.num_output_nodes
        node_pos_c = (c_double * (3 * n_out))()
        node_ori_c = (c_double * (9 * n_out))()
        node_vel_c = (c_float  * (6 * n_out))()
        node_acc_c = (c_float  * (6 * n_out))()
        reaction_c = (c_float  * 6)()
        outputs_c  = (c_float  * max(self._numChannels.value, 1))()

        self.BD_C_CalcOutput(
            byref(c_double(t)),                    # IN: time
            node_pos_c,                            # OUT: node positions
            node_ori_c,                            # OUT: node orientations
            node_vel_c,                            # OUT: node velocities
            node_acc_c,                            # OUT: node accelerations
            reaction_c,                            # OUT: root reaction force and moment
            outputs_c,                             # OUT: output channel values
            byref(self.error_status_c),            # OUT: ErrStat_C
            self.error_message_c                   # OUT: ErrMsg_C
        )

        self.check_error()

        for i in range(self._numChannels.value):
            output_channel_values[i] = outputs_c[i]

        return (
            np.array(node_pos_c, dtype=float).reshape(n_out, 3),
            np.array(node_ori_c, dtype=float).reshape(n_out, 3, 3),
            np.array(node_vel_c, dtype=float).reshape(n_out, 6),
            np.array(node_acc_c, dtype=float).reshape(n_out, 6),
            np.array(reaction_c, dtype=float),
        )

    # bd_packStates / bd_unpackStates ------------------------------------------------------------------------------------------------
    def bd_packStates(self, checkpoint_root):
        """Write the checkpoint file <checkpoint_root>.chkp with the complete state of the beam."""
        self.BD_C_PackStates(
            create_string_buffer(str(checkpoint_root).encode('utf-8'), self.IntfStrLen),   # IN: checkpoint root name
            byref(self.error_status_c),            # OUT: ErrStat_C
            self.error_message_c                   # OUT: ErrMsg_C
        )

        self.check_error()

    def bd_unpackStates(self, checkpoint_root):
        """Restore the state of the beam from <checkpoint_root>.chkp (call after bd_init with the same input file)."""
        self.BD_C_UnpackStates(
            create_string_buffer(str(checkpoint_root).encode('utf-8'), self.IntfStrLen),   # IN: checkpoint root name
            byref(self.error_status_c),            # OUT: ErrStat_C
            self.error_message_c                   # OUT: ErrMsg_C
        )

        self.check_error()

    # bd_end ------------------------------------------------------------------------------------------------------------
    def bd_end(self):

        if not self.ended:
            self.ended = True
            self.BD_C_End(
                byref(self.error_status_c),            # OUT: ErrStat_C
                self.error_message_c                   # OUT: ErrMsg_C
            )
            self.check_error()

    #===============================================================================
    # OTHER FUNCTIONS --------------------------------------------------------------------------------------------------

    @staticmethod
    def _to_c(values, c_type, length):
        """Copy values (any shape, flattened row by row) into a ctypes array of the given type and length."""
        flat = np.asarray(values, dtype=float).reshape(-1)
        if flat.size != length:
            raise ValueError(f"Expected {length} values, got {flat.size}")
        return (c_type * length)(*flat)

    # Error Handling Functions
    def check_error(self):
        if self.error_status_c.value == 0:
            return
        elif self.error_status_c.value < self.abort_error_level.value:
            print(f"BeamDyn error status: {self.error_levels[self.error_status_c.value]}: {self.error_message_c.value.decode('ascii', errors='replace')}")
        else:
            print(f"BeamDyn error status: {self.error_levels[self.error_status_c.value]}: {self.error_message_c.value.decode('ascii', errors='replace')}")
            self.bd_end()
            raise Exception("\nBeamDyn terminated prematurely.")

    # Mesh sizes
    @property
    def num_output_nodes(self):
        return self._numOutputNodes.value

    @property
    def num_point_load_nodes(self):
        return self._numPointLoadNodes.value

    @property
    def num_distr_load_nodes(self):
        return self._numDistrLoadNodes.value

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
#   For the regression testing to mirror the output from the BeamDyn Fortran
#   driver.  This may also have value for debugging the interfacing to BD.

class WriteOutChans():
    """
    This is only for testing purposes. Since we are not returning the
    output channels to anything, we will write them to file.  When coupled to
    another code, this data would be passed back for inclusion the any output
    file there.
    """
    def __init__(self,filename,chan_names,chan_units):
        chan_names = ['Time'] + list(chan_names)    # add time index header
        chan_units = ['(s)'] + list(chan_units)     # add time index unit
        self.OutFile=open(filename,'wt')            # open output file and write header info
        # write file header
        t_string=datetime.datetime.now()
        dt_string=datetime.date.today()
        self.OutFile.write(f"## This file was generated by BeamDyn c-bindings library on {dt_string.strftime('%b-%d-%Y')} at {t_string.strftime('%H:%M:%S')}\n")
        self.OutFile.write(f"## This file contains output channels requested from the OutList section of the input file\n")
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
        f_string = "{:15.6f}"+"{:25.8E}"*(l-1)
        for i in range(0,chan_data.shape[0]):
            self.OutFile.write(f_string.format(*chan_data[i,:]) + '\n')

    def end(self):
        if self.opened:
            self.OutFile.close()
            self.opened = False
