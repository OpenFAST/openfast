#-------------------------------------------------------------------------------
# LICENSING
#-------------------------------------------------------------------------------
# Copyright (C) 2021-present by National Renewable Energy Lab (NREL)
#
# This file is part of AeroDyn
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

#-------------------------------------------------------------------------------
# Overview
#-------------------------------------------------------------------------------
# This is the Python-C interface library for AeroDyn with InflowWind. This may
# be used directly with Python based codes to call and run AeroDyn and
# InflowWind together.

#--------------------------------------
# Key Features
#--------------------------------------
# - AeroDyn: Blade element momentum (BEM) theory for aerodynamic force
#   calculations (dynamic stall, unsteady aerodynamics, tower shadow, wind shear,
#   tip/hub losses)
# - InflowWind: Simulates complex inflow conditions (turbulence, wind shear,
#   gusts)

#--------------------------------------
# Usage
#--------------------------------------
# 1. Instantiate AeroDynInflowLib with shared library path
# 2. Initialize: adi_preinit() -> adi_setuprotor() -> adi_init()
# 3. Simulate: adi_setrotormotion() -> adi_updateStates() ->
#    adi_calcOutput() -> adi_getrotorloads()
# 4. adi_end() and handle errors via check_error()

# An example of using this library from Python is given in the accompanying
# Python driver program(s) in reg_tests/r-test/modules/aerodyn directory.
# Additional notes and information on the interfacing is included there.

#-------------------------------------------------------------------------------
# Imports
#-------------------------------------------------------------------------------
import os
from ctypes import (CDLL, POINTER, byref, c_char, c_char_p, c_double, c_float,
                    c_int, create_string_buffer)
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import numpy.typing as npt

#-------------------------------------------------------------------------------
# Helper functions and classes
#-------------------------------------------------------------------------------
def flatten_array(
    initial_count: int,
    current_count: int,
    array: npt.NDArray,
    array_name: str,
    elements_per_item: int,
    c_type: Any = c_float
) -> Any:
    """Flattens arrays for passing to C.

    This is a helper function to flatten arrays for passing to C. It is used to flatten
    positions, orientations, and velocities/accelerations.

    Args:
        initial_count: Initial number of items from setup
        current_count: Current number of items being processed
        array: Numpy array to flatten
        array_name: Descriptive name for error messages
        elements_per_item: Number of elements per item (e.g. 3 for positions, 9 for orientations)
        c_type: C type to convert to (default: c_float)

    Returns:
        Array: Flattened C array of size (elements_per_item * current_count)

    Raises:
        RuntimeError: If current_count differs from initial_count
    """
    if initial_count != current_count:
        error_msg = (
            f"The number of {array_name} points changed from initial value of"
            f" {initial_count}. This is not permitted during the simulation."
        )
        raise RuntimeError(f"Error in calling AeroDyn/InflowWind library: {error_msg}")

    c_array = (c_type * (elements_per_item * current_count))(*array.flatten())
    return c_array

def to_c_array(array: npt.NDArray, c_type: Any = c_float) -> Any:
    """Converts numpy array to C array of specified type.

    Args:
        array: Input numpy array
        c_type: C type to convert to (default: c_float)

    Returns:
        C-compatible array of the specified type
    """
    try:
        if isinstance(array, np.ndarray):
            flat_array = array.flatten()
            return (c_type * len(flat_array))(*flat_array)
        # If list/tuple, convert directly to C array
        return (c_type * len(array))(*array)
    except Exception as e:
        raise TypeError(f"Failed to convert to C array: {e}")

def to_c_string(input_array: List[str]) -> Tuple[bytes, int]:
    """Converts input string array into a null-separated byte string for use in C.

    Args:
        input_array: List of strings to join with null characters

    Returns:
        Tuple containing:
            - The encoded byte string
            - Length of the encoded string
    """
    encoded_string = '\x00'.join(input_array).encode('utf-8')
    return encoded_string, len(encoded_string)

@dataclass
class MotionData:
    """POD-style container for motion-related data i.e. state of a node."""
    position: npt.NDArray[np.float32]
    orientation: npt.NDArray[np.float64]
    velocity: npt.NDArray[np.float32]
    acceleration: npt.NDArray[np.float32]

#-------------------------------------------------------------------------------
# C-interface library class for AeroDyn x InflowWind
#-------------------------------------------------------------------------------
class AeroDynInflowLib(CDLL):
    """A Python interface to the AeroDyn/InflowWind library.

    This class provides a modern Python interface for calling and running AeroDyn
    and InflowWind together. It handles initialization, runtime operations, and cleanup
    of the underlying Fortran library.
    """

    #--------------------------------------
    # Error levels (from IfW)
    #--------------------------------------
    error_levels: Dict[int, str] = {
        0: "None",
        1: "Info",
        2: "Warning",
        3: "Severe Error",
        4: "Fatal Error"
    }

    #--------------------------------------
    # Constants
    #--------------------------------------
    # NOTE: The length of the error message in Fortran is determined by the
    #       ErrMsgLen variable in the NWTC_Base.f90 file. If ErrMsgLen is modified,
    #       the corresponding size here must also be updated to match.
    ERROR_MESSAGE_LENGTH: int = 1025
    DEFAULT_STRING_LENGTH: int = 1025
    CHANNEL_NAME_LENGTH: int = 20
    MAX_CHANNELS: int = 8000

    def __init__(self, library_path: Union[str, Path]) -> None:
        """Initializes the AeroDyn/InflowWind interface.

        Args:
            library_path: Path to the compiled library file (.dll, .so, or .dylib)

        Raises:
            FileNotFoundError: If the library file cannot be found
            OSError: If the library cannot be loaded
        """
        library_path = Path(library_path)
        if not library_path.exists():
            raise FileNotFoundError(f"Library not found at: {library_path}")

        super().__init__(str(library_path))
        self.library_path = library_path

        # Initialize the library interface
        self._initialize_routines()
        self.ended = False

        # Input file handling configuration
        self.aerodyn_inputs_passed_as_string: bool = True  # Pass input file as string
        self.inflow_inputs_passed_as_string: bool = True   # Pass input file as string

        # Error handling setup
        self.abort_error_level = 4
        self.error_status_c = c_int(0)
        self.error_message_c = create_string_buffer(self.ERROR_MESSAGE_LENGTH)

        # Channel information buffers
        self._channel_names_c = create_string_buffer(
            self.CHANNEL_NAME_LENGTH * self.MAX_CHANNELS
        )
        self._channel_units_c = create_string_buffer(
            self.CHANNEL_NAME_LENGTH * self.MAX_CHANNELS
        )

        #--------------------------------------
        # Flags
        #--------------------------------------
        # For flags: 0->false, 1->true
        self.store_hub_height_velocity = 1
        self.transpose_dcm = 1
        self.point_load_output = 1

        # 0->None, 1->Info, 2->Warning, 3->Severe Error, 4->Fatal Error
        self.debug_level = 0

        #--------------------------------------
        # VTK settings
        #--------------------------------------
        self.write_vtk = 0         # Default -> no vtk output
        self.vtk_type = 1          # Default -> surface meshes
        self.vtk_dt = 0.           # Default -> all
        self.vtk_nacelle_dimension = np.array(
            # Default nacelle dimension [x0, y0, z0, Lx, Ly, Lz] (m)
            [-2.5, -2.5, 0., 10., 5., 5.], dtype="float32"
        )
        self.vtk_hub_radius = 1.5  # Default hub radius for VTK surface rendering

        #--------------------------------------
        # Output file settings
        #--------------------------------------
        self.write_outputs = 0               # File format for writing outputs
        self.output_timestep = 0.            # Timestep for outputs to file
        self.output_root_name = "Output_ADIlib_default"
        self.output_vtk_dir = ""             # Set to specify a directory relative to input files

        #--------------------------------------
        # Interpolation settings
        #--------------------------------------
        # Options: 1 -> linear, 2 -> quadratic
        self.interpolation_order = 1         # Default -> linear interpolation

        #--------------------------------------
        # Time settings
        #--------------------------------------
        self.dt = 0.1               # Typical default for HD
        self.t_max = 600.           # Typical default for HD waves FFT
        self.total_time = 0.        # May be longer than t_max
        self.num_time_steps = 0

        #--------------------------------------
        # Channels and turbines
        #--------------------------------------
        self.num_channels = 0         # Number of channels returned
        self.num_turbines = 1

        #--------------------------------------
        # Initial position setup
        #--------------------------------------
        self.init_hub_pos = np.zeros(shape=(3), dtype=c_float)
        self.init_hub_orient = np.zeros(shape=(9), dtype=c_double)
        self.init_nacelle_pos = np.zeros(shape=(3), dtype=c_float)
        self.init_nacelle_orient = np.zeros(shape=(9), dtype=c_double)
        self.num_blades = 3
        self.init_root_pos = np.zeros(shape=(self.num_blades,3), dtype=c_float)
        self.init_root_orient = np.zeros(shape=(self.num_blades,9), dtype=c_double)

        #--------------------------------------
        # Structural Mesh
        #--------------------------------------
        self.num_mesh_pts = 1
        self.init_mesh_pos = np.zeros(shape=(self.num_mesh_pts,3), dtype=c_float)
        self.init_mesh_orient = np.zeros(shape=(self.num_mesh_pts,9), dtype=c_double)
        self.mesh_pt_to_blade_num = np.zeros(shape=(self.num_mesh_pts), dtype=c_int)

        #--------------------------------------
        # Environmental conditions
        #--------------------------------------
        self.gravity: float = 9.80665                # Gravitational acceleration (m/s^2)
        self.fluid_density: float = 1.225            # Air/fluid density (kg/m^3)
        self.kinematic_viscosity: float = 1.464E-05  # Kinematic viscosity (m^2/s)
        self.sound_speed: float = 335.               # Speed of sound (m/s)
        self.atmospheric_pressure: float = 103500.   # Atmospheric pressure (Pa)
        self.vapor_pressure: float = 1700.           # Vapor pressure (Pa)
        self.water_depth: float = 0.                 # Water depth (m)
        self.mean_sea_level_offset: float = 0.       # Mean sea level to still water level offset (m)

    def check_error(self) -> None:
        """Checks for and handles any errors from the Fortran library.

        Raises:
            RuntimeError: If a fatal error occurs in the Fortran code
        """
        # If the error status is 0, return
        if self.error_status_c.value == 0:
            return

        # Get the error level and error message
        error_level = self.error_levels.get(
            self.error_status_c.value,
            f"Unknown Error Level: {self.error_status_c.value}"
        )
        error_msg = self.error_message_c.value.decode('utf-8').strip()
        message = f"AeroDyn/InflowWind {error_level}: {error_msg}"

        # If the error level is fatal, call adi_end() and raise an error
        if self.error_status_c.value >= self.abort_error_level:
            try:
                self.adi_end()
            except Exception as e:
                message += f"\nAdditional error during cleanup: {e}"
            raise RuntimeError(message)
        else:
            print(message)

    def adi_preinit(self) -> None:
        """Pre-initializes the AeroDyn/InflowWind interface.

        Sets up initial arrays based on number of turbines.

        Raises:
            RuntimeError: If pre-initialization fails
        """
        self.ADI_C_PreInit(
            byref(c_int(self.num_turbines)),        # IN -> number of turbines
            byref(c_int(self.transpose_dcm)),       # IN -> transpose_dcm flag (0=false, 1=true)
            byref(c_int(self.point_load_output)),   # IN -> point_load_output flag (0=false, 1=true)
            byref(c_int(self.debug_level)),         # IN -> debug level (0=None to 4=Fatal)
            byref(self.error_status_c),             # OUT <- error status code
            self.error_message_c                    # OUT <- error message buffer
        )
        self.check_error()

    def adi_setuprotor(
        self,
        i_turbine: int,
        is_HAWT: int,
        turb_ref_pos: npt.NDArray[np.float32]
    ) -> None:
        """Sets up a single rotor with initial root/mesh information.

        Args:
            i_turbine: Turbine number
            is_HAWT: Flag indicating if turbine is horizontal axis (1) or not (0)
            turb_ref_pos: Reference position for the turbine [x, y, z]

        Raises:
            ValueError: If input arrays have incorrect dimensions
            RuntimeError: If rotor setup fails
        """
        self._init_num_mesh_pts = self.num_mesh_pts
        self._init_num_blades = self.num_blades
        turb_ref_pos_c = to_c_array(turb_ref_pos, c_float)

        # Validate the inputs
        self._validate_hub_root()
        self._validate_mesh()

        # Convert numpy arrays to C arrays with proper flattening
        init_arrays = self._prepare_init_arrays()

        self.ADI_C_SetupRotor(
            c_int(i_turbine),                     # IN -> current turbine number (0-based)
            c_int(is_HAWT),                       # IN -> 1: HAWT, 0: VAWT or cross-flow
            turb_ref_pos_c,                       # IN -> turbine reference position [x,y,z]
            init_arrays['hub_position_c'],        # IN -> initial hub position [x,y,z]
            init_arrays['hub_orientation_c'],     # IN -> initial hub orientation (DCM, flattened 3x3 matrix)
            init_arrays['nacelle_position_c'],    # IN -> initial nacelle position [x,y,z]
            init_arrays['nacelle_orientation_c'], # IN -> initial nacelle orientation (DCM, flattened 3x3 matrix)
            byref(c_int(self.num_blades)),        # IN -> number of blades
            init_arrays['root_position_c'],       # IN -> initial blade root positions (flattened array)
            init_arrays['root_orientation_c'],    # IN -> initial blade root orientations (flattened array)
            byref(c_int(self.num_mesh_pts)),      # IN -> number of structural mesh points
            init_arrays['mesh_position_c'],       # IN -> initial mesh point positions (flattened array)
            init_arrays['mesh_orientation_c'],    # IN -> initial mesh point orientations (flattened array)
            init_arrays['mesh_blade_num_c'],      # IN -> mapping of mesh points to blade numbers
            byref(self.error_status_c),           # OUT <- error status code
            self.error_message_c                  # OUT <- error message buffer
        )
        self.check_error()

    def adi_init(
        self,
        ad_input_string_array: List[str],
        ifw_input_string_array: List[str]
    ) -> None:
        """Initializes the AeroDyn/InflowWind simulation.

        Args:
            ad_input_string_array: List of strings containing AeroDyn input
            ifw_input_string_array: List of strings containing InflowWind input

        Raises:
            RuntimeError: If initialization fails
        """
        # Initialize channel counter
        self._num_channels_c = c_int(0)

        # Join input strings with null character separator
        ad_input_string, ad_input_string_length = to_c_string(ad_input_string_array)
        ifw_input_string, ifw_input_string_length = to_c_string(ifw_input_string_array)

        # Prepare output file paths
        output_file_root_name_c = create_string_buffer(
            self.output_root_name.ljust(self.DEFAULT_STRING_LENGTH).encode('utf-8')
        )
        vtk_output_dir_c = create_string_buffer(
            self.output_vtk_dir.ljust(self.DEFAULT_STRING_LENGTH).encode('utf-8')
        )

        # Convert VTK nacelle dimensions to C array
        vtk_nac_dimension_c = to_c_array(self.vtk_nacelle_dimension, c_float)

        self.ADI_C_Init(
            byref(c_int(self.aerodyn_inputs_passed_as_string)),  # IN -> AD input file is passed as string
            c_char_p(ad_input_string),                           # IN -> AD input file as string
            byref(c_int(ad_input_string_length)),                # IN -> AD input file string length
            byref(c_int(self.inflow_inputs_passed_as_string)),   # IN -> IfW input file is passed as string
            c_char_p(ifw_input_string),                          # IN -> IfW input file as string
            byref(c_int(ifw_input_string_length)),               # IN -> IfW input file string length
            output_file_root_name_c,                             # IN -> rootname for ADI file writing
            vtk_output_dir_c,                                    # IN -> directory for vtk output files
            byref(c_float(self.gravity)),                        # IN -> gravity
            byref(c_float(self.fluid_density)),                  # IN -> fluid density
            byref(c_float(self.kinematic_viscosity)),            # IN -> kinematic viscosity
            byref(c_float(self.sound_speed)),                    # IN -> speed of sound
            byref(c_float(self.atmospheric_pressure)),           # IN -> atmospheric pressure
            byref(c_float(self.vapor_pressure)),                 # IN -> vapor pressure
            byref(c_float(self.water_depth)),                    # IN -> water depth
            byref(c_float(self.mean_sea_level_offset)),          # IN -> MSL to SWL offset
            byref(c_int(self.interpolation_order)),              # IN -> interpolation order (1: linear, 2: quadratic)
            byref(c_double(self.dt)),                            # IN -> time step
            byref(c_double(self.t_max)),                         # IN -> maximum simulation time
            byref(c_int(self.store_hub_height_velocity)),        # IN -> store hub height velocity flag
            byref(c_int(self.write_vtk)),                        # IN -> write VTK flag
            byref(c_int(self.vtk_type)),                         # IN -> VTK write type
            byref(c_double(self.vtk_dt)),                        # IN -> VTK output time step
            vtk_nac_dimension_c,                                 # IN -> VTK nacelle dimensions
            byref(c_float(self.vtk_hub_radius)),                 # IN -> VTK hub radius
            byref(c_int(self.write_outputs)),                    # IN -> write outputs flag
            byref(c_double(self.output_timestep)),               # IN -> output time step
            byref(self._num_channels_c),                         # OUT <- number of channels
            self._channel_names_c,                               # OUT <- output channel names
            self._channel_units_c,                               # OUT <- output channel units
            byref(self.error_status_c),                          # OUT <- error status
            self.error_message_c                                 # OUT <- error message
        )
        self.check_error()

        # Store number of output channels
        self.num_channels = self._num_channels_c.value

    def adi_setrotormotion(
        self,
        i_turbine: int,
        hub: MotionData,
        nacelle: MotionData,
        root: MotionData,
        mesh: MotionData
    ) -> None:
        """Sets the rotor motion for simulation.
        Args:
            i_turbine: Turbine number
            hub: Hub motion data
            nacelle: Nacelle motion data
            root: Root motion data
            mesh: Mesh motion data

        Raises:
            ValueError: If motion data has incorrect dimensions
            RuntimeError: If setting rotor motion fails
        """
        # Validate inputs
        self._validate_motion_data(hub, "hub", single_pt=True)
        self._validate_motion_data(nacelle, "nacelle", single_pt=True)
        self._validate_motion_data(root, "root", num_pts=self._init_num_blades)
        self._validate_motion_data(mesh, "mesh", num_pts=self._init_num_mesh_pts)

        # Convert data to C arrays
        motion_arrays = self._prepare_motion_arrays(hub, nacelle, root, mesh)

        self.ADI_C_SetRotorMotion(
            c_int(i_turbine),                             # IN -> current turbine number (0-based)
            motion_arrays['hub_position_c'],              # IN -> hub positions
            motion_arrays['hub_orientation_c'],           # IN -> hub orientations
            motion_arrays['hub_velocity_c'],              # IN -> hub velocity [TVx, TVy, TVz, RVx, RVy, RVz]
            motion_arrays['hub_acceleration_c'],          # IN -> hub accelerations [TAx, TAy, TAz, RAx, RAy, RAz]
            motion_arrays['nacelle_position_c'],          # IN -> nacelle positions
            motion_arrays['nacelle_orientation_c'],       # IN -> nacelle orientations
            motion_arrays['nacelle_velocity_c'],          # IN -> nacelle velocity [TVx, TVy, TVz, RVx, RVy, RVz]
            motion_arrays['nacelle_acceleration_c'],      # IN -> nacelle accelerations [TAx, TAy, TAz, RAx, RAy, RAz]
            motion_arrays['root_position_c'],             # IN -> root positions
            motion_arrays['root_orientation_c'],          # IN -> root orientations (DCM)
            motion_arrays['root_velocity_c'],             # IN -> root velocities at desired positions
            motion_arrays['root_acceleration_c'],         # IN -> root accelerations at desired positions
            byref(c_int(self.num_mesh_pts)),              # IN -> number of attachment points expected (where motions are transferred into HD)
            motion_arrays['mesh_position_c'],             # IN -> mesh positions
            motion_arrays['mesh_orientation_c'],          # IN -> mesh orientations (DCM)
            motion_arrays['mesh_velocity_c'],             # IN -> mesh velocities at desired positions
            motion_arrays['mesh_acceleration_c'],         # IN -> mesh accelerations at desired positions
            byref(self.error_status_c),                   # OUT <- error status
            self.error_message_c                          # OUT <- error message
        )
        self.check_error()

    def adi_getrotorloads(
        self,
        i_turbine: int,
        mesh_forces_moments: npt.NDArray[np.float32],
        hub_height_velocity: Optional[npt.NDArray[np.float32]] = None
    ) -> None:
        """Get the rotor loads from the simulation.

        Args:
            i_turbine: Turbine number
            mesh_forces_moments: Array of mesh forces/moments [N, N-m]
            hub_height_velocity: Array of hub height velocity [m/s] (if storeHHVel is True)

        Raises:
            RuntimeError: If getting rotor loads fails
        """
        # Initialize C arrays for forces/moments and hub height velocity
        mesh_forces_moments_c = (c_float * (6 * self.num_mesh_pts))(0.)
        hub_height_velocity_c = (c_float * 3)(0.)

        # Call the C function to get rotor loads
        self.ADI_C_GetRotorLoads(
            c_int(i_turbine),                           # IN -> current turbine number
            byref(c_int(self.num_mesh_pts)),            # IN -> number of attachment points expected
            mesh_forces_moments_c,                      # OUT <- resulting forces/moments array
            hub_height_velocity_c,                      # OUT <- hub height velocity [Vx, Vy, Vz]
            byref(self.error_status_c),                 # OUT <- error status
            self.error_message_c                        # OUT <- error message
        )
        self.check_error()

        # Convert C arrays back to numpy arrays
        mesh_forces_moments[:, :] = np.reshape(mesh_forces_moments_c, (self.num_mesh_pts, 6))
        if hub_height_velocity is not None:
            hub_height_velocity[:] = np.reshape(hub_height_velocity_c, (3))

    def adi_getdiskavgvel(self, i_turbine: int, disk_avg_velocity: npt.NDArray[np.float32]) -> None:
        """Get the disk-averaged velocity.

        Args:
            i_turbine: Turbine number
            disk_avg_velocity: Array of disk-averaged velocities [m/s]

        Raises:
            RuntimeError: If getting disk average velocity fails
        """
        disk_avg_velocity_c = (c_float * 3)(0.)

        self.ADI_C_GetDiskAvgVel(
            c_int(i_turbine),                # IN -> current turbine number
            disk_avg_velocity_c,             # OUT <- disk-averaged velocity [Vx, Vy, Vz]
            byref(self.error_status_c),      # OUT <- error status
            self.error_message_c             # OUT <- error message
        )
        self.check_error()

        disk_avg_velocity[:] = np.reshape(disk_avg_velocity_c, (3))

    def adi_calcOutput(self, time: float, output_channel_values: npt.NDArray[np.float32]) -> None:
        """Calculate output values at the given time.

        Args:
            time: Current simulation time
            output_channel_values: Array to store calculated output values

        Raises:
            ValueError: If output_channel_values array has wrong size
            RuntimeError: If calculation fails
        """
        if output_channel_values.size != self.num_channels:
            raise ValueError(
                f"Output array must have size {self.num_channels}, "
                f"got {output_channel_values.size}"
            )

        output_channel_values_c = (c_float * self.num_channels)(0.)

        self.ADI_C_CalcOutput(
            byref(c_double(time)),           # IN -> current simulation time
            output_channel_values_c,         # OUT <- calculated output channel values
            byref(self.error_status_c),      # OUT <- error status
            self.error_message_c             # OUT <- error message
        )
        self.check_error()

        # Copy results back to numpy array
        output_channel_values[:] = np.reshape(output_channel_values_c, (self.num_channels))

    def adi_updateStates(self, time: float, time_next: float) -> None:
        """Update states from current time to next time step.

        Args:
            time: Current simulation time
            time_next: Next simulation time

        Raises:
            RuntimeError: If state update fails
        """
        self.ADI_C_UpdateStates(
            byref(c_double(time)),             # IN -> current simulation time, t
            byref(c_double(time_next)),        # IN -> next simulation time, t + dt
            byref(self.error_status_c),        # OUT <- error status
            self.error_message_c               # OUT <- error message
        )
        self.check_error()

    def adi_end(self) -> None:
        """Clean up and end the AeroDyn/InflowWind simulation.

        This method should be called when finishing the simulation to ensure
        proper cleanup of resources.

        Raises:
            RuntimeError: If cleanup fails
        """
        if not self.ended:
            self.ended = True
            self.ADI_C_End(
                byref(self.error_status_c),      # OUT <- error status
                self.error_message_c             # OUT <- error message
            )
            self.check_error()

    def get_output_info(self) -> Dict[str, List[str]]:
        """Get information about available output channels.

        Returns:
            Dictionary containing:
                - 'names': List of output channel names
                - 'units': List of corresponding units
        """
        return {
            'names': self.output_channel_names_c,
            'units': self.output_channel_units_c
        }

    @property
    def output_channel_names(self) -> List[str]:
        """Get the names of available output channels."""
        if not self._channel_names_c.value:
            return []
        names = self._channel_names_c.value.split()
        return [name.decode('utf-8') for name in names]

    @property
    def output_channel_units(self) -> List[str]:
        """Get the units of available output channels."""
        if not self._channel_units_c.value:
            return []
        units = self._channel_units_c.value.split()
        return [unit.decode('utf-8') for unit in units]

    def __enter__(self) -> 'AeroDynInflowLib':
        """Context manager entry.

        Returns:
            Self for use in with statement
        """
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Context manager exit.

        Ensures proper cleanup when used in a with statement.
        """
        self.adi_end()

    #-------------------------------------------------------------------------------
    # Private/internal methods
    #-------------------------------------------------------------------------------
    def _initialize_routines(self) -> None:
        """Initializes the Fortran routines i.e. the C-binding interfaces.

        Sets up the argument types and return types for all Fortran routines
        that will be called through the ctypes interface.
        """
        #--------------------------------------
        # ADI_C_PreInit
        #--------------------------------------
        self.ADI_C_PreInit.argtypes = [
            POINTER(c_int),                     # numTurbines
            POINTER(c_int),                     # transposeDCM
            POINTER(c_int),                     # pointLoadOutput
            POINTER(c_int),                     # debuglevel
            POINTER(c_int),                     # ErrStat_C
            POINTER(c_char)                     # ErrMsg_C
        ]
        self.ADI_C_PreInit.restype = c_int

        #--------------------------------------
        # ADI_C_SetupRotor
        #--------------------------------------
        self.ADI_C_SetupRotor.argtypes = [
            POINTER(c_int),                     # iturb
            POINTER(c_int),                     # isHAWT
            POINTER(c_float),                   # Turb_RefPos
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
            POINTER(c_int),                     # meshPtToBladeNum
            POINTER(c_int),                     # ErrStat_C
            POINTER(c_char)                     # ErrMsg_C
        ]
        self.ADI_C_SetupRotor.restype = c_int

        #--------------------------------------
        # ADI_C_Init
        #--------------------------------------
        self.ADI_C_Init.argtypes = [
            POINTER(c_int),                     # AD input file passed as string
            POINTER(c_char_p),                  # AD input file as string
            POINTER(c_int),                     # AD input file string length
            POINTER(c_int),                     # IfW input file passed as string
            POINTER(c_char_p),                  # IfW input file as string
            POINTER(c_int),                     # IfW input file string length
            POINTER(c_char),                    # OutRootName
            POINTER(c_char),                    # OutVTKdir
            POINTER(c_float),                   # gravity
            POINTER(c_float),                   # defFldDens
            POINTER(c_float),                   # defKinVisc
            POINTER(c_float),                   # defSpdSound
            POINTER(c_float),                   # defPatm
            POINTER(c_float),                   # defPvap
            POINTER(c_float),                   # WtrDpth
            POINTER(c_float),                   # MSL2SWL
            POINTER(c_int),                     # InterpOrder
            POINTER(c_double),                  # dt
            POINTER(c_double),                  # tmax
            POINTER(c_int),                     # storeHHVel
            POINTER(c_int),                     # WrVTK
            POINTER(c_int),                     # WrVTK_Type
            POINTER(c_double),                  # WrVTK_DT  -- 0 or negative to do every step
            POINTER(c_float),                   # VTKNacDim
            POINTER(c_float),                   # VTKHubRad
            POINTER(c_int),                     # wrOuts -- file format for writing outputs
            POINTER(c_double),                  # DT_Outs -- timestep for outputs to file
            POINTER(c_int),                     # number of channels
            POINTER(c_char),                    # output channel names
            POINTER(c_char),                    # output channel units
            POINTER(c_int),                     # ErrStat_C
            POINTER(c_char)                     # ErrMsg_C
        ]
        self.ADI_C_Init.restype = c_int

        #--------------------------------------
        # ADI_C_SetRotorMotion
        #--------------------------------------
        self.ADI_C_SetRotorMotion.argtypes = [
            POINTER(c_int),                     # iturb
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
        ]
        self.ADI_C_SetRotorMotion.restype = c_int

        #--------------------------------------
        # ADI_C_GetRotorLoads
        #--------------------------------------
        self.ADI_C_GetRotorLoads.argtypes = [
            POINTER(c_int),                     # iturb
            POINTER(c_int),                     # numMeshPts
            POINTER(c_float),                   # meshFrc -- mesh forces/moments in flat array of 6*numMeshPts
            POINTER(c_float),                   # hhVel -- wind speed at hub height in flat array of 3
            POINTER(c_int),                     # ErrStat_C
            POINTER(c_char)                     # ErrMsg_C
        ]
        self.ADI_C_GetRotorLoads.restype = c_int

        #--------------------------------------
        # ADI_C_GetDiskAvgVel
        #--------------------------------------
        self.ADI_C_GetDiskAvgVel.argtypes = [
            POINTER(c_int),                     # iturb
            POINTER(c_float),                   # Disk average vel vector
            POINTER(c_int),                     # ErrStat_C
            POINTER(c_char)                     # ErrMsg_C
        ]
        self.ADI_C_GetDiskAvgVel.restype = c_int

        #--------------------------------------
        # ADI_C_CalcOutput
        #--------------------------------------
        self.ADI_C_CalcOutput.argtypes = [
            POINTER(c_double),                  # Time_C
            POINTER(c_float),                   # Output Channel Values
            POINTER(c_int),                     # ErrStat_C
            POINTER(c_char)                     # ErrMsg_C
        ]
        self.ADI_C_CalcOutput.restype = c_int

        #--------------------------------------
        # ADI_C_UpdateStates
        #--------------------------------------
        self.ADI_C_UpdateStates.argtypes = [
            POINTER(c_double),                  # Time_C
            POINTER(c_double),                  # TimeNext_C
            POINTER(c_int),                     # ErrStat_C
            POINTER(c_char)                     # ErrMsg_C
        ]
        self.ADI_C_UpdateStates.restype = c_int

        #--------------------------------------
        # ADI_C_End
        #--------------------------------------
        self.ADI_C_End.argtypes = [
            POINTER(c_int),                     # ErrStat_C
            POINTER(c_char)                     # ErrMsg_C
        ]
        self.ADI_C_End.restype = c_int

    def _prepare_init_arrays(self) -> Dict[str, Any]:
        """Prepares C-compatible arrays for initialization.

        Returns:
            Dictionary containing all prepared C arrays
        """
        return {
            # Hub data
            'hub_position_c': to_c_array(self.init_hub_pos, c_float),
            'hub_orientation_c': to_c_array(self.init_hub_orient, c_double),
            # Nacelle data
            'nacelle_position_c': to_c_array(self.init_nacelle_pos, c_float),
            'nacelle_orientation_c': to_c_array(self.init_nacelle_orient, c_double),
            # Root data (per blade)
            'root_position_c': flatten_array(
                self._init_num_blades,
                self.num_blades,
                self.init_root_pos,
                'RootPos',
                elements_per_item=3,
                c_type=c_float
            ),
            'root_orientation_c': flatten_array(
                self._init_num_blades,
                self.num_blades,
                self.init_root_orient,
                'RootOrient',
                elements_per_item=9,
                c_type=c_double
            ),
            # Mesh data (structural points)
            'mesh_position_c': flatten_array(
                self._init_num_mesh_pts,
                self.num_mesh_pts,
                self.init_mesh_pos,
                'MeshPos',
                elements_per_item=3,
                c_type=c_float
            ),
            'mesh_orientation_c': flatten_array(
                self._init_num_mesh_pts,
                self.num_mesh_pts,
                self.init_mesh_orient,
                'MeshOrient',
                elements_per_item=9,
                c_type=c_double
            ),
            # Blade mapping for mesh points
            'mesh_blade_num_c': to_c_array(self.mesh_pt_to_blade_num, c_int)
        }

    def _validate_hub_root(self) -> None:
        """Validates hub and root configurations.

        Raises:
            ValueError: If any validation fails
        """
        if self.num_blades < 1:
            raise ValueError("Number of blades must be at least 1")
        if self.init_root_pos.shape != (self.num_blades, 3):
            raise ValueError(
                f"Root positions must have shape ({self.num_blades}, 3), "
                f"got {self.init_root_pos.shape}"
            )
        if self.init_root_orient.shape != (self.num_blades, 9):
            raise ValueError(
                f"Root orientations must have shape ({self.num_blades}, 9), "
                f"got {self.init_root_orient.shape}"
            )
        if self.init_hub_pos.shape != (3,):
            raise ValueError(
                f"Hub position must have shape (3,), got {self.init_hub_pos.shape}"
            )
        if self.init_hub_orient.shape != (9,):
            raise ValueError(
                f"Hub orientation must have shape (9,), got {self.init_hub_orient.shape}"
            )

    def _validate_mesh(self) -> None:
        """Validates mesh configurations.

        Raises:
            ValueError: If any validation fails
        """
        if self.init_mesh_pos.shape != (self.num_mesh_pts, 3):
            raise ValueError(
                f"Mesh positions must have shape ({self.num_mesh_pts}, 3), "
                f"got {self.init_mesh_pos.shape}"
            )
        if self.init_mesh_orient.shape != (self.num_mesh_pts, 9):
            raise ValueError(
                f"Mesh orientations must have shape ({self.num_mesh_pts}, 9), "
                f"got {self.init_mesh_orient.shape}"
            )
        if self.init_mesh_pos.shape[0] != self.init_mesh_orient.shape[0]:
            raise ValueError(
                "Inconsistent number of mesh points between position and orientation arrays"
            )

    def _validate_motion_data(
        self,
        motion: MotionData,
        name: str,
        *,
        single_pt: bool = False,
        num_pts: Optional[int] = None
    ) -> None:
        """Validates motion data dimensions.

        Args:
            motion: Motion data to validate
            name: Name of the component for error messages
            single_pt: Whether this is single-point data
            num_pts: Number of points expected (if not single_pt)

        Raises:
            ValueError: If dimensions are incorrect
        """
        expected_shape = (1,3) if single_pt else (num_pts, 3)
        expected_orient_shape = (1,9) if single_pt else (num_pts, 9)
        expected_vel_shape = (1,6) if single_pt else (num_pts, 6)

        if motion.position.shape != expected_shape:
            raise ValueError(
                f"{name} position must have shape {expected_shape}, "
                f"got {motion.position.shape}"
            )

        if motion.orientation.shape != expected_orient_shape:
            raise ValueError(
                f"{name} orientation must have shape {expected_orient_shape}, "
                f"got {motion.orientation.shape}"
            )

        if motion.velocity.shape != expected_vel_shape:
            raise ValueError(
                f"{name} velocity must have shape {expected_vel_shape}, "
                f"got {motion.velocity.shape}"
            )

        if motion.acceleration.shape != expected_vel_shape:
            raise ValueError(
                f"{name} acceleration must have shape {expected_vel_shape}, "
                f"got {motion.acceleration.shape}"
            )

    def _prepare_motion_arrays(
        self,
        hub: MotionData,
        nacelle: MotionData,
        root: MotionData,
        mesh: MotionData
    ) -> Dict[str, Any]:
        """Prepares C-compatible arrays for motion data.

        Args:
            hub: Hub motion data
            nacelle: Nacelle motion data
            root: Root motion data
            mesh: Mesh motion data

        Returns:
            Dictionary containing all prepared C arrays
        """
        return {
            # Hub data
            'hub_position_c': to_c_array(hub.position, c_float),
            'hub_orientation_c': to_c_array(hub.orientation, c_double),
            'hub_velocity_c': to_c_array(hub.velocity, c_float),
            'hub_acceleration_c': to_c_array(hub.acceleration, c_float),
            # Nacelle data
            'nacelle_position_c': to_c_array(nacelle.position, c_float),
            'nacelle_orientation_c': to_c_array(nacelle.orientation, c_double),
            'nacelle_velocity_c': to_c_array(nacelle.velocity, c_float),
            'nacelle_acceleration_c': to_c_array(nacelle.acceleration, c_float),
            # Root data (per blade)
            'root_position_c': flatten_array(
                self._init_num_blades,
                self.num_blades,
                root.position,
                'RootPos',
                elements_per_item=3,
                c_type=c_float
            ),
            'root_orientation_c': flatten_array(
                self._init_num_blades,
                self.num_blades,
                root.orientation,
                'RootOrient',
                elements_per_item=9,
                c_type=c_double
            ),
            'root_velocity_c': flatten_array(
                self._init_num_blades,
                self.num_blades,
                root.velocity,
                'RootVel',
                elements_per_item=6,
                c_type=c_float
            ),
            'root_acceleration_c': flatten_array(
                self._init_num_blades,
                self.num_blades,
                root.acceleration,
                'RootAcc',
                elements_per_item=6,
                c_type=c_float
            ),
            # Mesh data (structural points)
            'mesh_position_c': flatten_array(
                self._init_num_mesh_pts,
                self.num_mesh_pts,
                mesh.position,
                'MeshPos',
                elements_per_item=3,
                c_type=c_float
            ),
            'mesh_orientation_c': flatten_array(
                self._init_num_mesh_pts,
                self.num_mesh_pts,
                mesh.orientation,
                'MeshOrient',
                elements_per_item=9,
                c_type=c_double
            ),
            'mesh_velocity_c': flatten_array(
                self._init_num_mesh_pts,
                self.num_mesh_pts,
                mesh.velocity,
                'MeshVel',
                elements_per_item=6,
                c_type=c_float
            ),
            'mesh_acceleration_c': flatten_array(
                self._init_num_mesh_pts,
                self.num_mesh_pts,
                mesh.acceleration,
                'MeshAcc',
                elements_per_item=6,
                c_type=c_float
            )
        }

#-------------------------------------------------------------------------------
# Write output channels to a file
#-------------------------------------------------------------------------------
class WriteOutChans:
    """A helper class for writing output channels to file.

    This class writes simulation output channels to a text file in a tabular format.
    It's used for regression testing to mirror the output from AD15 and InflowWind
    from an OpenFAST simulation, and is valuable for debugging interfaces to the
    ADI_C_Binding library.

    When coupled to another code, this data would typically be passed back for inclusion
    in any output files there.

    Attributes:
        filename: Name of the output file
        opened: Boolean flag indicating if the output file is currently open
    """

    def __init__(self, filename: str, channel_names: List[str], channel_units: List[str]) -> None:
        channel_names.insert(0, 'Time')             # add time index header
        channel_units.insert(0, '(s)')              # add time index unit
        self.out_file = open(filename, 'wt')        # open output file and write header info
        # write file header
        t_string = datetime.now()
        dt_string = datetime.today()
        self.out_file.write(f"## This file was generated by AeroDyn_Inflow_Driver on {dt_string.strftime('%b-%d-%Y')} at {t_string.strftime('%H:%M:%S')}\n")
        self.out_file.write(f"## This file contains output channels requested from the OutList section of the AD15 and IfW input files")
        self.out_file.write(f"{filename}\n")
        self.out_file.write("#\n")
        self.out_file.write("#\n")
        self.out_file.write("#\n")
        self.out_file.write("#\n")
        l = len(channel_names)
        f_string = "{:^15s}"+"   {:^20s}  "*(l-1)
        self.out_file.write(f_string.format(*channel_names) + '\n')
        self.out_file.write(f_string.format(*channel_units) + '\n')
        self.opened = True

    def write(self, channel_data: npt.NDArray) -> None:
        """Write the channel data to the output file.

        Args:
            channel_data: Array of channel data with shape (n_timesteps, n_channels)
                          First column should be time values
        """
        time_format = "{:10.4f}"
        data_format = "{:25.7f}" * (channel_data.shape[1] - 1)
        format_str = time_format + data_format

        rows = [format_str.format(*row) for row in channel_data]
        self.out_file.write("\n".join(rows) + "\n")
        self.out_file.flush()

    def end(self) -> None:
        """Close the output file."""
        if self.opened:
            self.out_file.close()
            self.opened = False

    def __enter__(self) -> 'WriteOutChans':
        """Support for context manager protocol."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Ensure file is closed when exiting context."""
        self.end()

#-------------------------------------------------------------------------------
# Generate a debug file
#-------------------------------------------------------------------------------
class DriverDbg:
    """
    A helper class for debugging the ADI interface. This class writes out all the
    input positions/orientations, velocities, accelerations, and the resulting
    forces and moments at each input mesh point. If functioning correctly, this
    will be identical to the corresponding values in the AeroDyn/InflowWind output
    channels.

    NOTE: This may not output everything in the interface as updates have been made
    since writing this, but this routine was not updated accordingly.
    """

    def __init__(self, filename: str, num_mesh_pts: int) -> None:
        """Initializes the debugging class and open the output file."""
        self.filename = filename
        self.num_mesh_pts = num_mesh_pts
        self.opened = True

        with open(filename, 'wt') as self.debug_file:
            self._write_header()

        self.debug_file = open(filename, 'at') # switch to append mode

    def _write_header(self) -> None:
        """Writes the header information to the debug file."""
        # Build header components
        timestamp = datetime.now().strftime('%b-%d-%Y %H:%M:%S')
        header_lines = [
            f"## This file was generated by adi_c_lib on {timestamp}",
            f"## This file contains the resulting forces/moments at each of {self.num_mesh_pts} mesh points passed into the adi_c_lib",
            "#",
            "#",
            "#",
            "#"
        ]

        # Write column headers
        column_names = ["Time"]
        column_units = ["(s)"]
        for i in range(1, self.num_mesh_pts + 1):
            prefix = f"N{i:04d}_"
            # Position columns
            for suffix in ["x", "y", "z"]:
                column_names.append(f"{prefix}{suffix}")
                column_units.append("(m)")
            # Velocity columns
            for suffix in ["Vx", "Vy", "Vz"]:
                column_names.append(f"{prefix}{suffix}")
                column_units.append("(m/s)")
            # Angular velocity columns
            for suffix in ["RVx", "RVy", "RVz"]:
                column_names.append(f"{prefix}{suffix}")
                column_units.append("(rad/s)")
            # Acceleration columns
            for suffix in ["Ax", "Ay", "Az"]:
                column_names.append(f"{prefix}{suffix}")
                column_units.append("(m/s^2)")
            # Angular acceleration columns
            for suffix in ["RAx", "RAy", "RAz"]:
                column_names.append(f"{prefix}{suffix}")
                column_units.append("(rad/s^2)")
            # Force columns
            for suffix in ["Fx", "Fy", "Fz"]:
                column_names.append(f"{prefix}{suffix}")
                column_units.append("(N)")
            # Moment columns
            for suffix in ["Mx", "My", "Mz"]:
                column_names.append(f"{prefix}{suffix}")
                column_units.append("(N-m)")
        # Disk average velocity columns
        for suffix in ["DskAvgVx", "DskAvgVy", "DskAvgVz"]:
            column_names.append(suffix)
            column_units.append("(m/s)")

        f_string = "{:^25s}"
        header_lines.append("".join([f_string.format(name) for name in column_names]))
        header_lines.append("".join([f_string.format(unit) for name, unit in zip(column_names, column_units)]))

        self.debug_file.write("\n".join(header_lines) + "\n")

    def write(
        self,
        t: float,
        mesh_position: npt.NDArray,
        mesh_velocity: npt.NDArray,
        mesh_acceleration: npt.NDArray,
        mesh_forces_moments: npt.NDArray,
        disk_avg_velocity: npt.NDArray
    ) -> None:
        """Writes the current state to the debug file."""
        row_data = [f"{t:10.4f}"]

        for i in range(self.num_mesh_pts):
            row_data.extend([f"{val:25.7e}" for val in mesh_position[i, :]])
            row_data.extend([f"{val:25.7e}" for val in mesh_velocity[i, :]])
            row_data.extend([f"{val:25.7e}" for val in mesh_acceleration[i, :]])
            row_data.extend([f"{val:25.7e}" for val in mesh_forces_moments[i, :]])
        row_data.extend([f"{val:25.7e}" for val in disk_avg_velocity[:]])

        self.debug_file.write("".join(row_data) + "\n")
        self.debug_file.flush()

    def end(self) -> None:
        """Closes the debug file."""
        if self.opened:
            self.debug_file.close()
            self.opened = False
