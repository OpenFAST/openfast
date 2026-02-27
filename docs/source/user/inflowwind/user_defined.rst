.. _ifw_user_defined:

User-Defined Wind Fields
=========================

This section explains how to implement custom wind fields in InflowWind using ``WindType = 6``.

Overview
--------

The user-defined wind field feature allows developers to implement custom wind models by:

1. Defining a data structure to hold wind field parameters
2. Initializing that data structure from input files or parameters  
3. Implementing a function to return wind velocities at any position and time

This is useful for:

- Analytical wind models (e.g., vortex, wake models)
- Custom wind profiles not available in standard formats
- Coupling to external wind solvers
- Real-time wind measurements from sensors
- Research and development of new wind field representations

.. important::
   After modifying the registry files (``.txt`` files), you must rebuild the project 
   to regenerate the type definition files (``*_Types.f90``). The modifications to the 
   ``.txt`` files define the extended data structures, but they won't be available 
   until after regeneration.

Implementation Steps
--------------------

Step 1: Define Data Structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Edit ``modules/inflowwind/src/IfW_FlowField.txt`` and add fields to ``UserFieldType``:

.. code-block:: text

   typedef  ^              UserFieldType          ReKi                RefHeight           -     -         -     "reference height; used to center the wind"                   meters
   typedef  ^              ^                      IntKi               NumDataLines        -     0         -     "number of data lines (for time-varying user wind)"          -
   typedef  ^              ^                      DbKi                DTime               :     -         -     "time array for user-defined wind"                           seconds
   typedef  ^              ^                      ReKi                Data                ::    -         -     "user-defined wind data array [NumDataLines, NumDataColumns]" -
   typedef  ^              ^                      CHARACTER(1024)     FileName            -     -         -     "name of user wind file (if applicable)"                     -

Add any custom fields needed for your wind model implementation.

Step 2: Define Initialization Inputs  
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Edit ``modules/inflowwind/src/InflowWind_IO.txt`` and add fields to ``User_InitInputType``:

.. code-block:: text

   typedef  ^              User_InitInputType    CHARACTER(1024)         WindFileName            -     -     -     "name of file containing user-defined wind data (if applicable)" -
   typedef  ^              ^                     ReKi                    RefHt                   -     -     -     "reference height for user wind field"                        meters
   typedef  ^              ^                     IntKi                   NumDataColumns          -     0     -     "number of data columns in user wind file (if applicable)"    -

Add any parameters needed to initialize your wind model.

Step 3: Regenerate Type Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After modifying the registry files, rebuild the project to regenerate type definitions:

.. code-block:: bash

   cd build
   cmake .. -DGENERATE_TYPES=ON
   make

The build process automatically regenerates the ``*_Types.f90`` files from the ``.txt`` registry files.

Step 4: Implement Initialization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Edit ``modules/inflowwind/src/InflowWind_IO.f90`` and implement ``IfW_User_Init()``:

.. code-block:: fortran

   subroutine IfW_User_Init(InitInp, SumFileUnit, UF, FileDat, ErrStat, ErrMsg)
      ! Initialize UF%RefHeight, read data files, allocate arrays
      ! Set FileDat metadata (wind type, time range, spatial extent, etc.)
      ! Write summary information to SumFileUnit if > 0
   end subroutine

This routine:

- Reads any necessary input files specified in ``InitInp``
- Allocates and populates the ``UserFieldType`` (``UF``) data structure
- Sets appropriate metadata in the ``WindFileDat`` structure
- Writes initialization information to the summary file

Step 5: Implement Velocity Function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Edit ``modules/inflowwind/src/IfW_FlowField.f90`` and implement ``UserField_GetVel()``:

.. code-block:: fortran

   subroutine UserField_GetVel(UF, Time, Position, Velocity, ErrStat, ErrMsg)
      ! Use UF data to compute velocity at Position and Time
      ! Position(1) = X, Position(2) = Y, Position(3) = Z (meters)
      ! Return Velocity(1) = U, Velocity(2) = V, Velocity(3) = W (m/s)
   end subroutine

This function is called for each position where wind velocities are needed during simulation.

Coordinate Systems
------------------

Input Coordinates (Position)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **X**: Downstream direction (after rotation applied by InflowWind)
- **Y**: Lateral/crosswind direction  
- **Z**: Vertical direction (measured from ground, Z=0 is ground level)
- **Units**: meters

Output Velocities (Velocity)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **U**: Velocity component along X (positive = downwind)
- **V**: Velocity component along Y (positive = to the left when looking downwind)
- **W**: Velocity component along Z (positive = upward)
- **Units**: m/s

.. note::
   InflowWind handles the rotation between global coordinates and wind coordinates. 
   Your implementation should work in the wind coordinate system where X is aligned 
   with the mean wind direction.

Example Implementation
----------------------

Power-Law Wind Profile
~~~~~~~~~~~~~~~~~~~~~~

This example implements a simple power-law wind profile.

**Velocity Function** (in ``IfW_FlowField.f90``):

.. code-block:: fortran

   subroutine UserField_GetVel(UF, Time, Position, Velocity, ErrStat, ErrMsg)
      type(UserFieldType), intent(in)     :: UF
      real(DbKi), intent(in)              :: Time
      real(ReKi), intent(in)              :: Position(3)
      real(ReKi), intent(out)             :: Velocity(3)
      integer(IntKi), intent(out)         :: ErrStat
      character(*), intent(out)           :: ErrMsg
      
      real(ReKi)                          :: RefSpeed, Exponent, Height
      
      ErrStat = ErrID_None
      ErrMsg = ""
      
      ! Get reference speed and exponent from UF%Data
      RefSpeed = UF%Data(1, 1)  ! Reference wind speed (m/s)
      Exponent = UF%Data(1, 2)  ! Power law exponent
      Height = Position(3)       ! Height above ground
      
      ! Apply power law: U(z) = Uref * (z/zref)^alpha
      if (Height > 0.0_ReKi) then
         Velocity(1) = RefSpeed * (Height / UF%RefHeight)**Exponent
         Velocity(2) = 0.0_ReKi  ! No lateral wind
         Velocity(3) = 0.0_ReKi  ! No vertical wind
      else
         Velocity = 0.0_ReKi  ! Below ground
      end if
      
   end subroutine

**Initialization** (in ``InflowWind_IO.f90``):

.. code-block:: fortran

   subroutine IfW_User_Init(InitInp, SumFileUnit, UF, FileDat, ErrStat, ErrMsg)
      ! ... (declarations)
      
      ErrStat = ErrID_None
      ErrMsg = ""
      
      ! Set reference height
      UF%RefHeight = InitInp%RefHt
      
      ! Allocate data array for [RefSpeed, Exponent]
      UF%NumDataLines = 1
      call AllocAry(UF%Data, 1, 2, 'User wind data', TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return
      
      ! Set values (could read from file instead)
      UF%Data(1, 1) = 10.0_ReKi  ! 10 m/s reference speed
      UF%Data(1, 2) = 0.2_ReKi   ! Power law exponent
      
      ! Set metadata
      FileDat%WindType = 6
      FileDat%RefHt = UF%RefHeight
      FileDat%MWS = UF%Data(1, 1)
      FileDat%RefHt_Set = .true.
      ! ... (set other FileDat fields as needed)
      
   end subroutine

Common Use Cases
----------------

Steady Analytical Wind Field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Define wind as a function of position only (ignore ``Time`` parameter).

**Example**: Logarithmic wind profile, vortex wind field, uniform flow with shear.

Time-Varying Wind Field  
~~~~~~~~~~~~~~~~~~~~~~~

Store wind data in time series arrays and interpolate based on ``Time`` parameter.

**Example**: Measured wind data, prescribed wind transients, wake models.

Wind from External Solver
~~~~~~~~~~~~~~~~~~~~~~~~~~

Call external functions to get instantaneous wind fields.

**Example**: CFD coupling, external wake models, prescribed turbulence.

Real-Time Sensor Data
~~~~~~~~~~~~~~~~~~~~~

Load measured wind data from sensors and interpolate spatially/temporally.

**Example**: LIDAR measurements, met mast data, field measurements.

Limitations and Considerations
-------------------------------

Current Limitations
~~~~~~~~~~~~~~~~~~~

1. **No Acceleration Support**: User-defined wind fields do not currently support 
   acceleration calculations needed by some modules (e.g., MHK turbines).

2. **No Persistence**: Data must be recalculated if simulation is restarted.

Performance Considerations
~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``UserField_GetVel()`` is called for every point at every time step
- Implement efficiently; pre-compute values in ``IfW_User_Init()`` when possible
- Consider caching or interpolation strategies for complex calculations

Error Handling
~~~~~~~~~~~~~~

- Always validate input parameters in ``IfW_User_Init()``
- Check array bounds in ``UserField_GetVel()``
- Verify Position and Time values are within valid ranges
- Use ``SetErrStat()`` to report errors appropriately

Best Practices
--------------

1. **Start Simple**: Begin with analytical models before implementing complex wind fields

2. **Document Thoroughly**: Add detailed comments explaining your implementation and any file formats

3. **Use SI Units**: Always use meters, seconds, and m/s

4. **Pre-compute**: Calculate as much as possible during initialization rather than runtime

5. **Validate**: Test with known analytical solutions before using in production

6. **Handle Boundaries**: Implement appropriate behavior for points outside valid domain

7. **Report Metadata**: Properly populate ``WindFileDat`` with time range, spatial extent, etc.

File Locations
--------------

==========================================  ==========================================================
File                                        Purpose
==========================================  ==========================================================
``modules/inflowwind/src/``                 Source code directory
``IfW_FlowField.txt``                       Type definitions for flow field data structures
``InflowWind_IO.txt``                       Type definitions for initialization inputs
``IfW_FlowField.f90``                       Flow field implementation (``UserField_GetVel()``)
``InflowWind_IO.f90``                       Initialization implementation (``IfW_User_Init()``)
``IfW_FlowField_Types.f90``                 Auto-generated type definitions (regenerated from .txt)
``InflowWind_IO_Types.f90``                 Auto-generated type definitions (regenerated from .txt)
==========================================  ==========================================================

Additional Resources
--------------------

- See the original :download:`InflowWind Manual <InflowWind_Manual.pdf>` for general InflowWind information
- Review existing wind field implementations (Uniform, Grid3D) in the source code for reference
- Check :ref:`ifw_appendix` for example input files
- Refer to NWTC Library documentation for array allocation and error handling utilities
- See :ref:`ifw_angles` for information about wind coordinate systems and rotations
