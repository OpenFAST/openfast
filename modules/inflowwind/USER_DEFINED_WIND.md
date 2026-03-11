# User-Defined Wind Field Implementation Guide

This guide explains how to implement custom wind fields in InflowWind (WindType = 6).

## Overview

The user-defined wind field feature allows you to implement custom wind models by:
1. Defining a data structure to hold your wind field parameters
2. Initializing that data structure from input files or parameters
3. Implementing a function to return wind velocities at any position and time

**Important**: After modifying the registry files (`.txt` files), you must rebuild the project to regenerate the type definition files (`*_Types.f90`). The modifications to the `.txt` files in this commit define the extended data structures, but they won't be available until after regeneration.

## Quick Start

### Step 1: Define Your Data Structure

Edit `modules/inflowwind/src/IfW_FlowField.txt` and add fields to `UserFieldType`:

```
typedef  ^              UserFieldType          ReKi                RefHeight           -     -         -     "reference height; used to center the wind"                   meters
typedef  ^              ^                      IntKi               NumDataLines        -     0         -     "number of data lines (for time-varying user wind)"          -
typedef  ^              ^                      DbKi                DTime               :     -         -     "time array for user-defined wind"                           seconds
typedef  ^              ^                      ReKi                Data                ::    -         -     "user-defined wind data array [NumDataLines, NumDataColumns]" -
typedef  ^              ^                      CHARACTER(1024)     FileName            -     -         -     "name of user wind file (if applicable)"                     -
# Add your custom fields here
```

### Step 2: Define Initialization Inputs

Edit `modules/inflowwind/src/InflowWind_IO.txt` and add fields to `User_InitInputType`:

```
typedef  ^              User_InitInputType    CHARACTER(1024)         WindFileName            -     -     -     "name of file containing user-defined wind data (if applicable)" -
typedef  ^              ^                     ReKi                    RefHt                   -     -     -     "reference height for user wind field"                        meters
typedef  ^              ^                     IntKi                   NumDataColumns          -     0     -     "number of data columns in user wind file (if applicable)"    -
# Add your custom initialization parameters here
```

### Step 3: Regenerate Type Files

After modifying the registry files, regenerate the type files:

```bash
cd modules/inflowwind/src
# Run the registry generator (typically done during build)
# or rebuild the project which will regenerate types automatically
```

### Step 4: Implement Initialization

Edit `modules/inflowwind/src/InflowWind_IO.f90` and implement `IfW_User_Init()`:

```fortran
subroutine IfW_User_Init(InitInp, SumFileUnit, UF, FileDat, ErrStat, ErrMsg)
   ! Read input files
   ! Allocate arrays in UF
   ! Populate UF with your wind field data
   ! Set FileDat metadata
end subroutine
```

### Step 5: Implement Wind Velocity Function

Edit `modules/inflowwind/src/IfW_FlowField.f90` and implement `UserField_GetVel()`:

```fortran
subroutine UserField_GetVel(UF, Time, Position, Velocity, ErrStat, ErrMsg)
   ! Use UF data to compute velocity at Position and Time
   ! Position(1) = X, Position(2) = Y, Position(3) = Z
   ! Return Velocity(1) = U, Velocity(2) = V, Velocity(3) = W
end subroutine
```

## Example Implementation: Power-Law Wind Profile

Here's a simple example implementing a power-law wind profile:

### In UserField_GetVel():

```fortran
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
```

### In IfW_User_Init():

```fortran
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
   ! ... (set other FileDat fields)
   
end subroutine
```

## Coordinate Systems

### Input Coordinates (Position)
- **X**: Downstream direction (after rotation applied by InflowWind)
- **Y**: Lateral/crosswind direction  
- **Z**: Vertical direction (measured from ground, Z=0 is ground level)
- Units: meters

### Output Velocities (Velocity)
- **U**: Velocity component along X (positive = downwind)
- **V**: Velocity component along Y (positive = to the left when looking downwind)
- **W**: Velocity component along Z (positive = upward)
- Units: m/s

Note: InflowWind handles the rotation between global coordinates and wind coordinates. Your implementation should work in the wind coordinate system where X is aligned with the mean wind direction.

## Common Use Cases

### 1. Steady Analytical Wind Field
Define wind as a function of position only (ignore Time parameter).

### 2. Time-Varying Wind Field
Store wind data in time series arrays and interpolate based on Time parameter.

### 3. Wind from External Solver
Call external functions or read shared memory to get instantaneous wind fields.

### 4. Measured Wind Data
Load measured wind data from sensors and interpolate spatially/temporally.

## Limitations

1. **No Acceleration Support**: User-defined wind fields do not currently support acceleration calculations needed by some modules.

2. **Performance**: The UserField_GetVel() function is called for every point at every time step, so it should be efficient.

3. **No Built-in Interpolation**: You must implement any necessary spatial or temporal interpolation.

## Tips and Best Practices

1. **Error Handling**: Always check array bounds and validity of Position and Time values.

2. **Efficiency**: Pre-compute values in IfW_User_Init() rather than in UserField_GetVel().

3. **Testing**: Start with simple analytical models before implementing complex wind fields.

4. **Documentation**: Document your implementation and any file formats in comments.

5. **Units**: Always use SI units (meters, seconds, m/s).

## File Locations

- **Type Definitions**: `modules/inflowwind/src/IfW_FlowField.txt`, `InflowWind_IO.txt`
- **Initialization**: `modules/inflowwind/src/InflowWind_IO.f90` (IfW_User_Init)
- **Velocity Function**: `modules/inflowwind/src/IfW_FlowField.f90` (UserField_GetVel)

## Building

After modifying any `.txt` registry files, rebuild the project:

```bash
cd build
cmake .. -DGENERATE_TYPES=ON
make
```

The build process will automatically regenerate the type files from the `.txt` registry files.

## Further Information

- See existing wind field implementations (Uniform, Grid3D) for reference
- Check InflowWind documentation for information about wind coordinate systems
- Refer to NWTC Library documentation for array allocation and error handling utilities
