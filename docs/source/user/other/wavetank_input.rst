.. _WaveTank-Input:

Input File
----------

This document describes the WaveTank configuration input file (``wavetankconfig_input.txt``) used to set up and run the WaveTank model for marine hydrokinetic (MHK) turbine testing.

- The file is read by the WaveTank library during initialization.

Conventions and Units
---------------------

- SI units are used throughout: m, s, kg, N, Pa.
- Angles are in degrees unless otherwise specified.
- Rotational speed is in rpm where noted.
- Positions and heights are referenced to Mean Sea Level (MSL) unless otherwise noted.
- Input files for modules may be relative or absolute paths.

File and Simulation Control
---------------------------

OutRootName (string)
  Root name used when writing summary or other files.
  Example: ``FRM1Q_Floating_tank_test``.

DT (s)
  Nominal timestep for WaveTank internal scheduling. Currently unused/reserved.

TMax (s)
  Maximum simulation time for WaveTank internal scheduling. Currently unused/reserved.

MHK (switch)
  MHK turbine type switch:

  - 0: Not an MHK turbine
  - 1: Fixed MHK turbine
  - 2: Floating MHK turbine

  Only the floating option (2) is supported at present.

InterpOrd (-)
  Interpolation order for internal data interpolation. Currently unused/reserved.

DebugLevel (switch)
  Controls logging and visualization detail:

  - 0: none
  - 1: I/O summary
  - 2: + positions/orientations passed
  - 3: + input file
  - 4: + all meshes

.. note::
   Parameters marked “unused” are reserved for future development and are currently ignored by the code path.

Froude Scaling (disabled)
-------------------------

The following parameters may appear but are typically commented out. Froude scaling is not complete in the current code. Do not use unless explicitly enabled.

ScaleFact (-)
  Froude scaling factor λ = (full-size dimension) / (model-size dimension). Expected > 1 for scale-model testing.

DensFact (-)
  Density ratio ρ_full / ρ_model, used with Froude scaling of forces/moments.

Environment
-----------

Gravity (m/s^2)
  Gravitational acceleration.

WtrDens (kg/m^3)
  Water (working fluid) density.

WtrVisc (m^2/s)
  Kinematic viscosity of the working fluid.

SpdSound (m/s)
  Speed of sound in the working fluid.

Patm (Pa)
  Atmospheric pressure. Used for cavitation checks.

Pvap (Pa)
  Vapor pressure of the working fluid. Used for cavitation checks.

WtrDpth (m)
  Water depth.

MSL2SWL (m)
  Offset between still-water level (SWL) and mean sea level (MSL); positive upward.

Sea State
---------

SS_InputFile (string)
  Path to SeaState input file defining wave conditions. Ensure path is valid relative to the run directory or use an absolute path.

WaveTimeShift (s)
  Time shift applied to the SeaState wave time series to adjust phase and match tank conditions.

MoorDyn
-------

MD_InputFile (string)
  Path to MoorDyn input file defining mooring system properties and connections.

AeroDyn and InflowWind
----------------------

AD_InputFile (string)
  Path to AeroDyn input file defining aerodynamic model configuration (used for hydro/aero coupling as applicable in MHK context).

IfW_InputFile (string)
  Path to InflowWind input file defining inflow conditions for the rotor (e.g., currents or wind, depending on model setup).

Turbine Geometry and Reference Frames
-------------------------------------

NumBl (-)
  Number of blades on the rotor.

HubRad (m)
  Distance from the rotor apex to the blade root.

PreCone (deg)
  Blade cone angle.

OverHang (m)
  Distance from the yaw axis (tower centerline) to the rotor apex. Negative values indicate rotor apex aft of the yaw axis under the model’s convention.

ShftTilt (deg)
  Rotor shaft tilt angle.

Twr2Shft (m)
  Vertical distance from tower-top to the rotor shaft center (nacelle center). Negative values are below tower-top.

TowerHt (m)
  Height of the tower relative to MSL. Tower is vertically aligned with ``TowerBsPt`` (sloped towers not supported).

TowerBsPt (m, m, m)
  Tower base location relative to the platform reference position in x and y, and relative to MSL in z:

  - x: along surge axis
  - y: along sway axis
  - z: height relative to MSL

PtfmRefPos (m, m, m)
  Platform reference point position relative to MSL. All platform motions and loads connect at this point.

PtfmRefOrient (deg, deg, deg)
  Platform reference orientation given as Euler angles [roll, pitch, yaw].

Turbine Initial Conditions
--------------------------

RotSpeed (rpm)
  Initial rotational speed of the rotor (in rotor coordinates).

NacYaw (deg)
  Initial or fixed nacelle yaw angle.

BldPitch (deg)
  Initial blade 1 pitch angle. If a multi-blade model is used, blade pitch control typically applies per blade in other modules; here this initializes blade 1.

Azimuth (deg)
  Initial rotor azimuth angle.

Wave Buoy
---------

WaveBuoyLoc (m, m)
  Location of the wave elevation measurement buoy in the tank coordinate frame. SeaState data is returned at each timestep at this location.

Output
------

SendScreenToFile (flag)
  If true, send screen output to a file named ``<OutRootName>.screen.log``.

OutFile (switch)
  Controls tabular output of channels:

  - 0: no output file of channels
  - 1: output file in text format (written at default DT)

OutFmt (string)
  Format specifier for text tabular output channels (excluding the time channel). Uses a Fortran-like format string.
  Example: ``ES20.6E2``.

VTK Visualization Output
------------------------

WrVTK_Dir (string)
  Output directory for VTK visualization files.

WrVTK (switch)
  VTK visualization data output:

  - 0: none
  - 1: initialization data only
  - 2: animation
  - 3: mode shapes

WrVTK_type (switch)
  Type of VTK visualization data:

  - 1: surfaces
  - 2: basic meshes (lines/points)
  - 3: all meshes (debug)

.. note::
   Only lines/points may be supported in some builds. If surfaces are not
   supported, use ``WrVTK_type = 2`` to visualize line/point data.

WrVTK_DT (s)
  Timestep for writing VTK files.

VTKNacDim (m, m, m, m, m, m)
  Nacelle dimensions for VTK surface rendering in the format ``[x0, y0, z0, Lx, Ly, Lz]``:

  - ``x0, y0, z0``: nacelle origin offsets
  - ``Lx, Ly, Lz``: nacelle extents along x, y, z

Implementation Notes and Best Practices
---------------------------------------

- Only floating MHK (``MHK = 2``) is currently supported; other MHK modes will
  not perform as expected.
- Ensure external file paths (*SeaState*, *MoorDyn*, *AeroDyn*, *InflowWind*)
  are valid relative to the working directory or specify absolute paths.
- Coordinate conventions:

  - Positions and heights are referenced to MSL unless otherwise noted.
  - The platform reference point (``PtfmRefPos``) is the coupling point for
    motions and loads.
  - The tower base is defined relative to ``PtfmRefPos`` in x and y, and to MSL
    in z.

- Choose ``OutFmt`` to balance precision and file size. The example ``ES20.6E2``
  is suitable for scientific notation with fixed width.
