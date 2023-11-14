# Visual Studio builds for Windows
The following solution files are available for code development on Windows using the Intel Fortran compiler with Visual Studio.

- [OpenFAST](FAST/FAST.sln)
  This contains builds for both the command-line OpenFAST executable as well as the DLL for use with the OpenFAST-Simulink interface.
- [FAST.Farm](FAST-farm/FAST-Farm.sln)
  This contains the build configurations for FAST.Farm.
- Module-level drivers:
   - AeroDynamics and HydroDynamics:
     - [AeroDyn driver](AeroDyn/AeroDyn_Driver.sln)
     - [UnsteadyAero driver](UnsteadyAero/UnsteadyAero.sln)
     - [HydroDyn driver](HydroDyn/HydroDynDriver.sln)
   - Structural: 
     - [BeamDyn driver](BeamDyn/BeamDyn-w-registry.sln)
     - [SubDyn driver](SubDyn/SubDyn.sln)
   - Wind/Wave conditions:
      - [TurbSim](TurbSim/TurbSim.sln) Generates wind files
      - [InflowWind driver](InflowWind/InflowWind_driver.sln) Reads and interpolates existing wind files
      - [InflowWind c binding](InflowWind/InflowWind_c_binding.sln) Creates a library (DLL/so) of the InflowWind routines that can be called from a C interface. **TO DO: Combine this with InflowWind driver for easier maintenance**
      - [SeaState driver](SeaState/SeaStateDriver.sln) Waves and currents
- Other:
  - [Discon](Discon/Discon.sln) 
    This solution file contains all 3 controllers used in the OpenFAST r-test (with the NREL 5MW model).
    It also contains the controller used with the FAST.Farm super-controller.
  - [SC_DLL](SC_DLL.sln) This solution file builds a template supercontroller to be used with FAST.Farm.
  - [OpenFAST Registry](Registry/Registry.sln)
    The Registry project is included in almost every other solution file, so this solution file is only for debugging changes to the OpenFAST Registry.
