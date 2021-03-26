# Visual Studio builds for Windows
The following solution files are available for code development on Windows using the Intel Fortran compiler with Visual Studio.

- [OpenFAST](FAST/FAST.sln)
  This contains builds for both the command-line OpenFAST executable as well as the DLL for use with the OpenFAST-Simulink interface.
- Module-level drivers:
   - AeroDynamics:
     - [AeroDyn driver](AeroDyn/AeroDyn_Driver.sln)
     - [UnsteadyAero driver](UnsteadyAero/UnsteadyAero.sln)
   - Structural: 
     - [BeamDyn driver](BeamDyn/BeamDyn-w-registry.sln)
     - [SubDyn driver](SubDyn/SubDyn.sln)
   - Wind/Wave conditions
      - [TurbSim](TurbSim/TurbSim.sln) Generates wind files
      - [InflowWind driver](InflowWind/InflowWind_driver.sln) Reads and interpolates existing wind files
      - [HydroDyn driver](HydroDyn/HydroDynDriver.sln)
- Other:
  - [Discon](Discon/Discon.sln) This solution file contains all 3 controllers used in the OpenFAST r-test (with the NREL 5MW model).
  - [OpenFAST Registry](Registry/Registry.sln)
    The Registry project is included in almost every other solution file, so this solution file is only for debugging changes to the OpenFAST Registry.
