**Feature or improvement description**
Pull request to merge `dev` into `main` for release version X.Y.Z

See the milestone and project pages for additional information

    https://github.com/OpenFAST/openfast/milestone/XYZ

Test results, if applicable
See GitHub Actions

### Release checklist:
- [ ] Update the documentation version in docs/conf.py
- [ ] Update the versions in docs/source/user/api\_change.rst
- [ ] Update version info in openfast\_io/pyproject.toml (`openfast_io` package)
- [ ] Update version info in glue-codes/python/pyproject.toml (`pyOpenFAST` package for testing)
- [ ] Verify readthedocs builds correctly
- [ ] Create an annotated tag in OpenFAST during merge (mark as most recent if necessary)
- [ ] Create a merge commit in r-test and add a corresponding annotated tag
- [ ] Upload Docker image
- [ ] Compile executables for Windows builds
    - [ ] `AeroDisk_Driver_x64.exe`
    - [ ] `AeroDyn_Driver_x64.exe`
    - [ ] `AeroDyn_Driver_x64_OpenMP.exe`
    - [ ] `AeroDyn_Inflow_c_binding_x64.dll`
    - [ ] `AeroDyn_Inflow_c_binding_x64_OpenMP.dll`
    - [ ] `BeamDyn_Driver_x64.exe`
    - [ ] `DISCON.dll (x64)`
    - [ ] `DISCON_ITIBarge.dll (x64)`
    - [ ] `DISCON_OC3Hywind.dll (x64)`
    - [ ] `FAST.Farm_x64.exe`
    - [ ] `FAST.Farm_x64_OMP.exe`
    - [ ] `FAST_SFunc.mexw64`
    - [ ] `HydroDynDriver_x64.exe`
    - [ ] `HydroDyn_C_Binding_x64.dll`
    - [ ] `IinflowWind_c_binding_x64.dll`
    - [ ] `InflowWind_Driver_x64.exe`
    - [ ] `InflowWind_Driver_x64_OpenMP.exe`
    - [ ] `MoorDyn_Driver_x64.exe`
    - [ ] `MoorDyn_c_binding_x64.dll`
    - [ ] `OpenFAST-Simulink_x64.dll`
    - [ ] `openfast_x64.exe`
    - [ ] `SeaStateDriver_x64.exe`
    - [ ] `SeaState_c_binding_x64.dll`
    - [ ] `SimpleElastoDyn_x64.exe`
    - [ ] `SubDyn_x64.exe`
    - [ ] `Turbsim_x64.exe`
    - [ ] `UnsteadyAero_x64.exe`



# Release Overview
------



### Contribution Acknowledgements

### Statistics (since X.Y.Z)



# Changelog (from X.Y.Z)
------

## General

### Build systems

#### CMake


#### Visual Studio (Windows)



### Docker


### Documentation




## Solvers

### FAST.Farm


### OpenFAST


### OpenFAST interfaces

#### OpenFASTcpp


#### Simulink



## Modules

### Multiple


### AeroDisk


### AeroDyn

#### Unsteady



#### OLAF 


#### AeroDyn Driver / AeroDyn\_Inflow\_C\_Bindings interface



### ElastoDyn


### ExtInflow


### ExtLoads


### HydroDyn


### InflowWind


### MAP++


### MoorDyn


### NWTC-Library


### Registry


### SeaState


### Simplified-ElastoDyn


### SubDyn



## Testing and input file processing

### openfast_io


### GitHub actions


### Regression and Unit testing




## Input file changes


## Known issues



# Precompiled Windows Binaries
The binary files in this release were built with the Visual Studio solution files distributed with OpenFAST (not using cmake), using

- Intel Fortran Essentials 2025.3.0.333
- Microsoft Visual Studio 2022 Version 17.14.23.
- MATLAB 2025.2.999 (R2025b)
- Executables with `_OpenMP` or `_OMP` in the name are built with OpenMP libraries and linked with dynamic libraries.
   - You will need [this Intel Fortran redistributable package](https://registrationcenter-download.intel.com/akdlm/IRC_NAS/0dc56e76-d2c0-4bb8-9c83-c2ee3952b855/w_ifx_runtime_p_2025.2.1.1001.exe) installed to use these executables if you do not already have Intel Fortran OneAPI 2024 installed. See the installation instructions [here](https://software.intel.com/content/www/us/en/develop/articles/redistributable-libraries-for-intel-c-and-fortran-2022-compilers-for-windows.html).

**The other OpenFAST executables DO NOT require these redistributable libraries to be installed. Instead, they were built with static libraries.**
