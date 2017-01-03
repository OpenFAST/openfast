# OpenFAST

## CMake Build Instructions

```
git clone https://github.com/OpenFAST/OpenFAST.git
cd OpenFAST
mkdir build && cd build
cmake ../ 
make 
```

### Current CMake Options 

* `DOUBLE_PRECISION` - Enable/disable `-DDOUBLE_PRECISION` flag (Default: ON)
* `USE_DLL_INTERFACE` - Enable dynamic library loading capability (Default: ON)
* `CMAKE_BUILD_TYPE` - Release, Debug builds (Default: Release)
* `CMAKE_INSTALL_PREFIX` - Set desired installation directory
* `BUILD_SHARED_LIBS` - Enable/disable building shared libraries (Default: OFF)
