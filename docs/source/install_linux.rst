OpenFAST
========

CMake Build Instructions
------------------------

::

    git clone https://github.com/OpenFAST/OpenFAST.git
    cd OpenFAST
    mkdir build && cd build
    cmake ../ 
    make 


Current CMake Options
~~~~~~~~~~~~~~~~~~~~~

-  ``DOUBLE_PRECISION`` - Enable/disable ``-DDOUBLE_PRECISION`` flag
   (Default: ON)
-  ``USE_DLL_INTERFACE`` - Enable dynamic library loading capability
   (Default: ON)
-  ``CMAKE_BUILD_TYPE`` - Release, Debug builds (Default: Release)
-  ``CMAKE_INSTALL_PREFIX`` - Set desired installation directory
-  ``BUILD_SHARED_LIBS`` - Enable/disable building shared libraries
   (Default: OFF)
-  ``FPE_TRAP_ENABLED`` - Enable Floating Point Exception trap
-  ``BUILD_CPP_API`` - Enable C++ API

Dependencies
~~~~~~~~~~~~

OpenFAST depends on the ``LAPACK`` libraries provided through the variable ``BLASLIB``. When building the ``C++``API, OpenFAST also depends on `HDF5 <https://support.hdfgroup.org/HDF5/>` (provided by ``HDF5_ROOT``) and `yaml-cpp <https://github.com/jbeder/yaml-cpp>` (provided by ``YAML_ROOT``). We recommend installing OpenFAST using `spack <https://spack.readthedocs.io/en/latest>`. However, we also provide some sample scripts in ``share`` folder if you choose to install without ``spack``.

