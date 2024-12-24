# OpenFAST python readers/writers

This package is a python wrapper comprising readers and writers for converting OpenFAST files to/from python objects. It
was originally written for [WEIS](https://github.com/WISDEM/WEIS/tree/77a878d7989b8c1d07d2244135ccd308a193a924/weis/aeroelasticse) and has been ported over to OpenFAST to make it more widely accessible. 

## Installation
Installation with [Anaconda](https://www.anaconda.com) is the recommended approach because of the ability to create self-contained environments suitable for testing and analysis.

### Installation as a "library"

To use `openfast_io` as a library for incorporation into other scripts or tools, it is available via (assuming that you have already setup your python environment):

```shell
pip install openfast_io
```

### Installation as an editable library

These instructions are for interaction directly with the `openfast_io` source code.

0. Follow this step only if you have not cloned the OpenFAST repo.
    ```shell
    git clone https://github.com/OpenFAST/OpenFAST.git
    cd OpenFAST
    ```

1. Assuming you are within the OpenFAST directory.
    ```shell
    cd openfast_io
    pip install -e .
    ```

2. To test `openfast_io`, OpenFAST must be compiled within the build folder, then run:

    ```shell
    cd tests
    pytest test_of_io_pytest.py
    ```

### Extra options
[ROSCO](https://github.com/NREL/ROSCO) can be installed as an optional dependency. Run either
```shell
pip install openfast_io[rosco]
```

## Development and testing
To contribute to the development of `openfast_io`, install additioal depemndancies using:

```shell
pip install -e ".[all]"
```
