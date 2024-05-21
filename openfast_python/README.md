# OpenFAST python readers/writers

> [!CAUTION]
> The `openfast_io` package on PyPI is currently called `octue-openfast` but will soon be renamed.

This package is a python wrapper comprising readers and writers for converting OpenFAST files to/from python objects. It
was originally written for [WEIS](https://github.com/WISDEM/WEIS/tree/77a878d7989b8c1d07d2244135ccd308a193a924/weis/aeroelasticse) and has been ported over to OpenFAST to make it more widely accessible. 

## Installation
Run either
```shell
pip install openfast_io
```
or
```shell
poetry add openfast_io
```

### Extra options
[ROSCO](https://github.com/NREL/ROSCO) can be installed as an optional dependency. Run either
```shell
pip install openfast_io[rosco]
```
or
```shell
poetry add -E rosco openfast_io
```
