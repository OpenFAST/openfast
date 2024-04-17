# OpenFAST docker images

## Summary
The `Dockerfile` in this directory can be used to reliably build OpenFAST as a docker image that can be run locally and
in the cloud without much setup. By default, it's based on Ubuntu Jammy and is optimised in size and performance for 
production use. A multi-stage build is used, producing an Ubuntu image with just `libblas-dev`, `liblapack-dev`, `nano`
and `openfast` added. The image built by this `Dockerfile` can be customised at build time using build arguments.

## Image registry
Production images of OpenFAST for the `linux/amd64` platform are available on the 
[NREL docker hub](https://hub.docker.com/r/nrel/openfast).

## Build arguments
Provide any of the following build arguments to customise the image at build time. 

| Name            | Type    | Allowed values                                                                                                                                                                                                                                             | Default                                                                                     | Description                                               |
| --------------- | ------- |------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------|-----------------------------------------------------------|
| `BASE`          | String  | Any valid docker image URI that has the `apt` package manager installed.                                                                                                                                                                                   | `ubuntu:jammy`                                                                              | The docker image to base the OpenFAST image on.           |
| `CMAKE_OPTIONS` | String  | Any number of valid space-separated `cmake` options in the same format they're normally passed to `cmake` directly. See the options relevant to OpenFAST [here.](https://openfast.readthedocs.io/en/main/source/install/index.html#openfast-cmake-options) | `-DBUILD_TESTING=OFF -DBUILD_FASTFARM=ON -DDOUBLE_PRECISION=OFF -DCMAKE_BUILD_TYPE=RELEASE` | Options to control how CMake is used to build OpenFAST.   |
| `BUILD_CORES`   | Integer | Any integer greater than 0.                                                                                                                                                                                                                                | `4`                                                                                         | The number of cores to use to build OpenFAST with `make`. |
| `TIMEZONE`      | String  | Any [valid timezone](https://en.wikipedia.org/wiki/List_of_tz_database_time_zones).                                                                                                                                                                        | `UTC`                                                                                       | The timezone to use when running OpenFAST.                |

For example, to build OpenFAST v3.5.3 for the `linux/amd64` platform and set `CMAKE_OPTIONS` so the testing tree is built:

```shell
# Run from the root of this repository.
git checkout v3.5.3
docker build -f share/docker/Dockerfile -t openfast:3.5.3 --platform=linux/amd64 --build-arg=CMAKE_OPTIONS='-DBUILD_TESTING=ON' .
```

**NOTE:** This version of the `Dockerfile` is only available in v3.5.3 and up of this repository. To build earlier 
versions of OpenFAST, check out the code at that version and recreate the `Dockerfile` from v3.5.3 (or above) in the 
checked-out repository first.

## Building development images
Development images can be built from the production image as a base. Simply start a new `Dockerfile` with:

```dockerfile
FROM nrel/openfast:3.5.3
```

Images can be built for different platforms using the `--platform` option when building the image.
