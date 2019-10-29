openfast_dir=
yaml_install_dir=$openfast_dir/install/
hdf5_install_dir=$openfast_dir/install/

EXTRA_ARGS=$@

CC=mpicc CXX=mpic++ FC=gfortran cmake \
   -DCMAKE_INSTALL_PREFIX=$openfast_dir/install/ \
   -DCMAKE_BUILD_TYPE=RELEASE \
   -DBUILD_OPENFAST_CPP_API=ON \
   -DYAML_ROOT:PATH=$yaml_install_dir \
   -DHDF5_USE_STATIC_LIBRARIES=ON \
   -DHDF5_ROOT:PATH=$hdf5_install_dir \
   -DFPE_TRAP_ENABLED=OFF \
   $EXTRA_ARGS \
../ &> log.cmake

