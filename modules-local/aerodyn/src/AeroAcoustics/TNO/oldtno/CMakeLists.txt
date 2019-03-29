
add_library(TNO  int1.f90 int2.f90 TNOMods.f90 pressure.f90 qk61.f )

target_include_directories (TNO PUBLIC)
