
add_subdirectory(TNO)
add_subdirectory(TINoise)
###add_subdirectory(src/Xfoil)

# NOTE: generated as well in parent dir
generate_f90_types(../AirfoilInfo_Registry.txt AirfoilInfo_Types.f90 -noextrap)
generate_f90_types(../AeroAcoustics_Registry.txt AeroAcoustics_Types.f90)

add_library(AeroAcoustics  AeroAcoustics.f90 AeroAcoustics_IO.f90 AirfoilInfo_Types.f90   AeroAcoustics_Types.f90)

target_link_libraries(AeroAcoustics nwtclibs FullGuidati TNO )
###target_link_libraries(AeroAcoustics nwtclibs  RunXfoil  FullGuidati TNO )

