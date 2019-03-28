
add_subdirectory(src/TNO)
add_subdirectory(src/TINoise)
###add_subdirectory(src/Xfoil)


generate_f90_types(src/AeroAcoustics_Registry.txt AeroAcoustics_Types.f90)

add_library(AeroAcoustics  src/AeroAcoustics.f90 src/AeroAcoustics_IO.f90 AirfoilInfo_Types.f90   AeroAcoustics_Types.f90)

target_link_libraries(AeroAcoustics nwtclibs FullGuidati TNO )
###target_link_libraries(AeroAcoustics nwtclibs  RunXfoil  FullGuidati TNO )

