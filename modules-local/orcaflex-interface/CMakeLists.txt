
generate_f90_types(src/OrcaFlexInterface.txt OrcaFlexInterface_Types.f90)

if (ORCA_DLL_LOAD)
  add_definitions(-DLibLoad)
else (ORCA_DLL_LOAD)
  add_definitions(-DNO_LibLoad)
endif (ORCA_DLL_LOAD)

add_library(orcaflexlib
  src/OrcaFlexInterface.f90

  # Generated files
  OrcaFlexInterface_Types.f90
  )
target_link_libraries(orcaflexlib nwtclibs)

add_executable(orca_driver
  src/OrcaDriver_Subs.f90
  src/OrcaDriver_Types.f90
  src/OrcaDriver.f90)
target_link_libraries(orca_driver orcaflexlib)

install(TARGETS orcaflexlib orca_driver
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
