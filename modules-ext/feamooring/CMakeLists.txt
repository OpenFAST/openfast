
generate_f90_types(src/FEAM_Registry.txt FEAMooring_Types.f90)

set(FEAM_LIBS_SOURCES
  src/FEAM.f90

  # Generated files
  FEAMooring_Types.f90
  )

add_library(feamlib ${FEAM_LIBS_SOURCES})
target_link_libraries(feamlib nwtclibs)

add_executable(feam_driver
  src/FEAM_Driver.f90)
target_link_libraries(feam_driver feamlib)


install(TARGETS feamlib feam_driver
  RUNTIME DESTINATION lib
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
