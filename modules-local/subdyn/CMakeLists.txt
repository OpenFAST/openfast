
generate_f90_types(src/SubDyn_Registry.txt SubDyn_Types.f90)

set(SUBDYN_SOURCES
  src/SubDyn.f90
  src/SD_FEM.f90
  src/SubDyn_Output.f90
  src/qsort_c_module.f90

  SubDyn_Types.f90
  )

add_library(subdynlib ${SUBDYN_SOURCES})
target_link_libraries(subdynlib nwtclibs)

set(SUBDYN_DRIVER_SOURCES
  src/SubDyn_Driver.f90)

add_executable(subdyn ${SUBDYN_DRIVER_SOURCES})
target_link_libraries(subdyn subdynlib nwtclibs)

install(TARGETS subdynlib subdyn
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
