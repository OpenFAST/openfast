
generate_f90_types(src/MoorDyn_Registry.txt MoorDyn_Types.f90)

set(MOORDYN_LIBS_SOURCES
  src/MoorDyn.f90
  src/MoorDyn_IO.f90

  # Generated files
  MoorDyn_Types.f90)

add_library(moordynlib ${MOORDYN_LIBS_SOURCES})
target_link_libraries(moordynlib nwtclibs)

install(TARGETS moordynlib
  RUNTIME DESTINATION lib
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
