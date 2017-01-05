
generate_f90_types(src/Registry_IceDyn.txt IceDyn_Types.f90)

add_library(icedynlib
  src/IceDyn.f90

  # Generated files
  IceDyn_Types.f90)

target_link_libraries(icedynlib nwtclibs)

install(TARGETS icedynlib
  RUNTIME DESTINATION lib
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
