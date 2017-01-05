
generate_f90_types(src/OpenFOAM_Registry.txt OpenFOAM_Types.f90 -ccode)

add_library(openfoamtypeslib
  OpenFOAM_Types.f90)
target_link_libraries(openfoamtypeslib nwtclibs)

add_library(foamfastlib
  src/OpenFOAM.f90
  )
target_link_libraries(foamfastlib openfoamtypeslib openfast_prelib nwtclibs)

install(TARGETS openfoamtypeslib foamfastlib
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/OpenFOAM_Types.h
  DESTINATION include)
