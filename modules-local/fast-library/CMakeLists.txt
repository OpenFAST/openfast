
if (APPLE OR UNIX)
  if (FCNAME MATCHES "gfortran.*")
    add_definitions(-DIMPLICIT_DLLEXPORT)
  endif (FCNAME MATCHES "gfortran.*")
endif (APPLE OR UNIX)

generate_f90_types(src/FAST_Registry.txt FAST_Types.f90 -noextrap)

add_library(openfast_prelib
  FAST_Types.f90)
target_link_libraries(openfast_prelib
  nwtclibs
  ifwlib
  aerodynlib
  aerodyn14lib
  servodynlib
  elastodynlib
  beamdynlib
  subdynlib
  hydrodynlib
  orcaflexlib
  extptfm_mckflib
  openfoamtypeslib
  feamlib
  moordynlib
  icedynlib
  icefloelib
  maplib
  )

add_library(openfast_postlib
  src/FAST_Lin.f90
  src/FAST_Mods.f90
  src/FAST_Subs.f90
  src/FAST_Solver.f90
  )
target_link_libraries(openfast_postlib openfast_prelib foamfastlib)

add_library(openfastlib
  src/FAST_Library.f90)
target_link_libraries(openfastlib
  openfast_postlib openfast_prelib foamfastlib)

install(TARGETS openfastlib openfast_prelib openfast_postlib
  RUNTIME DESTINATION lib
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

install(FILES
  src/FAST_Library.h
  DESTINATION include)
