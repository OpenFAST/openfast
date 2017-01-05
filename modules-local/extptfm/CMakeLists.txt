
generate_f90_types(src/ExtPtfm_MCKF_Registry.txt ExtPtfm_MCKF_Types.f90)

add_library(extptfm_mckflib
  src/ExtPtfm_MCKF.f90

  # Generated files
  ExtPtfm_MCKF_Types.f90
  )

target_link_libraries(extptfm_mckflib nwtclibs)

install(TARGETS extptfm_mckflib
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
