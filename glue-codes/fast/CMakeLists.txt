
add_executable(openfast
  src/FAST_Prog.f90)
target_link_libraries(openfast openfast_postlib foamfastlib)

install(TARGETS openfast
  RUNTIME DESTINATION bin)
