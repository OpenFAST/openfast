
add_subdirectory(openfast)

if(BUILD_OPENFAST_CPP_API)
   add_subdirectory(openfast-cpp)
endif()

if(BUILD_FAST_FARM)
  add_subdirectory(fast-farm)
endif()

