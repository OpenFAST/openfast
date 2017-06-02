
add_subdirectory(fast)

if(BUILD_FAST_CPP_API)
   add_subdirectory(fast-cpp)
endif()

if(BUILD_FAST_FARM)
  add_subdirectory(fast-farm)
endif()

