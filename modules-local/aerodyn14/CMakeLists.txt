
generate_f90_types(src/Registry-AD14.txt AeroDyn14_Types.f90)
generate_f90_types(src/Registry-DWM.txt DWM_Types.f90)

set(AD14_LIBS_SOURCES
  src/AeroDyn14.f90
  src/AeroSubs.f90
  src/DWM.f90
  src/DWM_Wake_Sub_ver2.f90
  src/GenSubs.f90

  # Generated files
  AeroDyn14_Types.f90
  DWM_Types.f90
  )

add_library(aerodyn14lib ${AD14_LIBS_SOURCES})
target_link_libraries(aerodyn14lib nwtclibs)

# add_executable(aerodyn14
#   src/AD_Driver.f90)
# target_link_libraries(aerodyn14 aerodyn14lib)
# 
# 
# add_executable(dwm_driver_wind_farm
#   src/DWM_driver_wind_farm_mod.f90
#   src/DWM_driver_wind_farm_sub.f90
#   src/DWM_driver_wind_farm.f90
#   )
# target_link_libraries(dwm_driver_wind_farm aerodyn14lib)

install(TARGETS aerodyn14lib
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
