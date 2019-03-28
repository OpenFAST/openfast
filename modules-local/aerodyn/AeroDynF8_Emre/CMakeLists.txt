

add_subdirectory(AeroAcoustics)

generate_f90types(AeroAcoustics_Registry.txt AeroAcoustics_Types.f90)
generate_f90types(AirfoilInfo_Registry.txt AirfoilInfo_Types.f90 -noextrap)
generate_f90types(BEMT_Registry.txt BEMT_Types.f90)
generate_f90types(UnsteadyAero_Registry.txt UnsteadyAero_Types.f90)
generate_f90types(AeroDyn_Registry.txt AeroDyn_Types.f90)


add_library(aerodyn AeroAcoustics_Types.f90 AeroDyn_Types.f90 AeroDyn.f90 AeroDyn_IO.f90 AirfoilInfo.f90 AirfoilInfo_Types.f90 BEMT.f90 BEMTUncoupled.f90 BEMT_Registry.txt BEMT_Types.f90 UnsteadyAero.f90 UnsteadyAero_Types.f90 fmin_fcn.f90 mod_root1dim.f90)

target_link_libraries(aerodyn nwtclibs AeroAcoustics)

