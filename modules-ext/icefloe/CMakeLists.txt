
generate_f90_types(src/interfaces/FAST/IceFloe_FASTRegistry.inp
  IceFloe_Types.f90)


set(ICEFLOE_LIBS_SOURCES
  src/ranlux/RANLUX.f90

  src/icefloe/IceFlexBase.F90
  src/icefloe/IceFlexIEC.f90
  src/icefloe/IceFlexISO.f90
  src/icefloe/IceFloeBase.F90
  src/icefloe/coupledCrushing.F90
  src/icefloe/crushingIEC.F90
  src/icefloe/crushingISO.F90
  src/icefloe/iceInput.f90
  src/icefloe/iceLog.F90
  src/icefloe/intermittentCrushing.F90
  src/icefloe/lockInISO.F90
  src/icefloe/randomCrushing.F90

  src/interfaces/FAST/IceFloe.f90

  # Generated files
  IceFloe_Types.f90
  )

add_library(icefloelib ${ICEFLOE_LIBS_SOURCES})
target_link_libraries(icefloelib nwtclibs)

install(TARGETS icefloelib
  RUNTIME DESTINATION lib
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
