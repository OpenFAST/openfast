
set(NWTCLIBS_SOURCES
  src/ModMesh.f90
  src/ModMesh_Mapping.f90
  src/ModMesh_Types.f90
  src/NWTC_Base.f90
  src/NWTC_IO.f90
  src/NWTC_Library.f90
  src/NWTC_Num.f90
  src/SingPrec.f90
  src/NWTC_Library_Types.f90

  # NetLib sources
  src/NetLib/Dierckx_FitPack/dierckx_fitpack.f
  src/NetLib/fftpack/fftpack4.1.f
  src/NetLib/scalapack/dlasrt2.f
  src/NetLib/scalapack/slasrt2.f
  src/NetLib/Dierckx_FitPack/NWTC_FitPack.f90
  src/NetLib/fftpack/NWTC_FFTPACK.f90
  src/NetLib/lapack/NWTC_LAPACK.f90
  src/NetLib/scalapack/NWTC_ScaLAPACK.f90
  )

get_filename_component(FCNAME ${CMAKE_Fortran_COMPILER} NAME)

if (FCNAME MATCHES "gfortran.*")
  if (APPLE OR UNIX)
    set(NWTC_SYS_FILE src/SysGnuLinux.f90)
  elseif (WIN32)
    set(NWTC_SYS_FILE src/SysGnuWin.f90)
  endif (APPLE OR UNIX)
elseif (FCNAME MATCHES "ifort.*")
  if (APPLE OR UNIX)
    set(NWTC_SYS_FILE src/SysIFL.f90)
  elseif (WIN32)
    set(NWTC_SYS_FILE src/SysIFV.f90)
  endif (APPLE OR UNIX)
endif (FCNAME MATCHES "gfortran.*")

if (NWTC_SYS_FILE)
  message("-- Setting system file as: ${NWTC_SYS_FILE}")
  list(APPEND NWTCLIBS_SOURCES ${SYS_FILE})
else (NWTC_SYS_FILE)
  message(FATAL_ERROR "Cannot determine system file used with NWTC_Library")
endif (NWTC_SYS_FILE)

add_library(nwtclibs ${NWTCLIBS_SOURCES})
target_link_libraries(nwtclibs ${LAPACK_LIBRARIES} ${CMAKE_DL_LIBS})

install(TARGETS nwtclibs
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
