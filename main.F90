program main

  use amrex_base_module
  use read_amrex_subdomain_module

  implicit none

  character(1024) :: FileName
  integer(IntKi)  :: dims(3)
  real(ReKi)      :: origin(3)
  real(ReKi)      :: gridSpacing(3)
  real(c_double)  :: time

  call amrex_init()

  FileName = "my_subvol00006"
  call ReadAMReXHeader(FileName, dims, origin, gridSpacing, time)

  print*, "dims = " , dims(1) , dims(2) , dims(3)
  print*, "origin = " , origin(1) , origin(2) , origin(3)
  print*, "gridSpacing = ", gridSpacing(1), gridSpacing(2), gridSpacing(3)
  print*, "time = ", time

  call amrex_finalize()

end program main
