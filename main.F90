program main

  use amrex_base_module
  use read_amrex_subdomain_module

  implicit none

  character(1024) :: FileName
  integer(IntKi)  :: dims(3)
  real(ReKi)      :: origin(3)
  real(ReKi)      :: gridSpacing(3)
  real(c_double)  :: time

  integer :: idim, i, j, k
  integer(c_int) :: dlo(3), dhi(3)
  real(c_double), allocatable :: data(:,:,:,:)

  call amrex_init()

  FileName = "my_subvol00006"
  call ReadAMReXHeader(FileName, dims, origin, gridSpacing, time)

  print*, "dims = " , dims(1) , dims(2) , dims(3)
  print*, "origin = " , origin(1) , origin(2) , origin(3)
  print*, "gridSpacing = ", gridSpacing(1), gridSpacing(2), gridSpacing(3)
  print*, "time = ", time

  do idim = 1, 3
     dlo(idim) = dims(idim)/3
     dhi(idim) = dlo(idim) + dims(idim)/2
  end do

  allocate(data(3,dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3)))
  call ReadAMReXSubdomain(data, dlo, dhi)

  do k = dlo(3), dhi(3)
     do j = dlo(2), dhi(2)
        do i = dlo(1), dhi(1)
           print*, i,j,k, data(1,i,j,k), data(2,i,j,k), data(3,i,j,k)
        end do
     end do
  end do

  call amrex_finalize()

end program main
