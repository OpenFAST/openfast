module test_AMReX_reader
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use NWTC_Num

   use amrex_base_module
   use read_amrex_subdomain_module
 
   implicit none
 
   character(1024) :: FileName, ErrorMsg
   integer(c_int)  :: ErrorStat
   integer(c_int)  :: dims(3)
   real(c_double)  :: origin(3)
   real(c_double)  :: gridSpacing(3)
   real(c_double)  :: time
 
   integer :: idim, i, j, k
   integer(c_int) :: dlo(3), dhi(3)
   real(c_double), allocatable :: data(:,:,:,:)

contains

   subroutine run_test_AMReX_reader(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("AMReX_read_subvol", AMReX_read_subvol) &
                  ]
   end subroutine



   subroutine AMReX_read_subvol(error)
      type(error_type), allocatable, intent(out) :: error

      integer(IntKi)       :: ErrStat
      character(ErrMsgLen) :: ErrMsg
      character(1024)      :: testname

      ! initialize NWTC_Num constants
      call SetConstants()

      testname = "my_subvol00006"


      call amrex_init()
    
      FileName = "my_subvol00006"
      call ReadAMReXHeader(FileName, dims, origin, gridSpacing, time, ErrorStat, ErrorMsg)
    
      print*, "ErrorStat = ", ErrorStat
      if (ErrorStat .ne. 0) then
         print*, "ErrorMsg: ", ErrorMsg
      end if
      print*, "dims = " , dims(1) , dims(2) , dims(3)
      print*, "origin = " , origin(1) , origin(2) , origin(3)
      print*, "gridSpacing = ", gridSpacing(1), gridSpacing(2), gridSpacing(3)
      print*, "time = ", time

      ! Fail here: no file actually exists, so can't check anything yet
      ErrStat=ErrID_Fatal
      call check(error, ErrID_None, ErrStat); if (allocated(error)) return



      do idim = 1, 3
         dlo(idim) = dims(idim)/2
         dhi(idim) = dlo(idim) + dims(idim)/4
      end do
    
      allocate(data(3,dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3)))
      call ReadAMReXSubdomain(data, dlo, dhi, ErrorStat, ErrorMsg)
    
      print*, "ErrorStat = ", ErrorStat
      if (ErrorStat .ne. 0) then
         print*, "ErrorMsg: ", ErrorMsg
      end if
    
      do k = dlo(3), dhi(3)
         do j = dlo(2), dhi(2)
            do i = dlo(1), dhi(1)
               print*, i,j,k, data(1,i,j,k), data(2,i,j,k), data(3,i,j,k)
            end do
         end do
      end do
    
      call amrex_finalize()
   end subroutine

end module
