module test_AMReX_reader
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use NWTC_Library
   use amrex_utils
   use read_amrex_subdomain_module
 
   implicit none

contains

   subroutine run_test_AMReX_reader(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("AMReX_read_subvol_0", AMReX_read_subvol_0), &
                  new_unittest("AMReX_read_subvol_1", AMReX_read_subvol_1) &
                  ]
   end subroutine

   subroutine AMReX_read_subvol_0(error)
      
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter :: DirPath = "data/subvolmultiple_fewerpoints_0_00006"
      real(SiKi), allocatable :: data(:,:,:,:)
      integer(IntKi)          :: ErrStat
      character(ErrMsgLen)    :: ErrMsg
      integer(IntKi)          :: i, j, k
      integer(IntKi)          :: dims(3)
      real(ReKi)              :: time, origin(3), gridSpacing(3), gridSize(3), bounds(2,3)

      ! Read header
      call amrex_read_header(trim(DirPath), time, dims, gridSpacing, origin, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_None); if (allocated(error)) return
      if (ErrStat /= ErrID_None) print*, "ErrMsg = ", ErrMsg

      ! Calculate grid size
      gridSize = gridSpacing*real(dims - 1, c_double)

      ! Calculate the grid bounds
      bounds(1,:) = origin
      bounds(2,:) = origin + gridSize

      ! Check time
      call check(error, time, 0.6_c_double); if (allocated(error)) return

      ! Check dimensions
      call check(error, dims(1), 3_c_int); if (allocated(error)) return
      call check(error, dims(2), 4_c_int); if (allocated(error)) return
      call check(error, dims(3), 5_c_int); if (allocated(error)) return

      ! Check spacing
      call check(error, gridSpacing(1), 1.0_c_double); if (allocated(error)) return
      call check(error, gridSpacing(2), 1.0_c_double); if (allocated(error)) return
      call check(error, gridSpacing(3), 1.0_c_double); if (allocated(error)) return

      ! Check grid size
      call check(error, bounds(1,1), 2.0_c_double); if (allocated(error)) return
      call check(error, bounds(1,2), 3.0_c_double); if (allocated(error)) return
      call check(error, bounds(1,3), 4.0_c_double); if (allocated(error)) return

      ! Check lower bounds
      call check(error, bounds(1,1), 6.5_c_double); if (allocated(error)) return
      call check(error, bounds(1,2), 6.5_c_double); if (allocated(error)) return
      call check(error, bounds(1,3), 6.5_c_double); if (allocated(error)) return

      ! Check upper bounds
      call check(error, bounds(2,1), 8.5_c_double); if (allocated(error)) return
      call check(error, bounds(2,2), 9.5_c_double); if (allocated(error)) return
      call check(error, bounds(2,3), 10.5_c_double); if (allocated(error)) return

      ! Display grid properties
      print*, "dir         = ", trim(DirPath)
      print*, "time        = ", time
      print*, "dims        = ", dims
      print*, "origin      = ", origin
      print*, "gridSpacing = ", gridSpacing
      print*, "gridSize    = ", gridSize
      print*, "lb          = ", bounds(1,:)
      print*, "ub          = ", bounds(2,:)
      print*, "nPoints     = ", product(dims)

      ! Allocate data array and fill with value
      call AllocAry(data, 3, dims(1), dims(2), dims(3), "data", ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_None); if (allocated(error)) return
      data = 9999.0_c_float

      ! Read grid data
      call amrex_read_data(trim(DirPath), data, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_None); if (allocated(error)) return

      ! Print the grid data
      do k = 1, dims(3)
         do j = 1, dims(2)
            do i = 1, dims(1)
               print*, i,j,k, data(1,i,j,k), data(2,i,j,k), data(3,i,j,k)
            end do
         end do
      end do

   end subroutine

   subroutine AMReX_read_subvol_1(error)
      
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter :: DirPath = "data/subvolmultiple_1_00016"
      real(SiKi), allocatable :: data(:,:,:,:)
      integer(IntKi)          :: ErrStat
      character(ErrMsgLen)    :: ErrMsg
      integer(IntKi)          :: i, j, k
      integer(IntKi)          :: dims(3)
      real(ReKi)              :: time, origin(3), gridSpacing(3)

      ! Read header
      call amrex_read_header(trim(DirPath), time, dims, gridSpacing, origin, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_None); if (allocated(error)) return

      ! Check dimensions
      call check(error, dims(1), 24_c_int); if (allocated(error)) return
      call check(error, dims(2), 24_c_int); if (allocated(error)) return
      call check(error, dims(3), 8_c_int); if (allocated(error)) return

      ! Check spacing
      call check(error, gridSpacing(1), 1.0_c_double); if (allocated(error)) return
      call check(error, gridSpacing(2), 1.0_c_double); if (allocated(error)) return
      call check(error, gridSpacing(3), 1.0_c_double); if (allocated(error)) return

      ! Check origin
      call check(error, origin(1), 0.0_c_double); if (allocated(error)) return
      call check(error, origin(2), 0.0_c_double); if (allocated(error)) return
      call check(error, origin(3), 0.0_c_double); if (allocated(error)) return

      print*, "ErrStat = ", ErrStat
      if (ErrStat /= 0) then
         print*, "ErrMsg: ", ErrMsg
      end if

      print*, "dir = ", trim(DirPath)
      print*, "dims = " , dims(1) , dims(2) , dims(3)
      print*, "origin = " , origin(1) , origin(2) , origin(3)
      print*, "gridSpacing = ", gridSpacing(1), gridSpacing(2), gridSpacing(3)
      print*, "time = ", time

      ! Allocate data array and fill with value
      call AllocAry(data, 3, dims(1), dims(2), dims(3), "Data", ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_None); if (allocated(error)) return
      data = 9999.0_c_float

      ! Read data
      call amrex_read_data(trim(DirPath), data, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_None); if (allocated(error)) return

      ! do k = 1, dims(3)
      !    do j = 1, dims(2)
      !       do i = 1, dims(1)
      !          print*, i,j,k, data(1,i,j,k), data(2,i,j,k), data(3,i,j,k)
      !       end do
      !    end do
      ! end do

   end subroutine


end module
