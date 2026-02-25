module test_AMReX_reader
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use NWTC_Library
   use amrex_utils
 
   implicit none

contains

   subroutine run_test_AMReX_reader(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("AMReX_test_read_subvol_0", AMReX_test_read_subvol_0), &
                  new_unittest("AMReX_test_read_subvol_1", AMReX_test_read_subvol_1), &
                  new_unittest("AMReX_test_amrex_find_subvols_1", AMReX_test_amrex_find_subvols_1), &
                  new_unittest("AMReX_test_amrex_find_subvols_2", AMReX_test_amrex_find_subvols_2), &
                  new_unittest("AMReX_test_amrex_find_subvols_3", AMReX_test_amrex_find_subvols_3) &
                  ]
   end subroutine

   subroutine AMReX_test_read_subvol_0(error)
      
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
      call check(error, dims(1), 3_c_int, more="dims(1)"); if (allocated(error)) return
      call check(error, dims(2), 4_c_int, more="dims(2)"); if (allocated(error)) return
      call check(error, dims(3), 5_c_int, more="dims(3)"); if (allocated(error)) return

      ! Check spacing
      call check(error, gridSpacing(1), 1.0_c_double, more="gridSpacing(1)"); if (allocated(error)) return
      call check(error, gridSpacing(2), 1.0_c_double, more="gridSpacing(2)"); if (allocated(error)) return
      call check(error, gridSpacing(3), 1.0_c_double, more="gridSpacing(3)"); if (allocated(error)) return

      ! Check grid size
      call check(error, gridSize(1), 2.0_c_double, more="gridSize(1)"); if (allocated(error)) return
      call check(error, gridSize(2), 3.0_c_double, more="gridSize(2)"); if (allocated(error)) return
      call check(error, gridSize(3), 4.0_c_double, more="gridSize(3)"); if (allocated(error)) return

      ! Check lower bounds
      call check(error, bounds(1,1), 6.5_c_double, more="origin(1)"); if (allocated(error)) return
      call check(error, bounds(1,2), 6.5_c_double, more="origin(2)"); if (allocated(error)) return
      call check(error, bounds(1,3), 6.5_c_double, more="origin(3)"); if (allocated(error)) return

      ! Check upper bounds
      call check(error, bounds(2,1), 8.5_c_double, more="ub(1)"); if (allocated(error)) return
      call check(error, bounds(2,2), 9.5_c_double, more="ub(2)"); if (allocated(error)) return
      call check(error, bounds(2,3), 10.5_c_double, more="ub(3)"); if (allocated(error)) return

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

   subroutine AMReX_test_read_subvol_1(error)
      
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter :: DirPath = "data/subvolmultiple_1_00016"
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
      call check(error, time, 1.6_c_double); if (allocated(error)) return

      ! Check dimensions
      call check(error, dims(1), 2_c_int, more="dims(1)"); if (allocated(error)) return
      call check(error, dims(2), 3_c_int, more="dims(2)"); if (allocated(error)) return
      call check(error, dims(3), 5_c_int, more="dims(3)"); if (allocated(error)) return

      ! Check spacing
      call check(error, gridSpacing(1), 1.0_c_double, more="gridSpacing(1)"); if (allocated(error)) return
      call check(error, gridSpacing(2), 1.0_c_double, more="gridSpacing(2)"); if (allocated(error)) return
      call check(error, gridSpacing(3), 1.0_c_double, more="gridSpacing(3)"); if (allocated(error)) return

      ! Check grid size
      call check(error, gridSize(1), 1.0_c_double, more="gridSize(1)"); if (allocated(error)) return
      call check(error, gridSize(2), 2.0_c_double, more="gridSize(2)"); if (allocated(error)) return
      call check(error, gridSize(3), 4.0_c_double, more="gridSize(3)"); if (allocated(error)) return

      ! Check lower bounds
      call check(error, bounds(1,1), 6.5_c_double, more="origin(1)"); if (allocated(error)) return
      call check(error, bounds(1,2), 6.5_c_double, more="origin(2)"); if (allocated(error)) return
      call check(error, bounds(1,3), 6.5_c_double, more="origin(3)"); if (allocated(error)) return

      ! Check upper bounds
      call check(error, bounds(2,1), 7.5_c_double, more="ub(1)"); if (allocated(error)) return
      call check(error, bounds(2,2), 8.5_c_double, more="ub(2)"); if (allocated(error)) return
      call check(error, bounds(2,3), 10.5_c_double, more="ub(3)"); if (allocated(error)) return

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

   ! Test finding sub-volumes with start index of 0 and actual DT
   subroutine AMReX_test_amrex_find_subvols_1(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter    :: DirPath = "data/subvolmultiple"
      integer(IntKi), parameter  :: SubVol = 0
      real(DbKi), parameter      :: DT = 0.6_DbKi
      integer(IntKi), parameter  :: NumSteps = 5
      character(*), parameter    :: StartIndex = "00000"

      integer(IntKi)             :: FirstIndex, IndexDelta
      integer(IntKi)             :: ErrStat
      character(ErrMsgLen)       :: ErrMsg

      call amrex_find_subvols(DirPath, SubVol, DT, NumSteps, StartIndex, &
                              FirstIndex, IndexDelta, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_None, more="amrex_find_subvols: "//trim(ErrMsg)); if (allocated(error)) return

      call check(error, FirstIndex, 0); if (allocated(error)) return
      call check(error, IndexDelta, 6); if (allocated(error)) return

   end subroutine

   ! Test finding sub-volumes with nonzero start index and actual DT
   ! Sub-volume 1 has DT = 0.8 so this finds every directory starting at index 16
   subroutine AMReX_test_amrex_find_subvols_2(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter    :: DirPath = "data/subvolmultiple"
      integer(IntKi), parameter  :: SubVol = 1
      real(DbKi), parameter      :: DT = 0.8_DbKi
      integer(IntKi), parameter  :: NumSteps = 3
      character(*), parameter    :: StartIndex = "00016"

      integer(IntKi)             :: FirstIndex, IndexDelta
      integer(IntKi)             :: ErrStat
      character(ErrMsgLen)       :: ErrMsg

      call amrex_find_subvols(DirPath, SubVol, DT, NumSteps, StartIndex, &
                              FirstIndex, IndexDelta, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_None, more="amrex_find_subvols: "//trim(ErrMsg)); if (allocated(error)) return

      call check(error, FirstIndex, 16); if (allocated(error)) return
      call check(error, IndexDelta, 8); if (allocated(error)) return

   end subroutine

   ! Test finding sub-volumes with nonzero start index and larger DT
   ! Sub-volume 0 actual DT = 0.6 so this test skips every other directory
   ! starting at index 6
   subroutine AMReX_test_amrex_find_subvols_3(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter    :: DirPath = "data/subvolmultiple"
      integer(IntKi), parameter  :: SubVol = 0
      real(DbKi), parameter      :: DT = 1.2_DbKi
      integer(IntKi), parameter  :: NumSteps = 3
      character(*), parameter    :: StartIndex = "00006"

      integer(IntKi)             :: FirstIndex, IndexDelta
      integer(IntKi)             :: ErrStat
      character(ErrMsgLen)       :: ErrMsg

      call amrex_find_subvols(DirPath, SubVol, DT, NumSteps, StartIndex, &
                              FirstIndex, IndexDelta, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_None, more="amrex_find_subvols: "//trim(ErrMsg)); if (allocated(error)) return

      call check(error, FirstIndex, 6, more="FirstIndex"); if (allocated(error)) return
      call check(error, IndexDelta, 12, more="IndexDelta"); if (allocated(error)) return

   end subroutine

end module
