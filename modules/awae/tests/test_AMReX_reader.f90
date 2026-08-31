module test_AMReX_reader
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use NWTC_Library
   use amrex_utils
   use AWAE_Types, only: AWAE_ParameterType
   use AWAE_IO, only: ReadWindAMReX
 
   implicit none

contains

   subroutine run_test_AMReX_reader(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("AMReX_test_read_subvol_0", AMReX_test_read_subvol_0), &
                  new_unittest("AMReX_test_read_subvol_1", AMReX_test_read_subvol_1), &
                  new_unittest("AMReX_test_amrex_find_subvols_1", AMReX_test_amrex_find_subvols_1), &
                  new_unittest("AMReX_test_amrex_find_subvols_2", AMReX_test_amrex_find_subvols_2), &
                  new_unittest("AMReX_test_amrex_find_subvols_3", AMReX_test_amrex_find_subvols_3), &
                  new_unittest("AMReX_test_find_subvols_variable_stride", AMReX_test_find_subvols_variable_stride), &
                  new_unittest("AMReX_test_find_subvols_missing_step", AMReX_test_find_subvols_missing_step), &
                  new_unittest("AMReX_test_find_subvols_duplicate_step", AMReX_test_find_subvols_duplicate_step), &
                  new_unittest("AMReX_test_find_subvols_out_of_tolerance", AMReX_test_find_subvols_out_of_tolerance), &
                  new_unittest("AMReX_test_find_subvols_wide_index", AMReX_test_find_subvols_wide_index), &
                  new_unittest("AMReX_test_ReadWindAMReX_lookup", AMReX_test_ReadWindAMReX_lookup), &
                  new_unittest("AMReX_test_ReadWindAMReX_out_of_range", AMReX_test_ReadWindAMReX_out_of_range), &
                  new_unittest("AMReX_test_ReadWindAMReX_wide_index", AMReX_test_ReadWindAMReX_wide_index), &
                  new_unittest("AMReX_test_ReadWindAMReX_no_table", AMReX_test_ReadWindAMReX_no_table), &
                  new_unittest("AMReX_test_find_subvols_stale_beyond_window", AMReX_test_find_subvols_stale_beyond_window), &
                  new_unittest("AMReX_test_find_subvols_real_tiled", AMReX_test_find_subvols_real_tiled), &
                  new_unittest("AMReX_test_read_real_tiled", AMReX_test_read_real_tiled), &
                  new_unittest("AMReX_test_header_text_agrees_on_tiled", AMReX_test_header_text_agrees_on_tiled), &
                  new_unittest("AMReX_test_header_text_disagrees_on_subset", AMReX_test_header_text_disagrees_on_subset), &
                  new_unittest("AMReX_test_find_subvols_prefix_form", AMReX_test_find_subvols_prefix_form), &
                  new_unittest("AMReX_test_read_data_wrong_dims", AMReX_test_read_data_wrong_dims) &
                  ]
   end subroutine

   subroutine AMReX_test_read_subvol_0(error)
      
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter :: DirPath = "data/subvolmultiple_0_00006"
      real(SiKi), allocatable :: data(:,:,:,:)
      integer(IntKi)          :: ErrStat
      character(ErrMsgLen)    :: ErrMsg
      integer(IntKi)          :: i, j, k
      integer(IntKi)          :: dims(3)
      real(DbKi)              :: time
      real(ReKi)              :: origin(3), gridSpacing(3), gridSize(3), bounds(2,3)

      ! Read header
      call amrex_read_header(trim(DirPath), time, dims, gridSpacing, origin, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_None); if (allocated(error)) return
      if (ErrStat /= ErrID_None) print*, "ErrMsg = ", ErrMsg

      ! Calculate grid size
      gridSize = gridSpacing*real(dims - 1, ReKi)

      ! Calculate the grid bounds
      bounds(1,:) = origin
      bounds(2,:) = origin + gridSize

      ! Check time
      call check(error, time, 0.6_DbKi); if (allocated(error)) return

      ! Check dimensions
      call check(error, dims(1), 3_c_int, more="dims(1)"); if (allocated(error)) return
      call check(error, dims(2), 4_c_int, more="dims(2)"); if (allocated(error)) return
      call check(error, dims(3), 5_c_int, more="dims(3)"); if (allocated(error)) return

      ! Check spacing
      call check(error, gridSpacing(1), 1.0_ReKi, more="gridSpacing(1)"); if (allocated(error)) return
      call check(error, gridSpacing(2), 1.0_ReKi, more="gridSpacing(2)"); if (allocated(error)) return
      call check(error, gridSpacing(3), 1.0_ReKi, more="gridSpacing(3)"); if (allocated(error)) return

      ! Check grid size
      call check(error, gridSize(1), 2.0_ReKi, more="gridSize(1)"); if (allocated(error)) return
      call check(error, gridSize(2), 3.0_ReKi, more="gridSize(2)"); if (allocated(error)) return
      call check(error, gridSize(3), 4.0_ReKi, more="gridSize(3)"); if (allocated(error)) return

      ! Check lower bounds
      call check(error, bounds(1,1), 6.5_ReKi, more="origin(1)"); if (allocated(error)) return
      call check(error, bounds(1,2), 6.5_ReKi, more="origin(2)"); if (allocated(error)) return
      call check(error, bounds(1,3), 6.5_ReKi, more="origin(3)"); if (allocated(error)) return

      ! Check upper bounds
      call check(error, bounds(2,1), 8.5_ReKi, more="ub(1)"); if (allocated(error)) return
      call check(error, bounds(2,2), 9.5_ReKi, more="ub(2)"); if (allocated(error)) return
      call check(error, bounds(2,3), 10.5_ReKi, more="ub(3)"); if (allocated(error)) return

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
      real(DbKi)              :: time
      real(ReKi)              :: origin(3), gridSpacing(3), gridSize(3), bounds(2,3)

      ! Read header
      call amrex_read_header(trim(DirPath), time, dims, gridSpacing, origin, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_None); if (allocated(error)) return
      if (ErrStat /= ErrID_None) print*, "ErrMsg = ", ErrMsg

      ! Calculate grid size
      gridSize = gridSpacing*real(dims - 1, ReKi)

      ! Calculate the grid bounds
      bounds(1,:) = origin
      bounds(2,:) = origin + gridSize

      ! Check time
      call check(error, time, 1.6_DbKi); if (allocated(error)) return

      ! Check dimensions
      call check(error, dims(1), 2_c_int, more="dims(1)"); if (allocated(error)) return
      call check(error, dims(2), 3_c_int, more="dims(2)"); if (allocated(error)) return
      call check(error, dims(3), 5_c_int, more="dims(3)"); if (allocated(error)) return

      ! Check spacing
      call check(error, gridSpacing(1), 1.0_ReKi, more="gridSpacing(1)"); if (allocated(error)) return
      call check(error, gridSpacing(2), 1.0_ReKi, more="gridSpacing(2)"); if (allocated(error)) return
      call check(error, gridSpacing(3), 1.0_ReKi, more="gridSpacing(3)"); if (allocated(error)) return

      ! Check grid size
      call check(error, gridSize(1), 1.0_ReKi, more="gridSize(1)"); if (allocated(error)) return
      call check(error, gridSize(2), 2.0_ReKi, more="gridSize(2)"); if (allocated(error)) return
      call check(error, gridSize(3), 4.0_ReKi, more="gridSize(3)"); if (allocated(error)) return

      ! Check lower bounds
      call check(error, bounds(1,1), 6.5_ReKi, more="origin(1)"); if (allocated(error)) return
      call check(error, bounds(1,2), 6.5_ReKi, more="origin(2)"); if (allocated(error)) return
      call check(error, bounds(1,3), 6.5_ReKi, more="origin(3)"); if (allocated(error)) return

      ! Check upper bounds
      call check(error, bounds(2,1), 7.5_ReKi, more="ub(1)"); if (allocated(error)) return
      call check(error, bounds(2,2), 8.5_ReKi, more="ub(2)"); if (allocated(error)) return
      call check(error, bounds(2,3), 10.5_ReKi, more="ub(3)"); if (allocated(error)) return

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

      integer(IntKi), allocatable :: DirIndices(:)
      integer(IntKi), parameter  :: Expected(0:NumSteps-1) = [0, 6, 12, 18, 24]
      integer(IntKi)             :: ErrStat, i
      character(ErrMsgLen)       :: ErrMsg

      call amrex_find_subvols(DirPath, SubVol, DT, NumSteps, StartIndex, &
                              DirIndices, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_None, more="amrex_find_subvols: "//trim(ErrMsg)); if (allocated(error)) return

      call check(error, lbound(DirIndices,1), 0, more="lbound"); if (allocated(error)) return
      call check(error, ubound(DirIndices,1), NumSteps-1, more="ubound"); if (allocated(error)) return
      do i = 0, NumSteps-1
         call check(error, DirIndices(i), Expected(i), more="step "//trim(Num2LStr(i))); if (allocated(error)) return
      end do

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

      integer(IntKi), allocatable :: DirIndices(:)
      integer(IntKi), parameter  :: Expected(0:NumSteps-1) = [16, 24, 32]
      integer(IntKi)             :: ErrStat, i
      character(ErrMsgLen)       :: ErrMsg

      call amrex_find_subvols(DirPath, SubVol, DT, NumSteps, StartIndex, &
                              DirIndices, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_None, more="amrex_find_subvols: "//trim(ErrMsg)); if (allocated(error)) return

      do i = 0, NumSteps-1
         call check(error, DirIndices(i), Expected(i), more="step "//trim(Num2LStr(i))); if (allocated(error)) return
      end do

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

      integer(IntKi), allocatable :: DirIndices(:)
      integer(IntKi), parameter  :: Expected(0:NumSteps-1) = [6, 18, 30]
      integer(IntKi)             :: ErrStat, i
      character(ErrMsgLen)       :: ErrMsg

      call amrex_find_subvols(DirPath, SubVol, DT, NumSteps, StartIndex, &
                              DirIndices, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_None, more="amrex_find_subvols: "//trim(ErrMsg)); if (allocated(error)) return

      do i = 0, NumSteps-1
         call check(error, DirIndices(i), Expected(i), more="step "//trim(Num2LStr(i))); if (allocated(error)) return
      end do

   end subroutine

   ! Sub-volume directories written with a varying solver time step: the output interval is a
   ! uniform 0.1 s but the step index stride changes from 4 to 2 part way through, as happens when
   ! AMR-Wind transitions from time.initial_dt to fixed_dt. Header times carry the floating-point
   ! drift observed in real data. Before directories were matched by time this failed with
   ! "inconsistent delta between indices '00008' and '00010'".
   subroutine AMReX_test_find_subvols_variable_stride(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter    :: DirPath = "data/subvolvardt"
      integer(IntKi), parameter  :: SubVol = 0
      real(DbKi), parameter      :: DT = 0.1_DbKi
      integer(IntKi), parameter  :: NumSteps = 5
      character(*), parameter    :: StartIndex = "00000"

      integer(IntKi), allocatable :: DirIndices(:)
      integer(IntKi), parameter  :: Expected(0:NumSteps-1) = [0, 4, 8, 10, 12]
      integer(IntKi)             :: ErrStat, i
      character(ErrMsgLen)       :: ErrMsg

      call amrex_find_subvols(DirPath, SubVol, DT, NumSteps, StartIndex, &
                              DirIndices, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_None, more="amrex_find_subvols: "//trim(ErrMsg)); if (allocated(error)) return

      call check(error, lbound(DirIndices,1), 0, more="lbound"); if (allocated(error)) return
      do i = 0, NumSteps-1
         call check(error, DirIndices(i), Expected(i), more="step "//trim(Num2LStr(i))); if (allocated(error)) return
      end do

   end subroutine

   ! Asking for one step more than exists must name the missing step, not report a bare count
   subroutine AMReX_test_find_subvols_missing_step(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter    :: DirPath = "data/subvolvardt"
      integer(IntKi), parameter  :: SubVol = 0
      real(DbKi), parameter      :: DT = 0.1_DbKi
      integer(IntKi), parameter  :: NumSteps = 6
      character(*), parameter    :: StartIndex = "00000"

      integer(IntKi), allocatable :: DirIndices(:)
      integer(IntKi)             :: ErrStat
      character(ErrMsgLen)       :: ErrMsg

      call amrex_find_subvols(DirPath, SubVol, DT, NumSteps, StartIndex, &
                              DirIndices, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_Fatal, more="expected a fatal error"); if (allocated(error)) return
      call check(error, index(ErrMsg, "time step 5") > 0, .true., &
                 more="message should name the missing step: "//trim(ErrMsg)); if (allocated(error)) return

   end subroutine

   ! Two directories whose header times land on the same step is ambiguous and must be rejected
   subroutine AMReX_test_find_subvols_duplicate_step(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter    :: DirPath = "data/subvolvardt"
      integer(IntKi), parameter  :: SubVol = 1
      real(DbKi), parameter      :: DT = 0.1_DbKi
      integer(IntKi), parameter  :: NumSteps = 2
      character(*), parameter    :: StartIndex = "00000"

      integer(IntKi), allocatable :: DirIndices(:)
      integer(IntKi)             :: ErrStat
      character(ErrMsgLen)       :: ErrMsg

      call amrex_find_subvols(DirPath, SubVol, DT, NumSteps, StartIndex, &
                              DirIndices, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_Fatal, more="expected a fatal error"); if (allocated(error)) return
      call check(error, index(ErrMsg, "claim time step 1") > 0, .true., &
                 more="message should report the duplicated step: "//trim(ErrMsg)); if (allocated(error)) return

   end subroutine

   ! A directory sitting off the step grid by more than the tolerance should be reported as a
   ! near miss, so the user is pointed at the file rather than left with a bare count mismatch
   subroutine AMReX_test_find_subvols_out_of_tolerance(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter    :: DirPath = "data/subvolvardt"
      integer(IntKi), parameter  :: SubVol = 2
      real(DbKi), parameter      :: DT = 0.1_DbKi
      integer(IntKi), parameter  :: NumSteps = 2
      character(*), parameter    :: StartIndex = "00000"

      integer(IntKi), allocatable :: DirIndices(:)
      integer(IntKi)             :: ErrStat
      character(ErrMsgLen)       :: ErrMsg

      call amrex_find_subvols(DirPath, SubVol, DT, NumSteps, StartIndex, &
                              DirIndices, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_Fatal, more="expected a fatal error"); if (allocated(error)) return
      call check(error, index(ErrMsg, "outside") > 0, .true., &
                 more="message should flag the near miss: "//trim(ErrMsg)); if (allocated(error)) return

   end subroutine

   ! Directory indices that grow past the width of DirStartIndex. AMReX pads to a minimum width and
   ! widens beyond it, so a six-digit index follows a five-digit one. Comparing the suffixes as text
   ! silently drops those directories, because "100000" sorts before "99998".
   subroutine AMReX_test_find_subvols_wide_index(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter    :: DirPath = "data/subvolwide"
      integer(IntKi), parameter  :: SubVol = 0
      real(DbKi), parameter      :: DT = 0.1_DbKi
      integer(IntKi), parameter  :: NumSteps = 3
      character(*), parameter    :: StartIndex = "99998"

      integer(IntKi), allocatable :: DirIndices(:)
      integer(IntKi), parameter  :: Expected(0:NumSteps-1) = [99998, 100000, 100002]
      integer(IntKi)             :: ErrStat, i
      character(ErrMsgLen)       :: ErrMsg

      call amrex_find_subvols(DirPath, SubVol, DT, NumSteps, StartIndex, &
                              DirIndices, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_None, more="amrex_find_subvols: "//trim(ErrMsg)); if (allocated(error)) return

      do i = 0, NumSteps-1
         call check(error, DirIndices(i), Expected(i), more="step "//trim(Num2LStr(i))); if (allocated(error)) return
      end do

   end subroutine

   ! ---------------------------------------------------------------------------------------
   ! ReadWindAMReX: resolving a time step to a directory
   !
   ! These cover the contract between the index table built during initialization and the
   ! read itself: the table is indexed by the 0-based FAST.Farm time step, and the directory
   ! suffix is zero-padded to at least DirIndexLen characters and widened beyond it when the
   ! index needs more digits.
   ! ---------------------------------------------------------------------------------------

   ! Step n must resolve to the directory the table names for it, not to a computed index.
   ! Table [0, 6, 12] means step 1 is directory 00006; read it and compare against reading
   ! that directory directly.
   subroutine AMReX_test_ReadWindAMReX_lookup(error)
      type(error_type), allocatable, intent(out) :: error
      type(AWAE_ParameterType)   :: p
      real(SiKi), allocatable    :: viaTable(:,:,:,:), direct(:,:,:,:)
      integer(IntKi)             :: ErrStat
      character(ErrMsgLen)       :: ErrMsg

      ! Sub-volume 0 of the subvolmultiple fixture set is a 3x4x5 grid
      allocate(viaTable(3,3,4,5)); allocate(direct(3,3,4,5))

      p%WindFilePath = "data/subvolmultiple"
      p%DirIndexLen  = 5
      allocate(p%DirIndexLow(0:2))
      p%DirIndexLow = [0, 6, 12]

      call ReadWindAMReX(0, 1, p, viaTable, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_None, more="ReadWindAMReX: "//trim(ErrMsg)); if (allocated(error)) return

      call amrex_read_data("data/subvolmultiple_0_00006", direct, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_None, more="amrex_read_data: "//trim(ErrMsg)); if (allocated(error)) return

      call check(error, all(viaTable == direct), .true., &
                 more="step 1 did not resolve to directory 00006"); if (allocated(error)) return

   end subroutine

   ! A step outside the table is a hard error. It must never be clamped to the nearest
   ! entry, which would silently freeze the inflow for the rest of the run.
   subroutine AMReX_test_ReadWindAMReX_out_of_range(error)
      type(error_type), allocatable, intent(out) :: error
      type(AWAE_ParameterType)   :: p
      real(SiKi), allocatable    :: dat(:,:,:,:)
      integer(IntKi)             :: ErrStat
      character(ErrMsgLen)       :: ErrMsg

      allocate(dat(3,3,4,5))
      p%WindFilePath = "data/subvolmultiple"
      p%DirIndexLen  = 5
      allocate(p%DirIndexLow(0:2))
      p%DirIndexLow = [0, 6, 12]

      call ReadWindAMReX(0, 3, p, dat, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_Fatal, more="expected a fatal error"); if (allocated(error)) return
      call check(error, index(ErrMsg, "only steps 0 through 2") > 0, .true., &
                 more="message should report the available range: "//trim(ErrMsg)); if (allocated(error)) return

   end subroutine

   ! An index needing more digits than DirIndexLen must widen the field rather than overflow
   ! it. The subvolwide fixtures run 99998 -> 100000 -> 100002 with a five-character start
   ! index, so step 1 can only be read if the suffix widened to six digits; a fixed-width
   ! write would have produced '*****' and failed to open anything.
   subroutine AMReX_test_ReadWindAMReX_wide_index(error)
      type(error_type), allocatable, intent(out) :: error
      type(AWAE_ParameterType)   :: p
      real(SiKi), allocatable    :: viaTable(:,:,:,:), direct(:,:,:,:)
      integer(IntKi)             :: ErrStat
      character(ErrMsgLen)       :: ErrMsg

      allocate(viaTable(3,3,4,5)); allocate(direct(3,3,4,5))

      p%WindFilePath = "data/subvolwide"
      p%DirIndexLen  = 5
      allocate(p%DirIndexLow(0:2))
      p%DirIndexLow = [99998, 100000, 100002]

      call ReadWindAMReX(0, 1, p, viaTable, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_None, more="ReadWindAMReX: "//trim(ErrMsg)); if (allocated(error)) return

      call amrex_read_data("data/subvolwide_0_100000", direct, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_None, more="amrex_read_data: "//trim(ErrMsg)); if (allocated(error)) return

      call check(error, all(viaTable == direct), .true., &
                 more="six-digit index did not resolve to directory 100000"); if (allocated(error)) return

   end subroutine

   ! If initialization never populated the table, say so instead of reading whatever a
   ! zero-filled lookup would point at.
   subroutine AMReX_test_ReadWindAMReX_no_table(error)
      type(error_type), allocatable, intent(out) :: error
      type(AWAE_ParameterType)   :: p
      real(SiKi), allocatable    :: dat(:,:,:,:)
      integer(IntKi)             :: ErrStat
      character(ErrMsgLen)       :: ErrMsg

      allocate(dat(3,3,4,5))
      p%WindFilePath = "data/subvolmultiple"
      p%DirIndexLen  = 5

      call ReadWindAMReX(0, 0, p, dat, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_Fatal, more="expected a fatal error"); if (allocated(error)) return
      call check(error, index(ErrMsg, "never populated") > 0, .true., &
                 more="message should name the cause: "//trim(ErrMsg)); if (allocated(error)) return

   end subroutine

   ! ---------------------------------------------------------------------------------------
   ! Scan robustness and the fast header path
   ! ---------------------------------------------------------------------------------------

   ! A leftover directory from an earlier run can carry a time far past the window on a LOWER
   ! index than valid data (different time step, same index base). The ascending-index walk
   ! must not stop there while steps are still unclaimed: 00002 has t = 5.0 but 00004 and
   ! 00008 hold steps 1 and 2.
   subroutine AMReX_test_find_subvols_stale_beyond_window(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter    :: DirPath = "data/subvolstale"
      integer(IntKi), parameter  :: SubVol = 0
      real(DbKi), parameter      :: DT = 0.1_DbKi
      integer(IntKi), parameter  :: NumSteps = 3
      character(*), parameter    :: StartIndex = "00000"

      integer(IntKi), allocatable :: DirIndices(:)
      integer(IntKi), parameter  :: Expected(0:NumSteps-1) = [0, 4, 8]
      integer(IntKi)             :: ErrStat, i
      character(ErrMsgLen)       :: ErrMsg

      call amrex_find_subvols(DirPath, SubVol, DT, NumSteps, StartIndex, &
                              DirIndices, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_None, more="amrex_find_subvols: "//trim(ErrMsg)); if (allocated(error)) return
      do i = 0, NumSteps-1
         call check(error, DirIndices(i), Expected(i), more="step "//trim(Num2LStr(i))); if (allocated(error)) return
      end do

   end subroutine

   ! Real Kynema/AMR-Wind sub-volume output layout: one FAB tiling the domain box, three masked
   ! velocity components, headers exactly as written by the sampler. Unlike the hand-built
   ! subvolmultiple fixtures, the Header domain box equals the Cell_H box union here, so the
   ! text-parse fast path engages. The FAB payloads were replaced by a known field
   ! (u = 1+i, v = 2+j, w = 3+k in 0-based cell indices) so the read can be checked exactly.
   subroutine AMReX_test_find_subvols_real_tiled(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter    :: DirPath = "data/subvoltiled"
      integer(IntKi), parameter  :: SubVol = 1
      real(DbKi), parameter      :: DT = 0.1_DbKi
      integer(IntKi), parameter  :: NumSteps = 3
      character(*), parameter    :: StartIndex = "31220"

      integer(IntKi), allocatable :: DirIndices(:)
      integer(IntKi), parameter  :: Expected(0:NumSteps-1) = [31220, 31224, 31228]
      integer(IntKi)             :: ErrStat, i
      character(ErrMsgLen)       :: ErrMsg

      call amrex_find_subvols(DirPath, SubVol, DT, NumSteps, StartIndex, &
                              DirIndices, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_None, more="amrex_find_subvols: "//trim(ErrMsg)); if (allocated(error)) return
      do i = 0, NumSteps-1
         call check(error, DirIndices(i), Expected(i), more="step "//trim(Num2LStr(i))); if (allocated(error)) return
      end do

   end subroutine

   subroutine AMReX_test_read_real_tiled(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter :: DirPath = "data/subvoltiled_1_31220"
      real(SiKi), allocatable :: data(:,:,:,:)
      integer(IntKi)          :: ErrStat, dims(3)
      character(ErrMsgLen)    :: ErrMsg
      real(DbKi)              :: time
      real(ReKi)              :: origin(3), gridSpacing(3)

      call amrex_read_header(DirPath, time, dims, gridSpacing, origin, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_None, more="amrex_read_header: "//trim(ErrMsg)); if (allocated(error)) return
      call check(error, dims(1), 18, more="dims(1)"); if (allocated(error)) return
      call check(error, dims(2), 18, more="dims(2)"); if (allocated(error)) return
      call check(error, dims(3), 19, more="dims(3)"); if (allocated(error)) return
      call check(error, gridSpacing(1), 8.0_ReKi, thr=1.0e-6_ReKi, more="dx"); if (allocated(error)) return

      allocate(data(3, dims(1), dims(2), dims(3)))
      call amrex_read_data(DirPath, data, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_None, more="amrex_read_data: "//trim(ErrMsg)); if (allocated(error)) return
      ! Known payload: u = 1+i, v = 2+j, w = 3+k with 0-based cell indices, checked at the corners
      ! and an interior point. This also pins the component order and the x-fastest layout.
      call check(error, data(1,1,1,1),    1.0_SiKi, thr=1.0e-6_SiKi, more="u(0,0,0)");    if (allocated(error)) return
      call check(error, data(2,1,1,1),    2.0_SiKi, thr=1.0e-6_SiKi, more="v(0,0,0)");    if (allocated(error)) return
      call check(error, data(3,1,1,1),    3.0_SiKi, thr=1.0e-6_SiKi, more="w(0,0,0)");    if (allocated(error)) return
      call check(error, data(1,18,1,1),  18.0_SiKi, thr=1.0e-6_SiKi, more="u(17,0,0)");   if (allocated(error)) return
      call check(error, data(2,1,18,1),  19.0_SiKi, thr=1.0e-6_SiKi, more="v(0,17,0)");   if (allocated(error)) return
      call check(error, data(3,1,1,19),  21.0_SiKi, thr=1.0e-6_SiKi, more="w(0,0,18)");   if (allocated(error)) return
      call check(error, data(1,5,7,9),    5.0_SiKi, thr=1.0e-6_SiKi, more="u(4,6,8)");    if (allocated(error)) return
      call check(error, data(2,5,7,9),    8.0_SiKi, thr=1.0e-6_SiKi, more="v(4,6,8)");    if (allocated(error)) return
      call check(error, data(3,5,7,9),   11.0_SiKi, thr=1.0e-6_SiKi, more="w(4,6,8)");    if (allocated(error)) return

   end subroutine

   ! The fast path is only trusted after the text parse agrees with amrex_read_header on the
   ! starting directory. On real sub-volume output it must agree exactly.
   subroutine AMReX_test_header_text_agrees_on_tiled(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter :: DirPath = "data/subvoltiled_1_31220"
      logical                 :: ok
      integer(IntKi)          :: ErrStat, dimsA(3), dimsT(3), i
      character(ErrMsgLen)    :: ErrMsg
      real(DbKi)              :: timeA, timeT
      real(ReKi)              :: originA(3), dxA(3), originT(3), dxT(3)

      call amrex_read_header(DirPath, timeA, dimsA, dxA, originA, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_None, more="amrex_read_header: "//trim(ErrMsg)); if (allocated(error)) return
      call amrex_parse_header_text(DirPath, ok, timeT, dimsT, dxT, originT)
      call check(error, ok, .true., more="text parse should succeed on real output"); if (allocated(error)) return
      call check(error, timeT, timeA, thr=1.0e-9_DbKi, more="time"); if (allocated(error)) return
      do i = 1, 3
         call check(error, dimsT(i), dimsA(i), more="dims"); if (allocated(error)) return
         call check(error, dxT(i), dxA(i), thr=1.0e-6_ReKi, more="dx"); if (allocated(error)) return
         call check(error, originT(i), originA(i), thr=1.0e-3_ReKi, more="origin"); if (allocated(error)) return
      end do

   end subroutine

   ! The hand-built fixtures store a small box array inside a much larger domain box, so the
   ! text parse (domain) and amrex_read_header (box union) disagree on dims. That disagreement
   ! is exactly what makes the search fall back to the slow path for them -- record it.
   subroutine AMReX_test_header_text_disagrees_on_subset(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter :: DirPath = "data/subvolmultiple_0_00006"
      logical                 :: ok
      integer(IntKi)          :: ErrStat, dimsA(3), dimsT(3)
      character(ErrMsgLen)    :: ErrMsg
      real(DbKi)              :: timeA, timeT
      real(ReKi)              :: originA(3), dxA(3), originT(3), dxT(3)

      call amrex_read_header(DirPath, timeA, dimsA, dxA, originA, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_None, more="amrex_read_header: "//trim(ErrMsg)); if (allocated(error)) return
      call amrex_parse_header_text(DirPath, ok, timeT, dimsT, dxT, originT)
      call check(error, ok, .true., more="text parse should succeed"); if (allocated(error)) return
      call check(error, timeT, timeA, thr=1.0e-9_DbKi, more="time agrees"); if (allocated(error)) return
      call check(error, dimsT(1), 256, more="text dims are the domain box"); if (allocated(error)) return
      call check(error, dimsA(1), 3, more="reader dims are the box union"); if (allocated(error)) return

   end subroutine

   ! The scan lists the parent directory and matches its entries against the prefix. The entry
   ! paths the iterator hands back are not textually the prefix the caller passed -- a prefix with
   ! no directory of its own comes back as "./name", and a redundant separator is normalized away
   ! -- so the comparison has to be made on the final path component. Matching the whole path
   ! instead finds nothing beyond the starting directory even though every directory is present.
   subroutine AMReX_test_find_subvols_prefix_form(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter    :: DirPath = "data//subvolmultiple"
      integer(IntKi), parameter  :: SubVol = 0
      real(DbKi), parameter      :: DT = 0.6_DbKi
      integer(IntKi), parameter  :: NumSteps = 5
      character(*), parameter    :: StartIndex = "00000"

      integer(IntKi), allocatable :: DirIndices(:)
      integer(IntKi), parameter  :: Expected(0:NumSteps-1) = [0, 6, 12, 18, 24]
      integer(IntKi)             :: ErrStat, i
      character(ErrMsgLen)       :: ErrMsg

      call amrex_find_subvols(DirPath, SubVol, DT, NumSteps, StartIndex, &
                              DirIndices, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_None, more="amrex_find_subvols: "//trim(ErrMsg)); if (allocated(error)) return
      do i = 0, NumSteps-1
         call check(error, DirIndices(i), Expected(i), more="step "//trim(Num2LStr(i))); if (allocated(error)) return
      end do

   end subroutine

   ! amrex_read_data writes the destination by grid index, so a plotfile describing a different
   ! grid must be rejected rather than written past the end of the caller's array. The fixture
   ! grid is 3x4x5; ask for it into a 3x4x4 array.
   subroutine AMReX_test_read_data_wrong_dims(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter :: DirPath = "data/subvolmultiple_0_00006"
      real(SiKi), allocatable :: dat(:,:,:,:)
      integer(IntKi)          :: ErrStat
      character(ErrMsgLen)    :: ErrMsg

      allocate(dat(3,3,4,4))
      call amrex_read_data(DirPath, dat, ErrStat, ErrMsg)
      call check(error, ErrStat, ErrID_Fatal, more="expected a fatal error"); if (allocated(error)) return
      call check(error, index(ErrMsg, "do not match") > 0, .true., &
                 more="message should report the dimension mismatch: "//trim(ErrMsg)); if (allocated(error)) return

   end subroutine

end module
