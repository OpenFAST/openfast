module amrex_utils
use iso_c_binding
use NWTC_Library

#ifdef ENABLE_AMREX_LIB
use amrex_base_module, only: amrex_init, amrex_finalize
#endif

implicit none

private

interface
   subroutine amrex_read_header_c(dir_path, t, dims, dx, origin, err_stat, err_msg, err_msg_len) bind(c)
      import
      implicit none
      character(kind=c_char), intent(in)    :: dir_path(*)
      real(kind=c_double), intent(out)      :: t
      integer(kind=c_int), intent(out)      :: dims(3)
      real(kind=c_double), intent(out)      :: dx(3)
      real(kind=c_double), intent(out)      :: origin(3)
      integer(kind=c_int), intent(out)      :: err_stat
      character(kind=c_char), intent(out)   :: err_msg(*)
      integer(kind=c_int), intent(in)       :: err_msg_len
   end subroutine

   subroutine amrex_read_data_c(dir_path, data, dims_expected, err_stat, err_msg, err_msg_len) bind(c)
      import
      implicit none
      character(kind=c_char), intent(in)    :: dir_path(*)
      real(kind=c_float), intent(out)       :: data(*)
      integer(kind=c_int), intent(in)       :: dims_expected(3)
      integer(kind=c_int), intent(out)      :: err_stat
      character(kind=c_char), intent(out)   :: err_msg(*)
      integer(kind=c_int), intent(in)       :: err_msg_len
   end subroutine
   
   subroutine amrex_header_text_c(dir_path, t, dims, dx, origin, ok) bind(c)
      import
      implicit none
      character(kind=c_char), intent(in)    :: dir_path(*)
      real(kind=c_double), intent(out)      :: t
      integer(kind=c_int), intent(out)      :: dims(3)
      real(kind=c_double), intent(out)      :: dx(3)
      real(kind=c_double), intent(out)      :: origin(3)
      integer(kind=c_int), intent(out)      :: ok
   end subroutine

   subroutine amrex_find_subvols_c(dir_path, subvol, dt, num_step, start_index, dir_indices, &
                                   err_stat, err_msg, err_msg_len) bind(c)
      import
      implicit none
      character(kind=c_char), intent(in)    :: dir_path(*)
      integer(kind=c_int), intent(in)       :: subvol
      character(kind=c_char), intent(in)    :: start_index(*)
      real(kind=c_double), intent(in)       :: dt
      integer(kind=c_int), intent(in)       :: num_step
      integer(kind=c_int), intent(out)      :: dir_indices(*)
      integer(kind=c_int), intent(out)      :: err_stat
      character(kind=c_char), intent(out)   :: err_msg(*)
      integer(kind=c_int), intent(in)       :: err_msg_len
   end subroutine
end interface

public :: amrex_init, amrex_finalize
public :: amrex_read_header, amrex_read_data, amrex_find_subvols, amrex_parse_header_text

contains

#ifndef ENABLE_AMREX_LIB

subroutine amrex_init(arg_parmparse)
   logical, optional, intent(in) :: arg_parmparse
end subroutine

subroutine amrex_finalize()
end subroutine

#endif

! Read the header information for the AMReX grid and return it.
subroutine amrex_read_header(DirPath, time, nXYZ, dXYZ, oXYZ, ErrStat, ErrMsg)
   character(*), intent(in)      :: DirPath
   real(DbKi), intent(out)       :: Time
   Integer(IntKi), intent(out)   :: nXYZ(3)
   real(ReKi), intent(out)       :: dXYZ(3)
   real(ReKi), intent(out)       :: oXYZ(3)
   integer(IntKi), intent(out)   :: ErrStat
   character(*), intent(out)     :: ErrMsg

   character(c_char), allocatable :: dir_path(:)
   integer(c_int)    :: err_stat_c
   character(c_char) :: err_msg_c(ErrMsgLen)
   integer(c_int)    :: dims(3)
   real(c_double)    :: t, origin(3), gridSpacing(3)
   integer(IntKi)    :: i

   ErrStat = ErrID_None
   ErrMsg  = ""

#ifdef ENABLE_AMREX_LIB

   ! Convert directory path to C type
   allocate(dir_path(len_trim(DirPath) + 1))
   dir_path = transfer(trim(DirPath) // c_null_char, dir_path) 

   ! Call C++ function to read header
   call amrex_read_header_c(dir_path, t, dims, gridSpacing, origin, err_stat_c, err_msg_c, ErrMsgLen)

   ! Transfer outputs back to fortran types
   time = real(t, DbKi)
   nXYZ = int(dims, IntKi)
   dXYZ = real(gridSpacing, ReKi)
   oXYZ = real(origin, ReKi)
   ErrStat = int(err_stat_c, IntKi)
   ErrMsg = transfer(err_msg_c, ErrMsg)
   i = index(ErrMsg, c_null_char)
   if (i > 0) ErrMsg = ErrMsg(:i-1)
   
#else
   call SetErrStat(ErrID_Fatal, "AMReX library unavailable. Enable with -DAMREX_READER during compile with cmake on Linux, or change FAST.Farm Mod_AmbWind type", ErrStat, ErrMsg, "amrex_read_header")
#endif
end subroutine

! Read the XYZ velocity grid data into the FAST.Farm ambient wind data array [XYZ,NX,NY,NZ].
! This function cannot be called in parallel due to internal restrictions of the AMReX library.
subroutine amrex_read_data(DirPath, gridData, ErrStat, ErrMsg)
   character(*), intent(in)      :: DirPath
   real(SiKi), intent(out)       :: gridData(:,:,:,:)
   integer(IntKi), intent(out)   :: ErrStat
   character(*), intent(out)     :: ErrMsg

   character(*), parameter        :: RoutineName = 'amrex_read_data'
   character(c_char), allocatable :: dir_path(:)
   integer(c_int)    :: err_stat_c
   integer(c_int)    :: dims_expected(3)
   character(c_char) :: err_msg_c(ErrMsgLen)
   integer(IntKi)    :: i

   ErrStat = ErrID_None
   ErrMsg  = ""

#ifdef ENABLE_AMREX_LIB

   ! The C routine writes gridData by grid index, so it needs the extent of the destination to
   ! stay inside it. The first dimension is the three velocity components.
   if (size(gridData, 1) /= 3) then
      call SetErrStat(ErrID_Fatal, "gridData must have 3 velocity components, got "// &
                      trim(Num2LStr(size(gridData, 1))), ErrStat, ErrMsg, RoutineName)
      return
   end if
   dims_expected = int([size(gridData, 2), size(gridData, 3), size(gridData, 4)], c_int)

   ! Convert directory path to C type
   allocate(dir_path(len_trim(DirPath) + 1))
   dir_path = transfer(trim(DirPath) // c_null_char, dir_path) 

   ! Call C++ function to read the grid data
   call amrex_read_data_c(dir_path, gridData, dims_expected, err_stat_c, err_msg_c, ErrMsgLen)

   ! Transfer outputs back to fortran types
   ErrStat = int(err_stat_c, IntKi)
   ErrMsg = transfer(err_msg_c, ErrMsg)
   i = index(ErrMsg, c_null_char)
   if (i > 0) ErrMsg = ErrMsg(:i-1)

#else
   call SetErrStat(ErrID_Fatal, "AMReX library unavailable. Enable with -DAMREX_READER during compile with cmake on Linux, or change FAST.Farm Mod_AmbWind type", ErrStat, ErrMsg, "amrex_read_data")
#endif
end subroutine

! Search for AMReX directories based on given directory prefix, sub-volume
! number, time step, total number of steps, and the starting index string
! (e.g. `00000`). Returns DirIndices(0:NumStep-1), the directory index suffix
! to use for each time step. Directories are matched to time steps by the
! simulation time recorded in their Header; no constant stride between
! successive directory indices is assumed, so LES output written with a varying
! solver time step is supported. Also checks that every requested step is
! represented -- and, among the directories scanned, represented only once --
! and that the grid properties are consistent (size, origin, spacing).
subroutine amrex_find_subvols(DirPath, SubVol, DT, NumStep, StartIndex, &
                              DirIndices, ErrStat, ErrMsg)
   character(*), intent(in)                   :: DirPath
   integer(IntKi), intent(in)                 :: SubVol
   real(DbKi), intent(in)                     :: DT          ! Time step
   integer(IntKi), intent(in)                 :: NumStep     ! Number of steps
   character(*), intent(in)                   :: StartIndex
   integer(IntKi), allocatable, intent(out)   :: DirIndices(:)  ! (0:NumStep-1)
   integer(IntKi), intent(out)                :: ErrStat
   character(*), intent(out)                  :: ErrMsg

   character(*), parameter        :: RoutineName = 'amrex_find_subvols'
   character(c_char), allocatable :: dir_path(:)
   character(c_char), allocatable :: start_index(:)
   integer(c_int), allocatable    :: dir_indices_c(:)
   real(c_double)                 :: dt_c
   integer(c_int)                 :: num_step, subvol_c
   integer(c_int)                 :: err_stat_c
   character(c_char)              :: err_msg_c(ErrMsgLen)
   integer(IntKi)                 :: i, stat

   ErrStat = ErrID_None
   ErrMsg  = ""

#ifdef ENABLE_AMREX_LIB

   if (NumStep < 1) then
      call SetErrStat(ErrID_Fatal, "number of time steps must be at least 1", ErrStat, ErrMsg, RoutineName)
      return
   end if

   ! Convert directory path to C type
   allocate(dir_path(len_trim(DirPath) + 1))
   dir_path = transfer(trim(DirPath) // c_null_char, dir_path) 

   ! Convert start directory to C type
   allocate(start_index(len_trim(StartIndex) + 1))
   start_index = transfer(trim(StartIndex) // c_null_char, start_index) 
   
   ! Transfer inputs to C types
   subvol_c = int(SubVol, c_int)
   num_step = int(NumStep, c_int)
   dt_c = real(DT, c_double)

   ! Receive buffer for the C call. Allocated 1-based here and rebased to 0 on the way out,
   ! since callers index the result with the (0-based) FAST.Farm time step number.
   allocate(dir_indices_c(NumStep), stat=stat)
   if (stat /= 0) then
      call SetErrStat(ErrID_Fatal, "error allocating dir_indices_c", ErrStat, ErrMsg, RoutineName)
      return
   end if
   dir_indices_c = -1_c_int

   ! Call C++ routine to build the time step -> directory index table
   call amrex_find_subvols_c(dir_path, subvol_c, dt_c, num_step, start_index, dir_indices_c, &
                             err_stat_c, err_msg_c, ErrMsgLen)

   ! Transfer outputs to fortran types
   ErrStat = int(err_stat_c, IntKi)
   ErrMsg = transfer(err_msg_c, ErrMsg)
   i = index(ErrMsg, c_null_char)
   if (i > 0) ErrMsg = ErrMsg(:i-1)

   ! Leave DirIndices unallocated on failure so a caller that ignores ErrStat trips the
   ! allocated() guard in ReadWindAMReX rather than silently reading zeros.
   if (ErrStat >= AbortErrLev) return

   allocate(DirIndices(0:NumStep-1), stat=stat)
   if (stat /= 0) then
      call SetErrStat(ErrID_Fatal, "error allocating DirIndices", ErrStat, ErrMsg, RoutineName)
      return
   end if
   DirIndices = int(dir_indices_c, IntKi)

#else
   call SetErrStat(ErrID_Fatal, "AMReX library unavailable. Enable with -DAMREX_READER during compile with cmake on Linux, or change FAST.Farm Mod_AmbWind type", ErrStat, ErrMsg, RoutineName)
#endif
end subroutine

! Parse a plotfile Header directly, without AMReX, and return the grid information it holds.
! This is the fast path the sub-volume search uses once it has verified, on the starting
! directory, that the text agrees with amrex_read_header. Exposed so that agreement can be
! tested on real plotfiles. Ok is .false. if the header could not be parsed or describes a
! plotfile this reader does not accept.
subroutine amrex_parse_header_text(DirPath, Ok, Time, nXYZ, dXYZ, oXYZ)
   character(*), intent(in)      :: DirPath
   logical, intent(out)          :: Ok
   real(DbKi), intent(out)       :: Time
   integer(IntKi), intent(out)   :: nXYZ(3)
   real(ReKi), intent(out)       :: dXYZ(3)
   real(ReKi), intent(out)       :: oXYZ(3)

   character(c_char), allocatable :: dir_path(:)
   integer(c_int)    :: ok_c, dims(3)
   real(c_double)    :: t, origin(3), gridSpacing(3)

   Ok = .false.; Time = 0.0_DbKi; nXYZ = 0; dXYZ = 0.0_ReKi; oXYZ = 0.0_ReKi

#ifdef ENABLE_AMREX_LIB
   allocate(dir_path(len_trim(DirPath) + 1))
   dir_path = transfer(trim(DirPath) // c_null_char, dir_path)
   call amrex_header_text_c(dir_path, t, dims, gridSpacing, origin, ok_c)
   Ok = (ok_c /= 0)
   if (Ok) then
      Time = real(t, DbKi)
      nXYZ = int(dims, IntKi)
      dXYZ = real(gridSpacing, ReKi)
      oXYZ = real(origin, ReKi)
   end if
#endif
end subroutine

end module
