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

   subroutine amrex_read_data_c(dir_path, data, err_stat, err_msg, err_msg_len) bind(c)
      import
      implicit none
      character(kind=c_char), intent(in)    :: dir_path(*)
      real(kind=c_float), intent(out)       :: data(*)
      integer(kind=c_int), intent(out)      :: err_stat
      character(kind=c_char), intent(out)   :: err_msg(*)
      integer(kind=c_int), intent(in)       :: err_msg_len
   end subroutine
   
   subroutine amrex_find_subvols_c(dir_path, subvol, dt, num_step, start_index, first_index, &
                                   index_delta, err_stat, err_msg, err_msg_len) bind(c)
      import
      implicit none
      character(kind=c_char), intent(in)    :: dir_path(*)
      integer(kind=c_int), intent(in)       :: subvol
      character(kind=c_char), intent(in)    :: start_index(*)
      real(kind=c_double), intent(in)       :: dt
      integer(kind=c_int), intent(in)       :: num_step
      integer(kind=c_int), intent(out)      :: first_index
      integer(kind=c_int), intent(out)      :: index_delta
      integer(kind=c_int), intent(out)      :: err_stat
      character(kind=c_char), intent(out)   :: err_msg(*)
      integer(kind=c_int), intent(in)       :: err_msg_len
   end subroutine
end interface

public :: amrex_init, amrex_finalize
public :: amrex_read_header, amrex_read_data, amrex_find_subvols

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

#ifdef ENABLE_AMREX_LIB

   ! Convert directory path to C type
   dir_path = transfer(trim(DirPath) // c_null_char, dir_path) 

   ! Call C++ function to read header
   call amrex_read_header_c(dir_path, t, dims, gridSpacing, origin, err_stat_c, err_msg_c, ErrMsgLen)

   ! Transfer outputs back to fortran types
   time = transfer(t, time)
   nXYZ = transfer(dims, nXYZ)
   dXYZ = transfer(gridSpacing, dXYZ)
   oXYZ = transfer(origin, oXYZ)
   ErrStat = transfer(err_stat_c, ErrStat)
   ErrMsg = transfer(err_msg_c, ErrMsg)
   i = index(ErrMsg, c_null_char)
   if (i > 0) ErrMsg = ErrMsg(:i-1)
   
#else
   call SetErrStat(ErrID_Fatal, "Enable AMReX library with -DAMREX_READER", ErrStat, ErrMsg, "amrex_read_header")
#endif
end subroutine

! Read the XYZ velocity grid data into the FAST.Farm ambient wind data array [XYZ,NX,NY,NZ].
! This function cannot be called in parallel due to internal restrictions of the AMReX library.
subroutine amrex_read_data(DirPath, gridData, ErrStat, ErrMsg)
   character(*), intent(in)      :: DirPath
   real(SiKi), intent(out)       :: gridData(:,:,:,:)
   integer(IntKi), intent(out)   :: ErrStat
   character(*), intent(out)     :: ErrMsg

   character(c_char), allocatable :: dir_path(:)
   integer(c_int)    :: err_stat_c
   character(c_char) :: err_msg_c(ErrMsgLen)
   integer(IntKi)    :: i

#ifdef ENABLE_AMREX_LIB

   ! Convert directory path to C type
   dir_path = transfer(trim(DirPath) // c_null_char, dir_path) 

   ! Call C++ function to read header
   call amrex_read_data_c(dir_path, gridData, err_stat_c, err_msg_c, ErrMsgLen)

   ! Transfer outputs back to fortran types
   ErrStat = transfer(err_stat_c, ErrStat)
   ErrMsg = transfer(err_msg_c, ErrMsg)
   i = index(ErrMsg, c_null_char)
   if (i > 0) ErrMsg = ErrMsg(:i-1)

#else
   call SetErrStat(ErrID_Fatal, "Enable AMReX library with -DAMREX_READER", ErrStat, ErrMsg, "amrex_read_data")
#endif
end subroutine

! Search for AMReX directories based on given directory prefix, sub-volume 
! number, time step, total number of steps, and the starting index string 
! (e.g. `00000`). This function returns first index as a number, and the delta 
! between successive directory indices. It also checks that a sufficient number 
! of directories are available, that the data matches the requested time step, 
! and the grid properties are consistent (size, origin, spacing).
subroutine amrex_find_subvols(DirPath, SubVol, DT, NumStep, StartIndex, &
                              FirstIndex, IndexDelta, ErrStat, ErrMsg)
   character(*), intent(in)      :: DirPath
   integer(IntKi), intent(in)    :: SubVol
   real(DbKi), intent(in)        :: DT          ! Time step
   integer(IntKi), intent(in)    :: NumStep     ! Number of steps
   character(*), intent(in)      :: StartIndex
   integer(IntKi), intent(out)   :: FirstIndex, IndexDelta
   integer(IntKi), intent(out)   :: ErrStat
   character(*), intent(out)     :: ErrMsg

   character(c_char), allocatable :: dir_path(:)
   character(c_char), allocatable :: start_index(:)
   real(c_double)                 :: dt_c
   integer(c_int)                 :: num_step, first_index, index_delta, subvol_c
   integer(c_int)                 :: err_stat_c
   character(c_char)              :: err_msg_c(ErrMsgLen)
   integer(IntKi)                 :: i

#ifdef ENABLE_AMREX_LIB

   ! Convert directory path to C type
   dir_path = transfer(trim(DirPath) // c_null_char, dir_path) 

   ! Convert start directory to C type
   start_index = transfer(trim(StartIndex) // c_null_char, start_index) 
   
   ! Transfer inputs to C types
   subvol_c = transfer(SubVol, subvol_c)
   num_step = transfer(NumStep, num_step)
   dt_c = transfer(DT, dt_c)

   ! Call C++ routine to find first index and index delta values
   call amrex_find_subvols_c(dir_path, subvol_c, dt_c, num_step, start_index, first_index, &
                             index_delta, err_stat_c, err_msg_c, ErrMsgLen)

   ! Transfer outputs to fortran types
   FirstIndex = transfer(first_index, FirstIndex)
   IndexDelta = transfer(index_delta, IndexDelta)
   ErrStat = transfer(err_stat_c, ErrStat)
   ErrMsg = transfer(err_msg_c, ErrMsg)
   i = index(ErrMsg, c_null_char)
   if (i > 0) ErrMsg = ErrMsg(:i-1)

#else
   call SetErrStat(ErrID_Fatal, "Enable AMReX library with -DAMREX_READER", ErrStat, ErrMsg, "amrex_find_subvols")
#endif
end subroutine

end module
