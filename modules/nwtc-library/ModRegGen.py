
import textwrap
import itertools
from itertools import product

type_map = {
    'C1': 'character(*)',
    'L1': 'logical',
    'I4': 'integer(B4Ki)',
    'I8': 'integer(B8Ki)',
    'R4': 'real(R4Ki)',
    'R8': 'real(R8Ki)',
}

num_ranks = 5

module_header = '''
!STARTOFGENERATEDFILE 'ModReg.f90'
!
! WARNING This file is generated automatically by ModRegGen.py.
! Do not edit.  Your changes to this file will be lost.
!
! FAST Registry
!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2024  National Renewable Energy Laboratory
!
!    This file is part of the NWTC Subroutine Library.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!**********************************************************************************************************************************

!> This module contains routines for packing and unpacking data from a registry data file.
module ModReg
   use NWTC_Base
   implicit none

   private
   public :: RegFile
   public :: OpenRegFile, InitRegFile, CloseRegFile, RegCheckErr
   public :: RegPackBounds, RegUnpackBounds
   public :: RegPackPointer, RegUnpackPointer
   public :: RegPack, RegUnpack
   public :: RegPackAlloc, RegUnpackAlloc
   public :: RegPackPtr, RegUnpackPtr

   type :: RegFile
      integer(IntKi)             :: Unit
      integer(IntKi)             :: Offset
      type(c_ptr), allocatable   :: Pointers(:)
      integer(B8Ki)              :: NumData
      integer(B8Ki)              :: NumPointers
      integer(IntKi)             :: ErrStat = ErrID_Fatal
      character(ErrMsgLen)       :: ErrMsg = 'RegFile not initialized'
   end type
   {ifc_lines}

contains

   subroutine InitRegFile(RF, Unit, ErrStat, ErrMsg)
      type(RegFile), intent(inout)        :: RF
      integer(IntKi), intent(in)          :: Unit
      integer(IntKi), intent(out)         :: ErrStat
      character(ErrMsgLen), intent(out)   :: ErrMsg

      character(*), parameter             :: RoutineName = "InitRegFile"
      integer(B8Ki), parameter            :: NumPointersInit = 128
      integer(IntKi)                      :: stat

      ErrStat = ErrID_None
      ErrMsg = ""

      RF%ErrStat = ErrID_None
      RF%ErrMsg = ""
      RF%NumData = 0
      RF%NumPointers = 0
      RF%Unit = Unit

      ! Get current position in the file in case anything has been written to it
      inquire(Unit, POS=RF%Offset)

      ! Write invalid number of pointers at the beginning of file so we can
      ! check if the file if the file has been finalized and closed
      write (Unit, iostat=stat) -1_B8Ki
      if (stat /= 0) then
         ErrStat = ErrID_Fatal
         write (ErrMsg, *) 'InitRegFile: Unable to write offset at beginning of file'
         return
      end if

      ! If pointers have not been allocated, allocate with initial size
      if (.not. allocated(RF%Pointers)) then
         allocate (RF%Pointers(NumPointersInit), stat=stat)
         if (stat /= 0) then
            ErrStat = ErrID_Fatal
            write (ErrMsg, *) 'InitRegFile: Unable to init pointer index to with size of', NumPointersInit
            return
         end if
      end if

      ! Reset all pointers to null
      RF%Pointers = c_null_ptr
   end subroutine

   subroutine CloseRegFile(RF, ErrStat, ErrMsg)
      type(RegFile), intent(inout)        :: RF
      integer(IntKi), intent(out)         :: ErrStat
      character(ErrMsgLen), intent(out)   :: ErrMsg

      character(*), parameter             :: RoutineName = "CloseRegFile"
      integer(IntKi)                      :: stat

      ErrStat = ErrID_None
      ErrMsg = ""

      ! Check if there have been any errors while writing to the file
      if (RF%ErrStat /= ErrID_None) then
         call SetErrStat(RF%ErrStat, RF%ErrMsg, ErrStat, ErrMsg, RoutineName)
         return
      end if

      ! Write the actual number of pointers
      write (RF%Unit, POS=RF%Offset, iostat=stat) RF%NumPointers
      if (stat /= 0) then
         ErrStat = ErrID_Fatal
         write (ErrMsg, *) 'CloseRegFile: Unable to write offset at beginning of file'
         return
      end if

      ! Close the file
      close (RF%Unit)

      ! Deallocate pointer array
      if (allocated(RF%Pointers)) deallocate (RF%Pointers)
   end subroutine

   subroutine OpenRegFile(RF, Unit, ErrStat, ErrMsg)
      type(RegFile), intent(inout)        :: RF
      integer(IntKi), intent(in)          :: Unit
      integer(IntKi), intent(out)         :: ErrStat
      character(ErrMsgLen), intent(out)   :: ErrMsg

      character(*), parameter             :: RoutineName = "ReadRegFile"
      integer(IntKi)                      :: iostat

      ErrStat = ErrID_None
      ErrMsg = ''

      ! Save unit
      RF%Unit = Unit

      ! Read number of pointers
      read (Unit, iostat=iostat) RF%NumPointers
      if (iostat /= 0) then
         call SetErrStat(ErrID_Fatal, "Error reading number of pointers", ErrStat, ErrMsg, RoutineName)
         return
      end if

      ! If pointers are allocated, deallocate
      if (allocated(RF%Pointers)) deallocate (RF%Pointers)

      ! Allocate pointer index and initialize pointers to null
      allocate (RF%Pointers(1:RF%NumPointers), stat=ErrStat)
      RF%Pointers = c_null_ptr

      ! initialize the number of data
      RF%NumData = 0

      ! Clear error
      RF%ErrStat = ErrID_None
      RF%ErrMsg = ''
   end subroutine

   function RegCheckErr(RF, RoutineName) result(Err)
      type(RegFile), intent(inout)  :: RF
      character(*), intent(in)      :: RoutineName
      logical                       :: Err
      Err = RF%ErrStat /= ErrID_None
      if (Err) RF%ErrMsg = trim(RoutineName)//": "//trim(RF%ErrMsg)
   end function

   elemental function LogicalToByte(b) result(i)
      logical, intent(in) :: b
      integer(B1Ki)       :: i
      if (b) then
         i = 1_B1Ki
      else
         i = 0_B1Ki
      end if
   end function

   elemental function ByteToLogical(i) result(b)
      integer(B1Ki), intent(in)  :: i
      logical                    :: b
      if (i == 0) then
         b = .false.
      else
         b = .true.
      end if
   end function

   subroutine RegPackPointer(RF, Ptr, Found)
      type(RegFile), intent(inout)  :: RF
      type(c_ptr), intent(in)       :: Ptr
      logical, intent(out)          :: Found

      type(c_ptr), allocatable      :: PointersTmp(:)
      integer(B8Ki)                 :: NewSize
      integer(B8Ki)                 :: i

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Look for pointer in index, if found, pack pointer index and return
      do i = 1, RF%NumPointers
         if (c_associated(Ptr, RF%Pointers(i))) then
            call RegPack(RF, i)
            Found = .true.
            return
         end if
      end do

      ! Pointer was not found in index
      Found = .false.

      ! If pointer index is full, grow pointer index
      if (RF%NumPointers == size(RF%Pointers)) then
         NewSize = int(1.5_R8Ki*real(RF%NumPointers, R8Ki), B8Ki)
         call move_alloc(RF%Pointers, PointersTmp)
         allocate (RF%Pointers(NewSize), stat=RF%ErrStat)
         if (RF%ErrStat /= ErrID_None) then
            RF%ErrStat = ErrID_Fatal
            write (RF%ErrMsg, *) 'RegPackPointer: Unable to allocate pointer index to', NewSize, 'bytes'
            return
         end if
         RF%Pointers(1:size(PointersTmp)) = PointersTmp
         RF%Pointers(size(PointersTmp) + 1:) = c_null_ptr
      end if

      ! Increment number of pointers, add new pointer to index
      RF%NumPointers = RF%NumPointers + 1
      RF%Pointers(RF%NumPointers) = Ptr

      ! Pack pointer index
      call RegPack(RF, RF%NumPointers)
   end subroutine

   subroutine RegUnpackPointer(RF, Ptr, Idx)
      type(RegFile), intent(inout)  :: RF
      type(c_ptr), intent(out)      :: Ptr
      integer(B8Ki), intent(out)    :: Idx

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Unpack pointer index
      call RegUnpack(RF, Idx)

      ! Get pointer from index
      Ptr = RF%Pointers(Idx)
   end subroutine

   subroutine RegPackBounds(RF, R, LB, UB)
      type(RegFile), intent(inout)  :: RF
      integer(B4Ki), intent(in)     :: R
      integer(B4Ki), intent(in)     :: LB(:), UB(:)

      ! If has an error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Pack lower and upper bounds
      call RegPack(RF, LB(1:R))
      call RegPack(RF, UB(1:R))
      if (RegCheckErr(RF, "RegPackBounds")) return
   end subroutine

   subroutine RegUnpackBounds(RF, R, LB, UB)
      type(RegFile), intent(inout)  :: RF
      integer(B4Ki), intent(in)     :: R
      integer(B4Ki), intent(out)    :: LB(:), UB(:)

      ! If has an error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Unpack lower and upper bounds
      call RegUnpack(RF, LB(1:R))
      call RegUnpack(RF, UB(1:R))
      if (RegCheckErr(RF, "RegUnpackBounds")) return
   end subroutine

   function DataNumValid(RF) result(match)
      type(RegFile), intent(inout)  :: RF
      logical                       :: match
      integer(B8Ki)                 :: DataNum

      ! Increment the data number to be read
      RF%NumData = RF%NumData + 1

      ! Read the data number from the file
      read(RF%Unit) DataNum

      ! If data number from file does not match expected number, set match false
      ! and create error message; otherwise, set match to true
      if (DataNum /= RF%NumData) then
         match = .false.
         RF%ErrStat = ErrID_Fatal
         write(RF%ErrMsg, *) "Read data number", DataNum, "expected", RF%NumData
      else
         match = .true.
      end if
   end function
'''


def gen_pack(w, dt, decl, rank):
   dims = '' if rank == 0 else '('+','.join([':']*rank)+')'
   name = f'Pack_{dt}' if rank == 0 else f'Pack_{dt}_Rank{rank}'
   w.write(f'\n\n   subroutine {name}(RF, Data)')
   w.write(f'\n      type(RegFile), intent(inout)         :: RF')
   w.write(f'\n      {decl+", intent(in)":<35s}  :: Data{dims}')
   w.write(f'\n')
   w.write(f'\n      ! If error, return')
   w.write(f'\n      if (RF%ErrStat /= ErrID_None) return')
   w.write(f'\n')
   w.write(f'\n      ! Increment data number and write to file')
   w.write(f'\n      RF%NumData = RF%NumData + 1')
   w.write(f'\n      write(RF%Unit) RF%NumData')
   w.write(f'\n')
   w.write(f'\n      ! Write data to file')
   w.write(f'\n      write(RF%Unit) Data')
   w.write(f'\n   end subroutine')


def gen_unpack(w, dt, decl, rank):
   dims = '' if rank == 0 else '('+','.join([':']*rank)+')'
   name = f'Unpack_{dt}' if rank == 0 else f'Unpack_{dt}_Rank{rank}'
   w.write(f'\n')
   w.write(f'\n   subroutine {name}(RF, Data)')
   w.write(f'\n      type(RegFile), intent(inout)         :: RF')
   w.write(f'\n      {decl+", intent(out)":<35s}  :: Data{dims}')
   w.write(f'\n')
   w.write(f'\n      ! If error, return')
   w.write(f'\n      if (RF%ErrStat /= ErrID_None) return')
   w.write(f'\n')
   w.write(f'\n      ! Read data number, return if invalid')
   w.write(f'\n      if (.not. DataNumValid(RF)) return')
   w.write(f'\n')
   w.write(f'\n      ! Read data from file')
   w.write(f'\n      read(RF%Unit) Data')
   w.write(f'\n   end subroutine')

def gen_pack_alloc(w, dt, decl, rank):
   dims = '' if rank == 0 else '('+','.join([':']*rank)+')'
   name = f'PackAlloc_{dt}' + ("" if rank == 0 else f'_Rank{rank}')
   w.write(f'\n')
   w.write(f'\n   subroutine {name}(RF, Data)')
   w.write(f'\n      type(RegFile), intent(inout)         :: RF')
   w.write(f'\n      {decl+", allocatable, intent(in)":<35s}  :: Data{dims}')
   w.write(f'\n')
   w.write(f'\n      ! If error, return')
   w.write(f'\n      if (RF%ErrStat /= ErrID_None) return')
   w.write(f'\n')
   w.write(f'\n      ! Write if allocated')
   w.write(f'\n      call RegPack(RF, allocated(Data))')
   w.write(f'\n      if (RegCheckErr(RF, "{name}")) return')
   w.write(f'\n      if (.not. allocated(Data)) return')
   w.write(f'\n')
   if rank > 0:
      w.write(f'\n      ! Write array bounds')
      w.write(f'\n      call RegPackBounds(RF, {rank}, lbound(Data), ubound(Data))')
   w.write(f'\n')
   w.write(f'\n      ! Write data to file')
   w.write(f'\n      call RegPack(RF, Data)')
   w.write(f'\n      if (RegCheckErr(RF, "{name}")) return')
   w.write(f'\n   end subroutine')


def gen_unpack_alloc(w, dt, decl, rank):
   dims = '' if rank == 0 else '('+','.join([':']*rank)+')'
   name = f'UnpackAlloc_{dt}' + ("" if rank == 0 else f'_Rank{rank}')
   w.write(f'\n')
   w.write(f'\n   subroutine {name}(RF, Data)')
   w.write(f'\n      type(RegFile), intent(inout)         :: RF')
   w.write(f'\n      {decl+", allocatable, intent(out)":<35s}  :: Data{dims}')
   w.write(f'\n      integer(IntKi)                       :: stat')
   w.write(f'\n      logical                              :: IsAllocated')
   if rank > 0:
      w.write(f'\n      integer(B4Ki)                        :: LB({rank}), UB({rank})')
   w.write(f'\n')
   w.write(f'\n      ! If error, return')
   w.write(f'\n      if (RF%ErrStat /= ErrID_None) return')
   w.write(f'\n')
   w.write(f'\n      ! Deallocate if allocated')
   w.write(f'\n      if (allocated(Data)) deallocate(Data)')
   w.write(f'\n')
   w.write(f'\n      ! Read value to see if it was allocated, return if not')
   w.write(f'\n      call RegUnpack(RF, IsAllocated)')
   w.write(f'\n      if (RegCheckErr(RF, "{name}")) return')
   w.write(f'\n      if (.not. IsAllocated) return')
   w.write(f'\n')
   alloc_dims = ''
   if rank > 0:
      w.write(f'\n      ! Read array bounds')
      w.write(f'\n      call RegUnpackBounds(RF, {rank}, LB, UB)')
      w.write(f'\n      if (RegCheckErr(RF, "{name}")) return')
      alloc_dims = '(' + ','.join([f'LB({d+1}):UB({d+1})' for d in range(rank)]) + ')'
   w.write(f'\n')
   w.write(f'\n      ! Allocate data')
   w.write(f'\n      allocate(Data{alloc_dims}, stat=stat)')
   w.write(f'\n      if (stat /= 0) then')
   w.write(f'\n         RF%ErrStat = ErrID_Fatal')
   w.write(f'\n         RF%ErrMsg = "{name}: error allocating data"')
   w.write(f'\n         return')
   w.write(f'\n      end if')
   w.write(f'\n')
   w.write(f'\n      ! Read data')
   w.write(f'\n      call RegUnpack(RF, Data)')
   w.write(f'\n      if (RegCheckErr(RF, "{name}")) return')
   w.write(f'\n   end subroutine')


def gen_pack_ptr(w, dt, decl, rank):
   dims = '' if rank == 0 else '('+','.join([':']*rank)+')'
   name = f'PackPtr_{dt}'
   if rank > 0: name += f'_Rank{rank}'
   w.write(f'\n')
   w.write(f'\n   subroutine {name}(RF, Data)')
   w.write(f'\n      type(RegFile), intent(inout)         :: RF')
   w.write(f'\n      {decl+", pointer, intent(in)":<35s}  :: Data{dims}')
   w.write(f'\n      logical                              :: PtrInIndex')
   w.write(f'\n')
   w.write(f'\n      ! If error, return')
   w.write(f'\n      if (RF%ErrStat /= ErrID_None) return')
   w.write(f'\n')
   w.write(f'\n      ! Write if associated')
   w.write(f'\n      call RegPack(RF, associated(Data))')
   w.write(f'\n      if (RegCheckErr(RF, "{name}")) return')
   w.write(f'\n      if (.not. associated(Data)) return')
   if rank > 0:
      w.write(f'\n')
      w.write(f'\n      ! Write array bounds')
      w.write(f'\n      call RegPackBounds(RF, {rank}, lbound(Data), ubound(Data))')
   w.write(f'\n')
   w.write(f'\n      ! Write pointer info')
   w.write(f'\n      call RegPackPointer(RF, c_loc(Data), PtrInIndex)')
   w.write(f'\n      if (RegCheckErr(RF, "{name}")) return')
   w.write(f'\n      if (PtrInIndex) return')
   w.write(f'\n')
   w.write(f'\n      ! Write data to file')
   w.write(f'\n      call RegPack(RF, Data)')
   w.write(f'\n      if (RegCheckErr(RF, "{name}")) return')
   w.write(f'\n   end subroutine')

def gen_unpack_ptr(w, dt, decl, rank):
   dims = '' if rank == 0 else '('+','.join([':']*rank)+')'
   dt_size = int(dt[-1])
   name = f'UnpackPtr_{dt}' if rank == 0 else f'UnpackPtr_{dt}_Rank{rank}'
   w.write(f'\n')
   if rank == 0:
      w.write(f'\n   subroutine {name}(RF, Data)')
   else:
      w.write(f'\n   subroutine {name}(RF, Data, LB, UB)')
   w.write(f'\n      type(RegFile), intent(inout)          :: RF')
   w.write(f'\n      {decl+", pointer, intent(out)":<36s}  :: Data{dims}')
   if rank > 0:
      w.write(f'\n      integer(B4Ki), intent(out)            :: LB(:), UB(:)')
   w.write(f'\n      integer(IntKi)                        :: stat')
   w.write(f'\n      integer(B8Ki)                         :: PtrIdx')
   w.write(f'\n      logical                               :: IsAssociated')
   w.write(f'\n      type(c_ptr)                           :: Ptr')
   w.write(f'\n')
   w.write(f'\n      ! If error, return')
   w.write(f'\n      if (RF%ErrStat /= ErrID_None) return')
   w.write(f'\n')
   w.write(f'\n      ! If associated, deallocate and nullify')
   w.write(f'\n      if (associated(Data)) then')
   w.write(f'\n         deallocate(Data)')
   w.write(f'\n         nullify(Data)')
   w.write(f'\n      end if')
   w.write(f'\n')
   w.write(f'\n      ! Read value to see if it was associated, return if not')
   w.write(f'\n      call RegUnpack(RF, IsAssociated)')
   w.write(f'\n      if (RegCheckErr(RF, "{name}")) return')
   w.write(f'\n      if (.not. IsAssociated) return')
   if rank > 0:
      w.write(f'\n')
      w.write(f'\n      ! Read array bounds')
      w.write(f'\n      call RegUnpackBounds(RF, {rank}, LB, UB)')
      w.write(f'\n      if (RegCheckErr(RF, "{name}")) return')
   w.write(f'\n')
   w.write(f'\n      ! Unpack pointer inf')
   w.write(f'\n      call RegUnpackPointer(RF, Ptr, PtrIdx)')
   w.write(f'\n      if (RegCheckErr(RF, "{name}")) return')
   w.write(f'\n')
   w.write(f'\n      ! If pointer was in index, associate data with pointer, return')
   w.write(f'\n      if (c_associated(Ptr)) then')
   if rank == 0:
      alloc_dims = ''
      w.write(f'\n         call c_f_pointer(Ptr, Data)')
   else:
      alloc_dims = '(' + ','.join([f'LB({d+1}):UB({d+1})' for d in range(rank)]) + ')'
      remap_dims = ",".join([f'LB({d+1}):' for d in range(rank)])
      w.write(f'\n         call c_f_pointer(Ptr, Data, UB - LB)')   # Specify shape
      w.write(f'\n         Data({remap_dims}) => Data')        # Remap bounds
   w.write(f'\n         return')
   w.write(f'\n      end if')
   w.write(f'\n')
   w.write(f'\n      ! Allocate data')
   w.write(f'\n      allocate(Data{alloc_dims}, stat=stat)')
   w.write(f'\n      if (stat /= 0) then')
   w.write(f'\n         RF%ErrStat = ErrID_Fatal')
   w.write(f'\n         RF%ErrMsg = "{name}: error allocating data"')
   w.write(f'\n         return')
   w.write(f'\n      end if')
   w.write(f'\n')
   w.write(f'\n      ! Read data')
   w.write(f'\n      call RegUnpack(RF, Data)')
   w.write(f'\n      if (RegCheckErr(RF, "{name}")) return')
   w.write(f'\n   end subroutine')

# Registry interface
groups = ['Pack', 'Unpack', 'PackAlloc', 'UnpackAlloc', 'PackPtr', 'UnpackPtr']
ifc_lines = ''
ranks = [''] + [f'_Rank{r}' for r in range(1,num_ranks+1)]
for attr, punp in product([''], groups):
   ifc_lines += f'\n\n   interface Reg{punp}{attr}'
   funcs = [f'{punp}{attr}_{dt}{rank}'for dt, rank in product(type_map.keys(), ranks)]
   lines = textwrap.wrap('module procedure ' + ', '.join(funcs), 80,
                  initial_indent=" "*6, subsequent_indent=' '*9,break_long_words=False)
   ifc_lines += '\n' + ' &\n'.join(lines)
   ifc_lines += '\n   end interface'

with open('src/ModReg.f90', 'w') as w:
    w.write(module_header.format(ifc_lines=ifc_lines, maxrank=num_ranks))

    # Loop through data types and ranks
    for (dt,decl), rank in product(type_map.items(), range(num_ranks+1)):
        gen_pack(w, dt, decl, rank)
        gen_unpack(w, dt, decl, rank)
        gen_pack_alloc(w, dt, decl, rank)
        gen_unpack_alloc(w, dt, decl, rank)
        gen_pack_ptr(w, dt, decl, rank)
        gen_unpack_ptr(w, dt, decl, rank)

    w.write('\nend module')
