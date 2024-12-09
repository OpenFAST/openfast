
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
   

   interface RegPack
      module procedure Pack_C1, Pack_C1_Rank1, Pack_C1_Rank2, Pack_C1_Rank3, &
         Pack_C1_Rank4, Pack_C1_Rank5, Pack_L1, Pack_L1_Rank1, Pack_L1_Rank2, &
         Pack_L1_Rank3, Pack_L1_Rank4, Pack_L1_Rank5, Pack_I4, Pack_I4_Rank1, &
         Pack_I4_Rank2, Pack_I4_Rank3, Pack_I4_Rank4, Pack_I4_Rank5, Pack_I8, &
         Pack_I8_Rank1, Pack_I8_Rank2, Pack_I8_Rank3, Pack_I8_Rank4, &
         Pack_I8_Rank5, Pack_R4, Pack_R4_Rank1, Pack_R4_Rank2, Pack_R4_Rank3, &
         Pack_R4_Rank4, Pack_R4_Rank5, Pack_R8, Pack_R8_Rank1, Pack_R8_Rank2, &
         Pack_R8_Rank3, Pack_R8_Rank4, Pack_R8_Rank5
   end interface

   interface RegUnpack
      module procedure Unpack_C1, Unpack_C1_Rank1, Unpack_C1_Rank2, &
         Unpack_C1_Rank3, Unpack_C1_Rank4, Unpack_C1_Rank5, Unpack_L1, &
         Unpack_L1_Rank1, Unpack_L1_Rank2, Unpack_L1_Rank3, Unpack_L1_Rank4, &
         Unpack_L1_Rank5, Unpack_I4, Unpack_I4_Rank1, Unpack_I4_Rank2, &
         Unpack_I4_Rank3, Unpack_I4_Rank4, Unpack_I4_Rank5, Unpack_I8, &
         Unpack_I8_Rank1, Unpack_I8_Rank2, Unpack_I8_Rank3, Unpack_I8_Rank4, &
         Unpack_I8_Rank5, Unpack_R4, Unpack_R4_Rank1, Unpack_R4_Rank2, &
         Unpack_R4_Rank3, Unpack_R4_Rank4, Unpack_R4_Rank5, Unpack_R8, &
         Unpack_R8_Rank1, Unpack_R8_Rank2, Unpack_R8_Rank3, Unpack_R8_Rank4, &
         Unpack_R8_Rank5
   end interface

   interface RegPackAlloc
      module procedure PackAlloc_C1, PackAlloc_C1_Rank1, PackAlloc_C1_Rank2, &
         PackAlloc_C1_Rank3, PackAlloc_C1_Rank4, PackAlloc_C1_Rank5, &
         PackAlloc_L1, PackAlloc_L1_Rank1, PackAlloc_L1_Rank2, &
         PackAlloc_L1_Rank3, PackAlloc_L1_Rank4, PackAlloc_L1_Rank5, &
         PackAlloc_I4, PackAlloc_I4_Rank1, PackAlloc_I4_Rank2, &
         PackAlloc_I4_Rank3, PackAlloc_I4_Rank4, PackAlloc_I4_Rank5, &
         PackAlloc_I8, PackAlloc_I8_Rank1, PackAlloc_I8_Rank2, &
         PackAlloc_I8_Rank3, PackAlloc_I8_Rank4, PackAlloc_I8_Rank5, &
         PackAlloc_R4, PackAlloc_R4_Rank1, PackAlloc_R4_Rank2, &
         PackAlloc_R4_Rank3, PackAlloc_R4_Rank4, PackAlloc_R4_Rank5, &
         PackAlloc_R8, PackAlloc_R8_Rank1, PackAlloc_R8_Rank2, &
         PackAlloc_R8_Rank3, PackAlloc_R8_Rank4, PackAlloc_R8_Rank5
   end interface

   interface RegUnpackAlloc
      module procedure UnpackAlloc_C1, UnpackAlloc_C1_Rank1, &
         UnpackAlloc_C1_Rank2, UnpackAlloc_C1_Rank3, UnpackAlloc_C1_Rank4, &
         UnpackAlloc_C1_Rank5, UnpackAlloc_L1, UnpackAlloc_L1_Rank1, &
         UnpackAlloc_L1_Rank2, UnpackAlloc_L1_Rank3, UnpackAlloc_L1_Rank4, &
         UnpackAlloc_L1_Rank5, UnpackAlloc_I4, UnpackAlloc_I4_Rank1, &
         UnpackAlloc_I4_Rank2, UnpackAlloc_I4_Rank3, UnpackAlloc_I4_Rank4, &
         UnpackAlloc_I4_Rank5, UnpackAlloc_I8, UnpackAlloc_I8_Rank1, &
         UnpackAlloc_I8_Rank2, UnpackAlloc_I8_Rank3, UnpackAlloc_I8_Rank4, &
         UnpackAlloc_I8_Rank5, UnpackAlloc_R4, UnpackAlloc_R4_Rank1, &
         UnpackAlloc_R4_Rank2, UnpackAlloc_R4_Rank3, UnpackAlloc_R4_Rank4, &
         UnpackAlloc_R4_Rank5, UnpackAlloc_R8, UnpackAlloc_R8_Rank1, &
         UnpackAlloc_R8_Rank2, UnpackAlloc_R8_Rank3, UnpackAlloc_R8_Rank4, &
         UnpackAlloc_R8_Rank5
   end interface

   interface RegPackPtr
      module procedure PackPtr_C1, PackPtr_C1_Rank1, PackPtr_C1_Rank2, &
         PackPtr_C1_Rank3, PackPtr_C1_Rank4, PackPtr_C1_Rank5, PackPtr_L1, &
         PackPtr_L1_Rank1, PackPtr_L1_Rank2, PackPtr_L1_Rank3, PackPtr_L1_Rank4, &
         PackPtr_L1_Rank5, PackPtr_I4, PackPtr_I4_Rank1, PackPtr_I4_Rank2, &
         PackPtr_I4_Rank3, PackPtr_I4_Rank4, PackPtr_I4_Rank5, PackPtr_I8, &
         PackPtr_I8_Rank1, PackPtr_I8_Rank2, PackPtr_I8_Rank3, PackPtr_I8_Rank4, &
         PackPtr_I8_Rank5, PackPtr_R4, PackPtr_R4_Rank1, PackPtr_R4_Rank2, &
         PackPtr_R4_Rank3, PackPtr_R4_Rank4, PackPtr_R4_Rank5, PackPtr_R8, &
         PackPtr_R8_Rank1, PackPtr_R8_Rank2, PackPtr_R8_Rank3, PackPtr_R8_Rank4, &
         PackPtr_R8_Rank5
   end interface

   interface RegUnpackPtr
      module procedure UnpackPtr_C1, UnpackPtr_C1_Rank1, UnpackPtr_C1_Rank2, &
         UnpackPtr_C1_Rank3, UnpackPtr_C1_Rank4, UnpackPtr_C1_Rank5, &
         UnpackPtr_L1, UnpackPtr_L1_Rank1, UnpackPtr_L1_Rank2, &
         UnpackPtr_L1_Rank3, UnpackPtr_L1_Rank4, UnpackPtr_L1_Rank5, &
         UnpackPtr_I4, UnpackPtr_I4_Rank1, UnpackPtr_I4_Rank2, &
         UnpackPtr_I4_Rank3, UnpackPtr_I4_Rank4, UnpackPtr_I4_Rank5, &
         UnpackPtr_I8, UnpackPtr_I8_Rank1, UnpackPtr_I8_Rank2, &
         UnpackPtr_I8_Rank3, UnpackPtr_I8_Rank4, UnpackPtr_I8_Rank5, &
         UnpackPtr_R4, UnpackPtr_R4_Rank1, UnpackPtr_R4_Rank2, &
         UnpackPtr_R4_Rank3, UnpackPtr_R4_Rank4, UnpackPtr_R4_Rank5, &
         UnpackPtr_R8, UnpackPtr_R8_Rank1, UnpackPtr_R8_Rank2, &
         UnpackPtr_R8_Rank3, UnpackPtr_R8_Rank4, UnpackPtr_R8_Rank5
   end interface

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


   subroutine Pack_C1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      character(*), intent(in)             :: Data

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_C1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      character(*), intent(out)            :: Data

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_C1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      character(*), allocatable, intent(in)  :: Data

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_C1")) return
      if (.not. allocated(Data)) return


      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_C1")) return
   end subroutine

   subroutine UnpackAlloc_C1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      character(*), allocatable, intent(out)  :: Data
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_C1")) return
      if (.not. IsAllocated) return


      ! Allocate data
      allocate(Data, stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_C1: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_C1")) return
   end subroutine

   subroutine PackPtr_C1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      character(*), pointer, intent(in)    :: Data
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_C1")) return
      if (.not. associated(Data)) return

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_C1")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_C1")) return
   end subroutine

   subroutine UnpackPtr_C1(RF, Data)
      type(RegFile), intent(inout)          :: RF
      character(*), pointer, intent(out)    :: Data
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_C1")) return
      if (.not. IsAssociated) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_C1")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data)
         return
      end if

      ! Allocate data
      allocate(Data, stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_C1: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_C1")) return
   end subroutine

   subroutine Pack_C1_Rank1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      character(*), intent(in)             :: Data(:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_C1_Rank1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      character(*), intent(out)            :: Data(:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_C1_Rank1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      character(*), allocatable, intent(in)  :: Data(:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_C1_Rank1")) return
      if (.not. allocated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 1, lbound(Data), ubound(Data))

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_C1_Rank1")) return
   end subroutine

   subroutine UnpackAlloc_C1_Rank1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      character(*), allocatable, intent(out)  :: Data(:)
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated
      integer(B4Ki)                        :: LB(1), UB(1)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_C1_Rank1")) return
      if (.not. IsAllocated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 1, LB, UB)
      if (RegCheckErr(RF, "UnpackAlloc_C1_Rank1")) return

      ! Allocate data
      allocate(Data(LB(1):UB(1)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_C1_Rank1: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_C1_Rank1")) return
   end subroutine

   subroutine PackPtr_C1_Rank1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      character(*), pointer, intent(in)    :: Data(:)
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_C1_Rank1")) return
      if (.not. associated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 1, lbound(Data), ubound(Data))

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_C1_Rank1")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_C1_Rank1")) return
   end subroutine

   subroutine UnpackPtr_C1_Rank1(RF, Data, LB, UB)
      type(RegFile), intent(inout)          :: RF
      character(*), pointer, intent(out)    :: Data(:)
      integer(B4Ki), intent(out)            :: LB(:), UB(:)
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_C1_Rank1")) return
      if (.not. IsAssociated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 1, LB, UB)
      if (RegCheckErr(RF, "UnpackPtr_C1_Rank1")) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_C1_Rank1")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data, UB - LB)
         Data(LB(1):) => Data
         return
      end if

      ! Allocate data
      allocate(Data(LB(1):UB(1)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_C1_Rank1: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_C1_Rank1")) return
   end subroutine

   subroutine Pack_C1_Rank2(RF, Data)
      type(RegFile), intent(inout)         :: RF
      character(*), intent(in)             :: Data(:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_C1_Rank2(RF, Data)
      type(RegFile), intent(inout)         :: RF
      character(*), intent(out)            :: Data(:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_C1_Rank2(RF, Data)
      type(RegFile), intent(inout)         :: RF
      character(*), allocatable, intent(in)  :: Data(:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_C1_Rank2")) return
      if (.not. allocated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 2, lbound(Data), ubound(Data))

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_C1_Rank2")) return
   end subroutine

   subroutine UnpackAlloc_C1_Rank2(RF, Data)
      type(RegFile), intent(inout)         :: RF
      character(*), allocatable, intent(out)  :: Data(:,:)
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated
      integer(B4Ki)                        :: LB(2), UB(2)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_C1_Rank2")) return
      if (.not. IsAllocated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 2, LB, UB)
      if (RegCheckErr(RF, "UnpackAlloc_C1_Rank2")) return

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_C1_Rank2: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_C1_Rank2")) return
   end subroutine

   subroutine PackPtr_C1_Rank2(RF, Data)
      type(RegFile), intent(inout)         :: RF
      character(*), pointer, intent(in)    :: Data(:,:)
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_C1_Rank2")) return
      if (.not. associated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 2, lbound(Data), ubound(Data))

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_C1_Rank2")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_C1_Rank2")) return
   end subroutine

   subroutine UnpackPtr_C1_Rank2(RF, Data, LB, UB)
      type(RegFile), intent(inout)          :: RF
      character(*), pointer, intent(out)    :: Data(:,:)
      integer(B4Ki), intent(out)            :: LB(:), UB(:)
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_C1_Rank2")) return
      if (.not. IsAssociated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 2, LB, UB)
      if (RegCheckErr(RF, "UnpackPtr_C1_Rank2")) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_C1_Rank2")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data, UB - LB)
         Data(LB(1):,LB(2):) => Data
         return
      end if

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_C1_Rank2: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_C1_Rank2")) return
   end subroutine

   subroutine Pack_C1_Rank3(RF, Data)
      type(RegFile), intent(inout)         :: RF
      character(*), intent(in)             :: Data(:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_C1_Rank3(RF, Data)
      type(RegFile), intent(inout)         :: RF
      character(*), intent(out)            :: Data(:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_C1_Rank3(RF, Data)
      type(RegFile), intent(inout)         :: RF
      character(*), allocatable, intent(in)  :: Data(:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_C1_Rank3")) return
      if (.not. allocated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 3, lbound(Data), ubound(Data))

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_C1_Rank3")) return
   end subroutine

   subroutine UnpackAlloc_C1_Rank3(RF, Data)
      type(RegFile), intent(inout)         :: RF
      character(*), allocatable, intent(out)  :: Data(:,:,:)
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated
      integer(B4Ki)                        :: LB(3), UB(3)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_C1_Rank3")) return
      if (.not. IsAllocated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 3, LB, UB)
      if (RegCheckErr(RF, "UnpackAlloc_C1_Rank3")) return

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_C1_Rank3: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_C1_Rank3")) return
   end subroutine

   subroutine PackPtr_C1_Rank3(RF, Data)
      type(RegFile), intent(inout)         :: RF
      character(*), pointer, intent(in)    :: Data(:,:,:)
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_C1_Rank3")) return
      if (.not. associated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 3, lbound(Data), ubound(Data))

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_C1_Rank3")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_C1_Rank3")) return
   end subroutine

   subroutine UnpackPtr_C1_Rank3(RF, Data, LB, UB)
      type(RegFile), intent(inout)          :: RF
      character(*), pointer, intent(out)    :: Data(:,:,:)
      integer(B4Ki), intent(out)            :: LB(:), UB(:)
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_C1_Rank3")) return
      if (.not. IsAssociated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 3, LB, UB)
      if (RegCheckErr(RF, "UnpackPtr_C1_Rank3")) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_C1_Rank3")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data, UB - LB)
         Data(LB(1):,LB(2):,LB(3):) => Data
         return
      end if

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_C1_Rank3: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_C1_Rank3")) return
   end subroutine

   subroutine Pack_C1_Rank4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      character(*), intent(in)             :: Data(:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_C1_Rank4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      character(*), intent(out)            :: Data(:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_C1_Rank4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      character(*), allocatable, intent(in)  :: Data(:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_C1_Rank4")) return
      if (.not. allocated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 4, lbound(Data), ubound(Data))

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_C1_Rank4")) return
   end subroutine

   subroutine UnpackAlloc_C1_Rank4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      character(*), allocatable, intent(out)  :: Data(:,:,:,:)
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated
      integer(B4Ki)                        :: LB(4), UB(4)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_C1_Rank4")) return
      if (.not. IsAllocated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 4, LB, UB)
      if (RegCheckErr(RF, "UnpackAlloc_C1_Rank4")) return

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3),LB(4):UB(4)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_C1_Rank4: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_C1_Rank4")) return
   end subroutine

   subroutine PackPtr_C1_Rank4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      character(*), pointer, intent(in)    :: Data(:,:,:,:)
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_C1_Rank4")) return
      if (.not. associated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 4, lbound(Data), ubound(Data))

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_C1_Rank4")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_C1_Rank4")) return
   end subroutine

   subroutine UnpackPtr_C1_Rank4(RF, Data, LB, UB)
      type(RegFile), intent(inout)          :: RF
      character(*), pointer, intent(out)    :: Data(:,:,:,:)
      integer(B4Ki), intent(out)            :: LB(:), UB(:)
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_C1_Rank4")) return
      if (.not. IsAssociated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 4, LB, UB)
      if (RegCheckErr(RF, "UnpackPtr_C1_Rank4")) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_C1_Rank4")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data, UB - LB)
         Data(LB(1):,LB(2):,LB(3):,LB(4):) => Data
         return
      end if

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3),LB(4):UB(4)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_C1_Rank4: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_C1_Rank4")) return
   end subroutine

   subroutine Pack_C1_Rank5(RF, Data)
      type(RegFile), intent(inout)         :: RF
      character(*), intent(in)             :: Data(:,:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_C1_Rank5(RF, Data)
      type(RegFile), intent(inout)         :: RF
      character(*), intent(out)            :: Data(:,:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_C1_Rank5(RF, Data)
      type(RegFile), intent(inout)         :: RF
      character(*), allocatable, intent(in)  :: Data(:,:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_C1_Rank5")) return
      if (.not. allocated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 5, lbound(Data), ubound(Data))

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_C1_Rank5")) return
   end subroutine

   subroutine UnpackAlloc_C1_Rank5(RF, Data)
      type(RegFile), intent(inout)         :: RF
      character(*), allocatable, intent(out)  :: Data(:,:,:,:,:)
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated
      integer(B4Ki)                        :: LB(5), UB(5)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_C1_Rank5")) return
      if (.not. IsAllocated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 5, LB, UB)
      if (RegCheckErr(RF, "UnpackAlloc_C1_Rank5")) return

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3),LB(4):UB(4),LB(5):UB(5)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_C1_Rank5: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_C1_Rank5")) return
   end subroutine

   subroutine PackPtr_C1_Rank5(RF, Data)
      type(RegFile), intent(inout)         :: RF
      character(*), pointer, intent(in)    :: Data(:,:,:,:,:)
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_C1_Rank5")) return
      if (.not. associated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 5, lbound(Data), ubound(Data))

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_C1_Rank5")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_C1_Rank5")) return
   end subroutine

   subroutine UnpackPtr_C1_Rank5(RF, Data, LB, UB)
      type(RegFile), intent(inout)          :: RF
      character(*), pointer, intent(out)    :: Data(:,:,:,:,:)
      integer(B4Ki), intent(out)            :: LB(:), UB(:)
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_C1_Rank5")) return
      if (.not. IsAssociated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 5, LB, UB)
      if (RegCheckErr(RF, "UnpackPtr_C1_Rank5")) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_C1_Rank5")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data, UB - LB)
         Data(LB(1):,LB(2):,LB(3):,LB(4):,LB(5):) => Data
         return
      end if

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3),LB(4):UB(4),LB(5):UB(5)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_C1_Rank5: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_C1_Rank5")) return
   end subroutine

   subroutine Pack_L1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      logical, intent(in)                  :: Data

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_L1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      logical, intent(out)                 :: Data

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_L1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      logical, allocatable, intent(in)     :: Data

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_L1")) return
      if (.not. allocated(Data)) return


      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_L1")) return
   end subroutine

   subroutine UnpackAlloc_L1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      logical, allocatable, intent(out)    :: Data
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_L1")) return
      if (.not. IsAllocated) return


      ! Allocate data
      allocate(Data, stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_L1: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_L1")) return
   end subroutine

   subroutine PackPtr_L1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      logical, pointer, intent(in)         :: Data
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_L1")) return
      if (.not. associated(Data)) return

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_L1")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_L1")) return
   end subroutine

   subroutine UnpackPtr_L1(RF, Data)
      type(RegFile), intent(inout)          :: RF
      logical, pointer, intent(out)         :: Data
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_L1")) return
      if (.not. IsAssociated) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_L1")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data)
         return
      end if

      ! Allocate data
      allocate(Data, stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_L1: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_L1")) return
   end subroutine

   subroutine Pack_L1_Rank1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      logical, intent(in)                  :: Data(:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_L1_Rank1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      logical, intent(out)                 :: Data(:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_L1_Rank1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      logical, allocatable, intent(in)     :: Data(:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_L1_Rank1")) return
      if (.not. allocated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 1, lbound(Data), ubound(Data))

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_L1_Rank1")) return
   end subroutine

   subroutine UnpackAlloc_L1_Rank1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      logical, allocatable, intent(out)    :: Data(:)
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated
      integer(B4Ki)                        :: LB(1), UB(1)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_L1_Rank1")) return
      if (.not. IsAllocated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 1, LB, UB)
      if (RegCheckErr(RF, "UnpackAlloc_L1_Rank1")) return

      ! Allocate data
      allocate(Data(LB(1):UB(1)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_L1_Rank1: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_L1_Rank1")) return
   end subroutine

   subroutine PackPtr_L1_Rank1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      logical, pointer, intent(in)         :: Data(:)
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_L1_Rank1")) return
      if (.not. associated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 1, lbound(Data), ubound(Data))

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_L1_Rank1")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_L1_Rank1")) return
   end subroutine

   subroutine UnpackPtr_L1_Rank1(RF, Data, LB, UB)
      type(RegFile), intent(inout)          :: RF
      logical, pointer, intent(out)         :: Data(:)
      integer(B4Ki), intent(out)            :: LB(:), UB(:)
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_L1_Rank1")) return
      if (.not. IsAssociated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 1, LB, UB)
      if (RegCheckErr(RF, "UnpackPtr_L1_Rank1")) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_L1_Rank1")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data, UB - LB)
         Data(LB(1):) => Data
         return
      end if

      ! Allocate data
      allocate(Data(LB(1):UB(1)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_L1_Rank1: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_L1_Rank1")) return
   end subroutine

   subroutine Pack_L1_Rank2(RF, Data)
      type(RegFile), intent(inout)         :: RF
      logical, intent(in)                  :: Data(:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_L1_Rank2(RF, Data)
      type(RegFile), intent(inout)         :: RF
      logical, intent(out)                 :: Data(:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_L1_Rank2(RF, Data)
      type(RegFile), intent(inout)         :: RF
      logical, allocatable, intent(in)     :: Data(:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_L1_Rank2")) return
      if (.not. allocated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 2, lbound(Data), ubound(Data))

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_L1_Rank2")) return
   end subroutine

   subroutine UnpackAlloc_L1_Rank2(RF, Data)
      type(RegFile), intent(inout)         :: RF
      logical, allocatable, intent(out)    :: Data(:,:)
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated
      integer(B4Ki)                        :: LB(2), UB(2)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_L1_Rank2")) return
      if (.not. IsAllocated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 2, LB, UB)
      if (RegCheckErr(RF, "UnpackAlloc_L1_Rank2")) return

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_L1_Rank2: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_L1_Rank2")) return
   end subroutine

   subroutine PackPtr_L1_Rank2(RF, Data)
      type(RegFile), intent(inout)         :: RF
      logical, pointer, intent(in)         :: Data(:,:)
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_L1_Rank2")) return
      if (.not. associated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 2, lbound(Data), ubound(Data))

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_L1_Rank2")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_L1_Rank2")) return
   end subroutine

   subroutine UnpackPtr_L1_Rank2(RF, Data, LB, UB)
      type(RegFile), intent(inout)          :: RF
      logical, pointer, intent(out)         :: Data(:,:)
      integer(B4Ki), intent(out)            :: LB(:), UB(:)
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_L1_Rank2")) return
      if (.not. IsAssociated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 2, LB, UB)
      if (RegCheckErr(RF, "UnpackPtr_L1_Rank2")) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_L1_Rank2")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data, UB - LB)
         Data(LB(1):,LB(2):) => Data
         return
      end if

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_L1_Rank2: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_L1_Rank2")) return
   end subroutine

   subroutine Pack_L1_Rank3(RF, Data)
      type(RegFile), intent(inout)         :: RF
      logical, intent(in)                  :: Data(:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_L1_Rank3(RF, Data)
      type(RegFile), intent(inout)         :: RF
      logical, intent(out)                 :: Data(:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_L1_Rank3(RF, Data)
      type(RegFile), intent(inout)         :: RF
      logical, allocatable, intent(in)     :: Data(:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_L1_Rank3")) return
      if (.not. allocated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 3, lbound(Data), ubound(Data))

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_L1_Rank3")) return
   end subroutine

   subroutine UnpackAlloc_L1_Rank3(RF, Data)
      type(RegFile), intent(inout)         :: RF
      logical, allocatable, intent(out)    :: Data(:,:,:)
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated
      integer(B4Ki)                        :: LB(3), UB(3)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_L1_Rank3")) return
      if (.not. IsAllocated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 3, LB, UB)
      if (RegCheckErr(RF, "UnpackAlloc_L1_Rank3")) return

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_L1_Rank3: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_L1_Rank3")) return
   end subroutine

   subroutine PackPtr_L1_Rank3(RF, Data)
      type(RegFile), intent(inout)         :: RF
      logical, pointer, intent(in)         :: Data(:,:,:)
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_L1_Rank3")) return
      if (.not. associated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 3, lbound(Data), ubound(Data))

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_L1_Rank3")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_L1_Rank3")) return
   end subroutine

   subroutine UnpackPtr_L1_Rank3(RF, Data, LB, UB)
      type(RegFile), intent(inout)          :: RF
      logical, pointer, intent(out)         :: Data(:,:,:)
      integer(B4Ki), intent(out)            :: LB(:), UB(:)
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_L1_Rank3")) return
      if (.not. IsAssociated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 3, LB, UB)
      if (RegCheckErr(RF, "UnpackPtr_L1_Rank3")) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_L1_Rank3")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data, UB - LB)
         Data(LB(1):,LB(2):,LB(3):) => Data
         return
      end if

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_L1_Rank3: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_L1_Rank3")) return
   end subroutine

   subroutine Pack_L1_Rank4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      logical, intent(in)                  :: Data(:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_L1_Rank4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      logical, intent(out)                 :: Data(:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_L1_Rank4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      logical, allocatable, intent(in)     :: Data(:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_L1_Rank4")) return
      if (.not. allocated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 4, lbound(Data), ubound(Data))

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_L1_Rank4")) return
   end subroutine

   subroutine UnpackAlloc_L1_Rank4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      logical, allocatable, intent(out)    :: Data(:,:,:,:)
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated
      integer(B4Ki)                        :: LB(4), UB(4)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_L1_Rank4")) return
      if (.not. IsAllocated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 4, LB, UB)
      if (RegCheckErr(RF, "UnpackAlloc_L1_Rank4")) return

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3),LB(4):UB(4)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_L1_Rank4: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_L1_Rank4")) return
   end subroutine

   subroutine PackPtr_L1_Rank4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      logical, pointer, intent(in)         :: Data(:,:,:,:)
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_L1_Rank4")) return
      if (.not. associated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 4, lbound(Data), ubound(Data))

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_L1_Rank4")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_L1_Rank4")) return
   end subroutine

   subroutine UnpackPtr_L1_Rank4(RF, Data, LB, UB)
      type(RegFile), intent(inout)          :: RF
      logical, pointer, intent(out)         :: Data(:,:,:,:)
      integer(B4Ki), intent(out)            :: LB(:), UB(:)
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_L1_Rank4")) return
      if (.not. IsAssociated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 4, LB, UB)
      if (RegCheckErr(RF, "UnpackPtr_L1_Rank4")) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_L1_Rank4")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data, UB - LB)
         Data(LB(1):,LB(2):,LB(3):,LB(4):) => Data
         return
      end if

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3),LB(4):UB(4)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_L1_Rank4: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_L1_Rank4")) return
   end subroutine

   subroutine Pack_L1_Rank5(RF, Data)
      type(RegFile), intent(inout)         :: RF
      logical, intent(in)                  :: Data(:,:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_L1_Rank5(RF, Data)
      type(RegFile), intent(inout)         :: RF
      logical, intent(out)                 :: Data(:,:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_L1_Rank5(RF, Data)
      type(RegFile), intent(inout)         :: RF
      logical, allocatable, intent(in)     :: Data(:,:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_L1_Rank5")) return
      if (.not. allocated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 5, lbound(Data), ubound(Data))

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_L1_Rank5")) return
   end subroutine

   subroutine UnpackAlloc_L1_Rank5(RF, Data)
      type(RegFile), intent(inout)         :: RF
      logical, allocatable, intent(out)    :: Data(:,:,:,:,:)
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated
      integer(B4Ki)                        :: LB(5), UB(5)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_L1_Rank5")) return
      if (.not. IsAllocated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 5, LB, UB)
      if (RegCheckErr(RF, "UnpackAlloc_L1_Rank5")) return

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3),LB(4):UB(4),LB(5):UB(5)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_L1_Rank5: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_L1_Rank5")) return
   end subroutine

   subroutine PackPtr_L1_Rank5(RF, Data)
      type(RegFile), intent(inout)         :: RF
      logical, pointer, intent(in)         :: Data(:,:,:,:,:)
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_L1_Rank5")) return
      if (.not. associated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 5, lbound(Data), ubound(Data))

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_L1_Rank5")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_L1_Rank5")) return
   end subroutine

   subroutine UnpackPtr_L1_Rank5(RF, Data, LB, UB)
      type(RegFile), intent(inout)          :: RF
      logical, pointer, intent(out)         :: Data(:,:,:,:,:)
      integer(B4Ki), intent(out)            :: LB(:), UB(:)
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_L1_Rank5")) return
      if (.not. IsAssociated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 5, LB, UB)
      if (RegCheckErr(RF, "UnpackPtr_L1_Rank5")) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_L1_Rank5")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data, UB - LB)
         Data(LB(1):,LB(2):,LB(3):,LB(4):,LB(5):) => Data
         return
      end if

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3),LB(4):UB(4),LB(5):UB(5)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_L1_Rank5: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_L1_Rank5")) return
   end subroutine

   subroutine Pack_I4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B4Ki), intent(in)            :: Data

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_I4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B4Ki), intent(out)           :: Data

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_I4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B4Ki), allocatable, intent(in)  :: Data

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_I4")) return
      if (.not. allocated(Data)) return


      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_I4")) return
   end subroutine

   subroutine UnpackAlloc_I4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B4Ki), allocatable, intent(out)  :: Data
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_I4")) return
      if (.not. IsAllocated) return


      ! Allocate data
      allocate(Data, stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_I4: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_I4")) return
   end subroutine

   subroutine PackPtr_I4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B4Ki), pointer, intent(in)   :: Data
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_I4")) return
      if (.not. associated(Data)) return

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_I4")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_I4")) return
   end subroutine

   subroutine UnpackPtr_I4(RF, Data)
      type(RegFile), intent(inout)          :: RF
      integer(B4Ki), pointer, intent(out)   :: Data
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_I4")) return
      if (.not. IsAssociated) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_I4")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data)
         return
      end if

      ! Allocate data
      allocate(Data, stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_I4: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_I4")) return
   end subroutine

   subroutine Pack_I4_Rank1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B4Ki), intent(in)            :: Data(:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_I4_Rank1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B4Ki), intent(out)           :: Data(:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_I4_Rank1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B4Ki), allocatable, intent(in)  :: Data(:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_I4_Rank1")) return
      if (.not. allocated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 1, lbound(Data), ubound(Data))

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_I4_Rank1")) return
   end subroutine

   subroutine UnpackAlloc_I4_Rank1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B4Ki), allocatable, intent(out)  :: Data(:)
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated
      integer(B4Ki)                        :: LB(1), UB(1)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_I4_Rank1")) return
      if (.not. IsAllocated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 1, LB, UB)
      if (RegCheckErr(RF, "UnpackAlloc_I4_Rank1")) return

      ! Allocate data
      allocate(Data(LB(1):UB(1)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_I4_Rank1: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_I4_Rank1")) return
   end subroutine

   subroutine PackPtr_I4_Rank1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B4Ki), pointer, intent(in)   :: Data(:)
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_I4_Rank1")) return
      if (.not. associated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 1, lbound(Data), ubound(Data))

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_I4_Rank1")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_I4_Rank1")) return
   end subroutine

   subroutine UnpackPtr_I4_Rank1(RF, Data, LB, UB)
      type(RegFile), intent(inout)          :: RF
      integer(B4Ki), pointer, intent(out)   :: Data(:)
      integer(B4Ki), intent(out)            :: LB(:), UB(:)
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_I4_Rank1")) return
      if (.not. IsAssociated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 1, LB, UB)
      if (RegCheckErr(RF, "UnpackPtr_I4_Rank1")) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_I4_Rank1")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data, UB - LB)
         Data(LB(1):) => Data
         return
      end if

      ! Allocate data
      allocate(Data(LB(1):UB(1)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_I4_Rank1: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_I4_Rank1")) return
   end subroutine

   subroutine Pack_I4_Rank2(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B4Ki), intent(in)            :: Data(:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_I4_Rank2(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B4Ki), intent(out)           :: Data(:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_I4_Rank2(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B4Ki), allocatable, intent(in)  :: Data(:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_I4_Rank2")) return
      if (.not. allocated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 2, lbound(Data), ubound(Data))

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_I4_Rank2")) return
   end subroutine

   subroutine UnpackAlloc_I4_Rank2(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B4Ki), allocatable, intent(out)  :: Data(:,:)
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated
      integer(B4Ki)                        :: LB(2), UB(2)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_I4_Rank2")) return
      if (.not. IsAllocated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 2, LB, UB)
      if (RegCheckErr(RF, "UnpackAlloc_I4_Rank2")) return

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_I4_Rank2: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_I4_Rank2")) return
   end subroutine

   subroutine PackPtr_I4_Rank2(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B4Ki), pointer, intent(in)   :: Data(:,:)
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_I4_Rank2")) return
      if (.not. associated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 2, lbound(Data), ubound(Data))

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_I4_Rank2")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_I4_Rank2")) return
   end subroutine

   subroutine UnpackPtr_I4_Rank2(RF, Data, LB, UB)
      type(RegFile), intent(inout)          :: RF
      integer(B4Ki), pointer, intent(out)   :: Data(:,:)
      integer(B4Ki), intent(out)            :: LB(:), UB(:)
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_I4_Rank2")) return
      if (.not. IsAssociated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 2, LB, UB)
      if (RegCheckErr(RF, "UnpackPtr_I4_Rank2")) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_I4_Rank2")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data, UB - LB)
         Data(LB(1):,LB(2):) => Data
         return
      end if

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_I4_Rank2: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_I4_Rank2")) return
   end subroutine

   subroutine Pack_I4_Rank3(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B4Ki), intent(in)            :: Data(:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_I4_Rank3(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B4Ki), intent(out)           :: Data(:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_I4_Rank3(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B4Ki), allocatable, intent(in)  :: Data(:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_I4_Rank3")) return
      if (.not. allocated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 3, lbound(Data), ubound(Data))

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_I4_Rank3")) return
   end subroutine

   subroutine UnpackAlloc_I4_Rank3(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B4Ki), allocatable, intent(out)  :: Data(:,:,:)
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated
      integer(B4Ki)                        :: LB(3), UB(3)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_I4_Rank3")) return
      if (.not. IsAllocated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 3, LB, UB)
      if (RegCheckErr(RF, "UnpackAlloc_I4_Rank3")) return

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_I4_Rank3: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_I4_Rank3")) return
   end subroutine

   subroutine PackPtr_I4_Rank3(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B4Ki), pointer, intent(in)   :: Data(:,:,:)
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_I4_Rank3")) return
      if (.not. associated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 3, lbound(Data), ubound(Data))

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_I4_Rank3")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_I4_Rank3")) return
   end subroutine

   subroutine UnpackPtr_I4_Rank3(RF, Data, LB, UB)
      type(RegFile), intent(inout)          :: RF
      integer(B4Ki), pointer, intent(out)   :: Data(:,:,:)
      integer(B4Ki), intent(out)            :: LB(:), UB(:)
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_I4_Rank3")) return
      if (.not. IsAssociated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 3, LB, UB)
      if (RegCheckErr(RF, "UnpackPtr_I4_Rank3")) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_I4_Rank3")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data, UB - LB)
         Data(LB(1):,LB(2):,LB(3):) => Data
         return
      end if

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_I4_Rank3: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_I4_Rank3")) return
   end subroutine

   subroutine Pack_I4_Rank4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B4Ki), intent(in)            :: Data(:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_I4_Rank4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B4Ki), intent(out)           :: Data(:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_I4_Rank4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B4Ki), allocatable, intent(in)  :: Data(:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_I4_Rank4")) return
      if (.not. allocated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 4, lbound(Data), ubound(Data))

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_I4_Rank4")) return
   end subroutine

   subroutine UnpackAlloc_I4_Rank4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B4Ki), allocatable, intent(out)  :: Data(:,:,:,:)
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated
      integer(B4Ki)                        :: LB(4), UB(4)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_I4_Rank4")) return
      if (.not. IsAllocated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 4, LB, UB)
      if (RegCheckErr(RF, "UnpackAlloc_I4_Rank4")) return

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3),LB(4):UB(4)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_I4_Rank4: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_I4_Rank4")) return
   end subroutine

   subroutine PackPtr_I4_Rank4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B4Ki), pointer, intent(in)   :: Data(:,:,:,:)
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_I4_Rank4")) return
      if (.not. associated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 4, lbound(Data), ubound(Data))

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_I4_Rank4")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_I4_Rank4")) return
   end subroutine

   subroutine UnpackPtr_I4_Rank4(RF, Data, LB, UB)
      type(RegFile), intent(inout)          :: RF
      integer(B4Ki), pointer, intent(out)   :: Data(:,:,:,:)
      integer(B4Ki), intent(out)            :: LB(:), UB(:)
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_I4_Rank4")) return
      if (.not. IsAssociated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 4, LB, UB)
      if (RegCheckErr(RF, "UnpackPtr_I4_Rank4")) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_I4_Rank4")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data, UB - LB)
         Data(LB(1):,LB(2):,LB(3):,LB(4):) => Data
         return
      end if

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3),LB(4):UB(4)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_I4_Rank4: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_I4_Rank4")) return
   end subroutine

   subroutine Pack_I4_Rank5(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B4Ki), intent(in)            :: Data(:,:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_I4_Rank5(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B4Ki), intent(out)           :: Data(:,:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_I4_Rank5(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B4Ki), allocatable, intent(in)  :: Data(:,:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_I4_Rank5")) return
      if (.not. allocated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 5, lbound(Data), ubound(Data))

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_I4_Rank5")) return
   end subroutine

   subroutine UnpackAlloc_I4_Rank5(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B4Ki), allocatable, intent(out)  :: Data(:,:,:,:,:)
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated
      integer(B4Ki)                        :: LB(5), UB(5)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_I4_Rank5")) return
      if (.not. IsAllocated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 5, LB, UB)
      if (RegCheckErr(RF, "UnpackAlloc_I4_Rank5")) return

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3),LB(4):UB(4),LB(5):UB(5)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_I4_Rank5: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_I4_Rank5")) return
   end subroutine

   subroutine PackPtr_I4_Rank5(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B4Ki), pointer, intent(in)   :: Data(:,:,:,:,:)
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_I4_Rank5")) return
      if (.not. associated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 5, lbound(Data), ubound(Data))

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_I4_Rank5")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_I4_Rank5")) return
   end subroutine

   subroutine UnpackPtr_I4_Rank5(RF, Data, LB, UB)
      type(RegFile), intent(inout)          :: RF
      integer(B4Ki), pointer, intent(out)   :: Data(:,:,:,:,:)
      integer(B4Ki), intent(out)            :: LB(:), UB(:)
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_I4_Rank5")) return
      if (.not. IsAssociated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 5, LB, UB)
      if (RegCheckErr(RF, "UnpackPtr_I4_Rank5")) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_I4_Rank5")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data, UB - LB)
         Data(LB(1):,LB(2):,LB(3):,LB(4):,LB(5):) => Data
         return
      end if

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3),LB(4):UB(4),LB(5):UB(5)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_I4_Rank5: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_I4_Rank5")) return
   end subroutine

   subroutine Pack_I8(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B8Ki), intent(in)            :: Data

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_I8(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B8Ki), intent(out)           :: Data

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_I8(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B8Ki), allocatable, intent(in)  :: Data

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_I8")) return
      if (.not. allocated(Data)) return


      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_I8")) return
   end subroutine

   subroutine UnpackAlloc_I8(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B8Ki), allocatable, intent(out)  :: Data
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_I8")) return
      if (.not. IsAllocated) return


      ! Allocate data
      allocate(Data, stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_I8: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_I8")) return
   end subroutine

   subroutine PackPtr_I8(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B8Ki), pointer, intent(in)   :: Data
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_I8")) return
      if (.not. associated(Data)) return

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_I8")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_I8")) return
   end subroutine

   subroutine UnpackPtr_I8(RF, Data)
      type(RegFile), intent(inout)          :: RF
      integer(B8Ki), pointer, intent(out)   :: Data
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_I8")) return
      if (.not. IsAssociated) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_I8")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data)
         return
      end if

      ! Allocate data
      allocate(Data, stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_I8: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_I8")) return
   end subroutine

   subroutine Pack_I8_Rank1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B8Ki), intent(in)            :: Data(:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_I8_Rank1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B8Ki), intent(out)           :: Data(:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_I8_Rank1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B8Ki), allocatable, intent(in)  :: Data(:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_I8_Rank1")) return
      if (.not. allocated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 1, lbound(Data), ubound(Data))

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_I8_Rank1")) return
   end subroutine

   subroutine UnpackAlloc_I8_Rank1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B8Ki), allocatable, intent(out)  :: Data(:)
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated
      integer(B4Ki)                        :: LB(1), UB(1)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_I8_Rank1")) return
      if (.not. IsAllocated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 1, LB, UB)
      if (RegCheckErr(RF, "UnpackAlloc_I8_Rank1")) return

      ! Allocate data
      allocate(Data(LB(1):UB(1)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_I8_Rank1: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_I8_Rank1")) return
   end subroutine

   subroutine PackPtr_I8_Rank1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B8Ki), pointer, intent(in)   :: Data(:)
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_I8_Rank1")) return
      if (.not. associated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 1, lbound(Data), ubound(Data))

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_I8_Rank1")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_I8_Rank1")) return
   end subroutine

   subroutine UnpackPtr_I8_Rank1(RF, Data, LB, UB)
      type(RegFile), intent(inout)          :: RF
      integer(B8Ki), pointer, intent(out)   :: Data(:)
      integer(B4Ki), intent(out)            :: LB(:), UB(:)
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_I8_Rank1")) return
      if (.not. IsAssociated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 1, LB, UB)
      if (RegCheckErr(RF, "UnpackPtr_I8_Rank1")) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_I8_Rank1")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data, UB - LB)
         Data(LB(1):) => Data
         return
      end if

      ! Allocate data
      allocate(Data(LB(1):UB(1)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_I8_Rank1: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_I8_Rank1")) return
   end subroutine

   subroutine Pack_I8_Rank2(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B8Ki), intent(in)            :: Data(:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_I8_Rank2(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B8Ki), intent(out)           :: Data(:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_I8_Rank2(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B8Ki), allocatable, intent(in)  :: Data(:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_I8_Rank2")) return
      if (.not. allocated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 2, lbound(Data), ubound(Data))

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_I8_Rank2")) return
   end subroutine

   subroutine UnpackAlloc_I8_Rank2(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B8Ki), allocatable, intent(out)  :: Data(:,:)
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated
      integer(B4Ki)                        :: LB(2), UB(2)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_I8_Rank2")) return
      if (.not. IsAllocated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 2, LB, UB)
      if (RegCheckErr(RF, "UnpackAlloc_I8_Rank2")) return

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_I8_Rank2: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_I8_Rank2")) return
   end subroutine

   subroutine PackPtr_I8_Rank2(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B8Ki), pointer, intent(in)   :: Data(:,:)
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_I8_Rank2")) return
      if (.not. associated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 2, lbound(Data), ubound(Data))

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_I8_Rank2")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_I8_Rank2")) return
   end subroutine

   subroutine UnpackPtr_I8_Rank2(RF, Data, LB, UB)
      type(RegFile), intent(inout)          :: RF
      integer(B8Ki), pointer, intent(out)   :: Data(:,:)
      integer(B4Ki), intent(out)            :: LB(:), UB(:)
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_I8_Rank2")) return
      if (.not. IsAssociated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 2, LB, UB)
      if (RegCheckErr(RF, "UnpackPtr_I8_Rank2")) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_I8_Rank2")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data, UB - LB)
         Data(LB(1):,LB(2):) => Data
         return
      end if

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_I8_Rank2: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_I8_Rank2")) return
   end subroutine

   subroutine Pack_I8_Rank3(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B8Ki), intent(in)            :: Data(:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_I8_Rank3(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B8Ki), intent(out)           :: Data(:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_I8_Rank3(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B8Ki), allocatable, intent(in)  :: Data(:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_I8_Rank3")) return
      if (.not. allocated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 3, lbound(Data), ubound(Data))

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_I8_Rank3")) return
   end subroutine

   subroutine UnpackAlloc_I8_Rank3(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B8Ki), allocatable, intent(out)  :: Data(:,:,:)
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated
      integer(B4Ki)                        :: LB(3), UB(3)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_I8_Rank3")) return
      if (.not. IsAllocated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 3, LB, UB)
      if (RegCheckErr(RF, "UnpackAlloc_I8_Rank3")) return

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_I8_Rank3: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_I8_Rank3")) return
   end subroutine

   subroutine PackPtr_I8_Rank3(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B8Ki), pointer, intent(in)   :: Data(:,:,:)
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_I8_Rank3")) return
      if (.not. associated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 3, lbound(Data), ubound(Data))

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_I8_Rank3")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_I8_Rank3")) return
   end subroutine

   subroutine UnpackPtr_I8_Rank3(RF, Data, LB, UB)
      type(RegFile), intent(inout)          :: RF
      integer(B8Ki), pointer, intent(out)   :: Data(:,:,:)
      integer(B4Ki), intent(out)            :: LB(:), UB(:)
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_I8_Rank3")) return
      if (.not. IsAssociated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 3, LB, UB)
      if (RegCheckErr(RF, "UnpackPtr_I8_Rank3")) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_I8_Rank3")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data, UB - LB)
         Data(LB(1):,LB(2):,LB(3):) => Data
         return
      end if

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_I8_Rank3: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_I8_Rank3")) return
   end subroutine

   subroutine Pack_I8_Rank4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B8Ki), intent(in)            :: Data(:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_I8_Rank4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B8Ki), intent(out)           :: Data(:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_I8_Rank4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B8Ki), allocatable, intent(in)  :: Data(:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_I8_Rank4")) return
      if (.not. allocated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 4, lbound(Data), ubound(Data))

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_I8_Rank4")) return
   end subroutine

   subroutine UnpackAlloc_I8_Rank4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B8Ki), allocatable, intent(out)  :: Data(:,:,:,:)
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated
      integer(B4Ki)                        :: LB(4), UB(4)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_I8_Rank4")) return
      if (.not. IsAllocated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 4, LB, UB)
      if (RegCheckErr(RF, "UnpackAlloc_I8_Rank4")) return

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3),LB(4):UB(4)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_I8_Rank4: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_I8_Rank4")) return
   end subroutine

   subroutine PackPtr_I8_Rank4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B8Ki), pointer, intent(in)   :: Data(:,:,:,:)
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_I8_Rank4")) return
      if (.not. associated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 4, lbound(Data), ubound(Data))

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_I8_Rank4")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_I8_Rank4")) return
   end subroutine

   subroutine UnpackPtr_I8_Rank4(RF, Data, LB, UB)
      type(RegFile), intent(inout)          :: RF
      integer(B8Ki), pointer, intent(out)   :: Data(:,:,:,:)
      integer(B4Ki), intent(out)            :: LB(:), UB(:)
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_I8_Rank4")) return
      if (.not. IsAssociated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 4, LB, UB)
      if (RegCheckErr(RF, "UnpackPtr_I8_Rank4")) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_I8_Rank4")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data, UB - LB)
         Data(LB(1):,LB(2):,LB(3):,LB(4):) => Data
         return
      end if

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3),LB(4):UB(4)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_I8_Rank4: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_I8_Rank4")) return
   end subroutine

   subroutine Pack_I8_Rank5(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B8Ki), intent(in)            :: Data(:,:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_I8_Rank5(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B8Ki), intent(out)           :: Data(:,:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_I8_Rank5(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B8Ki), allocatable, intent(in)  :: Data(:,:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_I8_Rank5")) return
      if (.not. allocated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 5, lbound(Data), ubound(Data))

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_I8_Rank5")) return
   end subroutine

   subroutine UnpackAlloc_I8_Rank5(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B8Ki), allocatable, intent(out)  :: Data(:,:,:,:,:)
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated
      integer(B4Ki)                        :: LB(5), UB(5)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_I8_Rank5")) return
      if (.not. IsAllocated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 5, LB, UB)
      if (RegCheckErr(RF, "UnpackAlloc_I8_Rank5")) return

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3),LB(4):UB(4),LB(5):UB(5)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_I8_Rank5: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_I8_Rank5")) return
   end subroutine

   subroutine PackPtr_I8_Rank5(RF, Data)
      type(RegFile), intent(inout)         :: RF
      integer(B8Ki), pointer, intent(in)   :: Data(:,:,:,:,:)
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_I8_Rank5")) return
      if (.not. associated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 5, lbound(Data), ubound(Data))

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_I8_Rank5")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_I8_Rank5")) return
   end subroutine

   subroutine UnpackPtr_I8_Rank5(RF, Data, LB, UB)
      type(RegFile), intent(inout)          :: RF
      integer(B8Ki), pointer, intent(out)   :: Data(:,:,:,:,:)
      integer(B4Ki), intent(out)            :: LB(:), UB(:)
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_I8_Rank5")) return
      if (.not. IsAssociated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 5, LB, UB)
      if (RegCheckErr(RF, "UnpackPtr_I8_Rank5")) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_I8_Rank5")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data, UB - LB)
         Data(LB(1):,LB(2):,LB(3):,LB(4):,LB(5):) => Data
         return
      end if

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3),LB(4):UB(4),LB(5):UB(5)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_I8_Rank5: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_I8_Rank5")) return
   end subroutine

   subroutine Pack_R4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R4Ki), intent(in)               :: Data

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_R4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R4Ki), intent(out)              :: Data

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_R4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R4Ki), allocatable, intent(in)  :: Data

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_R4")) return
      if (.not. allocated(Data)) return


      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_R4")) return
   end subroutine

   subroutine UnpackAlloc_R4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R4Ki), allocatable, intent(out)  :: Data
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_R4")) return
      if (.not. IsAllocated) return


      ! Allocate data
      allocate(Data, stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_R4: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_R4")) return
   end subroutine

   subroutine PackPtr_R4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R4Ki), pointer, intent(in)      :: Data
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_R4")) return
      if (.not. associated(Data)) return

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_R4")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_R4")) return
   end subroutine

   subroutine UnpackPtr_R4(RF, Data)
      type(RegFile), intent(inout)          :: RF
      real(R4Ki), pointer, intent(out)      :: Data
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_R4")) return
      if (.not. IsAssociated) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_R4")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data)
         return
      end if

      ! Allocate data
      allocate(Data, stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_R4: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_R4")) return
   end subroutine

   subroutine Pack_R4_Rank1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R4Ki), intent(in)               :: Data(:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_R4_Rank1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R4Ki), intent(out)              :: Data(:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_R4_Rank1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R4Ki), allocatable, intent(in)  :: Data(:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_R4_Rank1")) return
      if (.not. allocated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 1, lbound(Data), ubound(Data))

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_R4_Rank1")) return
   end subroutine

   subroutine UnpackAlloc_R4_Rank1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R4Ki), allocatable, intent(out)  :: Data(:)
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated
      integer(B4Ki)                        :: LB(1), UB(1)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_R4_Rank1")) return
      if (.not. IsAllocated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 1, LB, UB)
      if (RegCheckErr(RF, "UnpackAlloc_R4_Rank1")) return

      ! Allocate data
      allocate(Data(LB(1):UB(1)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_R4_Rank1: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_R4_Rank1")) return
   end subroutine

   subroutine PackPtr_R4_Rank1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R4Ki), pointer, intent(in)      :: Data(:)
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_R4_Rank1")) return
      if (.not. associated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 1, lbound(Data), ubound(Data))

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_R4_Rank1")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_R4_Rank1")) return
   end subroutine

   subroutine UnpackPtr_R4_Rank1(RF, Data, LB, UB)
      type(RegFile), intent(inout)          :: RF
      real(R4Ki), pointer, intent(out)      :: Data(:)
      integer(B4Ki), intent(out)            :: LB(:), UB(:)
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_R4_Rank1")) return
      if (.not. IsAssociated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 1, LB, UB)
      if (RegCheckErr(RF, "UnpackPtr_R4_Rank1")) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_R4_Rank1")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data, UB - LB)
         Data(LB(1):) => Data
         return
      end if

      ! Allocate data
      allocate(Data(LB(1):UB(1)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_R4_Rank1: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_R4_Rank1")) return
   end subroutine

   subroutine Pack_R4_Rank2(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R4Ki), intent(in)               :: Data(:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_R4_Rank2(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R4Ki), intent(out)              :: Data(:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_R4_Rank2(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R4Ki), allocatable, intent(in)  :: Data(:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_R4_Rank2")) return
      if (.not. allocated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 2, lbound(Data), ubound(Data))

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_R4_Rank2")) return
   end subroutine

   subroutine UnpackAlloc_R4_Rank2(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R4Ki), allocatable, intent(out)  :: Data(:,:)
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated
      integer(B4Ki)                        :: LB(2), UB(2)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_R4_Rank2")) return
      if (.not. IsAllocated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 2, LB, UB)
      if (RegCheckErr(RF, "UnpackAlloc_R4_Rank2")) return

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_R4_Rank2: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_R4_Rank2")) return
   end subroutine

   subroutine PackPtr_R4_Rank2(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R4Ki), pointer, intent(in)      :: Data(:,:)
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_R4_Rank2")) return
      if (.not. associated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 2, lbound(Data), ubound(Data))

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_R4_Rank2")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_R4_Rank2")) return
   end subroutine

   subroutine UnpackPtr_R4_Rank2(RF, Data, LB, UB)
      type(RegFile), intent(inout)          :: RF
      real(R4Ki), pointer, intent(out)      :: Data(:,:)
      integer(B4Ki), intent(out)            :: LB(:), UB(:)
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_R4_Rank2")) return
      if (.not. IsAssociated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 2, LB, UB)
      if (RegCheckErr(RF, "UnpackPtr_R4_Rank2")) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_R4_Rank2")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data, UB - LB)
         Data(LB(1):,LB(2):) => Data
         return
      end if

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_R4_Rank2: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_R4_Rank2")) return
   end subroutine

   subroutine Pack_R4_Rank3(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R4Ki), intent(in)               :: Data(:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_R4_Rank3(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R4Ki), intent(out)              :: Data(:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_R4_Rank3(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R4Ki), allocatable, intent(in)  :: Data(:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_R4_Rank3")) return
      if (.not. allocated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 3, lbound(Data), ubound(Data))

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_R4_Rank3")) return
   end subroutine

   subroutine UnpackAlloc_R4_Rank3(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R4Ki), allocatable, intent(out)  :: Data(:,:,:)
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated
      integer(B4Ki)                        :: LB(3), UB(3)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_R4_Rank3")) return
      if (.not. IsAllocated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 3, LB, UB)
      if (RegCheckErr(RF, "UnpackAlloc_R4_Rank3")) return

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_R4_Rank3: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_R4_Rank3")) return
   end subroutine

   subroutine PackPtr_R4_Rank3(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R4Ki), pointer, intent(in)      :: Data(:,:,:)
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_R4_Rank3")) return
      if (.not. associated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 3, lbound(Data), ubound(Data))

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_R4_Rank3")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_R4_Rank3")) return
   end subroutine

   subroutine UnpackPtr_R4_Rank3(RF, Data, LB, UB)
      type(RegFile), intent(inout)          :: RF
      real(R4Ki), pointer, intent(out)      :: Data(:,:,:)
      integer(B4Ki), intent(out)            :: LB(:), UB(:)
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_R4_Rank3")) return
      if (.not. IsAssociated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 3, LB, UB)
      if (RegCheckErr(RF, "UnpackPtr_R4_Rank3")) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_R4_Rank3")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data, UB - LB)
         Data(LB(1):,LB(2):,LB(3):) => Data
         return
      end if

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_R4_Rank3: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_R4_Rank3")) return
   end subroutine

   subroutine Pack_R4_Rank4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R4Ki), intent(in)               :: Data(:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_R4_Rank4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R4Ki), intent(out)              :: Data(:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_R4_Rank4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R4Ki), allocatable, intent(in)  :: Data(:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_R4_Rank4")) return
      if (.not. allocated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 4, lbound(Data), ubound(Data))

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_R4_Rank4")) return
   end subroutine

   subroutine UnpackAlloc_R4_Rank4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R4Ki), allocatable, intent(out)  :: Data(:,:,:,:)
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated
      integer(B4Ki)                        :: LB(4), UB(4)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_R4_Rank4")) return
      if (.not. IsAllocated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 4, LB, UB)
      if (RegCheckErr(RF, "UnpackAlloc_R4_Rank4")) return

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3),LB(4):UB(4)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_R4_Rank4: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_R4_Rank4")) return
   end subroutine

   subroutine PackPtr_R4_Rank4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R4Ki), pointer, intent(in)      :: Data(:,:,:,:)
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_R4_Rank4")) return
      if (.not. associated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 4, lbound(Data), ubound(Data))

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_R4_Rank4")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_R4_Rank4")) return
   end subroutine

   subroutine UnpackPtr_R4_Rank4(RF, Data, LB, UB)
      type(RegFile), intent(inout)          :: RF
      real(R4Ki), pointer, intent(out)      :: Data(:,:,:,:)
      integer(B4Ki), intent(out)            :: LB(:), UB(:)
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_R4_Rank4")) return
      if (.not. IsAssociated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 4, LB, UB)
      if (RegCheckErr(RF, "UnpackPtr_R4_Rank4")) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_R4_Rank4")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data, UB - LB)
         Data(LB(1):,LB(2):,LB(3):,LB(4):) => Data
         return
      end if

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3),LB(4):UB(4)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_R4_Rank4: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_R4_Rank4")) return
   end subroutine

   subroutine Pack_R4_Rank5(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R4Ki), intent(in)               :: Data(:,:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_R4_Rank5(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R4Ki), intent(out)              :: Data(:,:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_R4_Rank5(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R4Ki), allocatable, intent(in)  :: Data(:,:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_R4_Rank5")) return
      if (.not. allocated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 5, lbound(Data), ubound(Data))

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_R4_Rank5")) return
   end subroutine

   subroutine UnpackAlloc_R4_Rank5(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R4Ki), allocatable, intent(out)  :: Data(:,:,:,:,:)
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated
      integer(B4Ki)                        :: LB(5), UB(5)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_R4_Rank5")) return
      if (.not. IsAllocated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 5, LB, UB)
      if (RegCheckErr(RF, "UnpackAlloc_R4_Rank5")) return

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3),LB(4):UB(4),LB(5):UB(5)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_R4_Rank5: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_R4_Rank5")) return
   end subroutine

   subroutine PackPtr_R4_Rank5(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R4Ki), pointer, intent(in)      :: Data(:,:,:,:,:)
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_R4_Rank5")) return
      if (.not. associated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 5, lbound(Data), ubound(Data))

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_R4_Rank5")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_R4_Rank5")) return
   end subroutine

   subroutine UnpackPtr_R4_Rank5(RF, Data, LB, UB)
      type(RegFile), intent(inout)          :: RF
      real(R4Ki), pointer, intent(out)      :: Data(:,:,:,:,:)
      integer(B4Ki), intent(out)            :: LB(:), UB(:)
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_R4_Rank5")) return
      if (.not. IsAssociated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 5, LB, UB)
      if (RegCheckErr(RF, "UnpackPtr_R4_Rank5")) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_R4_Rank5")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data, UB - LB)
         Data(LB(1):,LB(2):,LB(3):,LB(4):,LB(5):) => Data
         return
      end if

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3),LB(4):UB(4),LB(5):UB(5)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_R4_Rank5: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_R4_Rank5")) return
   end subroutine

   subroutine Pack_R8(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R8Ki), intent(in)               :: Data

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_R8(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R8Ki), intent(out)              :: Data

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_R8(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R8Ki), allocatable, intent(in)  :: Data

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_R8")) return
      if (.not. allocated(Data)) return


      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_R8")) return
   end subroutine

   subroutine UnpackAlloc_R8(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R8Ki), allocatable, intent(out)  :: Data
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_R8")) return
      if (.not. IsAllocated) return


      ! Allocate data
      allocate(Data, stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_R8: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_R8")) return
   end subroutine

   subroutine PackPtr_R8(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R8Ki), pointer, intent(in)      :: Data
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_R8")) return
      if (.not. associated(Data)) return

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_R8")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_R8")) return
   end subroutine

   subroutine UnpackPtr_R8(RF, Data)
      type(RegFile), intent(inout)          :: RF
      real(R8Ki), pointer, intent(out)      :: Data
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_R8")) return
      if (.not. IsAssociated) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_R8")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data)
         return
      end if

      ! Allocate data
      allocate(Data, stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_R8: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_R8")) return
   end subroutine

   subroutine Pack_R8_Rank1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R8Ki), intent(in)               :: Data(:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_R8_Rank1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R8Ki), intent(out)              :: Data(:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_R8_Rank1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R8Ki), allocatable, intent(in)  :: Data(:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_R8_Rank1")) return
      if (.not. allocated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 1, lbound(Data), ubound(Data))

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_R8_Rank1")) return
   end subroutine

   subroutine UnpackAlloc_R8_Rank1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R8Ki), allocatable, intent(out)  :: Data(:)
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated
      integer(B4Ki)                        :: LB(1), UB(1)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_R8_Rank1")) return
      if (.not. IsAllocated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 1, LB, UB)
      if (RegCheckErr(RF, "UnpackAlloc_R8_Rank1")) return

      ! Allocate data
      allocate(Data(LB(1):UB(1)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_R8_Rank1: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_R8_Rank1")) return
   end subroutine

   subroutine PackPtr_R8_Rank1(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R8Ki), pointer, intent(in)      :: Data(:)
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_R8_Rank1")) return
      if (.not. associated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 1, lbound(Data), ubound(Data))

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_R8_Rank1")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_R8_Rank1")) return
   end subroutine

   subroutine UnpackPtr_R8_Rank1(RF, Data, LB, UB)
      type(RegFile), intent(inout)          :: RF
      real(R8Ki), pointer, intent(out)      :: Data(:)
      integer(B4Ki), intent(out)            :: LB(:), UB(:)
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_R8_Rank1")) return
      if (.not. IsAssociated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 1, LB, UB)
      if (RegCheckErr(RF, "UnpackPtr_R8_Rank1")) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_R8_Rank1")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data, UB - LB)
         Data(LB(1):) => Data
         return
      end if

      ! Allocate data
      allocate(Data(LB(1):UB(1)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_R8_Rank1: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_R8_Rank1")) return
   end subroutine

   subroutine Pack_R8_Rank2(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R8Ki), intent(in)               :: Data(:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_R8_Rank2(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R8Ki), intent(out)              :: Data(:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_R8_Rank2(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R8Ki), allocatable, intent(in)  :: Data(:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_R8_Rank2")) return
      if (.not. allocated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 2, lbound(Data), ubound(Data))

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_R8_Rank2")) return
   end subroutine

   subroutine UnpackAlloc_R8_Rank2(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R8Ki), allocatable, intent(out)  :: Data(:,:)
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated
      integer(B4Ki)                        :: LB(2), UB(2)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_R8_Rank2")) return
      if (.not. IsAllocated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 2, LB, UB)
      if (RegCheckErr(RF, "UnpackAlloc_R8_Rank2")) return

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_R8_Rank2: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_R8_Rank2")) return
   end subroutine

   subroutine PackPtr_R8_Rank2(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R8Ki), pointer, intent(in)      :: Data(:,:)
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_R8_Rank2")) return
      if (.not. associated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 2, lbound(Data), ubound(Data))

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_R8_Rank2")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_R8_Rank2")) return
   end subroutine

   subroutine UnpackPtr_R8_Rank2(RF, Data, LB, UB)
      type(RegFile), intent(inout)          :: RF
      real(R8Ki), pointer, intent(out)      :: Data(:,:)
      integer(B4Ki), intent(out)            :: LB(:), UB(:)
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_R8_Rank2")) return
      if (.not. IsAssociated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 2, LB, UB)
      if (RegCheckErr(RF, "UnpackPtr_R8_Rank2")) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_R8_Rank2")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data, UB - LB)
         Data(LB(1):,LB(2):) => Data
         return
      end if

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_R8_Rank2: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_R8_Rank2")) return
   end subroutine

   subroutine Pack_R8_Rank3(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R8Ki), intent(in)               :: Data(:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_R8_Rank3(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R8Ki), intent(out)              :: Data(:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_R8_Rank3(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R8Ki), allocatable, intent(in)  :: Data(:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_R8_Rank3")) return
      if (.not. allocated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 3, lbound(Data), ubound(Data))

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_R8_Rank3")) return
   end subroutine

   subroutine UnpackAlloc_R8_Rank3(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R8Ki), allocatable, intent(out)  :: Data(:,:,:)
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated
      integer(B4Ki)                        :: LB(3), UB(3)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_R8_Rank3")) return
      if (.not. IsAllocated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 3, LB, UB)
      if (RegCheckErr(RF, "UnpackAlloc_R8_Rank3")) return

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_R8_Rank3: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_R8_Rank3")) return
   end subroutine

   subroutine PackPtr_R8_Rank3(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R8Ki), pointer, intent(in)      :: Data(:,:,:)
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_R8_Rank3")) return
      if (.not. associated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 3, lbound(Data), ubound(Data))

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_R8_Rank3")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_R8_Rank3")) return
   end subroutine

   subroutine UnpackPtr_R8_Rank3(RF, Data, LB, UB)
      type(RegFile), intent(inout)          :: RF
      real(R8Ki), pointer, intent(out)      :: Data(:,:,:)
      integer(B4Ki), intent(out)            :: LB(:), UB(:)
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_R8_Rank3")) return
      if (.not. IsAssociated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 3, LB, UB)
      if (RegCheckErr(RF, "UnpackPtr_R8_Rank3")) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_R8_Rank3")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data, UB - LB)
         Data(LB(1):,LB(2):,LB(3):) => Data
         return
      end if

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_R8_Rank3: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_R8_Rank3")) return
   end subroutine

   subroutine Pack_R8_Rank4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R8Ki), intent(in)               :: Data(:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_R8_Rank4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R8Ki), intent(out)              :: Data(:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_R8_Rank4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R8Ki), allocatable, intent(in)  :: Data(:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_R8_Rank4")) return
      if (.not. allocated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 4, lbound(Data), ubound(Data))

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_R8_Rank4")) return
   end subroutine

   subroutine UnpackAlloc_R8_Rank4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R8Ki), allocatable, intent(out)  :: Data(:,:,:,:)
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated
      integer(B4Ki)                        :: LB(4), UB(4)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_R8_Rank4")) return
      if (.not. IsAllocated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 4, LB, UB)
      if (RegCheckErr(RF, "UnpackAlloc_R8_Rank4")) return

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3),LB(4):UB(4)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_R8_Rank4: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_R8_Rank4")) return
   end subroutine

   subroutine PackPtr_R8_Rank4(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R8Ki), pointer, intent(in)      :: Data(:,:,:,:)
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_R8_Rank4")) return
      if (.not. associated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 4, lbound(Data), ubound(Data))

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_R8_Rank4")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_R8_Rank4")) return
   end subroutine

   subroutine UnpackPtr_R8_Rank4(RF, Data, LB, UB)
      type(RegFile), intent(inout)          :: RF
      real(R8Ki), pointer, intent(out)      :: Data(:,:,:,:)
      integer(B4Ki), intent(out)            :: LB(:), UB(:)
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_R8_Rank4")) return
      if (.not. IsAssociated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 4, LB, UB)
      if (RegCheckErr(RF, "UnpackPtr_R8_Rank4")) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_R8_Rank4")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data, UB - LB)
         Data(LB(1):,LB(2):,LB(3):,LB(4):) => Data
         return
      end if

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3),LB(4):UB(4)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_R8_Rank4: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_R8_Rank4")) return
   end subroutine

   subroutine Pack_R8_Rank5(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R8Ki), intent(in)               :: Data(:,:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Increment data number and write to file
      RF%NumData = RF%NumData + 1
      write(RF%Unit) RF%NumData

      ! Write data to file
      write(RF%Unit) Data
   end subroutine

   subroutine Unpack_R8_Rank5(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R8Ki), intent(out)              :: Data(:,:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Read data number, return if invalid
      if (.not. DataNumValid(RF)) return

      ! Read data from file
      read(RF%Unit) Data
   end subroutine

   subroutine PackAlloc_R8_Rank5(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R8Ki), allocatable, intent(in)  :: Data(:,:,:,:,:)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if allocated
      call RegPack(RF, allocated(Data))
      if (RegCheckErr(RF, "PackAlloc_R8_Rank5")) return
      if (.not. allocated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 5, lbound(Data), ubound(Data))

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackAlloc_R8_Rank5")) return
   end subroutine

   subroutine UnpackAlloc_R8_Rank5(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R8Ki), allocatable, intent(out)  :: Data(:,:,:,:,:)
      integer(IntKi)                       :: stat
      logical                              :: IsAllocated
      integer(B4Ki)                        :: LB(5), UB(5)

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Deallocate if allocated
      if (allocated(Data)) deallocate(Data)

      ! Read value to see if it was allocated, return if not
      call RegUnpack(RF, IsAllocated)
      if (RegCheckErr(RF, "UnpackAlloc_R8_Rank5")) return
      if (.not. IsAllocated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 5, LB, UB)
      if (RegCheckErr(RF, "UnpackAlloc_R8_Rank5")) return

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3),LB(4):UB(4),LB(5):UB(5)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackAlloc_R8_Rank5: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackAlloc_R8_Rank5")) return
   end subroutine

   subroutine PackPtr_R8_Rank5(RF, Data)
      type(RegFile), intent(inout)         :: RF
      real(R8Ki), pointer, intent(in)      :: Data(:,:,:,:,:)
      logical                              :: PtrInIndex

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! Write if associated
      call RegPack(RF, associated(Data))
      if (RegCheckErr(RF, "PackPtr_R8_Rank5")) return
      if (.not. associated(Data)) return

      ! Write array bounds
      call RegPackBounds(RF, 5, lbound(Data), ubound(Data))

      ! Write pointer info
      call RegPackPointer(RF, c_loc(Data), PtrInIndex)
      if (RegCheckErr(RF, "PackPtr_R8_Rank5")) return
      if (PtrInIndex) return

      ! Write data to file
      call RegPack(RF, Data)
      if (RegCheckErr(RF, "PackPtr_R8_Rank5")) return
   end subroutine

   subroutine UnpackPtr_R8_Rank5(RF, Data, LB, UB)
      type(RegFile), intent(inout)          :: RF
      real(R8Ki), pointer, intent(out)      :: Data(:,:,:,:,:)
      integer(B4Ki), intent(out)            :: LB(:), UB(:)
      integer(IntKi)                        :: stat
      integer(B8Ki)                         :: PtrIdx
      logical                               :: IsAssociated
      type(c_ptr)                           :: Ptr

      ! If error, return
      if (RF%ErrStat /= ErrID_None) return

      ! If associated, deallocate and nullify
      if (associated(Data)) then
         deallocate(Data)
         nullify(Data)
      end if

      ! Read value to see if it was associated, return if not
      call RegUnpack(RF, IsAssociated)
      if (RegCheckErr(RF, "UnpackPtr_R8_Rank5")) return
      if (.not. IsAssociated) return

      ! Read array bounds
      call RegUnpackBounds(RF, 5, LB, UB)
      if (RegCheckErr(RF, "UnpackPtr_R8_Rank5")) return

      ! Unpack pointer inf
      call RegUnpackPointer(RF, Ptr, PtrIdx)
      if (RegCheckErr(RF, "UnpackPtr_R8_Rank5")) return

      ! If pointer was in index, associate data with pointer, return
      if (c_associated(Ptr)) then
         call c_f_pointer(Ptr, Data, UB - LB)
         Data(LB(1):,LB(2):,LB(3):,LB(4):,LB(5):) => Data
         return
      end if

      ! Allocate data
      allocate(Data(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3),LB(4):UB(4),LB(5):UB(5)), stat=stat)
      if (stat /= 0) then
         RF%ErrStat = ErrID_Fatal
         RF%ErrMsg = "UnpackPtr_R8_Rank5: error allocating data"
         return
      end if

      ! Read data
      call RegUnpack(RF, Data)
      if (RegCheckErr(RF, "UnpackPtr_R8_Rank5")) return
   end subroutine
end module