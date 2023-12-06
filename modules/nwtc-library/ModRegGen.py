
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
module ModReg
   use NWTC_Base
   implicit none

   private
   public :: PackBuffer
   public :: WritePackBuffer, ReadPackBuffer, InitPackBuffer, DestroyPackBuffer, RegCheckErr
   public :: RegPack, RegPackBounds, RegPackPointer
   public :: RegUnpack, RegUnpackBounds, RegUnpackPointer

   type :: PackBuffer
      integer(B1Ki), allocatable  :: Bytes(:)
      integer(B8Ki)               :: NB
      type(c_ptr), allocatable    :: Pointers(:)
      integer(B8Ki)               :: NP
      integer(IntKi)              :: ErrStat = ErrID_Fatal
      character(ErrMsgLen)        :: ErrMsg = 'PackBuffer not initialized'
   end type
   {ifc_lines}

contains

   subroutine InitPackBuffer(Buf, ErrStat, ErrMsg)
      type(PackBuffer), intent(inout)   :: Buf
      integer(IntKi), intent(out)       :: ErrStat
      character(ErrMsgLen), intent(out) :: ErrMsg

      character(*), parameter           :: RoutineName = "InitPackBuffer"
      integer(B8Ki), parameter          :: NumPointersInit = 128
      integer(B8Ki), parameter          :: NumBytesInit = 1024
      integer(IntKi)                    :: stat

      ErrStat = ErrID_None
      ErrMsg = ""

      Buf%ErrStat = ErrID_None
      Buf%ErrMsg = ""
      Buf%NP = 0
      Buf%NB = 0
      
      ! If pointers have not been allocated, allocate with initial size
      if (.not. allocated(Buf%Pointers)) then
         allocate (Buf%Pointers(NumPointersInit), stat=stat)
         if (stat /= 0) then
            ErrStat = ErrID_Fatal
            write(ErrMsg,*) 'InitPackBuffer: Unable to init pointer index to with size of', NumPointersInit
            return
         end if
      end if
      
      ! Reset all pointers to null
      Buf%Pointers = c_null_ptr

      ! If byte array has not been allocated, allocate with initial size
      if (.not. allocated(Buf%Bytes)) then
         allocate (Buf%Bytes(NumBytesInit), stat=stat)
         if (stat /= 0) then
            ErrStat = ErrID_Fatal
            write(ErrMsg,*) 'Grow: Unable to init buffer to', NumBytesInit, 'bytes'
            return
         end if
      end if

   end subroutine

   subroutine DestroyPackBuffer(Buf, ErrStat, ErrMsg)
      type(PackBuffer), intent(inout)   :: Buf
      integer(IntKi), intent(out)       :: ErrStat
      character(ErrMsgLen), intent(out) :: ErrMsg

      character(*), parameter           :: RoutineName = "DestroyPackBuffer"

      ErrStat = ErrID_None
      ErrMsg = ""

      Buf%ErrStat = ErrID_None
      Buf%ErrMsg = ""
      Buf%NP = 0
      Buf%NB = 0
      
      if (allocated(Buf%Pointers)) deallocate (Buf%Pointers)
      if (allocated(Buf%Bytes   )) deallocate (Buf%Bytes)
   end subroutine

   subroutine WritePackBuffer(Buf, Unit, ErrStat, ErrMsg)
      type(PackBuffer), intent(inout)   :: Buf
      integer(IntKi), intent(in)        :: Unit
      integer(IntKi), intent(out)       :: ErrStat
      character(ErrMsgLen), intent(out) :: ErrMsg

      character(*), parameter           :: RoutineName = "WritePackBuffer"
      integer(IntKi)                    :: iostat

      ErrStat = ErrID_None
      ErrMsg = ''

      if (Buf%ErrStat /= ErrID_None) then
         call SetErrStat(Buf%ErrStat, Buf%ErrMsg, ErrStat, ErrMsg, 'Buf%WriteFile')
         return
      end if

      write(Unit, iostat=iostat) Buf%NP
      if (iostat /= 0) then
         call SetErrStat(ErrID_Fatal, "Error writing number of pointers", ErrStat, ErrMsg, RoutineName)
         return
      end if

      write(Unit, iostat=iostat) Buf%NB
      if (iostat /= 0) then
         call SetErrStat(ErrID_Fatal, "Error writing number of bytes", ErrStat, ErrMsg, RoutineName)
         return
      end if

      write(Unit, iostat=iostat) Buf%Bytes(1:Buf%NB)
      if (iostat /= 0) then
         call SetErrStat(ErrID_Fatal, "Error writing bytes", ErrStat, ErrMsg, RoutineName)
         return
      end if

   end subroutine

   subroutine ReadPackBuffer(Buf, Unit, ErrStat, ErrMsg)
      type(PackBuffer), intent(inout)   :: Buf
      integer(IntKi), intent(in)        :: Unit
      integer(IntKi), intent(out)       :: ErrStat
      character(ErrMsgLen), intent(out) :: ErrMsg

      character(*), parameter           :: RoutineName = "ReadPackBuffer"
      integer(IntKi)                    :: iostat

      ErrStat = ErrID_None
      ErrMsg = ''

      ! Read number of pointers
      read(Unit, iostat=iostat) Buf%NP
      if (iostat /= 0) then
         call SetErrStat(ErrID_Fatal, "Error reading number of pointers", ErrStat, ErrMsg, RoutineName)
         return
      end if

      ! If pointers are allocated, deallocate
      if (allocated(Buf%Pointers)) deallocate(Buf%Pointers)

      ! Allocate pointer index and initialize pointers to null
      allocate(Buf%Pointers(1:Buf%NP), stat=ErrStat)
      Buf%Pointers = c_null_ptr

      ! Read number of bytes in buffer
      read(Unit, iostat=iostat) Buf%NB
      if (iostat /= 0) then
         call SetErrStat(ErrID_Fatal, "Error reading number of bytes", ErrStat, ErrMsg, RoutineName)
         return
      end if

      ! If bytes are allocated, deallocate
      if (allocated(Buf%Bytes)) deallocate(Buf%Bytes)

      ! Allocate bytes
      allocate(Buf%Bytes(1:Buf%NB), stat=ErrStat)

      ! Read bytes
      read(Unit, iostat=iostat) Buf%Bytes
      if (iostat /= 0) then
         call SetErrStat(ErrID_Fatal, "Error reading bytes", ErrStat, ErrMsg, RoutineName)
         return
      end if

      ! Clear buffer error
      Buf%ErrStat = ErrID_None
      Buf%ErrMsg = ''

      ! Reset Number of bytes to be used by unpack routines
      Buf%NB = 0

   end subroutine

   function RegCheckErr(Buf, RoutineName) result(Err)
      type(PackBuffer), intent(inout)  :: Buf
      character(*), intent(in)         :: RoutineName
      logical                          :: Err
      Err = Buf%ErrStat /= ErrID_None
      if (Err) Buf%ErrMsg = trim(RoutineName)//": "//trim(Buf%ErrMsg)
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

   subroutine RegPackPointer(Buf, Ptr, Found)
      type(PackBuffer), intent(inout)   :: Buf
      type(c_ptr), intent(in)           :: Ptr
      logical, intent(out)              :: Found

      type(c_ptr), allocatable          :: PointersTmp(:)
      integer(B8Ki)                     :: NewSize
      integer(B8Ki)                     :: i

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Look for pointer in index, if found, pack pointer index and return
      do i = 1, Buf%NP
         if (c_associated(Ptr, Buf%Pointers(i))) then
            call RegPack(Buf, i)
            Found = .true.
            return
         end if
      end do

      ! Pointer was not found in index
      Found = .false.

      ! If pointer index is full, grow pointer index
      if (Buf%NP == size(Buf%Pointers)) then
         NewSize = int(1.5_R8Ki * real(Buf%NP, R8Ki), B8Ki)
         call move_alloc(Buf%Pointers, PointersTmp)
         allocate (Buf%Pointers(NewSize), stat=Buf%ErrStat)
         if (Buf%ErrStat /= ErrID_None) then
            Buf%ErrStat = ErrID_Fatal
            write(Buf%ErrMsg,*) 'RegPackPointer: Unable to allocate pointer index to', NewSize, 'bytes'
            return
         end if
         Buf%Pointers(1:size(PointersTmp)) = PointersTmp
         Buf%Pointers(size(PointersTmp)+1:) = c_null_ptr
      end if

      ! Increment number of pointers, add new pointer to index
      Buf%NP = Buf%NP + 1
      Buf%Pointers(Buf%NP) = Ptr

      ! Pack pointer index
      call RegPack(Buf, Buf%NP)

   end subroutine

   subroutine RegUnpackPointer(Buf, Ptr, Idx)
      type(PackBuffer), intent(inout)   :: Buf
      type(c_ptr), intent(out)          :: Ptr
      integer(B8Ki), intent(out)        :: Idx

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Unpack pointer index
      call RegUnpack(Buf, Idx)

      ! Get pointer from index
      Ptr = Buf%Pointers(Idx)

   end subroutine

   subroutine RegPackBounds(Buf, R, LB, UB)
      type(PackBuffer), intent(inout)  :: Buf
      integer(B4Ki), intent(in)        :: R
      integer(B8Ki), intent(in)        :: LB(:), UB(:)
      
      ! If buffer has an error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Pack lower and upper bounds
      call RegPack(Buf, LB(1:R))
      call RegPack(Buf, UB(1:R))
      if (RegCheckErr(Buf, "RegPackBounds")) return
   end subroutine

   subroutine RegUnpackBounds(Buf, R, LB, UB)
      type(PackBuffer), intent(inout)   :: Buf
      integer(B4Ki), intent(in)         :: R
      integer(B8Ki), intent(out)        :: LB(:), UB(:)

      ! If buffer has an error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Unpack lower and upper bounds
      call RegUnpack(Buf, LB(1:R))
      call RegUnpack(Buf, UB(1:R))
      if (RegCheckErr(Buf, "RegUnpackBounds")) return
   end subroutine

   subroutine GrowBuffer(Buf, N)
      type(PackBuffer), intent(inout)   :: Buf
      integer(B8Ki), intent(in)         :: N

      integer(B1Ki), allocatable        :: BytesTmp(:)
      integer(B8Ki)                     :: NewSize
      integer(IntKi)                    :: stat
      
      ! Return if there is a buffer error
      if (Buf%ErrStat /= ErrID_None) return

      ! If buffer can hold requested bytes, return
      if (size(Buf%Bytes) > Buf%NB + N) return

      ! Calculate new size
      NewSize = int(real(Buf%NB + N, R8Ki) * 1.8_R8Ki, B8Ki)

      ! Move allocation to temporary array and allocate buffer with new size
      call move_alloc(Buf%Bytes, BytesTmp)
      allocate (Buf%Bytes(NewSize), stat=stat)
      if (stat /= 0) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) 'Grow: Unable to grow buffer to', NewSize, 'bytes'
         return
      end if

      ! Copy contents of temporary bytes to buffer
      Buf%Bytes(1:size(BytesTmp)) = BytesTmp

   end subroutine
'''


def gen_pack(w, dt, decl, rank):
   dims = '' if rank == 0 else '('+','.join([':']*rank)+')'
   dt_size = int(dt[-1])
   name = f'Pack_{dt}' if rank == 0 else f'Pack_{dt}_Rank{rank}'
   w.write(f'\n\n   subroutine {name}(Buf, Data)')
   w.write(f'\n      type(PackBuffer), intent(inout)         :: Buf')
   w.write(f'\n      {decl+", intent(in)":<38s}  :: Data{dims}')
   w.write(f'\n      integer(B8Ki)                           :: DataSize')
   w.write(f'\n')
   w.write(f'\n      ! If buffer error, return')
   w.write(f'\n      if (Buf%ErrStat /= ErrID_None) return')
   w.write(f'\n')
   w.write(f'\n      ! Get size of data in bytes')
   if dt == 'C1' and rank == 0:
      w.write(f'\n      DataSize = len(Data)')
   elif dt == 'C1' and rank > 0:
      w.write(f'\n      DataSize = len(Data({",".join(["1"]*rank)}))*size(Data)')
   elif rank == 0:
      w.write(f'\n      DataSize = {dt_size}')
   elif dt_size == 1:
      w.write(f'\n      DataSize = size(Data)')
   else:
      w.write(f'\n      DataSize = {dt_size}*size(Data)')
   w.write(f'\n')
   w.write(f'\n      ! Grow buffer to accommodate Data')
   w.write(f'\n      call GrowBuffer(Buf, DataSize)')
   w.write(f'\n      if (RegCheckErr(Buf, "{name}")) return')
   w.write(f'\n')
   w.write(f'\n      ! Transfer data to buffer')
   if dt == 'L1':
      w.write(f'\n      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(LogicalToByte(Data), Buf%Bytes)')
   else:
      w.write(f'\n      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(Data, Buf%Bytes)')
   w.write(f'\n      Buf%NB = Buf%NB + DataSize')
   w.write(f'\n')
   w.write(f'\n   end subroutine')


def gen_unpack(w, dt, decl, rank):
   dims = '' if rank == 0 else '('+','.join([':']*rank)+')'
   dt_size = int(dt[-1])
   name = f'Unpack_{dt}' if rank == 0 else f'Unpack_{dt}_Rank{rank}'
   w.write(f'\n')
   w.write(f'\n   subroutine {name}(Buf, Data)')
   w.write(f'\n      type(PackBuffer), intent(inout)         :: Buf')
   w.write(f'\n      {decl+", intent(out)":<38s}  :: Data{dims}')
   w.write(f'\n      integer(B8Ki)                           :: DataSize')
   w.write(f'\n')
   w.write(f'\n      ! If buffer error, return')
   w.write(f'\n      if (Buf%ErrStat /= ErrID_None) return')
   w.write(f'\n')
   w.write(f'\n      ! Get size of data in bytes')
   if dt == 'C1' and rank == 0:
      w.write(f'\n      DataSize = len(Data)')
   elif dt == 'C1' and rank > 0:
      w.write(f'\n      DataSize = len(Data({",".join(["1"]*rank)}))*size(Data)')
   elif rank == 0:
      w.write(f'\n      DataSize = {dt_size}')
   elif dt_size == 1:
      w.write(f'\n      DataSize = size(Data)')
   else:
      w.write(f'\n      DataSize = {dt_size}*size(Data)')
   w.write(f'\n')
   w.write(f'\n      ! Check that buffer has sufficient bytes remaining')
   w.write(f'\n      if (size(Buf%Bytes) < Buf%NB + DataSize) then')
   w.write(f'\n         Buf%ErrStat = ErrID_Fatal')
   w.write(f'\n         write(Buf%ErrMsg,*) "{name}: buffer too small, requested", DataSize, "bytes"')
   w.write(f'\n         return')
   w.write(f'\n      end if')
   w.write(f'\n')
   w.write(f'\n      ! Transfer data from buffer')
   if dt == 'L1' and rank == 0:
      w.write(f'\n      Data = ByteToLogical(Buf%Bytes(Buf%NB+1))')
   elif dt == 'L1' and rank > 0:
      w.write(f'\n      Data = reshape(ByteToLogical(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize)), shape(Data))')
   elif rank == 0:
      w.write(f'\n      Data = transfer(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize), Data)')
   else:
      w.write(f'\n      Data = reshape(transfer(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize), Data), shape(Data))')
   w.write(f'\n      Buf%NB = Buf%NB + DataSize')
   w.write(f'\n')
   w.write(f'\n   end subroutine')

def gen_pack_alloc(w, dt, decl, rank):
   dims = '' if rank == 0 else '('+','.join([':']*rank)+')'
   name = f'PackAlloc_{dt}'
   if rank > 0: name += f'_Rank{rank}'
   w.write(f'\n')
   w.write(f'\n   subroutine {name}(Buf, Data)')
   w.write(f'\n      type(PackBuffer), intent(inout)         :: Buf')
   w.write(f'\n      {decl+", allocatable, intent(in)":<38s}  :: Data{dims}')
   w.write(f'\n')
   w.write(f'\n      ! If buffer error, return')
   w.write(f'\n      if (Buf%ErrStat /= ErrID_None) return')
   w.write(f'\n')
   w.write(f'\n      ! Write if allocated')
   w.write(f'\n      call RegPack(Buf, allocated(Data))')
   w.write(f'\n      if (RegCheckErr(Buf, "{name}")) return')
   w.write(f'\n      if (.not. allocated(Data)) return')
   w.write(f'\n')
   if rank > 0:
      w.write(f'\n      ! Write array bounds')
      w.write(f'\n      call RegPackBounds(Buf, {rank}, lbound(Data, kind=B8Ki), ubound(Data, kind=B8Ki))')
   w.write(f'\n')
   w.write(f'\n      ! Write data to buffer')
   w.write(f'\n      call RegPack(Buf, Data)')
   w.write(f'\n      if (RegCheckErr(Buf, "{name}")) return')
   w.write(f'\n')
   w.write(f'\n   end subroutine')


def gen_unpack_alloc(w, dt, decl, rank):
   dims = '' if rank == 0 else '('+','.join([':']*rank)+')'
   dt_size = int(dt[-1])
   name = f'UnpackAlloc_{dt}' if rank == 0 else f'UnpackAlloc_{dt}_Rank{rank}'
   w.write(f'\n')
   w.write(f'\n   subroutine {name}(Buf, Data)')
   w.write(f'\n      type(PackBuffer), intent(inout)         :: Buf')
   w.write(f'\n      {decl+", allocatable, intent(out)":<38s}  :: Data{dims}')
   w.write(f'\n      integer(IntKi)                          :: stat')
   w.write(f'\n      logical                               :: IsAllocated')
   if rank > 0:
      w.write(f'\n      integer(B8Ki)                          :: LB({rank}), UB({rank})')
   w.write(f'\n')
   w.write(f'\n      ! If buffer error, return')
   w.write(f'\n      if (Buf%ErrStat /= ErrID_None) return')
   w.write(f'\n')
   w.write(f'\n      ! Deallocate if allocated')
   w.write(f'\n      if (allocated(Data)) deallocate(Data)')
   w.write(f'\n')
   w.write(f'\n      ! Read value to see if it was allocated, return if not')
   w.write(f'\n      call RegUnpack(Buf, IsAllocated)')
   w.write(f'\n      if (RegCheckErr(Buf, "{name}")) return')
   w.write(f'\n      if (.not. IsAllocated) return')
   w.write(f'\n')
   alloc_dims = ''
   if rank > 0:
      w.write(f'\n      ! Read array bounds')
      w.write(f'\n      call RegUnpackBounds(Buf, {rank}, LB, UB)')
      w.write(f'\n      if (RegCheckErr(Buf, "{name}")) return')
      alloc_dims = '(' + ','.join([f'LB({d+1}):UB({d+1})' for d in range(rank)]) + ')'
   w.write(f'\n')
   w.write(f'\n      ! Allocate data')
   w.write(f'\n      allocate(Data{alloc_dims}, stat=stat)')
   w.write(f'\n      if (stat /= 0) then')
   w.write(f'\n         Buf%ErrStat = ErrID_Fatal')
   w.write(f'\n         Buf%ErrMsg = "{name}: error allocating data"')
   w.write(f'\n         return')
   w.write(f'\n      end if')
   w.write(f'\n')
   w.write(f'\n      ! Read data')
   w.write(f'\n      call RegUnpack(Buf, Data)')
   w.write(f'\n      if (RegCheckErr(Buf, "{name}")) return')
   w.write(f'\n')
   w.write(f'\n   end subroutine')


def gen_pack_ptr(w, dt, decl, rank):
   dims = '' if rank == 0 else '('+','.join([':']*rank)+')'
   name = f'PackPtr_{dt}'
   if rank > 0: name += f'_Rank{rank}'
   w.write(f'\n')
   w.write(f'\n   subroutine {name}(Buf, Data)')
   w.write(f'\n      type(PackBuffer), intent(inout)         :: Buf')
   w.write(f'\n      {decl+", pointer, intent(in)":<38s}  :: Data{dims}')
   w.write(f'\n      logical                               :: PtrInIndex')
   w.write(f'\n')
   w.write(f'\n      ! If buffer error, return')
   w.write(f'\n      if (Buf%ErrStat /= ErrID_None) return')
   w.write(f'\n')
   w.write(f'\n      ! Write if associated')
   w.write(f'\n      call RegPack(Buf, associated(Data))')
   w.write(f'\n      if (RegCheckErr(Buf, "{name}")) return')
   w.write(f'\n      if (.not. associated(Data)) return')
   if rank > 0:
      w.write(f'\n')
      w.write(f'\n      ! Write array bounds')
      w.write(f'\n      call RegPackBounds(Buf, {rank}, lbound(Data, kind=B8Ki), ubound(Data, kind=B8Ki))')
   w.write(f'\n')
   w.write(f'\n      ! Write pointer info')
   w.write(f'\n      call RegPackPointer(Buf, c_loc(Data), PtrInIndex)')
   w.write(f'\n      if (RegCheckErr(Buf, "{name}")) return')
   w.write(f'\n      if (PtrInIndex) return')
   w.write(f'\n')
   w.write(f'\n      ! Write data to buffer')
   w.write(f'\n      call RegPack(Buf, Data)')
   w.write(f'\n      if (RegCheckErr(Buf, "{name}")) return')
   w.write(f'\n')
   w.write(f'\n   end subroutine')

def gen_unpack_ptr(w, dt, decl, rank):
   dims = '' if rank == 0 else '('+','.join([':']*rank)+')'
   dt_size = int(dt[-1])
   name = f'UnpackPtr_{dt}' if rank == 0 else f'UnpackPtr_{dt}_Rank{rank}'
   w.write(f'\n')
   w.write(f'\n   subroutine {name}(Buf, Data)')
   w.write(f'\n      type(PackBuffer), intent(inout)         :: Buf')
   w.write(f'\n      {decl+", pointer, intent(out)":<38s}  :: Data{dims}')
   w.write(f'\n      integer(B8Ki)                           :: PtrIdx, stat')
   w.write(f'\n      logical                               :: IsAssociated')
   w.write(f'\n      type(c_ptr)                           :: Ptr')
   if rank > 0:
      w.write(f'\n      integer(B8Ki)                         :: LB({rank}), UB({rank})')
   w.write(f'\n')
   w.write(f'\n      ! If buffer error, return')
   w.write(f'\n      if (Buf%ErrStat /= ErrID_None) return')
   w.write(f'\n')
   w.write(f'\n      ! If associated, deallocate and nullify')
   w.write(f'\n      if (associated(Data)) then')
   w.write(f'\n         deallocate(Data)')
   w.write(f'\n         nullify(Data)')
   w.write(f'\n      end if')
   w.write(f'\n')
   w.write(f'\n      ! Read value to see if it was associated, return if not')
   w.write(f'\n      call RegUnpack(Buf, IsAssociated)')
   w.write(f'\n      if (RegCheckErr(Buf, "{name}")) return')
   w.write(f'\n      if (.not. IsAssociated) return')
   if rank > 0:
      w.write(f'\n')
      w.write(f'\n      ! Read array bounds')
      w.write(f'\n      call RegUnpackBounds(Buf, {rank}, LB, UB)')
   w.write(f'\n')
   w.write(f'\n      ! Unpack pointer inf')
   w.write(f'\n      call RegUnpackPointer(Buf, Ptr, PtrIdx)')
   w.write(f'\n      if (RegCheckErr(Buf, "{name}")) return')
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
   w.write(f'\n         Buf%ErrStat = ErrID_Fatal')
   w.write(f'\n         Buf%ErrMsg = "{name}: error allocating data"')
   w.write(f'\n         return')
   w.write(f'\n      end if')
   w.write(f'\n')
   w.write(f'\n      ! Read data')
   w.write(f'\n      call RegUnpack(Buf, Data)')
   w.write(f'\n      if (RegCheckErr(Buf, "{name}")) return')
   w.write(f'\n')
   w.write(f'\n   end subroutine')

# Registry interface
ifc_lines = ''
ranks = [''] + [f'_Rank{r}' for r in range(1,num_ranks+1)]
for attr, punp in product([''], ['Pack', 'Unpack']):
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

    w.write('\nend module')
