
module ModReg
   use NWTC_Base
   implicit none

   private
   public :: PackBuffer
   public :: WritePackBuffer, ReadPackBuffer, InitPackBuffer, RegCheckErr
   public :: RegPack, RegPackBounds, RegPackPointer
   public :: RegUnpack, RegUnpackBounds, RegUnpackPointer

   type :: PackBuffer
      integer(B1Ki), allocatable  :: Bytes(:)
      integer(IntKi)              :: NB
      type(c_ptr), allocatable    :: Pointers(:)
      integer(IntKi)              :: NP
      integer(IntKi)              :: ErrStat = ErrID_Fatal
      character(ErrMsgLen)        :: ErrMsg = 'PackBuffer not initialized'
   end type
   

   interface RegPack
      module procedure Pack_C1, Pack_C1_Rank1, Pack_C1_Rank2, Pack_C1_Rank3, &
         Pack_C1_Rank4, Pack_C1_Rank5, Pack_L1, Pack_L1_Rank1, Pack_L1_Rank2, &
         Pack_L1_Rank3, Pack_L1_Rank4, Pack_L1_Rank5, Pack_I4, Pack_I4_Rank1, &
         Pack_I4_Rank2, Pack_I4_Rank3, Pack_I4_Rank4, Pack_I4_Rank5, Pack_R4, &
         Pack_R4_Rank1, Pack_R4_Rank2, Pack_R4_Rank3, Pack_R4_Rank4, &
         Pack_R4_Rank5, Pack_R8, Pack_R8_Rank1, Pack_R8_Rank2, Pack_R8_Rank3, &
         Pack_R8_Rank4, Pack_R8_Rank5
   end interface

   interface RegUnpack
      module procedure Unpack_C1, Unpack_C1_Rank1, Unpack_C1_Rank2, &
         Unpack_C1_Rank3, Unpack_C1_Rank4, Unpack_C1_Rank5, Unpack_L1, &
         Unpack_L1_Rank1, Unpack_L1_Rank2, Unpack_L1_Rank3, Unpack_L1_Rank4, &
         Unpack_L1_Rank5, Unpack_I4, Unpack_I4_Rank1, Unpack_I4_Rank2, &
         Unpack_I4_Rank3, Unpack_I4_Rank4, Unpack_I4_Rank5, Unpack_R4, &
         Unpack_R4_Rank1, Unpack_R4_Rank2, Unpack_R4_Rank3, Unpack_R4_Rank4, &
         Unpack_R4_Rank5, Unpack_R8, Unpack_R8_Rank1, Unpack_R8_Rank2, &
         Unpack_R8_Rank3, Unpack_R8_Rank4, Unpack_R8_Rank5
   end interface

contains

   subroutine InitPackBuffer(Buf, ErrStat, ErrMsg)
      type(PackBuffer), intent(inout)   :: Buf
      integer(IntKi), intent(out)       :: ErrStat
      character(ErrMsgLen), intent(out) :: ErrMsg

      character(*), parameter           :: RoutineName = "InitPackBuffer"
      integer(IntKi), parameter         :: NumPointersInit = 128
      integer(IntKi), parameter         :: NumBytesInit = 1024
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
      integer(IntKi)                    :: NewSize
      integer(B4Ki)                     :: i

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
         NewSize = int(1.5_R4Ki * real(Buf%NP, R4Ki), IntKi)
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
      integer(B4Ki), intent(out)        :: Idx

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Unpack pointer index
      call RegUnpack(Buf, Idx)

      ! Get pointer from index
      Ptr = Buf%Pointers(Idx)

   end subroutine

   subroutine RegPackBounds(Buf, R, LB, UB)
      type(PackBuffer), intent(inout)  :: Buf
      integer(B4Ki), intent(in)        :: R, LB(:), UB(:)
      
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
      integer(B4Ki), intent(out)        :: LB(:), UB(:)

      ! If buffer has an error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Unpack lower and upper bounds
      call RegUnpack(Buf, LB(1:R))
      call RegUnpack(Buf, UB(1:R))
      if (RegCheckErr(Buf, "RegUnpackBounds")) return
   end subroutine

   subroutine GrowBuffer(Buf, N)
      type(PackBuffer), intent(inout)   :: Buf
      integer(B4Ki), intent(in)         :: N

      integer(B1Ki), allocatable        :: BytesTmp(:)
      integer(B4Ki)                     :: NewSize
      integer(IntKi)                    :: stat
      
      ! Return if there is a buffer error
      if (Buf%ErrStat /= ErrID_None) return

      ! If buffer can hold requested bytes, return
      if (size(Buf%Bytes) > Buf%NB + N) return

      ! Calculate new size
      NewSize = int(real(Buf%NB + N, R4Ki) * 1.8_R4Ki, IntKi)

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


   subroutine Pack_C1(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      character(*), intent(in)                :: Data
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = len(Data)

      ! Grow buffer to accommodate Data
      call GrowBuffer(Buf, DataSize)
      if (RegCheckErr(Buf, "Pack_C1")) return

      ! Transfer data to buffer
      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(Data, Buf%Bytes)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Unpack_C1(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      character(*), intent(out)               :: Data
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = len(Data)

      ! Check that buffer has sufficient bytes remaining
      if (size(Buf%Bytes) < Buf%NB + DataSize) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) "Unpack_C1: buffer too small, requested", DataSize, "bytes"
         return
      end if

      ! Transfer data from buffer
      Data = transfer(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize), Data)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Pack_C1_Rank1(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      character(*), intent(in)                :: Data(:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = len(Data(1))*size(Data)

      ! Grow buffer to accommodate Data
      call GrowBuffer(Buf, DataSize)
      if (RegCheckErr(Buf, "Pack_C1_Rank1")) return

      ! Transfer data to buffer
      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(Data, Buf%Bytes)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Unpack_C1_Rank1(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      character(*), intent(out)               :: Data(:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = len(Data(1))*size(Data)

      ! Check that buffer has sufficient bytes remaining
      if (size(Buf%Bytes) < Buf%NB + DataSize) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) "Unpack_C1_Rank1: buffer too small, requested", DataSize, "bytes"
         return
      end if

      ! Transfer data from buffer
      Data = reshape(transfer(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize), Data), shape(Data))
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Pack_C1_Rank2(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      character(*), intent(in)                :: Data(:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = len(Data(1,1))*size(Data)

      ! Grow buffer to accommodate Data
      call GrowBuffer(Buf, DataSize)
      if (RegCheckErr(Buf, "Pack_C1_Rank2")) return

      ! Transfer data to buffer
      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(Data, Buf%Bytes)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Unpack_C1_Rank2(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      character(*), intent(out)               :: Data(:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = len(Data(1,1))*size(Data)

      ! Check that buffer has sufficient bytes remaining
      if (size(Buf%Bytes) < Buf%NB + DataSize) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) "Unpack_C1_Rank2: buffer too small, requested", DataSize, "bytes"
         return
      end if

      ! Transfer data from buffer
      Data = reshape(transfer(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize), Data), shape(Data))
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Pack_C1_Rank3(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      character(*), intent(in)                :: Data(:,:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = len(Data(1,1,1))*size(Data)

      ! Grow buffer to accommodate Data
      call GrowBuffer(Buf, DataSize)
      if (RegCheckErr(Buf, "Pack_C1_Rank3")) return

      ! Transfer data to buffer
      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(Data, Buf%Bytes)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Unpack_C1_Rank3(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      character(*), intent(out)               :: Data(:,:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = len(Data(1,1,1))*size(Data)

      ! Check that buffer has sufficient bytes remaining
      if (size(Buf%Bytes) < Buf%NB + DataSize) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) "Unpack_C1_Rank3: buffer too small, requested", DataSize, "bytes"
         return
      end if

      ! Transfer data from buffer
      Data = reshape(transfer(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize), Data), shape(Data))
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Pack_C1_Rank4(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      character(*), intent(in)                :: Data(:,:,:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = len(Data(1,1,1,1))*size(Data)

      ! Grow buffer to accommodate Data
      call GrowBuffer(Buf, DataSize)
      if (RegCheckErr(Buf, "Pack_C1_Rank4")) return

      ! Transfer data to buffer
      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(Data, Buf%Bytes)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Unpack_C1_Rank4(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      character(*), intent(out)               :: Data(:,:,:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = len(Data(1,1,1,1))*size(Data)

      ! Check that buffer has sufficient bytes remaining
      if (size(Buf%Bytes) < Buf%NB + DataSize) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) "Unpack_C1_Rank4: buffer too small, requested", DataSize, "bytes"
         return
      end if

      ! Transfer data from buffer
      Data = reshape(transfer(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize), Data), shape(Data))
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Pack_C1_Rank5(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      character(*), intent(in)                :: Data(:,:,:,:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = len(Data(1,1,1,1,1))*size(Data)

      ! Grow buffer to accommodate Data
      call GrowBuffer(Buf, DataSize)
      if (RegCheckErr(Buf, "Pack_C1_Rank5")) return

      ! Transfer data to buffer
      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(Data, Buf%Bytes)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Unpack_C1_Rank5(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      character(*), intent(out)               :: Data(:,:,:,:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = len(Data(1,1,1,1,1))*size(Data)

      ! Check that buffer has sufficient bytes remaining
      if (size(Buf%Bytes) < Buf%NB + DataSize) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) "Unpack_C1_Rank5: buffer too small, requested", DataSize, "bytes"
         return
      end if

      ! Transfer data from buffer
      Data = reshape(transfer(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize), Data), shape(Data))
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Pack_L1(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      logical, intent(in)                     :: Data
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 1

      ! Grow buffer to accommodate Data
      call GrowBuffer(Buf, DataSize)
      if (RegCheckErr(Buf, "Pack_L1")) return

      ! Transfer data to buffer
      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(LogicalToByte(Data), Buf%Bytes)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Unpack_L1(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      logical, intent(out)                    :: Data
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 1

      ! Check that buffer has sufficient bytes remaining
      if (size(Buf%Bytes) < Buf%NB + DataSize) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) "Unpack_L1: buffer too small, requested", DataSize, "bytes"
         return
      end if

      ! Transfer data from buffer
      Data = ByteToLogical(Buf%Bytes(Buf%NB+1))
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Pack_L1_Rank1(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      logical, intent(in)                     :: Data(:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = size(Data)

      ! Grow buffer to accommodate Data
      call GrowBuffer(Buf, DataSize)
      if (RegCheckErr(Buf, "Pack_L1_Rank1")) return

      ! Transfer data to buffer
      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(LogicalToByte(Data), Buf%Bytes)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Unpack_L1_Rank1(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      logical, intent(out)                    :: Data(:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = size(Data)

      ! Check that buffer has sufficient bytes remaining
      if (size(Buf%Bytes) < Buf%NB + DataSize) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) "Unpack_L1_Rank1: buffer too small, requested", DataSize, "bytes"
         return
      end if

      ! Transfer data from buffer
      Data = reshape(ByteToLogical(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize)), shape(Data))
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Pack_L1_Rank2(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      logical, intent(in)                     :: Data(:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = size(Data)

      ! Grow buffer to accommodate Data
      call GrowBuffer(Buf, DataSize)
      if (RegCheckErr(Buf, "Pack_L1_Rank2")) return

      ! Transfer data to buffer
      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(LogicalToByte(Data), Buf%Bytes)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Unpack_L1_Rank2(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      logical, intent(out)                    :: Data(:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = size(Data)

      ! Check that buffer has sufficient bytes remaining
      if (size(Buf%Bytes) < Buf%NB + DataSize) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) "Unpack_L1_Rank2: buffer too small, requested", DataSize, "bytes"
         return
      end if

      ! Transfer data from buffer
      Data = reshape(ByteToLogical(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize)), shape(Data))
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Pack_L1_Rank3(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      logical, intent(in)                     :: Data(:,:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = size(Data)

      ! Grow buffer to accommodate Data
      call GrowBuffer(Buf, DataSize)
      if (RegCheckErr(Buf, "Pack_L1_Rank3")) return

      ! Transfer data to buffer
      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(LogicalToByte(Data), Buf%Bytes)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Unpack_L1_Rank3(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      logical, intent(out)                    :: Data(:,:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = size(Data)

      ! Check that buffer has sufficient bytes remaining
      if (size(Buf%Bytes) < Buf%NB + DataSize) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) "Unpack_L1_Rank3: buffer too small, requested", DataSize, "bytes"
         return
      end if

      ! Transfer data from buffer
      Data = reshape(ByteToLogical(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize)), shape(Data))
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Pack_L1_Rank4(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      logical, intent(in)                     :: Data(:,:,:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = size(Data)

      ! Grow buffer to accommodate Data
      call GrowBuffer(Buf, DataSize)
      if (RegCheckErr(Buf, "Pack_L1_Rank4")) return

      ! Transfer data to buffer
      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(LogicalToByte(Data), Buf%Bytes)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Unpack_L1_Rank4(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      logical, intent(out)                    :: Data(:,:,:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = size(Data)

      ! Check that buffer has sufficient bytes remaining
      if (size(Buf%Bytes) < Buf%NB + DataSize) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) "Unpack_L1_Rank4: buffer too small, requested", DataSize, "bytes"
         return
      end if

      ! Transfer data from buffer
      Data = reshape(ByteToLogical(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize)), shape(Data))
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Pack_L1_Rank5(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      logical, intent(in)                     :: Data(:,:,:,:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = size(Data)

      ! Grow buffer to accommodate Data
      call GrowBuffer(Buf, DataSize)
      if (RegCheckErr(Buf, "Pack_L1_Rank5")) return

      ! Transfer data to buffer
      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(LogicalToByte(Data), Buf%Bytes)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Unpack_L1_Rank5(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      logical, intent(out)                    :: Data(:,:,:,:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = size(Data)

      ! Check that buffer has sufficient bytes remaining
      if (size(Buf%Bytes) < Buf%NB + DataSize) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) "Unpack_L1_Rank5: buffer too small, requested", DataSize, "bytes"
         return
      end if

      ! Transfer data from buffer
      Data = reshape(ByteToLogical(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize)), shape(Data))
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Pack_I4(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      integer(B4Ki), intent(in)               :: Data
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 4

      ! Grow buffer to accommodate Data
      call GrowBuffer(Buf, DataSize)
      if (RegCheckErr(Buf, "Pack_I4")) return

      ! Transfer data to buffer
      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(Data, Buf%Bytes)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Unpack_I4(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      integer(B4Ki), intent(out)              :: Data
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 4

      ! Check that buffer has sufficient bytes remaining
      if (size(Buf%Bytes) < Buf%NB + DataSize) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) "Unpack_I4: buffer too small, requested", DataSize, "bytes"
         return
      end if

      ! Transfer data from buffer
      Data = transfer(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize), Data)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Pack_I4_Rank1(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      integer(B4Ki), intent(in)               :: Data(:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 4*size(Data)

      ! Grow buffer to accommodate Data
      call GrowBuffer(Buf, DataSize)
      if (RegCheckErr(Buf, "Pack_I4_Rank1")) return

      ! Transfer data to buffer
      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(Data, Buf%Bytes)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Unpack_I4_Rank1(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      integer(B4Ki), intent(out)              :: Data(:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 4*size(Data)

      ! Check that buffer has sufficient bytes remaining
      if (size(Buf%Bytes) < Buf%NB + DataSize) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) "Unpack_I4_Rank1: buffer too small, requested", DataSize, "bytes"
         return
      end if

      ! Transfer data from buffer
      Data = reshape(transfer(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize), Data), shape(Data))
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Pack_I4_Rank2(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      integer(B4Ki), intent(in)               :: Data(:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 4*size(Data)

      ! Grow buffer to accommodate Data
      call GrowBuffer(Buf, DataSize)
      if (RegCheckErr(Buf, "Pack_I4_Rank2")) return

      ! Transfer data to buffer
      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(Data, Buf%Bytes)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Unpack_I4_Rank2(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      integer(B4Ki), intent(out)              :: Data(:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 4*size(Data)

      ! Check that buffer has sufficient bytes remaining
      if (size(Buf%Bytes) < Buf%NB + DataSize) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) "Unpack_I4_Rank2: buffer too small, requested", DataSize, "bytes"
         return
      end if

      ! Transfer data from buffer
      Data = reshape(transfer(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize), Data), shape(Data))
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Pack_I4_Rank3(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      integer(B4Ki), intent(in)               :: Data(:,:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 4*size(Data)

      ! Grow buffer to accommodate Data
      call GrowBuffer(Buf, DataSize)
      if (RegCheckErr(Buf, "Pack_I4_Rank3")) return

      ! Transfer data to buffer
      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(Data, Buf%Bytes)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Unpack_I4_Rank3(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      integer(B4Ki), intent(out)              :: Data(:,:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 4*size(Data)

      ! Check that buffer has sufficient bytes remaining
      if (size(Buf%Bytes) < Buf%NB + DataSize) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) "Unpack_I4_Rank3: buffer too small, requested", DataSize, "bytes"
         return
      end if

      ! Transfer data from buffer
      Data = reshape(transfer(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize), Data), shape(Data))
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Pack_I4_Rank4(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      integer(B4Ki), intent(in)               :: Data(:,:,:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 4*size(Data)

      ! Grow buffer to accommodate Data
      call GrowBuffer(Buf, DataSize)
      if (RegCheckErr(Buf, "Pack_I4_Rank4")) return

      ! Transfer data to buffer
      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(Data, Buf%Bytes)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Unpack_I4_Rank4(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      integer(B4Ki), intent(out)              :: Data(:,:,:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 4*size(Data)

      ! Check that buffer has sufficient bytes remaining
      if (size(Buf%Bytes) < Buf%NB + DataSize) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) "Unpack_I4_Rank4: buffer too small, requested", DataSize, "bytes"
         return
      end if

      ! Transfer data from buffer
      Data = reshape(transfer(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize), Data), shape(Data))
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Pack_I4_Rank5(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      integer(B4Ki), intent(in)               :: Data(:,:,:,:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 4*size(Data)

      ! Grow buffer to accommodate Data
      call GrowBuffer(Buf, DataSize)
      if (RegCheckErr(Buf, "Pack_I4_Rank5")) return

      ! Transfer data to buffer
      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(Data, Buf%Bytes)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Unpack_I4_Rank5(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      integer(B4Ki), intent(out)              :: Data(:,:,:,:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 4*size(Data)

      ! Check that buffer has sufficient bytes remaining
      if (size(Buf%Bytes) < Buf%NB + DataSize) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) "Unpack_I4_Rank5: buffer too small, requested", DataSize, "bytes"
         return
      end if

      ! Transfer data from buffer
      Data = reshape(transfer(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize), Data), shape(Data))
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Pack_R4(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      real(R4Ki), intent(in)                  :: Data
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 4

      ! Grow buffer to accommodate Data
      call GrowBuffer(Buf, DataSize)
      if (RegCheckErr(Buf, "Pack_R4")) return

      ! Transfer data to buffer
      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(Data, Buf%Bytes)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Unpack_R4(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      real(R4Ki), intent(out)                 :: Data
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 4

      ! Check that buffer has sufficient bytes remaining
      if (size(Buf%Bytes) < Buf%NB + DataSize) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) "Unpack_R4: buffer too small, requested", DataSize, "bytes"
         return
      end if

      ! Transfer data from buffer
      Data = transfer(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize), Data)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Pack_R4_Rank1(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      real(R4Ki), intent(in)                  :: Data(:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 4*size(Data)

      ! Grow buffer to accommodate Data
      call GrowBuffer(Buf, DataSize)
      if (RegCheckErr(Buf, "Pack_R4_Rank1")) return

      ! Transfer data to buffer
      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(Data, Buf%Bytes)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Unpack_R4_Rank1(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      real(R4Ki), intent(out)                 :: Data(:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 4*size(Data)

      ! Check that buffer has sufficient bytes remaining
      if (size(Buf%Bytes) < Buf%NB + DataSize) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) "Unpack_R4_Rank1: buffer too small, requested", DataSize, "bytes"
         return
      end if

      ! Transfer data from buffer
      Data = reshape(transfer(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize), Data), shape(Data))
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Pack_R4_Rank2(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      real(R4Ki), intent(in)                  :: Data(:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 4*size(Data)

      ! Grow buffer to accommodate Data
      call GrowBuffer(Buf, DataSize)
      if (RegCheckErr(Buf, "Pack_R4_Rank2")) return

      ! Transfer data to buffer
      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(Data, Buf%Bytes)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Unpack_R4_Rank2(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      real(R4Ki), intent(out)                 :: Data(:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 4*size(Data)

      ! Check that buffer has sufficient bytes remaining
      if (size(Buf%Bytes) < Buf%NB + DataSize) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) "Unpack_R4_Rank2: buffer too small, requested", DataSize, "bytes"
         return
      end if

      ! Transfer data from buffer
      Data = reshape(transfer(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize), Data), shape(Data))
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Pack_R4_Rank3(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      real(R4Ki), intent(in)                  :: Data(:,:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 4*size(Data)

      ! Grow buffer to accommodate Data
      call GrowBuffer(Buf, DataSize)
      if (RegCheckErr(Buf, "Pack_R4_Rank3")) return

      ! Transfer data to buffer
      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(Data, Buf%Bytes)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Unpack_R4_Rank3(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      real(R4Ki), intent(out)                 :: Data(:,:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 4*size(Data)

      ! Check that buffer has sufficient bytes remaining
      if (size(Buf%Bytes) < Buf%NB + DataSize) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) "Unpack_R4_Rank3: buffer too small, requested", DataSize, "bytes"
         return
      end if

      ! Transfer data from buffer
      Data = reshape(transfer(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize), Data), shape(Data))
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Pack_R4_Rank4(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      real(R4Ki), intent(in)                  :: Data(:,:,:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 4*size(Data)

      ! Grow buffer to accommodate Data
      call GrowBuffer(Buf, DataSize)
      if (RegCheckErr(Buf, "Pack_R4_Rank4")) return

      ! Transfer data to buffer
      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(Data, Buf%Bytes)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Unpack_R4_Rank4(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      real(R4Ki), intent(out)                 :: Data(:,:,:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 4*size(Data)

      ! Check that buffer has sufficient bytes remaining
      if (size(Buf%Bytes) < Buf%NB + DataSize) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) "Unpack_R4_Rank4: buffer too small, requested", DataSize, "bytes"
         return
      end if

      ! Transfer data from buffer
      Data = reshape(transfer(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize), Data), shape(Data))
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Pack_R4_Rank5(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      real(R4Ki), intent(in)                  :: Data(:,:,:,:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 4*size(Data)

      ! Grow buffer to accommodate Data
      call GrowBuffer(Buf, DataSize)
      if (RegCheckErr(Buf, "Pack_R4_Rank5")) return

      ! Transfer data to buffer
      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(Data, Buf%Bytes)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Unpack_R4_Rank5(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      real(R4Ki), intent(out)                 :: Data(:,:,:,:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 4*size(Data)

      ! Check that buffer has sufficient bytes remaining
      if (size(Buf%Bytes) < Buf%NB + DataSize) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) "Unpack_R4_Rank5: buffer too small, requested", DataSize, "bytes"
         return
      end if

      ! Transfer data from buffer
      Data = reshape(transfer(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize), Data), shape(Data))
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Pack_R8(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      real(R8Ki), intent(in)                  :: Data
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 8

      ! Grow buffer to accommodate Data
      call GrowBuffer(Buf, DataSize)
      if (RegCheckErr(Buf, "Pack_R8")) return

      ! Transfer data to buffer
      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(Data, Buf%Bytes)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Unpack_R8(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      real(R8Ki), intent(out)                 :: Data
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 8

      ! Check that buffer has sufficient bytes remaining
      if (size(Buf%Bytes) < Buf%NB + DataSize) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) "Unpack_R8: buffer too small, requested", DataSize, "bytes"
         return
      end if

      ! Transfer data from buffer
      Data = transfer(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize), Data)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Pack_R8_Rank1(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      real(R8Ki), intent(in)                  :: Data(:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 8*size(Data)

      ! Grow buffer to accommodate Data
      call GrowBuffer(Buf, DataSize)
      if (RegCheckErr(Buf, "Pack_R8_Rank1")) return

      ! Transfer data to buffer
      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(Data, Buf%Bytes)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Unpack_R8_Rank1(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      real(R8Ki), intent(out)                 :: Data(:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 8*size(Data)

      ! Check that buffer has sufficient bytes remaining
      if (size(Buf%Bytes) < Buf%NB + DataSize) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) "Unpack_R8_Rank1: buffer too small, requested", DataSize, "bytes"
         return
      end if

      ! Transfer data from buffer
      Data = reshape(transfer(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize), Data), shape(Data))
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Pack_R8_Rank2(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      real(R8Ki), intent(in)                  :: Data(:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 8*size(Data)

      ! Grow buffer to accommodate Data
      call GrowBuffer(Buf, DataSize)
      if (RegCheckErr(Buf, "Pack_R8_Rank2")) return

      ! Transfer data to buffer
      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(Data, Buf%Bytes)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Unpack_R8_Rank2(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      real(R8Ki), intent(out)                 :: Data(:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 8*size(Data)

      ! Check that buffer has sufficient bytes remaining
      if (size(Buf%Bytes) < Buf%NB + DataSize) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) "Unpack_R8_Rank2: buffer too small, requested", DataSize, "bytes"
         return
      end if

      ! Transfer data from buffer
      Data = reshape(transfer(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize), Data), shape(Data))
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Pack_R8_Rank3(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      real(R8Ki), intent(in)                  :: Data(:,:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 8*size(Data)

      ! Grow buffer to accommodate Data
      call GrowBuffer(Buf, DataSize)
      if (RegCheckErr(Buf, "Pack_R8_Rank3")) return

      ! Transfer data to buffer
      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(Data, Buf%Bytes)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Unpack_R8_Rank3(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      real(R8Ki), intent(out)                 :: Data(:,:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 8*size(Data)

      ! Check that buffer has sufficient bytes remaining
      if (size(Buf%Bytes) < Buf%NB + DataSize) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) "Unpack_R8_Rank3: buffer too small, requested", DataSize, "bytes"
         return
      end if

      ! Transfer data from buffer
      Data = reshape(transfer(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize), Data), shape(Data))
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Pack_R8_Rank4(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      real(R8Ki), intent(in)                  :: Data(:,:,:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 8*size(Data)

      ! Grow buffer to accommodate Data
      call GrowBuffer(Buf, DataSize)
      if (RegCheckErr(Buf, "Pack_R8_Rank4")) return

      ! Transfer data to buffer
      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(Data, Buf%Bytes)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Unpack_R8_Rank4(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      real(R8Ki), intent(out)                 :: Data(:,:,:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 8*size(Data)

      ! Check that buffer has sufficient bytes remaining
      if (size(Buf%Bytes) < Buf%NB + DataSize) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) "Unpack_R8_Rank4: buffer too small, requested", DataSize, "bytes"
         return
      end if

      ! Transfer data from buffer
      Data = reshape(transfer(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize), Data), shape(Data))
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Pack_R8_Rank5(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      real(R8Ki), intent(in)                  :: Data(:,:,:,:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 8*size(Data)

      ! Grow buffer to accommodate Data
      call GrowBuffer(Buf, DataSize)
      if (RegCheckErr(Buf, "Pack_R8_Rank5")) return

      ! Transfer data to buffer
      Buf%Bytes(Buf%NB+1:Buf%NB+DataSize) = transfer(Data, Buf%Bytes)
      Buf%NB = Buf%NB + DataSize

   end subroutine

   subroutine Unpack_R8_Rank5(Buf, Data)
      type(PackBuffer), intent(inout)         :: Buf
      real(R8Ki), intent(out)                 :: Data(:,:,:,:,:)
      integer(IntKi)                          :: DataSize

      ! If buffer error, return
      if (Buf%ErrStat /= ErrID_None) return

      ! Get size of data in bytes
      DataSize = 8*size(Data)

      ! Check that buffer has sufficient bytes remaining
      if (size(Buf%Bytes) < Buf%NB + DataSize) then
         Buf%ErrStat = ErrID_Fatal
         write(Buf%ErrMsg,*) "Unpack_R8_Rank5: buffer too small, requested", DataSize, "bytes"
         return
      end if

      ! Transfer data from buffer
      Data = reshape(transfer(Buf%Bytes(Buf%NB+1:Buf%NB+DataSize), Data), shape(Data))
      Buf%NB = Buf%NB + DataSize

   end subroutine
end module