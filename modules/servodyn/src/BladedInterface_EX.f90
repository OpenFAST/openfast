!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013-2016  National Renewable Energy Laboratory
!
!    This file is part of FAST's Controls and Electrical Drive Module, "ServoDyn".
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
!
!**********************************************************************************************************************************
MODULE BladedInterface_EX

   USE NWTC_Library  
   USE ServoDyn_Types

   IMPLICIT                        NONE

   !> Extended AVRswap array --  This is an oversizing of the bladed 4.3 avrswap array to add an additional ~3000 communication channels
   !!
   !!    This is intended as a stop-gap to get a much larger set of channels integrated.  This is intended as a short term solution
   !!    until an entirely different architecture for the controller is adopted.  But given how many things get built as temporary
   !!    and never updated to a better solution, we are adding some flexibility to it (a temp-o-permanent solution).
   !!
   !!    Versioning:
   !!          It could be useful for the controller to know which version of the interface it is dealing with.  For example, if
   !!          the controller is expecting a block listed below, it would be useful for it to check version info for compatibility.
   !!          The controller could still talk to newer interfaces since all channels are set, but not to older interfaces (like
   !!          when a new block is added here, but the controller is of an older protocol, it should still work ok).
   !!       Major:   increment version to next integer when a block is added
   !!       Minor:   increment version by 0.001 for a single channel addition within an existing block
   !!       Sub-Major change: if the whole organization of the blocks is changed, then increment to next 0.010 integer multiple.
   !!
   !!    Extended array layout:
   !!
   !!       |  Index start |  Index end   | Channels in group  |  Description of array block                |  OpenFAST module      |
   !!       | :----------: | :----------: | :----------------: | :----------------------------------------: | :-------------------: |
   !!       |     1000     |              |           1        |  Version info for compatibility checks     |                       |
   !!       |     1001     |     2000     |        1000        |  Sensor inputs (non-lidar)                 |                       |
   !!       |     2001     |     2500     |         500        |  Lidar input measurements                  |  LidarSim             |
   !!       |     2501     |     2600     |         100        |  Lidar control signals                     |  LidarSim             |
   !!       |     2601     |     2800     |         200        |  Cable controls                            |  MoorDyn, SubDyn      |
   !!       |     2801     |     3000     |         200        |  Tuned Mass dampers control                |  StrucCtrl            |
   !!       |     3001     |     3300     |         300        |  Aero controls (individual flaps, etc)     |  AeroDyn              |
   !!       |              |              |                    |                                            |                       |
   !!
   integer(IntKi),   parameter   :: EXavrSWAP_Size = 3300      !< size of the avrSWAP array with the extended array sizing (increment this as new blocks ar added)
   real(ReKi),       parameter   :: EXavrSWAP_Ver  = 1.000     !< Version. increment minor for new signal addition, increment major for new block addition
   integer(IntKi),   parameter   :: CableCtrl_StartIdx = 2601  !< Starting index for the cable control
   integer(IntKi),   parameter   :: CableCtrl_MaxChan = 200    !< Maximum channels in cable control group


CONTAINS
!==================================================================================================================================
!> This routine sets the input and output necessary for the extended interface.
SUBROUTINE EXavrSWAP_Init( InitInp, u, p, y, dll_data, UnSum, ErrStat, ErrMsg)
   type(SrvD_InitInputType),     intent(in   )  :: InitInp        !< Input data for initialization routine
   type(SrvD_InputType),         intent(in   )  :: u              !< Inputs at t
   type(SrvD_ParameterType),     intent(in   )  :: p              !< Parameters
   type(SrvD_OutputType),        intent(inout)  :: y              !< Initial system outputs (outputs are not calculated)
   type(BladedDLLType),          intent(inout)  :: dll_data       !< data for the Bladed DLL
   integer(IntKi),               intent(in   )  :: UnSum          !< Unit number for summary file (>0 when active)
   integer(IntKi),               intent(  out)  :: ErrStat        !< Error status of the operation
   character(ErrMsgLen),         intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None
   character(*),                 parameter      :: RoutineName='EXavrSWAP_Init'
   character(1024),              allocatable    :: SumInfo(:)     ! Description strings for each avrSWAP record -- only used for summary file writing
   character(3),                 allocatable    :: DataFlow(:)    ! Direction of data flow -- only used for summary file writing
   character(64),                allocatable    :: Requestor(:)   ! Info on module requesting the channel
   integer(IntKi)                               :: ErrStat2
   character(ErrMsgLen)                         :: ErrMsg2

      ! Initialize ErrStat and ErrMsg
   ErrStat = ErrID_None
   ErrMsg  = ''



      ! Make sure we didn't get here by mistake
   if ( .not. p%EXavrSWAP ) return

   if (UnSum>0) then
      call AllocAry(SumInfo,size(dll_data%avrSwap),'SumInfo array for bladed interface',ErrStat2,ErrMsg2)
         if (Failed())  return
      SumInfo  = ''
      call AllocAry(Requestor,size(dll_data%avrSwap),'Requestor array for bladed interface',ErrStat2,ErrMsg2)
         if (Failed())  return
      Requestor  = ''
      call AllocAry(DataFlow,size(dll_data%avrSwap),'DataFlow array for bladed interface',ErrStat2,ErrMsg2)
         if (Failed())  return
      DataFlow  = ''
   endif


      ! Set the input and output arrays


      ! Initialize cable controls
   if (InitInp%NumCableControl > 0) then
      call InitCableCtrl()
        if (Failed())  return
   endif

      ! Write to summary file
   if (UnSum>0)   call WrBladedSumInfoToFile()

      ! Cleanup afterwards
   call Cleanup

contains
   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   end function Failed
   subroutine CleanUp()
      if (allocated(SumInfo))    deallocate(SumInfo)
      if (allocated(Requestor))  deallocate(Requestor)
      if (allocated(DataFlow))   deallocate(DataFlow)
   end subroutine CleanUp

   subroutine InitCableCtrl()
      integer(IntKi)    :: I,J   ! Generic counters
      if (.not. allocated(InitInp%CableControlRequestor)) then
         ErrStat2=ErrID_Fatal
         ErrMsg2='Cable control string array indicating which module requested cable controls is missing (CableControlRequestor)'
         if (Failed())  return
      endif
      if (size(InitInp%CableControlRequestor) /= InitInp%NumCableControl) then
         ErrStat2=ErrID_Fatal
         ErrMsg2='Size of cable control string array (CableControlRequestor) does not match the number of requested cable control channels.'
         if (Failed())  return
      endif
      if ( InitInp%NumCableControl*2_IntKi > CableCtrl_MaxChan ) then
         ErrStat2=ErrID_Fatal
         ErrMsg2='Maximum number of cable control channels exceeded:  requested '//trim(Num2LStr(InitInp%NumCableControl))// &
                  ' channel sets ('//trim(Num2LStr(InitInp%NumCableControl*2_IntKi))//' individual channels),'// &
                  ' but only '//trim(Num2LStr(CableCtrl_MaxChan))//' individual channels are available'
      endif
      call AllocAry( y%CableDeltaL, InitInp%NumCableControl, 'CableDeltaL', ErrStat2, ErrMsg2 )
         if (Failed())  return
      call AllocAry( y%CableDeltaLdot, InitInp%NumCableControl, 'CableDeltaLdot', ErrStat2, ErrMsg2 )
         if (Failed())  return
      ! Create info for summary file about channels
      if (UnSum > 0) then
         do I=1,InitInp%NumCableControl
            J=CableCtrl_StartIdx + (I*2)-1      ! Index into the full avrSWAP
            call WrSumInfoRcvd(J,  InitInp%CableControlRequestor(I),'Cable control channel group '//trim(Num2LStr(I))//' -- DeltaL')
            call WrSumInfoRcvd(J+1,InitInp%CableControlRequestor(I),'Cable control channel group '//trim(Num2LStr(I))//' -- DeltaLdot')
         enddo
      endif
   end subroutine InitCableCtrl

   subroutine WrBladedSumInfoToFile()
      integer(IntKi)    :: I  !< generic counter
      write(UnSum,'(A)') ''
      write(UnSum,'(A)') '   Legacy Bladed DLL interface with Extended avrSWAP'
      write(UnSum,'(A)') '          channel usage by SrvD:'
      write(UnSum,'(A)') ''
      write(UnSum,'(6x,8x,3x,A3,3x,A)') '-->','indicates from DLL to SrvD'
      write(UnSum,'(6x,8x,3x,A3,3x,A)') '-->','indicates from SrvD to DLL'
      write(UnSum,'(6x,A8,3x,A3,3x,A21,3x,A11)') 'Record #','   ','Requested by         ','Description'
      write(UnSum,'(6x,A8,3x,A3,3x,A21,3x,A11)') '--------','   ','---------------------','-----------'
      do I=1,size(SumInfo)
         if (len_trim(SumInfo(I)) > 0 )  write(UnSum,'(4x,I4,6x,A3,2x,A)')  I,DataFlow(I),trim(Requestor(I)),trim(SumInfo(I))
      enddo
   end subroutine WrBladedSumInfoToFile
   subroutine WrSumInfoSend(Record,Desc)
      integer(IntKi),   intent(in   )  :: Record
      character(*),     intent(in   )  :: Desc
      DataFlow(Record)  = '-->'
      SumInfo(Record)   = trim(Desc(1:min(len_trim(Desc),len(SumInfo(1)))))     ! prevent string length overrun
   end subroutine WrSumInfoSend
   subroutine WrSumInfoRcvd(Record,Rqst,Desc)
      integer(IntKi),   intent(in   )  :: Record
      character(*),     intent(in   )  :: Rqst
      character(*),     intent(in   )  :: Desc
      DataFlow(Record)  = '<--'
      Requestor(Record) = trim(Rqst(1:min(len_trim(Desc),len(Requestor(1)))))   ! prevent string length overrun
      SumInfo(Record)   = trim(Desc(1:min(len_trim(Desc),len(SumInfo(1)))))     ! prevent string length overrun
   end subroutine WrSumInfoRcvd
END SUBROUTINE

!==================================================================================================================================
!> This routine fills the extended avrSWAP
SUBROUTINE Fill_EXavrSWAP( t, u, p, dll_data )
   real(DbKi),                   intent(in   )  :: t           !< Current simulation time in seconds
   type(SrvD_InputType),         intent(in   )  :: u           !< Inputs at t
   type(SrvD_ParameterType),     intent(in   )  :: p           !< Parameters
   type(BladedDLLType),          intent(inout)  :: dll_data    !< data for the Bladed DLL

      ! local variables:
   integer(IntKi)                                 :: I           ! Loop counter

      ! Make sure we didn't get here by mistake
   if ( .not. p%EXavrSWAP ) return

      ! Set the version info on the extended avrswap if using
   dll_data%avrswap(1000) = real(EXavrSWAP_Ver,SiKi)

   call SetEXavrSWAP_Sensors()

   call SetEXavrSWAP_LidarSensors()
CONTAINS
   !> Set the sensor inputs for non-lidar channels
   !!    avrSWAP(1001:2000)
   subroutine SetEXavrSWAP_Sensors()
         ! in case something got set wrong, don't try to write beyond array
      if (size(dll_data%avrswap) < 2000 ) return
   end subroutine SetEXavrSWAP_Sensors

   !> Set the Lidar related sensor inputs
   !!    avrSWAP(2001:2500)
   subroutine SetEXavrSWAP_LidarSensors()
         ! in case something got set wrong, don't try to write beyond array
      if (size(dll_data%avrswap) < 2500 ) return
   end subroutine SetEXavrSWAP_LidarSensors

END SUBROUTINE Fill_EXavrSWAP


!==================================================================================================================================
!> This routine retrieves the extended avrSWAP
SUBROUTINE Retrieve_EXavrSWAP( p, dll_data, ErrStat, ErrMsg )
   type(SrvD_ParameterType),     intent(in   )  :: p           !< Parameters
   type(BladedDLLType),          intent(inout)  :: dll_data    !< data for the Bladed DLL
   integer(IntKi),               intent(  out)  :: ErrStat     !< Error status of the operation
   character(ErrMsgLen),         intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! local variables:
   integer(IntKi)                               :: K           ! Loop counter
   character(*),  parameter                     :: RoutineName = 'Retrieve_EXavrSWAP'
   
      ! Initialize ErrStat and ErrMsg
   ErrStat = ErrID_None
   ErrMsg  = ''   

   call Retrieve_EXavrSWAP_Lidar ()
   call Retrieve_EXavrSWAP_Cable ()
   call Retrieve_EXavrSWAP_StControls ()
   call Retrieve_EXavrSWAP_AeroControls ()
CONTAINS
   !> Controller signals to Lidar module
   !!    avrSWAP(2501:2600)
   subroutine Retrieve_EXavrSWAP_Lidar ()
      ! in case something got set wrong, don't try to read beyond array
      if (size(dll_data%avrswap) < 2000 ) return
   end subroutine Retrieve_EXavrSWAP_Lidar

   !> Controller signals to mooring systems (active mooring)
   !!    avrSWAP(2601:2800)
   subroutine Retrieve_EXavrSWAP_Cable ()
      ! in case something got set wrong, don't try to read beyond array
      if (size(dll_data%avrswap) < 2800 ) return
   end subroutine Retrieve_EXavrSWAP_Cable

   !> Controller signals to substructure controls (actuators in substructure)
   !!    avrSWAP(3001:3200)
   subroutine Retrieve_EXavrSWAP_StControls ()
      ! in case something got set wrong, don't try to read beyond array
      if (size(dll_data%avrswap) < 3200 ) return
   end subroutine Retrieve_EXavrSWAP_StControls

   !> Controller signals to active aero elements (aero actuators in blades)
   !!    avrSWAP(3201:3500)
   subroutine Retrieve_EXavrSWAP_AeroControls ()
      ! in case something got set wrong, don't try to read beyond array
      if (size(dll_data%avrswap) < 3500 ) return
   end subroutine Retrieve_EXavrSWAP_AeroControls

END SUBROUTINE Retrieve_EXavrSWAP

!==================================================================================================================================
END MODULE BladedInterface_EX
