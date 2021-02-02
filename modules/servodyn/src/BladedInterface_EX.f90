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
   integer(IntKi), parameter  :: EXavrSWAP_Size = 3300     !< size of the avrSWAP array with the extended array sizing (increment this as new blocks ar added)
   real(ReKi),     parameter  :: EXavrSWAP_Ver  = 1.000    !< Version. increment minor for new signal addition, increment major for new block addition 


CONTAINS


!==================================================================================================================================
!> This routine fills the extended avrSWAP
SUBROUTINE Fill_EXavrSWAP( t, u, p, dll_data )
   REAL(DbKi),                     INTENT(IN   )  :: t           !< Current simulation time in seconds
   TYPE(SrvD_InputType),           INTENT(IN   )  :: u           !< Inputs at t
   TYPE(SrvD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(BladedDLLType),            INTENT(INOUT)  :: dll_data    !< data for the Bladed DLL

      ! local variables:
   INTEGER(IntKi)                                 :: I           ! Loop counter

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
   TYPE(SrvD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(BladedDLLType),            INTENT(INOUT)  :: dll_data    !< data for the Bladed DLL
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! local variables:
   INTEGER(IntKi)                                 :: K           ! Loop counter
   CHARACTER(*), PARAMETER                        :: RoutineName = 'Retrieve_EXavrSWAP'
   
      
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
