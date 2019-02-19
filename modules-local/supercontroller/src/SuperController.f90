!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015  National Renewable Energy Laboratory
!
!    Lidar module, a submodule of InflowWind
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
MODULE SuperController

! This is a pseudo module used to couple FAST v8 with SuperController; it is considered part of the FAST glue code
   USE FAST_Types

   IMPLICIT NONE

   PRIVATE

   TYPE(ProgDesc), PARAMETER            :: SC_Ver = ProgDesc( 'SuperController Integration', '', '' )

  
! ==================================================================================================="

      
      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: Init_SC                           ! Initialization routine
   PUBLIC :: SC_SetInputs                      ! Glue-code routine to update inputs for SuperController
   PUBLIC :: SC_SetOutputs                     ! Glue-code routine to update inputs to turbine controller from SuperController
   
   
CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Init_SC( InitInp, SC, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(SC_InitInputType),        INTENT(IN   )  :: InitInp     ! Input data for initialization routine
   TYPE(SuperController_Data),             INTENT(INOUT)  :: SC        ! data for the SuperController integration module
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
                                          
      ! local variables                   
   INTEGER(IntKi)                                   :: ErrStat2    ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                             :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None
                                                    
   CHARACTER(*),   PARAMETER                        :: RoutineName = 'Init_SC'
                                          
      ! Initialize variables

   ErrStat = ErrID_None
   ErrMsg  = ""

   IF (InitInp%NumCtrl2SC > 0) THEN
      CALL AllocPAry( SC%u%toSC, InitInp%NumCtrl2SC, 'u%toSC', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF
   
   IF (ErrStat >= AbortErrLev) RETURN
   
      ! make sure the C versions are synced with these arrays
   if (InitInp%NumCtrl2SC > 0) then
      SC%u%c_obj%toSC_Len = InitInp%NumCtrl2SC
      SC%u%c_obj%toSC     = C_LOC( SC%u%toSC(1) )
      SC%u%toSC = 0.0_ReKi
   end if
      
      ! initialize the arrays:
      ! 
      
      !............................................................................................
      ! Define system output initializations (set up mesh) here:
      !............................................................................................   
   if (InitInp%NumSC2Ctrl > 0) then
      CALL AllocPAry( SC%y%fromSC, InitInp%NumSC2Ctrl, 'y%fromSC', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   end if
   
   IF (ErrStat >= AbortErrLev) RETURN
                        
      ! make sure the C versions are synced with these arrays
   if (InitInp%NumSC2Ctrl > 0) then
      SC%y%c_obj%fromSC_Len = InitInp%NumSC2Ctrl
      SC%y%c_obj%fromSC     = C_LOC( SC%y%fromSC(1) )
   end if
   

   if( (InitInp%NumSC2Ctrl > 0) .and. (InitInp%NumSC2Ctrl > 0)) then
      SC%p%scOn = .true. 
   else
      SC%p%scOn = .false.
   end if

   RETURN
   
END SUBROUTINE Init_SC
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SC_SetInputs(p_FAST, y_SrvD, SC, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(FAST_ParameterType),       INTENT(IN    )  :: p_FAST      ! Parameters for the glue code
   TYPE(SrvD_OutputType),          INTENT(IN)      :: y_SrvD      ! The outputs of the ServoDyn module (control)
   TYPE(SuperController_Data),            INTENT(INOUT)   :: SC        ! data for the SuperController integration module
   INTEGER(IntKi),                 INTENT(  OUT)   :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)   :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables
   INTEGER(IntKi)                                  :: ErrStat2    ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                            :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None
                                                   
   CHARACTER(*),   PARAMETER                       :: RoutineName = 'SC_SetInputs'
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
      ! set SuperController inputs
   if (p_FAST%CompServo == Module_SrvD) then
      if (allocated(y_SrvD%SuperController).and. associated(SC%u%toSC)) SC%u%toSC = y_SrvD%SuperController      
   end if
   
      
END SUBROUTINE SC_SetInputs
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SC_SetOutputs(p_FAST, u_SrvD, SC, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(FAST_ParameterType),       INTENT(IN    )  :: p_FAST      ! Parameters for the glue code
   TYPE(SrvD_InputType),          INTENT(INOUT)      :: u_SrvD      ! The inputs of the ServoDyn module (control)
   TYPE(SuperController_Data),            INTENT(IN)   :: SC        ! data for the SuperController integration module
   INTEGER(IntKi),                 INTENT(  OUT)   :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)   :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables
   INTEGER(IntKi)                                  :: ErrStat2    ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                            :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None
                                                   
   CHARACTER(*),   PARAMETER                       :: RoutineName = 'SC_SetOutputs'
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
      ! set SuperController inputs
   if (p_FAST%CompServo == Module_SrvD) then
      if (allocated(u_SrvD%SuperController).and. associated(SC%y%fromSC)) u_SrvD%SuperController = SC%y%fromSC
   end if
   
      
END SUBROUTINE SC_SetOutputs
!----------------------------------------------------------------------------------------------------------------------------------
END MODULE SuperController
!**********************************************************************************************************************************




