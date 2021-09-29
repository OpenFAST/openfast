!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015  National Renewable Energy Laboratory
!
!    SuperController DataExchange, a submodule of openfast
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
! File last committed: $Date: $
! (File) Revision #: $Rev: $
! URL: $HeadURL: $
!**********************************************************************************************************************************
MODULE SC_DataEx

! This is a pseudo module used to couple FAST v8 with SuperController; it is considered part of the FAST glue code
   USE FAST_Types
   USE SCDataEx_Types

   IMPLICIT NONE

   PRIVATE

   TYPE(ProgDesc), PARAMETER            :: SC_DX_Ver = ProgDesc( 'SuperController DataExchange', '', '' )

  
! ==================================================================================================="

      
      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: SC_DX_Init                           ! Initialization routine
   PUBLIC :: SC_DX_SetInputs                      ! Glue-code routine to update inputs for SuperController
   PUBLIC :: SC_DX_SetOutputs                     ! Glue-code routine to update inputs to turbine controller from SuperController
   
   
CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SC_DX_Init( NumSC2CtrlGlob, NumSC2Ctrl, NumCtrl2SC, SC_DX, ErrStat, ErrMsg )
!..................................................................................................................................
   INTEGER(IntKi),                  INTENT(IN   )  :: NumSC2CtrlGlob
   INTEGER(IntKi),                  INTENT(IN   )  :: NumSC2Ctrl
   INTEGER(IntKi),                  INTENT(IN   )  :: NumCtrl2SC
   TYPE(SCDataEx_Data),             INTENT(INOUT)  :: SC_DX        ! data for the SuperController integration module
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
                                          
      ! local variables                   
   INTEGER(IntKi)                                  :: ErrStat2    ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                            :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None
                                                   
   CHARACTER(*),   PARAMETER                       :: RoutineName = 'SC_DX_Init'
                                          
      ! Initialize variables

   ErrStat = ErrID_None
   ErrMsg  = ""

   IF (NumCtrl2SC > 0) THEN
      CALL AllocPAry( SC_DX%u%toSC, NumCtrl2SC, 'u%toSC', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF
   
   IF (ErrStat >= AbortErrLev) RETURN
   
      ! make sure the C versions are synced with these arrays
   if (NumCtrl2SC > 0) then
      SC_DX%u%c_obj%toSC_Len = NumCtrl2SC
      SC_DX%u%c_obj%toSC     = C_LOC( SC_DX%u%toSC(1) )
   else
      SC_DX%u%c_obj%toSC_Len = 0
      SC_DX%u%c_obj%toSC     = C_NULL_PTR
   end if
      
      
      !............................................................................................
      ! Define system output initializations (set up mesh) here:
      !............................................................................................   
   if (NumSC2CtrlGlob > 0) then
      CALL AllocPAry( SC_DX%y%fromSCglob, NumSC2CtrlGlob, 'y%fromSCglob', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   end if
   
   IF (ErrStat >= AbortErrLev) RETURN
                        
      ! make sure the C versions are synced with these arrays
   if (NumSC2CtrlGlob > 0) then
      SC_DX%y%c_obj%fromSCglob_Len = NumSC2CtrlGlob
      SC_DX%y%c_obj%fromSCglob     = C_LOC( SC_DX%y%fromSCglob(1) )
   else
      SC_DX%y%c_obj%fromSCglob_Len = 0
      SC_DX%y%c_obj%fromSCglob     = C_NULL_PTR
   end if
   
   if (NumSC2Ctrl > 0) then
      CALL AllocPAry( SC_DX%y%fromSC, NumSC2Ctrl, 'y%fromSC', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   end if
   
   IF (ErrStat >= AbortErrLev) RETURN
   
   ! make sure the C versions are synced with these arrays
   if (NumSC2Ctrl > 0) then
      SC_DX%y%c_obj%fromSC_Len = NumSC2Ctrl
      SC_DX%y%c_obj%fromSC     = C_LOC( SC_DX%y%fromSC(1) )
   else
      SC_DX%y%c_obj%fromSC_Len = 0
      SC_DX%y%c_obj%fromSC     = C_NULL_PTR
   end if

   if( (NumSC2CtrlGlob > 0) .or. (NumSC2Ctrl > 0) .or. (NumSC2Ctrl > 0)) then
      SC_DX%p%UseSC = .true. 
   else
      SC_DX%p%UseSC = .false.
   end if

   RETURN
   
END SUBROUTINE SC_DX_Init

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SC_DX_SetInputs(p_FAST, y_SrvD, SC_DX, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(FAST_ParameterType),       INTENT(IN    )  :: p_FAST      ! Parameters for the glue code
   TYPE(SrvD_OutputType),          INTENT(IN)      :: y_SrvD      ! The outputs of the ServoDyn module (control)
   TYPE(SCDataEx_Data),            INTENT(INOUT)   :: SC_DX       ! data for the SuperController integration module
   INTEGER(IntKi),                 INTENT(  OUT)   :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)   :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables
!   INTEGER(IntKi)                                  :: ErrStat2    ! temporary Error status of the operation
!   CHARACTER(ErrMsgLen)                            :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None
                                                   
   CHARACTER(*),   PARAMETER                       :: RoutineName = 'SC_DX_SetInputs'
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
      ! set SuperController inputs
   if (SC_DX%p%UseSC) then
      if (allocated(y_SrvD%toSC).and. associated(SC_DX%u%toSC)) SC_DX%u%toSC = y_SrvD%toSC     
   end if
   
      
END SUBROUTINE SC_DX_SetInputs
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SC_DX_SetOutputs(p_FAST, u_SrvD, SC_DX, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(FAST_ParameterType),       INTENT(IN    )  :: p_FAST      ! Parameters for the glue code
   TYPE(SrvD_InputType),           INTENT(INOUT)   :: u_SrvD      ! The inputs of the ServoDyn module (control)
   TYPE(SCDataEx_Data),            INTENT(IN   )   :: SC_DX       ! data for the SuperController integration module
   INTEGER(IntKi),                 INTENT(  OUT)   :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)   :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables
!   INTEGER(IntKi)                                  :: ErrStat2    ! temporary Error status of the operation
!   CHARACTER(ErrMsgLen)                            :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None
                                                   
   CHARACTER(*),   PARAMETER                       :: RoutineName = 'SC_DX_SetOutputs'
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
      ! set SuperController inputs
   if (SC_DX%p%UseSC) then
      if (allocated(u_SrvD%fromSC)    .and. associated(SC_DX%y%fromSC))     u_SrvD%fromSC     = SC_DX%y%fromSC
      if (allocated(u_SrvD%fromSCglob).and. associated(SC_DX%y%fromSCglob)) u_SrvD%fromSCglob = SC_DX%y%fromSCglob      
   end if
   
      
END SUBROUTINE SC_DX_SetOutputs
!----------------------------------------------------------------------------------------------------------------------------------
END MODULE SC_DataEx
!**********************************************************************************************************************************




