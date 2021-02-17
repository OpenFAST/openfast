!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2021 Nicole Mendoza
!
! This file is part of InflowWind.
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
MODULE InflowWindAPI

    USE ISO_C_BINDING
    USE InflowWind
    USE InflowWind_Subs
    USE InflowWind_Types
    USE NWTC_Library

IMPLICIT NONE

PUBLIC :: IFW_INIT_C
PUBLIC :: IFW_CALCOUTPUT_C
PUBLIC :: IFW_END_C

! Accessible to all routines inside module
TYPE(InflowWind_InputType)              :: InputGuess        !< An initial guess for the input; the input mesh must be defined, returned by Init
TYPE(InflowWind_InputType)              :: InputData         !< Created by IFW_CALCOUTPUT_C and used by IFW_END_C
TYPE(InflowWind_ParameterType)          :: p                 !< Parameters
TYPE(InflowWind_ContinuousStateType)    :: ContStates        !< Initial continuous states
TYPE(InflowWind_DiscreteStateType)      :: DiscStates        !< Initial discrete states
TYPE(InflowWind_ConstraintStateType)    :: ConstrStateGuess  !< Initial guess of the constraint states
TYPE(InflowWind_ConstraintStateType)    :: ConstrStates      !< Constraint states at Time
TYPE(InflowWind_OtherStateType)         :: OtherStates       !< Initial other/optimization states
TYPE(InflowWind_OutputType)             :: y                 !< Initial output (outputs are not calculated; only the output mesh is initialized)
TYPE(InflowWind_MiscVarType)            :: m                 !< Misc variables for optimization (not copied in glue code)
TYPE(InflowWind_InitOutputType)         :: InitOutData       !< Initial output data -- Names, units, and version info.

CONTAINS

!===============================================================================================================
!--------------------------------------------- IFW INIT --------------------------------------------------------
!===============================================================================================================
SUBROUTINE IFW_INIT_C(input_file_string_c_ptr,ErrStat_C,ErrMsg_C) BIND (C, NAME='IFW_INIT_C')

    TYPE(C_PTR)                   , INTENT(IN   )      :: input_file_string_c_ptr(*)
!    CHARACTER(KIND=C_CHAR),POINTER, INTENT(INOUT)      :: input_file_string_c(:)
    INTEGER(C_INT)                , INTENT(  OUT)      :: ErrStat_C
    CHARACTER(KIND=C_CHAR)        , INTENT(  OUT)      :: ErrMsg_C

    ! Local Variables    
    CHARACTER(1024), DIMENSION(55)   :: data
    REAL(DbKi)                       :: TimeInterval
    INTEGER                          :: ErrStat
    CHARACTER(ErrMsgLen)             :: ErrMsg
    TYPE(InflowWind_InitInputType)   :: InitInp           

! Convert Fortran character array to FileInfoType within InflowWind_InitInputType
CALL InitFileInfo(data, InitInp%PassedFileData, ErrStat, ErrMsg)           ! in, out (FileInfoType), out, out
PRINT*, ErrStat
PRINT*, "Done calling InitFileInfo ....."

InitInp%UseInputFile = .false.
InitInp%InputFileName = "passed_ifw_file" ! dummy
InitInp%RootName = "ifwRoot" ! used for making echo files
InitInp%NumWindPoints = 100 ! arbirtrary number - how many data points to interrogate for each time step, sets array sizes
TimeInterval = 0.01

! Pass FileInfoType into InflowWind_Init - only need InitInp and TimeInterval as inputs, the rest are set by InflowWind_Init
CALL InflowWind_Init( InitInp, InputGuess, p, ContStates, DiscStates, ConstrStateGuess, OtherStates, y, m, TimeInterval, InitOutData, ErrStat, ErrMsg )
PRINT*, ErrStat
PRINT*, "Done calling InflowWind_Init ....."

! Convert the outputs of InflowWind_Init from Fortran to C
PRINT*, ErrMsg
   if (ErrStat /= 0) then
      ErrStat_c = ErrID_Fatal
   else
      ErrStat_c = ErrID_None
   end if
ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )

CALL InflowWind_CopyInput(InputGuess, InputData, MESH_UPDATECOPY, ErrStat, ErrMsg )
CALL InflowWind_DestroyInput(InputGuess, ErrStat, ErrMsg ) 

CALL InflowWind_CopyConstrState(ConstrStateGuess, ConstrStates, MESH_UPDATECOPY, ErrStat, ErrMsg )
CALL InflowWind_DestroyConstrState(ConstrStateGuess, ErrStat, ErrMsg ) 

! NEED TO COPY WriteOutput CHANNEL INFO

PRINT*, "DONE WITH IFW_INIT_C! RETURNING"

END SUBROUTINE IFW_INIT_C

!===============================================================================================================
!--------------------------------------------- IFW CALCOUTPUT --------------------------------------------------
!===============================================================================================================

SUBROUTINE IFW_CALCOUTPUT_C(Time_C,ErrStat_C,ErrMsg_C) BIND (C, NAME='IFW_CALCOUTPUT_C')

REAL(C_DOUBLE)                , INTENT(IN   )      :: Time_C
INTEGER(C_INT)                , INTENT(  OUT)      :: ErrStat_C
CHARACTER(KIND=C_CHAR)        , INTENT(  OUT)      :: ErrMsg_C

! Local variables
REAL(DbKi)                :: Time
INTEGER                   :: ErrStat
CHARACTER(ErrMsgLen)      :: ErrMsg

! Need to convert C double to DbKi

! Need to add some code to prepare and process inputs

CALL InflowWind_CalcOutput( Time, InputData, p, ContStates, DiscStates, ConstrStates, OtherStates, y, m, ErrStat, ErrMsg )

! NEED TO COPY WriteOutput CHANNEL INFO
! NEED TO COPY SOME INFO OUT OF Y

! Convert the outputs of InflowWind_CalcOutput from Fortran to C
PRINT*, ErrMsg
   if (ErrStat /= 0) then
      ErrStat_C = ErrID_Fatal
   else
      ErrStat_C = ErrID_None
   end if
ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )

PRINT*, "DONE WITH IFW_CALCOUPUT_C!"

END SUBROUTINE IFW_CALCOUTPUT_C

!===============================================================================================================
!--------------------------------------------------- IFW END ---------------------------------------------------
!===============================================================================================================

SUBROUTINE IFW_END_C(ErrStat_C,ErrMsg_C) BIND (C, NAME='IFW_END_C')

INTEGER(C_INT)                , INTENT(  OUT)      :: ErrStat_C
CHARACTER(KIND=C_CHAR)        , INTENT(  OUT)      :: ErrMsg_C

! Local variables
INTEGER                          :: ErrStat
CHARACTER(ErrMsgLen)             :: ErrMsg

! Need InputData from ?
CALL InflowWind_End( InputData, p, ContStates, DiscStates, ConstrStateGuess, OtherStates, y, m, ErrStat, ErrMsg )

! Convert the outputs of InflowWind_End from Fortran to C
PRINT*, ErrMsg
   if (ErrStat /= 0) then
      ErrStat_C = ErrID_Fatal
   else
      ErrStat_C = ErrID_None
   end if
ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )

PRINT*, "DONE WITH IFW_END_C!"

END SUBROUTINE IFW_END_C

END MODULE