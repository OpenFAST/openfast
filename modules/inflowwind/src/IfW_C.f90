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
TYPE(InflowWind_InitInputType)          :: InitInp
TYPE(InflowWind_InitOutputType)         :: InitOutData       !< Initial output data -- Names, units, and version info.
TYPE(InflowWind_ParameterType)          :: p                 !< Parameters
TYPE(InflowWind_ContinuousStateType)    :: ContStates        !< Initial continuous states
TYPE(InflowWind_DiscreteStateType)      :: DiscStates        !< Initial discrete states
TYPE(InflowWind_ConstraintStateType)    :: ConstrStateGuess  !< Initial guess of the constraint states
TYPE(InflowWind_ConstraintStateType)    :: ConstrStates      !< Constraint states at Time
TYPE(InflowWind_OtherStateType)         :: OtherStates       !< Initial other/optimization states
TYPE(InflowWind_OutputType)             :: y                 !< Initial output (outputs are not calculated; only the output mesh is initialized)
TYPE(InflowWind_MiscVarType)            :: m                 !< Misc variables for optimization (not copied in glue code)

INTEGER, PARAMETER :: InputStringLength = 179                !< Fixed length for all lines of the string-based input file

CONTAINS

!===============================================================================================================
!--------------------------------------------- IFW INIT --------------------------------------------------------
!===============================================================================================================
SUBROUTINE IFW_INIT_C(InputFileStrings_C, InputFileStringLength_C, InputUniformStrings_C, InputUniformStringLength_C, NumWindPts_C, DT_C, NumChannels_C, OutputChannelNames_C, OutputChannelUnits_C, ErrStat_C, ErrMsg_C) BIND (C, NAME='IFW_INIT_C')

    TYPE(C_PTR)                                    , INTENT(IN   )   :: InputFileStrings_C
    INTEGER(C_INT)                                 , INTENT(IN   )   :: InputFileStringLength_C
    TYPE(C_PTR)                                    , INTENT(IN   )   :: InputUniformStrings_C
    INTEGER(C_INT)                                 , INTENT(IN   )   :: InputUniformStringLength_C
    INTEGER(C_INT)                                 , INTENT(IN   )   :: NumWindPts_C
    REAL(C_DOUBLE)                                 , INTENT(IN   )   :: DT_C
    INTEGER(C_INT)                                 , INTENT(  OUT)   :: NumChannels_C
    TYPE(C_PTR)                                    , INTENT(  OUT)   :: OutputChannelNames_C
    TYPE(C_PTR)                                    , INTENT(  OUT)   :: OutputChannelUnits_C
    INTEGER(C_INT)                                 , INTENT(  OUT)   :: ErrStat_C
    CHARACTER(KIND=C_CHAR)                         , INTENT(  OUT)   :: ErrMsg_C(1025) 

    ! Local Variables
    CHARACTER(InputStringLength), DIMENSION(InputFileStringLength_C) :: InputFileStrings
    CHARACTER(InputStringLength), DIMENSION(InputUniformStringLength_C) :: InputUniformStrings
    CHARACTER(kind=C_char, len=1), DIMENSION(:), POINTER             :: character_pointer
    CHARACTER(kind=C_char, len=1), DIMENSION(:), POINTER             :: character_pointer2
    CHARACTER, DIMENSION(InputStringLength)                          :: single_line_character_array
    CHARACTER, DIMENSION(InputStringLength)                          :: single_line_character_array2
    CHARACTER(InputStringLength)                                     :: single_line_chars
    CHARACTER(InputStringLength)                                     :: single_line_chars2

    CHARACTER(CHANLEN+1), ALLOCATABLE, TARGET                        :: tmp_OutputChannelNames_C(:)
    CHARACTER(CHANLEN+1), ALLOCATABLE, TARGET                        :: tmp_OutputChannelUnits_C(:)
    REAL(DbKi)                                                       :: TimeInterval
    INTEGER                                                          :: ErrStat
    CHARACTER(ErrMsgLen)                                             :: ErrMsg
    INTEGER                                                          :: I, J, K

   ! Convert the string-input from C-style character arrays (char**) to a Fortran-style array of characters.
   ! TODO: Add error checking
   CALL C_F_pointer(InputFileStrings_C, character_pointer, [InputFileStringLength_C * InputStringLength])
   DO i = 0, InputFileStringLength_C - 1
      single_line_character_array = character_pointer(i * InputStringLength + 1 : i * InputStringLength + InputStringLength)
      DO j = 1, InputStringLength
         single_line_chars(j:j) = single_line_character_array(j)
      END DO
      InputFileStrings(i + 1) = single_line_chars
   END DO

   CALL C_F_pointer(InputUniformStrings_C, character_pointer2, [InputUniformStringLength_C * InputStringLength])
   DO i = 0, InputUniformStringLength_C - 1
      single_line_character_array2 = character_pointer2(i * InputStringLength + 1 : i * InputStringLength + InputStringLength)
      DO j = 1, InputStringLength
         single_line_chars2(j:j) = single_line_character_array2(j)
      END DO
      InputUniformStrings(i + 1) = single_line_chars2
   END DO

   ! Store string-inputs as type FileInfoType within InflowWind_InitInputType
   CALL InitFileInfo(InputFileStrings, InitInp%PassedFileData, ErrStat, ErrMsg)           
   IF (ErrStat .NE. 0) THEN 
      PRINT *, "IFW_INIT_C: Failed to convert main input file to FileInfoType"
      PRINT *, ErrMsg
   END IF

   CALL InitFileInfo(InputUniformStrings, InitInp%WindType2Data, ErrStat, ErrMsg)        
   IF (ErrStat .NE. 0) THEN 
      PRINT *, "IFW_INIT_C: Failed to convert uniform input file to FileInfoType"
      PRINT *, ErrMsg
   END IF

   ! Set other inputs for calling InflowWind_Init
   InitInp%NumWindPoints         = NumWindPts_C              
   InitInp%InputFileName         = "passed_ifw_file"         ! dummy
   InitInp%RootName              = "ifwRoot"                 ! used for making echo files
   InitInp%UseInputFile          = .FALSE.
   InitInp%WindType2UseInputFile = .FALSE.                   
   TimeInterval                  = REAL(DT_C, DbKi)

   ! Call the main subroutine InflowWind_Init - only need InitInp and TimeInterval as inputs, the rest are set by InflowWind_Init
   CALL InflowWind_Init( InitInp, InputGuess, p, ContStates, DiscStates, ConstrStateGuess, OtherStates, y, m, TimeInterval, InitOutData, ErrStat, ErrMsg )
   IF (ErrStat .NE. 0) THEN 
      PRINT *, "IFW_INIT_C: Main InflowWind_Init subroutine failed!"
      PRINT *, ErrMsg
   ELSE
      PRINT*, "IFW_INIT_C: Successfully called InflowWind_Init ....."
   END IF

   ! Convert the outputs of InflowWind_Init from Fortran to C
   ALLOCATE(tmp_OutputChannelNames_C(size(InitOutData%WriteOutputHdr)))
   ALLOCATE(tmp_OutputChannelUnits_C(size(InitOutData%WriteOutputUnt)))
   NumChannels_C = size(InitOutData%WriteOutputHdr)

   DO I = 1,NumChannels_C
      tmp_OutputChannelNames_C(I) = TRANSFER(InitOutData%WriteOutputHdr(I)//C_NULL_CHAR, tmp_OutputChannelNames_C(I))
      tmp_OutputChannelUnits_C(I) = TRANSFER(InitOutData%WriteOutputUnt(I)//C_NULL_CHAR, tmp_OutputChannelUnits_C(I))
   END DO
   OutputChannelNames_C = C_LOC(tmp_OutputChannelNames_C)
   OutputChannelUnits_C = C_LOC(tmp_OutputChannelUnits_C)

   if (ErrStat /= 0) then
      ErrStat_C = ErrID_Fatal
   else
      ErrStat_C = ErrID_None
   end if
   ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )

   ! Clean up variables and set up for IFW_CALCOUTPUT_C
   CALL InflowWind_CopyInput(InputGuess, InputData, MESH_UPDATECOPY, ErrStat, ErrMsg )
   CALL InflowWind_DestroyInput(InputGuess, ErrStat, ErrMsg ) 

   CALL InflowWind_CopyConstrState(ConstrStateGuess, ConstrStates, MESH_UPDATECOPY, ErrStat, ErrMsg )
   CALL InflowWind_DestroyConstrState(ConstrStateGuess, ErrStat, ErrMsg ) 

   PRINT*, "DONE WITH IFW_INIT_C!"

END SUBROUTINE IFW_INIT_C

!===============================================================================================================
!--------------------------------------------- IFW CALCOUTPUT --------------------------------------------------
!===============================================================================================================

SUBROUTINE IFW_CALCOUTPUT_C(Time_C,Positions_C,Velocities_C,OutputChannelValues_C,ErrStat_C,ErrMsg_C) BIND (C, NAME='IFW_CALCOUTPUT_C')

REAL(C_DOUBLE)                , INTENT(IN   )      :: Time_C
REAL(C_FLOAT)                 , INTENT(IN   )      :: Positions_C(3*InitInp%NumWindPoints)
REAL(C_FLOAT)                 , INTENT(  OUT)      :: Velocities_C(3*InitInp%NumWindPoints)
REAL(C_FLOAT)                 , INTENT(  OUT)      :: OutputChannelValues_C(p%NumOuts)
INTEGER(C_INT)                , INTENT(  OUT)      :: ErrStat_C
CHARACTER(KIND=C_CHAR)        , INTENT(  OUT)      :: ErrMsg_C

! Local variables
REAL(DbKi)                                         :: Time
INTEGER                                            :: ErrStat
CHARACTER(ErrMsgLen)                               :: ErrMsg

! Convert the inputs from C to Fortran
Time = REAL(Time_C,DbKi)
InputData%PositionXYZ = reshape( real(Positions_C,ReKi), (/3, InitInp%NumWindPoints/) )

! Call the main subroutine InflowWind_CalcOutput to get the velocities
CALL InflowWind_CalcOutput( Time, InputData, p, ContStates, DiscStates, ConstrStates, OtherStates, y, m, ErrStat, ErrMsg )
IF (ErrStat .NE. 0) THEN
   PRINT *, "IFW_CALCOUTPUT_C: InflowWind_CalcOutput failed"
   PRINT *, ErrMsg
ELSE
   PRINT*, "IFW_CALCOUTPUT_C: Successfully called InflowWind_CalcOutput ....."
END IF

! Get velocities out of y and flatten them (still in same spot in memory)
Velocities_C = reshape( REAL(y%VelocityUVW, C_FLOAT), (/3*InitInp%NumWindPoints/) ) ! VelocityUVW is 2D array of ReKi (might need reshape or make into pointer); size [3,N]

! Get the output channel info out of y
OutputChannelValues_C = REAL(y%WriteOutput, C_FLOAT)

! Convert the outputs of InflowWind_CalcOutput from Fortran to C
if (ErrStat /= 0) then
   ErrStat_C = ErrID_Fatal
else
   ErrStat_C = ErrID_None
end if
ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )

PRINT*, "DONE WITH IFW_CALCOUTPUT_C!"

END SUBROUTINE IFW_CALCOUTPUT_C

!===============================================================================================================
!--------------------------------------------------- IFW END ---------------------------------------------------
!===============================================================================================================

SUBROUTINE IFW_END_C(ErrStat_C,ErrMsg_C) BIND (C, NAME='IFW_END_C')

INTEGER(C_INT)                , INTENT(  OUT)      :: ErrStat_C
CHARACTER(KIND=C_CHAR)        , INTENT(  OUT)      :: ErrMsg_C

! Local variables
INTEGER                                            :: ErrStat
CHARACTER(ErrMsgLen)                               :: ErrMsg

! Call the main subroutine InflowWind_End
CALL InflowWind_End( InputData, p, ContStates, DiscStates, ConstrStates, OtherStates, y, m, ErrStat, ErrMsg )
IF (ErrStat .NE. 0) THEN
   PRINT *, "IFW_END_C: InflowWind_End failed"
   PRINT *, ErrMsg
ELSE
   PRINT*, "IFW_END_C: Successfully called InflowWind_END ....."
END IF

! Convert the outputs of InflowWind_End from Fortran to C
if (ErrStat /= 0) then
   ErrStat_C = ErrID_Fatal
else
   ErrStat_C = ErrID_None
end if
ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )

PRINT*, "DONE WITH IFW_END_C!"

END SUBROUTINE IFW_END_C

END MODULE