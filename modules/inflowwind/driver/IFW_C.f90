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


INTEGER, PARAMETER :: IntfStrLen        = 1025               !< length of strings through the C interface
INTEGER, PARAMETER :: InputStringLength = 179                !< Fixed length for all lines of the string-based input file
INTEGER, PARAMETER :: InputFileLines    = 55                 !< Number of lines expected in the string-based input file array

CONTAINS

!===============================================================================================================
!--------------------------------------------- IFW INIT --------------------------------------------------------
!===============================================================================================================
SUBROUTINE IFW_INIT_C(InputFileStrings_C, InputUniformStrings_C, NumWindPts_C, DT_C, NumChannels_C, OutputChannelNames_C, OutputChannelUnits_C, ErrStat_C, ErrMsg_C) BIND (C, NAME='IFW_INIT_C')

    TYPE(C_PTR)                                    , INTENT(IN   )   :: InputFileStrings_C
    TYPE(C_PTR)                                    , INTENT(IN   )   :: InputUniformStrings_C
    INTEGER(C_INT)                                 , INTENT(IN   )   :: NumWindPts_C
    REAL(C_DOUBLE)                                 , INTENT(IN   )   :: DT_C
    TYPE(C_PTR)                                    , INTENT(  OUT)   :: OutputChannelNames_C
    TYPE(C_PTR)                                    , INTENT(  OUT)   :: OutputChannelUnits_C
    INTEGER(C_INT)                                 , INTENT(  OUT)   :: NumChannels_C
    INTEGER(C_INT)                                 , INTENT(  OUT)   :: ErrStat_C
    CHARACTER(KIND=C_CHAR)                         , INTENT(  OUT)   :: ErrMsg_C(IntfStrLen) 

    ! Local Variables
    CHARACTER(InputStringLength), DIMENSION(InputFileLines)          :: InputFileStrings
    CHARACTER(kind=C_char, len=1), DIMENSION(:), POINTER             :: character_pointer
    CHARACTER, DIMENSION(InputStringLength)                          :: single_line_character_array
    CHARACTER(InputStringLength)                                     :: single_line_chars
    CHARACTER(CHANLEN+1), ALLOCATABLE, TARGET                        :: tmp_OutputChannelNames_C(:)
    CHARACTER(CHANLEN+1), ALLOCATABLE, TARGET                        :: tmp_OutputChannelUnits_C(:)
    REAL(DbKi)                                                       :: TimeInterval
    INTEGER                                                          :: ErrStat
    CHARACTER(ErrMsgLen)                                             :: ErrMsg
    CHARACTER(LEN=CHANLEN+1),ALLOCATABLE                             :: FLAT_STRING
    INTEGER                                                          :: I, J, K

   ! Convert the string-input from C-style character arrays (char **) to a Fortran-style array of characters.
   ! This expects a fixed number of lines with a fixed length per line.
   ! TODO: Add error checking.
   CALL C_F_pointer(InputFileStrings_C, character_pointer, [InputFileLines * InputStringLength])
   DO i = 0, InputFileLines - 1
      single_line_character_array = character_pointer(i * InputStringLength + 1 : i * InputStringLength + InputStringLength)
      DO j = 1, InputStringLength
         single_line_chars(j:j) = single_line_character_array(j)
      END DO
      InputFileStrings(i + 1) = single_line_chars
   END DO

   ! @RAF: do we need to do the same treatment for the uniform wind input file? If so, it will have a variable number of lines (--> ALLOCATE?)

   ! Store string-inputs as type FileInfoType within InflowWind_InitInputType
   CALL InitFileInfo(InputFileStrings, InitInp%PassedFileData, ErrStat, ErrMsg)           ! in, out (FileInfoType), out, out
   ! Debugging
   ! CALL Print_FileInfo_Struct( CU, InitInp%PassedFileData )
   IF (ErrStat .NE. 0) PRINT *, "IFW_INIT_C: InitFileInfo failed"

   ! Set other inputs for calling InflowWind_Init
   InitInp%UseInputFile          = .false.
   InitInp%InputFileName         = "passed_ifw_file"         ! dummy
   InitInp%RootName              = "ifwRoot"                 ! used for making echo files
   InitInp%NumWindPoints         = NumWindPts_C              ! CHECK THIS!
   InitInp%WindType2UseInputFile = .TRUE.                    ! CHANGE TO FALSE ONCE GET UNIFORM INPUT FILE STRING WORKING!
   TimeInterval                  = REAL(DT_C, DbKi)

   ! Pass FileInfoType into InflowWind_Init - only need InitInp and TimeInterval as inputs, the rest are set by InflowWind_Init
   CALL InflowWind_Init( InitInp, InputGuess, p, ContStates, DiscStates, ConstrStateGuess, OtherStates, y, m, TimeInterval, InitOutData, ErrStat, ErrMsg )
   IF (ErrStat .NE. 0) PRINT *, "IFW_INIT_C: InflowWind_Init failed"
   PRINT*, "Done calling InflowWind_Init ....."

   ! Convert the outputs of InflowWind_Init from Fortran to C

! @ RAF: need some help here please

!DO i = 0, NumChannels_C - 1
!   single_line_character_array = character_pointer(i * InputStringLength + 1 : i * InputStringLength + InputStringLength)
!   DO j = 1, InputStringLength
!      single_line_chars(j:j) = single_line_character_array(j)
!   END DO
!   InputFileStrings(i + 1) = single_line_chars
!END DO
!CALL F_C_pointer(InputFileStrings_C, character_pointer, [InputFileLines * InputStringLength])

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

   PRINT*, "DONE WITH IFW_INIT_C! RETURNING"

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
PRINT *, "Time = ", Time
PRINT *, Positions_C
InputData%PositionXYZ = reshape( real(Positions_C,ReKi), (/3, InitInp%NumWindPoints/) )

! Call InflowWind_CalcOutput to get the velocities
CALL InflowWind_CalcOutput( Time, InputData, p, ContStates, DiscStates, ConstrStates, OtherStates, y, m, ErrStat, ErrMsg )
IF (ErrStat .NE. 0) PRINT *, "IFW_CALCOUTPUT_C: InflowWind_CalcOutput failed"
PRINT*, "Done calling InflowWind_CalcOutput ....."

! Get velocities out of y and flattens it (still in same spot in memory)
Velocities_C = reshape( REAL(y%VelocityUVW, C_FLOAT), (/3*InitInp%NumWindPoints/) ) ! VelocityUVW is 2D array of ReKi (might need reshape or make into pointer); size [3,N]

! NEED TO COPY WriteOutput CHANNEL INFO
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
INTEGER                          :: ErrStat
CHARACTER(ErrMsgLen)             :: ErrMsg

! Need InputData from ?
CALL InflowWind_End( InputData, p, ContStates, DiscStates, ConstrStates, OtherStates, y, m, ErrStat, ErrMsg )

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