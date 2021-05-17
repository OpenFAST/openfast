!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2021 Nicole Mendoza
!
! This file is part of HydroDyn.
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
MODULE HydroDynAPI

    USE ISO_C_BINDING
    USE HydroDyn
    USE HydroDyn_Types
    USE NWTC_Library

IMPLICIT NONE

PUBLIC :: HydroDyn_Init_C
PUBLIC :: HydroDyn_CalcOutput_C
PUBLIC :: HydroDyn_UpdateStates_C
PUBLIC :: HydroDyn_End_C

! Accessible to all routines inside module
TYPE(HydroDyn_InputType)              :: InputGuess        !< An initial guess for the input; the input mesh must be defined, returned by Init
TYPE(HydroDyn_InputType)              :: InputData         !< Created by IFW_CALCOUTPUT_C and used by IFW_END_C
TYPE(HydroDyn_InitInputType)          :: InitInp
TYPE(HydroDyn_InitOutputType)         :: InitOutData       !< Initial output data -- Names, units, and version info.
TYPE(HydroDyn_ParameterType)          :: p                 !< Parameters
TYPE(HydroDyn_ContinuousStateType)    :: ContStates        !< Initial continuous states
TYPE(HydroDyn_DiscreteStateType)      :: DiscStates        !< Initial discrete states
TYPE(HydroDyn_ConstraintStateType)    :: ConstrStateGuess  !< Initial guess of the constraint states
TYPE(HydroDyn_ConstraintStateType)    :: ConstrStates      !< Constraint states at Time
TYPE(HydroDyn_OtherStateType)         :: OtherStates       !< Initial other/optimization states
TYPE(HydroDyn_OutputType)             :: y                 !< Initial output (outputs are not calculated; only the output mesh is initialized)
TYPE(HydroDyn_MiscVarType)            :: m                 !< Misc variables for optimization (not copied in glue code)

!FIXME: we cannot have this fixed!!!!!
INTEGER, PARAMETER :: InputStringLength = 647                !< Fixed length for all lines of the string-based input file

CONTAINS

!===============================================================================================================
!--------------------------------------------- HydroDyn Init----------------------------------------------------
!===============================================================================================================
SUBROUTINE HydroDyn_Init_C(InputFileStrings_C, InputFileStringLength_C, DT_C, NumChannels_C, OutputChannelNames_C, OutputChannelUnits_C, ErrStat_C, ErrMsg_C) BIND (C, NAME='HydroDyn_Init_C')

    TYPE(C_PTR)                                    , INTENT(IN   )   :: InputFileStrings_C
    INTEGER(C_INT)                                 , INTENT(IN   )   :: InputFileStringLength_C
    REAL(C_DOUBLE)                                 , INTENT(IN   )   :: DT_C
    INTEGER(C_INT)                                 , INTENT(  OUT)   :: NumChannels_C
    TYPE(C_PTR)                                    , INTENT(  OUT)   :: OutputChannelNames_C
    TYPE(C_PTR)                                    , INTENT(  OUT)   :: OutputChannelUnits_C
    INTEGER(C_INT)                                 , INTENT(  OUT)   :: ErrStat_C
    CHARACTER(KIND=C_CHAR)                         , INTENT(  OUT)   :: ErrMsg_C(1025)

    ! Local Variables
    CHARACTER(InputStringLength), DIMENSION(InputFileStringLength_C) :: InputFileStrings
    CHARACTER(kind=C_char, len=1), DIMENSION(:), POINTER             :: character_pointer
    CHARACTER(kind=C_char, len=1), DIMENSION(:), POINTER             :: character_pointer2
    CHARACTER, DIMENSION(InputStringLength)                          :: single_line_character_array
    CHARACTER, DIMENSION(InputStringLength)                          :: single_line_character_array2
    CHARACTER(InputStringLength)                                     :: single_line_chars
    CHARACTER(InputStringLength)                                     :: single_line_chars2

    CHARACTER(CHANLEN+1), ALLOCATABLE, TARGET                        :: tmp_OutputChannelNames_C(:)
    CHARACTER(CHANLEN+1), ALLOCATABLE, TARGET                        :: tmp_OutputChannelUnits_C(:)
    REAL(DbKi)                                                       :: TimeInterval
    INTEGER                                                          :: ErrStat                    !< aggregated error message
    CHARACTER(ErrMsgLen)                                             :: ErrMsg                     !< aggregated error message
    INTEGER                                                          :: ErrStat2                   !< temporary error status  from a call
    CHARACTER(ErrMsgLen)                                             :: ErrMsg2                    !< temporary error message from a call
    INTEGER                                                          :: I, J, K
    character(*), parameter                                          :: RoutineName = 'HydroDyn_Init_C' !< for error handling

   ! Initialize error handling
   ErrStat  =  ErrID_None
   ErrMsg   =  ""

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

   ! Store string-inputs as type FileInfoType within HydroDyn_InitInputType
   CALL InitFileInfo(InputFileStrings, InitInp%PassedFileData, ErrStat2, ErrMsg2)
      if (Failed())  return

   ! Set other inputs for calling HydroDyn_Init
   InitInp%InputFile             = "passed_hd_file"         ! dummy
   InitInp%OutRootName           = "HDpassed"               ! used for making echo files
   InitInp%UseInputFile          = .FALSE.
   TimeInterval                  = REAL(DT_C, DbKi)
   InitInp%Linearize             = .FALSE.
   InitInp%Gravity               = 9.80665   ! Gravity (m/s^2)
   InitInp%defWtrDens            = 1025      ! Water density (kg/m^3)
   InitInp%defWtrDpth            = 200       ! Water depth (m)
   InitInp%defMSL2SWL            = 0         ! Offset between still-water level and mean sea level (m) [positive upward]
   InitInp%TMax                  = 600
   InitInp%HasIce                = .FALSE.
!   InitInp%WaveElevXY              ! Don't allocate this
   InitInp%PtfmLocationX         = 0.0_ReKi
   InitInp%PtfmLocationY         = 0.0_ReKi


   ! Call the main subroutine HydroDyn_Init - only need InitInp and TimeInterval as inputs, the rest are set by HydroDyn_Init
   CALL HydroDyn_Init( InitInp, InputGuess, p, ContStates, DiscStates, ConstrStateGuess, OtherStates, y, m, TimeInterval, InitOutData, ErrStat2, ErrMsg2 )
      if (Failed())  return

   ! Convert the outputs of HydroDyn_Init from Fortran to C
   ALLOCATE(tmp_OutputChannelNames_C(size(InitOutData%WriteOutputHdr)),STAT=ErrStat2)
      if (ErrStat2 /= 0) then
         ErrStat2 = ErrID_Fatal
         ErrMsg2  = "Could not allocate WriteOutputHdr array"
      endif
      if (Failed())  return
   ALLOCATE(tmp_OutputChannelUnits_C(size(InitOutData%WriteOutputUnt)),STAT=ErrStat2)
      if (ErrStat2 /= 0) then
         ErrStat2 = ErrID_Fatal
         ErrMsg2  = "Could not allocate WriteOutputUnt array"
      endif
      if (Failed())  return
   NumChannels_C = size(InitOutData%WriteOutputHdr)

   DO I = 1,NumChannels_C
      tmp_OutputChannelNames_C(I) = TRANSFER(InitOutData%WriteOutputHdr(I)//C_NULL_CHAR, tmp_OutputChannelNames_C(I))
      tmp_OutputChannelUnits_C(I) = TRANSFER(InitOutData%WriteOutputUnt(I)//C_NULL_CHAR, tmp_OutputChannelUnits_C(I))
   END DO
   OutputChannelNames_C = C_LOC(tmp_OutputChannelNames_C)
   OutputChannelUnits_C = C_LOC(tmp_OutputChannelUnits_C)

   ! Clean up variables and set up for HydroDyn_CalcOutput_C
   CALL HydroDyn_CopyInput(InputGuess, InputData, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
      if (Failed())  return
   CALL HydroDyn_CopyConstrState(ConstrStateGuess, ConstrStates, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
      if (Failed())  return

   call Cleanup()
   call SetErr()

CONTAINS
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed) then
         call Cleanup()
         call SetErr()
      endif
   end function Failed
   subroutine Cleanup()    ! NOTE: we are ignoring any error reporting from here
      CALL HydroDyn_DestroyInput(InputGuess, ErrStat2, ErrMsg2 )
      CALL HydroDyn_DestroyConstrState(ConstrStateGuess, ErrStat2, ErrMsg2 )
   end subroutine Cleanup
   subroutine SetErr()
      ErrStat_C = ErrStat     ! We will send back the same error status that is used in OpenFAST
      ErrMsg_C = TRANSFER( trim(ErrMsg)//C_NULL_CHAR, ErrMsg_C )
   end subroutine SetErr

END SUBROUTINE HydroDyn_Init_C

!===============================================================================================================
!--------------------------------------------- HydroDyn CalcOutput ---------------------------------------------
!===============================================================================================================

SUBROUTINE HydroDyn_CalcOutput_C(Time_C,OutputChannelValues_C,ErrStat_C,ErrMsg_C) BIND (C, NAME='HydroDyn_CalcOutput_C')
   REAL(C_DOUBLE)                , INTENT(IN   )      :: Time_C
   REAL(C_FLOAT)                 , INTENT(  OUT)      :: OutputChannelValues_C(p%NumOuts)
   INTEGER(C_INT)                , INTENT(  OUT)      :: ErrStat_C
   CHARACTER(KIND=C_CHAR)        , INTENT(  OUT)      :: ErrMsg_C

   ! Local variables
   REAL(DbKi)                                         :: Time
   INTEGER                                            :: ErrStat                          !< aggregated error message
   CHARACTER(ErrMsgLen)                               :: ErrMsg                           !< aggregated error message
   INTEGER                                            :: ErrStat2                         !< temporary error status  from a call
   CHARACTER(ErrMsgLen)                               :: ErrMsg2                          !< temporary error message from a call
   character(*), parameter                            :: RoutineName = 'HydroDyn_CalcOutput_C' !< for error handling

   ! Initialize error handling
   ErrStat  =  ErrID_None
   ErrMsg   =  ""

   ! Convert the inputs from C to Fortrn
   Time = REAL(Time_C,DbKi)
!   InputData%PositionXYZ = reshape( real(Positions_C,ReKi), (/3, InitInp%NumWindPoints/) )

!   ! Call the main subroutine HydroDyn_CalcOutput to get the velocities
!   CALL HydroDyn_CalcOutput( Time, InputData, p, ContStates, DiscStates, ConstrStates, OtherStates, y, m, ErrStat2, ErrMsg2 )
!      if (Failed())  return

!   ! Get velocities out of y and flatten them (still in same spot in memory)
!   Velocities_C = reshape( REAL(y%VelocityUVW, C_FLOAT), (/3*InitInp%NumWindPoints/) ) ! VelocityUVW is 2D array of ReKi (might need reshape or make into pointer); size [3,N]

!   ! Get the output channel info out of y
!   OutputChannelValues_C = REAL(y%WriteOutput, C_FLOAT)

   call SetErr()

CONTAINS
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call SetErr()
   end function Failed
   subroutine SetErr()
      ErrStat_C = ErrStat     ! We will send back the same error status that is used in OpenFAST
      ErrMsg_C = TRANSFER( trim(ErrMsg)//C_NULL_CHAR, ErrMsg_C )
   end subroutine SetErr
END SUBROUTINE HydroDyn_CalcOutput_C

!===============================================================================================================
!--------------------------------------------- HydroDyn UpdateStates ---------------------------------------------
!===============================================================================================================

SUBROUTINE HydroDyn_UpdateStates_C(Time_C,ErrStat_C,ErrMsg_C) BIND (C, NAME='HydroDyn_UpdateStates_C')
   REAL(C_DOUBLE)                , INTENT(IN   )      :: Time_C
   INTEGER(C_INT)                , INTENT(  OUT)      :: ErrStat_C
   CHARACTER(KIND=C_CHAR)        , INTENT(  OUT)      :: ErrMsg_C
   ! Local variables
   REAL(DbKi)                                         :: Time
   INTEGER                                            :: ErrStat                          !< aggregated error message
   CHARACTER(ErrMsgLen)                               :: ErrMsg                           !< aggregated error message
   INTEGER                                            :: ErrStat2                         !< temporary error status  from a call
   CHARACTER(ErrMsgLen)                               :: ErrMsg2                          !< temporary error message from a call
   character(*), parameter                            :: RoutineName = 'HydroDyn_UpdateStates_C' !< for error handling

   ! Initialize error handling
   ErrStat  =  ErrID_None
   ErrMsg   =  ""

   ! Convert the inputs from C to Fortran
   Time = REAL(Time_C,DbKi)

!   ! Call the main subroutine HydroDyn_UpdateStates to get the velocities
!   CALL HydroDyn_UpdateStates( Time, InputData, p, ContStates, DiscStates, ConstrStates, OtherStates, y, m, ErrStat2, ErrMsg2 )
!      if (Failed())  return

   call SetErr()

CONTAINS
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call SetErr()
   end function Failed
   subroutine SetErr()
      ErrStat_C = ErrStat     ! We will send back the same error status that is used in OpenFAST
      ErrMsg_C = TRANSFER( trim(ErrMsg)//C_NULL_CHAR, ErrMsg_C )
   end subroutine SetErr
END SUBROUTINE HydroDyn_UpdateStates_C

!===============================================================================================================
!--------------------------------------------------- HydroDyn End-----------------------------------------------
!===============================================================================================================

SUBROUTINE HydroDyn_End_C(ErrStat_C,ErrMsg_C) BIND (C, NAME='HydroDyn_End_C')

   INTEGER(C_INT)                , INTENT(  OUT)      :: ErrStat_C
   CHARACTER(KIND=C_CHAR)        , INTENT(  OUT)      :: ErrMsg_C

   ! Local variables
   INTEGER                                            :: ErrStat
   CHARACTER(ErrMsgLen)                               :: ErrMsg

   ! Call the main subroutine HydroDyn_End
   CALL HydroDyn_End( InputData, p, ContStates, DiscStates, ConstrStates, OtherStates, y, m, ErrStat, ErrMsg )

   call SetErr()

CONTAINS
   subroutine SetErr()
      ErrStat_C = ErrStat     ! We will send back the same error status that is used in OpenFAST
      ErrMsg_C = TRANSFER( trim(ErrMsg)//C_NULL_CHAR, ErrMsg_C )
   end subroutine SetErr
END SUBROUTINE HydroDyn_End_C

END MODULE HydroDynAPI
