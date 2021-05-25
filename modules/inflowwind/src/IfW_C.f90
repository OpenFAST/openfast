!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2021 National Renewable Energy Laboratory
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

CONTAINS

!===============================================================================================================
!--------------------------------------------- IFW INIT --------------------------------------------------------
!===============================================================================================================
SUBROUTINE IFW_INIT_C(InputFileString_C, InputFileStringLength_C, InputUniformString_C, InputUniformStringLength_C, NumWindPts_C, DT_C, NumChannels_C, OutputChannelNames_C, OutputChannelUnits_C, ErrStat_C, ErrMsg_C) BIND (C, NAME='IFW_INIT_C')

    TYPE(C_PTR)                                    , INTENT(IN   )   :: InputFileString_C
    INTEGER(C_INT)                                 , INTENT(IN   )   :: InputFileStringLength_C
    TYPE(C_PTR)                                    , INTENT(IN   )   :: InputUniformString_C
    INTEGER(C_INT)                                 , INTENT(IN   )   :: InputUniformStringLength_C
    INTEGER(C_INT)                                 , INTENT(IN   )   :: NumWindPts_C
    REAL(C_DOUBLE)                                 , INTENT(IN   )   :: DT_C
    INTEGER(C_INT)                                 , INTENT(  OUT)   :: NumChannels_C
    TYPE(C_PTR)                                    , INTENT(  OUT)   :: OutputChannelNames_C
    TYPE(C_PTR)                                    , INTENT(  OUT)   :: OutputChannelUnits_C
    INTEGER(C_INT)                                 , INTENT(  OUT)   :: ErrStat_C
    CHARACTER(KIND=C_CHAR)                         , INTENT(  OUT)   :: ErrMsg_C(1025)

    ! Local Variables
    CHARACTER(kind=C_char, len=InputFileStringLength_C),    POINTER  :: InputFileString            !< Input file as a single string with NULL chracter separating lines
    CHARACTER(kind=C_char, len=InputUniformStringLength_C), POINTER  :: UniformFileString          !< Input file as a single string with NULL chracter separating lines -- Uniform wind file

    CHARACTER(CHANLEN+1), ALLOCATABLE, TARGET                        :: tmp_OutputChannelNames_C(:)
    CHARACTER(CHANLEN+1), ALLOCATABLE, TARGET                        :: tmp_OutputChannelUnits_C(:)
    REAL(DbKi)                                                       :: TimeInterval
    INTEGER                                                          :: ErrStat                    !< aggregated error message
    CHARACTER(ErrMsgLen)                                             :: ErrMsg                     !< aggregated error message
    INTEGER                                                          :: ErrStat2                   !< temporary error status  from a call
    CHARACTER(ErrMsgLen)                                             :: ErrMsg2                    !< temporary error message from a call
    INTEGER                                                          :: I
    character(*), parameter                                          :: RoutineName = 'IFW_INIT_C' !< for error handling

   ! Initialize error handling
   ErrStat  =  ErrID_None
   ErrMsg   =  ""

   ! Get fortran pointer to C_NULL_CHAR deliniated input file as a string 
   CALL C_F_pointer(InputFileString_C, InputFileString)
   CALL C_F_pointer(InputUniformString_C, UniformFileString)

   ! Store string-inputs as type FileInfoType within InflowWind_InitInputType
   CALL InitFileInfo(InputFileString, InitInp%PassedFileData, ErrStat2, ErrMsg2);   if (Failed())  return
   InitInp%UseInputFile          = .FALSE.

   ! store Uniform File strings if they are non-zero sized
   if (len(UniformFileString) > 1) then 
      CALL InitFileInfo(UniformFileString, InitInp%WindType2Data, ErrStat2, ErrMsg2);  if (Failed())  return
      InitInp%WindType2UseInputFile = .FALSE.
   else  ! Default to reading from disk
      InitInp%WindType2UseInputFile = .TRUE.
   endif

   ! Set other inputs for calling InflowWind_Init
   InitInp%NumWindPoints         = NumWindPts_C              
   InitInp%InputFileName         = "passed_ifw_file"         ! dummy
   InitInp%RootName              = "ifwRoot"                 ! used for making echo files
   TimeInterval                  = REAL(DT_C, DbKi)

   ! Call the main subroutine InflowWind_Init - only need InitInp and TimeInterval as inputs, the rest are set by InflowWind_Init
   CALL InflowWind_Init( InitInp, InputGuess, p, ContStates, DiscStates, ConstrStateGuess, OtherStates, y, m, TimeInterval, InitOutData, ErrStat2, ErrMsg2 )
      if (Failed())  return

   ! Convert the outputs of InflowWind_Init from Fortran to C
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

   ! Clean up variables and set up for IFW_CALCOUTPUT_C
   CALL InflowWind_CopyInput(InputGuess, InputData, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
      if (Failed())  return
   CALL InflowWind_CopyConstrState(ConstrStateGuess, ConstrStates, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
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
      CALL InflowWind_DestroyInput(InputGuess, ErrStat2, ErrMsg2 )
      CALL InflowWind_DestroyConstrState(ConstrStateGuess, ErrStat2, ErrMsg2 )
   end subroutine Cleanup 
   subroutine SetErr()
      ! Make absolutely certain that we do not overrun the end of ErrMsg_C.  That is hard coded to 1025,
      ! but ErrMsgLen is set in the nwtc_library, and could change without updates here.  We don't want an
      ! inadvertant buffer overrun -- that can lead to bad things.
      integer(IntKi) :: CMsgLen
      ErrStat_C = ErrStat     ! We will send back the same error status that is used in OpenFAST
      CMsgLen = size(ErrMsg_C) - 1    ! Max length of ErrMsg_C without C_NULL_CHAR (probably 1024)
      if (ErrMsgLen > CMsgLen) then   ! If ErrMsgLen is > the space in ErrMsg_C, do not copy everything over
         ErrMsg_C = TRANSFER( trim(ErrMsg(1:CMsgLen))//C_NULL_CHAR, ErrMsg_C )
      else
         ErrMsg_C = TRANSFER( trim(ErrMsg)//C_NULL_CHAR, ErrMsg_C )
      endif
   end subroutine SetErr

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
   CHARACTER(KIND=C_CHAR)        , INTENT(  OUT)      :: ErrMsg_C(1025)

   ! Local variables
   REAL(DbKi)                                         :: Time
   INTEGER                                            :: ErrStat                          !< aggregated error message
   CHARACTER(ErrMsgLen)                               :: ErrMsg                           !< aggregated error message
   INTEGER                                            :: ErrStat2                         !< temporary error status  from a call
   CHARACTER(ErrMsgLen)                               :: ErrMsg2                          !< temporary error message from a call
   character(*), parameter                            :: RoutineName = 'IFW_CALCOUTPUT_C' !< for error handling

   ! Initialize error handling
   ErrStat  =  ErrID_None
   ErrMsg   =  ""

   ! Convert the inputs from C to Fortran
   Time = REAL(Time_C,DbKi)
   InputData%PositionXYZ = reshape( real(Positions_C,ReKi), (/3, InitInp%NumWindPoints/) )

   ! Call the main subroutine InflowWind_CalcOutput to get the velocities
   CALL InflowWind_CalcOutput( Time, InputData, p, ContStates, DiscStates, ConstrStates, OtherStates, y, m, ErrStat2, ErrMsg2 )
      if (Failed())  return

   ! Get velocities out of y and flatten them (still in same spot in memory)
   Velocities_C = reshape( REAL(y%VelocityUVW, C_FLOAT), (/3*InitInp%NumWindPoints/) ) ! VelocityUVW is 2D array of ReKi (might need reshape or make into pointer); size [3,N]

   ! Get the output channel info out of y
   OutputChannelValues_C = REAL(y%WriteOutput, C_FLOAT)

   call SetErr()

CONTAINS
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call SetErr()
   end function Failed
   subroutine SetErr()
      ! Make absolutely certain that we do not overrun the end of ErrMsg_C.  That is hard coded to 1025,
      ! but ErrMsgLen is set in the nwtc_library, and could change without updates here.  We don't want an
      ! inadvertant buffer overrun -- that can lead to bad things.
      integer(IntKi) :: CMsgLen
      ErrStat_C = ErrStat     ! We will send back the same error status that is used in OpenFAST
      CMsgLen = size(ErrMsg_C) - 1    ! Max length of ErrMsg_C without C_NULL_CHAR
      if (ErrMsgLen > CMsgLen) then   ! If ErrMsgLen is > the space in ErrMsg_C, do not copy everything over
         ErrMsg_C = TRANSFER( trim(ErrMsg(1:CMsgLen))//C_NULL_CHAR, ErrMsg_C )
      else
         ErrMsg_C = TRANSFER( trim(ErrMsg)//C_NULL_CHAR, ErrMsg_C )
      endif
   end subroutine SetErr
END SUBROUTINE IFW_CALCOUTPUT_C

!===============================================================================================================
!--------------------------------------------------- IFW END ---------------------------------------------------
!===============================================================================================================

SUBROUTINE IFW_END_C(ErrStat_C,ErrMsg_C) BIND (C, NAME='IFW_END_C')

   INTEGER(C_INT)                , INTENT(  OUT)      :: ErrStat_C
   CHARACTER(KIND=C_CHAR)        , INTENT(  OUT)      :: ErrMsg_C(1025)

   ! Local variables
   INTEGER                                            :: ErrStat
   CHARACTER(ErrMsgLen)                               :: ErrMsg

   ! Call the main subroutine InflowWind_End
   CALL InflowWind_End( InputData, p, ContStates, DiscStates, ConstrStates, OtherStates, y, m, ErrStat, ErrMsg )

   call SetErr()

CONTAINS
   subroutine SetErr()
      ! Make absolutely certain that we do not overrun the end of ErrMsg_C.  That is hard coded to 1025,
      ! but ErrMsgLen is set in the nwtc_library, and could change without updates here.  We don't want an
      ! inadvertant buffer overrun -- that can lead to bad things.
      integer(IntKi) :: CMsgLen
      ErrStat_C = ErrStat     ! We will send back the same error status that is used in OpenFAST
      CMsgLen = size(ErrMsg_C) - 1    ! Max length of ErrMsg_C without C_NULL_CHAR
      if (ErrMsgLen > CMsgLen) then   ! If ErrMsgLen is > the space in ErrMsg_C, do not copy everything over
         ErrMsg_C = TRANSFER( trim(ErrMsg(1:CMsgLen))//C_NULL_CHAR, ErrMsg_C )
      else
         ErrMsg_C = TRANSFER( trim(ErrMsg)//C_NULL_CHAR, ErrMsg_C )
      endif
   end subroutine SetErr
END SUBROUTINE IFW_END_C

END MODULE
