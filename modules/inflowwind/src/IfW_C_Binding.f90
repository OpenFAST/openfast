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
MODULE InflowWind_C_BINDING

    USE ISO_C_BINDING
    USE InflowWind
    USE InflowWind_Subs, only: MaxOutPts
    USE InflowWind_Types
    USE NWTC_Library

   IMPLICIT NONE

   PUBLIC :: IfW_C_Init
   PUBLIC :: IfW_C_CalcOutput
   PUBLIC :: IfW_C_End

   ! Accessible to all routines inside module
   TYPE(InflowWind_InputType)              :: InputData         !< Inputs to InflowWind
   TYPE(InflowWind_InitInputType)          :: InitInp
   TYPE(InflowWind_InitOutputType)         :: InitOutData       !< Initial output data -- Names, units, and version info.
   TYPE(InflowWind_ParameterType)          :: p                 !< Parameters
   TYPE(InflowWind_ContinuousStateType)    :: ContStates        !< Initial continuous states
   TYPE(InflowWind_DiscreteStateType)      :: DiscStates        !< Initial discrete states
   TYPE(InflowWind_ConstraintStateType)    :: ConstrStates      !< Constraint states at Time
   TYPE(InflowWind_OtherStateType)         :: OtherStates       !< Initial other/optimization states
   TYPE(InflowWind_OutputType)             :: y                 !< Initial output (outputs are not calculated; only the output mesh is initialized)
   TYPE(InflowWind_MiscVarType)            :: m                 !< Misc variables for optimization (not copied in glue code)

   !  This must exactly match the value in the Python interface. We are not using the variable 'ErrMsgLen'
   !  so that we avoid issues if ErrMsgLen changes in the NWTC Library. If the value of ErrMsgLen does change
   !  in the NWTC Library, ErrMsgLen_C (and the equivalent value in the Python interface) can be updated 
   !  to be equivalent to ErrMsgLen + 1, but the logic exists to correctly handle different lengths of the strings
   integer(IntKi),   parameter            :: ErrMsgLen_C=1025  ! Numerical equivalent of ErrMsgLen + 1

CONTAINS

!> This routine sets the error status in C_CHAR for export to calling code.
!! Make absolutely certain that we do not overrun the end of ErrMsg_C.  That is hard coded to 1025,
!! but ErrMsgLen is set in the nwtc_library, and could change without updates here.  We don't want an
!! inadvertant buffer overrun -- that can lead to bad things.
subroutine SetErr(ErrStat, ErrMsg, ErrStat_C, ErrMsg_C)
   integer,                intent(in   )  :: ErrStat                 !< aggregated error message (fortran type)
   character(ErrMsgLen),   intent(in   )  :: ErrMsg                  !< aggregated error message (fortran type)
   integer(c_int),         intent(  out)  :: ErrStat_C
   character(kind=c_char), intent(  out)  :: ErrMsg_C(ErrMsgLen_C)
   ErrStat_C = ErrStat     ! We will send back the same error status that is used in OpenFAST
   if (ErrMsgLen > ErrMsgLen_C-1) then   ! If ErrMsgLen is > the space in ErrMsg_C, do not copy everything over
      ErrMsg_C = TRANSFER( trim(ErrMsg(1:ErrMsgLen_C-1))//C_NULL_CHAR, ErrMsg_C )
   else
      ErrMsg_C = TRANSFER( trim(ErrMsg)//C_NULL_CHAR, ErrMsg_C )
   endif
end subroutine SetErr


!===============================================================================================================
!--------------------------------------------- IFW INIT --------------------------------------------------------
!===============================================================================================================
SUBROUTINE IfW_C_Init(InputFileString_C, InputFileStringLength_C, InputUniformString_C, InputUniformStringLength_C, NumWindPts_C, DT_C, NumChannels_C, OutputChannelNames_C, OutputChannelUnits_C, ErrStat_C, ErrMsg_C) BIND (C, NAME='IfW_C_Init')
   IMPLICIT NONE
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: IfW_C_Init
!GCC$ ATTRIBUTES DLLEXPORT :: IfW_C_Init
#endif
    TYPE(C_PTR)                                    , INTENT(IN   )   :: InputFileString_C
    INTEGER(C_INT)                                 , INTENT(IN   )   :: InputFileStringLength_C
    TYPE(C_PTR)                                    , INTENT(IN   )   :: InputUniformString_C
    INTEGER(C_INT)                                 , INTENT(IN   )   :: InputUniformStringLength_C
    INTEGER(C_INT)                                 , INTENT(IN   )   :: NumWindPts_C
    REAL(C_DOUBLE)                                 , INTENT(IN   )   :: DT_C
    INTEGER(C_INT)                                 , INTENT(  OUT)   :: NumChannels_C
    CHARACTER(KIND=C_CHAR)                         , INTENT(  OUT)   :: OutputChannelNames_C(ChanLen*MaxOutPts+1)
    CHARACTER(KIND=C_CHAR)                         , INTENT(  OUT)   :: OutputChannelUnits_C(ChanLen*MaxOutPts+1)
    INTEGER(C_INT)                                 , INTENT(  OUT)   :: ErrStat_C
    CHARACTER(KIND=C_CHAR)                         , INTENT(  OUT)   :: ErrMsg_C(ErrMsgLen_C)

    ! Local Variables
    CHARACTER(kind=C_char, len=InputFileStringLength_C),    POINTER  :: InputFileString            !< Input file as a single string with NULL chracter separating lines
    CHARACTER(kind=C_char, len=InputUniformStringLength_C), POINTER  :: UniformFileString          !< Input file as a single string with NULL chracter separating lines -- Uniform wind file

    REAL(DbKi)                                                       :: TimeInterval
    INTEGER                                                          :: ErrStat                    !< aggregated error message
    CHARACTER(ErrMsgLen)                                             :: ErrMsg                     !< aggregated error message
    INTEGER                                                          :: ErrStat2                   !< temporary error status  from a call
    CHARACTER(ErrMsgLen)                                             :: ErrMsg2                    !< temporary error message from a call
    INTEGER                                                          :: i,j,k
    character(*), parameter                                          :: RoutineName = 'IfW_C_Init' !< for error handling

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
   CALL InflowWind_Init( InitInp, InputData, p, ContStates, DiscStates, ConstrStates, OtherStates, y, m, TimeInterval, InitOutData, ErrStat2, ErrMsg2 )
      if (Failed()) return

   ! Number of channels
   NumChannels_C = size(InitOutData%WriteOutputHdr)

   ! transfer the output channel names and units to c_char arrays for returning
   k=1
   do i=1,NumChannels_C
      do j=1,ChanLen    ! max length of channel name.  Same for units
         OutputChannelNames_C(k)=InitOutData%WriteOutputHdr(i)(j:j)
         OutputChannelUnits_C(k)=InitOutData%WriteOutputUnt(i)(j:j)
         k=k+1
      enddo
   enddo

   ! null terminate the string
   OutputChannelNames_C(k) = C_NULL_CHAR
   OutputChannelUnits_C(k) = C_NULL_CHAR


   call Cleanup()
   call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)

CONTAINS
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed) then
         call Cleanup()
         call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
      endif
   end function Failed
   subroutine Cleanup()    ! NOTE: we are ignoring any error reporting from here
   end subroutine Cleanup 
END SUBROUTINE IfW_C_Init

!===============================================================================================================
!--------------------------------------------- IFW CALCOUTPUT --------------------------------------------------
!===============================================================================================================

SUBROUTINE IfW_C_CalcOutput(Time_C,Positions_C,Velocities_C,OutputChannelValues_C,ErrStat_C,ErrMsg_C) BIND (C, NAME='IfW_C_CalcOutput')
   IMPLICIT NONE
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: IfW_C_CalcOutput
!GCC$ ATTRIBUTES DLLEXPORT :: IfW_C_CalcOutput
#endif
   REAL(C_DOUBLE)                , INTENT(IN   )      :: Time_C
   REAL(C_FLOAT)                 , INTENT(IN   )      :: Positions_C(3*InitInp%NumWindPoints)
   REAL(C_FLOAT)                 , INTENT(  OUT)      :: Velocities_C(3*InitInp%NumWindPoints)
   REAL(C_FLOAT)                 , INTENT(  OUT)      :: OutputChannelValues_C(p%NumOuts)
   INTEGER(C_INT)                , INTENT(  OUT)      :: ErrStat_C
   CHARACTER(KIND=C_CHAR)        , INTENT(  OUT)      :: ErrMsg_C(ErrMsgLen_C)

   ! Local variables
   REAL(DbKi)                                         :: Time
   INTEGER                                            :: ErrStat                          !< aggregated error message
   CHARACTER(ErrMsgLen)                               :: ErrMsg                           !< aggregated error message
   INTEGER                                            :: ErrStat2                         !< temporary error status  from a call
   CHARACTER(ErrMsgLen)                               :: ErrMsg2                          !< temporary error message from a call
   character(*), parameter                            :: RoutineName = 'IfW_C_CalcOutput' !< for error handling

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

   call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)

CONTAINS
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
   end function Failed
END SUBROUTINE IfW_C_CalcOutput

!===============================================================================================================
!--------------------------------------------------- IFW END ---------------------------------------------------
!===============================================================================================================

SUBROUTINE IfW_C_End(ErrStat_C,ErrMsg_C) BIND (C, NAME='IfW_C_End')
   IMPLICIT NONE
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: IfW_C_End
!GCC$ ATTRIBUTES DLLEXPORT :: IfW_C_End
#endif
   INTEGER(C_INT)                , INTENT(  OUT)      :: ErrStat_C
   CHARACTER(KIND=C_CHAR)        , INTENT(  OUT)      :: ErrMsg_C(ErrMsgLen_C)

   ! Local variables
   INTEGER                                            :: ErrStat
   CHARACTER(ErrMsgLen)                               :: ErrMsg

   ! Call the main subroutine InflowWind_End
   CALL InflowWind_End( InputData, p, ContStates, DiscStates, ConstrStates, OtherStates, y, m, ErrStat, ErrMsg )

   call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)

END SUBROUTINE IfW_C_End

END MODULE
