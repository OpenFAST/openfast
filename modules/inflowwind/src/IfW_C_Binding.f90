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
   USE VersionInfo

   IMPLICIT NONE

   PUBLIC :: IfW_C_Init
   PUBLIC :: IfW_C_CalcOutput
   PUBLIC :: IfW_C_End

   !------------------------------------------------------------------------------------
   !  Version info for display
   type(ProgDesc), parameter              :: version   = ProgDesc( 'InflowWind library', '', '' )

   !------------------------------------------------------------------------------------
   !  Debugging: DebugLevel -- passed at PreInit
   !     0  - none
   !     1  - some summary info
   !     2  - above + all position/orientation info
   !     3  - above + input files (if direct passed)
   !     4  - above + meshes
   integer(IntKi)                         :: DebugLevel = 0

   !------------------------------------------------------------------------------------
   ! Primary IfW data derived types
   type(InflowWind_InputType)                    :: InputData         !< Inputs to InflowWind
   type(InflowWind_InitInputType)                :: InitInp
   type(InflowWind_InitOutputType)               :: InitOutData       !< Initial output data -- Names, units, and version info.
   type(InflowWind_ParameterType)                :: p                 !< Parameters
   type(InflowWind_ContinuousStateType)          :: ContStates        !< Initial continuous states
   type(InflowWind_DiscreteStateType)            :: DiscStates        !< Initial discrete states
   type(InflowWind_ConstraintStateType)          :: ConstrStates      !< Constraint states at Time
   type(InflowWind_OtherStateType)               :: OtherStates       !< Initial other/optimization states
   type(InflowWind_OutputType)                   :: y                 !< Initial output (outputs are not calculated; only the output mesh is initialized)
   type(InflowWind_MiscVarType)                  :: m                 !< Misc variables for optimization (not copied in glue code)

   !------------------------------------------------------------------------------------
   !  Error handling
   !     This must exactly match the value in the python-lib. If ErrMsgLen changes at
   !     some point in the nwtc-library, this should be updated, but the logic exists
   !     to correctly handle different lengths of the strings
   integer(IntKi),   parameter            :: ErrMsgLen_C = 1025
   integer(IntKi),   parameter            :: IntfStrLen  = 1025       ! length of other strings through the C interface



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
SUBROUTINE IfW_C_Init(IfWinputFilePassed, IfWinputFileString_C, IfWinputFileStringLength_C, OutRootName_C,  &
                     NumWindPts_C, DT_C, DebugLevel_in, NumChannels_C, OutputChannelNames_C, OutputChannelUnits_C,         &
                     ErrStat_C, ErrMsg_C) BIND (C, NAME='IfW_C_Init')
   IMPLICIT NONE
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: IfW_C_Init
!GCC$ ATTRIBUTES DLLEXPORT :: IfW_C_Init
#endif
   integer(c_int),            intent(in   )  :: IfWinputFilePassed                     !< Write VTK outputs [0: none, 1: init only, 2: animation]
   type(c_ptr),               intent(in   )  :: IfWinputFileString_C                   !< Input file as a single string with lines deliniated by C_NULL_CHAR
   integer(c_int),            intent(in   )  :: IfWinputFileStringLength_C             !< lenght of the input file string
   character(kind=c_char),    intent(in   )  :: OutRootName_C(IntfStrLen)              !< Root name to use for echo files and other
   integer(c_int),            intent(in   )  :: NumWindPts_C
   real(c_double),            intent(in   )  :: DT_C
   integer(c_int),            intent(in   )  :: DebugLevel_in
   integer(c_int),            intent(  out)  :: NumChannels_C
   character(kind=c_char),    intent(  out)  :: OutputChannelNames_C(ChanLen*MaxOutPts+1)
   character(kind=c_char),    intent(  out)  :: OutputChannelUnits_C(ChanLen*MaxOutPts+1)
   integer(c_int),            intent(  out)  :: ErrStat_C
   character(kind=c_char),    intent(  out)  :: ErrMsg_C(ErrMsgLen_C)

   ! local variables
   character(IntfStrLen)                                            :: OutRootName       !< Root name to use for echo files and other
   character(IntfStrLen)                                            :: TmpFileName                !< Temporary file name if not passing AD or IfW input file contents directly
   character(kind=c_char, len=IfWinputFileStringLength_C), pointer  :: IfWinputFileString         !< Input file as a single string with NULL chracter separating lines

   real(DbKi)                                                       :: TimeInterval
   integer                                                          :: ErrStat                    !< aggregated error message
   character(ErrMsgLen)                                             :: ErrMsg                     !< aggregated error message
   integer                                                          :: ErrStat2                   !< temporary error status  from a call
   character(ErrMsgLen)                                             :: ErrMsg2                    !< temporary error message from a call
   integer                                                          :: i,j,k
   character(*), parameter                                          :: RoutineName = 'IfW_C_Init' !< for error handling

   ! Initialize error handling
   ErrStat  =  ErrID_None
   ErrMsg   =  ""

   CALL NWTC_Init( ProgNameIn=version%Name )
   CALL DispCopyrightLicense( version%Name )
   CALL DispCompileRuntimeInfo( version%Name )


   ! interface debugging
   DebugLevel = int(DebugLevel_in,IntKi)

   ! Input files
   OutRootName = TRANSFER( OutRootName_C, OutRootName )
   i = INDEX(OutRootName,C_NULL_CHAR) - 1             ! if this has a c null character at the end...
   if ( i > 0 ) OutRootName = OutRootName(1:I)        ! remove it

   ! if non-zero, show all passed data here.  Then check valid values
   if (DebugLevel /= 0_IntKi) then
      call WrScr("   Interface debugging level "//trim(Num2Lstr(DebugLevel))//" requested.")
      call ShowPassedData()
   endif
   ! check valid debug level
   if (DebugLevel < 0_IntKi) then
      ErrStat2 = ErrID_Fatal
      ErrMsg2  = "Interface debug level must be 0 or greater"//NewLine// &
         "  0  - none"//NewLine// &
         "  1  - some summary info and variables passed through interface"//NewLine// &
         "  2  - above + all position/orientation info"//NewLine// &
         "  3  - above + input files (if direct passed)"//NewLine// &
         "  4  - above + meshes"
      if (Failed()) return;
   endif

   ! For debugging the interface:
   if (DebugLevel > 0) then
      call ShowPassedData()
   endif

   ! Get fortran pointer to C_NULL_CHAR deliniated input file as a string 
   CALL C_F_pointer(IfWinputFileString_C, IfWinputFileString)

   ! Format IfW input file contents
   if (IfWinputFilePassed==1_c_int) then
      InitInp%FilePassingMethod   = 1_IntKi                 ! Don't try to read an input -- use passed data instead (blades and AF tables not passed) using FileInfoType
      InitInp%InputFileName       = "passed_ifw_file"       ! not actually used
      call InitFileInfo(IfWinputFileString, InitInp%PassedFileInfo, ErrStat2, ErrMsg2); if (Failed())  return
   else
      InitInp%FilePassingMethod   = 0_IntKi                 ! Read input info from a primary input file
      i = min(IntfStrLen,IfWinputFileStringLength_C)
      TmpFileName = ''
      TmpFileName(1:i) = IfWinputFileString(1:i)
      i = INDEX(TmpFileName,C_NULL_CHAR) - 1                ! if this has a c null character at the end...
      if ( i > 0 ) TmpFileName = TmpFileName(1:I)           ! remove it
      InitInp%InputFileName  = TmpFileName
   endif

   ! For diagnostic purposes, the following can be used to display the contents
   ! of the InFileInfo data structure.
   !     CU is the screen -- system dependent.
   if (DebugLevel >= 3) then
      if (IfWinputFilePassed==1_c_int)    call Print_FileInfo_Struct( CU, InitInp%PassedFileInfo )
   endif

   ! Set other inputs for calling InflowWind_Init
   InitInp%NumWindPoints         = int(NumWindPts_C, IntKi)
   InitInp%RootName              = OutRootName        ! used for making echo files
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

   !> This subroutine prints out all the variables that are passed in.  Use this only
   !! for debugging the interface on the Fortran side.
   subroutine ShowPassedData()
      character(1) :: TmpFlag
      integer      :: i,j
      call WrSCr("")
      call WrScr("-----------------------------------------------------------")
      call WrScr("Interface debugging:  Variables passed in through interface")
      call WrScr("   ADI_C_Init")
      call WrScr("   --------------------------------------------------------")
      call WrScr("   FileInfo")
      TmpFlag="F";   if (IfWinputFilePassed==1_c_int) TmpFlag="T"
      call WrScr("       IfWinputFilePassed_C           "//TmpFlag )
      call WrScr("       IfWinputFileString_C (ptr addr)"//trim(Num2LStr(LOC(IfWinputFileString_C))) )
      call WrScr("       IfWinputFileStringLength_C     "//trim(Num2LStr( IfWinputFileStringLength_C )) )
      call WrScr("       OutRootName                    "//trim(OutRootName) )
      call WrScr("   Input variables")
      call WrScr("       NumWindPts_C                   "//trim(Num2LStr( NumWindPts_C)) )
      call WrScr("   Time variables")
      call WrScr("       DT_C                           "//trim(Num2LStr( DT_C          )) )
      call WrScr("-----------------------------------------------------------")
   end subroutine ShowPassedData

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
