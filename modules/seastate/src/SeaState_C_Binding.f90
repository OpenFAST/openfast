!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2025 National Renewable Energy Lab
!
! This file is part of SeaState.
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
MODULE SeaState_C_Binding

    USE ISO_C_BINDING
    USE SeaState
    USE SeaState_Types
    USE SeaState_Output
    USE NWTC_Library
    USE NWTC_C_Binding, ONLY: ErrMsgLen_C, IntfStrLen, SetErr, FileNameFromCString
    USE VersionInfo

    IMPLICIT NONE
    SAVE

    PUBLIC :: SeaSt_C_Init
    PUBLIC :: SeaSt_C_CalcOutput
    PUBLIC :: SeaSt_C_End

    !------------------------------------------------------------------------------------
    !  Version info for display
    type(ProgDesc), parameter              :: version   = SeaSt_ProgDesc

    !------------------------------------------------------------------------------------
    !  Debugging: DebugLevel -- passed at PreInit
    !     0  - none
    !     1  - some summary info
    !     2  - above + all position/orientation info
    !     3  - above + input files (if direct passed)
    !     4  - above + meshes
    integer(IntKi)                         :: DebugLevel = 0

    !------------------------------
    !  Primary derived types
    type(SeaSt_InputType)                    :: InputData         !< Inputs to SeaState
    type(SeaSt_InitInputType)                :: InitInp
    type(SeaSt_InitOutputType)               :: InitOutData       !< Initial output data -- Names, units, and version info.
    type(SeaSt_ParameterType)                :: p                 !< Parameters
    type(SeaSt_OutputType)                   :: y                 !< Initial output (outputs are not calculated; only the output mesh is initialized)
    type(SeaSt_MiscVarType)                  :: m                 !< Misc variables for optimization (not copied in glue code)

CONTAINS


subroutine SeaSt_C_Init(InputFile_c, OutRootName_c, Gravity_c, WtrDens_c, WtrDpth_c, MSL2SWL_c, NSteps_c, TimeInterval_c, WaveElevSeriesFlag_c, WrWvKinMod_c, NumChannels_C, OutputChannelNames_C, OutputChannelUnits_C, ErrStat_C, ErrMsg_C) BIND (C, NAME='SeaSt_C_Init')
implicit none
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: SeaSt_C_Init
!GCC$ ATTRIBUTES DLLEXPORT :: SeaSt_C_Init
#endif
    type(c_ptr),                intent(in   ) :: InputFile_c
    type(c_ptr),                intent(in   ) :: OutRootName_c
    real(c_float),              intent(in   ) :: Gravity_c
    real(c_float),              intent(in   ) :: WtrDens_c
    real(c_float),              intent(in   ) :: WtrDpth_c
    real(c_float),              intent(in   ) :: MSL2SWL_c
    integer(c_int),             intent(in   ) :: NSteps_c
    real(c_float),              intent(in   ) :: TimeInterval_c
    integer(c_int),             intent(in   ) :: WaveElevSeriesFlag_c
    integer(c_int),             intent(in   ) :: WrWvKinMod_c
    integer(c_int),             intent(  out) :: NumChannels_C
    character(kind=c_char),     intent(  out) :: OutputChannelNames_C(ChanLen*MaxOutPts+1)
    character(kind=c_char),     intent(  out) :: OutputChannelUnits_C(ChanLen*MaxOutPts+1)
    integer(c_int),             intent(  out) :: ErrStat_C
    character(kind=c_char),     intent(  out) :: ErrMsg_C(ErrMsgLen_C)

    ! Local variables
    character(kind=c_char, len=IntfStrLen), pointer :: InputFileString          !< Input file as a single string with NULL chracter separating lines
    character(kind=c_char, len=IntfStrLen), pointer :: OutputFileString          !< Input file as a single string with NULL chracter separating lines
    character(IntfStrLen)           :: InputFileName
    character(IntfStrLen)           :: OutRootName
    type(SeaSt_InputType)           :: u           !< An initial guess for the input; input mesh must be defined
    type(SeaSt_ContinuousStateType) :: x           !< Initial continuous states
    type(SeaSt_DiscreteStateType)   :: xd          !< Initial discrete states
    type(SeaSt_ConstraintStateType) :: z           !< Initial guess of the constraint states
    type(SeaSt_OtherStateType)      :: OtherState  !< Initial other states            
    real(DbKi)                      :: Interval    !< Coupling interval in seconds: the rate that 
                                                                   !!   (1) SeaSt_UpdateStates() is called in loose coupling &
                                                                   !!   (2) SeaSt_UpdateDiscState() is called in tight coupling.
                                                                   !!   Input is the suggested time from the glue code; 
                                                                   !!   Output is the actual coupling interval that will be used 
                                                                   !!   by the glue code.

    integer                    :: ErrStat                          !< aggregated error status
    character(ErrMsgLen)       :: ErrMsg                           !< aggregated error message
    integer                    :: ErrStat2                         !< temporary error status  from a call
    character(ErrMsgLen)       :: ErrMsg2                          !< temporary error message from a call
    integer                    :: i,j,k
    character(*), parameter    :: RoutineName = 'SeaSt_C_Init'  !< for error handling

    ! Initialize error handling
    ErrStat  =  ErrID_None
    ErrMsg   =  ""

    CALL NWTC_Init( ProgNameIn=version%Name )
    CALL DispCopyrightLicense( version%Name )
    CALL DispCompileRuntimeInfo( version%Name )

    ! interface debugging
    ! DebugLevel = int(DebugLevel_in,IntKi)

    ! Input files
    call c_f_pointer(InputFile_c, InputFileString)  ! Get a pointer to the input file string
    InputFileName = FileNameFromCString(InputFileString, IntfStrLen)  ! convert the input file name from c_char to fortran character

    call c_f_pointer(OutRootName_c, OutputFileString)  ! Get a pointer to the input file string
    OutRootName = FileNameFromCString(OutputFileString, IntfStrLen)  ! convert the input file name from c_char to fortran character

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

    ! Set other inputs for calling SeaSt_Init
    InitInp%InputFile    = InputFileName
    InitInp%UseInputFile = .TRUE. 
    InitInp%OutRootName  = OutRootName
    InitInp%Gravity      = Gravity_c
    InitInp%defWtrDens   = WtrDens_c
    InitInp%defWtrDpth   = WtrDpth_c
    InitInp%defMSL2SWL   = MSL2SWL_c
    InitInp%TMax         = (NSteps_c - 1) * TimeInterval_c   ! Using this to match the SeaState driver; could otherwise get TMax directly
    InitInp%WaveFieldMod = WaveElevSeriesFlag_c
    ! REAL(ReKi)  :: PtfmLocationX = 0.0_ReKi      !< Supplied by Driver:  X coordinate of platform location in the wave field [m]
    ! REAL(ReKi)  :: PtfmLocationY = 0.0_ReKi      !< Supplied by Driver:  Y coordinate of platform location in the wave field [m]
    InitInp%WrWvKinMod = WrWvKinMod_c
    ! LOGICAL  :: HasIce = .false.      !< Supplied by Driver:  Whether this simulation has ice loading (flag) [-]
    ! LOGICAL  :: Linearize = .FALSE.      !< Flag that tells this module if the glue code wants to linearize. [-]
    ! LOGICAL  :: SurfaceVis = .FALSE.      !< Turn on grid surface visualization outputs [-]
    ! INTEGER(IntKi)  :: SurfaceVisNx = 0      !< Number of points in X direction to output for visualization grid.  Use 0 or negative to set to SeaState resolution. [-]
    ! INTEGER(IntKi)  :: SurfaceVisNy = 0      !< Number of points in Y direction to output for visualization grid.  Use 0 or negative to set to SeaState resolution. [-]

    call SeaSt_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOutData, ErrStat2, ErrMsg2 )
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
        end do
    end do

    ! null terminate the string
    OutputChannelNames_C(k) = C_NULL_CHAR
    OutputChannelUnits_C(k) = C_NULL_CHAR

    call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)

contains
    logical function Failed()
        call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
        Failed = ErrStat >= AbortErrLev
        if (Failed) then
            call Cleanup()
            call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
        endif
    end function Failed

    subroutine Cleanup()    ! NOTE: we are ignoring any error reporting from here
    end subroutine Cleanup 

    subroutine ShowPassedData()
        ! character(1) :: TmpFlag
        ! integer      :: i,j
        call WrSCr("")
        call WrScr("-----------------------------------------------------------")
        call WrScr("Interface debugging:  Variables passed in through interface")
        call WrScr("   SeaSt_C_Init")
        call WrScr("   --------------------------------------------------------")
        call WrScr("   FileInfo")
        ! TmpFlag="F";   if (IfWinputFilePassed==1_c_int) TmpFlag="T"
        ! call WrScr("       IfWinputFilePassed_C           "//TmpFlag )
        ! call WrScr("       IfWinputFileString_C (ptr addr)"//trim(Num2LStr(LOC(IfWinputFileString_C))) )
        ! call WrScr("       IfWinputFileStringLength_C     "//trim(Num2LStr( IfWinputFileStringLength_C )) )
        ! call WrScr("       OutRootName                    "//trim(OutRootName) )
        ! call WrScr("   Input variables")
        ! call WrScr("       NumWindPts_C                   "//trim(Num2LStr( NumWindPts_C)) )
        ! call WrScr("   Time variables")
        ! call WrScr("       DT_C                           "//trim(Num2LStr( DT_C          )) )
        call WrScr("-----------------------------------------------------------")
    end subroutine ShowPassedData
end subroutine SeaSt_C_Init

subroutine SeaSt_C_CalcOutput(Time_C, OutputChannelValues_C, ErrStat_C, ErrMsg_C) BIND (C, NAME='SeaSt_C_CalcOutput')
implicit none
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: SeaSt_C_CalcOutput
!GCC$ ATTRIBUTES DLLEXPORT :: SeaSt_C_CalcOutput
#endif

    real(C_DOUBLE),             intent(in   ) :: Time_C
    real(C_FLOAT),              intent(  out) :: OutputChannelValues_C(p%NumOuts)
    integer(C_INT),             intent(  out) :: ErrStat_C
    character(KIND=C_CHAR),     intent(  out) :: ErrMsg_C(ErrMsgLen_C)

    ! Local variables
    type(SeaSt_InputType)           :: u           !< An initial guess for the input; input mesh must be defined
    type(SeaSt_ContinuousStateType) :: x           !< Initial continuous states
    type(SeaSt_DiscreteStateType)   :: xd          !< Initial discrete states
    type(SeaSt_ConstraintStateType) :: z           !< Initial guess of the constraint states
    type(SeaSt_OtherStateType)      :: OtherState  !< Initial other states            

    real(DbKi)                 :: Time
    integer                    :: ErrStat                          !< aggregated error status
    character(ErrMsgLen)       :: ErrMsg                           !< aggregated error message
    integer                    :: ErrStat2                         !< temporary error status  from a call
    character(ErrMsgLen)       :: ErrMsg2                          !< temporary error message from a call
    character(*), parameter    :: RoutineName = 'SeaSt_C_End'  !< for error handling

    ! Initialize error handling
    ErrStat  =  ErrID_None
    ErrMsg   =  ""

    ! Convert the inputs from C to Fortran
    Time = REAL(Time_C,DbKi)

    call SeaSt_CalcOutput( Time, u, p, x, xd, z, OtherState, y, m, ErrStat2, ErrMsg2 )
        if (Failed()) return

    ! Get the output channel info out of y
    OutputChannelValues_C = REAL(y%WriteOutput, C_FLOAT)

    call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)

contains
    logical function Failed()
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
        Failed = ErrStat >= AbortErrLev
        if (Failed) call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
    end function Failed
end subroutine

subroutine SeaSt_C_End(ErrStat_C,ErrMsg_C) BIND (C, NAME='SeaSt_C_End')
implicit none
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: SeaSt_C_End
!GCC$ ATTRIBUTES DLLEXPORT :: SeaSt_C_End
#endif
    integer(c_int),             intent(  out) :: ErrStat_C
    character(kind=c_char),     intent(  out) :: ErrMsg_C(ErrMsgLen_C)

    ! Local variables
    type(SeaSt_InputType)           :: u           !< An initial guess for the input; input mesh must be defined
    type(SeaSt_ContinuousStateType) :: x           !< Initial continuous states
    type(SeaSt_DiscreteStateType)   :: xd          !< Initial discrete states
    type(SeaSt_ConstraintStateType) :: z           !< Initial guess of the constraint states
    type(SeaSt_OtherStateType)      :: OtherState  !< Initial other states            

    integer                    :: ErrStat                          !< aggregated error status
    character(ErrMsgLen)       :: ErrMsg                           !< aggregated error message
    integer                    :: ErrStat2                         !< temporary error status  from a call
    character(ErrMsgLen)       :: ErrMsg2                          !< temporary error message from a call
    character(*), parameter    :: RoutineName = 'SeaSt_C_End'  !< for error handling

    ! Initialize error handling
    ErrStat  =  ErrID_None
    ErrMsg   =  ""
    call SeaSt_End(u, p, x, xd, z, OtherState, y, m, ErrStat2, ErrMsg2)

    call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
    call SetErr( ErrStat, ErrMsg, ErrStat_C, ErrMsg_C )

end subroutine

! subroutine get_wave_height(position)


! subroutine get_wave_field_pointer()
! pass back the internal pointer to the wave field to the calling code
! end subroutine

! subroutine set_flow_field_pointer()

! end subroutine

end module SeaState_C_Binding
