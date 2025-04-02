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
    TYPE(ProgDesc), PARAMETER              :: version   = SeaSt_ProgDesc

    !------------------------------------------------------------------------------------
    !  Debugging: DebugLevel -- passed at PreInit
    !     0  - none
    !     1  - some summary info
    !     2  - above + all position/orientation info
    !     3  - above + input files (if direct passed)
    !     4  - above + meshes
    INTEGER(IntKi)                         :: DebugLevel = 0

    !------------------------------
    !  Primary derived types
    TYPE(SeaSt_InputType)                    :: InputData         !< Inputs to SeaState
    TYPE(SeaSt_InitInputType)                :: InitInp
    TYPE(SeaSt_InitOutputType)               :: InitOutData       !< Initial output data -- Names, units, and version info.
    TYPE(SeaSt_ParameterType)                :: p                 !< Parameters
    TYPE(SeaSt_OutputType)                   :: y                 !< Initial output (outputs are not calculated; only the output mesh is initialized)
    TYPE(SeaSt_MiscVarType)                  :: m                 !< Misc variables for optimization (not copied in glue code)

CONTAINS


SUBROUTINE SeaSt_C_Init(InputFile_C, OutRootName_C, Gravity_C, WtrDens_C, WtrDpth_C, MSL2SWL_C, NSteps_C, TimeInterval_C, WaveElevSeriesFlag_C, WrWvKinMod_C, NumChannels_C, OutputChannelNames_C, OutputChannelUnits_C, ErrStat_C, ErrMsg_C) BIND (C, NAME='SeaSt_C_Init')
IMPLICIT NONE
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: SeaSt_C_Init
!GCC$ ATTRIBUTES DLLEXPORT :: SeaSt_C_Init
#endif
    TYPE(C_PTR),                INTENT(IN   ) :: InputFile_C
    TYPE(C_PTR),                INTENT(IN   ) :: OutRootName_C
    REAL(C_FLOAT),              INTENT(IN   ) :: Gravity_C
    REAL(C_FLOAT),              INTENT(IN   ) :: WtrDens_C
    REAL(C_FLOAT),              INTENT(IN   ) :: WtrDpth_C
    REAL(C_FLOAT),              INTENT(IN   ) :: MSL2SWL_C
    INTEGER(C_INT),             INTENT(IN   ) :: NSteps_C
    REAL(C_FLOAT),              INTENT(IN   ) :: TimeInterval_C
    INTEGER(C_INT),             INTENT(IN   ) :: WaveElevSeriesFlag_C
    INTEGER(C_INT),             INTENT(IN   ) :: WrWvKinMod_C
    INTEGER(C_INT),             INTENT(  OUT) :: NumChannels_C
    CHARACTER(KIND=C_CHAR),     INTENT(  OUT) :: OutputChannelNames_C(ChanLen*MaxOutPts+1)
    CHARACTER(KIND=C_CHAR),     INTENT(  OUT) :: OutputChannelUnits_C(ChanLen*MaxOutPts+1)
    INTEGER(C_INT),             INTENT(  OUT) :: ErrStat_C
    CHARACTER(KIND=C_CHAR),     INTENT(  OUT) :: ErrMsg_C(ErrMsgLen_C)

    ! Local variables
    CHARACTER(KIND=C_CHAR, len=IntfStrLen), POINTER :: InputFileString          !< Input file as a single string with NULL chracter separating lines
    CHARACTER(KIND=C_CHAR, len=IntfStrLen), POINTER :: OutputFileString          !< Input file as a single string with NULL chracter separating lines
    CHARACTER(IntfStrLen)           :: InputFileName
    CHARACTER(IntfStrLen)           :: OutRootName
    TYPE(SeaSt_InputType)           :: u           !< An initial guess for the input; input mesh must be defined
    TYPE(SeaSt_ContinuousStateType) :: x           !< Initial continuous states
    TYPE(SeaSt_DiscreteStateType)   :: xd          !< Initial discrete states
    TYPE(SeaSt_ConstraintStateType) :: z           !< Initial guess of the constraint states
    TYPE(SeaSt_OtherStateType)      :: OtherState  !< Initial other states            
    REAL(DbKi)                      :: Interval    !< Coupling interval in seconds: the rate that 
                                                                   !!   (1) SeaSt_UpdateStates() is called in loose coupling &
                                                                   !!   (2) SeaSt_UpdateDiscState() is called in tight coupling.
                                                                   !!   Input is the suggested time from the glue code; 
                                                                   !!   Output is the actual coupling interval that will be used 
                                                                   !!   by the glue code.

    INTEGER                    :: ErrStat                          !< aggregated error status
    CHARACTER(ErrMsgLen)       :: ErrMsg                           !< aggregated error message
    INTEGER                    :: ErrStat2                         !< temporary error status  from a call
    CHARACTER(ErrMsgLen)       :: ErrMsg2                          !< temporary error message from a call
    INTEGER                    :: i,j,k
    CHARACTER(*), PARAMETER    :: RoutineName = 'SeaSt_C_Init'  !< for error handling

    ! Initialize error handling
    ErrStat  =  ErrID_None
    ErrMsg   =  ""

    CALL NWTC_Init( ProgNameIn=version%Name )
    CALL DispCopyrightLicense( version%Name )
    CALL DispCompileRuntimeInfo( version%Name )

    ! interface debugging
    ! DebugLevel = int(DebugLevel_in,IntKi)

    ! Input files
    CALL C_F_POINTER(InputFile_C, InputFileString)  ! Get a pointer to the input file string
    InputFileName = FileNameFromCString(InputFileString, IntfStrLen)  ! convert the input file name from c_char to fortran character

    CALL C_F_POINTER(OutRootName_C, OutputFileString)  ! Get a pointer to the input file string
    OutRootName = FileNameFromCString(OutputFileString, IntfStrLen)  ! convert the input file name from c_char to fortran character

    ! if non-zero, show all passed data here.  Then check valid values
    IF (DebugLevel /= 0_IntKi) THEN
        CALL WrScr("   Interface debugging level "//trim(Num2Lstr(DebugLevel))//" requested.")
        CALL ShowPassedData()
    ENDIF
    ! check valid debug level
    IF (DebugLevel < 0_IntKi) THEN
        ErrStat2 = ErrID_Fatal
        ErrMsg2  = "Interface debug level must be 0 or greater"//NewLine// &
        "  0  - none"//NewLine// &
        "  1  - some summary info and variables passed through interface"//NewLine// &
        "  2  - above + all position/orientation info"//NewLine// &
        "  3  - above + input files (if direct passed)"//NewLine// &
        "  4  - above + meshes"
        IF (Failed()) RETURN;
    ENDIF

    ! For debugging the interface:
    IF (DebugLevel > 0) THEN
        CALL ShowPassedData()
    ENDIF

    ! Set other inputs for calling SeaSt_Init
    InitInp%InputFile    = InputFileName
    InitInp%UseInputFile = .TRUE. 
    InitInp%OutRootName  = OutRootName
    InitInp%Gravity      = Gravity_C
    InitInp%defWtrDens   = WtrDens_C
    InitInp%defWtrDpth   = WtrDpth_C
    InitInp%defMSL2SWL   = MSL2SWL_C
    InitInp%TMax         = (NSteps_C - 1) * TimeInterval_C   ! Using this to match the SeaState driver; could otherwise get TMax directly
    InitInp%WaveFieldMod = WaveElevSeriesFlag_C
    ! REAL(ReKi)  :: PtfmLocationX = 0.0_ReKi      !< Supplied by Driver:  X coordinate of platform location in the wave field [m]
    ! REAL(ReKi)  :: PtfmLocationY = 0.0_ReKi      !< Supplied by Driver:  Y coordinate of platform location in the wave field [m]
    InitInp%WrWvKinMod = WrWvKinMod_C
    ! LOGICAL  :: HasIce = .false.      !< Supplied by Driver:  Whether this simulation has ice loading (flag) [-]
    ! LOGICAL  :: Linearize = .FALSE.      !< Flag that tells this module if the glue code wants to linearize. [-]
    ! LOGICAL  :: SurfaceVis = .FALSE.      !< Turn on grid surface visualization outputs [-]
    ! INTEGER(IntKi)  :: SurfaceVisNx = 0      !< Number of points in X direction to output for visualization grid.  Use 0 or negative to set to SeaState resolution. [-]
    ! INTEGER(IntKi)  :: SurfaceVisNy = 0      !< Number of points in Y direction to output for visualization grid.  Use 0 or negative to set to SeaState resolution. [-]

    CALL SeaSt_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOutData, ErrStat2, ErrMsg2 )
        IF (Failed()) RETURN

    ! Number of channels
    NumChannels_C = size(InitOutData%WriteOutputHdr)

    ! transfer the output channel names and units to c_char arrays for returning
    k=1
    DO i=1,NumChannels_C
        DO j=1,ChanLen    ! max length of channel name.  Same for units
            OutputChannelNames_C(k)=InitOutData%WriteOutputHdr(i)(j:j)
            OutputChannelUnits_C(k)=InitOutData%WriteOutputUnt(i)(j:j)
            k=k+1
        ENDDO
    ENDDO

    ! null terminate the string
    OutputChannelNames_C(k) = C_NULL_CHAR
    OutputChannelUnits_C(k) = C_NULL_CHAR

    CALL SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)

CONTAINS
    LOGICAL FUNCTION Failed()
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
        Failed = ErrStat >= AbortErrLev
        IF (Failed) THEN
            CALL Cleanup()
            CALL SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
        ENDIF
    END FUNCTION Failed

    SUBROUTINE Cleanup()    ! NOTE: we are ignoring any error reporting from here
    END SUBROUTINE Cleanup 

    SUBROUTINE ShowPassedData()
        ! CHARACTER(1) :: TmpFlag
        ! integer      :: i,j
        CALL WrSCr("")
        CALL WrScr("-----------------------------------------------------------")
        CALL WrScr("Interface debugging:  Variables passed in through interface")
        CALL WrScr("   SeaSt_C_Init")
        CALL WrScr("   --------------------------------------------------------")
        CALL WrScr("   FileInfo")
        CALL WrScr("-----------------------------------------------------------")
    END SUBROUTINE ShowPassedData
END SUBROUTINE SeaSt_C_Init

SUBROUTINE SeaSt_C_CalcOutput(Time_C, OutputChannelValues_C, ErrStat_C, ErrMsg_C) BIND (C, NAME='SeaSt_C_CalcOutput')
IMPLICIT NONE
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: SeaSt_C_CalcOutput
!GCC$ ATTRIBUTES DLLEXPORT :: SeaSt_C_CalcOutput
#endif

    REAL(C_DOUBLE),             INTENT(IN   ) :: Time_C
    REAL(C_FLOAT),              INTENT(  OUT) :: OutputChannelValues_C(p%NumOuts)
    INTEGER(C_INT),             INTENT(  OUT) :: ErrStat_C
    CHARACTER(KIND=C_CHAR),     INTENT(  OUT) :: ErrMsg_C(ErrMsgLen_C)

    ! Local variables
    TYPE(SeaSt_InputType)           :: u           !< An initial guess for the input; input mesh must be defined
    TYPE(SeaSt_ContinuousStateType) :: x           !< Initial continuous states
    TYPE(SeaSt_DiscreteStateType)   :: xd          !< Initial discrete states
    TYPE(SeaSt_ConstraintStateType) :: z           !< Initial guess of the constraint states
    TYPE(SeaSt_OtherStateType)      :: OtherState  !< Initial other states            

    REAL(DbKi)                 :: Time
    INTEGER                    :: ErrStat                          !< aggregated error status
    CHARACTER(ErrMsgLen)       :: ErrMsg                           !< aggregated error message
    INTEGER                    :: ErrStat2                         !< temporary error status  from a call
    CHARACTER(ErrMsgLen)       :: ErrMsg2                          !< temporary error message from a call
    CHARACTER(*), PARAMETER    :: RoutineName = 'SeaSt_C_End'  !< for error handling

    ! Initialize error handling
    ErrStat  =  ErrID_None
    ErrMsg   =  ""

    ! Convert the inputs from C to Fortran
    Time = REAL(Time_C,DbKi)

    CALL SeaSt_CalcOutput( Time, u, p, x, xd, z, OtherState, y, m, ErrStat2, ErrMsg2 )
        IF (Failed()) RETURN

    ! Get the output channel info out of y
    OutputChannelValues_C = REAL(y%WriteOutput, C_FLOAT)

    CALL SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)

CONTAINS
    LOGICAL FUNCTION Failed()
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
        Failed = ErrStat >= AbortErrLev
        IF (Failed) CALL SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
    END FUNCTION Failed
END SUBROUTINE

SUBROUTINE SeaSt_C_End(ErrStat_C,ErrMsg_C) BIND (C, NAME='SeaSt_C_End')
IMPLICIT NONE
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: SeaSt_C_End
!GCC$ ATTRIBUTES DLLEXPORT :: SeaSt_C_End
#endif
    INTEGER(C_INT),             INTENT(  OUT) :: ErrStat_C
    CHARACTER(KIND=C_CHAR),     INTENT(  OUT) :: ErrMsg_C(ErrMsgLen_C)

    ! Local variables
    TYPE(SeaSt_InputType)           :: u           !< An initial guess for the input; input mesh must be defined
    TYPE(SeaSt_ContinuousStateType) :: x           !< Initial continuous states
    TYPE(SeaSt_DiscreteStateType)   :: xd          !< Initial discrete states
    TYPE(SeaSt_ConstraintStateType) :: z           !< Initial guess of the constraint states
    TYPE(SeaSt_OtherStateType)      :: OtherState  !< Initial other states            

    INTEGER                    :: ErrStat                          !< aggregated error status
    CHARACTER(ErrMsgLen)       :: ErrMsg                           !< aggregated error message
    INTEGER                    :: ErrStat2                         !< temporary error status  from a call
    CHARACTER(ErrMsgLen)       :: ErrMsg2                          !< temporary error message from a call
    CHARACTER(*), PARAMETER    :: RoutineName = 'SeaSt_C_End'  !< for error handling

    ! Initialize error handling
    ErrStat  =  ErrID_None
    ErrMsg   =  ""
    CALL SeaSt_End(u, p, x, xd, z, OtherState, y, m, ErrStat2, ErrMsg2)

    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
    CALL SetErr( ErrStat, ErrMsg, ErrStat_C, ErrMsg_C )

END SUBROUTINE

! SUBROUTINE get_wave_height(position)


! SUBROUTINE get_wave_field_pointer()
! pass back the internal pointer to the wave field to the calling code
! END SUBROUTINE

! SUBROUTINE set_flow_field_pointer()

! END SUBROUTINE

end module SeaState_C_Binding
