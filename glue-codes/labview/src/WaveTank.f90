
MODULE WaveTankTesting

    USE ISO_C_BINDING
    USE NWTC_Library
    ! USE Precision
    USE MoorDyn_C
    USE SeaState_C_Binding
    USE NWTC_C_Binding, ONLY: IntfStrLen, SetErr
    USE AeroDyn_Inflow_C_BINDING, ONLY: ADI_C_Init, MaxADIOutputs

    IMPLICIT NONE
    SAVE

    PUBLIC :: WaveTank_Init

    REAL(C_DOUBLE) :: dt_c = 0.01_C_DOUBLE  ! 100 hertz
    REAL(C_FLOAT) :: g_c = 9.8065_C_FLOAT
    REAL(C_FLOAT) :: rho_c = 1025.0_C_FLOAT
    REAL(C_FLOAT) :: depth_c = 200.0_C_FLOAT
    REAL(C_FLOAT), DIMENSION(6) :: ptfminit_c = 0.0_C_FLOAT
    INTEGER(C_INT) :: interporder_c = 2     ! 1: linear (uses two time steps) or 2: quadratic (uses three time steps)
    
    INTEGER(C_INT) :: N_CAMERA_POINTS
    
    INTEGER(C_INT) :: load_period = 20 ! seconds

CONTAINS

SUBROUTINE SetErrStat_C(ErrStatLocal, ErrMessLocal, ErrStatGlobal, ErrMessGlobal, RoutineName)

    INTEGER(C_INT),             INTENT(IN   )   :: ErrStatLocal     ! Error status of the operation
    CHARACTER(KIND=C_CHAR),     INTENT(IN   )   :: ErrMessLocal     ! Error message if ErrStat /= ErrID_None
    INTEGER(C_INT),             INTENT(INOUT)   :: ErrStatGlobal    ! Error status of the operation
    CHARACTER(KIND=C_CHAR),     INTENT(INOUT)   :: ErrMessGlobal    ! Error message if ErrStat /= ErrID_None
    CHARACTER(*),               INTENT(IN   )   :: RoutineName      ! Name of the routine error occurred in

    IF ( ErrStatLocal == ErrID_None ) RETURN

    IF (ErrStatGlobal /= ErrID_None) ErrMessGlobal = TRIM(ErrMessGlobal)//new_line('a')
    ErrMessGlobal = TRIM(ErrMessGlobal)//TRIM(RoutineName)//':'//TRIM(ErrMessLocal)
    ErrStatGlobal = MAX(ErrStatGlobal,ErrStatLocal)

END SUBROUTINE 

SUBROUTINE WaveTank_Init(   &
    MD_InputFile_C,         &
    SS_InputFile_C,         &
    AD_InputFile_C,         &
    IfW_InputFile_C,        &
    n_camera_points_C,      & 
    ErrStat_C,              &
    ErrMsg_C                &
) BIND (C, NAME='WaveTank_Init')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: WaveTank_Init
!GCC$ ATTRIBUTES DLLEXPORT :: WaveTank_Init
#endif

IMPLICIT NONE

    TYPE(C_PTR),        INTENT(IN   ) :: MD_InputFile_C
    TYPE(C_PTR),        INTENT(IN   ) :: SS_InputFile_C
    TYPE(C_PTR),        INTENT(IN   ) :: AD_InputFile_C
    TYPE(C_PTR),        INTENT(IN   ) :: IfW_InputFile_C
    INTEGER(C_INT),     INTENT(IN   ) :: n_camera_points_c
    INTEGER(C_INT),     INTENT(  OUT)         :: ErrStat_C
    CHARACTER(*, KIND=C_CHAR), INTENT(  OUT)  :: ErrMsg_C

    ! Local variables
    INTEGER(C_INT)                          :: ErrStat_C2
    CHARACTER(KIND=C_CHAR, LEN=ErrMsgLen_C) :: ErrMsg_C2
    INTEGER(C_INT)                          :: input_file_passed = 0  ! We're passing paths to input files rather than input strings for all modules

    ! SeaState variables
    CHARACTER(KIND=C_CHAR), TARGET :: SS_OutRootName(IntfStrLen)
    TYPE(C_PTR)    :: SS_OutRootName_C           ! TYPE(C_PTR),    intent(in   )
    REAL(C_FLOAT)  :: SS_Gravity_C               ! REAL(C_FLOAT),  intent(in   )
    REAL(C_FLOAT)  :: SS_WtrDens_C               ! REAL(C_FLOAT),  intent(in   )
    REAL(C_FLOAT)  :: SS_WtrDpth_C               ! REAL(C_FLOAT),  intent(in   )
    REAL(C_FLOAT)  :: SS_MSL2SWL_C               ! REAL(C_FLOAT),  intent(in   )
    INTEGER(C_INT) :: SS_NSteps_C                ! INTEGER(C_INT), intent(in   )
    REAL(C_FLOAT)  :: SS_TimeInterval_C          ! REAL(C_FLOAT),  intent(in   )
    INTEGER(C_INT) :: SS_WaveElevSeriesFlag_C    ! INTEGER(C_INT), intent(in   )
    INTEGER(C_INT) :: SS_WrWvKinMod_C            ! INTEGER(C_INT), intent(in   )
    INTEGER(C_INT) :: SS_NumChannels_C           ! INTEGER(C_INT), intent(  out)
    CHARACTER(KIND=C_CHAR) :: SS_OutputChannelNames_c(ChanLen*MaxOutPts+1)  ! CHARACTER(KIND=C_CHAR), intent(  out)
    CHARACTER(KIND=C_CHAR) :: SS_OutputChannelUnits_c(ChanLen*MaxOutPts+1)  ! CHARACTER(KIND=C_CHAR), intent(  out)

    ! MD variables
    ! CHARACTER(KIND=C_CHAR), TARGET :: MD_InputFileString(IntfStrLen)
    INTEGER(C_INT) :: MD_InputFilePassed            ! INTEGER(C_INT, INTENT(IN   ) :: InputFilePassed        !< Whether to load the file from the filesystem - 1: InputFileString_C contains the contents of the input file; otherwise, InputFileString_C contains the path to the input file
    TYPE(C_PTR)    :: MD_InputFileString_C          ! TYPE(C_PTR),   INTENT(IN   ) :: InputFileString_C        !< Input file as a single string with lines deliniated by C_NULL_CHAR
    INTEGER(C_INT) :: MD_InputFileStringLength_C    ! INTEGER(C_INT, INTENT(IN   ) :: InputFileStringLength_C  !< length of the input file string
    REAL(C_DOUBLE) :: MD_DT_C                       ! REAL(C_DOUBLE, INTENT(IN   ) :: DT_C
    REAL(C_FLOAT)  :: MD_G_C                        ! REAL(C_FLOAT), INTENT(IN   ) :: G_C
    REAL(C_FLOAT)  :: MD_RHO_C                      ! REAL(C_FLOAT), INTENT(IN   ) :: RHO_C
    REAL(C_FLOAT)  :: MD_DEPTH_C                    ! REAL(C_FLOAT), INTENT(IN   ) :: DEPTH_C
    REAL(C_FLOAT)  :: MD_PtfmInit_C(6)              ! REAL(C_FLOAT), INTENT(IN   ) :: PtfmInit_C(6)
    INTEGER(C_INT) :: MD_InterpOrder_C              ! INTEGER(C_INT, INTENT(IN   ) :: InterpOrder_C
    INTEGER(C_INT) :: MD_NumChannels_C              ! INTEGER(C_INT, INTENT(  OUT) :: NumChannels_C
    CHARACTER(KIND=C_CHAR) :: MD_OutputChannelNames_C(100000)  ! CHARACTER(KIND=C_CHAR)                         , INTENT(  OUT)   :: OutputChannelNames_C(100000)
    CHARACTER(KIND=C_CHAR) :: MD_OutputChannelUnits_C(100000)  ! CHARACTER(KIND=C_CHAR)                         , INTENT(  OUT)   :: OutputChannelUnits_C(100000)

    ! ADI variables
    INTEGER(C_INT) :: AD_InputFilePassed             ! intent(in   )  :: ADinputFilePassed                      !< Whether to load the file from the filesystem - 1: ADinputFileString_C contains the contents of the input file; otherwise, ADinputFileString_C contains the path to the input file
    TYPE(C_PTR)    :: AD_InputFileString_C           ! intent(in   )  :: ADinputFileString_C                    !< Input file as a single string with lines delineated by C_NULL_CHAR
    INTEGER(C_INT) :: AD_InputFileStringLength_C     ! intent(in   )  :: ADinputFileStringLength_C              !< length of the input file string
    INTEGER(C_INT) :: IfW_InputFilePassed            ! intent(in   )  :: IfWinputFilePassed                     !< Whether to load the file from the filesystem - 1: IfWinputFileString_C contains the contents of the input file; otherwise, IfWinputFileString_C contains the path to the input file
    TYPE(C_PTR)    :: IfW_InputFileString_C          ! intent(in   )  :: IfWinputFileString_C                   !< Input file as a single string with lines delineated by C_NULL_CHAR
    INTEGER(C_INT) :: IfW_InputFileStringLength_C    ! intent(in   )  :: IfWinputFileStringLength_C             !< length of the input file string
    CHARACTER(KIND=C_CHAR) :: ADI_OutRootName_C(IntfStrLen) ! intent(in   )  :: OutRootName_C(IntfStrLen)              !< Root name to use for echo files and other
    CHARACTER(KIND=C_CHAR) :: ADI_OutVTKDir_C(IntfStrLen)   ! intent(in   )  :: OutVTKDir_C(IntfStrLen)                !< Directory to put all vtk output
    ! Environmental
    REAL(C_FLOAT) :: ADI_gravity_C                   ! intent(in   )  :: gravity_C                              !< Gravitational acceleration (m/s^2)
    REAL(C_FLOAT) :: ADI_defFldDens_C                ! intent(in   )  :: defFldDens_C                           !< Air density (kg/m^3)
    REAL(C_FLOAT) :: ADI_defKinVisc_C                ! intent(in   )  :: defKinVisc_C                           !< Kinematic viscosity of working fluid (m^2/s)
    REAL(C_FLOAT) :: ADI_defSpdSound_C               ! intent(in   )  :: defSpdSound_C                          !< Speed of sound in working fluid (m/s)
    REAL(C_FLOAT) :: ADI_defPatm_C                   ! intent(in   )  :: defPatm_C                              !< Atmospheric pressure (Pa) [used only for an MHK turbine cavitation check]
    REAL(C_FLOAT) :: ADI_defPvap_C                   ! intent(in   )  :: defPvap_C                              !< Vapour pressure of working fluid (Pa) [used only for an MHK turbine cavitation check]
    REAL(C_FLOAT) :: ADI_WtrDpth_C                   ! intent(in   )  :: WtrDpth_C                              !< Water depth (m)
    REAL(C_FLOAT) :: ADI_MSL2SWL_C                   ! intent(in   )  :: MSL2SWL_C                              !< Offset between still-water level and mean sea level (m) [positive upward]
    ! Interpolation
    INTEGER(C_INT) :: ADI_InterpOrder_C              ! intent(in   )  :: InterpOrder_C                          !< Interpolation order to use (must be 1 or 2)
    ! Time
    REAL(C_DOUBLE) :: ADI_DT_C                       ! intent(in   )  :: DT_C                                   !< Timestep used with AD for stepping forward from t to t+dt.  Must be constant.
    REAL(C_DOUBLE) :: ADI_TMax_C                     ! intent(in   )  :: TMax_C                                 !< Maximum time for simulation
    ! Flags
    INTEGER(C_INT) :: ADI_storeHHVel                 ! intent(in   )  :: storeHHVel                             !< Store hub height time series from IfW
    ! VTK
    INTEGER(C_INT) :: ADI_WrVTK_in                   ! intent(in   )  :: WrVTK_in                               !< Write VTK outputs [0: none, 1: init only, 2: animation]
    INTEGER(C_INT) :: ADI_WrVTK_inType               ! intent(in   )  :: WrVTK_inType                           !< Write VTK outputs as [1: surface, 2: lines, 3: both]
    REAL(C_DOUBLE) :: ADI_WrVTK_inDT                 ! intent(in   )  :: WrVTK_inDT                             !< Timestep between VTK writes
    REAL(C_FLOAT)  :: ADI_VTKNacDim_in(6)            ! intent(in   )  :: VTKNacDim_in(6)                        !< Nacelle dimension passed in for VTK surface rendering [0,y0,z0,Lx,Ly,Lz] (m)
    REAL(C_FLOAT)  :: ADI_VTKHubrad_in               ! intent(in   )  :: VTKHubrad_in                           !< Hub radius for VTK surface rendering
    INTEGER(C_INT) :: ADI_wrOuts_C                   ! intent(in   )  :: wrOuts_C                               !< Write ADI output file
    REAL(C_DOUBLE) :: ADI_DT_Outs_C                  ! intent(in   )  :: DT_Outs_C                              !< Timestep to write output file from ADI
    ! Output
    INTEGER(C_INT) :: ADI_NumChannels_C              ! intent(  out)  :: NumChannels_C                          !< Number of output channels requested from the input file
    CHARACTER(KIND=C_CHAR) :: ADI_OutputChannelNames_C(ChanLen*MaxADIOutputs+1) ! intent(  out)  :: OutputChannelNames_C(ChanLen*MaxADIOutputs+1)    !< NOTE: if MaxADIOutputs is sufficiently large, we may overrun the buffer on the Python side.
    CHARACTER(KIND=C_CHAR) :: ADI_OutputChannelUnits_C(ChanLen*MaxADIOutputs+1) ! intent(  out)  :: OutputChannelUnits_C(ChanLen*MaxADIOutputs+1)

    ! Initialize error handling
    ErrStat_C  =  ErrID_None
    ErrMsg_C   =  ""

    ! N_CAMERA_POINTS = n_camera_points_c

    SS_OutRootName = "seastate"
    SS_OutRootName_C = c_loc(SS_OutRootName)
    SS_Gravity_C = 9.80665
    SS_WtrDens_C = 1025
    SS_WtrDpth_C = 200.0
    SS_MSL2SWL_C = 0.0
    SS_NSteps_C = 801
    SS_TimeInterval_C = 0.125
    SS_WaveElevSeriesFlag_C = 0
    SS_WrWvKinMod_C = 0

    CALL SeaSt_C_Init(              &    
        SS_InputFile_C,             & ! TYPE(C_PTR),                intent(in   ) :: InputFile_c
        SS_OutRootName_C,           & ! TYPE(C_PTR),                intent(in   ) :: OutRootName_c
        SS_Gravity_C,               & ! REAL(C_FLOAT),              intent(in   ) :: Gravity_c
        SS_WtrDens_C,               & ! REAL(C_FLOAT),              intent(in   ) :: WtrDens_c
        SS_WtrDpth_C,               & ! REAL(C_FLOAT),              intent(in   ) :: WtrDpth_c
        SS_MSL2SWL_C,               & ! REAL(C_FLOAT),              intent(in   ) :: MSL2SWL_c
        SS_NSteps_C,                & ! INTEGER(C_INT),             intent(in   ) :: NSteps_c
        SS_TimeInterval_C,          & ! REAL(C_FLOAT),              intent(in   ) :: TimeInterval_c
        SS_WaveElevSeriesFlag_C,    & ! INTEGER(C_INT),             intent(in   ) :: WaveElevSeriesFlag_c
        SS_WrWvKinMod_C,            & ! INTEGER(C_INT),             intent(in   ) :: WrWvKinMod_c
        SS_NumChannels_C,           & ! INTEGER(C_INT),             intent(  out) :: NumChannels_c
        SS_OutputChannelNames_C,    & ! CHARACTER(KIND=C_CHAR),     intent(  out) :: OutputChannelNames_c(ChanLen*MaxOutPts+1)
        SS_OutputChannelUnits_C,    & ! CHARACTER(KIND=C_CHAR),     intent(  out) :: OutputChannelUnits_c(ChanLen*MaxOutPts+1)
        ErrStat_C2,                 & ! INTEGER(C_INT),             intent(  out) :: ErrStat_c
        ErrMsg_C2                   & ! CHARACTER(KIND=C_CHAR),     intent(  out) :: ErrMsg_c(ErrMsgLen_c)
    )
    ! CALL SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'SeaSt_C_Init')

    ! SS_GetFlowFieldPtr(SSFFPTR)
    ! MD_SetFlowFieldPtr(SSFFPTR)

    MD_InputFilePassed = 0
    MD_InputFileString_C = MD_InputFile_C
    MD_InputFileStringLength_C = IntfStrLen
    MD_DT_C = 0.125
    MD_G_C = 9.80665
    MD_RHO_C = 1025
    MD_DEPTH_C = 200
    MD_PtfmInit_C = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
    MD_InterpOrder_C = 1  ! 1 - Linear, 2 - quadratic

    CALL MD_C_Init(                 &
        MD_InputFilePassed,         &
        MD_InputFileString_C,       &
        MD_InputFileStringLength_C, &
        MD_DT_C,                    &
        MD_G_C,                     &
        MD_RHO_C,                   &     ! Verify that this matches WtrDens in the pointer set above - this is checking that the pointer has values loaded by the init function
        MD_DEPTH_C,                 &
        MD_PtfmInit_C,              &
        MD_InterpOrder_C,           &
        MD_NumChannels_C,           &
        MD_OutputChannelNames_C,    &
        MD_OutputChannelUnits_C,    &
        ErrStat_C2,                 &
        ErrMsg_C2                   &
    )
    CALL SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'MD_C_Init')
    IF (ErrStat_C >= AbortErrLev) RETURN

    AD_InputFilePassed = 0
    AD_InputFileString_C = AD_InputFile_C
    AD_InputFileStringLength_C = IntfStrLen
    IfW_InputFilePassed = 0
    IfW_InputFileString_C = IfW_InputFile_C
    IfW_InputFileStringLength_C = IntfStrLen
    ADI_OutRootName_C(IntfStrLen) = "adi"
    ADI_OutVTKDir_C(IntfStrLen) = "vtk"
    ADI_gravity_C = 9.80655
    ADI_defFldDens_C = 1025
    ADI_defKinVisc_C = 1.06E-06
    ADI_defSpdSound_C = 1500
    ADI_defPatm_C = 101325
    ADI_defPvap_C = 2500
    ADI_WtrDpth_C = 200
    ADI_MSL2SWL_C = 0.0
    ADI_InterpOrder_C = 1
    ADI_DT_C = 0.125
    ADI_TMax_C = 100.0
    ADI_storeHHVel = 1
    ADI_WrVTK_in = 0
    ADI_WrVTK_inType = 3
    ADI_WrVTK_inDT = 0.125
    ADI_VTKNacDim_in = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
    ADI_VTKHubrad_in = 12.0
    ADI_wrOuts_C = 0
    ADI_DT_Outs_C = 0.125

    call ADI_C_Init(                 &
        AD_InputFilePassed,          &  ! INTEGER(C_INT),            intent(in   )  :: ADinputFilePassed                      !< Write VTK outputs [0: none, 1: init only, 2: animation]
        AD_InputFileString_C,        &  ! TYPE(C_PTR),               intent(in   )  :: ADinputFileString_C                    !< Input file as a single string with lines deliniated by C_NULL_CHAR
        AD_InputFileStringLength_C,  &  ! INTEGER(C_INT),            intent(in   )  :: ADinputFileStringLength_C              !< lenght of the input file string
        IfW_InputFilePassed,         &  ! INTEGER(C_INT),            intent(in   )  :: IfWinputFilePassed                     !< Write VTK outputs [0: none, 1: init only, 2: animation]
        IfW_InputFileString_C,       &  ! TYPE(C_PTR),               intent(in   )  :: IfWinputFileString_C                   !< Input file as a single string with lines deliniated by C_NULL_CHAR
        IfW_InputFileStringLength_C, &  ! INTEGER(C_INT),            intent(in   )  :: IfWinputFileStringLength_C             !< lenght of the input file string
        ADI_OutRootName_C,           &  ! CHARACTER(KIND=C_CHAR),    intent(in   )  :: OutRootName_C(IntfStrLen)              !< Root name to use for echo files and other
        ADI_OutVTKDir_C,             &  ! CHARACTER(KIND=C_CHAR),    intent(in   )  :: OutVTKDir_C(IntfStrLen)                !< Directory to put all vtk output
        ADI_gravity_C,               &  ! REAL(C_FLOAT),             intent(in   )  :: gravity_C                              !< Gravitational acceleration (m/s^2)
        ADI_defFldDens_C,            &  ! REAL(C_FLOAT),             intent(in   )  :: defFldDens_C                           !< Air density (kg/m^3)
        ADI_defKinVisc_C,            &  ! REAL(C_FLOAT),             intent(in   )  :: defKinVisc_C                           !< Kinematic viscosity of working fluid (m^2/s)
        ADI_defSpdSound_C,           &  ! REAL(C_FLOAT),             intent(in   )  :: defSpdSound_C                          !< Speed of sound in working fluid (m/s)
        ADI_defPatm_C,               &  ! REAL(C_FLOAT),             intent(in   )  :: defPatm_C                              !< Atmospheric pressure (Pa) [used only for an MHK turbine cavitation check]
        ADI_defPvap_C,               &  ! REAL(C_FLOAT),             intent(in   )  :: defPvap_C                              !< Vapour pressure of working fluid (Pa) [used only for an MHK turbine cavitation check]
        ADI_WtrDpth_C,               &  ! REAL(C_FLOAT),             intent(in   )  :: WtrDpth_C                              !< Water depth (m)
        ADI_MSL2SWL_C,               &  ! REAL(C_FLOAT),             intent(in   )  :: MSL2SWL_C                              !< Offset between still-water level and mean sea level (m) [positive upward]
        ADI_InterpOrder_C,           &  ! INTEGER(C_INT),            intent(in   )  :: InterpOrder_C                          !< Interpolation order to use (must be 1 or 2)
        ADI_DT_C,                    &  ! REAL(C_DOUBLE),            intent(in   )  :: DT_C                                   !< Timestep used with AD for stepping forward from t to t+dt.  Must be constant.
        ADI_TMax_C,                  &  ! REAL(C_DOUBLE),            intent(in   )  :: TMax_C                                 !< Maximum time for simulation
        ADI_storeHHVel,              &  ! INTEGER(C_INT),            intent(in   )  :: storeHHVel                             !< Store hub height time series from IfW
        ADI_WrVTK_in,                &  ! INTEGER(C_INT), intent(in   ) :: WrVTK_in        !< Write VTK outputs [0: none, 1: init only, 2: animation]
        ADI_WrVTK_inType,            &  ! INTEGER(C_INT), intent(in   ) :: WrVTK_inType    !< Write VTK outputs as [1: surface, 2: lines, 3: both]
        ADI_WrVTK_inDT,              &  ! REAL(C_DOUBLE), intent(in   ) :: WrVTK_inDT      !< Timestep between VTK writes
        ADI_VTKNacDim_in,            &  ! REAL(C_FLOAT),  intent(in   ) :: VTKNacDim_in(6) !< Nacelle dimension passed in for VTK surface rendering [0,y0,z0,Lx,Ly,Lz] (m)
        ADI_VTKHubRad_in,            &  ! REAL(C_FLOAT),  intent(in   ) :: VTKHubrad_in    !< Hub radius for VTK surface rendering
        ADI_wrOuts_C,                &
        ADI_DT_Outs_C,               &
        ADI_NumChannels_C,           &
        ADI_OutputChannelNames_C,    &
        ADI_OutputChannelUnits_C,    &
        ErrStat_C2,                  &
        ErrMsg_C2                    &
    )
    CALL SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'MD_C_Init')
    IF (ErrStat_C >= AbortErrLev) RETURN

! Set compiler flag for Labview
! Cmpl4LV   = .TRUE.

END SUBROUTINE WaveTank_Init

! delta_time,                 &
SUBROUTINE WaveTank_CalcOutput( &
    frame_number,               &
    positions_x,                &
    positions_y,                &
    positions_z,                &
    rotation_matrix,            &
    loads,                      &
    ErrStat_C,                  &
    ErrMsg_C                    &
) BIND (C, NAME='WaveTank_CalcOutput')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: WaveTank_CalcOutput
!GCC$ ATTRIBUTES DLLEXPORT :: WaveTank_CalcOutput
#endif

IMPLICIT NONE

! INTEGER(C_INT)                        :: delta_time
INTEGER(C_INT)                        :: frame_number
REAL(C_FLOAT),          INTENT(IN   ) :: positions_x(N_CAMERA_POINTS)
REAL(C_FLOAT),          INTENT(IN   ) :: positions_y(N_CAMERA_POINTS)
REAL(C_FLOAT),          INTENT(IN   ) :: positions_z(N_CAMERA_POINTS)
REAL(C_FLOAT),          INTENT(IN   ) :: rotation_matrix(9)
REAL(C_FLOAT),          INTENT(  OUT) :: loads(N_CAMERA_POINTS)
INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_C
CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_C(ErrMsgLen_C)

INTEGER :: i

IF ( MOD(frame_number / load_period, 2) == 0 ) THEN
    loads = -1.0
ELSE
    loads = 1.0
ENDIF

END SUBROUTINE

SUBROUTINE WaveTank_End() bind (C, NAME="WaveTank_End")


IMPLICIT NONE




END SUBROUTINE

end module WaveTankTesting
