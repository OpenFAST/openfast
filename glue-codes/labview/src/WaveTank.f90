
MODULE WaveTankTesting

    USE ISO_C_BINDING
    USE NWTC_Library
    USE NWTC_IO
    USE SeaState_C_Binding, ONLY: SeaSt_C_Init, SeaSt_C_CalcOutput, SeaSt_C_End, MaxOutPts, SeaSt_GetWaveFieldPointer_C
    USE SeaSt_WaveField_Types, ONLY: SeaSt_WaveFieldType
    USE AeroDyn_Inflow_C_BINDING, ONLY: ADI_C_PreInit, ADI_C_SetupRotor, ADI_C_Init, ADI_C_End, MaxADIOutputs, ADI_C_SetRotorMotion, ADI_C_UpdateStates, ADI_C_CalcOutput
    USE MoorDyn_C, ONLY: MD_C_Init, MD_C_End, MD_C_SetWaveFieldData, MD_C_UpdateStates, MD_C_CalcOutput
    USE NWTC_C_Binding, ONLY: IntfStrLen, SetErrStat_C, SetErrStat_F2C, ErrMsgLen_C, StringConvert_F2C, FileNameFromCString

    IMPLICIT NONE
    SAVE

    PUBLIC :: WaveTank_Init
    PUBLIC :: WaveTank_CalcOutput
    PUBLIC :: WaveTank_StepFunction
    PUBLIC :: WaveTank_End
    PUBLIC :: WaveTank_SetWaveFieldPointer
    PUBLIC :: WaveTank_NoOp
    PUBLIC :: WaveTank_Sizes

    INTEGER(C_INT) :: SS_NumChannels_C
    INTEGER(C_INT) :: MD_NumChannels_C
    INTEGER(C_INT) :: ADI_NumChannels_C

    REAL(C_FLOAT), DIMENSION(3,6) :: Positions = 0.0_C_FLOAT
    REAL(C_FLOAT), DIMENSION(2,6) :: Velocities = 0.0_C_FLOAT
    REAL(C_FLOAT), DIMENSION(1,6) :: Accelerations = 0.0_C_FLOAT

    TYPE WaveTank_InitInput
        ! SeaState variables
        TYPE(C_PTR)     :: SS_OutRootName_C
        REAL(C_FLOAT)   :: SS_Gravity_C
        REAL(C_FLOAT)   :: SS_WtrDens_C
        REAL(C_FLOAT)   :: SS_WtrDpth_C
        REAL(C_FLOAT)   :: SS_MSL2SWL_C
        INTEGER(C_INT)  :: SS_NSteps_C
        REAL(C_FLOAT)   :: SS_TimeInterval_C
        INTEGER(C_INT)  :: SS_WaveElevSeriesFlag_C
        INTEGER(C_INT)  :: SS_WrWvKinMod_C

        ! MD variables
        REAL(C_DOUBLE)  :: MD_DT_C
        REAL(C_FLOAT)   :: MD_G_C
        REAL(C_FLOAT)   :: MD_RHO_C
        REAL(C_FLOAT)   :: MD_DEPTH_C
        REAL(C_FLOAT)   :: MD_PtfmInit_C(6)
        INTEGER(C_INT)  :: MD_InterpOrder_C

        ! ADI variables
        ! Preinit
        INTEGER(C_INT)  :: NumTurbines_C
        INTEGER(C_INT)  :: TransposeDCM
        INTEGER(C_INT)  :: PointLoadOutput
        REAL(C_FLOAT)   :: ADI_gravity_C
        REAL(C_FLOAT)   :: ADI_defFldDens_C
        REAL(C_FLOAT)   :: ADI_defKinVisc_C
        REAL(C_FLOAT)   :: ADI_defSpdSound_C
        REAL(C_FLOAT)   :: ADI_defPatm_C
        REAL(C_FLOAT)   :: ADI_defPvap_C
        REAL(C_FLOAT)   :: ADI_WtrDpth_C
        REAL(C_FLOAT)   :: ADI_MSL2SWL_C
        INTEGER(C_INT)  :: MHK
        INTEGER(C_INT)  :: DebugLevel
        ! SetupRotor
        INTEGER(C_INT)  :: iWT_c                                    !< Wind turbine / rotor number
        INTEGER(C_INT)  :: TurbineIsHAWT_c                          !< true for HAWT, false for VAWT
        REAL(C_FLOAT)   :: TurbOrigin_C(3)                          !< turbine origin (tower base). Gets added to all meshes to shift turbine position.
        REAL(C_FLOAT)   :: HubPos_C( 3 )                            !< Hub position
        REAL(C_DOUBLE)  :: HubOri_C( 9 )                            !< Hub orientation
        REAL(C_FLOAT)   :: NacPos_C( 3 )                            !< Nacelle position
        REAL(C_DOUBLE)  :: NacOri_C( 9 )                            !< Nacelle orientation
        INTEGER(C_INT)  :: NumBlades_C
        REAL(C_FLOAT), DIMENSION(:), ALLOCATABLE :: BldRootPos_C    !< Blade root positions; 3xNumBlades_C
        REAL(C_DOUBLE), DIMENSION(:), ALLOCATABLE :: BldRootOri_C   !< Blade root orientations; 9xNumBlades_C
        ! Initial nodes
        INTEGER(C_INT)  :: NumMeshPts_C                             !< Number of mesh points we are transferring motions and outputting loads to
        REAL(C_FLOAT), DIMENSION(:), ALLOCATABLE :: InitMeshPos_C   !< A 3xNumMeshPts_C array [x,y,z]
        REAL(C_DOUBLE), DIMENSION(:), ALLOCATABLE :: InitMeshOri_C  !< A 9xNumMeshPts_C array [r11,r12,r13,r21,r22,r23,r31,r32,r33]
        INTEGER(C_INT), DIMENSION(:), ALLOCATABLE :: MeshPtToBladeNum_C !< A NumMeshPts_C array of blade numbers associated with each mesh point
        ! Init
        CHARACTER(KIND=C_CHAR) :: ADI_OutRootName_C(IntfStrLen)     !< Root name to use for echo files and other
        CHARACTER(KIND=C_CHAR) :: ADI_OutVTKDir_C(IntfStrLen)       !< Directory to put all vtk output
        ! Interpolation
        INTEGER(C_INT) :: ADI_InterpOrder_C                         !< Interpolation order to use (must be 1 or 2)
        ! Time
        REAL(C_DOUBLE) :: ADI_DT_C                                  !< Timestep used with AD for stepping forward from t to t+dt.  Must be constant.
        REAL(C_DOUBLE) :: ADI_TMax_C                                !< Maximum time for simulation
        ! Flags
        INTEGER(C_INT) :: ADI_storeHHVel                            !< Store hub height time series from IfW
        ! VTK
        INTEGER(C_INT) :: ADI_WrVTK_in
        INTEGER(C_INT) :: ADI_WrVTK_inType
        REAL(C_DOUBLE) :: ADI_WrVTK_inDT
        REAL(C_FLOAT)  :: ADI_VTKNacDim_in(6)
        REAL(C_FLOAT)  :: ADI_VTKHubrad_in
        INTEGER(C_INT) :: ADI_wrOuts_C
        REAL(C_DOUBLE) :: ADI_DT_Outs_C
    END TYPE WaveTank_InitInput

CONTAINS

SUBROUTINE ReadInput(InputFilePath, InitInp, ErrStat, ErrMsg)

    CHARACTER(*),               INTENT(IN   )   :: InputFilePath
    TYPE(WaveTank_InitInput),   INTENT(  OUT)   :: InitInp
    INTEGER(IntKi),             INTENT(  OUT)   :: ErrStat
    CHARACTER(*),               INTENT(  OUT)   :: ErrMsg

    ! Local variables  
    INTEGER :: UnIn = -1
    ! CHARACTER(1024)                                  :: Line                 ! String to temporarially hold value of read line   
    CHARACTER(1024), target                     :: TmpPath
    CHARACTER(1024), pointer                    :: TmpPointer
    ! CHARACTER(1024)                                  :: TmpFmt               ! Temporary storage for format statement
    CHARACTER(1024)                             :: FileName

    integer(IntKi)                              :: ErrStat2             ! local status of error message
    character(ErrMsgLen)                        :: ErrMsg2              ! local error message if errStat /= ErrID_None

    ErrStat = ErrID_None
    ErrMsg  = " "

    FileName = TRIM(InputFilePath)
    CALL GetNewUnit( UnIn )
    CALL OpenFInpFile( UnIn, FileName, ErrStat2, ErrMsg2)
    CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    IF (ErrStat >= AbortErrLev) RETURN

    CALL ReadCom( UnIn, FileName, 'Init comment', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadCom( UnIn, FileName, 'SeaState Init comment', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    IF (ErrStat >= AbortErrLev) RETURN

    CALL ReadVar( UnIn, FileName, TmpPath, 'SS_OutRootName_C', 'SS_OutRootName_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    TmpPointer => TmpPath
    InitInp%SS_OutRootName_C = C_LOC(TmpPointer)
    CALL ReadVar( UnIn, FileName, InitInp%SS_Gravity_C, 'SS_Gravity_C', 'SS_Gravity_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%SS_WtrDens_C, 'SS_WtrDens_C', 'SS_WtrDens_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%SS_WtrDpth_C, 'SS_WtrDpth_C', 'SS_WtrDpth_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%SS_MSL2SWL_C, 'SS_MSL2SWL_C', 'SS_MSL2SWL_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%SS_NSteps_C, 'SS_NSteps_C', 'SS_NSteps_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%SS_TimeInterval_C, 'SS_TimeInterval_C', 'SS_TimeInterval_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%SS_WaveElevSeriesFlag_C, 'SS_WaveElevSeriesFlag_C', 'SS_WaveElevSeriesFlag_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%SS_WrWvKinMod_C, 'SS_WrWvKinMod_C', 'SS_WrWvKinMod_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    IF (ErrStat >= AbortErrLev) RETURN

    CALL ReadCom( UnIn, FileName, 'MoorDyn Init comment', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%MD_DT_C, 'MD_DT_C', 'MD_DT_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%MD_G_C, 'MD_G_C', 'MD_G_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%MD_RHO_C, 'MD_RHO_C', 'MD_RHO_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%MD_DEPTH_C, 'MD_DEPTH_C', 'MD_DEPTH_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadAry( UnIn, FileName, InitInp%MD_PtfmInit_C, 6, 'MD_PtfmInit_C', 'MD_PtfmInit_C', ErrStat2,  ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%MD_InterpOrder_C, 'MD_InterpOrder_C', 'MD_InterpOrder_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    IF (ErrStat >= AbortErrLev) RETURN

    CALL ReadCom( UnIn, FileName, 'ADI Init comment', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%NumTurbines_C, 'NumTurbines_C', 'NumTurbines_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%TransposeDCM, 'TransposeDCM', 'TransposeDCM', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%PointLoadOutput, 'PointLoadOutput', 'PointLoadOutput', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%ADI_gravity_C, 'ADI_gravity_C', 'ADI_gravity_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%ADI_defFldDens_C, 'ADI_defFldDens_C', 'ADI_defFldDens_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%ADI_defKinVisc_C, 'ADI_defKinVisc_C', 'ADI_defKinVisc_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%ADI_defSpdSound_C, 'ADI_defSpdSound_C', 'ADI_defSpdSound_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%ADI_defPatm_C, 'ADI_defPatm_C', 'ADI_defPatm_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%ADI_defPvap_C, 'ADI_defPvap_C', 'ADI_defPvap_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%ADI_WtrDpth_C, 'ADI_WtrDpth_C', 'ADI_WtrDpth_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%ADI_MSL2SWL_C, 'ADI_MSL2SWL_C', 'ADI_MSL2SWL_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%MHK, 'MHK', 'MHK', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%DebugLevel, 'DebugLevel', 'DebugLevel', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    IF (ErrStat >= AbortErrLev) RETURN

    CALL ReadVar( UnIn, FileName, InitInp%iWT_c, 'iWT_c', 'iWT_c', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%TurbineIsHAWT_c, 'TurbineIsHAWT_c', 'TurbineIsHAWT_c', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadAry( UnIn, FileName, InitInp%TurbOrigin_C, 3, 'TurbOrigin_C', 'TurbOrigin_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadAry( UnIn, FileName, InitInp%HubPos_C, 3, 'HubPos_C', 'HubPos_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadAry( UnIn, FileName, InitInp%HubOri_C, 9, 'HubOri_C', 'HubOri_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadAry( UnIn, FileName, InitInp%NacPos_C, 3, 'NacPos_C', 'NacPos_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadAry( UnIn, FileName, InitInp%NacOri_C, 9, 'NacOri_C', 'NacOri_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%NumBlades_C, 'NumBlades_C', 'NumBlades_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL AllocAry(InitInp%BldRootPos_C, 3*InitInp%NumBlades_C, 'BldRootPos_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL AllocAry(InitInp%BldRootOri_C, 9*InitInp%NumBlades_C, 'BldRootPos_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadAry( UnIn, FileName, InitInp%BldRootPos_C, 3*InitInp%NumBlades_C, 'BldRootPos_C', 'BldRootPos_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadAry( UnIn, FileName, InitInp%BldRootOri_C, 9*InitInp%NumBlades_C, 'BldRootOri_C', 'BldRootOri_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%NumMeshPts_C, 'NumMeshPts_C', 'NumMeshPts_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL AllocAry(InitInp%InitMeshPos_C, 3*InitInp%NumMeshPts_C, 'InitMeshPos_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL AllocAry(InitInp%InitMeshOri_C, 9*InitInp%NumMeshPts_C, 'InitMeshOri_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL AllocAry(InitInp%MeshPtToBladeNum_C, InitInp%NumMeshPts_C, 'MeshPtToBladeNum_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadAry( UnIn, FileName, InitInp%InitMeshPos_C, 3*InitInp%NumMeshPts_C, 'InitMeshPos_C', 'InitMeshPos_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadAry( UnIn, FileName, InitInp%InitMeshOri_C, 9*InitInp%NumMeshPts_C, 'InitMeshOri_C', 'InitMeshOri_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadAry( UnIn, FileName, InitInp%MeshPtToBladeNum_C, InitInp%NumMeshPts_C, 'MeshPtToBladeNum_C', 'MeshPtToBladeNum_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    IF (ErrStat >= AbortErrLev) RETURN

    CALL ReadVar( UnIn, FileName, TmpPath, 'ADI_OutRootName_C', 'ADI_OutRootName_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL StringConvert_F2C(TmpPath, InitInp%ADI_OutRootName_C)
    CALL ReadVar( UnIn, FileName, TmpPath, 'ADI_OutVTKDir_C', 'ADI_OutVTKDir_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL StringConvert_F2C(TmpPath, InitInp%ADI_OutVTKDir_C)
    CALL ReadVar( UnIn, FileName, InitInp%ADI_InterpOrder_C, 'ADI_InterpOrder_C', 'ADI_InterpOrder_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%ADI_DT_C, 'ADI_DT_C', 'ADI_DT_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%ADI_TMax_C, 'ADI_TMax_C', 'ADI_TMax_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%ADI_storeHHVel, 'ADI_storeHHVel', 'ADI_storeHHVel', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%ADI_WrVTK_in, 'ADI_WrVTK_in', 'ADI_WrVTK_in', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%ADI_WrVTK_inType, 'ADI_WrVTK_inType', 'ADI_WrVTK_inType', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%ADI_WrVTK_inDT, 'ADI_WrVTK_inDT', 'ADI_WrVTK_inDT', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadAry( UnIn, FileName, InitInp%ADI_VTKNacDim_in, 6, 'ADI_VTKNacDim_in', 'ADI_VTKNacDim_in', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%ADI_VTKHubrad_in, 'ADI_VTKHubrad_in', 'ADI_VTKHubrad_in', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%ADI_wrOuts_C, 'ADI_wrOuts_C', 'ADI_wrOuts_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    CALL ReadVar( UnIn, FileName, InitInp%ADI_DT_Outs_C, 'ADI_DT_Outs_C', 'ADI_DT_Outs_C', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WaveTankTesting.ReadInput')
    IF (ErrStat >= AbortErrLev) RETURN

    if(UnIn>0) CLOSE( UnIn )

END SUBROUTINE

SUBROUTINE WaveTank_Init(   &
    WT_InputFile_C,         &
    MD_InputFile_C,         &
    SS_InputFile_C,         &
    AD_InputFile_C,         &
    IfW_InputFile_C,        &
    ErrStat_C,              &
    ErrMsg_C                &
) BIND (C, NAME='WaveTank_Init')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: WaveTank_Init
!GCC$ ATTRIBUTES DLLEXPORT :: WaveTank_Init
#endif

    TYPE(C_PTR),        INTENT(IN   ) :: WT_InputFile_C
    TYPE(C_PTR),        INTENT(IN   ) :: MD_InputFile_C
    TYPE(C_PTR),        INTENT(IN   ) :: SS_InputFile_C
    TYPE(C_PTR),        INTENT(IN   ) :: AD_InputFile_C
    TYPE(C_PTR),        INTENT(IN   ) :: IfW_InputFile_C
    INTEGER(C_INT),     INTENT(  OUT) :: ErrStat_C
    CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_C(ErrMsgLen_C)

    ! Local variables
    INTEGER(C_INT)                          :: ErrStat_C2
    CHARACTER(KIND=C_CHAR, LEN=ErrMsgLen_C) :: ErrMsg_C2
    INTEGER(IntKi)                          :: ErrStat_F2
    CHARACTER(ErrMsgLen)                    :: ErrMsg_F2
    TYPE(WaveTank_InitInput)                :: WT_InitInp
    CHARACTER(1024), POINTER                :: WT_InputFilePath

    ! The length of these arrays much match what is set in the corresponding C binding modules
    CHARACTER(KIND=C_CHAR) :: SS_OutputChannelNames_C(ChanLen*MaxOutPts+1)
    CHARACTER(KIND=C_CHAR) :: SS_OutputChannelUnits_C(ChanLen*MaxOutPts+1)
    CHARACTER(KIND=C_CHAR) :: MD_OutputChannelNames_C(100000)
    CHARACTER(KIND=C_CHAR) :: MD_OutputChannelUnits_C(100000)
    CHARACTER(KIND=C_CHAR) :: ADI_OutputChannelNames_C(ChanLen*MaxADIOutputs+1)
    CHARACTER(KIND=C_CHAR) :: ADI_OutputChannelUnits_C(ChanLen*MaxADIOutputs+1)

    ! Initialize error handling
    ErrStat_C = ErrID_None
    ErrMsg_C  = " "//C_NULL_CHAR

    CALL C_F_POINTER(WT_InputFile_C, WT_InputFilePath)
    call ReadInput(WT_InputFilePath, WT_InitInp, ErrStat_F2, ErrMsg_F2)
    CALL SetErrStat_F2C(ErrStat_F2, ErrMsg_F2, ErrStat_C, ErrMsg_C) !, 'WaveTank_Init')
    IF (ErrStat_C >= AbortErrLev) RETURN

    CALL SeaSt_C_Init(                          &    
        SS_InputFile_C,                         &
        WT_InitInp%SS_OutRootName_C,            &
        WT_InitInp%SS_Gravity_C,                &
        WT_InitInp%SS_WtrDens_C,                &
        WT_InitInp%SS_WtrDpth_C,                &
        WT_InitInp%SS_MSL2SWL_C,                &
        WT_InitInp%SS_NSteps_C,                 &
        WT_InitInp%SS_TimeInterval_C,           &
        WT_InitInp%SS_WaveElevSeriesFlag_C,     &
        WT_InitInp%SS_WrWvKinMod_C,             &
        SS_NumChannels_C,                       &
        SS_OutputChannelNames_C,                &
        SS_OutputChannelUnits_C,                &
        ErrStat_C2, ErrMsg_C2                   &
    )
    CALL SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'SeaSt_C_Init')
    IF (ErrStat_C >= AbortErrLev) RETURN

    ! Set the SeaState Wave Field pointer onto MoorDyn
    CALL WaveTank_SetWaveFieldPointer(ErrStat_C2, ErrMsg_C2)

    CALL MD_C_Init(                             &
        0,                                      &   !< InputFilePassed: 0 for file, 1 for string
        MD_InputFile_C,                         &
        IntfStrLen,                             &   !< InputFileStringLength_C
        WT_InitInp%MD_DT_C,                     &
        WT_InitInp%MD_G_C,                      &
        WT_InitInp%MD_RHO_C,                    &
        WT_InitInp%MD_DEPTH_C,                  &
        WT_InitInp%MD_PtfmInit_C,               &
        WT_InitInp%MD_InterpOrder_C,            &
        MD_NumChannels_C,                       &
        MD_OutputChannelNames_C,                &
        MD_OutputChannelUnits_C,                &
        ErrStat_C2,                             &
        ErrMsg_C2                               &
    )
    CALL SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'MD_C_Init')
    IF (ErrStat_C >= AbortErrLev) RETURN

    CALL ADI_C_PreInit(                         &
        WT_InitInp%NumTurbines_C,               &
        WT_InitInp%TransposeDCM,                &
        WT_InitInp%PointLoadOutput,             &
        WT_InitInp%ADI_gravity_C,               &
        WT_InitInp%ADI_defFldDens_C,            &
        WT_InitInp%ADI_defKinVisc_C,            &
        WT_InitInp%ADI_defSpdSound_C,           &
        WT_InitInp%ADI_defPatm_C,               &
        WT_InitInp%ADI_defPvap_C,               &
        WT_InitInp%ADI_WtrDpth_C,               &
        WT_InitInp%ADI_MSL2SWL_C,               &
        WT_InitInp%MHK,                         &
        WT_InitInp%DebugLevel,                  &
        ErrStat_C2, ErrMsg_C2                   &
    )
    CALL SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'ADI_C_PreInit')
    IF (ErrStat_C >= AbortErrLev) RETURN

    CALL ADI_C_SetupRotor(                      &
        WT_InitInp%iWT_c,                       &
        WT_InitInp%TurbineIsHAWT_c,             &
        WT_InitInp%TurbOrigin_C,                &
        WT_InitInp%HubPos_C,                    &
        WT_InitInp%HubOri_C,                    &
        WT_InitInp%NacPos_C,                    &
        WT_InitInp%NacOri_C,                    &
        WT_InitInp%NumBlades_C,                 &
        WT_InitInp%BldRootPos_C,                &
        WT_InitInp%BldRootOri_C,                &
        WT_InitInp%NumMeshPts_C,                &
        WT_InitInp%InitMeshPos_C,               &
        WT_InitInp%InitMeshOri_C,               &
        WT_InitInp%MeshPtToBladeNum_C,          &
        ErrStat_C2, ErrMsg_C2                   &
    )
    CALL SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'ADI_C_SetupRotor')
    IF (ErrStat_C >= AbortErrLev) RETURN

    CALL ADI_C_Init(                            &
        0,                                      &   !< ADinputFilePassed; 0 for file, 1 for string
        AD_InputFile_C,                         &   !< ADinputFileString_C; Input file as a single string with lines delineated by C_NULL_CHAR
        IntfStrLen,                             &   !< ADinputFileStringLength_C; length of the input file string
        0,                                      &   !< IfWinputFilePassed; 0 for file, 1 for string
        IfW_InputFile_C,                        &   !< IfWinputFileString_C; Input file as a single string with lines delineated by C_NULL_CHAR
        IntfStrLen,                             &   !< IfWinputFileStringLength_C; length of the input file string
        WT_InitInp%ADI_OutRootName_C,           &   !< Root name to use for echo files and other
        WT_InitInp%ADI_OutVTKDir_C,             &   !< Directory to put all vtk output
        WT_InitInp%ADI_InterpOrder_C,           &   !< Interpolation order to use (must be 1 or 2)
        WT_InitInp%ADI_DT_C,                    &   !< Timestep used with AD for stepping forward from t to t+dt.  Must be constant.
        WT_InitInp%ADI_TMax_C,                  &   !< Maximum time for simulation
        WT_InitInp%ADI_storeHHVel,              &   !< Store hub height time series from IfW
        WT_InitInp%ADI_WrVTK_in,                &   !< Write VTK outputs [0: none, 1: init only, 2: animation]
        WT_InitInp%ADI_WrVTK_inType,            &   !< Write VTK outputs as [1: surface, 2: lines, 3: both]
        WT_InitInp%ADI_WrVTK_inDT,              &   !< Timestep between VTK writes
        WT_InitInp%ADI_VTKNacDim_in,            &   !< Nacelle dimension passed in for VTK surface rendering [0,y0,z0,Lx,Ly,Lz] (m)
        WT_InitInp%ADI_VTKHubRad_in,            &   !< Hub radius for VTK surface rendering
        WT_InitInp%ADI_wrOuts_C,                &
        WT_InitInp%ADI_DT_Outs_C,               &
        ADI_NumChannels_C,           &
        ADI_OutputChannelNames_C,    &
        ADI_OutputChannelUnits_C,    &
        ErrStat_C2, ErrMsg_C2        &
    )
    CALL SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'ADI_C_Init')
    IF (ErrStat_C >= AbortErrLev) RETURN

END SUBROUTINE WaveTank_Init

SUBROUTINE WaveTank_CalcOutput( &
    time,                       &
    n_camera_points,            &
    positions_x,                &
    positions_y,                &
    positions_z,                &
    rotation_matrix,            &
    loads,                      &
    md_outputs,                 &
    adi_outputs,                &
    ErrStat_C,                  &
    ErrMsg_C                    &
) BIND (C, NAME='WaveTank_CalcOutput')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: WaveTank_CalcOutput
!GCC$ ATTRIBUTES DLLEXPORT :: WaveTank_CalcOutput
#endif

    REAL(C_DOUBLE),         INTENT(IN   ) :: time
    INTEGER(C_INT),         INTENT(IN   ) :: n_camera_points
    REAL(C_FLOAT),          INTENT(IN   ) :: positions_x(n_camera_points)
    REAL(C_FLOAT),          INTENT(IN   ) :: positions_y(n_camera_points)
    REAL(C_FLOAT),          INTENT(IN   ) :: positions_z(n_camera_points)
    REAL(C_FLOAT),          INTENT(IN   ) :: rotation_matrix(9)
    REAL(C_FLOAT),          INTENT(  OUT) :: loads(n_camera_points,6)
    REAL(C_FLOAT),          INTENT(  OUT) :: md_outputs(MD_NumChannels_C)
    REAL(C_FLOAT),          INTENT(  OUT) :: adi_outputs(ADI_NumChannels_C)
    INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_C
    CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_C(ErrMsgLen_C)

    ! Local variables
    INTEGER(C_INT)                          :: ErrStat_C2
    CHARACTER(KIND=C_CHAR, LEN=ErrMsgLen_C) :: ErrMsg_C2


    ! ! ADI
    ! ! SetRotorMotion
    integer(c_int)     :: ADI_iWT_c                         !< Wind turbine / rotor number
    real(c_float)      :: ADI_HubPos_C( 3 )                 !< Hub position
    real(c_double)     :: ADI_HubOri_C( 9 )                 !< Hub orientation
    real(c_float)      :: ADI_HubVel_C( 6 )                 !< Hub velocity
    real(c_float)      :: ADI_HubAcc_C( 6 )                 !< Hub acceleration
    real(c_float)      :: ADI_NacPos_C( 3 )                 !< Nacelle position
    real(c_double)     :: ADI_NacOri_C( 9 )                 !< Nacelle orientation
    real(c_float)      :: ADI_NacVel_C( 6 )                 !< Nacelle velocity
    real(c_float)      :: ADI_NacAcc_C( 6 )                 !< Nacelle acceleration
    integer(c_int), parameter :: NumBlades_C = 2                        !< Number of blades
    real(c_float)      :: ADI_BldRootPos_C( 3*NumBlades_C )   !< Blade root positions
    real(c_double)     :: ADI_BldRootOri_C( 9*NumBlades_C )   !< Blade root orientations
    real(c_float)      :: ADI_BldRootVel_C( 6*NumBlades_C )   !< Blade root velocities
    real(c_float)      :: ADI_BldRootAcc_C( 6*NumBlades_C )   !< Blade root accelerations
    ! Blade mesh nodes
    integer(c_int)     :: ADI_NumMeshPts_C                  !< Number of mesh points we are transfering motions to and output loads to
    real(c_float)      :: ADI_MeshPos_C( 3*n_camera_points )   !< A 3xNumMeshPts_C array [x,y,z]
    real(c_double)     :: ADI_MeshOri_C( 9*n_camera_points )   !< A 9xNumMeshPts_C array [r11,r12,r13,r21,r22,r23,r31,r32,r33]
    real(c_float)      :: ADI_MeshVel_C( 6*n_camera_points )   !< A 6xNumMeshPts_C array [x,y,z]
    real(c_float)      :: ADI_MeshAcc_C( 6*n_camera_points )   !< A 6xNumMeshPts_C array [x,y,z]


    ! Initialize error handling
    ErrStat_C = ErrID_None
    ErrMsg_C  = " "//C_NULL_CHAR

    ! Shift the positions and velocities over one index
    Positions(1,:) = Positions(2,:)
    Positions(2,:) = Positions(3,:)
    Velocities(1,:) = Velocities(2,:)

    ! Load the new positions
    Positions(3,:) = (/ positions_x, positions_y, positions_z, 0.0_C_FLOAT, 0.0_C_FLOAT, 0.0_C_FLOAT /)

    ! Calculate velocities and acceleration
    Velocities(1,:) = (/ (Positions(2,1:3) - Positions(1,1:3)) / real(0.1, c_float), 0.0_C_FLOAT, 0.0_C_FLOAT, 0.0_C_FLOAT /)
    Velocities(2,:) = (/ (Positions(3,1:3) - Positions(2,1:3)) / real(0.1, c_float), 0.0_C_FLOAT, 0.0_C_FLOAT, 0.0_C_FLOAT /)
    Accelerations(1,:) = (/ (Velocities(2,1:3) - Velocities(1,1:3)) / real(0.1, c_float), 0.0_C_FLOAT, 0.0_C_FLOAT, 0.0_C_FLOAT /)

    ! Get loads from MoorDyn
    CALL MD_C_UpdateStates(                 &
        time,                               &
        REAL(time + 0.1, C_DOUBLE),         &
        Positions(3,:),                     &
        Velocities(2,:),                    &
        Accelerations(1,:),                 &
        ErrStat_C2, ErrMsg_C2               &
    )
    CALL SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'MD_C_UpdateStates')
    IF (ErrStat_C >= AbortErrLev) RETURN

    CALL MD_C_CalcOutput(                   &
        time,                               &
        Positions(3,:),                     &
        Velocities(2,:),                    &
        Accelerations(1,:),                 &
        loads,                              &
        md_outputs,                         &
        ErrStat_C2, ErrMsg_C2               &
    )
    CALL SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'MD_C_CalcOutput')
    IF (ErrStat_C >= AbortErrLev) RETURN

    ! Get loads from ADI
    ADI_iWT_c = 1
    ADI_HubPos_C = Positions(3,1:3)
    ADI_HubOri_C = (/ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 /)
    ADI_HubVel_C = Velocities(2,:)
    ADI_HubAcc_C = Accelerations(1,:)
    ADI_NacPos_C = Positions(3,1:3)
    ADI_NacOri_C = (/ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 /)
    ADI_NacVel_C = Velocities(2,:)
    ADI_NacAcc_C = Accelerations(1,:)
    ADI_BldRootPos_C = Positions(3,:)
    ADI_BldRootOri_C = (/ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 /)
    ADI_BldRootVel_C = (/ Velocities(2,1:6), Velocities(2,1:6) /)
    ADI_BldRootAcc_C = (/ Accelerations(1,1:6), Accelerations(1,1:6) /)
    ADI_NumMeshPts_C = n_camera_points
    ADI_MeshPos_C = Positions(3,:)
    ADI_MeshOri_C = (/ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 /)
    ADI_MeshVel_C = Velocities(2,:)
    ADI_MeshAcc_C = Accelerations(1,:)

    CALL ADI_C_SetRotorMotion(              &
        ADI_iWT_c,                          & ! integer(c_int),  intent(in   )  :: iWT_c                         !< Wind turbine / rotor number
        ADI_HubPos_C,                       & ! real(c_float),   intent(in   )  :: HubPos_C( 3 )                 !< Hub position
        ADI_HubOri_C,                       & ! real(c_double),  intent(in   )  :: HubOri_C( 9 )                 !< Hub orientation
        ADI_HubVel_C,                       & ! real(c_float),   intent(in   )  :: HubVel_C( 6 )                 !< Hub velocity
        ADI_HubAcc_C,                       & ! real(c_float),   intent(in   )  :: HubAcc_C( 6 )                 !< Hub acceleration
        ADI_NacPos_C,                       & ! real(c_float),   intent(in   )  :: NacPos_C( 3 )                 !< Nacelle position
        ADI_NacOri_C,                       & ! real(c_double),  intent(in   )  :: NacOri_C( 9 )                 !< Nacelle orientation
        ADI_NacVel_C,                       & ! real(c_float),   intent(in   )  :: NacVel_C( 6 )                 !< Nacelle velocity
        ADI_NacAcc_C,                       & ! real(c_float),   intent(in   )  :: NacAcc_C( 6 )                 !< Nacelle acceleration
        ADI_BldRootPos_C,                   & ! real(c_float),   intent(in   )  :: BldRootPos_C( 3*Sim%WT(iWT_c)%NumBlades )   !< Blade root positions
        ADI_BldRootOri_C,                   & ! real(c_double),  intent(in   )  :: BldRootOri_C( 9*Sim%WT(iWT_c)%NumBlades )   !< Blade root orientations
        ADI_BldRootVel_C,                   & ! real(c_float),   intent(in   )  :: BldRootVel_C( 6*Sim%WT(iWT_c)%NumBlades )   !< Blade root velocities
        ADI_BldRootAcc_C,                   & ! real(c_float),   intent(in   )  :: BldRootAcc_C( 6*Sim%WT(iWT_c)%NumBlades )   !< Blade root accelerations
        ADI_NumMeshPts_C,                   & ! intent(in   )  :: NumMeshPts_C                  !< Number of mesh points we are transfering motions to and output loads to
        ADI_MeshPos_C,                      & ! intent(in   )  :: MeshPos_C( 3*NumMeshPts_C )   !< A 3xNumMeshPts_C array [x,y,z]
        ADI_MeshOri_C,                      & ! intent(in   )  :: MeshOri_C( 9*NumMeshPts_C )   !< A 9xNumMeshPts_C array [r11,r12,r13,r21,r22,r23,r31,r32,r33]
        ADI_MeshVel_C,                      & ! intent(in   )  :: MeshVel_C( 6*NumMeshPts_C )   !< A 6xNumMeshPts_C array [x,y,z]
        ADI_MeshAcc_C,                      & ! intent(in   )  :: MeshAcc_C( 6*NumMeshPts_C )   !< A 6xNumMeshPts_C array [x,y,z]
        ErrStat_C2, ErrMsg_C2               &
    )
    CALL SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'ADI_C_SetRotorMotion')
    IF (ErrStat_C >= AbortErrLev) RETURN

    CALL ADI_C_UpdateStates(                &
        time,                               &
        REAL(time + 0.1, C_DOUBLE),         &
        ErrStat_C2, ErrMsg_C2               &
    )
    CALL SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'ADI_C_UpdateStates')
    IF (ErrStat_C >= AbortErrLev) RETURN

    CALL ADI_C_CalcOutput(                  &
        time,                               &
        adi_outputs,                        &
        ErrStat_C2, ErrMsg_C2               &
    )
    CALL SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'ADI_C_CalcOutput')
    IF (ErrStat_C >= AbortErrLev) RETURN

    ! CALL ADI_C_GetRotorLoads(               &
    !     ErrStat_C, ErrMsg_C                 &
    ! )
    ! CALL SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'MD_C_Init')
    ! IF (ErrStat_C >= AbortErrLev) RETURN

END SUBROUTINE

SUBROUTINE WaveTank_StepFunction( &
    frame_number,               &
    n_camera_points,            &
    positions_x,                &
    positions_y,                &
    positions_z,                &
    rotation_matrix,            &
    loads,                      &
    ErrStat_C,                  &
    ErrMsg_C                    &
) BIND (C, NAME='WaveTank_StepFunction')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: WaveTank_StepFunction
!GCC$ ATTRIBUTES DLLEXPORT :: WaveTank_StepFunction
#endif

    ! INTEGER(C_INT)                        :: delta_time
    INTEGER(C_INT),         INTENT(IN   ) :: frame_number
    INTEGER(C_INT),         INTENT(IN   ) :: n_camera_points
    REAL(C_FLOAT),          INTENT(IN   ) :: positions_x(n_camera_points)
    REAL(C_FLOAT),          INTENT(IN   ) :: positions_y(n_camera_points)
    REAL(C_FLOAT),          INTENT(IN   ) :: positions_z(n_camera_points)
    REAL(C_FLOAT),          INTENT(IN   ) :: rotation_matrix(9)
    REAL(C_FLOAT),          INTENT(  OUT) :: loads(n_camera_points)
    INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_C
    CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_C(ErrMsgLen_C)

    INTEGER(C_INT) :: load_period = 20 ! seconds

    ! Initialize error handling
    ErrStat_C = ErrID_None
    ErrMsg_C  = " "//C_NULL_CHAR

    IF ( MOD(frame_number / load_period, 2) == 0 ) THEN
        loads = -1.0_C_FLOAT
    ELSE
        loads = 1.0_C_FLOAT
    ENDIF

END SUBROUTINE

SUBROUTINE WaveTank_End(ErrStat_C, ErrMsg_C) bind (C, NAME="WaveTank_End")
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: WaveTank_End
!GCC$ ATTRIBUTES DLLEXPORT :: WaveTank_End
#endif

    INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_C
    CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_C(ErrMsgLen_C)

    ! Local variables
    INTEGER(C_INT)                          :: ErrStat_C2
    CHARACTER(KIND=C_CHAR, LEN=ErrMsgLen_C) :: ErrMsg_C2

    CALL MD_C_END(ErrStat_C, ErrMsg_C)
    CALL SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'MD_C_END')
    IF (ErrStat_C >= AbortErrLev) RETURN

    CALL SeaSt_C_END(ErrStat_C, ErrMsg_C)
    CALL SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'SeaSt_C_END')
    IF (ErrStat_C >= AbortErrLev) RETURN

    CALL ADI_C_END(ErrStat_C, ErrMsg_C)
    CALL SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'ADI_C_END')
    IF (ErrStat_C >= AbortErrLev) RETURN

END SUBROUTINE

SUBROUTINE WaveTank_SetWaveFieldPointer(ErrStat_C, ErrMsg_C) bind (C, NAME="WaveTank_SetWaveFieldPointer")
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: WaveTank_SetWaveFieldPointer
!GCC$ ATTRIBUTES DLLEXPORT :: WaveTank_SetWaveFieldPointer
#endif

    INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_C
    CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_C(ErrMsgLen_C)

    ! Local variables
    INTEGER(C_INT)                          :: ErrStat_C2
    CHARACTER(KIND=C_CHAR, LEN=ErrMsgLen_C) :: ErrMsg_C2

    ! Set the SeaState FlowField pointer onto MoorDyn
    TYPE(C_PTR) :: WaveFieldPointer
    TYPE(SeaSt_WaveFieldType), POINTER :: WaveFieldPointer_F

    ! Initialize error handling
    ErrStat_C = ErrID_None
    ErrMsg_C  = " "//C_NULL_CHAR

    WaveFieldPointer = SeaSt_GetWaveFieldPointer_C()

    CALL C_F_POINTER(WaveFieldPointer, WaveFieldPointer_F)
    ! Verify that the data in the WaveField pointer has been set
    IF (WaveFieldPointer_F%WtrDpth == 0) THEN
        ErrStat_C2 = ErrID_Fatal
        ErrMsg_C2 = "SeaState WaveFieldPointer is WtrDpth is 0.0, so it it probably not initialized."
        CALL SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'WaveTank_SetWaveFieldPointer')
        RETURN
    END IF

    CALL MD_C_SetWaveFieldData(WaveFieldPointer)

END SUBROUTINE

SUBROUTINE WaveTank_NoOp(ErrStat_C, ErrMsg_C) bind (C, NAME="WaveTank_NoOp")
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: WaveTank_NoOp
!GCC$ ATTRIBUTES DLLEXPORT :: WaveTank_NoOp
#endif

    INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_C
    CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_C(ErrMsgLen_C)

    ! Local variables
    INTEGER(C_INT)                          :: ErrStat_C2
    CHARACTER(KIND=C_CHAR, LEN=ErrMsgLen_C) :: ErrMsg_C2

    ! Initialize error handling
    ErrStat_C = ErrID_None
    ErrMsg_C  = " "//C_NULL_CHAR

    ! No operation
    ErrStat_C2 = ErrID_Info
    ErrMsg_C2 = "Hi Stephen - No op here."

    CALL SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'WaveTank_NoOp')

END SUBROUTINE

SUBROUTINE WaveTank_Sizes(SS_NumOuts, MD_NumOuts, ADI_NumOuts) bind (C, NAME="WaveTank_Sizes")
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: WaveTank_Sizes
!GCC$ ATTRIBUTES DLLEXPORT :: WaveTank_Sizes
#endif
    INTEGER(C_INT),         INTENT(  OUT) :: SS_NumOuts
    INTEGER(C_INT),         INTENT(  OUT) :: MD_NumOuts
    INTEGER(C_INT),         INTENT(  OUT) :: ADI_NumOuts

    SS_NumOuts = SS_NumChannels_C
    MD_NumOuts = MD_NumChannels_C
    ADI_NumOuts = ADI_NumChannels_C

END SUBROUTINE

END MODULE WaveTankTesting
