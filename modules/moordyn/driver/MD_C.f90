!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2021 National Renewable Energy Laboratory
! Author: Nicole Mendoza
!
! This file is part of MoorDyn.
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
MODULE MoorDynAPI

    USE ISO_C_BINDING
    USE MoorDyn
    USE MoorDyn_Types
    USE NWTC_Library

IMPLICIT NONE

PUBLIC :: MD_INIT_C
PUBLIC :: MD_UPDATESTATES_C
PUBLIC :: MD_CALCOUTPUT_C
PUBLIC :: MD_END_C

! Global Variables
TYPE(MD_InitInputType)                  :: InitInp     !< Input data for initialization routine
TYPE(MD_InputType), ALLOCATABLE         :: u(:)        !< An initial guess for the input; input mesh must be defined
TYPE(MD_ParameterType)                  :: p           !< Parameters
TYPE(MD_ContinuousStateType)            :: x           !< Initial continuous states
TYPE(MD_DiscreteStateType)              :: xd          !< Initial discrete states
TYPE(MD_ConstraintStateType)            :: z           !< Initial guess of the constraint states
TYPE(MD_OtherStateType)                 :: other       !< Initial other states
TYPE(MD_OutputType)                     :: y           !< Initial system outputs (outputs are not calculated; only the output mesh is initialized)
TYPE(MD_MiscVarType)                    :: m           !< Initial misc/optimization variables
TYPE(MD_InitOutputType)                 :: InitOutData !< Output for initialization routine

INTEGER(IntKi)                          :: N_Global    !< Global timestep

CONTAINS

!===============================================================================================================
!---------------------------------------------- MD INIT --------------------------------------------------------
!===============================================================================================================
SUBROUTINE MD_INIT_C(InputFileString_C, InputFileStringLength_C, DT_C, G_C, RHO_C, DEPTH_C, PtfmInit_C, NumChannels_C, OutputChannelNames_C, OutputChannelUnits_C, ErrStat_C, ErrMsg_C) BIND (C, NAME='MD_INIT_C')

    !TEMPORARY hack until Waves handling is finalized
    USE WAVES, only: WaveGrid_n, WaveGrid_x0, WaveGrid_y0, WaveGrid_dx, WaveGrid_dy, WaveGrid_nx, WaveGrid_ny, WaveGrid_nz

    TYPE(C_PTR)                                    , INTENT(IN   )   :: InputFileString_C        !< Input file as a single string with lines deliniated by C_NULL_CHAR
    INTEGER(C_INT)                                 , INTENT(IN   )   :: InputFileStringLength_C  !< length of the input file string
    REAL(C_DOUBLE)                                 , INTENT(IN   )   :: DT_C
    REAL(C_DOUBLE)                                 , INTENT(IN   )   :: G_C
    REAL(C_DOUBLE)                                 , INTENT(IN   )   :: RHO_C
    REAL(C_DOUBLE)                                 , INTENT(IN   )   :: DEPTH_C
    REAL(C_FLOAT)                                  , INTENT(IN   )   :: PtfmInit_C(6)
    INTEGER(C_INT)                                 , INTENT(  OUT)   :: NumChannels_C
    TYPE(C_PTR)                                    , INTENT(  OUT)   :: OutputChannelNames_C
    TYPE(C_PTR)                                    , INTENT(  OUT)   :: OutputChannelUnits_C
    INTEGER(C_INT)                                 , INTENT(  OUT)   :: ErrStat_C
    CHARACTER(KIND=C_CHAR)                         , INTENT(  OUT)   :: ErrMsg_C(1025)

    ! Local Variables
    CHARACTER(KIND=C_char, LEN=InputFileStringLength_C), POINTER     :: InputFileString          !< Input file as a single string with NULL chracter separating lines
    CHARACTER(CHANLEN+1), ALLOCATABLE, TARGET                        :: tmp_OutputChannelNames_C(:)
    CHARACTER(CHANLEN+1), ALLOCATABLE, TARGET                        :: tmp_OutputChannelUnits_C(:)
    REAL(DbKi)                                                       :: DTcoupling
    INTEGER(IntKi)                                                   :: ErrStat, ErrStat2
    CHARACTER(ErrMsgLen)                                             :: ErrMsg, ErrMsg2
    INTEGER                                                          :: I, J

    ! NOTE: Wave info will be handled differently in the future.  So the following is a temporary hack until that is finalized
    ! Hard coded for 10 wave steps.  Doesn't actually matter since it will get zeroed
    INTEGER(IntKi)                                                   :: NStepWave = 10

    ! Initialize ErrStat
    ErrStat = ErrID_None
    ErrMsg  = ""

    ! Get fortran pointer to C_NULL_CHAR deliniated input file as a string 
    CALL C_F_pointer(InputFileString_C, InputFileString)

    ! Store string-inputs as type FileInfoType
    CALL InitFileInfo(InputFileString, InitInp%PassedPrimaryInputData, ErrStat, ErrMsg)           
    IF (ErrStat .GE. AbortErrLev) THEN 
       PRINT *, "MD_INIT_C: Failed to convert main input file string to FileInfoType"
       PRINT *, ErrMsg
       RETURN
    END IF

    ! Set other inputs for calling MD_Init
    DTcoupling               = REAL(DT_C, DbKi)
    InitInp%FileName         = 'notUsed'
    InitInp%RootName         = 'MDroot'
    InitInp%UsePrimaryInputFile = .FALSE.

    ! Environment variables -- These should be passed in from C.
    InitInp%g                = REAL(G_C, DbKi)
    InitInp%rhoW             = REAL(RHO_C, DbKi)
    InitInp%WtrDepth         = REAL(DEPTH_C, DbKi)

    ! Platform position (x,y,z,Rx,Ry,Rz) -- where rotations are small angle assumption in radians.
    ! This data is used to set the CoupledKinematics mesh that will be used at each timestep call
    CALL AllocAry (InitInp%PtfmInit, 6, 'InitInp%PtfmInit', ErrStat2, ErrMsg2 ); IF (Failed()) RETURN
    DO I = 1,6
        InitInp%PtfmInit(I)  = REAL(PtfmInit_C(I),ReKi)
    END DO

    ! Wave INformation - THIS IS A SHORT TERM HACK
    ! Fake wave info -- completely still, with no dynamic pressure terms
    ! Set wave info to zeros -- assume 10 timesteps for now (doesn't really matter since it isn't getting used)
    CALL AllocAry ( InitInp%WaveVel  ,NStepWave, WaveGrid_n, 3, 'InitInp%WaveVel' , ErrStat2, ErrMsg2 );    IF (Failed()) RETURN
    CALL AllocAry ( InitInp%WaveAcc  ,NStepWave, WaveGrid_n, 3, 'InitInp%WaveAcc' , ErrStat2, ErrMsg2 );    IF (Failed()) RETURN
    CALL AllocAry ( InitInp%WavePDyn ,NStepWave, WaveGrid_n,    'InitInp%WavePDyn', ErrStat2, ErrMsg2 );    IF (Failed()) RETURN
    CALL AllocAry ( InitInp%WaveElev ,NStepWave, WaveGrid_n,    'InitInp%WaveElev', ErrStat2, ErrMsg2 );    IF (Failed()) RETURN
    CALL AllocAry ( InitInp%WaveTime ,NStepWave,                'InitInp%WaveTime', ErrStat2, ErrMsg2 );    IF (Failed()) RETURN
    DO i=1,NStepWave
       InitInp%WaveTime(i) = DTcoupling * REAL(i-1, DbKi)
    END DO
    InitInp%WaveVel          = 0.0_ReKi
    InitInp%WaveAcc          = 0.0_ReKi
    InitInp%WavePDyn         = 0.0_ReKi
    InitInp%WaveElev         = 0.0_ReKi

    allocate(u(2), STAT=ErrStat)
      if (ErrStat .GE. AbortErrLev) then
         ErrStat = ErrID_Fatal
         ErrMsg  = "MD_INIT_C: Could not allocate input"
         RETURN
      endif

    ! Call the main subroutine MD_Init
    CALL MD_Init(InitInp, u(1), p, x, xd, z, other, y, m, DTcoupling, InitOutData, ErrStat, ErrMsg)
    !FIXME: this may catch messages labelled as Info as fatal errors.  You probably don't want that.
    IF (ErrStat .GE. AbortErrLev) THEN
        PRINT *, "MD_INIT_C: Main MD_Init subroutine failed!"
        PRINT *, ErrMsg
        RETURN
    END IF

    ! Convert the outputs of MD_Init from Fortran to C
    ALLOCATE(tmp_OutputChannelNames_C(size(InitOutData%writeOutputHdr)))
    ALLOCATE(tmp_OutputChannelUnits_C(size(InitOutData%writeOutputUnt)))
    NumChannels_C = size(InitOutData%writeOutputHdr)
    PRINT *, 'MD_INIT_C: The number of output channels is ', NumChannels_C

    DO I = 1,NumChannels_C
        tmp_OutputChannelNames_C(I) = TRANSFER(InitOutData%writeOutputHdr(I)//C_NULL_CHAR, tmp_OutputChannelNames_C(I))
        tmp_OutputChannelUnits_C(I) = TRANSFER(InitOutData%writeOutputUnt(I)//C_NULL_CHAR, tmp_OutputChannelUnits_C(I))
    END DO
    OutputChannelNames_C = C_LOC(tmp_OutputChannelNames_C)
    OutputChannelUnits_C = C_LOC(tmp_OutputChannelUnits_C)

    ! Clean up variables and set up for IFW_CALCOUTPUT_C
    CALL MD_DestroyInitInput( InitInp, ErrStat, ErrMsg )
    IF (ErrStat .GE. AbortErrLev) THEN
        PRINT *, ErrMsg
        RETURN
    END IF

    CALL MD_DestroyInitOutput( InitOutData, ErrStat, ErrMsg )
    IF (ErrStat .GE. AbortErrLev) THEN
        PRINT *, ErrMsg
        RETURN
    END IF

    IF (ErrStat /= 0) THEN
        ErrStat_C = ErrID_Fatal
    ELSE
        ErrStat_C = ErrID_None
    END IF
    ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )

    PRINT*, "DONE WITH MD_INIT_C!"

CONTAINS

    SUBROUTINE Cleanup()
        ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )
        ErrStat_C = ErrStat
    END SUBROUTINE Cleanup

    logical FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'MD_Init_C')
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
    END FUNCTION Failed

END SUBROUTINE MD_INIT_C

!===============================================================================================================
!---------------------------------------------- MD UPDATE STATES -----------------------------------------------
!===============================================================================================================
SUBROUTINE MD_UPDATESTATES_C(TIME_C, TIME2_C, ErrStat_C, ErrMsg_C) BIND (C, NAME='MD_UPDATESTATES_C')

    REAL(C_DOUBLE)                                 , INTENT(IN   )   :: TIME_C, TIME2_C
    INTEGER(C_INT)                                 , INTENT(  OUT)   :: ErrStat_C
    CHARACTER(KIND=C_CHAR)                         , INTENT(  OUT)   :: ErrMsg_C

    ! Local Variables
    REAL(DbKi)                                                       :: t_array(2)
    INTEGER(IntKi)                                                   :: ErrStat, J, InterpOrder
    CHARACTER(ErrMsgLen)                                             :: ErrMsg

    ! Set up inputs to MD_UpdateStates
    t_array(1)  = REAL(TIME_C, DbKi)          ! t
    t_array(2)  = REAL(TIME2_C, DbKi)         ! t + dt
    N_Global    = 0_IntKi                     ! MOORDYN IS NOT CURRENTLY USING, BUT MAY CHANGE IN THE FUTURE
    InterpOrder = 1                           ! InterpOrder must be 1 (linear) or 2 (quadratic)

    allocate(u(InterpOrder+1), STAT=ErrStat)
    IF (ErrStat .GE. AbortErrLev) then
         ErrStat = ErrID_Fatal
         ErrMsg  = "MD_UPDATESTATES_C: Could not allocate input"
         RETURN
    END IF

    ! Call the main subroutine MD_UpdateStates
    CALL MD_UpdateStates( t_array(1), N_Global, u, t_array, p, x, xd, z, other, m, ErrStat, ErrMsg)
    IF (ErrStat .GE. AbortErrLev) THEN
        PRINT *, "MD_UPDATESTATES_C: Main MD_calcOutput subroutine failed!"
        PRINT *, ErrMsg
        RETURN
    END IF

    ! Convert the outputs of MD_calcOutput back to C
    IF (ErrStat /= 0) THEN
        ErrStat_C = ErrID_Fatal
    ELSE
        ErrStat_C = ErrID_None
    END IF
    ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )

    PRINT*, "DONE WITH MD_UPDATESTATES_C!"

END SUBROUTINE MD_UPDATESTATES_C

!===============================================================================================================
!---------------------------------------------- MD CALC OUTPUT -------------------------------------------------
!===============================================================================================================
SUBROUTINE MD_CALCOUTPUT_C(Time_C, ErrStat_C, ErrMsg_C) BIND (C, NAME='MD_CALCOUTPUT_C')

    REAL(C_DOUBLE)                                 , INTENT(IN   )   :: Time_C
    INTEGER(C_INT)                                 , INTENT(  OUT)   :: ErrStat_C
    CHARACTER(KIND=C_CHAR)                         , INTENT(  OUT)   :: ErrMsg_C

    ! Local Variables
    REAL(DbKi)                                                       :: t
    INTEGER(IntKi)                                                   :: ErrStat
    CHARACTER(ErrMsgLen)                                             :: ErrMsg

    ! Set up inputs to MD_CalcOutput
    t = REAL(Time_C, DbKi)
    allocate(u(2), STAT=ErrStat)
    if (ErrStat .GE. AbortErrLev) then
        ErrStat = ErrID_Fatal
        ErrMsg  = "MD_CALCOUTPUT_C: Could not allocate input"
        RETURN
    end if

    ! Call the main subroutine MD_CalcOutput
    CALL MD_CalcOutput( t, u(1), p, x, xd, z, other, y, m, ErrStat, ErrMsg )
    IF (ErrStat .GE. AbortErrLev) THEN
        PRINT *, "MD_CALCOUTPUT_C: Main MD_calcOutput subroutine failed!"
        PRINT *, ErrMsg
        RETURN
    END IF

    ! Convert the outputs of MD_calcOutput back to C
    IF (ErrStat /= 0) THEN
        ErrStat_C = ErrID_Fatal
    ELSE
        ErrStat_C = ErrID_None
    END IF
    ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )

    PRINT*, "DONE WITH MD_CALCOUTPUT_C!"

END SUBROUTINE MD_CALCOUTPUT_C

!===============================================================================================================
!----------------------------------------------- MD END --------------------------------------------------------
!===============================================================================================================
SUBROUTINE MD_END_C(ErrStat_C,ErrMsg_C) BIND (C, NAME='MD_END_C')

    INTEGER(C_INT)                , INTENT(  OUT)      :: ErrStat_C
    CHARACTER(KIND=C_CHAR)        , INTENT(  OUT)      :: ErrMsg_C

    ! Local variables
    INTEGER(IntKi)                                     :: ErrStat
    CHARACTER(ErrMsgLen)                               :: ErrMsg

    ! Set up inputs for MD_End
    allocate(u(2), STAT=ErrStat)
    if (ErrStat .GE. AbortErrLev) then
       ErrStat = ErrID_Fatal
       ErrMsg  = "MD_END_C: Could not allocate input"
       RETURN
    END IF

    ! Call the main subroutine MD_End
    CALL MD_End(u(1), p, x, xd, z, other, y, m, ErrStat , ErrMsg)
    IF (ErrStat .GE. AbortErrLev) THEN
        PRINT *, "MD_END_C: MD_End failed!"
        PRINT *, ErrMsg
        RETURN
    END IF

    ! Convert the outputs of MD_End from Fortran to C
    IF (ErrStat /= 0) THEN
        ErrStat_C = ErrID_Fatal
    ELSE
        ErrStat_C = ErrID_None
    END IF
    ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )

    PRINT*, "DONE WITH MD_END_C!"

END SUBROUTINE MD_END_C

END MODULE
