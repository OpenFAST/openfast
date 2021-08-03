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
TYPE(MD_InitInputType)                  :: InitInp             !< Input data for initialization routine
TYPE(MD_InputType), ALLOCATABLE         :: u(:)                !< An initial guess for the input; input mesh must be defined
TYPE(MD_ParameterType)                  :: p                   !< Parameters
TYPE(MD_ContinuousStateType)            :: x                   !< Initial continuous states
TYPE(MD_DiscreteStateType)              :: xd                  !< Initial discrete states
TYPE(MD_ConstraintStateType)            :: z                   !< Initial guess of the constraint states
TYPE(MD_OtherStateType)                 :: other               !< Initial other states
TYPE(MD_OutputType)                     :: y                   !< Initial system outputs (outputs are not calculated; only the output mesh is initialized)
TYPE(MD_MiscVarType)                    :: m                   !< Initial misc/optimization variables
TYPE(MD_InitOutputType)                 :: InitOutData         !< Output for initialization routine

! Time tracking
INTEGER(IntKi)                          :: N_Global=0          !< Global timestep. MOORDYN IS NOT CURRENTLY USING, BUT MAY CHANGE IN THE FUTURE
INTEGER(IntKi)                          :: InterpOrder=1       !< Interpolation order: must be 1 (linear) or 2 (quadratic)

!------------------------------------------------------------------------------------
!  Meshes for motions and loads
TYPE(MeshType)                          :: MD_MotionMesh       !< mesh for motions of external nodes
TYPE(MeshType)                          :: MD_LoadMesh         !< mesh for loads  for external nodes
TYPE(MeshType)                          :: MD_LoadMesh_tmp     !< mesh for loads  for external nodes 
TYPE(MeshMapType)                       :: Map_Motion_2_MD_WB  !< Mesh mapping between input motion mesh and WAMIT body(ies) mesh
TYPE(MeshMapType)                       :: Map_MD_WB_2_Load    !< Mesh mapping between MD output WAMIT body loads mesh and external nodes mesh

!  Motions input (so that we don't have to reallocate all the time)
REAL(ReKi)                              :: tmpPositions(6,1)    !< temp array.  Probably don't need this, but makes conversion from C clearer.
REAL(ReKi)                              :: tmpVelocities(6,1)   !< temp array.  Probably don't need this, but makes conversion from C clearer.
REAL(ReKi)                              :: tmpForces(6,1)       !< temp array.  Probably don't need this, but makes conversion to   C clearer.

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
    CHARACTER(KIND=C_CHAR)                         , INTENT(  OUT)   :: OutputChannelNames_C(100000)
    CHARACTER(KIND=C_CHAR)                         , INTENT(  OUT)   :: OutputChannelUnits_C(100000)
    INTEGER(C_INT)                                 , INTENT(  OUT)   :: ErrStat_C
    CHARACTER(KIND=C_CHAR)                         , INTENT(  OUT)   :: ErrMsg_C(1025)

    ! Local Variables
    CHARACTER(KIND=C_char, LEN=InputFileStringLength_C), POINTER     :: InputFileString          !< Input file as a single string with NULL chracter separating lines
    REAL(DbKi)                                                       :: DTcoupling
    INTEGER(IntKi)                                                   :: ErrStat, ErrStat2
    CHARACTER(ErrMsgLen)                                             :: ErrMsg, ErrMsg2
    INTEGER                                                          :: I, J, K

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
    CALL AllocAry (InitInp%PtfmInit, 6, 'InitInp%PtfmInit', ErrStat, ErrMsg ); IF (Failed()) RETURN
    DO I = 1,6
        InitInp%PtfmInit(I)  = REAL(PtfmInit_C(I),ReKi)
    END DO

    ! Wave Information - THIS IS A SHORT TERM HACK
    ! Fake wave info -- completely still, with no dynamic pressure terms
    ! Set wave info to zeros -- assume 10 timesteps for now (doesn't really matter since it isn't getting used)
    CALL AllocAry ( InitInp%WaveVel  ,NStepWave, WaveGrid_n, 3, 'InitInp%WaveVel' , ErrStat, ErrMsg );    IF (Failed()) RETURN
    CALL AllocAry ( InitInp%WaveAcc  ,NStepWave, WaveGrid_n, 3, 'InitInp%WaveAcc' , ErrStat, ErrMsg );    IF (Failed()) RETURN
    CALL AllocAry ( InitInp%WavePDyn ,NStepWave, WaveGrid_n,    'InitInp%WavePDyn', ErrStat, ErrMsg );    IF (Failed()) RETURN
    CALL AllocAry ( InitInp%WaveElev ,NStepWave, WaveGrid_n,    'InitInp%WaveElev', ErrStat, ErrMsg );    IF (Failed()) RETURN
    CALL AllocAry ( InitInp%WaveTime ,NStepWave,                'InitInp%WaveTime', ErrStat, ErrMsg );    IF (Failed()) RETURN
    DO i=1,NStepWave
       InitInp%WaveTime(i) = DTcoupling * REAL(i-1, DbKi)
    END DO
    InitInp%WaveVel          = 0.0_ReKi
    InitInp%WaveAcc          = 0.0_ReKi
    InitInp%WavePDyn         = 0.0_ReKi
    InitInp%WaveElev         = 0.0_ReKi

    ALLOCATE(u(2), STAT=ErrStat)
    IF (ErrStat .GE. AbortErrLev) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = "MD_INIT_C: Could not allocate input"
         RETURN
    END IF
    
    !-------------------------------------------------
    ! Call the main subroutine MD_Init
    !-------------------------------------------------
    CALL MD_Init(InitInp, u(1), p, x, xd, z, other, y, m, DTcoupling, InitOutData, ErrStat, ErrMsg)
    IF (ErrStat .GE. AbortErrLev) THEN
        PRINT *, "MD_INIT_C: Main MD_Init subroutine failed!"
        PRINT *, ErrMsg
        RETURN
    END IF

    !-------------------------------------------------
    !  Set output channel information for driver code
    !-------------------------------------------------

    ! Number of channels
    NumChannels_C = size(InitOutData%WriteOutputHdr)

    ! Transfer the output channel names and units to c_char arrays for returning
    k=1
    DO i=1,NumChannels_C
        DO j=1,ChanLen    ! max length of channel name.  Same for units
            OutputChannelNames_C(k)=InitOutData%WriteOutputHdr(i)(j:j)
            OutputChannelUnits_C(k)=InitOutData%WriteOutputUnt(i)(j:j)
            k=k+1
        END DO
    END DO

    ! null terminate the string
    OutputChannelNames_C(k) = C_NULL_CHAR
    OutputChannelUnits_C(k) = C_NULL_CHAR

    !-------------------------------------------------
    ! Clean up variables and set up for MD_CALCOUTPUT_C
    !------------------------------------------------- 
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

    LOGICAL FUNCTION Failed()
    CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'MD_Init_C')
    Failed =  ErrStat >= AbortErrLev
    IF (Failed) THEN
        ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )
        ErrStat_C = ErrStat
    END IF
    END FUNCTION Failed

END SUBROUTINE MD_INIT_C

!===============================================================================================================
!---------------------------------------------- MD UPDATE STATES -----------------------------------------------
!===============================================================================================================
SUBROUTINE MD_UPDATESTATES_C(TIME_C, TIME2_C, POSITIONS_C, VELOCITIES_C, ErrStat_C, ErrMsg_C) BIND (C, NAME='MD_UPDATESTATES_C')

    REAL(C_DOUBLE)                                 , INTENT(IN   )   :: TIME_C, TIME2_C
    REAL(C_FLOAT)                                  , INTENT(IN   )   :: POSITIONS_C(1,6)
    REAL(C_FLOAT)                                  , INTENT(IN   )   :: VELOCITIES_C(1,6)
    INTEGER(C_INT)                                 , INTENT(  OUT)   :: ErrStat_C
    CHARACTER(KIND=C_CHAR)                         , INTENT(  OUT)   :: ErrMsg_C

    ! Local Variables
    REAL(DbKi)                                                       :: t_array(2)
    INTEGER(IntKi)                                                   :: ErrStat, ErrStat2, J
    CHARACTER(ErrMsgLen)                                             :: ErrMsg, ErrMsg2

    ! Set up inputs to MD_UpdateStates
    t_array(1)  = REAL(TIME_C, DbKi)          ! t
    t_array(2)  = REAL(TIME2_C, DbKi)         ! t + dt

    ALLOCATE(u(InterpOrder+1), STAT=ErrStat)
    IF (ErrStat .GE. AbortErrLev) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = "MD_UPDATESTATES_C: Could not allocate input"
         RETURN
    END IF

    ! Reshape position and velocity (from a row vector to a column vector)
    DO J = 1,6
        tmpPositions(J,1)  = REAL(POSITIONS_C(1,J),ReKi)
        tmpVelocities(J,1) = REAL(VELOCITIES_C(1,J),ReKi)
    END DO

    ! Transfer motions to input meshes
   call Set_MotionMesh()                                       ! update motion mesh with input motion arrays
   call MD_SetInputMotion( u(1), ErrStat2, ErrMsg2 )  ! transfer input motion mesh to u(1) meshes
      if (Failed())  return

    !-------------------------------------------------
    ! Call the main subroutine MD_UpdateStates
    !-------------------------------------------------
    CALL MD_UpdateStates( t_array(1), N_Global, u, t_array, p, x, xd, z, other, m, ErrStat, ErrMsg)
    IF (ErrStat .GE. AbortErrLev) THEN
        PRINT *, "MD_UPDATESTATES_C: Main MD_UpdateStates subroutine failed!"
        PRINT *, ErrMsg
        RETURN
    END IF

    !-------------------------------------------------
    ! Convert the outputs of MD_UpdateStates back to C
    !-------------------------------------------------
    IF (ErrStat .GE. 4) THEN
        ErrStat_C = ErrID_Fatal
    ELSE
        ErrStat_C = ErrID_None
    END IF
    ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )

    PRINT*, "DONE WITH MD_UPDATESTATES_C!"

CONTAINS

    LOGICAL FUNCTION Failed()
    CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'MD_Init_C')
    Failed =  ErrStat >= AbortErrLev
    IF (Failed) THEN
        ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )
        ErrStat_C = ErrStat
    END IF
    END FUNCTION Failed

END SUBROUTINE MD_UPDATESTATES_C

!===============================================================================================================
!---------------------------------------------- MD CALC OUTPUT -------------------------------------------------
!===============================================================================================================
SUBROUTINE MD_CALCOUTPUT_C(Time_C, POSITIONS_C, VELOCITIES_C, FORCES_C, OUTPUTS_C, ErrStat_C, ErrMsg_C) BIND (C, NAME='MD_CALCOUTPUT_C')

    REAL(C_DOUBLE)                                 , INTENT(IN   )   :: Time_C
    REAL(C_FLOAT)                                  , INTENT(IN   )   :: POSITIONS_C(1,6)
    REAL(C_FLOAT)                                  , INTENT(IN   )   :: VELOCITIES_C(1,6)
    REAL(C_FLOAT)                                  , INTENT(  OUT)   :: FORCES_C(1,6)
    REAL(C_FLOAT)                                  , INTENT(  OUT)   :: OUTPUTS_C(p%NumOuts)
    INTEGER(C_INT)                                 , INTENT(  OUT)   :: ErrStat_C
    CHARACTER(KIND=C_CHAR)                         , INTENT(  OUT)   :: ErrMsg_C

    ! Local Variables
    REAL(DbKi)                                                       :: t
    INTEGER(IntKi)                                                   :: ErrStat, ErrStat2, J
    CHARACTER(ErrMsgLen)                                             :: ErrMsg, ErrMsg2

    PRINT*, 'inside MD_CALCOUTPUT_C'

    ! Set up inputs to MD_CalcOutput
    t = REAL(Time_C, DbKi)
    ALLOCATE(u(2), STAT=ErrStat)
    if (ErrStat .GE. AbortErrLev) then
        ErrStat = ErrID_Fatal
        ErrMsg  = "MD_CALCOUTPUT_C: Could not allocate input"
        RETURN
    end if

    ! Reshape position and velocity (from row vector to a column vector)
    DO J = 1,6
        tmpPositions(J,1)    = REAL(POSITIONS_C(1,J),ReKi)
        tmpVelocities(J,1)   = REAL(VELOCITIES_C(1,J),ReKi)
    END DO

    ! Transfer motions to input meshes
    CALL Set_MotionMesh()                            ! update motion mesh with input motion arrays
    CALL MD_SetInputMotion( u(1), ErrStat, ErrMsg )  ! transfer input motion mesh to u(1) meshes
    IF (Failed())  RETURN

    !-------------------------------------------------
    ! Call the main subroutine MD_CalcOutput
    !-------------------------------------------------
    CALL MD_CalcOutput( t, u(1), p, x, xd, z, other, y, m, ErrStat, ErrMsg )
    IF (ErrStat .GE. AbortErrLev) THEN
        PRINT *, "MD_CALCOUTPUT_C: Main MD_calcOutput subroutine failed!"
        PRINT *, ErrMsg
        RETURN
    END IF

    !-------------------------------------------------
    ! Convert the outputs of MD_calcOutput back to C
    !-------------------------------------------------
    ! Transfer resulting load meshes to intermediate mesh
    CALL MD_TransferLoads( u(1), y, ErrStat, ErrMsg )
    if (Failed())  return

    ! Set output force/moment array
    CALL Set_OutputLoadArray( )
    ! Reshape for return
    DO J = 1,6
        FORCES_C(1,J) = REAL(tmpForces(J,1), c_float)
    END DO

    OUTPUTS_C = REAL(y%WriteOutput, C_FLOAT)

    IF (ErrStat .GE. 4) THEN
        ErrStat_C = ErrID_Fatal
    ELSE
        ErrStat_C = ErrID_None
    END IF
    ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )

    PRINT*, "DONE WITH MD_CALCOUTPUT_C!"

CONTAINS

    LOGICAL FUNCTION Failed()
    CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'MD_Init_C')
    Failed =  ErrStat >= AbortErrLev
    IF (Failed) THEN
        ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )
        ErrStat_C = ErrStat
    END IF
    END FUNCTION Failed

END SUBROUTINE MD_CALCOUTPUT_C

!===============================================================================================================
!----------------------------------------------- MD END --------------------------------------------------------
!===============================================================================================================
SUBROUTINE MD_END_C(ErrStat_C,ErrMsg_C) BIND (C, NAME='MD_END_C')

    INTEGER(C_INT)                , INTENT(  OUT)      :: ErrStat_C
    CHARACTER(KIND=C_CHAR)        , INTENT(  OUT)      :: ErrMsg_C

    ! Local variables
    INTEGER(IntKi)                                     :: ErrStat, ErrStat2
    CHARACTER(ErrMsgLen)                               :: ErrMsg, ErrMsg2

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
    IF (ErrStat .GE. 4) THEN
        ErrStat_C = ErrID_Fatal
    ELSE
        ErrStat_C = ErrID_None
    END IF
    ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )

    PRINT*, "DONE WITH MD_END_C!"

CONTAINS

    LOGICAL FUNCTION Failed()
    CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'MD_Init_C')
    Failed =  ErrStat >= AbortErrLev
    IF (Failed) THEN
        ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )
        ErrStat_C = ErrStat
    END IF
    END FUNCTION Failed

END SUBROUTINE MD_END_C

!===============================================================================================================
!----------------------------------------- ADDITIONAL SUBROUTINES ----------------------------------------------
!===============================================================================================================

! SUBROUTINE SetMotionLoadsInterfaceMeshes
!---------------------------------------------------------------------------------------------------------------
!! This subroutine sets the interface meshes to map to the input motions to the MD
!! meshes for WAMIT.  This subroutine also sets the meshes for the output loads.
SUBROUTINE SetMotionLoadsInterfaceMeshes(ErrStat,ErrMsg)
    integer(IntKi),         intent(  out)  :: ErrStat    !< temporary error status
    character(ErrMsgLen),   intent(  out)  :: ErrMsg     !< temporary error message
    real(ReKi)                             :: InitPos(3)
    real(R8Ki)                             :: theta(3)
    real(R8Ki)                             :: Orient(3,3)
    !-------------------------------------------------------------
    ! Set the interface  meshes for motion inputs and loads output
    !-------------------------------------------------------------
    ! Motion mesh
    !     This point mesh may contain more than one point. Mapping will be used to map
    !     this to the input meshes for WAMIT and Morison.
    call MeshCreate(  MD_MotionMesh                       ,  &
                      IOS              = COMPONENT_INPUT  ,  &
                      Nnodes           = 1                ,  &
                      ErrStat          = ErrStat         ,  &
                      ErrMess          = ErrMsg          ,  &
                      TranslationDisp  = .TRUE.,    Orientation = .TRUE., &
                      TranslationVel   = .TRUE.,    RotationVel = .TRUE., &
                      TranslationAcc   = .TRUE.,    RotationAcc = .TRUE.  )
    if (ErrStat >= AbortErrLev) RETURN

    ! initial position and orientation of node
    InitPos  = tmpPositions(1:3,1)
    theta    = real(tmpPositions(4:6,1),DbKi)    ! convert ReKi to DbKi to avoid roundoff
    Orient   = EulerConstruct( theta )
    CALL MeshPositionNode(  MD_MotionMesh            , &
                            1                        , &
                            InitPos                  , &  ! position
                            ErrStat, ErrMsg        , &
                            Orient                     )  ! orientation
    IF (ErrStat >= AbortErrLev) RETURN
     
    CALL MeshConstructElement ( MD_MotionMesh, ELEMENT_POINT, ErrStat, ErrMsg, 1 )
    IF (ErrStat >= AbortErrLev) RETURN

    CALL MeshCommit ( MD_MotionMesh, ErrStat, ErrMsg )
    IF (ErrStat >= AbortErrLev) RETURN

    MD_MotionMesh%RemapFlag  = .TRUE.

    ! For checking the mesh, uncomment this.
    !     note: CU is is output unit (platform dependent).
    !call MeshPrintInfo( CU, MD_MotionMesh )
 
    !-------------------------------------------------------------
    ! Loads mesh
    !     This point mesh may contain more than one point. Mapping will be used to map
    !     the loads from output meshes for WAMIT.
    ! Output mesh for loads at each WAMIT body
    CALL MeshCopy( SrcMesh  = MD_MotionMesh      ,&
                   DestMesh = MD_LoadMesh        ,&
                   CtrlCode = MESH_SIBLING       ,&
                   IOS      = COMPONENT_OUTPUT   ,&
                   ErrStat  = ErrStat           ,&
                   ErrMess  = ErrMsg            ,&
                   Force    = .TRUE.             ,&
                   Moment   = .TRUE.             )
    IF (ErrStat >= AbortErrLev) return
    
    MD_LoadMesh%RemapFlag  = .TRUE.

    ! For checking the mesh, uncomment this.
    !     note: CU is is output unit (platform dependent).
    !call MeshPrintInfo( CU, MD_LoadMesh )

    !-------------------------------------------------------------
    ! Loads mesh
    !     This point mesh may contain more than one point. Mapping will be used to map
    !     the loads from output meshes for WAMIT.
    ! Output mesh for loads at each WAMIT body
    ! CALL MeshCopy( SrcMesh  = MD_LoadMesh        ,&
    !                DestMesh = MD_LoadMesh_tmp    ,&
    !                CtrlCode = MESH_COUSIN        ,&
    !                IOS      = COMPONENT_OUTPUT   ,&
    !                ErrStat  = ErrStat           ,&
    !                ErrMess  = ErrMsg            ,&
    !                Force    = .TRUE.             ,&
    !                Moment   = .TRUE.             )
    !    if (ErrStat >= AbortErrLev) return
    
    ! MD_LoadMesh_tmp%RemapFlag  = .TRUE.

    ! For checking the mesh, uncomment this.
    !     note: CU is is output unit (platform dependent).
    !call MeshPrintInfo( CU, MD_LoadMesh_tmp )

    !-------------------------------------------------------------
    ! Set the mapping meshes
    ! WAMIT - floating bodies using potential flow
    IF ( u(1)%CoupledKinematics%Committed ) THEN      ! input motions
       CALL MeshMapCreate( MD_MotionMesh, u(1)%CoupledKinematics, Map_Motion_2_MD_WB, ErrStat, ErrMsg )
          IF (ErrStat >= AbortErrLev) RETURN
    END IF
    IF (    y%CoupledLoads%Committed ) THEN           ! output loads
       CALL MeshMapCreate( y%CoupledLoads, MD_LoadMesh, Map_MD_WB_2_Load, ErrStat, ErrMsg )
          IF (ErrStat >= AbortErrLev) RETURN
        END IF

END SUBROUTINE SetMotionLoadsInterfaceMeshes

! SUBROUTINE Set_MotionMesh
!---------------------------------------------------------------------------------------------------------------
!> This routine is operating on module level data, hence few inputs
SUBROUTINE Set_MotionMesh()
    real(R8Ki)                                :: theta(3)
    real(R8Ki)                                :: Orient(3,3)
    ! Set mesh corresponding to input motions
       theta    = real(tmpPositions(4:6,1),DbKi)    ! convert ReKi to DbKi to avoid roundoff
       Orient   = EulerConstruct( theta )
       MD_MotionMesh%TranslationDisp(1:3,1) = tmpPositions(1:3,1) - MD_MotionMesh%Position(1:3,1)  ! relative displacement only
       MD_MotionMesh%Orientation(1:3,1:3,1) = Orient
       MD_MotionMesh%TranslationVel( 1:3,1) = tmpVelocities(1:3,1)
       MD_MotionMesh%RotationVel(    1:3,1) = tmpVelocities(4:6,1)
END SUBROUTINE Set_MotionMesh

! SUBROUTINE MD_SetInputMotion
!---------------------------------------------------------------------------------------------------------------
!> Map the motion of the intermediate input mesh over to the input meshes
!! This routine is operating on module level data, hence few inputs
SUBROUTINE MD_SetInputMotion( u_local, ErrStat, ErrMsg )
    TYPE(MD_InputType),        INTENT(inout)  :: u_local
    INTEGER(IntKi),            INTENT(  out)  :: ErrStat
    CHARACTER(ErrMsgLen),      INTENT(  out)  :: ErrMsg
    !  WAMIT mesh
    IF ( u_local%CoupledKinematics%Committed ) then
        CALL Transfer_Point_to_Point( MD_MotionMesh, u_local%CoupledKinematics, Map_Motion_2_MD_WB, ErrStat, ErrMsg )
        IF (ErrStat >= AbortErrLev) RETURN
    END IF

END SUBROUTINE MD_SetInputMotion

! SUBROUTINE MD_TransferLoads
!---------------------------------------------------------------------------------------------------------------
!> Map the loads of the output meshes to the intermediate output mesh.  Since
!! we are mapping two meshes over to a single one, we use an intermediate
!! temporary mesh -- prevents accidental overwrite of WAMIT loads on MD_LoadMesh with the 
!! mapping of the Morison loads. This routine is operating on module level data, hence few inputs
SUBROUTINE MD_TransferLoads( u_local, y_local, ErrStat, ErrMsg )
    type(MD_InputType),        intent(in   )  :: u_local           ! Only one input (probably at T)
    type(MD_OutputType),       intent(in   )  :: y_local     ! Only one input (probably at T)
    integer(IntKi),            intent(  out)  :: ErrStat
    character(ErrMsgLen),      intent(  out)  :: ErrMsg
 
    MD_LoadMesh%Force    = 0.0_ReKi
    MD_LoadMesh%Moment   = 0.0_ReKi
 
    !  WAMIT mesh
    IF ( y_local%CoupledLoads%Committed ) THEN
        MD_LoadMesh_tmp%Force    = 0.0_ReKi
        MD_LoadMesh_tmp%Moment   = 0.0_ReKi
        CALL Transfer_Point_to_Point( y_local%CoupledLoads, MD_LoadMesh_tmp, Map_MD_WB_2_Load, ErrStat, ErrMsg, u_local%CoupledKinematics, MD_MotionMesh )
        IF (ErrStat >= AbortErrLev)  RETURN
        MD_LoadMesh%Force    = MD_LoadMesh%Force  + MD_LoadMesh_tmp%Force
        MD_LoadMesh%Moment   = MD_LoadMesh%Moment + MD_LoadMesh_tmp%Moment
    END IF
END SUBROUTINE MD_TransferLoads

! SUBROUTINE Set_OutputLoadArray
!---------------------------------------------------------------------------------------------------------------
!> Transfer the loads from the load mesh to the temporary array for output
!! This routine is operating on module level data, hence few inputs
SUBROUTINE Set_OutputLoadArray()
    ! Set mesh corresponding to input motions
       tmpForces(1:3,1)   = MD_LoadMesh%Force (1:3,1)
       tmpForces(4:6,1)   = MD_LoadMesh%Moment(1:3,1)
END SUBROUTINE Set_OutputLoadArray

END MODULE
